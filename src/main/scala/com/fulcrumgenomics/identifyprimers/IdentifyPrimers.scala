/*
 * The MIT License
 *
 * Copyright (c) 2018 Fulcrum Genomics LLC
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

package com.fulcrumgenomics.identifyprimers

import java.io.File
import java.nio.file.{Files, Paths}
import java.util.concurrent.ForkJoinPool
import java.util.concurrent.atomic.AtomicLong

import com.fulcrumgenomics.FgBioDef.{FilePath, PathToBam}
import com.fulcrumgenomics.alignment.{Aligner, Alignment, AsyncAligner, AsyncKswAligner, AsyncScalaAligner, Mode => AlignmentMode}
import com.fulcrumgenomics.bam.api.{SamOrder, SamRecord, SamSource, SamWriter}
import com.fulcrumgenomics.bam.{Bams, Template}
import com.fulcrumgenomics.cmdline.{ClpGroups, FgBioTool}
import com.fulcrumgenomics.commons.CommonsDef._
import com.fulcrumgenomics.commons.io.PathUtil
import com.fulcrumgenomics.commons.util.{LazyLogging, SimpleCounter}
import com.fulcrumgenomics.sopt.cmdline.ValidationException
import com.fulcrumgenomics.sopt.{arg, clp}
import com.fulcrumgenomics.util.{Io, Metric, ProgressLogger}
import htsjdk.samtools.util.SequenceUtil

import scala.collection.mutable.ListBuffer
import scala.concurrent.duration.Duration
import scala.concurrent.{Await, ExecutionContext, Future}


@clp(group=ClpGroups.SamOrBam, description=
  """
    |Identifies primers that generate the reads in the input BAM file.  Takes data that is generated by single-plex or
    |multiplex targeted PCR amplification and identifies the most likely primer for both the R1 and R2 ends of each
    |template.
    |
    |## Primers
    |
    |The input primers file must be tab-separated, and contain a header, and then one row per primer:
    |
    |  * `pair_id` - the unique identifier for the primer pair
    |  * `primer_id` - the unique identifier for the primer in a primer pair
    |  * `sequence` - the DNA sequence of the primer as ordered, including degenerate bases.
    |  * `strand` -  if the primer maps to the positive strand of the genome or R/r/- if the primer maps to the negative
    |                strand
    |  * `ref_name` - the chromosome/contig
    |  * `start` - the start position (1-based)
    |  * `end` - the end position (1-based, inclusive)
    |
    |`pair_id` and `primer_id` may not have commas.
    |
    |The primer file must contain headers.  If mapping information is not available for a primer then
    |ref_name/start/end should be left empty.  E.g.:
    |
    |```
    |pair_id    primer_id sequence strand chr  start   end
    |mapped1    m1_fwd    GATTACA  +      chr1 1010873 1010894
    |mapped1    m1_rev    ACATTAG  -      chr1 1011118 1011137
    |unmapped1  u1_fwd    ACGTGTA  +
    |unmapped1  u1_rev    GGATACA  -
    |```
    |
    |Use of the --multi-primer-pairs option allows the same pair_id to be used for > 2 primers. In this case any
    |combination of primers with the same pair_id will be deemed valid (i.e. a valid primer pair).  This is useful when
    |searching for primer matches on the 3' end of the read.
    |
    |## Primer Matching
    |
    |A match is looked for in the following order:
    |1. Using alignment position.  Reads whose starts overlap primer positions are then aligned (ungapped) to the primer
    |   sequence, and only reads that match the primer sequence within acceptable bounds are accepted.
    |2. Using ungapped alignment search.  The start of the read is compared against all possible primer sequences with
    |   an ungapped alignment and a match triggered if the alignment is sufficiently high quality.
    |3. Using gapped alignment search.  The start of the read is compared against all possible primer sequences using a
    |   gapped alignment and a match triggered if the alignment is sufficiently high quality.
    |
    |Alignments are accepted/rejected based on the following criteria:
    |* Ungapped alignments are accepted if the alignment has too many mismatches (see `--max-mismatch-rate`).
    |* Gapped alignments are accepted if the alignment has too low of a score (see `--min-alignment-score-rate`).  The
    |  alignment score rate is calculated as the alignment score divided by the primer length.
    |
    |## Matching Primers on the 3' End
    |
    |The `--three-prime` option can be used to also search the 3' end of every read for a primer sequence.
    |
    |If the read itself has a primer match on the 5' end, all other primers on the opposite strand in with the same
    |`pair_id` are compared.  Additionally, the primer match for the 5' end of its mate is compared (if present). If
    |neither yields a primer, then all possible primer sequences are compared.
    |
    |### Unmapped Data or Primers without Mappings
    |
    |If no reference is given or mapping information is not available, matching using mapping locations is skipped.
    |
    |### Speeding Up Gapped Alignment
    |
    |Install [the ksw executable](https://github.com/nh13/ksw) manually or via conda (`conda install -c bioconda ksw`)
    |and supply the path to the `ksw` executable via `--ksw`.  This may help when many gapped alignments need to be
    |performed.
    |
    |## SAM Tags
    |
    |Each read will be annotated with SAM tags based on results of primer matching.
    |
    |The `f5`, `r5` describe the primer match for the forward and reverse strand primers matching the 5' end of the read
    |respectively.  The `f3`, `r3` describe the primer match for the forward and reverse strand primers matching the 3'
    |end of the read respectively (when `--three-prime` is used).  If no match was found, then the tag is set "none",
    |otherwise, a comma-delimited set of values are given as follows
    |
    |  * `pair_id` - the unique identifier for the primer pair for the given primer match
    |  * `primer_id` - the unique identifier for the primer in a primer pair for the given primer match
    |  * `<ref_name>:<start>-<end>` - the primer chromosome/contig, start (1-based), and end (1-based, inclusive)
    |  * `strand` - '+' if the forward primer, '-' otherwise
    |  * `read_num` - the read number of the match (1 or 2).
    |  * `match_type` - how the match was found; valid values are 'location', 'gapped', or 'ungapped'
    |  * `match_type_info` - additional information based on the `match_type`, containing one or more values
    |
    |The `match_type_info` for each `match_type` is as follow:
    |  * location:
    |    * `num_mismatches` - the number of mismatches between the primer and read sequence
    |  * ungapped:
    |    * `num_mismatches` - the number of mismatches between the primer and read sequence
    |    * `next_best` - the number of mismatches in the next best ungapped alignment (or "na" if none was found)
    |  * gapped:
    |    * `score` - the alignment score
    |    * `next_best` - the alignment score of then next best alignment
    |    * `start` - the offset from the start of the read where the primer match starts
    |    * `end` - the offset from the start of the read where the primer match ends
    |
    |Note: if primer matches are found on both the 5' and 3' ends of a read, the `pair_id` will be different if they
    |belong to different primer pairs.
  """)
class IdentifyPrimers
(@arg(flag='i', doc="Input BAM file.")  val input: PathToBam,
 @arg(flag='p', doc="File containing information about the primer pairs.") val primerPairs: FilePath,
 @arg(flag='o', doc="Output BAM file.") val output: PathToBam,
 @arg(flag='m', doc="Path prefix for the metrics files.") val metrics: PathPrefix,
 @arg(flag='S', doc="Match to primer locations +/- this many bases.") val slop: Int = 5,
 @arg(flag='3', doc="The # of extra bases in addition to the primer to search for primers on the 3' end (-1 to disable).") val threePrime: Int = -1,
 @arg(          doc="Maximum per-base mismatch rate for a primer to be considered a match.") val maxMismatchRate: Double = 0.05,
 @arg(          doc="The minimum per-base alignment score rate for gapped alignment for a primer to be considered a match. " +
   "This is the total score divided by the primer length.") val minAlignmentScoreRate: Double = 0.01,
 @arg(          doc="The match score to use for aligning to primer sequences (must be >= 0).") val matchScore: Int = 1,
 @arg(          doc="The mismatch score to use for aligning to primer sequences (must be <= 0).") val mismatchScore: Int = -4,
 @arg(          doc="The gap open score to use for aligning to primer sequences (must be <= 0).") val gapOpen: Int = -6,
 @arg(          doc="The gap extension score to use for aligning to primer sequences (must be <= 0).") val gapExtend: Int = -1,
 @arg(flag='t', doc="The number of threads to use.") val threads: Int = 1,
 @arg(          doc="Skip gapped-alignment matching") val skipGappedAlignment: Boolean = false,
 @arg(          doc="Path to the ksw aligner.") val ksw: Option[String] = None,
 @arg(flag='k', doc="Skip ungapped alignment if no kmer of this length is common between any primer and a given read.") val ungappedKmerLength: Option[Int] = None,
 @arg(flag='K', doc="Skip gapped alignment if no kmer of this length is common between any primer and a given read.") val gappedKmerLength: Option[Int] = None,
 @arg(          doc="Allow multiple primers on the same strand to have the same `pair_id`.") val multiPrimerPairs: Boolean = false
) extends FgBioTool with LazyLogging {

  /** The number of templates to process at a time per thread. */
  private val templatesPerThread: Int = 5000

  /** Finds the ksw executable, which may be on the system path. */
  private val kswExecutable: Option[FilePath] = this.ksw.map(Paths.get(_)).flatMap {
    case p if Files.exists(p)                  => Some(p)
    case name if name.contains(File.separator) =>
      throw new ValidationException(s"Path to the ksw executable does not exist: $ksw")
    case name                                  =>
      // try the system path
      val path = System.getenv("PATH")
      validate(path != null, s"Searching for the ksw executable '$ksw' on the PATH, but the PATH environment variable was not defined.")
      path.split(File.pathSeparatorChar)
        .view
        .map(PathUtil.pathTo(_))
        .map(p => p.resolve(name))
        .find(ex => Files.exists(ex))
        .orElse {
          throw new ValidationException(s"Is the path to the ksw executable mis-typed? Could not find ksw executable $ksw in PATH: $path")
        }
  }

  Io.assertReadable(input)
  Io.assertCanWriteFile(output)
  Io.assertCanWriteFile(metrics)
  kswExecutable.foreach(Io.assertReadable)

  validate(maxMismatchRate >= 0, "--max-mismatches must be >= 0")
  validate(matchScore >= 0,       "--match-score must be >= 0")
  validate(mismatchScore <= 0,    "--mismatch-score must be <= 0")
  validate(gapOpen <= 0,          "--gap-open must be <= 0")
  validate(gapExtend <= 0,        "--gap-extend must be <= 0")
  validate(minAlignmentScoreRate >= 0, "--min-alignment-score-rate must be >= 0")

  private val (locationBasedMatcher, ungappedBasedMatcher, gappedBasedMatcher, pairIdToPrimers) = {
    val primers  = Primer.read(this.primerPairs, multiPrimerPairs = multiPrimerPairs).toIndexedSeq
    logger.info(f"Read in ${primers.length}%,d primers.")
    val location = new LocationBasedPrimerMatcher(primers, slop, maxMismatchRate)
    val ungapped = new UngappedAlignmentBasedPrimerMatcher(primers, slop, maxMismatchRate, ungappedKmerLength)
    val gapped   = {
      // NB: `GappedAlignmentBasedPrimerMatcher` requires an aligner.  Since we do not call `.find()` on it, it will
      // never be used, so we give it null.
      new GappedAlignmentBasedPrimerMatcher(primers, slop, null, minAlignmentScoreRate, gappedKmerLength)
    }
    (location, ungapped, gapped, primers.groupBy(_.pair_id))
  }

  val numAlignments: AtomicLong = new AtomicLong(0)

  private val PrimerInfoForward5PrimeTag: String   = "f5"
  private val PrimerInfoReverse5PrimeTag: String   = "r5"
  private val PrimerInfoForward3PrimeTag: String   = "f3"
  private val PrimerInfoReverse3PrimeTag: String   = "r3"
  private val NoPrimerMatchInfo: String            = "none"
  private val comments: Seq[String]                = {
    val allTags = Seq(PrimerInfoForward5PrimeTag, PrimerInfoReverse5PrimeTag, PrimerInfoForward3PrimeTag, PrimerInfoReverse3PrimeTag).mkString("/")
    Seq(
      s"The $allTags tags store the primer match metadata for the forward and reverse strand respectively.",
      s"The $allTags tags are formatted as follow: <pair_id>,<primer_id>,<ref_name>:<start>-<end>,<strand>,<read-num>,<match-offset>,<match-length>,<match-type>,<match-type-info>.",
      s"The match-type is 'location', 'gapped', or 'ungapped' based on if the match was found using the location, mismatch-based (ungapped) alignment, or gapped-alignment."
    )
  }

  override def execute(): Unit = {
    this.kswExecutable match {
      case Some(p) => logger.info(s"Using ksw aligner: $p")
      case None    => logger.info("Using the internal aligner; install ksw for faster operation.")
    }
    val metricCounter = new SimpleCounter[TemplateTypeMetric]()

    // Group the reads by template
    val in: SamSource = SamSource(this.input)
    val iterator: Iterator[Template] = Bams.templateIterator(in)

    // NB: Add comments explaining tags in the output writer
    val (out: SamWriter, programGroupId: Option[String]) = {
      val header = in.header.clone()
      comments.foreach { comment => header.addComment(comment) }
      // FIXME: https://github.com/fulcrumgenomics/fgbio/issues/439
      //val pgId   = this.toolInfo.map { info => info.applyTo(header) }
      val writer = SamWriter(path = output, header = header, sort = SamOrder(in.header))
      //(writer, pgId)
      (writer, None) // FIXME
    }

    // We take the input records, and batch them into `majorBatchSize` number of records.  For each major-batch, we
    // split those into `templatesPerThread` sub-batches.  Each sub-batch is processed by a single thread.  We process
    // each major-batch serially, meaning we wait for one major-batch to complete before moving onto the next one.  This
    // is so we don't have to read all records into memory to parallelize.  We set the `majorBatchSize` to have more
    // sub-batches than just `templatesPerThread * threads` so that we can more efficiently utilize the available threads.
    val majorBatchSize  = templatesPerThread * threads * 4
    val readingProgress = ProgressLogger(this.logger, verb="read", unit=5e4.toInt)

    // Batch templates to process them in individual threads.
    val outputIterator: Iterator[Seq[Template]] = if (threads > 1) {
      logger.info(f"Batching $majorBatchSize%,d templates with $templatesPerThread%,d templates per thread.")
      val pool = new ForkJoinPool(threads, ForkJoinPool.defaultForkJoinWorkerThreadFactory, null, true)
      implicit val ec: ExecutionContext = ExecutionContext.fromExecutor(pool)
      import com.fulcrumgenomics.commons.CommonsDef.ParSupport
      // Developer Note: Iterator does not support parallel operations, so we need to group together records into a
      // [[List]] or [[Seq]].  A fixed number of records are grouped to reduce memory overhead.
      iterator
        .grouped(majorBatchSize)
        .flatMap { templates =>
          templates
            .grouped(templatesPerThread)
            .toStream
            .parWith(threads, fifo = false)
            .map { templates => processBatch(templates, metricCounter, readingProgress) }
            .seq // ensures that the only parallelism is processBatch
            .toIterator
        }
    }
    else {
      import scala.concurrent.ExecutionContext.Implicits.global
      logger.info(f"Batching $templatesPerThread%,d templates.")
      iterator
        .grouped(templatesPerThread)
        .map { templates => processBatch(templates, metricCounter, readingProgress) }
    }

    // Write the results
    val writingProgress = ProgressLogger(this.logger, "written", unit=1e6.toInt)
    outputIterator
      .flatMap(_.flatMap(_.allReads))
      .foreach { rec =>
        programGroupId.foreach { pgId => rec("PG") = pgId }
        writingProgress.record(rec)
        out += rec
      }

    val rate = numAlignments.get() / readingProgress.getElapsedSeconds.toDouble
    logger.info(f"Performed ${numAlignments.get()}%,d gapped alignments in total ($rate%,.2f alignments/second).")
    logger.info(f"Wrote ${readingProgress.getCount}%,d records.")

    in.safelyClose()
    out.close()

    // Create the [[TemplateTypeMetric]]
    val templateTypeMetrics = TemplateTypeMetric.metricsFrom(metricCounter)

    // Write metrics
    Metric.write(PathUtil.pathTo(this.metrics + ".detailed.txt"), templateTypeMetrics)
    Metric.write(PathUtil.pathTo(this.metrics + ".summary.txt"), IdentifyPrimersMetric(templateTypeMetrics))
  }

  /** Creates a new [[com.fulcrumgenomics.alignment.AsyncAligner]]. */
  private def newAligner[T <: AsyncAligner](mode: AlignmentMode = AlignmentMode.Glocal): AsyncAligner = {
    this.kswExecutable match {
      case Some(executable) => new AsyncKswAligner(executable, matchScore, mismatchScore, gapOpen, gapExtend, mode)
      case None             => new AsyncScalaAligner(Aligner(matchScore, mismatchScore, gapOpen, gapExtend, mode))
    }
  }

  /** Processes a single batch of templates.
    *
    * Attempts to match primers based on location, then ungapped alignment.  If no match was found, attempts to perform
    * gapped alignment, which is performed asynchronously.
    *
    * Returns the batch of templates, with reads annotated with primer matches if found.
    * */
  def processBatch(templates: Seq[Template],
                   metricCounter: SimpleCounter[TemplateTypeMetric],
                   progress: ProgressLogger)
                  (implicit e :ExecutionContext): Seq[Template] = {
    val fivePrimeAligner: AsyncAligner = newAligner()
    val threePrimerAligner: AsyncAligner = newAligner(AlignmentMode.Local) // NB: we run this in local mode to get partial matches

    val futures: Seq[Future[TemplateTypeMetric]] = templates.map { template =>
      // Get futures for performing primer matching on the 5' end of the reads
      val r1FivePrimeFuture = template.r1.map { r => toPrimeMatchFuture(r, fivePrimeAligner) }.getOrElse(Future.successful(None))
      val r2PrimerMatchFuture = template.r2.map { r => toPrimeMatchFuture(r, fivePrimeAligner) }.getOrElse(Future.successful(None))
      // Return a future with when all primer matching is complete
      r1FivePrimeFuture.flatMap { r1FivePrime =>
        r2PrimerMatchFuture.flatMap { r2FivePrime =>
          // Get futures for performing primer matching on the 3' end of the reads, which depend on the results of 5' primer matching
          val r1ThreePrimeFuture = template.r1.map { r1 => matchThreePrime(r1, r1FivePrime, r2FivePrime, threePrimerAligner) }.getOrElse(Future.successful(None))
          val r2ThreePrimeFuture = template.r2.map { r2 => matchThreePrime(r2, r2FivePrime, r1FivePrime, threePrimerAligner) }.getOrElse(Future.successful(None))
          r1ThreePrimeFuture.flatMap { r1ThreePrime =>
            r2ThreePrimeFuture.map { r2ThreePrime =>
              // All primer matching is complete, annotate the reads and update the counts

              val templateType = (template.r1, template.r2) match {
                case (Some(r1), Some(r2)) =>
                  addTagsToPair(r1, r2, r1FivePrime, r2FivePrime, r1ThreePrime, r2ThreePrime)
                  TemplateType(r1, Some(r2))
                case (Some(r1), None) =>
                  require(!r1.paired, s"Found paired read but missing R2 for ${r1.name}")
                  addTagsToFragment(r1, r1FivePrime, r1ThreePrime)
                  TemplateType(r1, None)
                case _ =>
                  throw new IllegalStateException(s"Template did not have an R1: ${template.name}")
              }
              TemplateTypeMetric(templateType, template.r1.get.isFrPair, r1FivePrime, r2FivePrime)
            }
          }
        }
      }
    }

    // Wait for primer matching to complete on all templates
    val metrics = Await.result(Future.sequence(futures), Duration.Inf)

    // Update counters and progress
    numAlignments.addAndGet(fivePrimeAligner.numAligned)
    metricCounter.synchronized {
      metrics.foreach { m => metricCounter.count(m) }
      templates.flatMap(_.allReads).foreach { rec =>
        if (progress.record(rec)) {
          val rate = numAlignments.get() / progress.getElapsedSeconds.toDouble
          logger.info(f"Performed ${numAlignments.get()}%,d gapped alignments so far ($rate%,.2f alignments/second).")
        }
      }
    }

    templates
  }


  /** Matches the read based on location, then ungapped alignment.  If no match was found, return a list of alignment
    * tasks. */
  private def toPrimeMatchFuture(rec: SamRecord, aligner: AsyncAligner)(implicit e :ExecutionContext): Future[Option[PrimerMatch]] = {
    locationBasedMatcher.find(rec).orElse { ungappedBasedMatcher.find(rec) } match {
      case Some(pm) => Future.successful(Some(pm))
      case None     =>
        val futures: Seq[Future[(Primer, Alignment)]] = {
          gappedBasedMatcher.toAlignmentTasks(rec).map { case (primer, target) =>
            aligner.align(primer.basesInSequencingOrder, target).map { alignment => (primer, alignment) }
          }
        }
        finalizedGappedAlignment(futures)
    }
  }

  /** Returns a future that when completed, gives the best primer match from the given alignments, otherwise None.
    *
    * The alignments are filtered by the minimum alignment score rate and to ensure that the end of the alignment is
    * within the slop distance of the end of the read.
    * */
  private def finalizedGappedAlignment(futures: Seq[Future[(Primer, Alignment)]])
                                      (implicit e :ExecutionContext): Future[Option[PrimerMatch]] = {
    Future.sequence(futures).map { primerAndAlignments: Seq[(Primer, Alignment)] =>

      // Get the best two alignments, sorted by score descending.
      gappedBasedMatcher.getBestGappedAlignment(primerAndAlignments).filter { gappedPrimerMatch =>
        // Filter by minAlignmentScoreRate
        // TODO: we could ensure that the end of the alignment is within the slop distance
        val alignmentScoreRate = gappedPrimerMatch.score / gappedPrimerMatch.primer.length.toDouble
        alignmentScoreRate >= minAlignmentScoreRate
      }
    }
  }

  /** Returns a future that performs matching on the 3' end of the read.
    *
    * The list of primers against which to match are determined as follows
    * 1. Use the primer found for the mate if given.
    * 2. Use the opposite primers in the same primer pair set if a 5' match was found.
    * 3. Use all primers if (1 & 2) are not valid.
    *
    * If the read was paired, then use the primer found for the mate.  Otherwise, if no 5' match, use all the primers,
    * otherwise, get the primers from the same primer pair set and return the ones on the opposite strand.
    * */
  private def matchThreePrime(rec: SamRecord,
                              fivePrimeMatch: Option[PrimerMatch],
                              otherReadFivePrimerMatch: Option[PrimerMatch],
                              aligner: AsyncAligner)
                             (implicit e :ExecutionContext): Future[Option[PrimerMatch]] = if (threePrime < 0) Future.successful(None) else {
    import PrimerMatcher.SequencingOrderBases

    val primersToCheck: Seq[Primer] = {
      val buffer = ListBuffer[Primer]()

      // Use the paired primer(s)
      fivePrimeMatch.foreach { primerMatch =>
        pairIdToPrimers(primerMatch.primer.pair_id)
          .filter(_.positive_strand == primerMatch.primer.negativeStrand)
          .foreach { primer => buffer.append(primer) }
      }

      // Use the primer on the mate's 5' end
      otherReadFivePrimerMatch.foreach { primerMatch => buffer.append(primerMatch.primer) }

      // If not primers matches were found on 5' end of the read, or on the mate, then try all primers
      if (buffer.nonEmpty) buffer.toList else this.gappedBasedMatcher.primers
    }

    // for 3' matching, it is the same as 5' matching, just that we match the end of the reads
    def toTarget(targetLength: Int): Array[Byte] = {
      val bases = rec.basesInSequencingOrder.clone()
      SequenceUtil.reverseComplement(bases)
      bases
    }
    // use a little cache so we don't create to many objects
    val targetCache = scala.collection.mutable.HashMap[Int, Array[Byte]]()
    val futures = primersToCheck.map { primer =>
      val targetLength = math.min(primer.sequence.length + threePrime + slop, rec.length)
      val target       = targetCache.getOrElseUpdate(targetLength, toTarget(targetLength))
      aligner.align(primer.basesInSequencingOrder, target).map { alignment => (primer, alignment) }
    }
    // wait for it...
    finalizedGappedAlignment(futures)
  }

  /** A little cse class to store some information about primer matches on a given record, used in [[addTagsToPair]]. */
  private case class PrimerMatches(rec: SamRecord, fivePrimeMatch: Option[PrimerMatch], threePrimeMatch: Option[PrimerMatch], positiveStrand: Boolean)

  /** Adds tags to paired end reads based on the primer matching results. */
  private def addTagsToPair(r1: SamRecord,
                            r2: SamRecord,
                            r1FivePrimeMatch: Option[PrimerMatch],
                            r2FivePrimeMatch: Option[PrimerMatch],
                            r1ThreePrimeMatch: Option[PrimerMatch],
                            r2ThreePrimeMatch: Option[PrimerMatch]): Unit = {
    // Get which rec is on the forward strand, and which is reverse.  First looks at the primer match for r1, then the
    // primer match for r2, otherwise r1 is forward.
    val (forwardPrimerMatch, reversePrimerMatch) = {
      val r1PrimerMatch = PrimerMatches(r1, r1FivePrimeMatch, r1ThreePrimeMatch, positiveStrand=true)
      val r2PrimerMatch = PrimerMatches(r2, r2FivePrimeMatch, r2ThreePrimeMatch, positiveStrand=false)
      (r1FivePrimeMatch, r2FivePrimeMatch) match {
        case (Some(m), _)    =>
          if (m.primer.positive_strand) (r1PrimerMatch, r2PrimerMatch) else  (r2PrimerMatch.copy(positiveStrand=true), r1PrimerMatch.copy(positiveStrand=false))
        case (None, Some(m)) =>
          if (m.primer.negativeStrand) (r1PrimerMatch, r2PrimerMatch) else  (r2PrimerMatch.copy(positiveStrand=true), r1PrimerMatch.copy(positiveStrand=false))
        case _               => (r1PrimerMatch, r2PrimerMatch)
      }
    }
    require(forwardPrimerMatch.positiveStrand && !reversePrimerMatch.positiveStrand)

    // Set the primer match info to place in the forward/reverse primer match tags
    def tagit(fwdTag: String, revTag: String, f: PrimerMatches => Option[PrimerMatch]): Unit = {
      val forwardInfo = f(forwardPrimerMatch).map(_.info(forwardPrimerMatch.rec)).getOrElse(NoPrimerMatchInfo)
      val reverseInfo = f(reversePrimerMatch).map(_.info(reversePrimerMatch.rec)).getOrElse(NoPrimerMatchInfo)
      Seq(r1, r2).foreach { rec =>
        rec(fwdTag) = forwardInfo
        rec(revTag) = reverseInfo
      }
    }
    tagit(PrimerInfoForward5PrimeTag, PrimerInfoReverse5PrimeTag, _.fivePrimeMatch)
    if (threePrime >= 0) tagit(PrimerInfoForward3PrimeTag, PrimerInfoReverse3PrimeTag, _.threePrimeMatch)
  }

  /** Adds tags to a fragment read based on the primer matching results. */
  private def addTagsToFragment(frag: SamRecord,
                                fragFivePrimeMatch: Option[PrimerMatch],
                                fragThreePrimeMatch: Option[PrimerMatch]): Unit = {
    val forwardInfo = fragFivePrimeMatch.map(_.info(frag)).getOrElse(NoPrimerMatchInfo)

    frag(PrimerInfoForward5PrimeTag) = forwardInfo
    frag(PrimerInfoReverse5PrimeTag) = NoPrimerMatchInfo

    if (threePrime >= 0) {
      frag(PrimerInfoForward5PrimeTag) = fragThreePrimeMatch.map(_.info(frag)).getOrElse(NoPrimerMatchInfo)
      frag(PrimerInfoReverse5PrimeTag) = NoPrimerMatchInfo
    }
  }
}
