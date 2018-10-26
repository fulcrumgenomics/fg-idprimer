[![Build Status](https://travis-ci.org/fulcrumgenomics/fg-idprimer.svg?branch=master)](https://travis-ci.org/fulcrumgenomics/fg-idprimer)
[![codecov](https://codecov.io/gh/fulcrumgenomics/fg-idprimer/branch/master/graph/badge.svg)](https://codecov.io/gh/fulcrumgenomics/fg-idprimer)
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/fc4f5fe8dbe34bf784114435b202fab4)](https://www.codacy.com/app/contact_32/fg-idprimer?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=fulcrumgenomics/fg-idprimer&amp;utm_campaign=Badge_Grade)
[![Javadocs](http://javadoc.io/badge/com.fulcrumgenomics/fg-idprimer_2.12.svg)](http://javadoc.io/doc/com.fulcrumgenomics/fg-idprimer_2.12)
[![License](http://img.shields.io/badge/license-MIT-blue.svg)](https://github.com/fulcrumgenomics/fg-idprimer/blob/master/LICENSE)
[![Language](http://img.shields.io/badge/language-scala-brightgreen.svg)](http://www.scala-lang.org/)

fg-idprimer
====

A set of tools to identify primers in on Next Generation Sequencing.  

The main tool to identify primers is `IdentifyPrimers`.
See the [Example Analysis](#example-analysis) section below.

<!---toc start-->
  * [Building](#building)
  * [Command line](#command-line)
  * [Contributing](#contributing)
  * [Authors](#authors)
  * [License](#license)
  * [Example Analysis](#example-analysis)

<!---toc end-->

## Building 
### Cloning the Repository

`git clone https://github.com/fulcrumgenomics/fg-idprimer.git`

### Running the build
fg-idprimer is built using [sbt](http://www.scala-sbt.org/).

Use ```sbt assembly``` to build an executable jar in ```target/scala-2.12/```.

Tests may be run with ```sbt test```. `R` and `ggplot2` are test dependencies.

Java SE 8 is required.

Python 3.6.4+ is required to run the various python scripts.

## Command line

`java -jar target/scala-2.12/fg-idprimer-<version>.jar` to see the commands supported.  Use `java -jar target/scala-2.12/fg-idprimer-<version>.jar <command>` to see the help message for a particular command.

## Contributing

Contributions are welcome and encouraged.
We will do our best to provide an initial response to any pull request or issue within one-week.
For urgent matters, please contact us directly.

## Authors

* [Tim Fennell](https://github.com/tfenne) (maintainer)
* [Nils Homer](https://github.com/nh13) (maintainer)

## License

`fg-idprimer` is open source software released under the [MIT License](https://github.com/fulcrumgenomics/fg-idprimer/blob/master/LICENSE).


## Primer File

The input primers file must be tab-separated, and contain a header, and then one row per primer:

* `pair_id` - the unique identifier for the primer pair
* `primer_id` - the unique identifier for the primer in a primer pair
* `sequence` - the DNA sequence of the primer as ordered, including degenerate bases.
* `strand` -  if the primer maps to the positive strand of the genome or R/r/- if the primer maps to the negative strand
* `ref_name` - the chromosome/contig
* `start` - the start position (1-based)
* `end` - the end position (1-based, inclusive)

`pair_id` and `primer_id` may not have commas.

The primer file must contain headers.  If mapping information is not available for a primer then
ref_name/start/end should be left empty.  E.g.:

```
pair_id    primer_id sequence strand chr  start   end
mapped1    m1_fwd    GATTACA  +      chr1 1010873 1010894
mapped1    m1_rev    ACATTAG  -      chr1 1011118 1011137
unmapped1  u1_fwd    ACGTGTA  +
unmapped1  u1_rev    GGATACA  -
```

## Example Analysis

### Preparing the List of Primers

The (`scripts`)[https://github.com/fulcrumgenomics/fg-idprimer/tree/master/scripts] folder contains two useful scripts for reformatting various primer format files to 
the format expected by `IdentifyPrimers`.
Requirements for these scripts can be found in the same folder.


* detailed usage: `python scripts/reformat_mapped_primers.py --help`
* recommended command line:

```
python scripts/reformat_mapped_primers \
    -i <in.primers.from.idt.txt> \
    -o <out.primers.for.fg-idprimer.txt>
```

For primers without mapping information, instead use `scripts/reformat_unmapepd_primers.py`

### Identify Primers

To identify primers from reads in a BAM file, use the `fg-idprimer IdentifyPrimers` tool.

* detailed usage `fg-idprimer IdentifyPrimers --help`
* recommended command line:

```
java -Xmx8G -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 \
    -jar <path/to/fg-idprimer.jar> \
    IdentifyPrimers \
    -i <in.bam> \
    -p <in.primers.for.fg-idprimer.txt> \
    -o <out/basename>.bam \
    -m <out/basename>.metrics \
    -t <number-of-threads> \
    -k 6 \
    -K 10 \
    --ignore-primer-strand false
```

### Categorizing Primers

How primer matches are classified is dependent on the question being posed and the experiment being performed.
The (`scripts`)[https://github.com/fulcrumgenomics/fg-idprimer/tree/master/scripts] folder contains an example 
script to post-process the BAM file produced by the `fg-idprimer IdentifyPrimers` tool.
Requirements for these scripts can be found in the same folder.

* detailed usage `python scripts/post_process_bam.py --help`
* recommended command line:

```
python scripts/post_process_bam.py \
    -i <out/basename>.bam \
    -o <out/basename>.post_processed.bam \
    --metrics <out/basename>.post_processed.txt \
	--min-insert-length 50 \
	--max-insert-length 250
```

This script outputs two summary files:

1) `<out/basename>.matches.txt` give the count of primer pair match types found by the script.

The following pertains to paired end reads, where primers were matched on the 5' end of the reads.
For a given paired end read, the primer pair match types are:

| Name | Description |
| --- | --- |
| `Canonical` | the primers are from the same pair (i.e. same `pair_id`) but on opposite strands |
| `SelfDimer` | the primers are the same primer! |
| `CrossDimer` | from different pairs (i.e. different `pair_id`s) but outside the expected insert range |
| `NonCanonical` | the primers are (a) from the same pair (i.e. same `pair_id`) and on the same strands, or (b) from different pairs (i.e. different `pair_id`s) but within the expected insert range |
| `Single` | a single primer match on one end of the paired end read |
| `NoMatch` | no primer matches |

2) `<out/basename>.dimers.txt` gives the count of dimers found by the script.

| Name | Primer Pair Match Types |
| --- | --- |
| `Dimer` | `SelfDimer` and `CrossDimer` |
| `NotADimer`| `Canonical` |
| `MaybeDimer` | `NonCanonical` and `Single` |
| `NoInformation` | `NoMatch` |

