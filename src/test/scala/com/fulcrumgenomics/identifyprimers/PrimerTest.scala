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

import com.fulcrumgenomics.commons.io.Io
import com.fulcrumgenomics.testing.UnitSpec

class PrimerTest extends UnitSpec {

  private val pairOneForward1 = new Primer("pair1", "primer1", "GATTACA", true, "chr1", 1, 7)
  private val pairOneReverse1 = new Primer("pair1", "primer2", "GATTACA", false, "chr1", 1, 57)
  private val pairOneForward2 = new Primer("pair1", "primer3", "GATTACA", true, "chr1", 1, 67)
  private val pairTwoForward1 = new Primer("pair2", "primer4", "GATTACA", true, "chr1", 1, 7)
  private val pairTwoReverse1 = new Primer("pair2", "primer5", "GATTACA", false, "chr1", 1, 7)

  "Primer.validatePrimers" should "with multiPrimerPairs==true should validate that there is at least one forward and one reverse primer for each pair" in {
    an[Exception] should be thrownBy Primer.validatePrimers(Seq(pairOneForward1), multiPrimerPairs=true) // NOK
    an[Exception] should be thrownBy Primer.validatePrimers(Seq(pairOneReverse1), multiPrimerPairs=true) // NOK
    an[Exception] should be thrownBy Primer.validatePrimers(Seq(pairOneForward1, pairOneForward2), multiPrimerPairs=true) // NOK
    Primer.validatePrimers(Seq(pairOneForward1, pairOneForward2, pairOneReverse1), multiPrimerPairs=true) // OK

  }

  it should "with multiPrimerPairs==false validate that two primers exist for each pair, and that one is forward and the other is reverse" in {
    an[Exception] should be thrownBy Primer.validatePrimers(Seq(pairOneForward1), multiPrimerPairs=false) // NOK
    an[Exception] should be thrownBy Primer.validatePrimers(Seq(pairOneReverse1), multiPrimerPairs=false) // NOK
    an[Exception] should be thrownBy Primer.validatePrimers(Seq(pairOneForward1, pairOneForward2, pairOneReverse1), multiPrimerPairs=false) // NOK
    Primer.validatePrimers(Seq(pairOneForward1, pairOneReverse1), multiPrimerPairs=true) // OK
    Primer.validatePrimers(Seq(pairOneForward2, pairOneReverse1), multiPrimerPairs=true) // OK
  }

  it should "validate we do not have the same refName/start/end for primers" in {
    an[Exception] should be thrownBy Primer.validatePrimers(Seq(pairOneForward1, pairOneReverse1, pairTwoReverse1), multiPrimerPairs=false) // NOK
    Primer.validatePrimers(Seq(pairOneForward1, pairOneReverse1), multiPrimerPairs=false) // OK
  }

  it should "pass validation when primers do not have mapping information" in {
    val primers = Seq(pairOneForward1, pairOneReverse1, pairTwoForward1, pairTwoReverse1).map { primer =>
      primer.copy(ref_name = "")
    }
    Primer.validatePrimers(primers, multiPrimerPairs = false)
  }


  "Primer.strandToForward" should "return true if the positive strand, false if the negative strand" in {
    Seq("+", "f", "fwd", "for", "forward", "positive", "true").foreach { strand => Primer.isPositiveStrand(strand) shouldBe true }
    Seq("-", "r", "rev", "reverse", "negative", "false").foreach { strand => Primer.isPositiveStrand(strand) shouldBe false }
  }

  "Primer" should "round trip primers" in {
    val primers = Seq(pairOneForward1, pairOneReverse1, pairOneForward2)
    val out     = makeTempFile("primers.", ".tab")
    Primer.write(out, primers)
    Primer.read(out, multiPrimerPairs = true) should contain theSameElementsInOrderAs primers
  }

  it should "read a primer file without mapping information" in {
    val out = makeTempFile("primers.", ".tab")
    val lines = Seq(
      "pair_id\tprimer_id\tsequence\tstrand",
      "1\t1+\tGATTACA\t+",
      "1\t1-\tGATTACA\t-"
    )
    Io.writeLines(out, lines)
    val primers = Primer.read(out)
    primers.length shouldBe 2
    primers(0) shouldBe Primer("1", "1+", "GATTACA", positive_strand=true)
    primers(1) shouldBe Primer("1", "1-", "GATTACA", positive_strand=false)
  }

  "Primer.mapped" should "return if the primer is mapped to the ref" in {
    val primer = new Primer("pair1", "primer1", "GATTACA", true)
    primer.mapped shouldBe false
    primer.copy(ref_name="chr1", start=1, end=2).mapped shouldBe true
  }

  "Primer" should "not contain invalid bases" in {
    new Primer("pair1", "primer1", "GARTTACA", true) // OK
    an[Exception] should be thrownBy new Primer("pair1", "primer1", "GA_TTACA", true) // OK
  }
}