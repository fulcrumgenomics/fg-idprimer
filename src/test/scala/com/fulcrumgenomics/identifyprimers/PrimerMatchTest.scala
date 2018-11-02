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

import com.fulcrumgenomics.testing.{SamBuilder, UnitSpec}
import org.scalatest.OptionValues

class PrimerMatchTest extends UnitSpec with OptionValues {
  import PrimerMatch._

  private val primer = new Primer("pair1", "primer1+", "GATTACA", true, "chr1", 1, 7)

  val primerMatches = Seq(
    (LocationBasedPrimerMatch(primer, 2), LocationName),
    (UngappedAlignmentPrimerMatch(primer, 2, Int.MaxValue), UngappedName),
    (UngappedAlignmentPrimerMatch(primer, 2, 3), UngappedName),
    (GappedAlignmentPrimerMatch(primer, 3, 2, 0, 0), GappedName)
  )

  "PrimerMatch.toName" should "give the name of the primer match" in {
    // provide the class
    PrimerMatch.toName[LocationBasedPrimerMatch] shouldBe LocationName
    PrimerMatch.toName[UngappedAlignmentPrimerMatch] shouldBe UngappedName
    PrimerMatch.toName[GappedAlignmentPrimerMatch] shouldBe GappedName

    // with an instance
    primerMatches.foreach { case (primerMatch, expectedName) =>
      PrimerMatch.toName(Some(primerMatch)) shouldBe expectedName
    }
    PrimerMatch.toName(None) shouldBe NoPrimerMatchName
  }

  "PrimerMatch.info" should "give information about the match" in {
    val rec = new SamBuilder().addFrag().value

    // first five fields should all be the same for all info
    primerMatches.map(_._1).map { pm =>
      pm.info(rec).split(InfoDelimiter).toList.take(5).mkString(InfoDelimiter)
    }.distinct should have length 1

    primerMatches.foreach { case (pm, expectedName) =>
      val fields = pm.info(rec).split(InfoDelimiter).toIndexedSeq
      // sixth field should be the match name
      fields(5) shouldBe expectedName

      // the rest of the fields should depend on the type
      val rest = fields.drop(6)
      pm match {
        case m: LocationBasedPrimerMatch     => rest should contain theSameElementsInOrderAs Seq(m.numMismatches).map(_.toString)
        case m: UngappedAlignmentPrimerMatch => rest should contain theSameElementsInOrderAs Seq(m.numMismatches, if (m.nextNumMismatches == Int.MaxValue) "na" else m.nextNumMismatches).map(_.toString)
        case m: GappedAlignmentPrimerMatch   => rest should contain theSameElementsInOrderAs Seq(f"${m.score}%.2f", f"${m.secondBestScore}%.2f", m.start, m.end).map(_.toString)
      }
    }
  }
}
