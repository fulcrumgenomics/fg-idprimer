/*
 * The MIT License
 *
 * Copyright (c) 2018 Fulcrum Genomics
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
 *
 */

package com.fulcrumgenomics.cmdline

/**
  * Main program for fgbio that loads everything up and runs the appropriate sub-command
  */
object FgIdPrimersMain {
  /** The main method */
  def main(args: Array[String]): Unit = new FgIdPrimersMain().makeItSoAndExit(args)

  /**
    * Exception class intended to be used by [[FgIdPrimersMain]] and [[FgBioTool]] to communicate
    * non-exceptional(!) failures when running a tool.
    */
  case class FailureException private[cmdline] (exit:Int = 1, message:Option[String] = None) extends RuntimeException
}

class FgIdPrimersMain extends FgBioMain {

  /** The packages we wish to include in our command line **/
  override protected def packageList: List[String] = List[String]("com.fulcrumgenomics.identifyprimers")
}
