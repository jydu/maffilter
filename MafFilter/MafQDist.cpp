//
// File: MafQDist.cpp
// Created by: Julien Dutheil
// Created on: Jul 21 2010
//

/*
Copyright or © or Copr. CNRS

This software is a computer program whose purpose is to test the
homogeneity of the substitution process of a given alignment.

This software is governed by the CeCILL  license under French law and
abiding by the rules of distribution of free software.  You can  use, 
modify and/ or redistribute the software under the terms of the CeCILL
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info". 

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability. 

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or 
data to be ensured and,  more generally, to use and operate it in the 
same conditions as regards security. 

The fact that you are presently reading this means that you have had
knowledge of the CeCILL license and that you accept its terms.
*/

// From the STL:
#include <iostream>
#include <iomanip>
#include <string>
#include <memory>
using namespace std;

//From boost:
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/filter/zlib.hpp>
using namespace boost::iostreams;

// From bpp-core:
#include <Bpp/App/BppApplication.h>

// From bpp-seq:
#include <Bpp/Seq/Io/Maf.all>
#include <Bpp/Seq/SequenceWithQuality.h>

using namespace bpp;

void help()
{

}

int main(int args, char** argv)
{
  cout << "******************************************************************" << endl;
  cout << "*           MAF Quality Distribution, version 0.1.0              *" << endl;
  cout << "* Author: J. Dutheil                        Created on  27/04/10 *" << endl;
  cout << "*                                           Last Modif. 27/04/10 *" << endl;
  cout << "******************************************************************" << endl;
  cout << endl;

  if (args == 1)
  {
    help();
    exit(0);
  }
  
  try
  {
    BppApplication mafqdist(args, argv, "MafQDist");
    mafqdist.startTimer();

    string inputFile = ApplicationTools::getAFilePath("input.file", mafqdist.getParams(), true, true);
    string compress = ApplicationTools::getStringParameter("compression", mafqdist.getParams(), "none");

    filtering_istream stream;
    if (compress == "none") {
    } else if (compress == "gzip") {
      stream.push(gzip_decompressor());
    } else if (compress == "zip") {
      stream.push(zlib_decompressor());
    } else if (compress == "bzip2") {
      stream.push(bzip2_decompressor());
    } else
      throw Exception("Bad compression format: " + compress);
    stream.push(file_source(inputFile));
    
    string species = ApplicationTools::getStringParameter("species", mafqdist.getParams(), "none");
    if (species == "none") {
      throw Exception("A candidate species must be specified.");
    }

    string outputFile = ApplicationTools::getAFilePath("output.file", mafqdist.getParams(), true, false);

    ofstream out(outputFile.c_str(), ios::out);
    MafParser parser(&stream);

    out << "Block\tc0\tc1\tc2\tc3\tc4\tc5\tc6\tc7\tc8\tc9\tcF\tGaps\tUnknown" << endl;
    unsigned int blockCount = 0;
    while (MafBlock* block = parser.nextBlock())
    {
      blockCount++;
      cout << "."; cout.flush();
      for (unsigned int i = 0; i < block->getAlignment().getNumberOfSequences(); ++i) {
        string name = block->getAlignment().getSequence(i).getName();
        size_t pos = name.find(".");
        if (pos != string::npos)
          name = name.substr(0, pos);
        if (name == species) {
          const MafSequence* mafSeq      = &dynamic_cast<const MafSequence&>(block->getAlignment().getSequence(i));
          const SequenceQuality* qscores = &dynamic_cast<const SequenceQuality&>(mafSeq->getAnnotation(SequenceQuality::QUALITY_SCORE));
          map<int, size_t> counts = VectorTools::countValues(qscores->getScores());
          out << blockCount;
          for (int j = 0; j <= 10; ++j)
            out << "\t" << counts[j];
          out << "\t" << counts[-1];
          out << "\t" << counts[-2] << endl;
          break;
        }
      }
      delete block;
    }
    cout << endl;

    out.close();

    mafqdist.done();
  }
  catch (exception& e)
  {
    cout << e.what() << endl;
    exit(-1);
  }

  return (0);
}

