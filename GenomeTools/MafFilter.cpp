//
// File: MafFilter.cpp
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

// From Utils:
#include <Utils/BppApplication.h>
#include <Utils/KeyvalTools.h>
#include <Utils/StringTokenizer.h>

// From SeqLib:
#include <Seq/MafAlignmentParser.h>
#include <Seq/SequenceWithQuality.h>

using namespace bpp;

void help()
{

}

void getList(const string& desc, vector<string>& list) throw (Exception)
{
  if (desc == "none")
    throw Exception("You must specify at least one species name.");
  if (desc[0] == '(') {
    StringTokenizer st(desc.substr(1, desc.size() - 2), ",");
    while (st.hasMoreToken())
      list.push_back(st.nextToken());
  } else {
    list.push_back(desc);
  }
}

int main(int args, char** argv)
{
  cout << "******************************************************************" << endl;
  cout << "*                  MAF Filter, version 0.1.0                     *" << endl;
  cout << "* Author: J. Dutheil                        Created on  10/09/10 *" << endl;
  cout << "*                                           Last Modif. 10/09/10 *" << endl;
  cout << "******************************************************************" << endl;
  cout << endl;

  if (args == 1)
  {
    help();
    exit(0);
  }
  
  try
  {
    BppApplication maffilter(args, argv, "MafFilter");
    maffilter.startTimer();

    string inputFile = ApplicationTools::getAFilePath("input.file", maffilter.getParams(), true, true);
    string compress = ApplicationTools::getStringParameter("input.file.compression", maffilter.getParams(), "none");

    filtering_istream stream;
    if (compress == "none") {
    } else if (compress == "gzip") {
      stream.push(gzip_decompressor());
    } else if (compress == "zip") {
      stream.push(zlib_decompressor());
    } else if (compress == "bzip2") {
      stream.push(bzip2_decompressor());
    } else
      throw Exception("Bad input incompression format: " + compress);
    stream.push(file_source(inputFile));
    
    string logFile = ApplicationTools::getAFilePath("output.log", maffilter.getParams(), true, false);

    StlOutputStream log(auto_ptr<ostream>(new ofstream(logFile.c_str(), ios::out)));
    MafAlignmentParser parser(&stream, true);

    vector<string> actions = ApplicationTools::getVectorParameter<string>("maf.filter", maffilter.getParams(), ',', "", "", false, false);
    MafIterator* currentIterator = &parser;
    vector<MafIterator*> its;
    for (unsigned int a = 0; a < actions.size(); a++) {
      string cmdName;
      map<string, string> cmdArgs;
      KeyvalTools::parseProcedure(actions[a], cmdName, cmdArgs);
      ApplicationTools::displayResult("Adding filter", cmdName);
      

      // +-----------------+
      // | Sequence subset |
      // +-----------------+
      if (cmdName == "Subset") {
        string speciesList = ApplicationTools::getStringParameter("species", cmdArgs, "none");
        bool strict = ApplicationTools::getBooleanParameter("strict", cmdArgs, false);
        ApplicationTools::displayBooleanResult("All species should be in output blocks", strict);
        bool rmdupl = ApplicationTools::getBooleanParameter("rm.duplicates", cmdArgs, false);
        ApplicationTools::displayBooleanResult("Species should be present only once", rmdupl);
        vector<string> species;
        getList(speciesList, species);
        SequenceFilterMafIterator* iterator = new SequenceFilterMafIterator(currentIterator, species, strict, rmdupl);
        iterator->setLogStream(&log);
        currentIterator = iterator;
        its.push_back(iterator);
      }


      // +---------------+
      // | Block merging |
      // +---------------+
      if (cmdName == "Merge") {
        string speciesList = ApplicationTools::getStringParameter("species", cmdArgs, "none");
        vector<string> species;
        getList(speciesList, species);
        BlockMergerMafIterator* iterator = new BlockMergerMafIterator(currentIterator, species);
        iterator->setLogStream(&log);
        string ignoreChrList = ApplicationTools::getStringParameter("ignore.chr", cmdArgs, "none");
        if (ignoreChrList != "none") {
          if (ignoreChrList[0] == '(') {
            StringTokenizer st(ignoreChrList.substr(1, ignoreChrList.size() - 2), ",");
            while (st.hasMoreToken())
              iterator->ignoreChromosome(st.nextToken());
          } else {
            iterator->ignoreChromosome(ignoreChrList);
          }
        }
        currentIterator = iterator;
        its.push_back(iterator);
      }


      // +--------------------+
      // | Full gap filtering |
      // +--------------------+
      if (cmdName == "XFullGap") {
        string speciesList = ApplicationTools::getStringParameter("species", cmdArgs, "none");
        vector<string> species;
        getList(speciesList, species);
        FullGapFilterMafIterator* iterator = new FullGapFilterMafIterator(currentIterator, species);
        iterator->setLogStream(&log);
        currentIterator = iterator;
        its.push_back(iterator);
      }


      // +---------------------+
      // | Alignment filtering |
      // +---------------------+
      if (cmdName == "AlnFilter") {
        string speciesList = ApplicationTools::getStringParameter("species", cmdArgs, "none");
        vector<string> species;
        getList(speciesList, species);
        unsigned int ws = ApplicationTools::getParameter<unsigned int>("window.size", cmdArgs, 10);
        unsigned int st = ApplicationTools::getParameter<unsigned int>("window.step", cmdArgs, 5);
        unsigned int gm = ApplicationTools::getParameter<unsigned int>("max.gap", cmdArgs, 0);
        ApplicationTools::displayResult("Window size", ws);
        ApplicationTools::displayResult("Window step", st);
        ApplicationTools::displayResult("Max. gaps allowed in Window", gm);
        AlignmentFilterMafIterator* iterator = new AlignmentFilterMafIterator(currentIterator, species, ws, st, gm);
        iterator->setLogStream(&log);
        currentIterator = iterator;
        its.push_back(iterator);
      }


      // +--------+
      // | Output |
      // +--------+
      if (cmdName == "Output") {
        string outputFile = ApplicationTools::getAFilePath("file", cmdArgs, true, false);
        compress = ApplicationTools::getStringParameter("compression", cmdArgs, "none");
        filtering_ostream* out = new filtering_ostream;
        if (compress == "none") {
        } else if (compress == "gzip") {
          out->push(gzip_compressor());
        } else if (compress == "zip") {
          out->push(zlib_compressor());
        } else if (compress == "bzip2") {
          out->push(bzip2_compressor());
        } else
          throw Exception("Bad output compression format: " + compress);
        out->push(file_sink(outputFile));
        ApplicationTools::displayResult("File compression", compress);
        bool mask = ApplicationTools::getBooleanParameter("mask", cmdArgs, true);
        ApplicationTools::displayBooleanResult("Output mask", mask);
        OutputMafIterator* iterator = new OutputMafIterator(currentIterator, out, mask); //NB: there is a memory leak here because the stream is never deleted... TODO
        currentIterator = iterator;
        its.push_back(iterator);
      }

    }

    //Now loop over the last iterator and that's it!
    while (MafBlock* block = currentIterator->nextBlock())
    {
      cout << "."; cout.flush();
      delete block;
    }

    maffilter.done();
  }
  catch (exception& e)
  {
    cout << e.what() << endl;
    exit(-1);
  }

  return (0);
}

