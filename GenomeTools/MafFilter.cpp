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

// From bpp-core:
#include <Bpp/App/BppApplication.h>
#include <Bpp/Text/KeyvalTools.h>
#include <Bpp/Text/StringTokenizer.h>

// From bpp-seq:
#include <Bpp/Seq/Io/Maf.all>
#include <Bpp/Seq/SequenceWithQuality.h>
#include <Bpp/Seq/Feature/Gff/GffFeatureReader.h>
#include <Bpp/Seq/Io/Clustal.h>

using namespace bpp;

void help()
{

}

//void getList(const string& desc, vector<string>& list) throw (Exception)
//{
//  if (desc == "none")
//    throw Exception("You must specify at least one element.");
//  if (desc[0] == '(') {
//    StringTokenizer st(desc.substr(1, desc.size() - 2), ",");
//    while (st.hasMoreToken())
//      list.push_back(st.nextToken());
//  } else {
//    list.push_back(desc);
//  }
//}

int main(int args, char** argv)
{
  cout << "******************************************************************" << endl;
  cout << "*                  MAF Filter, version 0.1.0                     *" << endl;
  cout << "* Author: J. Dutheil                        Created on  10/09/10 *" << endl;
  cout << "*                                           Last Modif. 27/06/12 *" << endl;
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
    vector<filtering_ostream*> ostreams;
    for (unsigned int a = 0; a < actions.size(); a++) {
      string cmdName;
      map<string, string> cmdArgs;
      KeyvalTools::parseProcedure(actions[a], cmdName, cmdArgs);
      ApplicationTools::displayResult("Adding filter", cmdName);
      

      // +-----------------+
      // | Sequence subset |
      // +-----------------+
      if (cmdName == "Subset") {
        //string speciesList = ApplicationTools::getStringParameter("species", cmdArgs, "none");
        bool strict = ApplicationTools::getBooleanParameter("strict", cmdArgs, false);
        ApplicationTools::displayBooleanResult("All species should be in output blocks", strict);
        bool rmdupl = ApplicationTools::getBooleanParameter("rm.duplicates", cmdArgs, false);
        ApplicationTools::displayBooleanResult("Species should be present only once", rmdupl);
        vector<string> species = ApplicationTools::getVectorParameter<string>("species", cmdArgs, ',', "");
        if (species.size() == 0)
          throw Exception("At least one species should be provided for command 'Subset'.");
        //getList(speciesList, species);
        SequenceFilterMafIterator* iterator = new SequenceFilterMafIterator(currentIterator, species, strict, rmdupl);
        iterator->setLogStream(&log);
        currentIterator = iterator;
        its.push_back(iterator);
      }


      // +---------------+
      // | Block merging |
      // +---------------+
      if (cmdName == "Merge") {
        //string speciesList = ApplicationTools::getStringParameter("species", cmdArgs, "none");
        //vector<string> species;
        //getList(speciesList, species);
        vector<string> species = ApplicationTools::getVectorParameter<string>("species", cmdArgs, ',', "");
        if (species.size() == 0)
          throw Exception("At least one species should be provided for command 'Merge'.");

        unsigned int distMax = ApplicationTools::getParameter<unsigned int>("dist.max", cmdArgs, 0);
        ApplicationTools::displayResult("Maximum distance allowed", distMax);
        BlockMergerMafIterator* iterator = new BlockMergerMafIterator(currentIterator, species, distMax);
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
        bool verbose = ApplicationTools::getBooleanParameter("verbose", cmdArgs, true);
        //string speciesList = ApplicationTools::getStringParameter("species", cmdArgs, "none");
        //vector<string> species;
        //getList(speciesList, species);
        vector<string> species = ApplicationTools::getVectorParameter<string>("species", cmdArgs, ',', "");
        if (species.size() == 0)
          throw Exception("At least one species should be provided for command 'XFullGap'.");
        FullGapFilterMafIterator* iterator = new FullGapFilterMafIterator(currentIterator, species);
        iterator->setLogStream(&log);
        iterator->verbose(verbose);
        currentIterator = iterator;
        its.push_back(iterator);
      }


      // +---------------------+
      // | Alignment filtering |
      // +---------------------+
      if (cmdName == "AlnFilter") {
        bool verbose = ApplicationTools::getBooleanParameter("verbose", cmdArgs, true);
        //string speciesList = ApplicationTools::getStringParameter("species", cmdArgs, "none");
        //vector<string> species;
        //getList(speciesList, species);
        vector<string> species = ApplicationTools::getVectorParameter<string>("species", cmdArgs, ',', "");
        if (species.size() == 0)
          throw Exception("At least one species should be provided for command 'AlnFilter'.");
        unsigned int ws = ApplicationTools::getParameter<unsigned int>("window.size", cmdArgs, 10);
        unsigned int st = ApplicationTools::getParameter<unsigned int>("window.step", cmdArgs, 5);
        unsigned int gm = ApplicationTools::getParameter<unsigned int>("max.gap", cmdArgs, 0);
        string outputFile = ApplicationTools::getAFilePath("file", cmdArgs, false, false);
        bool trash = outputFile == "none";
        ApplicationTools::displayResult("Window size", ws);
        ApplicationTools::displayResult("Window step", st);
        ApplicationTools::displayResult("Max. gaps allowed in Window", gm);
        ApplicationTools::displayBooleanResult("Output removed blocks", !trash);
        AlignmentFilterMafIterator* iterator = new AlignmentFilterMafIterator(currentIterator, species, ws, st, gm, !trash);
        iterator->setLogStream(&log);
        iterator->verbose(verbose);
        its.push_back(iterator);

        if (!trash) {
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
          ostreams.push_back(out);
          ApplicationTools::displayResult("File compression for removed blocks", compress);

          //Now build an adaptor for retrieving the trashed blocks:
          TrashIteratorAdapter* trashIt = new TrashIteratorAdapter(iterator);
          //Add an output iterator:
          OutputMafIterator* outIt = new OutputMafIterator(trashIt, out);
          //And then synchronize the two iterators:
          MafIteratorSynchronizer* syncIt = new MafIteratorSynchronizer(iterator, outIt);
          //Returns last iterator:
          currentIterator = syncIt;
          //Keep track of all those iterators:
          its.push_back(trashIt);
          its.push_back(syncIt);
        } else {
          //We only get the remaining blocks here:
          currentIterator = iterator;
        }
      }


      // +-----------------------+
      // | Alignment filtering 2 |
      // +-----------------------+
      if (cmdName == "AlnFilter2") {
        bool verbose = ApplicationTools::getBooleanParameter("verbose", cmdArgs, true);
        //string speciesList = ApplicationTools::getStringParameter("species", cmdArgs, "none");
        //vector<string> species;
        //getList(speciesList, species);
        vector<string> species = ApplicationTools::getVectorParameter<string>("species", cmdArgs, ',', "");
        if (species.size() == 0)
          throw Exception("At least one species should be provided for command 'AlnFilter2'.");
        unsigned int ws = ApplicationTools::getParameter<unsigned int>("window.size", cmdArgs, 10);
        unsigned int st = ApplicationTools::getParameter<unsigned int>("window.step", cmdArgs, 5);
        unsigned int gm = ApplicationTools::getParameter<unsigned int>("max.gap", cmdArgs, 0);
        unsigned int pm = ApplicationTools::getParameter<unsigned int>("max.pos", cmdArgs, 0);
        string outputFile = ApplicationTools::getAFilePath("file", cmdArgs, false, false);
        bool trash = outputFile == "none";
        ApplicationTools::displayResult("Window size", ws);
        ApplicationTools::displayResult("Window step", st);
        ApplicationTools::displayResult("Max. gaps allowed per position", gm);
        ApplicationTools::displayResult("Max. gap positions allowed", pm);
        ApplicationTools::displayBooleanResult("Output removed blocks", !trash);
        AlignmentFilter2MafIterator* iterator = new AlignmentFilter2MafIterator(currentIterator, species, ws, st, gm, pm, !trash);
        iterator->setLogStream(&log);
        iterator->verbose(verbose);
        its.push_back(iterator);

        if (!trash) {
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
          ostreams.push_back(out);
          ApplicationTools::displayResult("File compression for removed blocks", compress);

          //Now build an adaptor for retrieving the trashed blocks:
          TrashIteratorAdapter* trashIt = new TrashIteratorAdapter(iterator);
          //Add an output iterator:
          OutputMafIterator* outIt = new OutputMafIterator(trashIt, out);
          //And then synchronize the two iterators:
          MafIteratorSynchronizer* syncIt = new MafIteratorSynchronizer(iterator, outIt);
          //Returns last iterator:
          currentIterator = syncIt;
          //Keep track of all those iterators:
          its.push_back(trashIt);
          its.push_back(syncIt);
        } else {
          //We only get the remaining blocks here:
          currentIterator = iterator;
        }
      }



      // +----------------+
      // | Mask filtering |
      // +----------------+
      if (cmdName == "MaskFilter") {
        bool verbose = ApplicationTools::getBooleanParameter("verbose", cmdArgs, true);
        //string speciesList = ApplicationTools::getStringParameter("species", cmdArgs, "none");
        //vector<string> species;
        //getList(speciesList, species);
        vector<string> species = ApplicationTools::getVectorParameter<string>("species", cmdArgs, ',', "");
        if (species.size() == 0)
          throw Exception("At least one species should be provided for command 'MaskFilter'.");
        unsigned int ws = ApplicationTools::getParameter<unsigned int>("window.size", cmdArgs, 10);
        unsigned int st = ApplicationTools::getParameter<unsigned int>("window.step", cmdArgs, 5);
        unsigned int mm = ApplicationTools::getParameter<unsigned int>("max.masked", cmdArgs, 0);
        string outputFile = ApplicationTools::getAFilePath("file", cmdArgs, false, false);
        bool trash = outputFile == "none";
        ApplicationTools::displayResult("Window size", ws);
        ApplicationTools::displayResult("Window step", st);
        ApplicationTools::displayResult("Max. masked sites allowed in Window", mm);
        ApplicationTools::displayBooleanResult("Output removed blocks", !trash);
        MaskFilterMafIterator* iterator = new MaskFilterMafIterator(currentIterator, species, ws, st, mm, !trash);
        iterator->setLogStream(&log);
        iterator->verbose(verbose);
        its.push_back(iterator);

        if (!trash) {
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
          ostreams.push_back(out);
          ApplicationTools::displayResult("File compression for removed blocks", compress);

          //Now build an adaptor for retrieving the trashed blocks:
          TrashIteratorAdapter* trashIt = new TrashIteratorAdapter(iterator);
          //Add an output iterator:
          OutputMafIterator* outIt = new OutputMafIterator(trashIt, out);
          //And then synchronize the two iterators:
          MafIteratorSynchronizer* syncIt = new MafIteratorSynchronizer(iterator, outIt);
          //Returns last iterator:
          currentIterator = syncIt;
          //Keep track of all those iterators:
          its.push_back(trashIt);
          its.push_back(syncIt);
        } else {
          //We only get the remaining blocks here:
          currentIterator = iterator;
        }
      }


      // +-------------------+
      // | Quality filtering |
      // +-------------------+
      if (cmdName == "QualFilter") {
        bool verbose = ApplicationTools::getBooleanParameter("verbose", cmdArgs, true);
        //string speciesList = ApplicationTools::getStringParameter("species", cmdArgs, "none");
        //vector<string> species;
        //getList(speciesList, species);
        vector<string> species = ApplicationTools::getVectorParameter<string>("species", cmdArgs, ',', "");
        if (species.size() == 0)
          throw Exception("At least one species should be provided for command 'QualFilter'.");
        unsigned int ws = ApplicationTools::getParameter<unsigned int>("window.size", cmdArgs, 10);
        unsigned int st = ApplicationTools::getParameter<unsigned int>("window.step", cmdArgs, 5);
        double       mq = ApplicationTools::getDoubleParameter("min.qual", cmdArgs, 0);
        string outputFile = ApplicationTools::getAFilePath("file", cmdArgs, false, false);
        bool trash = outputFile == "none";
        ApplicationTools::displayResult("Window size", ws);
        ApplicationTools::displayResult("Window step", st);
        ApplicationTools::displayResult("Min. average quality allowed in Window", mq);
        ApplicationTools::displayBooleanResult("Output removed blocks", !trash);
        QualityFilterMafIterator* iterator = new QualityFilterMafIterator(currentIterator, species, ws, st, mq, !trash);
        iterator->setLogStream(&log);
        iterator->verbose(verbose);
        its.push_back(iterator);

        if (!trash) {
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
          ostreams.push_back(out);
          ApplicationTools::displayResult("File compression for removed blocks", compress);

          //Now build an adaptor for retrieving the trashed blocks:
          TrashIteratorAdapter* trashIt = new TrashIteratorAdapter(iterator);
          //Add an output iterator:
          OutputMafIterator* outIt = new OutputMafIterator(trashIt, out);
          //And then synchronize the two iterators:
          MafIteratorSynchronizer* syncIt = new MafIteratorSynchronizer(iterator, outIt);
          //Returns last iterator:
          currentIterator = syncIt;
          //Keep track of all those iterators:
          its.push_back(trashIt);
          its.push_back(syncIt);
        } else {
          //We only get the remaining blocks here:
          currentIterator = iterator;
        }
      }


      
      // +-------------------------+
      // | Feature-based filtering |
      // +-------------------------+
      if (cmdName == "FeatureFilter") {
        bool verbose = ApplicationTools::getBooleanParameter("verbose", cmdArgs, true);
        string refSpecies = ApplicationTools::getStringParameter("ref_species", cmdArgs, "none");
        string featureFile = ApplicationTools::getAFilePath("feature.file", cmdArgs, false, false);
        string featureFormat = ApplicationTools::getStringParameter("feature.format", cmdArgs, "GFF");
        string outputFile = ApplicationTools::getAFilePath("file", cmdArgs, false, false);
        bool trash = outputFile == "none";
        ApplicationTools::displayResult("Features to remove", featureFile + " (" + featureFormat + ")");
        ApplicationTools::displayResult("features are for species", refSpecies);
        ApplicationTools::displayBooleanResult("Output removed blocks", !trash);
        if (featureFormat != "GFF")
          throw Exception("Sorry, but so far, only GFF features are supported :(");
        ifstream file(featureFile.c_str(), ios::in);
        SequenceFeatureSet featuresSet;
        GffFeatureReader reader(file);
        reader.getAllFeatures(featuresSet);
        FeatureFilterMafIterator* iterator = new FeatureFilterMafIterator(currentIterator, refSpecies, featuresSet, !trash);
        iterator->setLogStream(&log);
        iterator->verbose(verbose);
        its.push_back(iterator);

        if (!trash) {
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
          ostreams.push_back(out);
          ApplicationTools::displayResult("File compression for removed blocks", compress);

          //Now build an adaptor for retrieving the trashed blocks:
          TrashIteratorAdapter* trashIt = new TrashIteratorAdapter(iterator);
          //Add an output iterator:
          OutputMafIterator* outIt = new OutputMafIterator(trashIt, out);
          //And then synchronize the two iterators:
          MafIteratorSynchronizer* syncIt = new MafIteratorSynchronizer(iterator, outIt);
          //Returns last iterator:
          currentIterator = syncIt;
          //Keep track of all those iterators:
          its.push_back(trashIt);
          its.push_back(syncIt);
        } else {
          //We only get the remaining blocks here:
          currentIterator = iterator;
        }
      }



      // +----------------------+
      // | Block size filtering |
      // +----------------------+
      if (cmdName == "MinBlockSize") {
        unsigned int minSize = ApplicationTools::getParameter<unsigned int>("min.size", cmdArgs, 0);
        ApplicationTools::displayResult("Minimum block size required", minSize);
        BlockSizeMafIterator* iterator = new BlockSizeMafIterator(currentIterator, minSize);
        iterator->setLogStream(&log);
        currentIterator = iterator;
        its.push_back(iterator);
      }


      // +----------------------+
      // | Chromosome filtering |
      // +----------------------+
      if (cmdName == "SelectChr") {
        string ref = ApplicationTools::getStringParameter("reference", cmdArgs, "");
        ApplicationTools::displayResult("Reference species", ref);
        string chr = ApplicationTools::getStringParameter("chromosome", cmdArgs, "");
        ApplicationTools::displayResult("Chromosome", chr);
        ChromosomeMafIterator* iterator = new ChromosomeMafIterator(currentIterator, ref, chr);
        iterator->setLogStream(&log);
        currentIterator = iterator;
        its.push_back(iterator);
      }


      // +---------------------+
      // | Duplicate filtering |
      // +---------------------+
      if (cmdName == "DuplicateFilter") {
        string ref = ApplicationTools::getStringParameter("reference", cmdArgs, "");
        ApplicationTools::displayResult("Reference species", ref);
        DuplicateFilterMafIterator* iterator = new DuplicateFilterMafIterator(currentIterator, ref);
        iterator->setLogStream(&log);
        currentIterator = iterator;
        its.push_back(iterator);
      }


      // +---------------------+
      // | Sequence statistics |
      // +---------------------+
      //if (cmdName == "SequenceStatistics") {
      //  string speciesList = ApplicationTools::getStringParameter("species", cmdArgs, "none");
      //  vector<string> species;
      //  getList(speciesList, species);
      //  string outputFile = ApplicationTools::getAFilePath("file", cmdArgs, true, false);
      //  ApplicationTools::displayResult("Output file", outputFile);
      //  auto_ptr<ostream> ofs(new ofstream(outputFile.c_str(), ios::out));
      //  StlOutputStream* output = new StlOutputStream(ofs);
      //  SequenceStatisticsMafIterator* iterator = new SequenceStatisticsMafIterator(currentIterator, species, output);
      //  currentIterator = iterator;
      //  its.push_back(iterator);
      //}
      if (cmdName == "SequenceStatistics") {
        vector<string> statisticsDesc = ApplicationTools::getVectorParameter<string>("statistics", cmdArgs, ',', "", "", false, true);
        //Parse all statistics:
        vector<MafStatistics*> statistics;
        //Get output file:
        string outputFile = ApplicationTools::getAFilePath("file", cmdArgs, true, false);
        ApplicationTools::displayResult("Output file", outputFile);
        auto_ptr<ostream> ofs(new ofstream(outputFile.c_str(), ios::out));
        StlOutputStream* output = new StlOutputStream(ofs);
        SequenceStatisticsMafIterator* iterator = new SequenceStatisticsMafIterator(currentIterator, statistics);
        CsvStatisticsOutputIterationListener* listener = new CsvStatisticsOutputIterationListener(iterator, output);
        iterator->addIterationListener(listener);
        currentIterator = iterator;
        its.push_back(iterator);
      }


      // +------------------------------+
      // | Pairwise sequence statistics |
      // +------------------------------+
      //if (cmdName == "PairwiseSequenceStatistics") {
      //  string species1 = ApplicationTools::getStringParameter("species1", cmdArgs, "none");
      //  string species2 = ApplicationTools::getStringParameter("species2", cmdArgs, "none");
      //  string outputFile = ApplicationTools::getAFilePath("file", cmdArgs, true, false);
      //  ApplicationTools::displayResult("Output file", outputFile);
      //  auto_ptr<ostream> ofs(new ofstream(outputFile.c_str(), ios::out));
      //  StlOutputStream* output = new StlOutputStream(ofs);
      //  PairwiseSequenceStatisticsMafIterator* iterator = new PairwiseSequenceStatisticsMafIterator(currentIterator, species1, species2, output);
      //  currentIterator = iterator;
      //  its.push_back(iterator);
      //}


      
      // +--------------------+
      // | Feature extraction |
      // +--------------------+
      if (cmdName == "ExtractFeature") {
        bool verbose         = ApplicationTools::getBooleanParameter("verbose", cmdArgs, true);
        bool ignoreStrand    = ApplicationTools::getBooleanParameter("ignore_strand", cmdArgs, false);
        string refSpecies    = ApplicationTools::getStringParameter("ref_species", cmdArgs, "none");
        string featureFile   = ApplicationTools::getAFilePath("feature.file", cmdArgs, false, false);
        string featureFormat = ApplicationTools::getStringParameter("feature.format", cmdArgs, "GFF");
        ApplicationTools::displayResult("Features to extract", featureFile + " (" + featureFormat + ")");
        ApplicationTools::displayResult("Features are for species", refSpecies);
        ApplicationTools::displayBooleanResult("Features are strand-aware", !ignoreStrand);
        //string featuresList = ApplicationTools::getStringParameter("feature.select", cmdArgs, "none");
        if (featureFormat != "GFF")
          throw Exception("Sorry, but so far, only GFF features are supported :(");
        ifstream file(featureFile.c_str(), ios::in);
        SequenceFeatureSet featuresSet;
        GffFeatureReader reader(file);
        reader.getAllFeatures(featuresSet);
        
        vector<string> features = ApplicationTools::getVectorParameter<string>("feature.select", cmdArgs, ',', "", "", false, true);
        if (features.size() == 0)
          throw Exception("Error, please set some features to extract!");
        if (features[0] == "all") { 
          set<string> tmp = featuresSet.getTypes();
          features = vector<string>(tmp.begin(), tmp.end());
        }
        if (features.size() == 0)
          throw Exception("Error, no feature to extract!");
        string lst = features[0];
        for (size_t i = 1; i < features.size(); ++i)
          lst += ", " + features[i];
        ApplicationTools::displayResult("Features to extract", lst);
        auto_ptr<SequenceFeatureSet> featuresSet2(featuresSet.getSubsetForType(features));

        FeatureExtractor* iterator = new FeatureExtractor(currentIterator, refSpecies, *featuresSet2, ignoreStrand);
        iterator->setLogStream(&log);
        iterator->verbose(verbose);
        its.push_back(iterator);

        currentIterator = iterator;
      }



      // +--------+
      // | Output |
      // +--------+
      if (cmdName == "Output") {
        string outputFile = ApplicationTools::getAFilePath("file", cmdArgs, true, false);
        compress = ApplicationTools::getStringParameter("compression", cmdArgs, "none");
        ApplicationTools::displayResult("Output file", outputFile);
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
        ostreams.push_back(out);
        ApplicationTools::displayResult("File compression", compress);
        bool mask = ApplicationTools::getBooleanParameter("mask", cmdArgs, true);
        ApplicationTools::displayBooleanResult("Output mask", mask);
        OutputMafIterator* iterator = new OutputMafIterator(currentIterator, out, mask); //NB: there is a memory leak here because the stream is never deleted... TODO
        currentIterator = iterator;
        its.push_back(iterator);
      }



      // +-------------------+
      // | Output alignments |
      // +-------------------+
      if (cmdName == "OutputAlignments") {
        string outputFile = ApplicationTools::getAFilePath("file", cmdArgs, true, false);
        compress = ApplicationTools::getStringParameter("compression", cmdArgs, "none");
        ApplicationTools::displayResult("Output alignment file", outputFile);
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
        ostreams.push_back(out);
        ApplicationTools::displayResult("File compression", compress);
        bool mask = ApplicationTools::getBooleanParameter("mask", cmdArgs, true);
        ApplicationTools::displayBooleanResult("Output mask", mask);
        OutputAlignmentMafIterator* iterator = new OutputAlignmentMafIterator(currentIterator, out, mask); //NB: there is a memory leak here because the stream is never deleted... TODO
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

    //Flush all streams:
    for (size_t i = 0; i < ostreams.size(); ++i) {
      close(*ostreams[i]);
    }

    //Clean memory:
    for (size_t i = 0; i < its.size(); ++i) {
      delete its[i];
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

