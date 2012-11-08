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

#include "MafIterators.h"

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
#include <Bpp/Seq/SequenceWithQuality.h>
#include <Bpp/Seq/Io/BppOAlignmentWriterFormat.h>
#include <Bpp/Seq/Container/SiteContainerTools.h>

// From bpp-seq-omics:
#include <Bpp/Seq/Io/Maf.all>
#include <Bpp/Seq/Feature/Gff/GffFeatureReader.h>
#include <Bpp/Seq/Feature/Gtf/GtfFeatureReader.h>

// From bpp-phyl:
#include <Bpp/Phyl/Distance.all>

using namespace bpp;

void help()
{

}

int main(int args, char** argv)
{
  cout << "******************************************************************" << endl;
  cout << "*                  MAF Filter, version 1.0.0                     *" << endl;
  cout << "* Author: J. Dutheil                        Created on  10/09/10 *" << endl;
  cout << "*                                           Last Modif. 08/11/12 *" << endl;
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
    string inputFormat = ApplicationTools::getStringParameter("input.format", maffilter.getParams(), "Maf", "", true, false);
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

    MafIterator* currentIterator;
    if (inputFormat == "Maf") {
      currentIterator = new MafAlignmentParser(&stream, true);
    } else if (inputFormat == "Fasta") {
      ISequenceStream* seqStream = new Fasta();
      currentIterator = new SequenceStreamToMafIterator(seqStream, &stream);
    } else throw Exception("Unvalid input format: " + inputFormat);

    vector<string> actions = ApplicationTools::getVectorParameter<string>("maf.filter", maffilter.getParams(), ',', "", "", false, false);
    vector<MafIterator*> its;
    its.push_back(currentIterator);
    vector<filtering_ostream*> ostreams;
    for (unsigned int a = 0; a < actions.size(); a++) {
      string cmdName;
      map<string, string> cmdArgs;
      KeyvalTools::parseProcedure(actions[a], cmdName, cmdArgs);
      (*ApplicationTools::message << "-------------------------------------------------------------------").endLine();
      ApplicationTools::displayResult("Adding filter", cmdName);
      
      bool verbose = ApplicationTools::getBooleanParameter("verbose", cmdArgs, true, "", true, false);
      ApplicationTools::displayBooleanResult("-- Verbose", verbose);

      // +-----------------+
      // | Sequence subset |
      // +-----------------+
      if (cmdName == "Subset") {
        //string speciesList = ApplicationTools::getStringParameter("species", cmdArgs, "none");
        bool strict = ApplicationTools::getBooleanParameter("strict", cmdArgs, false);
        ApplicationTools::displayBooleanResult("-- All species should be in output blocks", strict);
        bool keep = ApplicationTools::getBooleanParameter("keep", cmdArgs, false);
        ApplicationTools::displayBooleanResult("-- Sequences not in the list will be kept", keep);
        if (cmdArgs.find("rm.duplicates") != cmdArgs.end()) {
          throw Exception("rm.duplicates argument in Subset is deprecated: use remove_duplicates instead.");
        }
        bool rmdupl = ApplicationTools::getBooleanParameter("remove_duplicates", cmdArgs, false);
        ApplicationTools::displayBooleanResult("-- Species should be present only once", rmdupl);
        vector<string> species = ApplicationTools::getVectorParameter<string>("species", cmdArgs, ',', "");
        if (species.size() == 0)
          throw Exception("At least one species should be provided for command 'Subset'.");
        //getList(speciesList, species);
        SequenceFilterMafIterator* iterator = new SequenceFilterMafIterator(currentIterator, species, strict, keep, rmdupl);
        iterator->setLogStream(&log);
        iterator->verbose(verbose);
        currentIterator = iterator;
        its.push_back(iterator);
      }


      // +---------------+
      // | Block merging |
      // +---------------+
      if (cmdName == "Merge") {
        vector<string> species = ApplicationTools::getVectorParameter<string>("species", cmdArgs, ',', "");
        if (species.size() == 0)
          throw Exception("At least one species should be provided for command 'Merge'.");

        if (cmdArgs.find("dist.max") != cmdArgs.end()) {
          throw Exception("dist.max argument in Merge is deprecated: use dist_max instead.");
        }
        unsigned int distMax = ApplicationTools::getParameter<unsigned int>("dist_max", cmdArgs, 0);
        ApplicationTools::displayResult("-- Maximum distance allowed", distMax);
        BlockMergerMafIterator* iterator = new BlockMergerMafIterator(currentIterator, species, distMax);
        iterator->setLogStream(&log);
        iterator->verbose(verbose);
        if (cmdArgs.find("ignore.chr") != cmdArgs.end()) {
          throw Exception("ignore.chr argument in Merge is deprecated: use ignore_chr instead.");
        }
        string ignoreChrList = ApplicationTools::getStringParameter("ignore_chr", cmdArgs, "none");
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


      // +---------------------+
      // | Block concatenation |
      // +---------------------+
      if (cmdName == "Concatenate") {
        unsigned int minimumSize = ApplicationTools::getParameter<unsigned int>("minimum_size", cmdArgs, 0);
        ApplicationTools::displayResult("-- Minimum final block size", minimumSize);
        ConcatenateMafIterator* iterator = new ConcatenateMafIterator(currentIterator, minimumSize);
        iterator->setLogStream(&log);
        iterator->verbose(verbose);
        currentIterator = iterator;
        its.push_back(iterator);
      }


      // +--------------------+
      // | Full gap filtering |
      // +--------------------+
      if (cmdName == "XFullGap") {
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
        vector<string> species = ApplicationTools::getVectorParameter<string>("species", cmdArgs, ',', "");
        if (species.size() == 0)
          throw Exception("At least one species should be provided for command 'AlnFilter'.");
        unsigned int ws = ApplicationTools::getParameter<unsigned int>("window.size", cmdArgs, 10);
        unsigned int st = ApplicationTools::getParameter<unsigned int>("window.step", cmdArgs, 5);
        unsigned int gm = ApplicationTools::getParameter<unsigned int>("max.gap", cmdArgs, 0);
        double em       = ApplicationTools::getParameter<double>("max.ent", cmdArgs, 0); //Default means no entropy threshold
        bool missingAsGap = ApplicationTools::getParameter<bool>("missing_as_gap", cmdArgs, false);
        string outputFile = ApplicationTools::getAFilePath("file", cmdArgs, false, false);
        bool trash = outputFile == "none";
        ApplicationTools::displayResult("-- Window size", ws);
        ApplicationTools::displayResult("-- Window step", st);
        ApplicationTools::displayResult("-- Max. gaps allowed in Window", gm);
        ApplicationTools::displayResult("-- Max. total entropy in Window", em);
        ApplicationTools::displayBooleanResult("-- Missing sequence replaced by gaps", missingAsGap);
        ApplicationTools::displayBooleanResult("-- Output removed blocks", !trash);
        AlignmentFilterMafIterator* iterator = new AlignmentFilterMafIterator(currentIterator, species, ws, st, gm, em, !trash, missingAsGap);
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
          ApplicationTools::displayResult("-- File compression for removed blocks", compress);

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
        vector<string> species = ApplicationTools::getVectorParameter<string>("species", cmdArgs, ',', "");
        if (species.size() == 0)
          throw Exception("At least one species should be provided for command 'AlnFilter2'.");
        unsigned int ws = ApplicationTools::getParameter<unsigned int>("window.size", cmdArgs, 10);
        unsigned int st = ApplicationTools::getParameter<unsigned int>("window.step", cmdArgs, 5);
        unsigned int gm = ApplicationTools::getParameter<unsigned int>("max.gap", cmdArgs, 0);
        unsigned int pm = ApplicationTools::getParameter<unsigned int>("max.pos", cmdArgs, 0);
        bool missingAsGap = ApplicationTools::getParameter<bool>("missing_as_gap", cmdArgs, false);
        string outputFile = ApplicationTools::getAFilePath("file", cmdArgs, false, false);
        bool trash = outputFile == "none";
        ApplicationTools::displayResult("-- Window size", ws);
        ApplicationTools::displayResult("-- Window step", st);
        ApplicationTools::displayResult("-- Max. gaps allowed per position", gm);
        ApplicationTools::displayResult("-- Max. gap positions allowed", pm);
        ApplicationTools::displayBooleanResult("-- Missing sequence replaced by gaps", missingAsGap);
        ApplicationTools::displayBooleanResult("-- Output removed blocks", !trash);
        AlignmentFilter2MafIterator* iterator = new AlignmentFilter2MafIterator(currentIterator, species, ws, st, gm, pm, !trash, missingAsGap);
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
          ApplicationTools::displayResult("-- File compression for removed blocks", compress);

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
      // | Entropy filtering |
      // +-------------------+
      if (cmdName == "EntropyFilter") {
        vector<string> species = ApplicationTools::getVectorParameter<string>("species", cmdArgs, ',', "");
        if (species.size() == 0)
          throw Exception("At least one species should be provided for command 'AlnFilter2'.");
        unsigned int ws = ApplicationTools::getParameter<unsigned int>("window.size", cmdArgs, 10);
        unsigned int st = ApplicationTools::getParameter<unsigned int>("window.step", cmdArgs, 5);
        double       em = ApplicationTools::getParameter<double>      ("max.ent", cmdArgs, 0);
        unsigned int pm = ApplicationTools::getParameter<unsigned int>("max.pos", cmdArgs, 0);
        bool missingAsGap = ApplicationTools::getParameter<bool>("missing_as_gap", cmdArgs, false);
        bool ignoreGaps   = ApplicationTools::getParameter<bool>("ignore_gaps", cmdArgs, false);
        string outputFile = ApplicationTools::getAFilePath("file", cmdArgs, false, false);
        bool trash = outputFile == "none";
        ApplicationTools::displayResult("-- Window size", ws);
        ApplicationTools::displayResult("-- Window step", st);
        ApplicationTools::displayResult("-- Max. entropy allowed per position", em);
        ApplicationTools::displayResult("-- Max. high entropy positions allowed", pm);
        ApplicationTools::displayBooleanResult("-- Missing sequence replaced by gaps", missingAsGap);
        ApplicationTools::displayBooleanResult("-- Gaps should be ignored", ignoreGaps);
        ApplicationTools::displayBooleanResult("-- Output removed blocks", !trash);
        if (ignoreGaps && missingAsGap)
          throw Exception("Error, incompatible options ingore_gaps=yes and missing_as_gap=yes.");
        EntropyFilterMafIterator* iterator = new EntropyFilterMafIterator(currentIterator, species, ws, st, em, pm, !trash, missingAsGap, ignoreGaps);
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
          ApplicationTools::displayResult("-- File compression for removed blocks", compress);

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
        vector<string> species = ApplicationTools::getVectorParameter<string>("species", cmdArgs, ',', "");
        if (species.size() == 0)
          throw Exception("At least one species should be provided for command 'MaskFilter'.");
        unsigned int ws = ApplicationTools::getParameter<unsigned int>("window.size", cmdArgs, 10);
        unsigned int st = ApplicationTools::getParameter<unsigned int>("window.step", cmdArgs, 5);
        unsigned int mm = ApplicationTools::getParameter<unsigned int>("max.masked", cmdArgs, 0);
        string outputFile = ApplicationTools::getAFilePath("file", cmdArgs, false, false);
        bool trash = outputFile == "none";
        ApplicationTools::displayResult("-- Window size", ws);
        ApplicationTools::displayResult("-- Window step", st);
        ApplicationTools::displayResult("-- Max. masked sites allowed in Window", mm);
        ApplicationTools::displayBooleanResult("-- Output removed blocks", !trash);
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
          ApplicationTools::displayResult("-- File compression for removed blocks", compress);

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
        vector<string> species = ApplicationTools::getVectorParameter<string>("species", cmdArgs, ',', "");
        if (species.size() == 0)
          throw Exception("At least one species should be provided for command 'QualFilter'.");
        unsigned int ws = ApplicationTools::getParameter<unsigned int>("window.size", cmdArgs, 10);
        unsigned int st = ApplicationTools::getParameter<unsigned int>("window.step", cmdArgs, 5);
        double       mq = ApplicationTools::getDoubleParameter("min.qual", cmdArgs, 0);
        string outputFile = ApplicationTools::getAFilePath("file", cmdArgs, false, false);
        bool trash = outputFile == "none";
        ApplicationTools::displayResult("-- Window size", ws);
        ApplicationTools::displayResult("-- Window step", st);
        ApplicationTools::displayResult("-- Min. average quality allowed in Window", mq);
        ApplicationTools::displayBooleanResult("-- Output removed blocks", !trash);
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
          ApplicationTools::displayResult("-- File compression for removed blocks", compress);

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
        string refSpecies = ApplicationTools::getStringParameter("ref_species", cmdArgs, "none");
        string featureFile = ApplicationTools::getAFilePath("feature.file", cmdArgs, false, false);
        string featureFormat = ApplicationTools::getStringParameter("feature.format", cmdArgs, "GFF");
        string outputFile = ApplicationTools::getAFilePath("file", cmdArgs, false, false);
        bool trash = outputFile == "none";
        ApplicationTools::displayResult("-- Features to remove", featureFile + " (" + featureFormat + ")");
        ApplicationTools::displayResult("-- Features are for species", refSpecies);
        ApplicationTools::displayBooleanResult("-- Output removed blocks", !trash);
        compress = ApplicationTools::getStringParameter("feature.file.compression", cmdArgs, "none");
        filtering_istream featureStream;
        if (compress == "none") {
        } else if (compress == "gzip") {
          featureStream.push(gzip_decompressor());
        } else if (compress == "zip") {
          featureStream.push(zlib_decompressor());
        } else if (compress == "bzip2") {
          featureStream.push(bzip2_decompressor());
        } else
          throw Exception("Bad input incompression format: " + compress);
        featureStream.push(file_source(featureFile));
        SequenceFeatureSet featuresSet;
        if (featureFormat == "GFF") {
          GffFeatureReader reader(featureStream);
          reader.getAllFeatures(featuresSet);
        } else if (featureFormat == "GTF") {
          GtfFeatureReader reader(featureStream);
          reader.getAllFeatures(featuresSet);
        } else
          throw Exception("Unsupported feature format: " + featureFormat);
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
          ApplicationTools::displayResult("-- File compression for removed blocks", compress);

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



      // +------------------------+
      // | Block length filtering |
      // +------------------------+
      if (cmdName == "MinBlockLength") {
        if (cmdArgs.find("min.length") != cmdArgs.end()) {
          throw Exception("min.length argument in MinBlockLength is deprecated: use min_length instead.");
        }
        unsigned int minLength = ApplicationTools::getParameter<unsigned int>("min_length", cmdArgs, 0);
        ApplicationTools::displayResult("-- Minimum block length required", minLength);
        BlockLengthMafIterator* iterator = new BlockLengthMafIterator(currentIterator, minLength);
        iterator->setLogStream(&log);
        currentIterator = iterator;
        its.push_back(iterator);
      }


      // +----------------------+
      // | Block size filtering |
      // +----------------------+
      if (cmdName == "MinBlockSize") {
        if (cmdArgs.find("min.size") != cmdArgs.end()) {
          throw Exception("min.size argument in MinBlockSize is deprecated: use min_size instead.");
        }
        unsigned int minSize = ApplicationTools::getParameter<unsigned int>("min_size", cmdArgs, 0);
        ApplicationTools::displayResult("-- Minimum block size required", minSize);
        if (minSize > 5)
          ApplicationTools::displayWarning("!! Warning, in previous version of maffilter BlockLength was named BlockSize... Check!");
        BlockSizeMafIterator* iterator = new BlockSizeMafIterator(currentIterator, minSize);
        iterator->setLogStream(&log);
        currentIterator = iterator;
        its.push_back(iterator);
      }


      // +----------------------+
      // | Chromosome filtering |
      // +----------------------+
      if (cmdName == "SelectChr") {
        if (cmdArgs.find("reference") != cmdArgs.end()) {
          throw Exception("reference argument in SelectChr is deprecated: use ref_species instead.");
        }
        string ref = ApplicationTools::getStringParameter("ref_species", cmdArgs, "");
        ApplicationTools::displayResult("-- Reference species", ref);
        string chr = ApplicationTools::getStringParameter("chromosome", cmdArgs, "");
        ApplicationTools::displayResult("-- Chromosome", chr);
        ChromosomeMafIterator* iterator = new ChromosomeMafIterator(currentIterator, ref, chr);
        iterator->setLogStream(&log);
        iterator->verbose(verbose);
        currentIterator = iterator;
        its.push_back(iterator);
      }


      // +---------------------+
      // | Duplicate filtering |
      // +---------------------+
      //Nb: this is kind of deprecated, should be done better by looking at partial overlap.
      //could be useful for debugging though. We do not report it in the documentation for now.
      if (cmdName == "DuplicateFilter") {
        string ref = ApplicationTools::getStringParameter("reference", cmdArgs, "");
        ApplicationTools::displayResult("-- Reference species", ref);
        DuplicateFilterMafIterator* iterator = new DuplicateFilterMafIterator(currentIterator, ref);
        iterator->setLogStream(&log);
        iterator->verbose(verbose);
        currentIterator = iterator;
        its.push_back(iterator);
      }


      // +---------------------+
      // | Sequence statistics |
      // +---------------------+
      if (cmdName == "SequenceStatistics") {
        vector<string> statisticsDesc = ApplicationTools::getVectorParameter<string>("statistics", cmdArgs, ',', "", "", false, true);
        
        //Parse all statistics:
        vector<MafStatistics*> statistics;
        for (size_t i = 0; i < statisticsDesc.size(); ++i) {
          string statName;
          map<string, string> statArgs;
          KeyvalTools::parseProcedure(statisticsDesc[i], statName, statArgs);
          MafStatistics* mafStat = 0;
          string statDesc = "";
          if (statName == "BlockSize") {
            mafStat = new BlockSizeMafStatistics();
          } else if (statName == "BlockLength") {
            mafStat = new BlockLengthMafStatistics();
          } else if (statName == "SequenceLength") {
            string sp = ApplicationTools::getStringParameter("species", statArgs, "");
            mafStat = new SequenceLengthMafStatistics(sp);
          } else if (statName == "AlnScore") {
            mafStat = new AlignmentScoreMafStatistics();
          } else if (statName == "BlockCounts") {
            mafStat = new CharacterCountsMafStatistics(&AlphabetTools::DNA_ALPHABET);
          } else if (statName == "PairwiseDivergence") {
            string sp1 = ApplicationTools::getStringParameter("species1", statArgs, "");
            string sp2 = ApplicationTools::getStringParameter("species2", statArgs, "");
            mafStat = new PairwiseDivergenceMafStatistics(sp1, sp2);
          } else if (statName == "SiteFrequencySpectrum") {
            vector<double> bounds  = ApplicationTools::getVectorParameter<double>("bounds", statArgs, ',', "", "", false, true);
            vector<string> ingroup = ApplicationTools::getVectorParameter<string>("ingroup", statArgs, ',', "", "", false, true);
            mafStat = new SiteFrequencySpectrumMafStatistics(&AlphabetTools::DNA_ALPHABET, bounds, ingroup); 
          } else if (statName == "SiteStatistics") {
            vector<string> species = ApplicationTools::getVectorParameter<string>("species", statArgs, ',', "", "", false, true);
            mafStat = new SiteMafStatistics(species); 
          } else if (statName == "CountClusters") {
            string treeProperty = ApplicationTools::getStringParameter("tree", statArgs, "none");
            double threshold = ApplicationTools::getDoubleParameter("threshold", statArgs, 0);
            mafStat = new CountClustersMafStatistics(treeProperty, threshold);
            statDesc = " / " + treeProperty;
           } else {
            throw Exception("Unknown statistic: " + statName);
          }
          statistics.push_back(mafStat);
          ApplicationTools::displayResult("-- Adding statistic", mafStat->getFullName() + " <" + mafStat->getShortName() + ">" + statDesc);
        }

        //Get output file:
        string outputFile = ApplicationTools::getAFilePath("file", cmdArgs, true, false);
        ApplicationTools::displayResult("-- Output file", outputFile);
        auto_ptr<ostream> ofs(new ofstream(outputFile.c_str(), ios::out));
        StlOutputStream* output = new StlOutputStream(ofs);
        SequenceStatisticsMafIterator* iterator = new SequenceStatisticsMafIterator(currentIterator, statistics);
        
        if (cmdArgs.find("reference") != cmdArgs.end()) {
          throw Exception("reference argument in SequenceStatistics is deprecated: use ref_species instead.");
        }
        string ref = ApplicationTools::getStringParameter("ref_species", cmdArgs, "none");
        ApplicationTools::displayResult("-- Reference species", ref);
        CsvStatisticsOutputIterationListener* listener = new CsvStatisticsOutputIterationListener(iterator, ref, output);
        
        iterator->addIterationListener(listener);
        currentIterator = iterator;
        iterator->verbose(verbose);
        its.push_back(iterator);
      }

      
      // +--------------------+
      // | Feature extraction |
      // +--------------------+
      if (cmdName == "ExtractFeature") {
        bool ignoreStrand    = ApplicationTools::getBooleanParameter("ignore_strand", cmdArgs, false);
        string refSpecies    = ApplicationTools::getStringParameter("ref_species", cmdArgs, "none");
        string featureFile   = ApplicationTools::getAFilePath("feature.file", cmdArgs, false, false);
        string featureFormat = ApplicationTools::getStringParameter("feature.format", cmdArgs, "GFF");
        ApplicationTools::displayResult("-- Features to extract", featureFile + " (" + featureFormat + ")");
        ApplicationTools::displayResult("-- Features are for species", refSpecies);
        ApplicationTools::displayBooleanResult("-- Features are strand-aware", !ignoreStrand);
        //string featuresList = ApplicationTools::getStringParameter("feature.select", cmdArgs, "none");
        compress = ApplicationTools::getStringParameter("feature.file.compression", cmdArgs, "none");
        filtering_istream featureStream;
        if (compress == "none") {
        } else if (compress == "gzip") {
          featureStream.push(gzip_decompressor());
        } else if (compress == "zip") {
          featureStream.push(zlib_decompressor());
        } else if (compress == "bzip2") {
          featureStream.push(bzip2_decompressor());
        } else
          throw Exception("Bad input incompression format: " + compress);
        featureStream.push(file_source(featureFile));

        SequenceFeatureSet featuresSet;
        if (featureFormat == "GFF") {
          GffFeatureReader reader(featureStream);
          reader.getAllFeatures(featuresSet);
        } else if (featureFormat == "GTF") {
          GtfFeatureReader reader(featureStream);
          reader.getAllFeatures(featuresSet);
        } else
          throw Exception("Unsupported feature format: " + featureFormat);
        
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
        ApplicationTools::displayResult("-- Features to extract", lst);
        auto_ptr<SequenceFeatureSet> featuresSet2(featuresSet.getSubsetForType(features));

        FeatureExtractor* iterator = new FeatureExtractor(currentIterator, refSpecies, *featuresSet2, ignoreStrand);
        iterator->setLogStream(&log);
        iterator->verbose(verbose);
        its.push_back(iterator);

        currentIterator = iterator;
      }



      // +------------------+
      // | Window splitting |
      // +------------------+
      if (cmdName == "WindowSplit") {
        if (cmdArgs.find("preferred.size") != cmdArgs.end()) {
          throw Exception("preferred.size argument in WindowSplit is deprecated: use preferred_size instead.");
        }
        unsigned int preferredSize = ApplicationTools::getParameter<unsigned int>("preferred_size", cmdArgs, 0);
        ApplicationTools::displayResult("-- Preferred size", preferredSize);
        string splitOptionStr = ApplicationTools::getStringParameter("align", cmdArgs, "center");
        short splitOption;
        if (splitOptionStr == "ragged_left")
          splitOption = WindowSplitMafIterator::RAGGED_LEFT;
        else if (splitOptionStr == "ragged_right")
          splitOption = WindowSplitMafIterator::RAGGED_RIGHT;
        else if (splitOptionStr == "center")
          splitOption = WindowSplitMafIterator::CENTER;
        else if (splitOptionStr == "adjust")
          splitOption = WindowSplitMafIterator::ADJUST;
        else throw Exception("Unvalid alignment option for WindowSplit: " + splitOptionStr);
        ApplicationTools::displayResult("-- Alignment option", splitOptionStr);

        WindowSplitMafIterator* iterator = new WindowSplitMafIterator(currentIterator, preferredSize, splitOption);
        iterator->setLogStream(&log);
        currentIterator = iterator;
        its.push_back(iterator);
      }



      // +---------------------+
      // | Distance estimation |
      // +---------------------+
      if (cmdName == "DistanceEstimation") {
        string distMethod = ApplicationTools::getStringParameter("method", cmdArgs, "count");
        ApplicationTools::displayResult("-- Method", distMethod);
        if (distMethod == "count") {
          string gapOption = ApplicationTools::getStringParameter("gap_option", cmdArgs, "no_gap");
          if (gapOption == "all") {
            gapOption = SiteContainerTools::SIMILARITY_ALL;
          } else if (gapOption == "no_gap") {
            gapOption = SiteContainerTools::SIMILARITY_NOGAP;
          } else if (gapOption == "no_full_gap") {
            gapOption = SiteContainerTools::SIMILARITY_NOFULLGAP;
          } else if (gapOption == "no_double_gap") {
            gapOption = SiteContainerTools::SIMILARITY_NODOUBLEGAP;
          } else {
            throw Exception("Unrecognized gap option, should be either 'all', 'no_full_gap', 'no_double_gap' or 'no_gap'.");
          }
          ApplicationTools::displayResult("-- Gap option", gapOption);
          bool unresolvedAsGap = ApplicationTools::getBooleanParameter("unresolved_as_gap", cmdArgs, "no");
          ApplicationTools::displayBooleanResult("-- Unresolved as gaps", unresolvedAsGap);

          CountDistanceEstimationMafIterator* iterator = new CountDistanceEstimationMafIterator(currentIterator, gapOption, unresolvedAsGap);
          ApplicationTools::displayResult("-- Block-wise matrices are registered as", iterator->getPropertyName());
          iterator->setLogStream(&log);
          currentIterator = iterator;
          its.push_back(iterator);
        } else {
          throw Exception("Unknown distance method: " + distMethod);
        }
      }



      // +--------------------------+
      // | Phylogeny reconstruction |
      // +--------------------------+
      if (cmdName == "DistanceBasedPhylogeny") {
        string distMethodName = ApplicationTools::getStringParameter("method", cmdArgs, "bionj");
        string distProperty = ApplicationTools::getStringParameter("dist_mat", cmdArgs, "none");
        DistanceMethod* distMethod = 0;
        if (distMethodName == "upgma") {
          distMethod = new PGMA(false);
        } else if (distMethodName == "wpgma") {
          distMethod = new PGMA(true);
        } else if (distMethodName == "nj") {
          distMethod = new NeighborJoining(false, false);
        } else if (distMethodName == "bionj") {
          distMethod = new BioNJ(false, false);
        } else {
          throw Exception("Unknown distance-based phylogenetic method: " + distMethodName); 
        }
        distMethod->setVerbose(false);
        ApplicationTools::displayResult("-- Reading distance matrix from", distProperty);
        ApplicationTools::displayResult("-- Build distance tree using", distMethodName);

        DistanceBasedPhylogenyReconstructionMafIterator* iterator = new DistanceBasedPhylogenyReconstructionMafIterator(currentIterator, distMethod, distProperty);
        ApplicationTools::displayResult("-- Writing block-wise trees to", iterator->getPropertyName());
        iterator->setLogStream(&log);
        currentIterator = iterator;
        its.push_back(iterator);
      }



      // +-------------------+
      // | Phylogeny rooting |
      // +-------------------+
      if (cmdName == "NewOutgroup") {
        string treePropertyInput = ApplicationTools::getStringParameter("tree_input", cmdArgs, "none");
        string treePropertyOutput = ApplicationTools::getStringParameter("tree_output", cmdArgs, "none");
        string outgroup = ApplicationTools::getStringParameter("outgroup", cmdArgs, "none");
        ApplicationTools::displayResult("-- Reading tree from", treePropertyInput);
        ApplicationTools::displayResult("-- Rerooting according to species", outgroup);
        NewOutgroupMafIterator* iterator = new NewOutgroupMafIterator(currentIterator, treePropertyInput, treePropertyOutput, outgroup);
        ApplicationTools::displayResult("-- Writing tree to", treePropertyOutput);
        iterator->setLogStream(&log);
        currentIterator = iterator;
        its.push_back(iterator);
      }



      // +------------------------+
      // | Phylogeny drop species |
      // +------------------------+
      if (cmdName == "DropSpecies") {
        string treePropertyInput = ApplicationTools::getStringParameter("tree_input", cmdArgs, "none");
        string treePropertyOutput = ApplicationTools::getStringParameter("tree_output", cmdArgs, "none");
        string species = ApplicationTools::getStringParameter("species", cmdArgs, "none");
        ApplicationTools::displayResult("-- Reading tree from", treePropertyInput);
        ApplicationTools::displayResult("-- Removing leaves from species", species);
        DropSpeciesMafIterator* iterator = new DropSpeciesMafIterator(currentIterator, treePropertyInput, treePropertyOutput, species);
        ApplicationTools::displayResult("-- Writing tree to", treePropertyOutput);
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
        ApplicationTools::displayResult("-- Output file", outputFile);
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
        ApplicationTools::displayResult("-- File compression", compress);
        bool mask = ApplicationTools::getBooleanParameter("mask", cmdArgs, true);
        ApplicationTools::displayBooleanResult("-- Output mask", mask);
        OutputMafIterator* iterator = new OutputMafIterator(currentIterator, out, mask);
        currentIterator = iterator;
        its.push_back(iterator);
      }



      // +-------------------+
      // | Output alignments |
      // +-------------------+
      if (cmdName == "OutputAlignments") {
        string outputFile = ApplicationTools::getAFilePath("file", cmdArgs, true, false);
        bool multipleFiles = (outputFile.find("%i") != string::npos);
        ApplicationTools::displayResult("-- Output alignment file" + string(multipleFiles ? "s" : ""), outputFile);
        bool mask = ApplicationTools::getBooleanParameter("mask", cmdArgs, true);
        ApplicationTools::displayBooleanResult("-- Output mask", mask);
        
        OutputAlignmentMafIterator* iterator; 
        BppOAlignmentWriterFormat bppoWriter;
        string description = ApplicationTools::getStringParameter("format", cmdArgs, "Clustal");
        map<string, string> fargs;
        OAlignment* oAln = bppoWriter.read(description, fargs, true);
        if (multipleFiles) {
          iterator = new OutputAlignmentMafIterator(currentIterator, outputFile, oAln, mask);
        } else {
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
          ApplicationTools::displayResult("-- File compression", compress);
          iterator = new OutputAlignmentMafIterator(currentIterator, out, oAln, mask);
        }
        currentIterator = iterator;
        its.push_back(iterator);
      }



      // +--------------+
      // | Output trees |
      // +--------------+
      if (cmdName == "OutputTrees") {
        string outputFile = ApplicationTools::getAFilePath("file", cmdArgs, true, false);
        compress = ApplicationTools::getStringParameter("compression", cmdArgs, "none");
        ApplicationTools::displayResult("-- Output tree file", outputFile);
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
        ApplicationTools::displayResult("-- File compression", compress);
        string treeProperty = ApplicationTools::getStringParameter("tree", cmdArgs, "none");
        ApplicationTools::displayResult("-- Tree to write", treeProperty);
        OutputTreeMafIterator* iterator = new OutputTreeMafIterator(currentIterator, out, treeProperty);
        currentIterator = iterator;
        its.push_back(iterator);
      }

    }

    //Now loop over the last iterator and that's it!
    size_t blockCounter = 0;
    while (MafBlock* block = currentIterator->nextBlock())
    {
      ApplicationTools::displayUnlimitedGauge(blockCounter++, "Parsing...");
      delete block;
    }
    ApplicationTools::message->endLine();

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

