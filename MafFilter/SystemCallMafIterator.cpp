//
// File: SystemCallMafIterator.cpp
// Authors: Julien Dutheil
// Created: Tue Sep 29 2015
//

/*
Copyright or © or Copr. Julien Y. Dutheil, (2015)

This file is part of MafFilter.

MafFilter is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MafFilter is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MafFilter.  If not, see <https://www.gnu.org/licenses/>.
*/

#include "SystemCallMafIterator.h"

#include <Bpp/Seq/Container/VectorSiteContainer.h>
#include <Bpp/Seq/Container/SiteContainerTools.h>
#include <Bpp/Seq/SequenceWithQuality.h>

using namespace bpp;
using namespace std;

unique_ptr<MafBlock> SystemCallMafIterator::analyseCurrentBlock_() {
  currentBlock_ = iterator_->nextBlock();
  if (! currentBlock_)
    return 0;
  auto aln = currentBlock_->getAlignment();
  
  //We translate sequence names to avoid compatibility issues
  vector<string> names(aln->getNumberOfSequences());
  for (size_t i = 0; i < names.size(); ++i) {
    names[i] = "seq" + TextTools::toString(i);
  }
  aln->setSequenceNames(names, true);
  
  //Write sequences to file:
  alnWriter_->writeAlignment(inputFile_, *aln, true);

  //Call the external program:
  int rc = system(call_.c_str());
  if (rc) throw Exception("SystemCallMafIterator::analyseCurrentBlock_(). System call exited with non-zero status.");
  
  //Read the results:
  auto result = alnReader_->readAlignment(outputFile_, AlphabetTools::DNA_ALPHABET);

  //Convert and assign the realigned sequences:
  vector<unique_ptr<MafSequence>> tmp;
  for (size_t i = 0; i < currentBlock_->getNumberOfSequences(); ++i) {
    unique_ptr<MafSequence> mseq(currentBlock_->sequence(i).cloneMeta());
    //NB: we discard any putative score associated to this sequence.
    string seqKey = "seq" + TextTools::toString(i);
    mseq->setContent(dynamic_cast<const Sequence&>(result->sequence(seqKey)).toString()); //NB shall we use getContent here?
    tmp.push_back(std::move(mseq));
  }
  currentBlock_->clear();
  for (size_t i = 0; i < tmp.size(); ++i) {
    currentBlock_->addSequence(tmp[i]);
  }

  //Perform a Head-or-Tail test to compute alignment scores (optional):
  if (hotTest_) {
    //Reverse the alignment:
    auto rev = make_unique<AlignedSequenceContainer>(AlphabetTools::DNA_ALPHABET);
    for (size_t i = 0; i < aln->getNumberOfSequences(); ++i) {
      auto tmpSeq = SequenceTools::getInvert(aln->sequence(i));
      auto invSeq = unique_ptr<Sequence>(dynamic_cast<Sequence*>(tmpSeq.release()));
      rev->addSequence(invSeq->getName(), invSeq);
    }
    //Write reversed sequences to file (so far we use the same file as for the non-reversed alignment):
    alnWriter_->writeAlignment(inputFile_, *rev, true);

    //Call the external program again:
    rc = system(call_.c_str());
    if (rc) throw Exception("SystemCallMafIterator::analyseCurrentBlock_(). System call exited with non-zero status.");
    
    //Read the results:
    auto result2 = alnReader_->readAlignment(outputFile_, AlphabetTools::DNA_ALPHABET);

    //Re-reverse the results, and make sure they are in the same order as the non-inverted results:
    auto result2rev = make_unique<VectorSiteContainer>(result->getAlphabet());
    vector<string> seqNames = result->getSequenceNames();
    for (const auto& i : seqNames) {
      auto tmpSeq = SequenceTools::getInvert(result2->sequence(i));
      auto invSeq = unique_ptr<Sequence>(dynamic_cast<Sequence*>(tmpSeq.release()));
      result2rev->addSequence(invSeq->getName(), invSeq);
    }

    //Now we compare the two alignments:
    // Build alignment indexes:
    RowMatrix<size_t> indexTest, indexRef;
    SiteContainerTools::getSequencePositions(*result2rev, indexTest);
    SiteContainerTools::getSequencePositions(*result, indexRef);

    // Now build scores:
    vector<double> sps = SiteContainerTools::getSumOfPairsScores(indexRef, indexTest, 1);
    double nbPos = 0.;
    double sumPos = 0.;
    for (size_t i = 0; i < sps.size(); ++i) {
      if (sps[i] >= 0.) {
	nbPos ++;
	sumPos += sps[i];
      }
    }

    // Assign global score:
    double s = nbPos > 0 ? sumPos / nbPos : 0.;
    currentBlock_->setScore(s);

    // Record and assign site-specific scores:
    if (refSeq_ != "" && currentBlock_->hasSequenceForSpecies(refSeq_)) 
    {
      auto seqQual = make_shared<SequenceQuality>(sps.size());
      for (size_t i = 0; i < sps.size(); ++i) {
	if (sps[i] < 0) {
	  seqQual->setScore(i, -1);
	}
	else if (sps[i] == 1.)
	{
	  seqQual->setScore(i, 10);
	}
	else
	{
	  seqQual->setScore(i, static_cast<int>(floor(sps[i] * 10)));
	}
      }
      currentBlock_->addAnnotationToSequenceForSpecies(refSeq_, seqQual); 
    }
  }

  //Done:
  return std::move(currentBlock_);
}

