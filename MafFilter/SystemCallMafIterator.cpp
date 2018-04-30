//
// File: SystemCallMafIterator.cpp
// Authors: Julien Dutheil
// Created: Tue Sep 29 2015
//

/*
Copyright or © or Copr. Julien Y. Dutheil, (2015)

This software is a computer program whose purpose is to provide classes
for sequences analysis.

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

#include "SystemCallMafIterator.h"

using namespace bpp;
using namespace std;

MafBlock* SystemCallMafIterator::analyseCurrentBlock_() {
  currentBlock_ = iterator_->nextBlock();
  if (! currentBlock_)
    return 0;
  unique_ptr<AlignedSequenceContainer> aln(currentBlock_->getAlignment().clone());
  
  //We translate sequence names to avoid compatibility issues
  vector<string> names(aln->getNumberOfSequences());
  for (size_t i = 0; i < names.size(); ++i) {
    names[i] = "seq" + TextTools::toString(i);
  }
  aln->setSequencesNames(names);
  
  //Write sequences to file:
  alnWriter_->writeAlignment(inputFile_, *aln, true);

  //Call the external program:
  int rc = system(call_.c_str());
  if (rc) throw Exception("SystemCallMafIterator::analyseCurrentBlock_(). System call exited with non-zero status.");

  //Then read and assign the realigned sequences:
  unique_ptr<SiteContainer> result(alnReader_->readAlignment(outputFile_, &AlphabetTools::DNA_ALPHABET));
  vector<MafSequence*> tmp;
  for (size_t i = 0; i < currentBlock_->getNumberOfSequences(); ++i) {
    MafSequence* mseq = currentBlock_->getSequence(i).cloneMeta();
    //NB: we discard any putative score associated to this sequence.
    string name = "seq" + TextTools::toString(i);
    mseq->setContent(dynamic_cast<const BasicSequence&>(result->getSequence(name)).toString()); //NB shall we use getContent here?
    tmp.push_back(mseq);
  }
  currentBlock_->getAlignment().clear();
  for (size_t i = 0; i < tmp.size(); ++i) {
    currentBlock_->getAlignment().addSequence(*tmp[i], false);
    delete tmp[i];
  }

  //Done:
  return currentBlock_;
}

