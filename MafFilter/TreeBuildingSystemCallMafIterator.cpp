//
// File: TreeBuildingSystemCallMafIterator.cpp
// Authors: Julien Dutheil
// Created: Sat Jun 18 2016
//

/*
Copyright or © or Copr. Julien Y. Dutheil, (2016)

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

#include "TreeBuildingSystemCallMafIterator.h"

// From bpp-phyl
#include "Bpp/Phyl/TreeTemplate.h"

using namespace bpp;
using namespace std;

MafBlock* TreeBuildingSystemCallMafIterator::analyseCurrentBlock_() throw (Exception) {
  currentBlock_ = iterator_->nextBlock();
  if (! currentBlock_)
    return 0;
  unique_ptr<AlignedSequenceContainer> aln(currentBlock_->getAlignment().clone());
  
  //We translate sequence names to avoid compatibility issues
  vector<string> names(aln->getNumberOfSequences());
  map<string, string> nameIndex;
  for (size_t i = 0; i < names.size(); ++i) {
    names[i] = "seq" + TextTools::toString(i);
    nameIndex[names[i]] = currentBlock_->getSequence(i).getName();
  }
  aln->setSequencesNames(names);
  
  //Write sequences to file:
  alnWriter_->writeAlignment(inputFile_, *aln, true);

  //Call the external program:
  int rc = system(call_.c_str());
  if (rc) throw Exception("TreeBuildingSystemCallMafIterator::analyseCurrentBlock_(). System call exited with non-zero status.");

  //Then read the generated tree and assign sequence names:
  unique_ptr< Tree > result(treeReader_->read(outputFile_));
  unique_ptr< TreeTemplate<Node> > tree(new TreeTemplate<Node>(*result));
  vector<Node*> leaves = tree->getLeaves();
  for (size_t i = 0; i < leaves.size(); ++i) {
    leaves[i]->setName(nameIndex[leaves[i]->getName()]);
  }
  currentBlock_->setProperty(propertyName_, tree.release());

  //Done:
  return currentBlock_;
}

