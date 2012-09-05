//
// File: MafIterators.cpp
// Created by: Julien Dutheil
// Created on: Jul 24 2012
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

#include "MafIterators.h"

#include <Bpp/Seq/Container/SiteContainerTools.h>

DistanceMatrix* CountDistanceEstimationMafIterator::estimateDistanceMatrixForBlock(const MafBlock& block)
{
  DistanceMatrix* dist = SiteContainerTools::computeSimilarityMatrix(block.getAlignment(), true, gapOption_, unresolvedAsGap_);
  return dist;
}

Tree* DistanceBasedPhylogenyReconstructionMafIterator::buildTreeForBlock(const MafBlock& block) throw (Exception)
{
  //First get the distance matrix for this block:
  if (!block.hasProperty(distanceProperty_))
    throw Exception("DistanceBasedPhylogenyReconstructionMafIterator::buildTreeForBlock. No property available for " + distanceProperty_);
  try {
    const DistanceMatrix& dist = dynamic_cast<const DistanceMatrix&>(block.getProperty(distanceProperty_));
    builder_->setDistanceMatrix(dist);
    builder_->computeTree(false);
    Tree* tree = builder_->getTree();
    if (!tree)
      throw Exception("DistanceBasedPhylogenyReconstructionMafIterator::buildTreeForBlock. Tree reconstruction failed!");
    return tree;
  } catch (bad_cast& e) {
    throw Exception("DistanceBasedPhylogenyReconstructionMafIterator::buildTreeForBlock. A property was found for '" + distanceProperty_ + "' but does not appear to contain a distance matrix.");
  }
}

void OutputTreeMafIterator::writeBlock(std::ostream& out, const MafBlock& block) const
{
  //First get the tree for this block:
  if (!block.hasProperty(treeProperty_))
    throw Exception("OutputTreeAlignmentMafIterator::writeBlock. No property available for " + treeProperty_);
  try {
    const Tree& tree = dynamic_cast<const Tree&>(block.getProperty(treeProperty_));
    writer_.write(tree, out);
  } catch (Exception& e) {
    throw Exception("OutputTreeAlignmentMafIterator::writeBlock. A property was found for '" + treeProperty_ + "' but does not appear to contain a phylogenetic tree.");
  }

}

unsigned int CountClustersMafStatistics::getNumberOfClusters_(const Node* node, map<const Node*, double>& heights)
{
  unsigned int nClust = 0;
  double h = heights[node];
  if (h < threshold_) {
    nClust++;
  } else {
    for (unsigned int i = 0; i < node->getNumberOfSons(); ++i) {
      nClust += getNumberOfClusters_((*node)[i], heights);
    }
  }
  return nClust;
}

void CountClustersMafStatistics::compute(const MafBlock& block)
{
  if (!block.hasProperty(treeProperty_))
    throw Exception("CountClustersMafStatistics::compute. No property available for " + treeProperty_);
  try {
    TreeTemplate<Node> tree(dynamic_cast<const Tree&>(block.getProperty(treeProperty_)));
    if (!tree.isRooted())
      throw Exception("CountClustersMafStatistics::compute. Cluster count only works with a rooted tree.");
    //Compute all tree heights:
    map<const Node*, double> heights;
    TreeTemplateTools::getHeights(*tree.getRootNode(), heights);
    unsigned int nClust = getNumberOfClusters_(tree.getRootNode(), heights);
    result_.setValue("NbClusters", nClust);
  } catch (bad_cast& e) {
    throw Exception("CountClustersMafStatistics::compute. A property was found for '" + treeProperty_ + "' but does not appear to contain a phylogenetic tree.");
  }
}


MafBlock* TreeManipulationMafIterator::analyseCurrentBlock_() throw (Exception)
{
  currentBlock_ = iterator_->nextBlock();
  if (currentBlock_) {
    if (!currentBlock_->hasProperty(treePropertyRead_))
      throw Exception("TreeManipulationMafIterator::analyseCurrentBlock_(). No property available for " + treePropertyRead_);
    try {
      TreeTemplate<Node>* tree = new TreeTemplate<Node>(dynamic_cast<const Tree&>(currentBlock_->getProperty(treePropertyRead_)));
      manipulateTree_(tree);
      currentBlock_->setProperty(treePropertyWrite_, tree);
    } catch (Exception& e) {
      throw Exception("TreeManipulationMafIterator::analyseCurrentBlock_(). A property was found for '" + treePropertyRead_ + "' but does not appear to contain a phylogenetic tree.");
    }
  }
  return currentBlock_;
}


void NewOutgroupMafIterator::manipulateTree_(TreeTemplate<Node>* tree) throw (Exception)
{
  vector<Node*> leaves = tree->getLeaves();
  Node* outgroup = 0;
  bool outgroupFound = false;
  for (size_t i = 0; i < leaves.size() && !outgroupFound; ++i) {
    string species, chr;
    MafSequence::splitNameIntoSpeciesAndChromosome(leaves[i]->getName(), species, chr);
    if (species == outgroupSpecies_) {
      outgroup = leaves[i];
      outgroupFound = true;
    }
  }
  if (!outgroupFound)
    throw Exception("NewOutgroupTreeMafIterator::analyseCurrentBlock_(). No ougroup species was found in the attached tree.");
  tree->newOutGroup(outgroup);
}


void DropSpeciesMafIterator::manipulateTree_(TreeTemplate<Node>* tree) throw (Exception)
{
  vector<Node*> leaves = tree->getLeaves();
  for (size_t i = 0; i < leaves.size(); ++i) {
    string species, chr;
    MafSequence::splitNameIntoSpeciesAndChromosome(leaves[i]->getName(), species, chr);
    if (species == species_) {
      TreeTemplateTools::dropSubtree(*tree, leaves[i]);
    }
  }
}

