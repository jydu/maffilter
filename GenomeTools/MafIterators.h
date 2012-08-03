//
// File: MafIterators.h
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

#include <Bpp/Seq/Io/Maf/MafIterator.h>
#include <Bpp/Seq/Container/SiteContainer.h>
#include <Bpp/Phyl/Tree.h>
#include <Bpp/Phyl/Distance/DistanceEstimation.h>
#include <Bpp/Phyl/Distance/PGMA.h>
#include <Bpp/Phyl/Io/Newick.h>

using namespace bpp;

/**
 * @brief Partial implementation for distance estimation iterator.
 *
 * This iterator calls a distance reconstruction method (to be implemented by the derivated class)
 * and store the resulting distance matrix as an associated block property for the block,
 * before forwarding it.
 */
class AbstractDistanceEstimationMafIterator:
  public AbstractFilterMafIterator
{
  public:
    AbstractDistanceEstimationMafIterator(MafIterator* iterator):
      AbstractFilterMafIterator(iterator)
    {}

  private:
    MafBlock* analyseCurrentBlock_() throw (Exception)
    {
      MafBlock* block = iterator_->nextBlock();
      if (!block) return 0;
      DistanceMatrix* dist = estimateDistanceMatrixForBlock(*block);
      block->setProperty(getPropertyName(), dist);
      return block;
    }

  public:
    virtual std::string getPropertyName() const = 0;
    virtual DistanceMatrix* estimateDistanceMatrixForBlock(const MafBlock& block) = 0;

};


/**
 * @brief Compute A simple distance using observed counts.
 */
class CountDistanceEstimationMafIterator:
  public AbstractDistanceEstimationMafIterator
{
  private:
    std::string gapOption_;
    bool unresolvedAsGap_;

  public:
    /**
     * @brief Build a new distance estimation maf iterator, based on the SiteContainerTools::computeSimilarityMatrix method.
     *
     * @see SiteContainerTools
     * @param gapOption How to deal with gaps. Option forawarded to the computeSimilarityMatrix method.
     * @param unresolvedAsGap Tell if unresolved characters should be considered as gaps. Option forawarded to the computeSimilarityMatrix method.
     */
    CountDistanceEstimationMafIterator(MafIterator* iterator, const std::string& gapOption, bool unresolvedAsGap):
      AbstractDistanceEstimationMafIterator(iterator),
      gapOption_(gapOption), unresolvedAsGap_(unresolvedAsGap)
    {}
    
  public:
    std::string getPropertyName() const { return "CountDistance"; }
    DistanceMatrix* estimateDistanceMatrixForBlock(const MafBlock& block);

};



/**
 * @brief Partial implementation for phylogeny reconstruction iterator.
 *
 * This iterator calls a tree reconstruction method (to be implemented by the derivated class)
 * and store the resulting tree as an associated block property for the block,
 * before forwarding it.
 */
class AbstractPhylogenyReconstructionMafIterator:
  public AbstractFilterMafIterator
{
  public:
    AbstractPhylogenyReconstructionMafIterator(MafIterator* iterator):
      AbstractFilterMafIterator(iterator)
    {}

  virtual ~AbstractPhylogenyReconstructionMafIterator() {}

  private:
    MafBlock* analyseCurrentBlock_() throw (Exception)
    {
      MafBlock* block = iterator_->nextBlock();
      if (!block) return 0;
      Tree* tree = buildTreeForBlock(*block);
      block->setProperty(getPropertyName(), tree);
      return block;
    }

  public:
    virtual std::string getPropertyName() const = 0;
    virtual Tree* buildTreeForBlock(const MafBlock& block) = 0;

};



/**
 * @brief Implementation for distance-based phylogeny reconstruction iterator.
 */
class DistanceBasedPhylogenyReconstructionMafIterator:
  public AbstractPhylogenyReconstructionMafIterator
{
  private:
    std::string distanceProperty_;
    std::auto_ptr<DistanceMethod> builder_;
  
  public:
    DistanceBasedPhylogenyReconstructionMafIterator(MafIterator* iterator, DistanceMethod* method, const std::string& property):
      AbstractPhylogenyReconstructionMafIterator(iterator),
      distanceProperty_(property), builder_(method)
    {}

  private:
    DistanceBasedPhylogenyReconstructionMafIterator(const DistanceBasedPhylogenyReconstructionMafIterator& it):
      AbstractPhylogenyReconstructionMafIterator(0),
      distanceProperty_(it.distanceProperty_),
      builder_()
    {}

    DistanceBasedPhylogenyReconstructionMafIterator& operator=(const DistanceBasedPhylogenyReconstructionMafIterator& it)
    {
      distanceProperty_ = it.distanceProperty_;
      return *this;
    }

  public:
    void setDistanceProperty(const std::string& property) { distanceProperty_ = property; }
    const std::string& getDistanceProperty() const { return distanceProperty_; }

    std::string getPropertyName() const { return builder_->getName(); }
    Tree* buildTreeForBlock(const MafBlock& block) throw (Exception);
};


/**
 * @brief Count the number of sequence cluster, given a certain threshold.
 * This requires that a phylogenetic tree was previously computed.
 */
class CountClustersMafStatistics:
  public AbstractMafStatisticsSimple
{
  private:
    std::string treeProperty_;
    double threshold_;
  
  public:
    CountClustersMafStatistics(const std::string& property, double threshold):
      AbstractMafStatisticsSimple("NbClusters"),
      treeProperty_(property), threshold_(threshold)
    {}

  public:
    void setTreeProperty(const std::string& property) { treeProperty_ = property; }
    const std::string& getTreeProperty() const { return treeProperty_; }

    std::string getShortName() const { return "CountClusters(" + TextTools::toString(threshold_) + ")"; }
    std::string getFullName() const { return "Number of sequence clusters with divergence <= " + TextTools::toString(threshold_) + "."; }
    void compute(const MafBlock& block);

    std::string getPropertyName() const { return "CountClusters"; }

  private:
    unsigned int getNumberOfClusters_(const Node* node, map<const Node*, double>& heights);
};



/**
 * @brief This iterator root associated trees according to an outgroup sequence.
 */
class TreeManipulationMafIterator:
  public AbstractFilterMafIterator
{
  private:
    std::string treePropertyRead_;
    std::string treePropertyWrite_;

  public:
    //Write can be the same as read.
    TreeManipulationMafIterator(MafIterator* iterator, const std::string& treePropertyRead, const std::string& treePropertyWrite) :
      AbstractFilterMafIterator(iterator), treePropertyRead_(treePropertyRead), treePropertyWrite_(treePropertyWrite)
    {}

  public:
    MafBlock* analyseCurrentBlock_() throw (Exception);

  protected:
    virtual void manipulateTree_(TreeTemplate<Node>* tree) throw (Exception) = 0;

};



/**
 * @brief This iterator root associated trees according to an outgroup sequence.
 */
class NewOutgroupMafIterator:
  public TreeManipulationMafIterator
{
  private:
    std::string outgroupSpecies_;

  public:
    //Write can be the same as read.
    NewOutgroupMafIterator(MafIterator* iterator, const std::string& treePropertyRead, const std::string& treePropertyWrite, const std::string& outgroupSpecies) :
      TreeManipulationMafIterator(iterator, treePropertyRead, treePropertyWrite), outgroupSpecies_(outgroupSpecies)
    {}

  private:
    void manipulateTree_(TreeTemplate<Node>* tree) throw (Exception);

};



/**
 * @brief This iterator removes leaves of a certain species in an attached tree.
 */
class DropSpeciesMafIterator:
  public TreeManipulationMafIterator
{
  private:
    std::string species_;

  public:
    //Write can be the same as read.
    DropSpeciesMafIterator(MafIterator* iterator, const std::string& treePropertyRead, const std::string& treePropertyWrite, const std::string& species) :
      TreeManipulationMafIterator(iterator, treePropertyRead, treePropertyWrite), species_(species)
    {}

  private:
    void manipulateTree_(TreeTemplate<Node>* tree) throw (Exception);

};



/**
 * @brief This iterator print an attached tree to a newick file.
 */
class OutputTreeMafIterator:
  public AbstractFilterMafIterator
{
  private:
    std::ostream* output_;
    std::string treeProperty_;
    Newick writer_;

  public:
    OutputTreeMafIterator(MafIterator* iterator, std::ostream* out, const std::string treeProperty) :
      AbstractFilterMafIterator(iterator), output_(out), treeProperty_(treeProperty), writer_()
    {}

  private:
    OutputTreeMafIterator(const OutputTreeMafIterator& iterator) :
      AbstractFilterMafIterator(0),
      output_(iterator.output_),
      treeProperty_(iterator.treeProperty_),
      writer_()
    {}
    
    OutputTreeMafIterator& operator=(const OutputTreeMafIterator& iterator)
    {
      output_ = iterator.output_;
      treeProperty_ = iterator.treeProperty_;
      writer_ = iterator.writer_;
      return *this;
    }


  public:
    MafBlock* analyseCurrentBlock_() throw (Exception) {
      currentBlock_ = iterator_->nextBlock();
      if (output_ && currentBlock_)
        writeBlock(*output_, *currentBlock_);
      return currentBlock_;
    }

  private:
    void writeBlock(std::ostream& out, const MafBlock& block) const;
};



