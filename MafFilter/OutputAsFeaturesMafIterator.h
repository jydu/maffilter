//
// File: OutputasFeaturesMafIterator.h
// Authors: Julien Dutheil
// Created: Fri Mar 22 2013
//

/*
Copyright or © or Copr. Bio++ Development Team, (2013)

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

#ifndef _OUTPUTASFEATURESMAFITERATOR_H_
#define _OUTPUTASFEATURESMAFITERATOR_H_

#include <Bpp/Seq/Io/Maf/MafIterator.h>

//From the STL:
#include <iostream>
#include <string>
#include <deque>

namespace bpp {

/**
 * @brief This iterator print as features the sequences for a given species.
 */
class OutputAsFeaturesMafIterator:
  public AbstractFilterMafIterator
{
  private:
    std::ostream* output_;
    std::string species_;

  public:
    OutputAsFeaturesMafIterator(MafIterator* iterator, std::ostream* out, const std::string& species) :
      AbstractFilterMafIterator(iterator), output_(out), species_(species)
    {
      if (output_)
        writeHeader(*output_);
    }

  private:
    OutputAsFeaturesMafIterator(const OutputAsFeaturesMafIterator& iterator) :
      AbstractFilterMafIterator(0),
      output_(iterator.output_),
      species_(iterator.species_)
    {}
    
    OutputAsFeaturesMafIterator& operator=(const OutputAsFeaturesMafIterator& iterator)
    {
      output_ = iterator.output_;
      species_ = iterator.species_;
      return *this;
    }


  public:
    MafBlock* analyseCurrentBlock_() {
      currentBlock_ = iterator_->nextBlock();
      if (output_ && currentBlock_)
        writeBlock(*output_, *currentBlock_);
      return currentBlock_;
    }

  private:
    void writeHeader(std::ostream& out) const;
    void writeBlock(std::ostream& out, const MafBlock& block) const;
};

} //end of namespace bpp.

#endif //_OUTPUTASFEATURESMAFITERATOR_H_

