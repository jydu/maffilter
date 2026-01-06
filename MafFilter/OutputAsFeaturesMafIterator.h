//
// File: OutputasFeaturesMafIterator.h
// Authors: Julien Dutheil
// Created: Fri Mar 22 2013
//

// Copyright or © or Copr. Julien Y. Dutheil, (2013)
// SPDX-FileCopyrightText: 2026 Julien Y. Dutheil <jy.dutheil@gmail.com>
//
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef _OUTPUTASFEATURESMAFITERATOR_H_
#define _OUTPUTASFEATURESMAFITERATOR_H_

#include <Bpp/Seq/Io/Maf/AbstractMafIterator.h>

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
    std::shared_ptr<std::ostream> output_;
    std::string species_;

  public:
    OutputAsFeaturesMafIterator(
	std::shared_ptr<MafIteratorInterface> iterator,
       	std::shared_ptr<std::ostream> out,
       	const std::string& species) :
      AbstractFilterMafIterator(iterator),
      output_(out),
      species_(species)
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
    std::unique_ptr<MafBlock> analyseCurrentBlock_() {
      currentBlock_ = iterator_->nextBlock();
      if (output_ && currentBlock_)
        writeBlock(*output_, *currentBlock_);
      return std::move(currentBlock_);
    }

  private:
    void writeHeader(std::ostream& out) const;
    void writeBlock(std::ostream& out, const MafBlock& block) const;
};

} //end of namespace bpp.

#endif //_OUTPUTASFEATURESMAFITERATOR_H_

