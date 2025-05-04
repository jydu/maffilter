//
// File: SystemCallMafIterator.h
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

#ifndef _SYSTEMCALLMAFITERATOR_H_
#define _SYSTEMCALLMAFITERATOR_H_

#include <Bpp/Seq/Io/Maf/AbstractMafIterator.h>
#include <Bpp/Seq/Io/ISequence.h>
#include <Bpp/Seq/Io/OSequence.h>

//From the STL:
#include <iostream>
#include <string>
#include <memory>

namespace bpp {

/**
 * @brief This iterator calls an external program on each block.
 */
class SystemCallMafIterator:
  public AbstractFilterMafIterator
{
  private:
    std::unique_ptr<OAlignment> alnWriter_;
    std::string inputFile_;
    std::unique_ptr<IAlignment> alnReader_;
    std::string outputFile_;
    std::string call_;
    bool hotTest_;
    std::string refSeq_; //Sequence where to store alignment scores. If empty, do no store scores.

  public:
    SystemCallMafIterator(
        std::shared_ptr<MafIteratorInterface> iterator,
        std::unique_ptr<OAlignment> alnWriter,
        const std::string& inputFile,
        std::unique_ptr<IAlignment> alnReader,
        const std::string& outputFile,
        const std::string& callCmd,
        bool hotTest = false,
	const std::string refSeq = "") :
      AbstractFilterMafIterator(iterator),
          alnWriter_(std::move(alnWriter)),
          inputFile_(inputFile),
          alnReader_(std::move(alnReader)),
          outputFile_(outputFile),
          call_(callCmd),
          hotTest_(hotTest),
	  refSeq_(refSeq)
    {}

  private:
    SystemCallMafIterator(const SystemCallMafIterator& iterator) = delete;
          
    SystemCallMafIterator& operator=(const SystemCallMafIterator& iterator) = delete;
    
  public:
    std::unique_ptr<MafBlock> analyseCurrentBlock_();

};

} // namespace bpp.

#endif //_SYSTEMCALLMAFITERATOR_H_

