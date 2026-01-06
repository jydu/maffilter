//
// File: TreeBuildingSystemCallMafIterator.h
// Authors: Julien Dutheil
// Created: Sat Jun 18 2016
//

// Copyright or © or Copr. Julien Y. Dutheil, (2016)
// SPDX-FileCopyrightText: 2026 Julien Y. Dutheil <jy.dutheil@gmail.com>
//
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef _TREEBUILDINGSYSTEMCALLMAFITERATOR_H_
#define _TREEBUILDINGSYSTEMCALLMAFITERATOR_H_

#include <Bpp/Seq/Io/Maf/AbstractMafIterator.h>
#include <Bpp/Phyl/Io/IoTree.h>
#include <Bpp/Seq/Io/OSequence.h>

//From the STL:
#include <iostream>
#include <string>
#include <memory>

namespace bpp {

/**
 * @brief This iterator calls an external program on each block.
 */
class TreeBuildingSystemCallMafIterator:
  public AbstractFilterMafIterator
{
  private:
    std::unique_ptr<OAlignment> alnWriter_;
    std::string inputFile_;
    std::unique_ptr<ITree> treeReader_;
    std::string outputFile_;
    std::string call_;
    std::string propertyName_;

  public:
    TreeBuildingSystemCallMafIterator(
        std::shared_ptr<MafIteratorInterface> iterator,
        std::unique_ptr<OAlignment> alnWriter,
        const std::string& inputFile,
        std::unique_ptr<ITree> treeReader,
        const std::string& outputFile,
        const std::string& callCmd,
        const std::string& propertyName) :
      AbstractFilterMafIterator(iterator),
          alnWriter_(std::move(alnWriter)),
          inputFile_(inputFile),
          treeReader_(std::move(treeReader)),
          outputFile_(outputFile),
          call_(callCmd),
          propertyName_(propertyName)
    {}

  private:
    TreeBuildingSystemCallMafIterator(const TreeBuildingSystemCallMafIterator& iterator) :
      AbstractFilterMafIterator(0),
      alnWriter_(),
      inputFile_(iterator.inputFile_),
      treeReader_(),
      outputFile_(iterator.outputFile_),
      call_(iterator.call_),
      propertyName_(iterator.propertyName_)
    {}
    
    TreeBuildingSystemCallMafIterator& operator=(const TreeBuildingSystemCallMafIterator& iterator)
    {
      inputFile_ = iterator.inputFile_;
      outputFile_ = iterator.outputFile_;
      call_ = iterator.call_;
      propertyName_ = iterator.propertyName_;
      return *this;
    }


  public:
    std::unique_ptr<MafBlock> analyseCurrentBlock_() override;

};

} // namespace bpp.

#endif //_TREEBUILDINGSYSTEMCALLMAFITERATOR_H_

