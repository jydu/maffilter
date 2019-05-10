//
// File: SystemCallMafIterator.h
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

#ifndef _SYSTEMCALLMAFITERATOR_H_
#define _SYSTEMCALLMAFITERATOR_H_

#include <Bpp/Seq/Io/Maf/MafIterator.h>
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

  public:
    SystemCallMafIterator(
        MafIterator* iterator,
        OAlignment* alnWriter,
        const std::string& inputFile,
        IAlignment* alnReader,
        const std::string& outputFile,
        const std::string& callCmd,
        bool hotTest = false) :
      AbstractFilterMafIterator(iterator),
          alnWriter_(alnWriter),
          inputFile_(inputFile),
          alnReader_(alnReader),
          outputFile_(outputFile),
          call_(callCmd),
          hotTest_(hotTest)
    {}

  private:
    SystemCallMafIterator(const SystemCallMafIterator& iterator) :
      AbstractFilterMafIterator(0),
      alnWriter_(),
      inputFile_(iterator.inputFile_),
      alnReader_(),
      outputFile_(iterator.outputFile_),
      call_(iterator.call_),
      hotTest_(iterator.hotTest_)
    {}
    
    SystemCallMafIterator& operator=(const SystemCallMafIterator& iterator)
    {
      inputFile_ = iterator.inputFile_;
      outputFile_ = iterator.outputFile_;
      call_ = iterator.call_;
      hotTest_ = iterator.hotTest_;
      return *this;
    }


  public:
    MafBlock* analyseCurrentBlock_();

};

} // namespace bpp.

#endif //_SYSTEMCALLMAFITERATOR_H_

