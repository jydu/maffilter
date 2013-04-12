//
// File: OutputasFeaturesMafIterator.cpp
// Authors: Julien Dutheil
// Created: Fri Mar 22 2013
//

/*
Copyright or © or Copr. Bio++ Development Team, (2010)

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

#include "OutputAsFeaturesMafIterator.h"

using namespace bpp;

//From the STL:
#include <string>
#include <vector>

using namespace std;

void OutputAsFeaturesMafIterator::writeHeader(std::ostream& out) const
{
  out << "##Maf2Annotations" << endl;
  out << "Chr\tStart\tStop\tStrand\tSrcSize" << endl;
  //We probably here ultimately want to write in GFF or sthg like this?
}

void OutputAsFeaturesMafIterator::writeBlock(std::ostream& out, const MafBlock& block) const
{
  if (block.hasSequenceForSpecies(species_)) {
    vector<const MafSequence*> sequences = block.getSequencesForSpecies(species_);
    for (size_t i = 0; i < sequences.size(); ++i) {
      out << sequences[i]->getChromosome() << "\t" << sequences[i]->start() << "\t" << sequences[i]->stop() << "\t" << sequences[i]->getStrand() << "\t" << sequences[i]->getSrcSize() << endl; 
    }
  }
}

