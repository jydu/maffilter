//
// File: OutputasFeaturesMafIterator.cpp
// Authors: Julien Dutheil
// Created: Fri Mar 22 2013
//

// Copyright or © or Copr. Julien Y. Dutheil, (2013)
// SPDX-FileCopyrightText: 2026 Julien Y. Dutheil <jy.dutheil@gmail.com>
//
// SPDX-License-Identifier: GPL-3.0-or-later

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

