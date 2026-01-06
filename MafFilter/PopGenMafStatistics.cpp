//
// File: PopGenMafSttistics.cpp
// Authors: Julien Dutheil
// Created: Wed Jun 26 2019
//

// Copyright or © or Copr. Julien Y. Dutheil, (2019)
// SPDX-FileCopyrightText: 2026 Julien Y. Dutheil <jy.dutheil@gmail.com>
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include "PopGenMafStatistics.h"

#include <Bpp/PopGen/PolymorphismSequenceContainer.h>
#include <Bpp/PopGen/SequenceStatistics.h>
#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Numeric/Random/RandomTools.h>

using namespace std;
using namespace bpp;

vector<string> FstMafStatistics::FstMafStatistics::getSupportedTags() const
{
  vector<string> tags;
  tags.push_back("Fst");
  if (minNbPermutations_ > 0 || maxNbPermutations_ > 0) {
    tags.push_back("Fst.p-value");
    tags.push_back("Fst.nb-perm");
  }
  return tags;
}

void FstMafStatistics::compute(const MafBlock& block)
{
  auto aln = block.getAlignment();
  PolymorphismSequenceContainer poly(*aln);
  vector<string> names = block.getSpeciesList();
  poly.setSequenceNames(names, true); //Have to be unique, no paralog allowed, otherwise exception thrown.
  for (auto it = pop1_.begin(); it != pop1_.end(); ++it) {
    poly.setGroupId(*it, 1);
  }
  for (auto it = pop2_.begin(); it != pop2_.end(); ++it) {
    poly.setGroupId(*it, 2);
  }
  double fst = SequenceStatistics::fstHudson92(poly, 1, 2);
  result_.setValue("Fst", fst);
  if (minNbPermutations_ > 0 || maxNbPermutations_ > 0) {
    double nbTests = 0;
    double nbPermutations = 0;
    while(nbPermutations < minNbPermutations_ || (nbTests == 0 && nbPermutations < maxNbPermutations_)) {
      ++nbPermutations;
      if (verbose_) {
        ApplicationTools::displayGauge(static_cast<size_t>(nbPermutations), static_cast<size_t>(maxNbPermutations_), '=', "Compute Fst on permutations");
      }
      vector<string> individuals = pop1_;
      individuals.insert(individuals.end(), pop2_.begin(), pop2_.end());
      shuffle(individuals.begin(), individuals.end(), RandomTools::DEFAULT_GENERATOR);
      for (size_t j = 0; j < pop1_.size(); ++j) {
        poly.setGroupId(individuals[j], 1);
      }
      for (size_t j = pop1_.size(); j < pop1_.size() + pop2_.size(); ++j) {
        poly.setGroupId(individuals[j], 2);
      }
      double randFst = SequenceStatistics::fstHudson92(poly, 1, 2);
      if (randFst >= fst) ++nbTests;
    }
    double pv = (nbTests + 1) / static_cast<double>(nbPermutations);
    result_.setValue("Fst.p-value", pv);
    result_.setValue("Fst.nb-perm", nbPermutations);
  }
}


