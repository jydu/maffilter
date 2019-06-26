//
// File: PopGenMafSttistics.cpp
// Authors: Julien Dutheil
// Created: Wed Jun 26 2019
//

/*
Copyright or © or Copr. Julien Y. Dutheil, (2019)

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

#include "PopGenMafStatistics.h"

#include <Bpp/PopGen/PolymorphismSequenceContainer.h>
#include <Bpp/PopGen/SequenceStatistics.h>

using namespace std;
using namespace bpp;

void FstMafStatistics::compute(const MafBlock& block)
{
  PolymorphismSequenceContainer poly(block.getAlignment());
  for (auto it = pop1_.begin(); it != pop1_.end(); ++it) {
    poly.setGroupId(*it, 1);
  }
  for (auto it = pop2_.begin(); it != pop2_.end(); ++it) {
    poly.setGroupId(*it, 2);
  }
  double fst = SequenceStatistics::fstHudson92(poly, 1, 2);
  result_.setValue("Fst", fst);
}


