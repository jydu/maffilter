//
// File: PopGenMafSttistics.h
// Authors: Julien Dutheil
// Created: Tue Sep 29 2015
//

// Copyright or © or Copr. Julien Y. Dutheil, (2019)
// SPDX-FileCopyrightText: 2026 Julien Y. Dutheil <jy.dutheil@gmail.com>
//
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef _POPGENMAFSTATISTICS_H_
#define _POPGENMAFSTATISTICS_H_

#include <Bpp/Seq/Io/Maf/MafStatistics.h>

//From the STL:
#include <map>
#include <string>

namespace bpp {

/**
 * @brief Compute Fst value between two populations.
 *
 * Hudson's 92 estimator is used.
 */
class FstMafStatistics:
  public AbstractMafStatistics
{

  private:
    std::vector<std::string> pop1_, pop2_;
    unsigned int minNbPermutations_;
    unsigned int maxNbPermutations_;
    bool verbose_;

  public:
    FstMafStatistics(std::vector<std::string>& pop1, std::vector<std::string>& pop2,
        unsigned int minNbPermutations = 0,
        unsigned int maxNbPermutations = 0,
        bool verbose = true):
      AbstractMafStatistics(), pop1_(pop1), pop2_(pop2),
      minNbPermutations_(minNbPermutations),
      maxNbPermutations_(maxNbPermutations),
      verbose_(verbose)
    {}

    virtual ~FstMafStatistics() {}

  public:
    std::string getShortName() const { return "FstStatistics"; }
    std::string getFullName() const { return "Fst statistics."; }
    void compute(const MafBlock& block);
    std::vector<std::string> getSupportedTags() const;
};


} // end of namespace bpp

#endif //_POPGENMAFSTATISTICS_H_

