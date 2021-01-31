/// 786

/// This file is subject to the terms and conditions defined in
/// file 'LICENSE', which is part of this source code package.

/// Author: inumanag

/******************************************************************************/

#pragma once

/******************************************************************************/

#include <deque>
#include <string>

#include "align.h"
#include "common.h"
#include "fasta.h"

/******************************************************************************/

struct Hit {
  std::shared_ptr<Sequence> query;
  int query_start, query_end; // query range

  std::shared_ptr<Sequence> ref;
  int ref_start, ref_end; // reference range

  int jaccard; // coordinates of seed matches
  std::string name, comment;

  Alignment aln;

public:
  static Hit from_bed(const std::string &bed, std::string *cigar = nullptr);
  static Hit from_bed(const std::string &bed, std::shared_ptr<Sequence> query,
                      std::shared_ptr<Sequence> ref);
  static Hit from_wgac(const std::string &bed);

public:
  std::string to_bed(bool do_rc = true, bool with_cigar = true,
                     const FastaReference *fr = nullptr) const;
  bool operator<(const Hit &h) const {
    return std::tie(query_start, query_end, ref_start, ref_end) <
           std::tie(h.query_start, h.query_end, h.ref_start, h.ref_end);
  }

public:
  void extend(const double factor, const int max_extend);
};

/******************************************************************************/

void update_from_alignment(Hit &h);
