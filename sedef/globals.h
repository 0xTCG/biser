/// 786

/// This file is subject to the terms and conditions defined in
/// file 'LICENSE', which is part of this source code package.

/// Author: inumanag

/******************************************************************************/

#pragma once

/******************************************************************************/

#include "common.h"

/******************************************************************************/

const size_t KB = 1000;
const size_t MB = 1000 * KB;
const size_t GB = 1000 * MB;

/******************************************************************************/

namespace Globals {
struct Search {
  // search_main.cc
  static int KMER_SIZE;
  static int WINDOW_SIZE;
  static int MIN_UPPERCASE;

  // search_main.cc
  static double MAX_ERROR;
  static double MAX_EDIT_ERROR;
  static double GAP_FREQUENCY;
  static int MIN_READ_SIZE;

  // search.cc
  static const int MAX_SD_SIZE = 1 * 1024 * 1024; /// 1MB at most
};

struct Hash {
  // hash.cc
  static constexpr double INDEX_CUTOFF = 0.001;
};

struct Align { // Full SD alignment (via KSW)
  // align_main.cc + align.cc
  static int MATCH;
  static int MISMATCH;
  static int GAP_OPEN;
  static int GAP_EXTEND;

  // align.cc
  static const int MAX_KSW_SEQ_LEN = 60 * KB;
};

struct Extend { // Extension of initial SDs
  // align_main.cc
  // Formula for left and right extension:
  //		let width = max(query_len, ref_len)
  // 	extend min(MAX_EXTEND, width * EXTEND_RATIO) bases
  static double RATIO;
  static int MAX_EXTEND;
  // Max dist of two hits to be merged
  static int MERGE_DIST;
};

struct Chain {
  // chain.cc
  // Minimum number of uppercase letters in the chain
  static const int MIN_UPPERCASE_MATCH = 90;
  static const int MATCH_CHAIN_SCORE = 4;

  // Minimum chain gap (depends on user parameters)
  static int MAX_CHAIN_GAP;

  // refine.cc
  struct Refine { // SD alignment via approximate global alignment (after
                  // chaining)
    static constexpr double MATCH = 10;
    static constexpr double MISMATCH = 1;
    static constexpr double GAP = 0.5;
    static constexpr double GAPOPEN = 100; // try to approximate WGAC
    static const int MIN_READ = 900;       // Minimal refined read size
    static const int SIDE_ALIGN = 500;
    static const int MAX_GAP = 10 * KB; // Max gap during refining process
  };
};

struct Stats {
  // stats_main.cc
  // Maximal gap size (in percent of the alignment length) that should remain in
  // alignment
  static int MAX_OK_GAP;
  // Minimal size of the SDs after splitting (1KB by default)
  static int MIN_SPLIT_SIZE; // TODO <-- change to 900
  static int MIN_UPPERCASE;
  static double MAX_SCALED_ERROR;

  // Minimum number of consecutive Ns needed to call assembly gap
  static const int MIN_ASSEMBLY_GAP_SIZE = 100;
  static const int BIG_OVERLAP_THRESHOLD = 100;
};

struct Internal {
  static const bool DoUppercase = true;
  static const bool DoUppercaseSeeds = true;
  static const bool DoQgram = true;
};
}; // namespace Globals
