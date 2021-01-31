/// 786

/// This file is subject to the terms and conditions defined in
/// file 'LICENSE', which is part of this source code package.

/// Author: inumanag

/******************************************************************************/

#pragma once

/******************************************************************************/

#include <iostream>
#include <list>
#include <unordered_map>
#include <vector>

/******************************************************************************/

struct Hash {
  enum Status { HAS_UPPERCASE, ALL_LOWERCASE, HAS_N };
  uint32_t hash;
  Status status; // 2 if N, 1 if lowercase, 0 otherwise
};

struct Minimizer {
  Hash hash;
  int loc;
};

namespace std {
template <> struct hash<Hash> {
  size_t operator()(const Hash &h) const {
    return std::hash<uint32_t>()(h.hash) ^ std::hash<char>()((char)h.status);
  }
};
} // namespace std

/******************************************************************************/

struct Sequence {
  std::string name;
  std::string seq;
  bool is_rc;

  Sequence(const std::string &name, const std::string &seq, bool is_rc = false);
};

struct Index {
  const int kmer_size;
  const int window_size;
  unsigned int threshold;

  std::shared_ptr<Sequence> seq;

  // (hash, loci), sorted by loci
  std::vector<Minimizer> minimizers;
  // hash -> list of locations
  std::unordered_map<Hash, std::list<int>> index;

public:
  Index(std::shared_ptr<Sequence> seq, int kmer_size, int window_size,
        bool separate_lowercase = true);

  // Find first minimizer at loci p
  int find_minimizers(int p) const;
};

#include "extern/ostream.h"

std::ostream &operator<<(std::ostream &os, const Hash &dt);

bool operator<(const Hash &x, const Hash &y);
bool operator<=(const Hash &x, const Hash &y);
bool operator==(const Hash &x, const Hash &y);
bool operator!=(const Hash &x, const Hash &y);
bool operator==(const Minimizer &x, const Minimizer &y);
