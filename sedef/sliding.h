/// 786

/// This file is subject to the terms and conditions defined in
/// file 'LICENSE', which is part of this source code package.

/// Author: inumanag

/******************************************************************************/

#pragma once

/******************************************************************************/

#include <map>
#include <queue>
#include <set>
#include <string>
#include <vector>

#include "hash.h"

/******************************************************************************/

struct SlidingMap {
  std::map<Hash, char> storage;
  typename std::map<Hash, char>::iterator boundary; // this is inclusive!
  int query_size;
  int intersection;
  double limit;
  int kmer_size;

  std::unordered_map<int, int> estimate_memoize;

private:
  SlidingMap() = default;

public:
  SlidingMap(int kmer_size);
  SlidingMap(const SlidingMap &other);
  SlidingMap(SlidingMap &&other);
  SlidingMap &operator=(SlidingMap other);

public:
  int jaccard();

  // static SlidingMap fromMap(const SlidingMap &m);
  std::string print_it(const std::map<Hash, char>::iterator &boundary) const;

  // when adding, if same hash is found: try to match earliest one if added by
  // another set (i.e. try increase jaccard) when removing, if same hash is
  // found: try to remove the latest one (try to preserve jaccard)

  bool add(const Hash &h, int BIT, int FULL = 3);
  bool remove(const Hash &h, int BIT, int FULL = 3);

  void add_to_query(const Hash &h);
  void remove_from_query(const Hash &h);
  void add_to_reference(const Hash &h);
  void remove_from_reference(const Hash &h);

  friend void swap(SlidingMap &first, SlidingMap &second) // nothrow
  {
    using std::swap;

    swap(first.storage, second.storage);
    swap(first.boundary, second.boundary);
    swap(first.query_size, second.query_size);
    swap(first.intersection, second.intersection);
    swap(first.limit, second.limit);
    swap(first.kmer_size, second.kmer_size);
  }
};
