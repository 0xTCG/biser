/// 786

/// This file is subject to the terms and conditions defined in
/// file 'LICENSE', which is part of this source code package.

/// Author: inumanag
/// Fast winnowing algorithm:
//  https://people.cs.uct.ac.za/~ksmith/articles/sliding_window_minimum.html

/******************************************************************************/

#include <algorithm>
#include <deque>
#include <string>
#include <vector>

#include "common.h"
#include "hash.h"

using namespace std;

/******************************************************************************/

ostream &operator<<(ostream &os, const Hash &dt) {
  os << int(dt.status) << '.' << std::hex << dt.hash;
  return os;
}

bool operator<(const Hash &x, const Hash &y) {
  return tie(x.status, x.hash) < tie(y.status, y.hash);
}

bool operator==(const Hash &x, const Hash &y) {
  return tie(x.status, x.hash) == tie(y.status, y.hash);
}

bool operator!=(const Hash &x, const Hash &y) {
  return tie(x.status, x.hash) != tie(y.status, y.hash);
}

bool operator<=(const Hash &x, const Hash &y) { return (x < y) || (x == y); }

bool operator==(const Minimizer &x, const Minimizer &y) {
  return tie(x.loc, x.hash) == tie(y.loc, y.hash);
}

bool operator<(const Minimizer &x, const Minimizer &y) {
  return tie(x.loc, x.hash) < tie(y.loc, y.hash);
}

/******************************************************************************/

vector<Minimizer> get_minimizers(const string &s, int kmer_size,
                                 const int window_size,
                                 bool separate_lowercase) {
  vector<Minimizer> minimizers;
  minimizers.reserve((2 * s.size()) / window_size);
  deque<Minimizer> window;

  const uint32_t MASK = (1 << (2 * kmer_size)) - 1;
  uint32_t h = 0;
  int last_n = -kmer_size - window_size;
  int last_u = last_n;
  for (int i = 0; i < s.size(); i++) {
    if (toupper(s[i]) == 'N') {
      last_n = i;
    } else if (isupper(s[i])) {
      last_u = i;
    }

    h = ((h << 2) | hash_dna(s[i])) & MASK;

    if (i < kmer_size - 1)
      continue;

    Hash hh{h, last_n >= (i - kmer_size + 1)
                   ? Hash::Status::HAS_N
                   : (last_u >= (i - kmer_size + 1))
                         ? Hash::Status::HAS_UPPERCASE
                         : Hash::Status::ALL_LOWERCASE};
    if (!separate_lowercase && hh.status == Hash::Status::ALL_LOWERCASE) {
      hh.status = Hash::Status::HAS_UPPERCASE;
    }
    while (!window.empty() && !(window.back().hash < hh)) {
      window.pop_back();
    }
    while (!window.empty() &&
           window.back().loc < (i - kmer_size + 1) - window_size) {
      window.pop_front();
    }
    window.push_back({hh, i - kmer_size + 1});

    if (i - kmer_size + 1 < window_size)
      continue;
    if (!minimizers.size() || !(window.front() == minimizers.back())) {
      minimizers.push_back(window.front());
    }
  }
  return minimizers;
}

/******************************************************************************/

Sequence::Sequence(const string &name, const string &seq, bool is_rc)
    : name(name), seq(seq), is_rc(is_rc) {
  if (is_rc) {
    this->seq = rc(seq);
  }
}

/******************************************************************************/

Index::Index(shared_ptr<Sequence> seq, int kmer_size, int window_size,
             bool separate_lowercase)
    : seq(seq), kmer_size(kmer_size), window_size(window_size) {
  assert(kmer_size <= 16);
  minimizers =
      get_minimizers(seq->seq, kmer_size, window_size, separate_lowercase);

  for (auto &i : minimizers) {
    index[i.hash].push_back(i.loc);
  }

  int ignore = (minimizers.size() * Globals::Hash::INDEX_CUTOFF) / 100.0;

  map<int, int> hist;
  for (auto &i : index) {
    hist[i.second.size()] += 1;
  }
  int sum = 0;
  threshold = 1 << 31;
  int j = 0;
  for (auto i = hist.rbegin(); i != hist.rend(); i++, j++) {
    sum += i->second;
    if (sum <= ignore) {
      threshold = i->first;
    } else {
      break;
    }
  }
}

int Index::find_minimizers(int p) const {
  int lo = 0, hi = minimizers.size() - 1, mid;
  while (lo <= hi) {
    mid = lo + (hi - lo) / 2;
    if (minimizers[mid].loc >= p && (!mid || minimizers[mid - 1].loc < p))
      break;
    if (minimizers[mid].loc < p) {
      lo = mid + 1;
    } else {
      hi = mid;
    }
  }
  assert(minimizers[mid].loc >= p || mid == minimizers.size() - 1);
  assert(!mid || minimizers[mid - 1].loc < p);
  if (minimizers[mid].loc < p) {
    mid++; // the last one--- no solution
  }
  return mid;
}
