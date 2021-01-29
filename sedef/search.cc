/// 786

/// This file is subject to the terms and conditions defined in
/// file 'LICENSE', which is part of this source code package.

/// Author: inumanag
/// Based on
/// http://www.biorxiv.org/content/biorxiv/early/2017/03/24/103812.full.pdf

/******************************************************************************/

#include <chrono>
#include <cmath>
#include <list>
#include <queue>
#include <string>
#include <vector>

#include "chain.h"
#include "common.h"
#include "filter.h"
#include "search.h"
#include "sliding.h"

using namespace std;

/******************************************************************************/

/* extern */ int64_t TOTAL_ATTEMPTED = 0;
/* extern */ int64_t JACCARD_FAILED = 0;
/* extern */ int64_t INTERVAL_FAILED = 0;

/******************************************************************************/

bool is_overlap(Tree &tree, int pf_pos, int pf_end, int pfp_pos, int pfp_end) {
  assert(pf_pos <= pf_end);
  assert(pfp_pos <= pfp_end);

  auto pf = tree.find(pf_pos);
  if (pf == tree.end())
    return false;

  auto pfp = pf->second.find(pfp_pos);
  if (pfp == pf->second.end())
    return false;

  for (auto &it : pfp->second) {
    int sA = it.first.lower(), eA = it.first.upper();
    int sB = it.second.lower(), eB = it.second.upper();

    // 1. total overlap
    if (pf_pos >= sA && pf_end <= eA && pfp_pos >= sB && pfp_end <= eB)
      return true;

    // 2. do not check overlaps with too small intervals
    if (min(eA - sA, eB - sB) < Globals::Search::MIN_READ_SIZE * 1.5)
      continue;

    // 3. if partial overlap, it must be greater than RIGHT_ALLOWANCE (i.e. do
    // not count small overlaps as hits)
    // ----------+
    //     +------------
    //     |*****| the space inside must be at least RIGHT_ALLOWANCE to trigger
    //     overlap
    const int RIGHT_ALLOWANCE =
        Globals::Search::MIN_READ_SIZE; /// TODO more mathy formulation
    if (eA - pf_pos >= RIGHT_ALLOWANCE && eB - pfp_pos >= RIGHT_ALLOWANCE)
      return true;
  }
  return false;
}

auto parse_hits(vector<Hit> &hits) {
  vector<Hit> hits_real;
  for (auto &h : hits) { // TODO fix brute force
    bool add = true;
    for (auto &ph : hits)
      if ((&h - &hits[0]) != (&ph - &hits[0])) {
        if (h.ref_start >= ph.ref_start && h.ref_end <= ph.ref_end &&
            h.query_start >= ph.query_start &&
            h.query_end <= ph.query_end) { // if full match
          add = false;
          break;
        }
      }
    if (add) {
      hits_real.push_back(h);
    }
  }
  return hits_real;
}

/******************************************************************************/

Hit extend(SlidingMap &winnow, shared_ptr<Index> query_hash, int query_start,
           int query_end, int query_winnow_start, int query_winnow_end,
           shared_ptr<Index> ref_hash, int ref_start, int ref_end,
           int ref_winnow_start, int ref_winnow_end, bool same_genome) {
  assert(query_start < query_hash->seq->seq.size());
  assert(ref_start < ref_hash->seq->seq.size());
  assert(query_end <= query_hash->seq->seq.size());
  assert(ref_end <= ref_hash->seq->seq.size());

  // Always extend to the boundary of the next winnow
  auto do_extend_query_right = [&]() {
    if (query_winnow_end >= query_hash->minimizers.size())
      return false;
    winnow.add_to_query(query_hash->minimizers[query_winnow_end++].hash);
    query_end = query_winnow_end < query_hash->minimizers.size()
                    ? query_hash->minimizers[query_winnow_end].loc
                    : query_hash->seq->seq.size();
    return true;
  };
  auto undo_extend_query_right = [&]() {
    winnow.remove_from_query(query_hash->minimizers[--query_winnow_end].hash);
    query_end = query_hash->minimizers[query_winnow_end].loc;
  };

  auto do_extend_ref_right = [&]() {
    if (ref_winnow_end >= ref_hash->minimizers.size())
      return false;
    winnow.add_to_reference(ref_hash->minimizers[ref_winnow_end++].hash);
    ref_end = ref_winnow_end < ref_hash->minimizers.size()
                  ? ref_hash->minimizers[ref_winnow_end].loc
                  : ref_hash->seq->seq.size();
    return true;
  };
  auto undo_extend_ref_right = [&]() {
    winnow.remove_from_reference(ref_hash->minimizers[--ref_winnow_end].hash);
    ref_end = ref_hash->minimizers[ref_winnow_end].loc;
  };

  auto do_extend_both_right = [&]() {
    if (ref_winnow_end >= ref_hash->minimizers.size() ||
        query_winnow_end >= query_hash->minimizers.size())
      return false;
    bool r = do_extend_query_right();
    r &= do_extend_ref_right();
    return r;
  };
  auto undo_extend_both_right = [&]() {
    undo_extend_ref_right();
    undo_extend_query_right();
  };

  auto do_extend_query_left = [&]() {
    if (!query_winnow_start)
      return false;
    winnow.add_to_query(query_hash->minimizers[--query_winnow_start].hash);
    query_start = query_winnow_start
                      ? query_hash->minimizers[query_winnow_start - 1].loc + 1
                      : 0;
    return true;
  };
  auto undo_extend_query_left = [&]() {
    query_start = query_hash->minimizers[query_winnow_start].loc + 1;
    winnow.remove_from_query(query_hash->minimizers[query_winnow_start++].hash);
  };

  auto do_extend_ref_left = [&]() {
    if (!ref_winnow_start)
      return false;
    winnow.add_to_reference(ref_hash->minimizers[--ref_winnow_start].hash);
    ref_start = ref_winnow_start
                    ? ref_hash->minimizers[ref_winnow_start - 1].loc + 1
                    : 0;
    return true;
  };
  auto undo_extend_ref_left = [&]() {
    ref_start = ref_hash->minimizers[ref_winnow_start].loc + 1;
    winnow.remove_from_reference(ref_hash->minimizers[ref_winnow_start++].hash);
  };

  auto do_extend_both_left = [&]() {
    if (!query_winnow_start || !ref_winnow_start)
      return false;
    bool r = do_extend_query_left();
    r &= do_extend_ref_left();
    return r;
  };
  auto undo_extend_both_left = [&]() {
    undo_extend_ref_left();
    undo_extend_query_left();
  };

  auto do_extend_both_both = [&]() {
    if (!query_winnow_start || !ref_winnow_start)
      return false;
    if (ref_winnow_end >= ref_hash->minimizers.size() ||
        query_winnow_end >= query_hash->minimizers.size())
      return false;
    bool r = do_extend_both_left();
    r &= do_extend_both_right();
    return r;
  };
  auto undo_extend_both_both = [&]() {
    undo_extend_both_right();
    undo_extend_both_left();
  };

#define p(q) make_pair(do_extend_##q, undo_extend_##q)
  auto extensions = vector<pair<function<bool(void)>, function<void(void)>>>{
      p(both_both), p(both_right), p(both_left)};
#undef p

  // First extend to the boundaries
  query_start = query_winnow_start
                    ? query_hash->minimizers[query_winnow_start - 1].loc + 1
                    : 0;
  query_end = query_winnow_end < query_hash->minimizers.size()
                  ? query_hash->minimizers[query_winnow_end].loc
                  : query_hash->seq->seq.size();
  ref_start =
      ref_winnow_start ? ref_hash->minimizers[ref_winnow_start - 1].loc + 1 : 0;
  ref_end = ref_winnow_end < ref_hash->minimizers.size()
                ? ref_hash->minimizers[ref_winnow_end].loc
                : ref_hash->seq->seq.size();

  const double MAX_GAP_ERROR =
      Globals::Search::MAX_ERROR - Globals::Search::MAX_EDIT_ERROR;
  for (;;) {
    int max_match =
        min(Globals::Search::MAX_SD_SIZE,
            same_genome
                ? int((1.0 / MAX_GAP_ERROR + .5) * abs(query_start - ref_start))
                : Globals::Search::MAX_SD_SIZE);
    int aln_len = max(query_end - query_start, ref_end - ref_start);
    int seq_len = min(query_end - query_start, ref_end - ref_start);
    if (aln_len > max_match ||
        pct(seq_len, aln_len) < 100 * (1 - 2 * MAX_GAP_ERROR)) {
      break;
    }

    if (same_genome) {
      int overlap = query_end - ref_start;
      if (overlap > 0 &&
          pct(overlap, ref_end - ref_start) > 100 * Globals::Search::MAX_ERROR)
        break;
    }

    bool extended = false;
    for (auto &fn : extensions) {
      if (!fn.first())
        continue;
      if (winnow.jaccard() >= 0) {
        extended = true;
        break;
      } else {
        fn.second();
      }
    }
    if (!extended)
      break;
  }

  return Hit{
      query_hash->seq, query_start,      query_end, ref_hash->seq, ref_start,
      ref_end,         winnow.jaccard(), "",        "OK",          {}};
}

/******************************************************************************/

vector<Hit> search_in_reference_interval(
    int query_start, int query_winnow_start, int query_winnow_end,
    shared_ptr<Index> query_hash, shared_ptr<Index> ref_hash, Tree &tree,
    bool same_genome, int init_len, bool allow_extend, bool report_fails,
    SlidingMap winnow, int t_start, int t_end) {
  assert(t_start <= t_end);
  assert(t_start >= 0);
  assert(winnow.query_size > 0);

  TOTAL_ATTEMPTED++;

  int ref_start = t_start,
      ref_end = min(t_start + init_len, (int)ref_hash->seq->seq.size());
  int ref_winnow_start = ref_hash->find_minimizers(ref_start);
  assert(ref_winnow_start < ref_hash->minimizers.size());

  int ref_winnow_end =
      ref_winnow_start; // winnow is W(query) ; extend it to W(query) | W(ref)
  for (; ref_winnow_end < ref_hash->minimizers.size() &&
         ref_hash->minimizers[ref_winnow_end].loc < ref_end;
       ref_winnow_end++) {
    winnow.add_to_reference(ref_hash->minimizers[ref_winnow_end].hash);
  }

  // Roll until we find best inital match
  // TODO: optimize
  SlidingMap best_winnow(winnow);
  int best_ref_start = ref_start, best_ref_end = ref_end;
  int best_ref_winnow_start = ref_winnow_start,
      best_ref_winnow_end = ref_winnow_end;
  while (ref_start < t_end && ref_end < ref_hash->seq->seq.size()) {
    if (ref_winnow_start < ref_hash->minimizers.size() &&
        ref_hash->minimizers[ref_winnow_start].loc < ref_start + 1) {
      winnow.remove_from_reference(
          ref_hash->minimizers[ref_winnow_start++].hash);
    }
    if (ref_winnow_end < ref_hash->minimizers.size() &&
        ref_hash->minimizers[ref_winnow_end].loc == ref_end) {
      winnow.add_to_reference(ref_hash->minimizers[ref_winnow_end++].hash);
    }
    if (winnow.jaccard() > best_winnow.jaccard()) {
      best_ref_start = ref_start;
      best_ref_end = ref_end;
      best_ref_winnow_start = ref_winnow_start;
      best_ref_winnow_end = ref_winnow_end;
      best_winnow = winnow;
    }
    ref_start++;
    ref_end++;
    if (ref_end == ref_hash->seq->seq.size())
      break;
  }
  // END TODO

  vector<Hit> hits;

  if (best_winnow.jaccard() < 0) {
    JACCARD_FAILED++;
    if (report_fails)
      hits.push_back({query_hash->seq,
                      query_start,
                      query_start + init_len,
                      ref_hash->seq,
                      best_ref_start,
                      best_ref_end,
                      best_winnow.jaccard(),
                      "",
                      fmt::format("jaccard: {} < {}",
                                  best_winnow.limit + best_winnow.jaccard(),
                                  best_winnow.limit),
                      {}});
  } else if (allow_extend) {
    if (!is_overlap(tree, query_start, query_start + init_len, best_ref_start,
                    best_ref_end)) {
      auto f = filter(query_hash->seq->seq, query_start, query_start + init_len,
                      ref_hash->seq->seq, ref_start, ref_end);
      if (!f.first) {
        if (report_fails)
          hits.push_back({query_hash->seq,
                          query_start,
                          query_start + init_len,
                          ref_hash->seq,
                          ref_start,
                          ref_end,
                          0,
                          "",
                          f.second,
                          {}});
      } else {
        Hit h = extend(best_winnow, query_hash, query_start,
                       query_start + init_len, query_winnow_start,
                       query_winnow_end, ref_hash, best_ref_start, best_ref_end,
                       best_ref_winnow_start, best_ref_winnow_end, same_genome);
        f = filter(query_hash->seq->seq, h.query_start, h.query_end,
                   ref_hash->seq->seq, h.ref_start, h.ref_end);
        if (!f.first) {
          if (report_fails) {
            h.comment = f.second;
            hits.push_back(h);
          }
        } else {
          hits.push_back(h);
          auto a = Interval(h.query_start, h.query_end);
          auto b = Interval(h.ref_start, h.ref_end);
          tree += make_pair(a, Subtree({b, {make_pair(a, b)}}));
        }
      }
    } else {
      INTERVAL_FAILED++;
    }
  } else {
    auto f = filter(query_hash->seq->seq, query_start, query_start + init_len,
                    ref_hash->seq->seq, best_ref_start, best_ref_end);
    if (f.first || report_fails) {
      hits.push_back({query_hash->seq,
                      query_start,
                      query_start + init_len,
                      ref_hash->seq,
                      best_ref_start,
                      best_ref_end,
                      best_winnow.jaccard(),
                      "",
                      f.second == "" ? "OK_INIT" : f.second,
                      {}});
    }
  }

  return hits;
}

/******************************************************************************/

vector<Hit> search(int query_winnow_start, shared_ptr<Index> query_hash,
                   shared_ptr<Index> ref_hash, Tree &tree,
                   const bool same_genome, const int init_len,
                   const bool allow_extend, const bool report_fails) {
  if (query_winnow_start >= query_hash->minimizers.size())
    return {};

  int query_start = query_hash->minimizers[query_winnow_start].loc;
  if (query_start + init_len > query_hash->seq->seq.size())
    return {};

  assert(query_hash->kmer_size == ref_hash->kmer_size);
  SlidingMap init_winnow(query_hash->kmer_size);
  set<int> candidates_prel;
  int query_winnow_end = query_winnow_start;
  for (; query_winnow_end < query_hash->minimizers.size() &&
         query_hash->minimizers[query_winnow_end].loc - query_start <= init_len;
       query_winnow_end++) {
    auto &h = query_hash->minimizers[query_winnow_end].hash;
    init_winnow.add_to_query(h);

    if (Globals::Internal::DoUppercaseSeeds &&
        h.status != Hash::Status::HAS_UPPERCASE) // use only hashes with
                                                 // uppercase character!
      continue;
    auto pf = tree.find(query_hash->minimizers[query_winnow_end].loc);
    auto ptr = ref_hash->index.find(h);
    if (ptr == ref_hash->index.end() ||
        ptr->second.size() >= ref_hash->threshold) {
      continue;
    } else
      for (auto pos : ptr->second) {
        if (!same_genome ||
            pos >=
                query_start + init_len) { // Make sure to have at least read_len
                                          // spacing if reference = query
          if (pf == tree.end() || pf->second.find(pos) == pf->second.end()) {
            candidates_prel.insert(pos);
          }
        }
      }
  }
  if (!init_winnow.query_size)
    return {};

  vector<pair<int, int>> T;
  vector<int> candidates(candidates_prel.begin(), candidates_prel.end());
  for (int i = 0; i <= (int)candidates.size() - init_winnow.limit; i++) {
    int j = i + (init_winnow.limit - 1);
    if (candidates[j] - candidates[i] <= init_len) {
      int x = max(0, candidates[j] - init_len + 1), y = candidates[i] + 1;
      if (T.size() && x < T.back().second) {
        T.back().second = max(T.back().second, y);
      } else {
        T.push_back({x, y});
      }
    }
  }

  vector<Hit> hits;
  for (auto &t : T) {
    if (same_genome) {
      t.first = max(t.first, query_start + init_len);
    }
    if (t.first > t.second)
      continue;
    auto h = search_in_reference_interval(
        query_start, query_winnow_start, query_winnow_end, query_hash, ref_hash,
        tree, same_genome, init_len, allow_extend, report_fails, init_winnow,
        t.first, t.second);
    for (auto &hh : h)
      hits.push_back(hh);
  }

  tree -= Interval(0, query_start - Globals::Search::MIN_READ_SIZE);
  return parse_hits(hits);
}
