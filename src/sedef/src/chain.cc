/// 786

/// This file is subject to the terms and conditions defined in
/// file 'LICENSE', which is part of this source code package.

/// Author: inumanag

/******************************************************************************/

#include <fstream>

#include "align.h"
#include "chain.h"
#include "common.h"
#include "fasta.h"
#include "refine.h"
#include "search.h"
#include "segment.h"

using namespace std;

/******************************************************************************/

vector<Anchor> generate_anchors(const string &query, const string &ref,
                                const Hit &orig, const int kmer_size) {
  const uint32_t MASK = (1 << (2 * kmer_size)) - 1;

  unordered_map<uint32_t, list<int>> ref_hashes;
  int last_n = -kmer_size;
  uint32_t h = 0;
  for (int i = 0; i < ref.size(); i++) {
    if (toupper(ref[i]) == 'N')
      last_n = i;
    h = ((h << 2) | hash_dna(ref[i])) & MASK;
    if (i < kmer_size - 1)
      continue;
    if (last_n >= (i - kmer_size + 1))
      continue;
    ref_hashes[h].push_back(i - kmer_size + 1);
  }

  vector<int> slide(query.size() + ref.size(), -1);
  vector<Anchor> anchors;
  dprn("got {} hashes", ref_hashes.size());

  bool same_chr = orig.query->name == orig.ref->name &&
                  orig.query->is_rc == orig.ref->is_rc;

  last_n = -kmer_size, h = 0;
  int w = 0;
  for (int i = 0; i < query.size(); i++) {
    if (toupper(query[i]) == 'N')
      last_n = i;
    h = ((h << 2) | hash_dna(query[i])) & MASK;
    if (i < kmer_size - 1)
      continue;
    if (last_n >= (i - kmer_size + 1))
      continue;

    auto it = ref_hashes.find(h);
    if (it == ref_hashes.end() || it->second.size() >= 1000)
      continue;

    int q = i - kmer_size + 1;
    int off = query.size();
    for (int r : it->second) {
      if (same_chr &&
          abs(orig.ref_start + r - (orig.query_start + q)) <= kmer_size)
        continue;
      int d = off + r - q;
      assert(d >= 0 && d < slide.size());
      if (q >= slide[d]) {
        assert(r >= d - off + slide[d]);
        bool has_u = 0;
        int len;
        for (len = 0; q + len < query.size() && r + len < ref.size(); len++) {
          if (toupper(query[q + len]) == 'N' || toupper(ref[r + len]) == 'N') {
            assert(len >= kmer_size);
            break;
          }
          if (toupper(query[q + len]) != toupper(ref[r + len])) {
            break;
          }
          has_u += bool(isupper(query[q + len]) || isupper(ref[r + len]));
        }
        if (len >= kmer_size) {
          if (anchors.size() >= (1 << 20) &&
              anchors.size() == anchors.capacity()) {
            anchors.reserve(anchors.size() * 1.5);
          }
          anchors.emplace_back(Anchor{q, r, len, has_u});
          slide[d] = q + len;
        }
      } else {
        assert(slide[d] >= q + kmer_size); // subset unless it has N!!!
        assert(d - off + slide[d] >= r + kmer_size);
      }
    }
  }
  return anchors;
}

auto chain_anchors(vector<Anchor> &anchors) {
  auto T = cur_time();

  struct Coor {
    pair<int, int> x;
    int score, pos;
    bool operator<(const Coor &a) const { return x < a.x; }
  };
  vector<Coor> xs, ys; // QUERY; pos in ANCHORS
  xs.reserve(2 * anchors.size());
  ys.reserve(anchors.size());

  int max_q = 0, max_r = 0, l = 0;
  for (int i = 0; i < anchors.size(); i++) {
    l++;
    auto &a = anchors[i];
    xs.push_back({{a.q, i}, SegmentTree<Coor>::MIN, i});
    xs.push_back({{a.q + a.l, i}, SegmentTree<Coor>::MIN, i});
    ys.push_back({{a.r + a.l - 1, i}, SegmentTree<Coor>::MIN, i});

    assert(a.l);
    max_q = max(max_q, a.q + a.l);
    max_r = max(max_r, a.r + a.l);
  }
  dprn("-- anchors to dp: {:n}", l);

  sort(xs.begin(), xs.end());
  SegmentTree<Coor> tree(ys);

  vector<int> prev(anchors.size(), -1);
  vector<pair<int, int>> dp(anchors.size());
  for (int i = 0; i < dp.size(); i++) {
    dp[i] = {0, i};
  }
  int deactivate_bound = 0;
  for (auto &x : xs) {
    int i = x.x.second;
    auto &a = anchors[i];
    if (x.x.first == a.q) {
      while (deactivate_bound < (&x - &xs[0])) {
        int t = xs[deactivate_bound].x.second; // index
        if (xs[deactivate_bound].x.first ==
            anchors[t].q + anchors[t].l) { // end point
          if (a.q - (anchors[t].q + anchors[t].l) <=
              Globals::Chain::MAX_CHAIN_GAP)
            break;
          tree.deactivate({anchors[t].r + anchors[t].l - 1, t});
        }
        deactivate_bound++;
      }

      assert(a.has_u <= a.l);
      int w = Globals::Chain::MATCH_CHAIN_SCORE * a.has_u +
              (Globals::Chain::MATCH_CHAIN_SCORE / 2) * (a.l - a.has_u);
      int j = tree.rmq({a.r - Globals::Chain::MAX_CHAIN_GAP, 0},
                       {a.r - 1, anchors.size()});
      if (j != -1 && ys[j].score != SegmentTree<Coor>::MIN) {
        j = ys[j].pos;
        auto &p = anchors[j];
        assert(a.q >= p.q + p.l);
        assert(a.r >= p.r + p.l);
        int gap = (a.q - (p.q + p.l) + a.r - (p.r + p.l));
        if (w + dp[j].first - gap > 0) {
          dp[i].first = w + dp[j].first - gap;
          prev[i] = j;
        } else {
          dp[i].first = w;
        }
      } else {
        dp[i].first = w;
      }
    } else {
      int gap = (max_q + 1 - (a.q + a.l) + max_r + 1 - (a.r + a.l));
      tree.activate({a.r + a.l - 1, i}, dp[i].first - gap);
    }
  }
  sort(dp.begin(), dp.end(), greater<pair<int, int>>());

  vector<int> path;
  path.reserve(anchors.size());
  vector<pair<int, bool>> boundaries{{0, 0}};
  vector<char> used(anchors.size(), 0);
  for (auto &m : dp) {
    int maxi = m.second;
    if (used[maxi])
      continue;
    int has_u = 0;
    while (maxi != -1 && !used[maxi]) {
      path.push_back(maxi);
      has_u += anchors[maxi].has_u;
      used[maxi] = true;
      maxi = prev[maxi];
    }
    boundaries.push_back({path.size(), has_u});
  }
  return make_pair(path, boundaries);
}

/******************************************************************************/

vector<Hit> fast_align(const string &query, const string &ref, const Hit &orig,
                       int kmer_size) {
  auto T = cur_time();
  dprn("-- aligning query {:n} --> ref {:n}", query.size(), ref.size());
  auto query_ptr = make_shared<Sequence>("QRY", query);
  auto ref_ptr = make_shared<Sequence>("REF", ref);

  /// 1. Generate the list of hits (small anchors) inside the dot graph
  auto anchors = generate_anchors(query, ref, orig, kmer_size);
  dprn("-- got {} anchors in {} s", anchors.size(), elapsed(T));
  T = cur_time();

  /// 2. Run DP on the anchors and collect all different anchors
  vector<Hit> hits;
  vector<vector<int>> guides;
  auto chains_init = chain_anchors(anchors);
  auto &bounds = chains_init.second;
  auto &chain = chains_init.first;
  for (int bi = 1; bi < bounds.size(); bi++) {
    bool has_u = bounds[bi].second;
    int be = bounds[bi].first;
    int bs = bounds[bi - 1].first;
    int up = bounds[bi].second;

    int qlo = anchors[chain[be - 1]].q,
        qhi = anchors[chain[bs]].q + anchors[chain[bs]].l;
    int rlo = anchors[chain[be - 1]].r,
        rhi = anchors[chain[bs]].r + anchors[chain[bs]].l;

    // check error
    int span = max(rhi - rlo, qhi - qlo);
    if ((!has_u || span < Globals::Chain::MIN_UPPERCASE_MATCH) &&
        span <
            Globals::Search::MIN_READ_SIZE * (1 - Globals::Search::MAX_ERROR)) {
      continue;
    }

    assert(qhi <= query.size());
    assert(rhi <= ref.size());

    Hit a{query_ptr, qlo, qhi, ref_ptr, rlo, rhi, up};
    guides.push_back(vector<int>());
    for (int bi = be - 1; bi >= bs; bi--) {
      guides.back().emplace_back(chain[bi]);
    }
    hits.push_back(a);
  }
  dprn(":: elapsed/dp = {}s", elapsed(T));
  T = cur_time();

  /// 3. Perform the full alignment
  vector<Hit> new_hits;
  for (auto &h : hits) {
    h.aln = Alignment(query, ref, anchors, guides[&h - &hits[0]]);
    update_from_alignment(h);
  }
  dprn(":: elapsed/alignment = {}s", elapsed(T));
  T = cur_time();

  /// 3. Refine these chains
  refine_chains(hits, query, ref, orig);
  dprn(":: elapsed/refinement = {}s", elapsed(T));
  T = cur_time();

  return hits;
}
