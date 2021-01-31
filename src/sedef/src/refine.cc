/// 786

/// This file is subject to the terms and conditions defined in
/// file 'LICENSE', which is part of this source code package.

/// Author: inumanag

/******************************************************************************/

#include "refine.h"
#include "align.h"
#include "chain.h"
#include "common.h"
#include "fasta.h"
#include "hit.h"
#include "search.h"
#include "segment.h"

using namespace std;

/******************************************************************************/

void refine_chains(vector<Hit> &anchors, const string &qseq, const string &rseq,
                   const Hit &orig) {
  dprn(":: taking {} anchors for refinement", anchors.size());

  sort(anchors.begin(), anchors.end());

  bool same_chr = orig.query->name == orig.ref->name &&
                  orig.query->is_rc == orig.ref->is_rc;
  vector<int> score;
  for (auto &a : anchors) {
    score.push_back(+Globals::Chain::Refine::MATCH * a.aln.matches() -
                    Globals::Chain::Refine::MISMATCH * a.aln.mismatches() -
                    Globals::Chain::Refine::GAP * a.aln.gap_bases());
  }

  vector<int> dp(anchors.size(), 0);
  vector<int> prev(anchors.size(), -1);
  set<pair<int, int>, greater<pair<int, int>>> maxes;
  for (int ai = 0; ai < anchors.size(); ai++) {
    if (same_chr) {
      auto &c = anchors[ai];
      int qlo = c.query_start, qhi = c.query_end;
      int rlo = c.ref_start, rhi = c.ref_end;
      int qo = max(0, min(orig.query_start + qhi, orig.ref_start + rhi) -
                          max(orig.query_start + qlo, orig.ref_start + rlo));
      if ((rhi - rlo) - qo < Globals::Chain::Refine::SIDE_ALIGN &&
          (qhi - qlo) - qo <
              Globals::Chain::Refine::SIDE_ALIGN) { // no gap between
        continue;
      }
    }

    dp[ai] = score[ai];
    for (int aj = ai - 1; aj >= 0; aj--) {
      auto &c = anchors[ai];
      auto &p = anchors[aj];

      int cqs = c.query_start;
      if (cqs < p.query_end) {
        cqs = p.query_end;
      }
      int crs = c.ref_start;
      if (crs < p.ref_end) {
        crs = p.ref_end;
      }

      if (p.query_end >= c.query_end || p.ref_end >= c.ref_end)
        continue;
      if (p.ref_start >= c.ref_start)
        continue;

      int ma = max(cqs - p.query_end, crs - p.ref_end);
      int mi = min(cqs - p.query_end, crs - p.ref_end);

      if (ma >= Globals::Chain::Refine::MAX_GAP)
        continue;

      if (same_chr) {
        int qlo = p.query_end, qhi = cqs;
        int rlo = p.ref_end, rhi = crs;
        int qo = max(0, min(orig.query_start + qhi, orig.ref_start + rhi) -
                            max(orig.query_start + qlo, orig.ref_start + rlo));
        if (qo >= 1) { // no gap between
          continue;
        }
      }

      int mis = Globals::Chain::Refine::MISMATCH * mi,
          gap = Globals::Chain::Refine::GAPOPEN +
                Globals::Chain::Refine::GAP * (ma - mi);
      int sco = dp[aj] + score[ai] - mis - gap;
      if (sco >= dp[ai]) {
        dp[ai] = sco;
        prev[ai] = aj;
      }
    }
    maxes.insert({dp[ai], ai});
  }

  vector<bool> used(anchors.size(), 0);
  vector<deque<int>> paths;
  vector<Hit> hits;
  for (auto &m : maxes) {
    if (m.first == 0)
      break;
    int maxi = m.second;
    if (used[maxi])
      continue;
    paths.push_back(deque<int>());
    int hasu = 0;
    while (maxi != -1 && (!used[maxi])) {
      paths.back().push_front(maxi);
      hasu += anchors[maxi].jaccard;
      used[maxi] = true;
      maxi = prev[maxi];
    }

    int qlo = anchors[paths.back().front()].query_start,
        qhi = anchors[paths.back().back()].query_end;
    int rlo = anchors[paths.back().front()].ref_start,
        rhi = anchors[paths.back().back()].ref_end;

    int est_size = anchors[paths.back()[0]].aln.span();
    for (int i = 1; i < paths.back().size(); i++) {
      est_size += anchors[paths.back()[i]].aln.span();
      est_size += max(anchors[paths.back()[i]].query_start -
                          anchors[paths.back()[i - 1]].query_end,
                      anchors[paths.back()[i]].ref_start -
                          anchors[paths.back()[i - 1]].ref_end);
    }
    dprn("-- chain: (est:{}, size:{}) {}..{} --> {}..{}  ## {}..{} --> {}..{} ",
         est_size, paths.back().size(), qlo, qhi, rlo, rhi,
         qlo + orig.query_start, qhi + orig.query_start, rlo + orig.ref_start,
         rhi + orig.ref_start);
    for (auto p : paths.back()) {
      auto &y = anchors[p];
      dprn("    {}..{}->{}..{}", y.query_start, y.query_end, y.ref_start,
           y.ref_end);
    }

    if (est_size <
        Globals::Chain::Refine::MIN_READ - Globals::Chain::Refine::SIDE_ALIGN) {
      dprn("est size failed");
      continue;
    }

    bool overlap = 0;
    for (auto &h : hits) {
      int qo = max(0, min(qhi, h.query_end) - max(qlo, h.query_start));
      int ro = max(0, min(rhi, h.ref_end) - max(rlo, h.ref_start));

      if (qhi - qlo - qo < Globals::Chain::Refine::SIDE_ALIGN &&
          rhi - rlo - ro < Globals::Chain::Refine::SIDE_ALIGN) {
        dprn("between overlap failed");
        overlap = 1;
        break;
      }
    }
    if (overlap)
      continue;

    auto hit =
        Hit{anchors.front().query, qlo, qhi, anchors.front().ref, rlo, rhi};

    vector<Hit> guide;
    Hit *prev = &anchors[paths.back()[0]];
    for (int pi = 1; pi < paths.back().size(); pi++) {
      auto &cur = anchors[paths.back()[pi]];
      if (cur.query_start < prev->query_end || cur.ref_start < prev->ref_end) {
        prev->aln.merge(cur.aln, qseq, rseq);
        update_from_alignment(*prev);
      } else {
        guide.push_back(*prev);
        prev = &cur;
      }
    }
    guide.push_back(*prev);

    hit.aln = Alignment(hit.query->seq, hit.ref->seq, guide,
                        Globals::Chain::Refine::SIDE_ALIGN);
    update_from_alignment(hit);
    if (hit.aln.span() >= Globals::Chain::Refine::MIN_READ) {
      dprn("IN!");
      hits.push_back(hit);
    } else {
      dprn("failed final size ");
    }
  }

  anchors = hits;
}
