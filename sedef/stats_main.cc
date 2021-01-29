/// 786

/// This file is subject to the terms and conditions defined in
/// file 'LICENSE', which is part of this source code package.

/// Authors: alimg, inumanag

/******************************************************************************/

#include <algorithm>
#include <bitset>
#include <chrono>
#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <boost/dynamic_bitset.hpp>

#include "align.h"
#include "common.h"
#include "extern/argh.h"
#include "fasta.h"
#include "hit.h"
#include "merge.h"

using namespace std;

/******************************************************************************/

bool subhit(const Hit &hin, int start, int end, Hit &h) {
  dprn(">subhitting {} {} vs {} {} ", start, end, 0, hin.aln.alignment.size());
  if (end >= hin.aln.alignment.size())
    end = hin.aln.alignment.size();
  if (start >= end)
    return false;
  h = hin;
  int sa = 0, la = 0, sb = 0, lb = 0;
  for (int i = 0; i < end; i++) {
    if (h.aln.align_a[i] != '-') {
      if (i < start) {
        sa++;
      } else {
        la++;
      }
    }
    if (h.aln.align_b[i] != '-') {
      if (i < start) {
        sb++;
      } else {
        lb++;
      }
    }
  }
  h.aln.align_a = h.aln.align_a.substr(start, end - start);
  h.aln.alignment = h.aln.alignment.substr(start, end - start);
  h.aln.align_b = h.aln.align_b.substr(start, end - start);

  h.aln.a = h.aln.a.substr(sa, la);
  h.aln.start_a = 0;
  h.aln.end_a = la;

  h.aln.b = h.aln.b.substr(sb, lb);
  h.aln.start_b = 0;
  h.aln.end_b = lb;

  h.aln.cigar_from_alignment();
  h.aln.trim_back();
  h.aln.trim_front();

  h.query_start += sa;
  h.query_end = h.query_start + la;
  assert(!h.query->is_rc);
  if (h.ref->is_rc) {
    dprn("applying RC");
    h.ref_start = h.ref_end - (lb + sb);
    h.ref_end = h.ref_end - sb;
  } else {
    h.ref_start += sb;
    h.ref_end = h.ref_start + lb;
  }
  return true;
}

vector<Hit> gap_split(Hit h) {
  struct Gap {
    int start_a, start_b, len_a, len_b;
    int start, len;
  };
  vector<Gap> gaps;
  Gap g{h.aln.start_a, h.aln.start_b, 0, 0, 0, 0};
  for (auto &c : h.aln.cigar) {
    if (c.second && c.first != 'M') {
      if (c.first != 'D')
        g.len_a = 0, g.len_b = c.second;
      else
        g.len_b = 0, g.len_a = c.second;
      g.len = c.second;
      gaps.push_back(g);
    }
    if (c.first != 'D')
      g.start_b += c.second;
    if (c.first != 'I')
      g.start_a += c.second;
    g.start += c.second;
  }
  sort(gaps.begin(), gaps.end(),
       [](const Gap &a, const Gap &b) { return a.len > b.len; });

  Hit hh;
  vector<Hit> hits;
  double gap_score = h.aln.gap_error();
  if (Globals::Stats::MAX_OK_GAP > -1)
    for (auto &g : gaps) {
      dprn("--> {} :: a:{} b:{} ... a:{} b:{}", g.len,
           g.start_a - h.aln.start_a, g.start_b - h.aln.start_b,
           h.aln.end_a - (g.start_a + g.len_a),
           h.aln.end_b - (g.start_b + g.len_b));
      if (g.start_a - h.aln.start_a < Globals::Stats::MIN_SPLIT_SIZE ||
          g.start_b - h.aln.start_b < Globals::Stats::MIN_SPLIT_SIZE)
        continue;
      if (h.aln.end_a - (g.start_a + g.len_a) <
              Globals::Stats::MIN_SPLIT_SIZE ||
          h.aln.end_b - (g.start_b + g.len_b) < Globals::Stats::MIN_SPLIT_SIZE)
        continue;

      double g_score = pct(g.len, h.aln.error.matches + h.aln.error.gap_bases +
                                      h.aln.error.mismatches);

      dprn("{} ~ {}", g_score, g.len);
      if (g_score >= Globals::Stats::MAX_OK_GAP) {
        dprn(":: Found gap of size {} and score {}\n:: {}\n:: {}\n:: {}", g.len,
             g_score, h.aln.align_a.substr(g.start, g.len),
             h.aln.alignment.substr(g.start, g.len),
             h.aln.align_b.substr(g.start, g.len));

        bool x = subhit(h, 0, g.start, hh);
        assert(x);
        auto r = gap_split(hh);
        for (auto &hx : r)
          hits.push_back(hx);

        x = subhit(h, g.start + g.len, h.aln.alignment.size(), hh);
        assert(x);
        r = gap_split(hh);
        for (auto &hx : r)
          hits.push_back(hx);

        return hits;
      }
    }
  if (!hits.size())
    hits.push_back(h);
  return hits;
}

vector<Hit> split_alignment(Hit h) {
  vector<Hit> hits;

  // Find stretch of Ns
  int prev_an = 0, prev_bn = 0;
  int hit_begin = 0;
  Hit hh;
  for (int i = 0; i < h.aln.alignment.size(); i++) {
    if (toupper(h.aln.align_a[i]) == 'N') {
      prev_an++;
    } else {
      if (prev_an >= Globals::Stats::MIN_ASSEMBLY_GAP_SIZE) {
        dprn(":: Found assembly gap of size {}\n:: {}\n:: {}\n:: {}", prev_an,
             h.aln.align_a.substr(i - prev_an, prev_an),
             h.aln.alignment.substr(i - prev_an, prev_an),
             h.aln.align_b.substr(i - prev_an, prev_an));
        if (subhit(h, hit_begin, i - prev_an, hh))
          hits.push_back(hh);
        hit_begin = i;
      }
      prev_an = 0;
    }
    if (toupper(h.aln.align_b[i]) == 'N') {
      prev_bn++;
    } else {
      if (prev_bn >= Globals::Stats::MIN_ASSEMBLY_GAP_SIZE) {
        dprn(":: Found assembly gap of size {}\n:: {}\n:: {}\n:: {}", prev_bn,
             h.aln.align_a.substr(i - prev_bn, prev_bn),
             h.aln.alignment.substr(i - prev_bn, prev_bn),
             h.aln.align_b.substr(i - prev_bn, prev_bn));
        if (subhit(h, hit_begin, i - prev_bn, hh))
          hits.push_back(hh);
        hit_begin = i;
      }
      prev_bn = 0;
    }
  }
  if (!hit_begin)
    hits.push_back(h);
  else if (subhit(h, hit_begin, h.aln.alignment.size(), hh))
    hits.push_back(hh);

  // Find gaps
  vector<Hit> hits_final;
  for (auto &h : hits) {
    auto hh = gap_split(h);
    // dprn("first round: {} ret", hh.size());
    for (auto &hhh : hh) {
      hits_final.push_back(hhh);
    }
  }
  return hits_final;
}

void process(Hit hs, string cigar, FastaReference &fr) {
  string fa = fr.get_sequence(hs.query->name, hs.query_start, &hs.query_end);
  string fb = fr.get_sequence(hs.ref->name, hs.ref_start, &hs.ref_end);
  assert(!hs.query->is_rc);
  if (hs.query->is_rc) {
    fa = rc(fa);
  }
  if (hs.ref->is_rc) {
    fb = rc(fb);
  }
  assert(cigar.size());
  hs.aln = Alignment(fa, fb, cigar);

  auto hs_split = split_alignment(hs);
  for (auto &h : hs_split) {
    if (h.aln.alignment.size() < Globals::Chain::Refine::MIN_READ)
      continue;

    int align_length = h.aln.span();
    int indel_a = 0;
    int indel_b = 0;
    int alignB = 0;
    int matchB = 0;
    int mismatchB = 0;
    int transitionsB = 0;
    int transversionsB = 0;

    int uppercaseA = 0;
    int uppercaseB = 0;
    int uppercaseMatches = 0;

    for (int i = 0; i < align_length; i++) {
      char a = toupper(h.aln.align_a[i]);
      char b = toupper(h.aln.align_b[i]);
      indel_a += a == '-';
      indel_b += b == '-';
      matchB += a != '-' && a == b;
      uppercaseA +=
          (h.aln.align_a[i] != '-' && toupper(h.aln.align_a[i]) != 'N' &&
           isupper(h.aln.align_a[i]));
      uppercaseB +=
          (h.aln.align_b[i] != '-' && toupper(h.aln.align_b[i]) != 'N' &&
           isupper(h.aln.align_b[i]));
      if (a != '-' && b != '-') {
        alignB += 1;
        if (a != b) {
          mismatchB += 1;
          if (a == 'A' || a == 'G') {
            transitionsB += b == 'A' || b == 'G';
            transversionsB += !(b == 'A' || b == 'G');
          } else {
            transitionsB += b == 'C' || b == 'T';
            transversionsB += !(b == 'C' || b == 'T');
          }
        } else if (isupper(h.aln.align_a[i]) && isupper(h.aln.align_b[i])) {
          uppercaseMatches++;
        }
      }
    }

    double fracMatch = double(matchB) / (alignB);
    double fracMatchIndel = double(matchB) / (align_length);

    double jcp = double(mismatchB) / (alignB);
    double jcK = -0.75 * log(1.0 - 4.0 / 3 * jcp);

    double p = double(transitionsB) / (alignB);
    double q = double(transversionsB) / (alignB);
    double w1 = 1.0 / (1 - 2.0 * p - q);
    double w2 = 1.0 / (1 - 2.0 * q);
    double k2K = 0.5 * log(w1) + 0.25 * log(w2);

    // TODO handle this smarter
    bool same_chr =
        h.query->name == h.ref->name && h.query->is_rc == h.ref->is_rc;
    int overlap = !same_chr ? 0
                            : max(0, min(h.query_end, h.ref_end) -
                                         max(h.query_start, h.ref_start));
    bool too_big_overlap = (h.query_end - h.query_start - overlap) <
                               Globals::Stats::BIG_OVERLAP_THRESHOLD ||
                           (h.ref_end - h.ref_start - overlap) <
                               Globals::Stats::BIG_OVERLAP_THRESHOLD;
    too_big_overlap &= same_chr;

    double errorScaled =
        (h.aln.gaps() + h.aln.mismatches()) /
        double(h.aln.gaps() + h.aln.mismatches() + h.aln.matches());

    // Split large gaps?
    if (uppercaseA >= Globals::Stats::MIN_UPPERCASE &&
        uppercaseB >= Globals::Stats::MIN_UPPERCASE && !too_big_overlap &&
        errorScaled <= Globals::Stats::MAX_SCALED_ERROR &&
        uppercaseMatches >= Globals::Stats::MIN_UPPERCASE) {
      // string l = h.to_bed(false);

      h.name = "S";
      h.comment = "";
#pragma omp critical
      prn("{}\t"
          "{}\t{}\t"
          "{}\t{}\t{}\t"
          "{}\t{}\t"
          "{}\t{}\t"
          "{}\t{}\t"
          "{}\t"
          "{}\t{}\t{}\t"
          "{}\t{}\t{}\t"
          "{}\t"
          "{}\t{}",
          h.to_bed(false, false, &fr),              // 1-13
          indel_a, indel_b,                         // 15-16
          alignB, matchB, mismatchB,                // 17-19
          transitionsB, transversionsB,             // 20-21
          fracMatch, fracMatchIndel,                // 22-23
          jcK, k2K,                                 // 24-25
          h.aln.gaps(),                             // 26
          uppercaseA, uppercaseB, uppercaseMatches, // 27-29
          h.aln.matches(), h.aln.mismatches(), h.aln.gaps(),
          h.aln.gap_bases(),                    // 30-33
          h.aln.cigar_string(), 1 - errorScaled // 34-35
      );
    }
  }
}

void stats(const string &ref_path, const string &bed_path) {
  FastaReference fr(ref_path);
  ifstream fin(bed_path.c_str());
  if (!fin.is_open()) {
    throw fmt::format("BED file {} does not exist", bed_path);
  }

  string s;
  vector<pair<Hit, string>> hits;
  while (getline(fin, s)) {
    string cigar;
    Hit h = Hit::from_bed(s, &cigar);

    assert(h.ref != nullptr);
    assert(h.query != nullptr);
    if (tie(h.query->name, h.query_start, h.query_end) >
        tie(h.ref->name, h.ref_start, h.ref_end)) {
      swap(h.query->name, h.ref->name);
      swap(h.query_start, h.ref_start);
      swap(h.query_end, h.ref_end);
      for (auto &c : cigar) {
        if (c == 'I')
          c = 'D';
        else if (c == 'D')
          c = 'I';
      }
    }
    hits.push_back({h, cigar});
  }
  sort(hits.begin(), hits.end(),
       [](const pair<Hit, string> &a, const pair<Hit, string> &b) {
         return tie(a.first.ref->is_rc, a.first.query->name, a.first.ref->name,
                    a.first.query_start, a.first.ref_start) <
                tie(b.first.ref->is_rc, b.first.query->name, b.first.ref->name,
                    b.first.query_start, b.first.ref_start);
       });

  int hit_count = 0, out_count = 0;
  string prev;

  prn("#chr1\tstart1\tend1\tchr2\tstart2\tend2\tname\tscore\tstrand1\tstrand2\t"
      "max_len\taln_len\tcomment\t"                 // 1-13
      "indel_a\tindel_b\talnB\tmatchB\tmismatchB\t" // 14-19
      "transitionsB\ttransversions\tfracMatch\tfracMatchIndel\tjck\tk2K\t" // 20-25
      "aln_gaps\tuppercaseA\tuppercaseB\tuppercaseMatches\t"   // 26-29
      "aln_matches\taln_mismatches\taln_gaps\taln_gap_bases\t" // 30-33
      "cigar\tfilter_score"                                    // 34
  );
#pragma omp parallel for
  for (auto hsi = 0; hsi < hits.size(); hsi++) {
    process(hits[hsi].first, hits[hsi].second, fr);
#pragma omp critical
    eprnn("\rProcessed hit {:n} out of {:n}", ++hit_count, hits.size());
  }
  eprn("\rProcessed hit {:n} out of {:n}... done!", hit_count, hits.size());
}

/******************************************************************************/

void get_differences(const string &ref_path, const string &bed_path,
                     const string &wgac_path) {
  map<string, boost::dynamic_bitset<>> sedef;
  map<string, boost::dynamic_bitset<>> wgac;

  FastaReference fr(ref_path);

  string s;
  ifstream fin(bed_path);
  int q = 0, w = 0;
  while (getline(fin, s)) {
    if (s[0] == '#')
      continue;
    string cigar;
    Hit h = Hit::from_bed(s, &cigar);

    string fa = fr.get_sequence(h.query->name, h.query_start, &h.query_end);
    string fb = fr.get_sequence(h.ref->name, h.ref_start, &h.ref_end);
    int qa = 0;
    for (auto f : fa)
      if (isupper(f))
        qa++;
    int qb = 0;
    for (auto f : fb)
      if (isupper(f))
        qb++;
    if (qa < 100 || qb < 100) {
      w++;
      continue;
    }
    q++;

    auto c1 = fmt::format("{}", h.query->name, "+-"[h.query->is_rc]);
    auto c2 = fmt::format("{}", h.ref->name, "+-"[h.ref->is_rc]);
    if (sedef.find(c1) == sedef.end()) {
      sedef[c1] = boost::dynamic_bitset<>(250 * MB);
    }
    if (sedef.find(c2) == sedef.end()) {
      sedef[c2] = boost::dynamic_bitset<>(250 * MB);
    }
    for (int i = h.query_start; i < h.query_end; i++)
      sedef[c1].set(i);
    for (int i = h.ref_start; i < h.ref_end; i++)
      sedef[c2].set(i);
  }

  eprn("SEDEF reading done (in {:n}, miss {:n})!", q, w);

  ifstream fiw(wgac_path);
  getline(fiw, s);
  unordered_set<string> seen;
  while (getline(fiw, s)) {
    Hit h = Hit::from_wgac(s);
    auto c1 = fmt::format("{}", h.query->name, "+-"[h.query->is_rc]);
    auto c2 = fmt::format("{}", h.ref->name, "+-"[h.ref->is_rc]);
    if (c1.size() > 6 || c2.size() > 6)
      continue;

    if (seen.find(h.name) == seen.end()) {
      seen.insert(h.name);
      if (wgac.find(c1) == wgac.end()) {
        wgac[c1] = boost::dynamic_bitset<>(250 * MB);
      }
      if (wgac.find(c2) == wgac.end()) {
        wgac[c2] = boost::dynamic_bitset<>(250 * MB);
      }
      for (int i = h.query_start; i < h.query_end; i++)
        wgac[c1].set(i);
      for (int i = h.ref_start; i < h.ref_end; i++)
        wgac[c2].set(i);
    }
  }

  eprn("WGAC reading done!");

  int intersect = 0, wgac_only = 0, wgac_span = 0, sedef_only = 0,
      sedef_span = 0;

  int sedef_extra_upper = 0;
  int miss_upper = 0;

  for (auto &p : sedef) {
    auto &s = p.second;
    auto &w = wgac[p.first];

    auto seq = fr.get_sequence(p.first);

    for (int i = 0; i < seq.size(); i++) {
      if ((s[i] & (~w[i])) && isupper(seq[i]) && toupper(seq[i]) != 'N') {
        sedef_extra_upper++;
      }
      if ((w[i] & (~s[i])) && isupper(seq[i]) && toupper(seq[i]) != 'N') {
        miss_upper++;
      }
    }

    intersect += (s & w).count();
    wgac_only += (w & (~s)).count();
    sedef_only += (s & (~w)).count();
    sedef_span += s.count();
    wgac_span += w.count();
  }

  eprn("SEDEF: spans              {:12n}\n"
       "       unique             {:12n}\n"
       "       unique (uppercase) {:12n}\n"
       "       misses             {:12n}\n"
       "       misses (uppercase) {:12n}\n"
       "WGAC:  spans              {:12n}\n"
       "       intersects         {:12n}",
       sedef_span, sedef_only, sedef_extra_upper, wgac_only, miss_upper,
       wgac_span, intersect);
}

/******************************************************************************/

void stats_main(int argc, char **argv) {
  using namespace Globals;
  argh::parser cmdl(argc, argv, argh::parser::PREFER_PARAM_FOR_UNREG_OPTION);

  cmdl("max-ok-gap", Stats::MAX_OK_GAP) >> Stats::MAX_OK_GAP;
  cmdl("min-split", Stats::MIN_SPLIT_SIZE) >> Stats::MIN_SPLIT_SIZE;
  cmdl("uppercase", Stats::MIN_UPPERCASE) >> Stats::MIN_UPPERCASE;
  cmdl("max-error", Stats::MAX_SCALED_ERROR) >> Stats::MAX_SCALED_ERROR;

  if (!cmdl(2)) {
    throw fmt::format("Not enough arguments to stats");
  }

  string command = cmdl[0];
  if (command == "generate") {
    stats(cmdl[1], cmdl[2]);
  } else if (command == "diff") {
    if (!cmdl(3)) {
      throw fmt::format("Not enough arguments to stats");
    }
    get_differences(cmdl[1], cmdl[2], cmdl[3]);
  } else {
    throw fmt::format("Unknown stats command");
  }
}
