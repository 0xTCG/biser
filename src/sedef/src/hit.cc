/// 786

/// This file is subject to the terms and conditions defined in
/// file 'LICENSE', which is part of this source code package.

/// Author: inumanag

/******************************************************************************/

#include <algorithm>
#include <chrono>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include <glob.h>

#include "align.h"
#include "common.h"
#include "fasta.h"
#include "hit.h"

using namespace std;

/******************************************************************************/

Hit Hit::from_bed(const string &bed, std::string *cigar) {
  auto ss = split(bed, '\t');
  assert(ss.size() >= 10);

  Hit h{make_shared<Sequence>(ss[0], "", ss[8][0] != '+'),
        0,
        0,
        make_shared<Sequence>(ss[3], "", ss[9][0] != '+'),
        0,
        0,
        0,
        "",
        "",
        {}};

  h.query_start = atoi(ss[1].c_str());
  h.query_end = atoi(ss[2].c_str());

  h.ref_start = atoi(ss[4].c_str());
  h.ref_end = atoi(ss[5].c_str());

  h.name = ss[6];
  if (ss.size() >= 15) {
    h.comment = ss[14];
  }
  if (ss.size() >= 14) {
    h.jaccard = atoi(ss[13].c_str());
  }

  if (ss.size() >= 13 && cigar != nullptr) {
    *cigar = ss[12];
  }

  return h;
}

Hit Hit::from_bed(const string &bed, shared_ptr<Sequence> query,
                  shared_ptr<Sequence> ref) {
  auto ss = split(bed, '\t');
  assert(ss.size() >= 10);

  Hit h{query, 0, 0, ref, 0, 0, 0, "", "", {}};

  h.query_start = atoi(ss[1].c_str());
  h.query_end = atoi(ss[2].c_str());

  h.ref_start = atoi(ss[4].c_str());
  h.ref_end = atoi(ss[5].c_str());

  assert(h.query->is_rc == (ss[8][0] != '+'));
  assert(h.ref->is_rc == (ss[9][0] != '+'));

  assert(!h.query->is_rc);
  if (ref->is_rc) {
    swap(h.ref_start, h.ref_end);
    h.ref_start = ref->seq.size() - h.ref_start + 1;
    h.ref_end = ref->seq.size() - h.ref_end + 1;
  }

  h.name = ss[6];
  if (ss.size() >= 14) {
    h.comment = ss[13];
  }
  if (ss.size() >= 13) {
    h.aln = Alignment(query->seq, ref->seq, ss[12]);
  }

  return h;
}

Hit Hit::from_wgac(const string &bed) {
  auto ss = split(bed, '\t');
  assert(ss.size() >= 27);

  Hit h{make_shared<Sequence>(ss[0], "", false),
        atoi(ss[1].c_str()),
        atoi(ss[2].c_str()),
        make_shared<Sequence>(ss[6], "", ss[5][0] != '+'),
        atoi(ss[7].c_str()),
        atoi(ss[8].c_str()),
        0,
        ss[16],
        fmt::format("err={:.1f}", 100 - 100 * atof(ss[26].c_str())),
        {}};

  assert(h.ref->is_rc == (ss[5][0] != '+'));
  assert(!h.query->is_rc);

  return h;
}

/******************************************************************************/

int get_position(const vector<pair<size_t, string>> &ar, size_t ppos) {
  auto lb = lower_bound(ar.begin(), ar.end(), make_pair(ppos, string()));
  if (lb == ar.end())
    return ar.size() - 1;
  else if (lb->first == ppos)
    return lb - ar.begin();
  else {
    assert(lb != ar.begin());
    return lb - ar.begin() - 1;
  }
}

string Hit::to_bed(bool do_rc, bool with_cigar,
                   const FastaReference *fr) const {
  assert(!query->is_rc);

  string qn = query->name;
  int qs = query_start, qe = query_end;
  string rn = ref->name;
  int rs = do_rc && ref->is_rc ? ref->seq.size() - ref_end + 1 : ref_start;
  int re = do_rc && ref->is_rc ? ref->seq.size() - ref_start + 1 : ref_end;

  if (fr && fr->translation_index.size()) {
    const auto p = fr->translation_index.find(qn);
    assert(p != fr->translation_index.end());

    int pos = get_position(p->second, qs);

    // if (pos != p->second.size() - 1 && p->second[pos + 1].first <= qe) {
    // 	eprn("{} {} {}", qn, qs, qe);
    // 	eprn("{} {} ... {} {}",
    // 	 p->second[pos].first, p->second[pos].second,
    // 	 p->second[pos + 1].first, p->second[pos + 1].second);
    // }

    qn = p->second[pos].second;
    qs -= p->second[pos].first;
    assert(pos == p->second.size() - 1 || p->second[pos + 1].first > qe);
    qe -= p->second[pos].first;
  }
  if (fr && fr->translation_index.size()) {
    const auto p = fr->translation_index.find(rn);
    assert(p != fr->translation_index.end());
    int pos = get_position(p->second, rs);

    rn = p->second[pos].second;
    rs -= p->second[pos].first;
    assert(pos == p->second.size() - 1 || p->second[pos + 1].first > re);
    re -= p->second[pos].first;
  }

  return fmt::format("{}\t{}\t{}\t" // QUERY 0 1 3
                     "{}\t{}\t{}\t" // REF   3 4 5
                     "{}\t{}\t"     // NAME 6 SCORE 7
                     "{}\t{}\t"     // STRAND 8 STRAND 9
                     "{}\t{}\t"     // MAXLEN 10 ALNLEN 11
                     "{}{}",        // CIGAR 12 COMMENT 13
                     qn, qs, qe, rn, rs, re, name,
                     aln.span() ? fmt::format("{:.1f}", aln.total_error()) : "",
                     query->is_rc ? "-" : "+", ref->is_rc ? "-" : "+",
                     // Optional fields
                     // - Max. span
                     // - Aln. length
                     // - CIGAR
                     // - Jaccard similarity
                     // - Reason
                     max(query_end - query_start, ref_end - ref_start),
                     aln.span(), with_cigar ? aln.cigar_string() + "\t" : "",
                     fmt::format("{}{}",
                                 aln.span() ? fmt::format("m={:.1f};g={:.1f}",
                                                          aln.mismatch_error(),
                                                          aln.gap_error())
                                            : "",
                                 comment.size() ? ";" + comment : ""));
}

/******************************************************************************/

void Hit::extend(const double factor, const int max_extend) {
  int w = max(query_end - query_start, ref_end - ref_start);
  w = min(max_extend, int(factor * w));
  query_start = max(0, query_start - w);
  query_end += w;
  ref_start = max(0, ref_start - w);
  ref_end += w;
}

/******************************************************************************/

void update_from_alignment(Hit &h) {
  h.query_start = h.aln.start_a;
  h.query_end = h.aln.end_a;
  h.ref_start = h.aln.start_b;
  h.ref_end = h.aln.end_b;
}
