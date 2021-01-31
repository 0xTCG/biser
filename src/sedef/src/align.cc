/// 786

/// This file is subject to the terms and conditions defined in
/// file 'LICENSE', which is part of this source code package.

/// Author: inumanag

/******************************************************************************/

#include <algorithm>
#include <chrono>
#include <fstream>
#include <glob.h>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "align.h"
#include "common.h"
#include "extern/ksw2.h"
#include "fasta.h"
#include "hit.h"

using namespace std;

/******************************************************************************/

inline bool ceq(char a, char b) {
  if (a == '-' || b == '-')
    return false;
  if (toupper(a) == 'N' || toupper(b) == 'N')
    return false;
  return (toupper(a) == toupper(b));
}

/******************************************************************************/

auto align_helper(const string &qseq, const string &tseq, int sc_mch,
                  int sc_mis, int gapo, int gape, int bandwidth) {
  int8_t a = (int8_t)sc_mch,
         b = sc_mis < 0 ? (int8_t)sc_mis : (int8_t)(-sc_mis); // a>0 and b<0
  int8_t mat[25] = {a, b, b, b, 0, b, a, b, b, 0, b, b, a,
                    b, 0, b, b, b, a, 0, 0, 0, 0, 0, 0};
  deque<pair<char, int>> cigar;
  for (int SP = 0; SP < min(tseq.size(), qseq.size());
       SP += Globals::Align::MAX_KSW_SEQ_LEN) {
    ksw_extz_t ez;
    ksw_extz2_sse(0,
                  min(Globals::Align::MAX_KSW_SEQ_LEN, (int)(qseq.size() - SP)),
                  (const uint8_t *)(qseq.c_str() + SP),
                  min(Globals::Align::MAX_KSW_SEQ_LEN, (int)(tseq.size() - SP)),
                  (const uint8_t *)(tseq.c_str() + SP), 5, mat, // M; MxM matrix
                  gapo, gape, bandwidth,
                  -1, // band width; off-diagonal drop-off to stop extension (-1
                      // to disable)
                  0, &ez);
    for (int i = 0; i < ez.n_cigar; i++) {
      int idx = ez.cigar[i] & 0xf;
      int len = ez.cigar[i] >> 4;
      if (idx < 3) {
        cigar.push_back({"MDI"[idx], len});
      }
    }
    free(ez.cigar);
  }
  return cigar;
}

/******************************************************************************/

Alignment::Alignment() : start_a(0), end_a(0), start_b(0), end_b(0) {
  error = {0, 0, 0, 0};
}

Alignment::Alignment(const string &fa, const string &fb)
    : chr_a("A"), start_a(0), end_a(fa.size()), chr_b("B"), start_b(0),
      end_b(fb.size()), a(fa), b(fb) {
  string xa = fa, xb = fb;

  transform(xa.begin(), xa.end(), xa.begin(), align_dna);
  transform(xb.begin(), xb.end(), xb.begin(), align_dna);

  cigar =
      align_helper(xa, xb, Globals::Align::MATCH, Globals::Align::MISMATCH,
                   -Globals::Align::GAP_OPEN, -Globals::Align::GAP_EXTEND, -1);
  populate_nice_alignment();
}

Alignment::Alignment(const string &fa, const string &fb,
                     const string &cigar_str)
    : chr_a("A"), start_a(0), end_a(fa.size()), chr_b("B"), start_b(0),
      end_b(fb.size()), a(fa), b(fb) {
  for (int ci = 0, num = 0; ci < cigar_str.size(); ci++) {
    if (isdigit(cigar_str[ci])) {
      num = 10 * num + (cigar_str[ci] - '0');
    } else if (cigar_str[ci] == ';') {
      continue;
    } else {
      cigar.push_back({cigar_str[ci], num});
      num = 0;
    }
  }
  populate_nice_alignment();
}

Alignment::Alignment(const string &qstr, const string &rstr,
                     const vector<Hit> &guide, const int side) {
  auto prev = guide.begin();
  *this = prev->aln;
  for (auto cur = next(prev); cur != guide.end(); cur++) {
    int qs = cur->query_start, qe = cur->query_end;
    int qps = prev->query_start, qpe = prev->query_end;

    int rs = cur->ref_start, re = cur->ref_end;
    int rps = prev->ref_start, rpe = prev->ref_end;

    assert(qpe <= qs);
    assert(rpe <= rs);

    end_a = qe;
    end_b = re;
    a += qstr.substr(qpe, qe - qpe);
    b += rstr.substr(rpe, re - rpe);

    int qgap = qs - qpe, rgap = rs - rpe;
    if (qgap && rgap) {
      if (qgap <= 1000 && rgap <= 1000) { // "close" hits
        Alignment gap(qstr.substr(qpe, qgap), rstr.substr(rpe, rgap));
        append_cigar(gap.cigar);
      } else { // assume only one part is the gap
        int ma = max(qgap, rgap);
        int mi = min(qgap, rgap);
        Alignment ma1(qstr.substr(qpe, mi), rstr.substr(rpe, mi));
        ma1.cigar.push_back({qgap == mi ? 'I' : 'D', ma - mi});
        Alignment ma2(qstr.substr(qs - mi, mi), rstr.substr(rs - mi, mi));
        ma2.cigar.push_front({qgap == mi ? 'I' : 'D', ma - mi});
        append_cigar(ma2.total_error() < ma2.total_error() ? ma2.cigar
                                                           : ma1.cigar);
      }
    } else if (qgap) {
      append_cigar({{'D', qgap}});
    } else if (rgap) {
      append_cigar({{'I', rgap}});
    }
    append_cigar(cur->aln.cigar);
    prev = cur;
  }

  int qlo = start_a, qhi = end_a;
  int rlo = start_b, rhi = end_b;
  assert(a == qstr.substr(qlo, qhi - qlo));
  assert(b == rstr.substr(rlo, rhi - rlo));

  if (side) {
    int qlo_n = max(0, qlo - side);
    int rlo_n = max(0, rlo - side);
    if (qlo - qlo_n && rlo - rlo_n) {
      Alignment gap(qstr.substr(qlo_n, qlo - qlo_n),
                    rstr.substr(rlo_n, rlo - rlo_n));
      gap.trim_front();

      qlo_n = qlo - (gap.end_a - gap.start_a);
      rlo_n = rlo - (gap.end_b - gap.start_b);
      prepend_cigar(gap.cigar);
      a = qstr.substr(qlo_n, qlo - qlo_n) + a;
      b = rstr.substr(rlo_n, rlo - rlo_n) + b;
      start_a = qlo = qlo_n;
      start_b = rlo = rlo_n;
    }

    int qhi_n = min(qhi + side, (int)qstr.size());
    int rhi_n = min(rhi + side, (int)rstr.size());
    if (qhi_n - qhi && rhi_n - rhi) {
      Alignment gap(qstr.substr(qhi, qhi_n - qhi),
                    rstr.substr(rhi, rhi_n - rhi));
      gap.trim_back();

      qhi_n = qhi + gap.end_a;
      rhi_n = rhi + gap.end_b;
      append_cigar(gap.cigar);
      a += qstr.substr(qhi, qhi_n - qhi);
      b += rstr.substr(rhi, rhi_n - rhi);
      end_a = qhi = qhi_n;
      end_b = rhi = rhi_n;
    }
  }

  assert(qlo >= 0);
  assert(rlo >= 0);
  assert(qhi <= qstr.size());
  assert(rhi <= rstr.size());
  assert(a == qstr.substr(qlo, qhi - qlo));
  assert(b == rstr.substr(rlo, rhi - rlo));

  populate_nice_alignment();
}

Alignment::Alignment(const string &qstr, const string &rstr,
                     const vector<Anchor> &guide, const vector<int> &guide_idx)
    : chr_a("A"), chr_b("B") {
  if (guide_idx.size() == 0) {
    *this = Alignment();
    return;
  }

  auto prev = guide_idx.begin();
  start_a = guide[*prev].q;
  end_a = guide[*prev].q + guide[*prev].l;
  start_b = guide[*prev].r;
  end_b = guide[*prev].r + guide[*prev].l;
  a = qstr.substr(start_a, end_a - start_a);
  b = rstr.substr(start_b, end_b - start_b);
  cigar = {{'M', end_a - start_a}};
  assert(end_a - start_a == end_b - start_b);

  for (auto cur = next(prev); cur != guide_idx.end(); cur++) {
    int qs = guide[*cur].q, qe = guide[*cur].q + guide[*cur].l;
    int qps = guide[*prev].q, qpe = guide[*prev].q + guide[*prev].l;

    int rs = guide[*cur].r, re = guide[*cur].r + guide[*cur].l;
    int rps = guide[*prev].r, rpe = guide[*prev].r + guide[*prev].l;

    assert(qpe <= qs);
    assert(rpe <= rs);

    end_a = qe;
    end_b = re;
    a += qstr.substr(qpe, qe - qpe);
    b += rstr.substr(rpe, re - rpe);

    int qgap = qs - qpe, rgap = rs - rpe;
    if (qgap && rgap) {
      if (qgap <= 1000 && rgap <= 1000) { // "close" hits
        Alignment gap(qstr.substr(qpe, qgap), rstr.substr(rpe, rgap));
        append_cigar(gap.cigar);
      } else { // assume only one part is the gap
        int ma = max(qgap, rgap);
        int mi = min(qgap, rgap);
        Alignment ma1(qstr.substr(qpe, mi), rstr.substr(rpe, mi));
        ma1.cigar.push_back({qgap == mi ? 'I' : 'D', ma - mi});
        Alignment ma2(qstr.substr(qs - mi, mi), rstr.substr(rs - mi, mi));
        ma2.cigar.push_front({qgap == mi ? 'I' : 'D', ma - mi});
        append_cigar(ma2.total_error() < ma2.total_error() ? ma2.cigar
                                                           : ma1.cigar);
      }
    } else if (qgap) {
      append_cigar({{'D', qgap}});
    } else if (rgap) {
      append_cigar({{'I', rgap}});
    }
    // eprn("");
    assert(qe - qs == re - rs);
    append_cigar({{'M', qe - qs}});
    prev = cur;
  }

  int qlo = start_a, qhi = end_a;
  int rlo = start_b, rhi = end_b;
  assert(a == qstr.substr(qlo, qhi - qlo));
  assert(b == rstr.substr(rlo, rhi - rlo));
  assert(qlo >= 0);
  assert(rlo >= 0);
  assert(qhi <= qstr.size());
  assert(rhi <= rstr.size());
  assert(a == qstr.substr(qlo, qhi - qlo));
  assert(b == rstr.substr(rlo, rhi - rlo));

  populate_nice_alignment();
}

/******************************************************************************/

void Alignment::populate_nice_alignment() {
  align_a = "";
  align_b = "";
  alignment = "";
  int ia = 0, ib = 0;
  for (auto &c : cigar) {
    for (int i = 0; i < c.second; i++) {
      assert(c.first != 'M' || ia < a.size());
      assert(c.first != 'M' || ib < b.size());
      if (c.first == 'M' && ceq(a[ia], b[ib])) {
        alignment += "|";
      } else {
        alignment += "*";
      }
      if (c.first != 'D')
        align_b += b[ib++];
      else
        align_b += "-";
      if (c.first != 'I')
        align_a += a[ia++];
      else
        align_a += "-";
    }
  }

  error = AlignmentError{0, 0, 0, 0};
  for (auto &c : cigar) {
    if (c.first != 'M') {
      error.gaps++;
      error.gap_bases += c.second;
    }
  }
  for (int i = 0; i < alignment.size(); i++) {
    if (align_a[i] != '-' && align_b[i] != '-') {
      if (ceq(align_a[i], align_b[i])) {
        error.matches++;
      } else {
        error.mismatches++;
      }
    }
  }
}

void Alignment::trim() {
  while (cigar.size()) {
    if (cigar[0].first == 'D') {
      a = a.substr(cigar[0].second);
      start_a += cigar[0].second;
      cigar.pop_front();
    } else if (cigar[0].first == 'I') {
      b = b.substr(cigar[0].second);
      start_b += cigar[0].second;
      cigar.pop_front();
    } else if (cigar.back().first == 'D') {
      end_a -= cigar.back().second;
      a = a.substr(0, a.size() - cigar.back().second);
      cigar.pop_back();
    } else if (cigar.back().first == 'I') {
      end_b -= cigar.back().second;
      b = b.substr(0, b.size() - cigar.back().second);
      cigar.pop_back();
    } else {
      break;
    }
  }

  populate_nice_alignment();
}

void Alignment::trim_front() // ABCD -> --CD
{
  int max_score = 0, max_i = a.size();
  int score = 0;
  for (int i = alignment.size() - 1; i >= 0; i--) {
    if (alignment[i] == '|') {
      score += Globals::Align::MATCH;
    } else {
      if (align_a[i] != '-' && align_b[i] != '-') {
        score += Globals::Align::MISMATCH;
      } else {
        if (i == alignment.size() - 1 ||
            (align_a[i] == '-' && align_a[i + 1] != '-') ||
            (align_b[i] == '-' && align_b[i + 1] != '-')) {
          score += Globals::Align::GAP_OPEN;
        }
        score += Globals::Align::GAP_EXTEND;
      }
    }
    if (score >= max_score) {
      max_score = score, max_i = i;
    }
  }
  if (max_i == a.size()) {
    a = "";
    b = "";
    start_a = end_a;
    start_b = end_b;
    cigar.clear();
    return;
  }
  for (int ci = 0, cur_len = 0; ci < cigar.size(); ci++) {
    if (cigar[ci].second + cur_len > max_i) {
      assert(cigar[ci].first == 'M');
      int need = max_i - cur_len;
      cigar[ci].second -= need;
      for (int cj = 0; cj < ci; cj++)
        cigar.pop_front();
      start_a += need;
      start_b += need;
      break;
    }
    cur_len += cigar[ci].second;
    if (cigar[ci].first == 'M') {
      start_a += cigar[ci].second;
      start_b += cigar[ci].second;
    } else if (cigar[ci].first == 'I') {
      start_b += cigar[ci].second;
    } else {
      start_a += cigar[ci].second;
    }
  }
  a = a.substr(start_a, end_a - start_a);
  b = b.substr(start_b, end_b - start_b);
  populate_nice_alignment();
}

void Alignment::trim_back() // ABCD -> AB--
{
  int max_score = 0, max_i = -1;
  int score = 0;
  for (int i = 0; i < alignment.size(); i++) {
    if (alignment[i] == '|') {
      score += Globals::Align::MATCH;
    } else {
      if (align_a[i] != '-' && align_b[i] != '-') {
        score += Globals::Align::MISMATCH;
      } else {
        if (i == 0 || (align_a[i] == '-' && align_a[i - 1] != '-') ||
            (align_b[i] == '-' && align_b[i - 1] != '-')) {
          score += Globals::Align::GAP_OPEN;
        }
        score += Globals::Align::GAP_EXTEND;
      }
    }
    if (score >= max_score) {
      max_score = score, max_i = i;
    }
  }
  if (max_i == -1) {
    a = "";
    b = "";
    end_a = start_a;
    end_b = start_b;
    cigar.clear();
    return;
  }
  max_i++;
  end_a = start_a, end_b = start_b;
  for (int ci = 0, cur_len = 0; ci < cigar.size(); ci++) {
    if (cigar[ci].second + cur_len >= max_i) {
      assert(cigar[ci].first == 'M');
      int need = max_i - cur_len;
      cigar[ci].second = need;
      while (cigar.size() - 1 > ci)
        cigar.pop_back();
      end_a += need;
      end_b += need;
      break;
    }
    cur_len += cigar[ci].second;
    if (cigar[ci].first == 'M') {
      end_a += cigar[ci].second;
      end_b += cigar[ci].second;
    } else if (cigar[ci].first == 'I') {
      end_b += cigar[ci].second;
    } else {
      end_a += cigar[ci].second;
    }
  }
  a = a.substr(start_a, end_a - start_a);
  b = b.substr(start_b, end_b - start_b);
  populate_nice_alignment();
}

void Alignment::prepend_cigar(const deque<pair<char, int>> &app) {
  if (!app.size())
    return;
  if (cigar.size() && cigar.front().first == app.back().first) {
    cigar.front().second += app.back().second;
    cigar.insert(cigar.begin(), app.begin(), app.begin() + (app.size() - 1));
  } else {
    cigar.insert(cigar.begin(), app.begin(), app.end());
  }
}

void Alignment::append_cigar(const deque<pair<char, int>> &app) {
  if (!app.size())
    return;
  if (cigar.size() && cigar.back().first == app.front().first) {
    cigar.back().second += app.front().second;
    cigar.insert(cigar.end(), next(app.begin()), app.end());
  } else {
    cigar.insert(cigar.end(), app.begin(), app.end());
  }
}

void Alignment::cigar_from_alignment() {
  cigar.clear();
  int sz = 0;
  char op = 0, top;
  for (int i = 0; i < alignment.size(); i++) {
    if (align_a[i] == '-') {
      top = 'I';
    } else if (align_b[i] == '-') {
      top = 'D';
    } else {
      top = 'M';
    }

    if (op != top) {
      if (op)
        cigar.push_back(make_pair(op, sz));
      op = top, sz = 0;
    }
    sz++;
  }
  cigar.push_back(make_pair(op, sz));
}

/******************************************************************************/

void Alignment::merge(Alignment &cur, const string &qstr, const string &rstr) {
  assert(cur.start_a < end_a || cur.start_b < end_b);
  assert(end_a <= cur.end_a);
  assert(end_b <= cur.end_b);

  int trim = end_a - cur.start_a;
  int q = 0, r = 0, i = 0;
  for (i = alignment.size() - 1; i >= 0 && q < trim; i--) {
    if (align_a[i] != '-')
      q++;
    if (align_b[i] != '-')
      r++;
  }
  align_a = align_a.substr(0, i + 1);
  alignment = alignment.substr(0, i + 1);
  align_b = align_b.substr(0, i + 1);
  end_a = start_a + a.size() - q;
  end_b = start_b + b.size() - r;
  a = a.substr(0, a.size() - q);
  b = b.substr(0, b.size() - r);

  q = 0, r = 0, i = 0;
  for (i = 0; i < cur.alignment.size() && q < trim; i++) {
    if (cur.align_a[i] != '-')
      q++;
    if (cur.align_b[i] != '-')
      r++;
  }
  cur.align_a = cur.align_a.substr(i);
  cur.alignment = cur.alignment.substr(i);
  cur.align_b = cur.align_b.substr(i);
  cur.start_a += q;
  cur.start_b += r;
  cur.a = cur.a.substr(q);
  cur.b = cur.b.substr(r);

  trim = end_b - cur.start_b;
  q = 0, r = 0, i = 0;
  for (i = alignment.size() - 1; i >= 0 && r < trim; i--) {
    if (align_a[i] != '-')
      q++;
    if (align_b[i] != '-')
      r++;
  }
  align_a = align_a.substr(0, i + 1);
  alignment = alignment.substr(0, i + 1);
  align_b = align_b.substr(0, i + 1);
  end_a = start_a + a.size() - q;
  end_b = start_b + b.size() - r;
  a = a.substr(0, a.size() - q);
  b = b.substr(0, b.size() - r);

  q = 0, r = 0;
  for (i = 0; i < cur.alignment.size() && r < trim; i++) {
    if (cur.align_a[i] != '-')
      q++;
    if (cur.align_b[i] != '-')
      r++;
  }
  cur.align_a = cur.align_a.substr(i);
  cur.alignment = cur.alignment.substr(i);
  cur.align_b = cur.align_b.substr(i);
  cur.start_a += q;
  cur.start_b += r;
  cur.a = cur.a.substr(q);
  cur.b = cur.b.substr(r);

  cigar_from_alignment();
  cur.cigar_from_alignment();

  assert(start_a <= cur.start_a);
  assert(start_b <= cur.start_b);
  assert(end_a <= cur.start_a);
  assert(end_b <= cur.start_b);
  int qgap = cur.start_a - end_a;
  int rgap = cur.start_b - end_b;
  if (qgap && rgap) {
    if (qgap <= 1000 && rgap <= 1000) { // "close" hits
      Alignment gap(qstr.substr(end_a, qgap), rstr.substr(end_b, rgap));
      append_cigar(gap.cigar);
    } else {                    // assume only one part is the gap
      int ma = max(qgap, rgap); // gap
      int mi = min(qgap, rgap); // mismatch
      Alignment ma1(qstr.substr(end_a, mi), rstr.substr(end_b, mi));
      ma1.cigar.push_back({qgap == mi ? 'I' : 'D', ma - mi});
      Alignment ma2(qstr.substr(cur.start_a - mi, mi),
                    rstr.substr(cur.start_b - mi, mi));
      ma2.cigar.push_front({qgap == mi ? 'I' : 'D', ma - mi});
      append_cigar(ma2.total_error() < ma2.total_error() ? ma2.cigar
                                                         : ma1.cigar);
    }
  } else if (qgap) {
    append_cigar({{'D', qgap}});
  } else if (rgap) {
    append_cigar({{'I', rgap}});
  }

  a += qstr.substr(end_a, qgap) + cur.a;
  b += rstr.substr(end_b, rgap) + cur.b;
  assert(cur.end_a >= end_a);
  assert(cur.end_b >= end_b);
  end_a = cur.end_a;
  end_b = cur.end_b;
  append_cigar(cur.cigar);
  populate_nice_alignment();
}

/******************************************************************************/

string Alignment::cigar_string() const {
  string res;
  for (auto &p : cigar)
    if (p.second) {
      res += fmt::format("{}{}", p.second, p.first);
    }
  return res;
}

void Alignment::swap() {
  std::swap(a, b);
  std::swap(chr_a, chr_b);
  std::swap(start_a, start_b);
  std::swap(end_a, end_b);
  for (auto &p : cigar)
    if (p.second) {
      if (p.first == 'I')
        p.first = 'D';
      else if (p.first == 'D')
        p.first = 'I';
    }
  populate_nice_alignment();
}

string Alignment::print(int width, bool only_alignment) const {
  assert(alignment.size());

  string res;
  int qa = start_a, sa = 0;
  int qb = start_b, sb = 0;

  if (width == -1) {
    width = alignment.size();
  }

  if (!only_alignment) {
    res += fmt::format(
        "       A: {:>9}..{:<9} (len {:7})    Gaps:       {:5} = {:.0f}% ({})\n"
        "       B: {:>9}..{:<9} (len {:7})    Mismatches: {:5} = {:.0f}%\n"
        "   CIGAR: {}\n",
        start_a, end_a, end_a - start_a, error.gap_bases, gap_error(),
        error.gaps, start_b, end_b, end_b - start_b, error.mismatches,
        mismatch_error(), cigar_string());
  }
  for (int i = 0; i < alignment.size(); i += width) {
    if (only_alignment) {
      res += fmt::format("{}\n{}\n{}\n\n", align_a.substr(i, width),
                         alignment.substr(i, width), align_b.substr(i, width));
    } else {
      res += fmt::format("   {:10}: {} {}\n   {:10}  {} {}\n   {:10}: {} {}\n",
                         qa, align_a.substr(i, width), sa, "",
                         alignment.substr(i, width),
                         i + align_a.substr(i, width).size(), qb,
                         align_b.substr(i, width), sb);
    }
    for (auto c : align_a.substr(i, width))
      if (c != '-')
        qa++, sa++;
    for (auto c : align_b.substr(i, width))
      if (c != '-')
        qb++, sb++;
  }
  return res;
}
