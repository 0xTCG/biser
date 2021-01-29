/// 786

/// This file is subject to the terms and conditions defined in
/// file 'LICENSE', which is part of this source code package.

/// Author: inumanag

/******************************************************************************/

#include <cmath>
#include <list>
#include <queue>
#include <string>
#include <vector>

#include "common.h"
#include "filter.h"

using namespace std;

/******************************************************************************/

/* extern */ int64_t QGRAM_NORMAL_FAILED = 0;
/* extern */ int64_t OTHER_FAILED = 0;

/******************************************************************************/

inline int min_qgram(int l, int q) {
  return l * (1 -
              (Globals::Search::MAX_ERROR - Globals::Search::MAX_EDIT_ERROR) -
              q * Globals::Search::MAX_EDIT_ERROR) -
         (Globals::Search::GAP_FREQUENCY * l + 1) * (q - 1);
}

/******************************************************************************/

pair<bool, string> uppercase_filter(const string &q, int q_pos, int q_len,
                                    const string &r, int r_pos, int r_len) {
  int q_up = 0;
  for (int i = 0; i < q_len; i++)
    q_up += (bool)isupper(q[q_pos + i]);
  int r_up = 0;
  for (int i = 0; i < r_len; i++)
    r_up += (bool)isupper(r[r_pos + i]);

  if (q_up < Globals::Search::MIN_UPPERCASE ||
      r_up < Globals::Search::MIN_UPPERCASE) {
    OTHER_FAILED++;
    return {false, fmt::format("upper ({}, {}) < {}", q_up, r_up,
                               Globals::Search::MIN_UPPERCASE)};
  }
  return {true, ""};
}

pair<bool, string> qgram_filter(const string &q, int q_pos, int q_len,
                                const string &r, int r_pos, int r_len) {
  int maxlen = max(q_len, r_len);
  int QG = 5;
  uint64_t QSZ = (1 << (2 * QG));
  uint64_t MASK = QSZ - 1;

  vector<int> qgram_p(QSZ, 0), qgram_r(QSZ, 0);

  int minqg = min_qgram(maxlen, QG);
  assert(minqg >= 10);

  assert(q_pos + q_len <= q.size());
  assert(r_pos + r_len <= r.size());

  for (uint64_t qi = q_pos, qgram = 0; qi < q_pos + q_len; qi++) {
    qgram = ((qgram << 2) | hash_dna(q[qi])) & MASK;
    if (qi - q_pos >= QG - 1)
      qgram_p[qgram] += 1;
  }
  for (uint64_t qi = r_pos, qgram = 0; qi < r_pos + r_len; qi++) {
    qgram = ((qgram << 2) | hash_dna(r[qi])) & MASK;
    if (qi - r_pos >= QG - 1)
      qgram_r[qgram] += 1;
  }
  int dist = 0;
  for (uint64_t qi = 0; qi < QSZ; qi++) {
    dist += min(qgram_p[qi], qgram_r[qi]);
    qgram_p[qi] = qgram_r[qi] = 0;
  }

  if (dist < minqg) {
    QGRAM_NORMAL_FAILED++;
    return {false, fmt::format("q-grams {} < {}", dist, minqg)};
  }
  return {true, ""};
}

/******************************************************************************/

pair<bool, string> filter(const string &q, int q_pos, int q_end,
                          const string &r, int r_pos, int r_end) {
  if (Globals::Internal::DoUppercase) {
    auto f = uppercase_filter(q, q_pos, q_end - q_pos, r, r_pos, r_end - r_pos);
    if (!f.first)
      return f;
  }

  if (Globals::Internal::DoQgram) {
    auto f = qgram_filter(q, q_pos, q_end - q_pos, r, r_pos, r_end - r_pos);
    if (!f.first)
      return f;
  }

  return {true, ""};
}
