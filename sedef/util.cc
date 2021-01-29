/// 786

/// This file is subject to the terms and conditions defined in
/// file 'LICENSE', which is part of this source code package.

/// Author: inumanag

/******************************************************************************/

#include <sstream>
#include <time.h>
#include <unordered_map>

#include <sys/stat.h>
#include <sys/types.h>

#include <boost/math/distributions/binomial.hpp>
#include <boost/math/tools/roots.hpp>

#include "common.h"

using namespace std;

/******************************************************************************/

mode_t stat_file(const string &path) {
  struct stat path_stat;
  int s = stat(path.c_str(), &path_stat);
  assert(s == 0);
  return path_stat.st_mode;
}

vector<string> split(const string &s, char delim) {
  vector<string> elems;
  stringstream ss(s);
  string item;
  while (getline(ss, item, delim)) {
    elems.push_back(item);
  }
  return elems;
}

string rc(const string &s) {
  auto r = s;
  reverse(r.begin(), r.end());
  transform(r.begin(), r.end(), r.begin(), rev_dna);
  return r;
}

/******************************************************************************/

double tau(double edit_error, int kmer_size) {
  const double ERROR_RATIO =
      (Globals::Search::MAX_ERROR - Globals::Search::MAX_EDIT_ERROR) /
      Globals::Search::MAX_EDIT_ERROR;
  double gap_error = std::min(1.0, ERROR_RATIO * edit_error);
  double a = (1 - gap_error) / (1 + gap_error);
  double b = 1 / (2 * std::exp(kmer_size * edit_error) - 1);
  return a * b;
}

double solve_inverse_jaccard(int j, int kmer_size) {
  if (j == 0)
    return 1;
  if (j == 1)
    return 0;
  return boost::math::tools::newton_raphson_iterate(
      [j, kmer_size](double d) {
        const double ERROR_RATIO =
            (Globals::Search::MAX_ERROR - Globals::Search::MAX_EDIT_ERROR) /
            Globals::Search::MAX_EDIT_ERROR;
        double E = exp(d * kmer_size);
        return make_tuple(((1 - d * ERROR_RATIO) / (1 + d * ERROR_RATIO)) *
                                  (1.0 / (2 * E - 1)) -
                              j,
                          2 *
                              (-kmer_size * E + ERROR_RATIO -
                               2 * ERROR_RATIO * E +
                               E * kmer_size * pow(d * ERROR_RATIO, 2)) /
                              pow((2 * E - 1) * (1 + d * ERROR_RATIO), 2));
      },
      0.10, 0.0, 1.0, numeric_limits<double>::digits);
}

int relaxed_jaccard_estimate(int s, int kmer_size,
                             unordered_map<int, int> &mm) {
  double result = -1;
  auto it = mm.find(s);
  if (it != mm.end())
    result = it->second;
  if (result != -1)
    return result;

  using namespace boost::math;
  const double CI = 0.75;
  const double Q2 = (1.0 - CI) / 2; // one side interval probability

  result = ceil(s * tau(Globals::Search::MAX_EDIT_ERROR, kmer_size));
  for (; result >= 0; result--) {
    double d =
        solve_inverse_jaccard(result / s, kmer_size); // returns edit error
    double x = quantile(
        complement(binomial(s, tau(d, kmer_size)), Q2)); // inverse binomial
    double low_d = solve_inverse_jaccard(x / s, kmer_size);
    if (100 * (1 - low_d) < Globals::Search::MAX_EDIT_ERROR) {
      result++;
      break;
    }
  }
  result = max(result, 0.0);
  mm[s] = result;
  return result;
}
