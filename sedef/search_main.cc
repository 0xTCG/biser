/// 786

/// This file is subject to the terms and conditions defined in
/// file 'LICENSE', which is part of this source code package.

/// Author: inumanag

/******************************************************************************/

#include <algorithm>
#include <chrono>
#include <fstream>
#include <iostream>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

// #include <boost/asio/io_service.hpp>
// #include <boost/bind.hpp>
// #include <boost/thread/thread.hpp>

#include "common.h"
#include "extern/argh.h"
#include "search.h"
#include "search_main.h"

using namespace std;

/******************************************************************************/

extern int64_t TOTAL_ATTEMPTED;
extern int64_t JACCARD_FAILED;
extern int64_t QGRAM_NORMAL_FAILED;
extern int64_t OTHER_FAILED;
extern int64_t INTERVAL_FAILED;

/******************************************************************************/

template <typename T>
int initial_search(shared_ptr<Index> query_hash, shared_ptr<Index> ref_hash,
                   bool is_same_genome, T print_function,
                   bool show_progress = true) {
  Tree tree;
  int total = 0, track = 0;
  int next_to_attain = 0;

  const int TRACK_PROGRESS = 10000;
  for (int qi = 0; qi < query_hash->minimizers.size(); qi++) {
    auto &qm = query_hash->minimizers[qi];

    if (show_progress && qm.loc / TRACK_PROGRESS != track) {
      eprnn("\r |>{}<| {:.1f}% (loci={:n} hits={:n})",
            string(int(pct(qm.loc, query_hash->seq->seq.size()) / 2) + 1, '-'),
            pct(qm.loc, query_hash->seq->seq.size()), qm.loc, total);
      track = qm.loc / TRACK_PROGRESS;
    }

    if (qm.loc < next_to_attain)
      continue;

    if (Globals::Internal::DoUppercaseSeeds &&
        qm.hash.status != Hash::Status::HAS_UPPERCASE)
      continue;

    auto hits = search(qi, query_hash, ref_hash, tree, is_same_genome,
                       Globals::Search::MIN_READ_SIZE, true, false);
    int min_len = query_hash->seq->seq.size();
    for (auto &pp : hits) {
      min_len = min(min_len, pp.query_end - pp.query_start);
      print_function(pp);
    }
    total += hits.size();

    next_to_attain = (min_len >= Globals::Search::MIN_READ_SIZE
                          ? qm.loc + (Globals::Search::MIN_READ_SIZE *
                                      Globals::Search::MAX_ERROR) /
                                         2
                          : qm.loc);
  }
  return total;
}

/******************************************************************************/

string v2s(const vector<string> &v) {
  string r = "";
  for (auto &s : v)
    r += s + ", ";
  return r.size() ? r.substr(0, r.size() - 2) : r;
};

vector<vector<string>> generate_translation(const string &ref_path,
                                            bool print) {
  FastaReference fr(ref_path);
  eprn("Translating {}...", ref_path);

  vector<pair<size_t, string>> vv;
  for (const auto &rf : fr.index) {
    vv.push_back({rf.second.length, rf.second.name});
  }
  sort(vv.begin(), vv.end(), std::greater<pair<size_t, string>>());

  vector<vector<string>> ref;
  int cur_size = 0;
  const int MAX_SIZE = 100 * MB;
  for (auto &v : vv) {
    if (!ref.size() || cur_size + v.first > MAX_SIZE) {
      ref.push_back({v.second});
      cur_size = v.first;
    } else {
      ref.back().push_back(v.second);
      cur_size += v.first;
    }
  }
  if (print)
    for (int i = 0; i < ref.size(); i++)
      eprn(" [Translate] {} -> {}", i, v2s(ref[i]));
  return ref;
}

void search_single(const string &ref_path, const string &query_chr,
                   const string &ref_chr, bool is_ref_complement, int kmer_size,
                   int window_size, bool transform) {
  eprn("        Parameters: READ_SIZE      = {}\n"
       "                    MAX_ERROR      = {:.2f} ({:.2f} EDIT + {:.2f} GAP; "
       "GAPFREQ={:.3f})",
       Globals::Search::MIN_READ_SIZE, Globals::Search::MAX_ERROR,
       Globals::Search::MAX_EDIT_ERROR,
       Globals::Search::MAX_ERROR - Globals::Search::MAX_EDIT_ERROR,
       Globals::Search::GAP_FREQUENCY);

  eprn("Reverse complement: {}", is_ref_complement);
  eprn("k-mer size:         {}", kmer_size);
  eprn("Window size:        {}", window_size);
  eprn("");

  auto T = cur_time();

  FastaReference fr(ref_path);

  vector<string> qr, rr;
  if (!transform) {
    qr.push_back(query_chr);
    rr.push_back(ref_chr);
  } else {
    auto ref = generate_translation(ref_path);
    qr = ref[std::stoi(query_chr)];
    rr = ref[std::stoi(ref_chr)];
  }

  // q < r ?
  int total = 0;

  map<pair<string, bool>, shared_ptr<Index>> indices;
  for (auto &r : rr) {
    string ref = fr.get_sequence(r);
    indices[{r, is_ref_complement}] =
        make_shared<Index>(make_shared<Sequence>(r, ref, is_ref_complement),
                           kmer_size, window_size);
  }
  for (auto &r : qr) {
    if (indices.find({r, false}) == indices.end()) {
      string query = fr.get_sequence(r);
      indices[{r, false}] = make_shared<Index>(make_shared<Sequence>(r, query),
                                               kmer_size, window_size);
    }
  }
  eprn("Building index took {:.1f}s", elapsed(T)), T = cur_time();

  for (auto &r : rr) {
    auto ref_hash = indices[{r, is_ref_complement}];
    for (auto &q : qr) {
      auto query_hash = indices[{q, false}];
      bool is_same_genome = (q == r) && !is_ref_complement;
      total += initial_search(
          query_hash, ref_hash, is_same_genome,
          [](Hit &h) { prn("{}", h.to_bed()); },
          qr.size() == 1 && rr.size() == 1);
      // if (!(qr.size() == 1 && rr.size() == 1))
      //   eprnn("\r{:10} / {:10} ({:10})", ri * rr.size() + qi,
      //         rr.size() * rr.size(), ri);
    }
  }

  eprn("Total:           = {:10n}", total);
  eprn("Fails: attempts  = {:10n}\n"
       "       Jaccard   = {:10n}\n"
       "       interval  = {:10n}\n"
       "       lowercase = {:10n}\n"
       "       q-grams   = {:10n}",
       TOTAL_ATTEMPTED, JACCARD_FAILED, INTERVAL_FAILED, OTHER_FAILED,
       QGRAM_NORMAL_FAILED);

  exit(0);
}

/******************************************************************************/

void search_main(int argc, char **argv) {
  using namespace Globals;
  argh::parser cmdl;
  cmdl.add_params({
      "-k", "--kmer", "-w", "--window", "-u", "--uppercase", "-e", "--error",
      "-E", "--edit-error", "-g", "--gap-freq", "-l", "--min-read-size",
      // "-t", "--threads", "-m", "--max-parallel-size"
  });
  cmdl.parse(argc, argv);

  cmdl({"-k", "--kmer"}, Search::KMER_SIZE) >> Search::KMER_SIZE;
  cmdl({"-w", "--window"}, Search::WINDOW_SIZE) >> Search::WINDOW_SIZE;
  cmdl({"-u", "--uppercase"}, Search::MIN_UPPERCASE) >> Search::MIN_UPPERCASE;
  cmdl({"-e", "--error"}, Search::MAX_ERROR) >> Search::MAX_ERROR;
  cmdl({"-E", "--edit-error"}, Search::MAX_EDIT_ERROR) >>
      Search::MAX_EDIT_ERROR;
  cmdl({"-g", "--gap-freq"}, Search::GAP_FREQUENCY) >> Search::GAP_FREQUENCY;

  // int threads = -1;
  // cmdl({"-t", "--threads"}, -1) >> threads;
  // size_t max_size = 0;
  // cmdl({"-m", "--max-parallel-size"}, 0) >> max_size;

  Search::MIN_READ_SIZE = KB * (1 - Search::MAX_ERROR); // 700 by default

  bool transform = cmdl[{"-t", "transform"}];
  bool is_complement = cmdl[{"-r", "reverse"}];
  if (!cmdl(2))
    throw fmt::format("Not enough arguments to search");

  search_single(cmdl[0], cmdl[1], cmdl[2], is_complement, Search::KMER_SIZE,
                Search::WINDOW_SIZE, transform);
}

void trans_main(int argc, char **argv) {
  argh::parser cmdl;
  cmdl.parse(argc, argv);
  if (!cmdl(0))
    throw fmt::format("Not enough arguments to translate");
  auto t = generate_translation(cmdl[0], true);
  prn("{}", t.size());
}
