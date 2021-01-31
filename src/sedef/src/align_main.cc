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
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include "align.h"
#include "align_main.h"
#include "chain.h"
#include "common.h"
#include "extern/argh.h"
#include "fasta.h"
#include "hit.h"
#include "merge.h"
#include "search_main.h"
#include "util.h"

using namespace std;

/******************************************************************************/

void bucket_alignments_extern(const string &bed_path, int nbins,
                              string output_dir, bool extend,
                              const string &reference) {

  auto ref = generate_translation(reference);
  map<string, int> lookup;
  for (int i = 0; i < ref.size(); i++)
    for (auto &j : ref[i])
      lookup[j] = i;

  vector<string> files;
  if (S_ISREG(stat_file(bed_path))) {
    files.push_back(bed_path);
  } else if (S_ISDIR(stat_file(bed_path))) {
    glob_t glob_result;
    glob((bed_path + "/*.bed").c_str(), GLOB_TILDE, NULL, &glob_result);
    for (int i = 0; i < glob_result.gl_pathc; i++) {
      string f = glob_result.gl_pathv[i];
      if (S_ISREG(stat_file(f))) {
        files.push_back(f);
      }
    }
  } else {
    throw fmt::format("Path {} is neither file nor directory", bed_path);
  }

  vector<Hit> hits;
  map<string, FILE *> tmp_bins;
  map<string, int> lens;
  int ix = 0, total_nhits = 0;
  for (auto &file : files) {
    ifstream fin(file.c_str());
    if (!fin.is_open()) {
      throw fmt::format("BED file {} does not exist", bed_path);
    }

    int nhits = 0;
    string s;
    while (getline(fin, s)) {
      Hit h = Hit::from_bed(s);
      if (extend) {
        h.extend(Globals::Extend::RATIO, Globals::Extend::MAX_EXTEND);
      }
      assert(h.ref != nullptr);
      assert(h.query != nullptr);
      if (tie(h.query->name, h.query_start, h.query_end) >
          tie(h.ref->name, h.ref_start, h.ref_end)) {
        swap(h.query->name, h.ref->name);
        swap(h.query_start, h.ref_start);
        swap(h.query_end, h.ref_end);
      }
      // Bucket reads to allow external sorting
      string fno =
          output_dir + fmt::format("/tmp_{}_{}.tmp", lookup[h.query->name],
                                   lookup[h.ref->name]);
      auto it = tmp_bins.find(fno);
      if (it == tmp_bins.end()) {
        tmp_bins[fno] = fopen(fno.c_str(), "w");
        it = tmp_bins.find(fno);
        assert(it != tmp_bins.end());
      }
      FILE *fo = it->second;
      fputs(h.to_bed(false).c_str(), fo);
      fputs("\n", fo);
      lens[fno]++;

      // hits.push_back(h);
      nhits++;
      total_nhits++;
    }
    eprnn("\rRead {:10} alignments in {}         ", nhits, file);
  }
  eprn("\nRead total {} alignments", total_nhits);

  int max_complexity = 0;
  map<int, int> complexity;
  for (auto bin : tmp_bins) {
    fclose(bin.second);
    ifstream fin(bin.first.c_str());
    eprnn("\rProcessing bucket {}...", bin.first);
    vector<Hit> hits;
    hits.reserve(lens[bin.first]);
    string s;
    while (getline(fin, s)) {
      Hit h = Hit::from_bed(s);
      hits.push_back(h);
    }
    fin.close();

    if (extend) {
      hits = merge(hits, Globals::Extend::MERGE_DIST);
      eprnn("after merging remaining {} alignments       ", hits.size());
    }
    for (auto &h : hits) {
      int c = (int)sqrt(double(h.query_end - h.query_start) *
                        double(h.ref_end - h.ref_start));
      max_complexity = max(max_complexity, c);
      complexity[c / 1000]++;
    }

    FILE *fo = fopen(bin.first.c_str(), "w");
    for (auto &h : hits) {
      fputs(h.to_bed(false).c_str(), fo);
      fputs("\n", fo);
    }
    fclose(fo);
  }
  eprn("\nFinished with sorting");

  auto next_bin = vector<int>(1, 0);
  for (int c = 1; c <= max_complexity / 1000; c++)
    next_bin.push_back((next_bin[c - 1] + complexity[c - 1]) % nbins);

  const int BUFF_SZ = 1000;
  vector<vector<Hit>> buffer(nbins);

  vector<ofstream> fout;
  for (int b = 0; b < nbins; b++) {
    string of = output_dir + fmt::format("/bucket_{:04d}", b);
    fout.push_back(ofstream(of.c_str()));
    if (!fout.back().is_open()) {
      throw fmt::format("Cannot open file {} for writing", of);
    }
  }

  for (auto bin : tmp_bins) {
    ifstream fin(bin.first.c_str());
    eprnn("\rProcessing bucket {}...", bin.first);
    string s;
    while (getline(fin, s)) {
      Hit h = Hit::from_bed(s);
      int complexity = sqrt(double(h.query_end - h.query_start) *
                            double(h.ref_end - h.ref_start));
      complexity /= 1000;

      int bin = next_bin[complexity];
      next_bin[complexity] = (next_bin[complexity] + 1) % nbins;

      if (h.query->is_rc) {
        swap(h.query, h.ref);
        swap(h.query_start, h.ref_start);
        swap(h.query_end, h.ref_end);
      }
      buffer[bin].push_back(h);
      if (buffer[bin].size() == BUFF_SZ) {
        for (auto &h : buffer[bin])
          fout[bin] << h.to_bed(false) << endl;
        buffer[bin].clear();
      }
    }
    fin.close();
  }

  for (int b = 0; b < nbins; b++) {
    for (auto &h : buffer[b])
      fout[b] << h.to_bed(false) << endl;
    fout[b].close();
  }
  for (auto &s : tmp_bins)
    unlink(s.first.c_str());
}

auto bucket_alignments(const string &bed_path, int nbins, string output_dir,
                       bool extend) {
  vector<string> files;
  if (S_ISREG(stat_file(bed_path))) {
    files.push_back(bed_path);
  } else if (S_ISDIR(stat_file(bed_path))) {
    glob_t glob_result;
    glob((bed_path + "/*.bed").c_str(), GLOB_TILDE, NULL, &glob_result);
    for (int i = 0; i < glob_result.gl_pathc; i++) {
      string f = glob_result.gl_pathv[i];
      if (S_ISREG(stat_file(f))) {
        files.push_back(f);
      }
    }
  } else {
    throw fmt::format("Path {} is neither file nor directory", bed_path);
  }

  vector<Hit> hits;
  for (auto &file : files) {
    ifstream fin(file.c_str());
    if (!fin.is_open()) {
      throw fmt::format("BED file {} does not exist", bed_path);
    }

    int nhits = 0;
    string s;
    while (getline(fin, s)) {
      Hit h = Hit::from_bed(s);
      if (extend) {
        h.extend(Globals::Extend::RATIO, Globals::Extend::MAX_EXTEND);
      }
      hits.push_back(h);
      nhits++;
    }
    eprn("Read {} alignments in {}", nhits, file);
  }

  eprn("Read total {} alignments", hits.size());
  if (extend) {
    hits = merge(hits, Globals::Extend::MERGE_DIST);
    eprn("After merging remaining {} alignments", hits.size());
  }
  int max_complexity = 0;
  for (auto &h : hits) {
    max_complexity =
        max(max_complexity, (int)sqrt(double(h.query_end - h.query_start) *
                                      double(h.ref_end - h.ref_start)));
  }
  vector<vector<Hit>> bins(max_complexity / 1000 + 1);
  for (auto &h : hits) {
    int complexity = sqrt(double(h.query_end - h.query_start) *
                          double(h.ref_end - h.ref_start));
    assert(complexity / 1000 < bins.size());
    bins[complexity / 1000].push_back(h);
  }

  vector<vector<Hit>> results(nbins);
  int bc = 0;
  for (auto &bin : bins) {
    for (auto &hit : bin) {
      results[bc].push_back(hit);
      bc = (bc + 1) % nbins;
    }
  }

  if (output_dir != "") {
    int count = 0;
    for (auto &bin : results) {
      string of = output_dir + fmt::format("/bucket_{:04d}", count++);
      ofstream fout(of.c_str());
      if (!fout.is_open()) {
        throw fmt::format("Cannot open file {} for writing", of);
      }
      for (auto &h : bin) {
        fout << h.to_bed(false) << endl;
      }
      fout.close();
      eprn("Wrote {} alignments in {}", bin.size(), of);
    }
  }

  return results;
}

void generate_alignments(const string &ref_path, const string &bed_path,
                         int kmer_size, string ref_path2 = "") {
  auto T = cur_time();

  auto schedule = bucket_alignments(bed_path, 1, "", false);
  if (ref_path2 == "")
  {
    ref_path2 = ref_path;
  }
  FastaReference fr(ref_path);
  FastaReference fr2(ref_path2);


  int lines = 0, total = 0;
  for (auto &s : schedule)
    total += s.size();

  eprn("Using k-mer size {}", kmer_size);

  int total_written = 0;
  for (int i = 0; i < schedule.size(); i++) {
    for (auto &h : schedule[i]) {
      lines++;
      // auto name1 = h.query->name
      
      string delimiter = "#";
      string name1 = h.query->name.substr(h.query->name.find_last_of(delimiter) + 1, h.query->name.size());  //   h.query->name.substr(0, h.query->name.find(delimiter));
      string name2 = h.ref->name.substr(h.ref->name.find_last_of(delimiter) + 1, h.ref->name.size());

      string fa = fr.get_sequence(name1, h.query_start, &h.query_end);
      string fb = fr2.get_sequence(name2, h.ref_start, &h.ref_end);
      if (h.ref->is_rc)
        fb = rc(fb);

      eprnn("\r Processing {} out of {} ({:.1f}%, len {:10n} to {:10n})", lines,
            total, pct(lines, total), fa.size(), fb.size());

      // eprn("{}", h.to_bed(0));

      auto alns = fast_align(fa, fb, h, kmer_size);
      for (auto &hh : alns) {
        hh.query_start += h.query_start;
        hh.query_end += h.query_start;
        if (h.ref->is_rc) {
          swap(hh.ref_start, hh.ref_end);
          hh.ref_start = h.ref_end - hh.ref_start;
          hh.ref_end = h.ref_end - hh.ref_end;
          hh.ref->is_rc = true;
        } else {
          hh.ref_start += h.ref_start;
          hh.ref_end += h.ref_start;
        }
        hh.query->name = h.query->name;
        hh.ref->name = h.ref->name;
        hh.ref->name = h.ref->name;
        total_written++;
        prn("{}\t{}", hh.to_bed(false), h.to_bed(0));
      }
    }
  }

  eprn("\nFinished BED {} in {}s ({} lines, generated {} hits)", bed_path,
       elapsed(T), lines, total_written);
}

/******************************************************************************/

void align_main(int argc, char **argv) {
  using namespace Globals;
  argh::parser cmdl(argc, argv, argh::parser::PREFER_PARAM_FOR_UNREG_OPTION);

  cmdl("match", Align::MATCH) >> Align::MATCH;
  cmdl("mismatch", Align::MISMATCH) >> Align::MISMATCH;
  cmdl("gap-open", Align::GAP_OPEN) >> Align::GAP_OPEN;
  cmdl("gap-extend", Align::GAP_EXTEND) >> Align::GAP_EXTEND;

  cmdl("extend-ratio", Extend::RATIO) >> Extend::RATIO;
  cmdl("max-extend", Extend::MAX_EXTEND) >> Extend::MAX_EXTEND;
  cmdl("merge-dist", Extend::MERGE_DIST) >> Extend::MERGE_DIST;

  if (!cmdl(2)) {
    throw fmt::format("Not enough arguments to align");
  }

  string command = cmdl[0];
  if (command == "bucket") {
    int nbins;
    if (!(cmdl({"-n", "--bins"}) >> nbins)) {
      throw fmt::format("Must provide number of bins (--bins)");
    }
    bucket_alignments_extern(cmdl[1], nbins, cmdl[2], true, cmdl[3]);
  } else if (command == "generate") {
    int kmer_size;
    if (!(cmdl({"-k", "--kmer"}) >> kmer_size)) {
      throw fmt::format("Must provide k-mer size (--kmer)");
    }
    auto ref_str = cmdl[1];
    auto query_str = cmdl[1];

    if (cmdl(3)) {
      query_str = cmdl[3];
      // eprn("HERE {}, {}", query_str, ref_str);

    }
    generate_alignments(ref_str, cmdl[2], kmer_size, query_str);
  } else {
    throw fmt::format("Unknown align command");
  }
}
