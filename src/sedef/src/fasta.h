/// 786

/// This file is subject to the terms and conditions defined in
/// file 'LICENSE', which is part of this source code package.

/// Author: inumanag
/// Adapted from Bedtools source code
/// https://github.com/arq5x/bedtools/tree/master/src/fastaFromBed

/******************************************************************************/

#pragma once

/******************************************************************************/

#include <cstdio>
#include <map>
#include <string>
#include <vector>

/******************************************************************************/

struct FastaIndexEntry {
  std::string name;
  int length;
  long long offset;
  int line_blen, line_len;

  FastaIndexEntry(std::string name, int length, long long offset, int line_blen,
                  int line_len)
      : name(name), length(length), offset(offset), line_blen(line_blen),
        line_len(line_len) {}
};

/******************************************************************************/

struct FastaIndex : public std::map<std::string, FastaIndexEntry> {
  std::vector<std::string> sequenceNames;

  FastaIndex(){};
  FastaIndex(const std::string &fname);
  FastaIndexEntry entry(const std::string &key);
};

/******************************************************************************/

struct FastaReference {
  FILE *file;
  void *filemm;
  size_t filesize;
  FastaIndex index;

  std::map<std::string, std::vector<std::pair<size_t, std::string>>>
      translation_index;

  FastaReference(std::string filename);
  ~FastaReference();
  std::string get_sequence(std::string seqname, int start = 0, int *end = 0);
};
