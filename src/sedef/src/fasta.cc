/// 786

/// This file is subject to the terms and conditions defined in
/// file 'LICENSE', which is part of this source code package.

/// Author: inumanag
/// Adapted from Bedtools source code
/// https://github.com/arq5x/bedtools/tree/master/src/fastaFromBed

/******************************************************************************/

#include <algorithm>
#include <fstream>
#include <string>
#include <sys/mman.h>
#include <sys/stat.h>

#include "common.h"
#include "fasta.h"

using namespace std;

/******************************************************************************/

FastaIndex::FastaIndex(const string &fname) {
  string line;
  long long linenum = 0;
  ifstream indexFile;
  indexFile.open(fname.c_str(), ifstream::in);
  if (indexFile.is_open()) {
    while (getline(indexFile, line)) {
      ++linenum;
      // the fai format defined in samtools is tab-delimited, every line being:
      // fai->name[i], (int)x.len, (long long)x.offset, (int)x.line_blen,
      // (int)x.line_len
      vector<string> fields = split(line, '\t');
      if (fields.size() == 5) { // if we don't get enough fields then there is a
                                // problem with the file
        // note that fields[0] is the sequence name
        char *end;
        string name = split(fields[0], ' ').at(0); // key by first token of name
        sequenceNames.push_back(name);
        this->insert(make_pair(
            name,
            FastaIndexEntry(fields[0], atoi(fields[1].c_str()),
                            strtoll(fields[2].c_str(), &end, 10),
                            atoi(fields[3].c_str()), atoi(fields[4].c_str()))));
      } else {
        throw fmt::format("Index file {} is malformed at line {}", fname,
                          linenum);
      }
    }
    indexFile.close();
  } else {
    throw fmt::format("Index file {} does not exist", fname);
  }
}

/******************************************************************************/

FastaIndexEntry FastaIndex::entry(const string &name) {
  FastaIndex::iterator e = this->find(name);
  if (e == this->end()) {
    throw fmt::format("Chromosome {} does not exist", name);
  } else {
    return e->second;
  }
}

/******************************************************************************/

FastaReference::FastaReference(string filename) {
  if (!(file = fopen(filename.c_str(), "r"))) {
    throw fmt::format("Cannot open file {}", filename);
  }

  struct stat stFileInfo;
  string indexFileName = filename + ".fai";
  // if we can find an index file, use it
  if (stat(indexFileName.c_str(), &stFileInfo) == 0) {
    // check if the index file is older than the FASTA file
    struct stat index_attrib, fasta_attrib;
    stat(indexFileName.c_str(), &index_attrib);
    stat(filename.c_str(), &fasta_attrib);
    if (fasta_attrib.st_mtime > index_attrib.st_mtime) {
      eprn("Warning: the index file is older than the FASTA file");
    }
    index = FastaIndex(indexFileName);
  }
  int fd = fileno(file);
  struct stat sb;
  if (fstat(fd, &sb) == -1) {
    throw fmt::format("Cannot stat file {}", filename);
  }
  filesize = sb.st_size;
  // map the whole file
  filemm = mmap(NULL, filesize, PROT_READ, MAP_SHARED, fd, 0);
}

FastaReference::~FastaReference(void) {
  fclose(file);
  munmap(filemm, filesize);
}

string FastaReference::get_sequence(string seqname, int start, int *end) {
  FastaIndexEntry entry = index.entry(seqname);

  if (start < 0) {
    start = 0;
  }
  int length;
  if (end == nullptr || *end > entry.length) {
    length = entry.length - start;
    if (end != nullptr) {
      *end = entry.length;
    }
  } else {
    length = *end - start;
  }

  // we have to handle newlines
  // approach: count newlines before start
  //           count newlines by end of read
  //             subtracting newlines before start find count of embedded
  //             newlines
  int newlines_before = start > 0 ? (start - 1) / entry.line_blen : 0;
  int newlines_by_end = (start + length - 1) / entry.line_blen;
  int newlines_inside = newlines_by_end - newlines_before;
  int seqlen = length + newlines_inside;
  char *seq = (char *)calloc(seqlen + 1, sizeof(char));

  memcpy(seq, (char *)filemm + entry.offset + newlines_before + start, seqlen);
  seq[seqlen] = '\0';
  char *pbegin = seq;
  char *pend = seq + (seqlen / sizeof(char));
  pend = std::remove(pbegin, pend, '\n');
  pend = std::remove(pbegin, pend, '\0');
  string s = seq;
  free(seq);
  s.resize((pend - pbegin) / sizeof(char));
  return s;
}
