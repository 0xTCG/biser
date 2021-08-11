#include <algorithm>
#include <boost/heap/fibonacci_heap.hpp>
#include <cassert>
#include <chrono>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <list>
#include <map>
#include <set>
#include <sstream>
#include <stdio.h>
#include <string>
#include <sys/resource.h>
#include <tuple>
#include <unistd.h>
#include <vector>
using namespace std;
using namespace std::chrono;

#define E(s, ...) fprintf(stderr, ". " s "\n", __VA_ARGS__)
#define TST(x)                                                                 \
  {                                                                            \
    fprintf(stderr, "[stage] " x ": %.1f (%.1f Mb)\n",                         \
            duration_cast<milliseconds>(high_resolution_clock::now() - T)      \
                    .count() /                                                 \
                1000.0,                                                        \
            getPeakRSS() / 1024.0 / 1024.0);                                   \
    T = high_resolution_clock::now();                                          \
  }

size_t getPeakRSS() {
  struct rusage rusage;
  getrusage(RUSAGE_SELF, &rusage);
  return (size_t)(rusage.ru_maxrss);
}

struct Node {
  Node *parent;
  int size;
  Node *find() {
    auto n = this;
    while (n->parent != n) {
      n->parent = n->parent->parent;
      n = n->parent;
    }
    return n;
  }
  Node *merge(Node *y) {
    auto x = this->find();
    assert(x && x->parent == x);
    y = y->find();
    assert(y && y->parent == y);
    if (x != y) {
      if (x->size < y->size)
        swap(x, y);
      y->parent = x;
      x->size += y->size;
    }
    return x;
  }
};

int main(int argc, char **argv) {
  ios_base::sync_with_stdio(false);
  auto T = high_resolution_clock::now();
  string path = argv[1];

  bool swap_ID = argc == 3;
  if (swap_ID)
    E("swapping I and D in CIGAR... %s", "");

  struct Block {
    int64_t a, b;
    int32_t len, sd_id;
    char edge;
  };
  vector<Block> blocks;
  ifstream fin(path);
  int sds = 0;
  map<string, int> chrs;
  for (string s; getline(fin, s); sds++) {
    istringstream ss(s);
    string chra, chrb, cigar, da, db, _;
    int64_t sa, ea, sb, eb;
    ss >> chra >> sa >> ea >> chrb >> sb >> eb >> _ >> da >> db >> _ >> _ >>
        cigar;
    if (swap_ID) {
      swap(chra, chrb);
      swap(da, db);
      swap(sa, sb);
      swap(ea, eb);
    }
  again:
    // if (db != "+") continue;
    // if (chra != "hg19#chr1") continue;
    string ca = chra + da, cb = chrb + db;
    if (!chrs[ca])
      chrs[ca] = chrs.size();
    if (!chrs[cb])
      chrs[cb] = chrs.size();
    int64_t cia = chrs[ca], cib = chrs[cb];
    auto ssa = sa, ssb = sb;
    for (int ci = 0, num = 0, co = 0; ci < cigar.size(); ci++) {
      if (isdigit(cigar[ci])) {
        num = 10 * num + cigar[ci] - '0';
      } else {
        if (cigar[ci] == 'M') {
          char edge = 0;
          if (sa == ssa || sb == ssb)
            edge |= 1;
          if (sa + num == ea || sb + num == eb)
            edge |= 2;
          blocks.push_back(
              {(cia << 32) + sa, (cib << 32) + sb, num, sds, edge});
          sa += num;
          sb += num;
        } else if (cigar[ci] == 'I') {
          sb += num;
        } else if (cigar[ci] == 'D') {
          sa += num;
        }
        co++;
        num = 0;
      }
    }
    // assert(sa == ea && sb == eb);
    if (db[0] != '+') {
      swap(da, db);
      sa = ssa, sb = ssb;
      sds++;
      goto again;
    }
  }
  TST("load");
  E("%d SDs, %lu blocks", sds, blocks.size());

  struct BlockInt {
    int64_t start;
    int len, block;
    bool operator<(const BlockInt &y) const {
      return make_pair(start, len) < make_pair(y.start, y.len);
    }
  };
  vector<BlockInt> block_intervals;
  vector<int64_t> pos;
  vector<char> is_edge;
  {
    for (auto &b : blocks) {
      block_intervals.push_back({b.a, b.len, int(&b - &blocks[0])});
      block_intervals.push_back({b.b, b.len, int(&b - &blocks[0])});
    }
    sort(block_intervals.begin(), block_intervals.end());
    int64_t last = 0;
    for (auto &bb : block_intervals) {
      for (last = max(bb.start, last); last <= bb.start + bb.len; last++)
        pos.push_back(last);
    }
    is_edge.resize(pos.size());
    for (auto &b : blocks) {
      if (b.edge & 1) {
        auto p = std::lower_bound(pos.begin(), pos.end(), b.a) - pos.begin();
        is_edge[p] |= 1;
        p = std::lower_bound(pos.begin(), pos.end(), b.b) - pos.begin();
        is_edge[p] |= 1;
      }
      if (b.edge & 2) {
        auto p = std::lower_bound(pos.begin(), pos.end(), b.a + b.len - 1) -
                 pos.begin();
        is_edge[p] |= 2;
        p = std::lower_bound(pos.begin(), pos.end(), b.b + b.len - 1) -
            pos.begin();
        is_edge[p] |= 2;
      }
    }
    for (auto &b : block_intervals)
      b.start = &(*lower_bound(pos.begin(), pos.end(), b.start)) - &pos[0];
    TST("pos");
    E("pos: %lu, bi = %lu", pos.size(), block_intervals.size());
  }
  vector<Node> nodes(pos.size());
  for (auto &n : nodes) {
    n.parent = &n;
    n.size = 0;
  }
  TST("nodes");

  auto progress = [](int c, int n) {
    if (c % (n / 1000) == 0)
      fprintf(stderr, "\r. progress: %.1f", double(c) * 100 / n);
  };

  const int D = 200;
  int rounded = 0;
  auto round = [&](int64_t x, vector<size_t> &activated) {
    auto pi = lower_bound(pos.begin(), pos.end(), x);
    assert(pi != pos.end() && x == *pi);
    size_t i = &(*pi) - &pos[0];
    for (auto ii = max(int64_t(0), int64_t(i) - D); ii < min(i + D, pos.size()); ii++)
      if (i != ii && labs(pos[ii] - x) < D && nodes[ii].size) {
        i = ii;
        rounded++;
        break;
      }
    if (!nodes[i].size) {
      nodes[i].size = 1;
      activated.push_back(i); // activate node
    }
    return nodes[i].find();
  };
  vector<size_t> activated_old;
  for (auto &b : blocks) {
    for (int i = 0; i < 2; i++) {
      int64_t x = !i ? b.a : b.a + b.len;
      int64_t y = !i ? b.b : b.b + b.len;
      auto nx = round(x, activated_old), ny = round(y, activated_old);
      nx->merge(ny);
    }
  }
  E("start with %lu dots, rounded %d dots", activated_old.size(), rounded);
  TST("stage_1");
  for (int iter = 2;; iter++) {
    rounded = 0;
    vector<size_t> activated_new;
    sort(activated_old.begin(), activated_old.end());
    for (auto &intv : block_intervals) {
      auto image = blocks[intv.block].a != pos[intv.start]
                       ? blocks[intv.block].a
                       : blocks[intv.block].b;
      auto i =
          lower_bound(activated_old.begin(), activated_old.end(), intv.start);
      for (; i != activated_old.end() && pos[*i] - pos[intv.start] <= intv.len;
           i++) {
        auto nx = round(pos[*i], activated_new),
             ny = round(image + pos[*i] - pos[intv.start], activated_new);
        nx->merge(ny);
      }
    }
    fprintf(stderr, "\r%3d, activated %8lu, rounded %8d ...", iter,
            activated_new.size(), rounded);
    if (!activated_new.size())
      break;
    else
      activated_old = activated_new;
  }
  E("%40s", "done");
  TST("stage_2");

  const int64_t mask = ((1llu << 32) - 1);
  map<int, string> rev_chrs;
  for (auto &c : chrs)
    rev_chrs[c.second] = c.first;


  struct ElemSD {
    int id;
    set<int> sds;
    vector<pair<int64_t, int>> copies;
  };
  vector<ElemSD> elems;
  vector<tuple<int64_t, int, int>> all_elems;
  map<pair<int64_t, int64_t>, int> elementary_ids;
  for (size_t i = 0, prev = -1; i < nodes.size(); i++)
    if (nodes[i].size) {
      if (prev != -1) {
        int64_t ca = pos[prev] >> 32, a = pos[prev] & mask;
        int64_t cb = pos[i] >> 32, b = pos[i] & mask;
        if (ca != cb)
          goto end;
        int64_t miss = 0;
        for (auto j = prev + 1; j <= i; j++)
          miss += pos[j] - pos[j - 1] - 1;
        if (b - a >= 1 && double(miss) / (b - a) < 0.1) {
          auto p = make_pair(nodes[prev].find() - &nodes[0],
                             nodes[i].find() - &nodes[0]);
          if (elementary_ids.find(p) == elementary_ids.end()) {
            int id = elems.size();
            elems.push_back(ElemSD{id, {}, {}});
            elementary_ids[p] = id;
          }
          int id = elementary_ids[p];
          elems[id].copies.push_back({pos[prev], pos[i] - pos[prev]});
          all_elems.push_back({pos[prev], pos[i] - pos[prev], id});
        }
      }
    end:
      prev = i;
    }
  sort(all_elems.begin(), all_elems.end());
  E("%lu elementary SDs, %lu copies", elems.size(), all_elems.size());

  for (auto &i : all_elems) {
    int64_t ca = get<0>(i) >> 32, a = get<0>(i) & mask;
    printf("%d\t%s\t%lld\t%lld\t%d\n", get<2>(i), rev_chrs[ca].c_str(), a,
           a + get<1>(i), get<1>(i));
  }
  TST("elementary");

  auto set_cover = [](vector<ElemSD> &elems, vector<set<int>> &sds) {
    boost::heap::fibonacci_heap<pair<int64_t, int>> heap;
    vector<boost::heap::fibonacci_heap<pair<int64_t, int>>::handle_type>
        heap_handles;
    for (auto &e : elems) {
      heap_handles.push_back(heap.push({e.sds.size(), e.id}));
    }
    // do the greedy algorithm
    vector<char> elems_visited(elems.size(), 0);
    vector<int> cover;
    while (!heap.empty()) {
      auto h = heap.top();
      heap.pop();
      printf("[core] %d (%lld SDs)\n", h.second, h.first);
      elems_visited[h.second] = true;
      if (h.first <= 1)
        break;
      cover.push_back(h.second);
      for (auto sd : elems[h.second].sds)
        for (auto elem : sds[sd])
          if (!elems_visited[elem]) {
            auto it = elems[elem].sds.find(sd);
            if (it != elems[elem].sds.end()) {
              elems[elem].sds.erase(it);
              (*heap_handles[elem]).first -= 1;
              heap.decrease(heap_handles[elem]);
            }
          }
    }
    return cover;
  };

  vector<set<int>> sd2elem(sds);
  for (auto &b : blocks) {
    auto it = lower_bound(all_elems.begin(), all_elems.end(),
                          make_tuple(b.a - D, 0, 0));
    for (; it != all_elems.end() && get<0>(*it) <= b.a + b.len; it++) {
      int olp =
          min((int64_t)min(get<1>(*it), b.len),
              min(get<0>(*it) + get<1>(*it) - b.a, b.a + b.len - get<0>(*it)));
      if (olp >= D / 2) {
        auto eid = get<2>(*it);
        sd2elem[b.sd_id].insert(eid);
        elems[eid].sds.insert(b.sd_id);
      }
    }
  }
  auto cores = set_cover(elems, sd2elem);
  TST("cores");

  return 0;
}
