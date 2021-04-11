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
    // if (db != "+") continue;
    // if (chra != "hg19#chr1") continue;
    string ca = chra + da, cb = chrb + db;
    if (!chrs[ca])
      chrs[ca] = chrs.size();
    if (!chrs[cb])
      chrs[cb] = chrs.size();
  }
  TST("load");

  struct ElemSD {
    int id;
    set<int> sds;
    vector<pair<int64_t, int>> copies;
  };
  
  vector<ElemSD> elems;
  vector<tuple<int64_t, int, int>> all_elems;
  map<pair<int64_t, int64_t>, int> elementary_ids;
  for (size_t i = 0, prev = -1; i < nodes.size(); i++)
    elems.push_back(ElemSD{i, {}, {}});
    
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
//   for (auto &b : blocks) {
//     auto it = lower_bound(all_elems.begin(), all_elems.end(),
//                           make_tuple(b.a - D, 0, 0));
//     for (; it != all_elems.end() && get<0>(*it) <= b.a + b.len; it++) {
//       int olp =
//           min((int64_t)min(get<1>(*it), b.len),
//               min(get<0>(*it) + get<1>(*it) - b.a, b.a + b.len - get<0>(*it)));
//       if (olp >= D / 2) {
//         auto eid = get<2>(*it);
//         sd2elem[b.sd_id].insert(eid);
//         elems[eid].sds.insert(b.sd_id);
//       }
//     }
//   }
  auto cores = set_cover(elems, sd2elem);
  TST("cores");

  return 0;
}
