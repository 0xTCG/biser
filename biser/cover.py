import sys
import ncls


def greedy_set_cover(subsets, parent_set):
  # https://stackoverflow.com/questions/21973126/set-cover-or-hitting-set-numpy-least-element-combinations-to-make-up-full-set

  import heapq

  parent_set = set(parent_set)
  max = len(parent_set)
  # create the initial heap. Note 'subsets' can be unsorted,
  # so this is independent of whether remove_redunant_subsets is used.
  heap = []
  for se, s in subsets.items():
    # Python's heapq lets you pop the *smallest* value, so we
    # want to use max-len(s) as a score, not len(s).
    # len(heap) is just proving a unique number to each subset,
    # used to tiebreak equal scores.
    heapq.heappush(heap, [max-len(s), se, s])
  results = []
  result_set = set()
  while result_set < parent_set:
    best = []
    unused = []
    while heap:
      score, count, s = heapq.heappop(heap)
      if not best:
        best = [max-len(s - result_set), count, s]
        continue
      if score >= best[0]:
        # because subset scores only get worse as the resultset
        # gets bigger, we know that the rest of the heap cannot beat
        # the best score. So push the subset back on the heap, and
        # stop this iteration.
        heapq.heappush(heap, [score, count, s])
        break
      score = max-len(s - result_set)
      if score >= best[0]:
        unused.append([score, count, s])
      else:
        unused.append(best)
        best = [score, count, s]
    add_set = best[2]
    results.append((best[1], add_set))
    result_set.update(add_set)
    while unused:
      heapq.heappush(heap, unused.pop())
  return results


def cover(bed, elems):
  trees = {}
  elem_sd = {}
  sd_elem = {}
  elem_ids = {}
  lines = []
  with open(elems) as f:
    for lp in f:
      l = lp.strip().split('\t')
      l[7] = l[7] if l[7] == '+' else '-'
      key = (l[0], l[1], l[7])
      if key not in trees:
        trees[key] = [], [], []
      trees[key][0].append(int(l[2]))
      trees[key][1].append(int(l[3]))
      if l[4] not in elem_ids:
        elem_ids[l[4]] = len(elem_ids)
      trees[key][2].append(elem_ids[l[4]])
      elem_sd[elem_ids[l[4]]] = set()
      lines.append(l)
      # Seq TODO: *key w/o parens does not work
      # ints[(*key, int(l[2]), int(l[3]))] = l[4]

  for i in trees:
    trees[i] = ncls.NCLS(*trees[i])

  with open(bed) as f:
    li = 0
    for lp in f:
      l = lp.strip().split('\t')
      sp = l[6].split(':')
      k1, i1 = (sp[0], l[0], l[8]), (int(l[1]), int(l[2]))
      k2, i2 = (sp[1], l[3], l[9]), (int(l[4]), int(l[5]))

      se = set()
      for k, i in ((k1, i1), (k2, i2)):
        if k in trees:
          # IntervalTree is broken (line 224: z.k is too large)
          for it in trees[k].find_overlap(*i):
            elem_sd[it[2]].add(li)
            se.add(it[2])
      if len(se) > 0:
        sd_elem[li] = se
        li += 1
    # print(f'{len(elem_sd)} elementary SDs')

    cores = set()
    for eid, sds in greedy_set_cover(elem_sd, sd_elem):
      if len(sds) < 2:
        continue
      cores.add(eid)
    # print(f'{len(cores)} core SDs')

    with open(elems, 'w') as fo:
      for l in lines:
        if elem_ids[l[4]] in cores and len(l) <= 8:
          l.append('CORE')
        print(*l, sep='\t', file=fo)
