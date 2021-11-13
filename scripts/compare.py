#%%
import sys
import itertools
import pathlib
import pyranges as pr
import pandas as pd
import tabulate

# sys.argv[1:] = ["../_a", "../_data/_hg19_search_2"]
ranges = {}
for fn in sys.argv[1:]:
  d = {"Chromosome": [],  "Start": [], "End": [], "Line": []}
  stm = pathlib.Path(fn).stem
  with open(fn) as f:
    for li, l in enumerate(f):
      if 'sdquest'in stm:
        l = l.split()
        l = l[1:4] + l[5:8] + [l[4], l[8]]
      elif 'wgac' in stm:
        l = l.split('\t')
        l = l[:3] + l[6:9] + ['+', l[5]]
      else:
        l = l.split('\t')
        l = l[:6] + l[8:10]
      for i in range(6, 8):
        if l[i] == '_': l[i] = '-'
      c1, s1, e1, st1 = l[0].split('#')[-1], int(l[1]), int(l[2]), l[6]
      c2, s2, e2, st2 = l[3].split('#')[-1], int(l[4]), int(l[5]), l[7]
      if (c1, s1, e1) >= (c2, s2, e2):
        c1, s1, e1, c2, s2, e2 = c2, s2, e2, c1, s1, e1
      d["Chromosome"] += [c1, c2]
      d["Start"] += [s1, s2]
      d["End"] += [e1, e2]
      # d["Strand"] += [st1, st2]
      d["Line"] += [f'{stm}.{li}'] * 2
      # if st1 != '+' or st2 != '+':
      #   cs += [c1, c2]
      #   st += [s1, s2]
      #   ed += [e1, e2]
      #   sd += [st2, st1]
  ranges[stm] = pr.from_dict(d)


#%%
ovl = pr.count_overlaps(ranges)
do = ovl.df
do["Size"] = do["End"] - do["Start"]
do

#%%
totals = {s: 0 for s in ranges}
rows = []
for i in range(len(ranges)):
  for c in itertools.combinations(ranges.keys(), i + 1):
    rule = do["Size"] > 0
    for r in ranges:
      rule &= (do[r] > 0) if r in c else (do[r] == 0)
    sz = do[rule]["Size"].sum()
    for r in c:
      totals[r] += sz
    rows.append(['*' if r in c else ' ' for r in ranges] + [sz])
    # name = ", ".join(c)
    # print(f'{name:50} -> {sz:15,}')
print(tabulate.tabulate(
  sorted(rows, key=lambda x: ''.join(x[:-1][::-1])),
  headers=[r.split('.')[-1][:3] for r in ranges] + ['Coverage'],
  floatfmt=',.0f'
))
print('='* 38)
for t, sz in totals.items():
  print(f'{t:20} -> {sz:14,}')

# pd.set_option("display.max_rows", None)
# print(do[(do._a > 0) & (do._hg19_search_2 == 0)].sort_values(by='Size').tail())

#%%
# f = pr.count_overlaps({"x": rx, "y": ry}).df
# q = f[(f["y"] == 0) & (f["x"] != 0)]
# q["dst"] = q["End"] - q["Start"]
# q.sort_values(by="dst")
