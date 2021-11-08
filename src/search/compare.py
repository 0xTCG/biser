#%%
import sys
import pyranges as pr
import pandas as pd
import ncls

chrs, start, end = [], [], []
with open(sys.argv[1]) as f:
    for l in f:
        l = l.split()
        chrs.append(l[0] + l[7])
        start.append(int(l[1]))
        end.append(int(l[2]))
        chrs.append(l[3] + l[8])
        start.append(int(l[4]))
        end.append(int(l[5]))
rx = pr.PyRanges(chromosomes=chrs, starts=start, ends=end)

chrs, start, end = [], [], []
with open(sys.argv[2]) as f:
    for l in f:
        l = l.split()
        chrs.append(l[0] + l[7])
        start.append(int(l[1]))
        end.append(int(l[2]))
        chrs.append(l[3] + l[8])
        start.append(int(l[4]))
        end.append(int(l[5]))
ry = pr.PyRanges(chromosomes=chrs, starts=start, ends=end)
print('done')

d = pr.count_overlaps({"x": rx, "y": ry})
print('done2')

xuniq, yuniq, shared = 0, 0, 0
for _, r in d.df.iterrows():
    span = r["End"] - r["Start"]
    if r["x"] == 0 and r["y"] == 0:
        continue
    elif r["x"] == 0:
        yuniq += span
    elif r["y"] == 0:
        xuniq += span
    else:
        shared += span
print(f"x={xuniq :,}, y={yuniq :,}, xy={shared :,}")

#%%
# pd.set_option("display.max_rows", None)
# f = pr.count_overlaps({"x": rx, "y": ry}).df
# q = f[(f["y"] == 0) & (f["x"] != 0)]
# q["dst"] = q["End"] - q["Start"]
# q.sort_values(by="dst")
