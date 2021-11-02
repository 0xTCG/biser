#%%
import pyranges as pr
import pandas as pd

chrs, start, end = [], [], []
with open("_x") as f:
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
with open("_y") as f:
    for l in f:
        l = l.split()
        chrs.append(l[0] + l[7])
        start.append(int(l[1]))
        end.append(int(l[2]))
        chrs.append(l[3] + l[8])
        start.append(int(l[4]))
        end.append(int(l[5]))
ry = pr.PyRanges(chromosomes=chrs, starts=start, ends=end)

xuniq, yuniq, shared = 0, 0, 0
for _, r in pr.count_overlaps({"x": rx, "y": ry}).df.iterrows():
    span = r["End"] - r["Start"]
    if r["x"] == 0 and r["y"] == 0:
        continue
    elif r["x"] == 0:
        yuniq += span
    elif r["y"] == 0:
        xuniq += span
    else:
        shared += span
print(f"x={xuniq / 1e6:,.2f}, y={yuniq / 1e6:,.2f}, xy={shared / 1e6:,.2f}")

#%%
# pd.set_option("display.max_rows", None)
# f = pr.count_overlaps({"x": rx, "y": ry}).df
# q = f[(f["y"] == 0) & (f["x"] != 0)]
# q["dst"] = q["End"] - q["Start"]
# q.sort_values(by="dst")
