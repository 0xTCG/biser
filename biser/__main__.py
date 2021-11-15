# 786

import multiprocessing as mp
import tqdm
import sys
import subprocess
import tempfile
import os
import glob
import time
import contextlib
from pathlib import Path


@contextlib.contextmanager
def timing(title, results=None, force=False):
  t = time.time()
  yield
  t = int(time.time() - t)
  if results:
    st = int(sum(i[0][0] for i in results if i[0]))
    print(f'{title}: {t//60:02d}:{t%60:02d}s (single: {st//60:02d}:{st%60:02d}s)')
  elif force:
    print(f'{title}: {t//60:02d}:{t%60:02d}s')


def progress(*args, **kwargs):
  return tqdm.tqdm(
    *args, **kwargs, bar_format='{l_bar}{bar:40}| {n_fmt}/{total_fmt}'
  )


def run_biser(*args):
  path = f'{os.path.dirname(__file__)}/seq/biser'
  t = time.time()
  o = subprocess.run(
    [path, *args],
    env={
      "LD_LIBRARY_PATH": "/Users/inumanag/Projekti/seq/devel/build_release",
      "OMP_NUM_THREADS": "1"},
    capture_output=True
  )
  t = time.time() - t
  l = ['[p] ' + ' '.join(o.args)]
  l += [f'[o] {l}' for l in o.stdout.decode('ascii').strip().split('\n')]
  l += [f'[e] {l}' for l in o.stderr.decode('ascii').strip().split('\n')]
  if o.returncode != 0:
    err = '\n'.join(l)
    raise RuntimeError(f"BISER failed:\n{err}")
  return t, l  # time, output

# def decompose(file):
#   out = f"{file}_dec"
#   run_biser("decompose", file, '-o', out)
#   return out

def biser_search(args):
  if len(args) == 4:
    tmp, genome, species, chr = args
    out = f"{tmp}/search/{species}_{chr}.bed"
    o = run_biser("search", genome, "-c", chr, "-o", out)
    return (o, species, chr, out)
  else:
    tmp, genome1, chr1, genome2, intv2, species1, species2 = args
    out = f"{tmp}/cross_search/{species1}_{chr1}_{species2}.bed"
    o = run_biser("search", genome1, chr1, genome2, intv2, "-o", out)
    return (o, species1, chr1, species2, out)


def search(tmp, genomes, threads): # genomes is {sp = Path(g).stem -> PATH}
  os.makedirs(f'{tmp}/search', exist_ok=True)

  jobs = []
  for sp, genome in genomes.items():
    with open(f'{genome}.fai') as f:
      chrs = [l.split()[0] for l in f]
      jobs += [(tmp, genome, sp, c) for c in chrs if '_' not in c and c != 'chrM']

  results = []
  with timing('Search', results), mp.Pool(threads) as pool:
    results[:] = list(progress(pool.imap(biser_search, jobs), total=len(jobs)))

  with open(f'{tmp}/search.log', 'w') as fo:
    for (_, log), sp, ch, _ in results:
      for l in log:
        print(f'{sp}.{ch}: {l}', file=fo)
  return results


def biser_align(args):
  species, bed, genomes = args
  out = f"{bed}.align"
  o = run_biser("align", *genomes, bed, '-o', out)
  return (o, species, bed, out)


def align(tmp, genomes, threads, search, nbuckets=50):
  os.makedirs(f'{tmp}/align', exist_ok=True)

  jobs = []
  with timing('Spread'):
    hits = {sp: [] for sp in genomes}
    for _, sp, _, out in search:
      with open(out) as f:
        for l in f:
          l = l.strip()
          s = l.split()
          span = max(int(s[2]) - int(s[1]), int(s[5]) - int(s[4]))
          hits[sp].append((span, l))
    for h in hits.values():
      h.sort()
    buckets = {sp: [[] for i in range(nbuckets)] for sp in genomes}
    for sp, hs in hits.items():
      for i, (_, h) in enumerate(hs):
        buckets[sp][i % nbuckets].append(h)

    for sp in buckets:
      for i, b in enumerate(buckets[sp]):
        bed = f'{tmp}/align/{sp}.{i:05d}.bed'
        with open(bed, 'w') as fo:
          for l in b:
            print(l, file=fo)
        jobs.append((sp, bed, [genomes[sp]]))

  results = []
  with timing('Align', results), mp.Pool(threads) as pool:
    results[:] = list(progress(pool.imap(biser_align, jobs), total=len(jobs)))

  with open(f'{tmp}/align.log', 'w') as fo:
    for (_, log), sp, bed, _ in results:
      for l in log:
        print(f'{sp}.{bed.split(".")[-2]}: {l}', file=fo)
  return results


if __name__ == '__main__':
  tmp = sys.argv[1]
  os.makedirs(tmp, exist_ok=True)

  threads = int(sys.argv[2])
  final = sys.argv[3]
  genomes = {Path(path).stem: path for path in sys.argv[4:]}

  print('Genomes:', list(genomes))
  # with timing("Total", force=True):
    # s = search(tmp, genomes, threads)
    # a = align(tmp, genomes, threads, s)

  # with tempfile.TemporaryDirectory(prefix='biser') as tmp:

  a = []
  for g in glob.glob(f'{tmp}/align/*.bed'):
    a.append((None, Path(g).stem.split('.')[0], g, f'{g}.align'))

  with timing('Cross-search'):
    print(len(a), 'results')
    os.makedirs(f'{tmp}/cross_search', exist_ok=True)
    beds = {}
    for _, sp, _, out in a:
      beds.setdefault(sp, []).append(out)
    for sp in beds:
      sp_bed = f'{tmp}/cross_search/{sp}.bed'
      with open(sp_bed, 'w') as fo:
        for b in beds[sp]:
          with open(b) as f:
            for l in f:
              print(l.strip(), file=fo)
      run_biser('extract', sp_bed, '-o', f'{sp_bed}.regions.txt')

    chrs = {}
    for sp, genome in genomes.items():
      chrs[sp] = []
      with open(f'{genome}.fai') as f:
        chrs[sp] = [l.split()[0] for l in f]
        chrs[sp] = [c for c in chrs[sp] if '_' not in c and c != 'chrM']

    jobs = [
      (tmp, genomes[g1], c, genomes[g2], f'{tmp}/cross_search/{g2}.bed.regions.txt', g1, g2)
      for g1 in genomes for g2 in genomes if g1 < g2 for c in chrs[g1]
    ]
    results = []
    with timing('Search', results), mp.Pool(threads) as pool:
      results[:] = list(progress(pool.imap(biser_search, jobs), total=len(jobs)))

  with timing('Cross-align'):
    os.makedirs(f'{tmp}/cross_align', exist_ok=True)

    jobs = []
    hits = {}
    for _, sp1, ch1, sp2, out in results:
      with open(out) as f:
        for l in f:
          l = l.strip()
          s = l.split()
          span = max(int(s[2]) - int(s[1]), int(s[5]) - int(s[4]))
          hits.setdefault((sp1, sp2), []).append((span, l))
    for h in hits.values():
      h.sort()
    nbuckets = 50
    buckets = {sp: [[] for i in range(nbuckets)] for sp in hits}
    for sp, hs in hits.items():
      for i, (_, h) in enumerate(hs):
        buckets[sp][i % nbuckets].append(h)

    for sp1, sp2 in buckets:
      for i, b in enumerate(buckets[sp1, sp2]):
        bed = f'{tmp}/cross_align/{sp1}.{sp2}.{i:05d}.bed'
        with open(bed, 'w') as fo:
          for l in b:
            print(l, file=fo)
        jobs.append(((sp1, sp2), bed, [genomes[sp1], genomes[sp2]]))

    results = []
    with timing('Align', results), mp.Pool(threads) as pool:
      results[:] = list(progress(pool.imap(biser_align, jobs), total=len(jobs)))

    with open(f'{tmp}/cross-align.log', 'w') as fo:
      for (_, log), (sp1, sp2), bed, _ in results:
        for l in log:
          print(f'{sp1}.{sp2}.{bed.split(".")[-2]}: {l}', file=fo)

  #   with open(final, 'w') as fo:
  #     for _, r in results:
  #       with open(r) as f:
  #         for l in f:
  #           print(l, end='', file=fo)

  # with timing('Decomposition', results):
  #   output = f'{tmp}/clusters'
  #   os.makedirs(output, exist_ok=True)
  #   run_biser('cluster', final, path, '-o', output)
  #   clusters = []
  #   for dirpath, _, files in os.walk(output):
  #     for f in files:
  #       if f.endswith('.fa'):
  #         clusters.append(('decompose', os.path.abspath(os.path.join(dirpath, f))))
  #   with mp.Pool(threads) as pool:
  #     results[:] = list(progress(pool.imap(star, clusters), total=len(clusters)))
  #   with open(f'{final}.elem.txt', 'w') as fo:
  #     for _, r in results:
  #       with open(r) as f:
  #         for l in f:
  #           print(l, end='', file=fo)


