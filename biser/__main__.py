# 786

import multiprocessing as mp
import tqdm
import sys
import subprocess as sp
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
    st = int(sum(i[0][0] for i in results))
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
  o = sp.run(
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
  tmp, genome, species, chr = args
  out = f"{tmp}/search/{species}_{chr}.bed"
  o = run_biser("search", genome, "-c", chr, "-o", out)
  return (o, species, chr, out)


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
  species, bed, genome = args
  out = f"{bed}.align"
  o = run_biser("align", genome, bed, '-o', out)
  return (o, species, bed, out)


def align(tmp, genomes, threads, search, nbuckets=50):
  os.makedirs(f'{tmp}/search/align', exist_ok=True)

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
        bed = f'{tmp}/search/align/{sp}.{i:05d}.bed'
        with open(bed, 'w') as fo:
          for l in b:
            print(l, file=fo)
        jobs.append((sp, bed, genomes[sp]))

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
  with timing("Total", force=True):
    s = search(tmp, genomes, threads)
    a = align(tmp, genomes, threads, s)

  # with tempfile.TemporaryDirectory(prefix='biser') as tmp:



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


