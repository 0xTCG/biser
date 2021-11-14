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


@contextlib.contextmanager
def timing(title, results=None):
  t = time.time()
  yield
  t = int(time.time() - t)
  if results:
    st = int(sum(i for i, _ in results))
    print(f'{title}: {t//60:02d}:{t%60:02d}s (single: {st//60:02d}:{st%60:02d}s)')


def progress(*args, **kwargs):
  return tqdm.tqdm(
    *args, **kwargs, bar_format='{l_bar}{bar:40}| {n_fmt}/{total_fmt}'
  )


def run_biser(*args):
  path = f'{os.path.dirname(__file__)}/seq/biser'
  o = sp.run(
    [path, *args],
    env={"LD_LIBRARY_PATH": "/Users/inumanag/Projekti/seq/devel/build_release"},
    capture_output=True
  )
  if o.returncode != 0:
    print(o.stdout.decode('ascii'))
    print(o.stderr.decode('ascii'))
  return o


def search(tmp, path, chr):
  out = f"{tmp}/search_{chr}.bed"
  # out = f"{tmp}/search"
  run_biser("search", path, "-c", chr, "-o", out)
  return out


def align(path, file):
  out = f"{file}_align"
  run_biser("align", path, file, '-o', out)
  return out


def decompose(file):
  out = f"{file}_dec"
  run_biser("decompose", file, '-o', out)
  return out


def star(args):
  t = time.time()
  o = None
  if args[0] == 'search':
    o = search(*args[1:])
  if args[0] == 'align':
    o = align(*args[1:])
  if args[0] == 'decompose':
    o = decompose(*args[1:])
  return time.time() - t, o


if __name__ == '__main__':
  path = sys.argv[1]
  threads = int(sys.argv[2])
  final = sys.argv[3]
  nbuckets = 100
  tmp = sys.argv[4]
  # with tempfile.TemporaryDirectory(prefix='biser') as tmp:

  print(f'Temp: {tmp} on {threads} threads')
  results = []
  with open(f'{path}.fai') as f:
    chrs = [l.split()[0] for l in f]
    chrs = [('search', tmp, path, c) for c in chrs if '_' not in c and c != 'chrM']
  with timing('Search', results), mp.Pool(threads) as pool:
    results[:] = list(progress(pool.imap(star, chrs), total=len(chrs)))

  with timing('Spread'):
    hits = []
    for _, r in results:
      with open(r) as f:
        for l in f:
          l = l.strip()
          s = l.split()
          span = max(int(s[2]) - int(s[1]), int(s[5]) - int(s[4]))
          hits.append((span, l))
    hits.sort()
    buckets = [[] for i in range(nbuckets)]
    for i, (_, h) in enumerate(hits):
      buckets[i % nbuckets].append(h)
    aligns = []
    for i, b in enumerate(buckets):
      out = f'{tmp}/bucket{i}.bed'
      with open(out, 'w') as fo:
        for l in b:
          print(l, file=fo)
      aligns.append(('align', path, out))

  with timing('Align', results), mp.Pool(threads) as pool:
    results[:] = list(progress(pool.imap(star, aligns), total=len(aligns)))
    with open(final, 'w') as fo:
      for _, r in results:
        with open(r) as f:
          for l in f:
            print(l, end='', file=fo)

  with timing('Decomposition', results):
    output = f'{tmp}/clusters'
    os.makedirs(output, exist_ok=True)
    run_biser('cluster', final, path, '-o', output)
    clusters = []
    for dirpath, _, files in os.walk(output):
      for f in files:
        if f.endswith('.fa'):
          clusters.append(('decompose', os.path.abspath(os.path.join(dirpath, f))))
    with mp.Pool(threads) as pool:
      results[:] = list(progress(pool.imap(star, clusters), total=len(clusters)))
    with open(f'{final}.elem.txt', 'w') as fo:
      for _, r in results:
        with open(r) as f:
          for l in f:
            print(l, end='', file=fo)


