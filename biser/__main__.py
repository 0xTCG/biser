# 786

import multiprocessing as mp
import tqdm
import sys
import subprocess
import tempfile
import argparse
import shutil
import os
import glob
import time
import contextlib
from pathlib import Path

from .cover import cover


@contextlib.contextmanager
def timing(title, results=None, force=False, indent=2):
  t = time.time()
  yield
  t = int(time.time() - t)
  ind = " " * indent
  if results:
    st = int(sum(i[0][0] for i in results if i[0]))
    print(f'{ind}{title}: {t//60:02d}:{t%60:02d}s (single: {st//60:02d}:{st%60:02d}s)')
  elif force:
    print(f'{ind}{title}: {t//60:02d}:{t%60:02d}s')


def progress(*args, **kwargs):
  return tqdm.tqdm(*args, **kwargs, bar_format='{l_bar}{bar:40}| {n_fmt}/{total_fmt}')


def run_biser(*args):
  root = os.path.dirname(__file__)
  root = '/home/vagrant/.pyenv/versions/3.7.13/lib/python3.7/site-packages/biser'
  path = f'{root}/exe/biser.exe'
  t = time.time()
  o = subprocess.run(
    [path, *args],
    env={"OMP_NUM_THREADS": "1"},
    stdout=subprocess.PIPE, 
    stderr=subprocess.PIPE,
  )
  t = time.time() - t
  l = ['[p] ' + ' '.join(o.args)]
  l += [f'[o] {l}' for l in o.stdout.decode('ascii').strip().split('\n')]
  l += [f'[e] {l}' for l in o.stderr.decode('ascii').strip().split('\n')]
  # print([path, *args], o.returncode)
  if o.returncode != 0:
    err = '\n'.join(l)
    raise RuntimeError(f"BISER failed:\n{err}")
  return t, l  # time, output


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


def search(tmp, genomes, threads, contigs = False):
  os.makedirs(f'{tmp}/search', exist_ok=True)
  print('1. Putative SD detection')

  jobs = []
  for sp, genome in genomes.items():
    with open(f'{genome}.fai') as f:
      chrs = [l.split()[0] for l in f]
      jobs += [
        (tmp, genome, sp, c) 
        for c in chrs 
        if contigs or ('_' not in c and c != 'chrM')
      ]

  results = []

  if jobs:
    with timing('Search', results), mp.Pool(threads) as pool:
      results[:] = list(progress(pool.imap(biser_search, jobs), total=len(jobs)))

    with open(f'{tmp}/search.log', 'w') as fo:
      for (_, log), sp, ch, _ in results:
        for l in log:
          print(f'{sp}.{ch}: {l}', file=fo)
  else:
    print(
      "   No chromosomes found. Try running with --keep-contigs "
      "if you are using scaffolds."
    )
  return results


def biser_align(args):
  species, bed, genomes = args
  out = f"{bed}.align"
  o = run_biser("align", *genomes, bed, '-o', out)
  return (o, species, bed, out)


def align(tmp, genomes, threads, search, nbuckets=50):
  os.makedirs(f'{tmp}/align', exist_ok=True)
  print('2. Putative SD alignment')

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
      os.remove(bed)
  return results


def biser_decompose(file):
  out = f"{file}_dec"
  o = run_biser("decompose", file, '-o', out)
  return o, out


def cross_biser(tmp, genomes, threads, alignments, contigs = False):
  with timing('Cross-search'):
    os.makedirs(f'{tmp}/cross_search', exist_ok=True)
    print('1. Cross-genome putative SD detection')

    beds = {}
    for _, sp, _, out in alignments:
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
        chrs[sp] = [
          c 
          for c in chrs[sp] 
          if contigs or ('_' not in c and c != 'chrM')
        ]

    jobs = [
      (tmp, genomes[g1], c, genomes[g2], f'{tmp}/cross_search/{g2}.bed.regions.txt', g1, g2)
      for g1 in genomes for g2 in genomes if g1 < g2 for c in chrs[g1]
    ]
    results = []
    with timing('Search', results), mp.Pool(threads) as pool:
      results[:] = list(progress(pool.imap(biser_search, jobs), total=len(jobs)))

  with timing('Cross-align'):
    os.makedirs(f'{tmp}/cross_align', exist_ok=True)
    print('2. Cross-genome putative SD alignment')

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
        os.remove(bed)

  return results


def decompose(tmp, genomes, threads, final):
  print('3. SD decomposition')

  results = []
  with timing('Decomposition', results):
    output = f'{tmp}/clusters'
    os.makedirs(output, exist_ok=True)
    o = run_biser('cluster', final, *list(genomes.values()), '-o', output)
    clusters = []
    for dirpath, _, files in os.walk(output):
      for f in files:
        if f.endswith('.fa'):
          clusters.append(os.path.abspath(os.path.join(dirpath, f)))
    if clusters:
      with mp.Pool(threads) as pool:
        results[:] = list(progress(pool.imap(biser_decompose, clusters), total=len(clusters)))
      with open(f'{tmp}/decompose.log', 'w') as fo:
        for l in o[1]:
          print(f'cluster: {l}', file=fo)
        for (_, log), o in results:
          for l in log:
            print(f'{o}: {l}', file=fo)
      with open(f'{final}.elem.txt', 'w') as fo:
        for _, r in results:
          with open(r) as f:
            for l in f:
              print(l, end='', file=fo)
      cover(final, f'{final}.elem.txt')
    else:
      print("   No SDs found to decompose.")

def biser_mask(args):
  tmp, species, genome = args
  out = f'{tmp}/genomes/{species}.fa'
  o = run_biser('mask', genome, "-o", out)
  subprocess.check_call(['samtools', 'faidx', out])
  return o, species, out


def main(argv):
  if len(argv) > 0 and argv[0] == 'test':
    _, l = run_biser('hello')
    if len(l) < 2 or l[1] != '[o] BISER v1.0':
      print('Error: unexpected respose from BISER')
      print('\n'.join(l))
      sys.exit(1)
    sys.exit(0)

  parser = argparse.ArgumentParser(
    prog="biser",
    description="Segmental duplication detection tool",
  )
  parser.add_argument(
    "--temp", "-T", help="Temporary directory location",
  )
  parser.add_argument(
    "--threads", "-t", type=int, default=1, help="Number of threads",
  )
  parser.add_argument(
    "genomes", nargs="+", help="Indexed genomes in FASTA format."
  )
  parser.add_argument(
    "--output", "-o", required=True, help="Indexed genomes in FASTA format."
  )
  parser.add_argument(
    "--hard", "-H", action="store_true", help="Are input genomes already hard-masked?"
  )
  parser.add_argument(
    "--keep-contigs", action="store_true",
    help=(
      "Do not ignore contigs, unplaced sequences, alternate alleles, patch chromosomes "
      "and mitochondrion sequences "
      "(i.e., chrM and chromosomes whose name contains underscore). "
      "Enable this when running BISER on scaffolds and custom assemblies."
    )
  )
  parser.add_argument(
    "--keep-temp", "-k", action="store_true",
    help="Keep temporary directory after the execution. Useful for debugging."
  )
  parser.add_argument(
    "--no-decomposition", action="store_true",
    help="Skip SD decomposition step."
  )
  args = parser.parse_args(argv)

  try:
    threads = args.threads

    genomes = {Path(path).stem: os.path.abspath(path) for path in args.genomes}
    with timing("BISER", force=True, indent=0):
      tmp = tempfile.mkdtemp(prefix='biser', dir=args.temp)
      print(f'Running BISER on {len(genomes)} genome(s): {", ".join(genomes)}')
      if args.keep_temp:
        print(f'Temporary directory: {tmp}')
      print()

      if not args.hard:
        orig_genomes = genomes.copy()
        results = []
        print('0. Hard-masking genomes')
        with timing("Hard-masking", results):
          os.makedirs(f'{tmp}/genomes')
          jobs = [(tmp, *g) for g in genomes.items()]
          with mp.Pool(threads) as pool:
            results[:] = list(progress(pool.imap(biser_mask, jobs), total=len(jobs)))
          for _, g, out in results:
            genomes[g] = out
        print()

      r = search(tmp, genomes, threads, args.keep_contigs)
      print()

      r = align(tmp, genomes, threads, r)
      print()

      if len(genomes) > 1:
        r = cross_biser(tmp, genomes, threads, r, args.keep_contigs)
        print()

      files = list(glob.glob(f'{tmp}/align/*.align'))
      if len(genomes) > 1:
        files += list(glob.glob(f'{tmp}/cross_align/*.align'))

      final = f'{tmp}/final.bed'
      with open(final, 'w') as fo:
        for fl in files:
          with open(fl) as f:
            for l in f:
              print(l, end='', file=fo)

      if not args.no_decomposition:
        decompose(tmp, genomes, threads, final)

      if not args.hard:
        run_biser(
          'translate', 
          '-o', args.output, 
          os.path.abspath(final), 
          *list(orig_genomes.values())
        )
      else:
        shutil.copy(final, args.output)
        if not args.no_decompostion:
          shutil.copy(f'{final}.elem.txt', f'{args.output}.elem.txt')

      if not args.keep_temp:
        shutil.rmtree(tmp, ignore_errors=True)

      print(f'Done! Results are available in {args.output}')
  except Exception as e:
    print(f'ERROR:', e)
    sys.exit(1)


def console():
  main(sys.argv[1:])


if __name__ == '__main__':
  console()
