<p align="center">
  <img src="https://emojipedia-us.s3.dualstack.us-west-1.amazonaws.com/thumbs/320/emojipedia/181/oyster_1f9aa.png" height=100 />
</p>

# BISER

[BISER](https://en.wiktionary.org/wiki/biser#Serbo-Croatian) (🦪🔮; Brisk Inference of Segmental duplication Evolutionary stRucture) is
a fast tool for detecting and decomposing segmental duplications (SDs) in a single genome
or multiple genomes.
BISER is [SEDEF](https://github.com/vpc-ccg/sedef)'s successor.

## Instalation

BISER needs Python 3.7+ and Samtools to run.

To install BISER, just run:
```bash
pip install biser
```

If you wish to build BISER from source, you will also need
[Seq programming language](https://docs.seq-lang.org/intro.html#install).
To install BISER from source, run:
```bash
pip install git+https://github.com/0xTCG/biser.git
```

## Usage
### Single genome

To find SDs in a single genome, just run:
```bash
biser -o <output> -t <threads> <genome.fa>
```

BISER will also produce a file called `output.elem` that will contain the elementary SD
decomposition of the found SDs.

All genomes should be indexed beforehand with `samtools faidx genome.fa`.

> ⚠️: BISER requires a soft-masked or a hard-masked genome assemblies for
> the optimal performance.
> Check for the presence of lowercase bases in your genome; if you have them,
> you are good to go.

> ⚠️: If you are experiences crashes on Linux machines (especially in cluster environments),
> try setting --gc-heap 1G (or higher).

### Multiple genomes

To find SDs in multiple genomes, just run:
```bash
biser -o <output> -j <jobs> <genome1.fa> <genome2.fa> ...
```

### Other options

```
Usage: biser [-h] [--temp TEMP] [--threads THREADS] --output OUTPUT [--hard]
             [--keep-contigs] [--keep-temp] [--no-decomposition]
             genomes [genomes ...]

Positional arguments:
  genomes               Indexed genomes in FASTA format.

Optional arguments:
  -h, --help            show this help message and exit
  --temp TEMP, -T TEMP  Temporary directory location
  --threads THREADS, -t THREADS
                        Number of threads
  --output OUTPUT, -o OUTPUT
                        Indexed genomes in FASTA format.
  --hard, -H            Are input genomes already hard-masked?
  --keep-contigs        Do not ignore contigs, unplaced sequences, alternate
                        alleles, patch chromosomes and mitochondrion sequences
                        (i.e., chrM and chromosomes whose name contains
                        underscore). Enable this when running BISER on
                        scaffolds and custom assemblies.
  --keep-temp, -k       Keep temporary directory after the execution. Useful
                        for debugging.
  --resume RESUME       Resume the previously interrupted run (that was run
                        with --keep-temp; needs the temp directory for
                        resume).
  --no-decomposition    Skip SD decomposition step.
  --max-error MAX_ERROR
                        Maximum SD error (large gaps includes).
  --max-edit-error MAX_EDIT_ERROR
                        Maximum SD edit error (large gaps NOT included).
  --max-chromosome-size MAX_CHROMOSOME_SIZE
                        Maximum chromosome size.
  --kmer-size KMER_SIZE
                        Search k-mer size.
  --winnow-size WINNOW_SIZE
                        Search winnow size.
  --version, -v         show program's version number and exit
  --gc-heap GC_HEAP     Set GC_INITIAL_HEAP_SIZE.
```

### Output format

The output follows the [BEDPE file format](https://bedtools.readthedocs.io/en/latest/content/general-usage.html#bedpe-format).

The first six (6) fields are the standard BEDPE fields describing the coordinates of SD mates:
- `chr1`, `start1` and `end1`
- `chr2`, `start2` and `end2` (both intervals are semi-open and 0-indexed).

Other fields are as follows:

| Field              | Description |
|--------------------|--------------------|
| `reference`        | Reference genome names of the first and the second mate, separated by `:`. |
| `score`            | Total alignment error (0--100\%): the number of mismatches and indels divided by the total alignment span.  |
| `strand1`          | Strand (`+` or `-`) of the first SD mate. |
| `strand2`          | Strand (`+` or `-`) of the second SD mate. |
| `max_len`          | Length of the longer mate. |
| `aln_len`          | Alignment span (mate length with gaps included) |
| `cigar`            | [CIGAR string](https://samtools.github.io/hts-specs/SAMv1.pdf) that describes the alignment |
| `optional`         | Optional fields in the format `NAME=VALUE;...`. Currently contains the mismatch rate (starts with `X=`) and the gap rate (starts with `ID=`). |

In addition to BEDPE output, BISER might also output the decomposition file (with the `.elem` extension) as well.
This file contains the list of core SD regions in the analyzed reference genomes.
The format of decomposition file is as follows:

| Field              | Description |
|--------------------|--------------------|
| `reference`        | Reference genome name. |
| `start`            | Start position of the core region (0-indexed). |
| `end`              | End position of the core region. |
| `id`               | Core region. Note that many regions share the same core ID because core regions are duplicated across the genome(s). |
| `len`              | Length of the core region. |
| `score`            | Core region score (internal use only). |
| `strand`           | Strand (`+` or `-`) of the core region. |

## Paper & Simulations

BISER was published in the [Algorithms for Molecular Biology](https://link.springer.com/article/10.1186/s13015-022-00210-2) and was presented at the [WABI 2021](https://drops.dagstuhl.de/opus/volltexte/2021/14368/pdf/LIPIcs-WABI-2021-15.pdf).

Please cite:
> Išerić, H., Alkan, C., Hach, F. et al. Fast characterization of segmental duplication structure in multiple genome assemblies. Algorithms Mol Biol 17, 4 (2022). https://doi.org/10.1186/s13015-022-00210-2

BibTeX entry:
```
@article{ivseric2022fast,
  title={Fast characterization of segmental duplication structure in multiple genome assemblies},
  author={I{\v{s}}eri{\'c}, Hamza and Alkan, Can and Hach, Faraz and Numanagi{\'c}, Ibrahim},
  journal={Algorithms for Molecular Biology},
  volume={17},
  number={1},
  pages={1--15},
  year={2022},
  publisher={Springer}
}
```

Paper simulations are available in [paper](paper/) directory.

## Changelog

- **BISER v1.4** (Mar 2023):
  - Change of alignment refinement heuristics (should be faster now)
    - Note: SDs generated with v1.4 might be slightly different
      than those generated by the earlier version
  - Switch to Codon
  - Minor bugfixes

## Contact

Please reach out to [Ibrahim Numanagić](mailto:inumanag_at_uvic_dot_ca).
