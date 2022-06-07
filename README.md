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

> ⚠️ BISER requires a soft-masked or a hard-masked genome assemblies for 
> the optimal performance. 
> Check for the presence of lowercase bases in your genome; if you have them,
> you are good to go.

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
  -h, --help                  Show the help message and exit.
  --temp TEMP, -T TEMP        Temporary directory location. Must exist beforehand.
  --threads N, -t N           Number of cores to use.
  --output OUTPUT, -o OUTPUT  Output filename.
  --hard, -H                  Pass if the input genomes are already hard-masked.
  --keep-contigs              Do not ignore contigs, unplaced sequences, alternate
                              alleles, patch chromosomes and mitochondrion sequences
                              (e.g., chrM and chromosomes whose name contains underscore). 
                              Enable this when running BISER on scaffolds and custom assemblies.
  --keep-temp, -k             Keep temporary directory after the execution. Useful for debugging.
  --no-decomposition          Skip the SD decomposition step.
```

### Output format

The output follows the [BEDPE file format](https://bedtools.readthedocs.io/en/latest/content/general-usage.html#bedpe-format).

The first six (6) fields are the standard BEDPE fields describing the coordinates of SD mates:
- `chr1`, `start1` and `end1`
- `chr2`, `start2` and `end2`

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
| `start`            | Start position of the core region. |
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

> We will provide more experimental details soon.

## Contact

Please reach out to [Ibrahim Numanagić](mailto:inumanag_at_uvic_dot_ca).
