<p align="center">
  <img src="https://emojipedia-us.s3.dualstack.us-west-1.amazonaws.com/thumbs/320/emojipedia/181/oyster_1f9aa.png" height=100 />
</p>

# BISER

BISER (🦪🔮; Brisk Inference of Segmental duplication Evolutionary stRucture) is 
a fast tool for detecting and decomposing segmental duplications in genome assemblies.


## Instalation

BISER requires [Seq programming language](https://docs.seq-lang.org/intro.html#install),
Python 3.7+ and Samtools to run.

To install BISER, install Seq and then run:
```bash
pip install git+https://github.com/0xTCG/biser.git
```

## Usage
### Single genome

To find SDs in a single genome, just run:
```bash
biser -o <output> -t <threads> <genome.fa> 
```

BISER will also produce `output.elem.txt` file that contains the elementary SD
decomposition of the final SD set.

All genomes should be indexed beforehand with `samtools faidx genome.fa`.

> ⚠️ BISER requires soft-masked or hard-masked genome assemblies for 
> the optimal performance. 
> Check for the presence of lowercase bases in your genome; if you have them,
> you are good to go.

### Multiple genomes

To find SDs in multiple genomes, just run:
```bash
biser -o <output> -j <jobs> <genome1.fa> <genome2.fa> ...
```

### Other options

Run `biser -h` to list available options.

## Paper & Simulations

BISER was published in WABI 2021 proceedings and [is available here](https://drops.dagstuhl.de/opus/volltexte/2021/14368/pdf/LIPIcs-WABI-2021-15.pdf).

Paper simulations are available in [paper](paper/) directory.

> We will provide more experimental details soon.

## Contact

Please reach out to [Ibrahim Numanagić](mailto:inumanag_at_uvic_dot_ca).
