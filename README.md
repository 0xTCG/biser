# BISER: BriskInference of Segmental duplication Evolutionary stRucture

Biser is a fast tool for:
* Detecting SDs in a genome
* Finding core SDs that drive the evolutionary process

Needed instalations:
* SEDEF
* Seq programming language


## Instalation
For compailing SEDEF, run these simple commands:
```bash
cd src/sedef
make -j release
```

By default, SEDEF uses Intel C++ compiler. If you are using g++, build with:
```bash
make -j release CXX=g++
```

If you are using Clang on macOS, compile as 
```bash
brew install libomp
make -j release OPENMP="-Xpreprocessor -fopenmp" CXX=clang++
```

> You need at least g++ 5.1.0 (C++14) to compile SEDEF. Clang should work fine as well.

SEDEF requires Boost libraries in order to compile. In case you installed Boost in a non-standard directory, you can still compile as follows:
```bash
CPATH={path_to_boost} make -j release
```

For installing Seq consult:
https://docs.seq-lang.org/intro.html#install

For compiling script for finding core SDs run:
```bash
g++ -O3 -o uf uf.cpp
```
Same note if you are using Boost in a non-standard directory, you should specify it with 
`CPATH={path_to_boost}`


## Usage
For finding SDs and decomposing them in elementary SDs in one specie, go to `src` directory and run:
```bash
./biser_1s.sh  -o <output> -j <jobs> <hard_masked_genome> 
```
For example, to run hg19_hard_masked.fa on 8 cores type:
```bash
./biser_1s.sh -o biser_out_hg19 -j 8 hg19_hard_masked.fa
```
Other possible parameters you can define:
* `-p` - change padding (default 5000)
* `-d` - change to dynamic gap
* `-l` - if you want to use filtering
* `-g` - if you want to use winnowing
* `f` - force, if output folder already exists

## Output
Output is in specified `output` folder. `{output}/final.bed` contains all SDs and `{output}/elementaries.txt` contains all elementaries in format:
* `elementary_id`, `chromosome`, `begin`, `end`,  `length`
* at the end are lines in format: `[core]`, `elementary_id`, `# of SDs that contain this core`

