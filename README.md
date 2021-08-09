# BISER: Brisk Inference of Segmental duplication Evolutionary stRucture

Biser is a fast tool for detecting and decomposing SDs in a one or multiple genomes.

Needed installations:
* SEDEF
* Seq programming language


## Instalation

BISER's depends on Seq, and to install it consult:
https://docs.seq-lang.org/intro.html#install

For installign BISER, run the following:
```bash
make
```

BISER requires Boost libraries in order to compile. In case you installed Boost in a non-standard directory, you can still compile as follows:
```bash
CPATH={path_to_boost} make
```

## Usage
### One species
For finding SDs and decomposing them into elementary SDs in one species, go to `src` directory and run:
```bash
./biser_1s.sh  -o <output> -j <jobs> <hard_masked_genome> 
```
For example, to run hg19_hard_masked.fa on 8 cores, type:
```bash
./biser_1s.sh -o biser_out_hg19 -j 8 hg19_hard_masked.fa
```

### Multiple species
For finding SDs and decomposing them into elementary SDs in multiple species, go to `src` directory and run:
```bash
./biser_ms.sh  -o <output> -j <jobs> <destionation_where_hard_masked_genomes_are> 
```
For example, to run hg19_hard_masked.fa on 8 cores, type:
```bash
./biser_ms.sh -o biser_out_ms -j 8 data/genomes/
```


All genome names must be in the following format `{species_name}_hard_50.fa`

Optional parameters you can set (for both, one and multiple species):
* `-p <padding value>` - change padding (default 5000)
* `-d <1/0>` - change to dynamic gap (max SD length set to 1Mbp)
* `-l <1/0>` - if you want to use filtering
* `-g <1/0>` - if you want to use winnowing
* `-f` - force, if output folder already exists

## Output
Output is generated in the specified `output` folder. `{output}/final.bed` contains all SDs and `{output}/elementaries.txt` contains all elementaries in format:
* `elementary_id`, `chromosome`, `begin`, `end`,  `length`
* at the end are lines in format: `[core]`, `elementary_id`, `# of SDs that contain this core`

For multiple species, chromosome names are written in the following format `{species_name}#{chromosome}`

## Simulations
Path for SD detection simulations:
https://github.com/0xTCG/biser/blob/master/src/plot_simulations.ipynb   


Path for core detection simulations:
https://github.com/0xTCG/biser/blob/master/src/simulate_cores.ipynb   

