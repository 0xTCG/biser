#!/bin/bash


path="$1"
jobs="15"
TIME="/cvmfs/soft.computecanada.ca/gentoo/2020/usr/bin/time"

output="$2"
rm -rf ${output}

mkdir ${output}
mkdir ${output}/colors_elems
mkdir ${output}/colors_logs
module load samtools

for i in ${path}/*; do

    f="$(basename -- $i)"
    if [[ $i == *.fa ]]; then
        samtools faidx $i
        echo "${TIME} -f'TIMING %e %M' ./decompose ${i} ${output}/colors_elems/ >${output}/colors_logs/${f}_1 2>${output}/colors_logs/${f}_2"

    fi
done | ${TIME} -f'Decomposition time: %E' parallel --will-cite -j 15


# ./decompose2.sh /home/hiseric1/new_sedef/biser/src/results/5_all/same/hg19_hg19/colors_fas/ results/hg19_5_2/
# ./decompose2.sh results/test1 results/hg19_5_2/
