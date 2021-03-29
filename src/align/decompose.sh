#!/bin/bash
path="$1"
jobs=10
TIME="/cvmfs/soft.computecanada.ca/gentoo/2020/usr/bin/time"

output="$2"
mkdir ${output}
mkdir ${output}/colors_elems
mkdir ${output}/colors_logs

for i in ${path}/*; do
    # echo ${i}
    # echo `wc -l ${i}`
    # /cvmfs/soft.computecanada.ca/gentoo/2020/usr/bin/time -f'TIMING %e %M' seqc decompose.seq >l1
    if [[ $i == *.fa ]]
    then
        # echo samtools faidx ${i}
        f="$(basename -- $i)"
        # echo "/cvmfs/soft.computecanada.ca/gentoo/2020/usr/bin/time -f'TIMING %e %M' seqc decompose.seq ${i} results/hg19_5/colors_elems/ >results/hg19_5/colors_logs/${f}_1 2>results/hg19_5/colors_logs/${f}_1"
        /cvmfs/soft.computecanada.ca/gentoo/2020/usr/bin/time -f'TIMING %e %M' seqc run decompose.seq ${i} ${output}/colors_elems/ >${output}/colors_logs/${f}_1 2>${output}/colors_logs/${f}_2

    fi
    # /cvmfs/soft.computecanada.ca/gentoo/2020/usr/bin/time -f'TIMING %e %M' seqc decompose.seq ${i} results/hg19_5/colors_elems/ >results/l1
    # break
done | ${TIME} -f'Decomposition time: %E' parallel --will-cite -j ${jobs} --joblog "${path}/decompose.joblog"


# ./decompose.sh /home/hiseric1/new_sedef/biser/src/results/5_all/same/hg19_hg19/colors_fas/ results/hg19_5_2/
