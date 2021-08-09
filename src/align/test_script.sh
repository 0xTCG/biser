#!/bin/bash

genome="/home/hiseric1/new_sedef/seq/search/data/genomes/hg19_hard_50.fa"
TIME=/usr/bin/time

if [[ ! -f "${TIME}" ]]; then
	TIME="/cvmfs/soft.computecanada.ca/gentoo/2020/usr/bin/time"
fi

p="5000"
d='0'
g='0'
l='1'
jobs=6

module load samtools




for i in 500,100 700,100 750,150 750,200 800,300 1000,300; do
    IFS=","; set -- $i; 
    r=$1
    q=$2
    output="results/hg19_1_${r}_${q}"
    echo ${r} ${q}
    # echo "./biser_1s.sh ${genome} -r ${r} -q ${q} -p ${p} -d ${d} -g ${g} -l ${l} -o ${output} -j 8 -f"
    ./biser_1s.sh ${genome} -r ${r} -q ${q} -p ${p} -d ${d} -g ${g} -l ${l} -o ${output} -j 8 -f
    # -r ${r} -q ${q}
    # ./biser_1s.sh ${genome} -p ${p} -d ${d} -g ${g} -l ${l} -o ${output} -r ${r} -q ${q} -j 8 -f
	# ${TIME} -f'TIMING %e %M' seqc biser_search.seq -p ${p} -d ${d} -g ${g} -f ${l} -k 14 -w 16 -r ${r} -q ${q} -o ${output}/seeds ${input} ${input} chr1

	# echo "${TIME} -f'TIMING %e %M' seqc biser_search.seq -p ${p} -d ${d} -g ${g} -f ${l} -k 14 -w 16 -r ${r} -q ${q} -o ${output}/seeds ${input} ${input} $i >${output}/log/seeds/${i}_.log 2>${output}/log/seeds/${i}_2.log"
    # break

done #| ${TIME} -f'All time: %E' parallel --will-cite -j ${jobs} --joblog "results/test.joblog"