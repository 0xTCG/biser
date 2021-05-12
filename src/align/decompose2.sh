#!/bin/bash
# shopt -s nullglob
# numfiles=(results/test1/*)
# numfiles=${#numfiles[@]}
# echo  ${numfiles}



# --

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
    # echo ${i}
    # echo `wc -l ${i}`
    # /cvmfs/soft.computecanada.ca/gentoo/2020/usr/bin/time -f'TIMING %e %M' seqc decompose.seq >l1
    
    # echo samtools faidx ${i}
    f="$(basename -- $i)"
    # shopt -s nullglob
    # numfiles=(${i}/*.fa)
    # numfiles=${#numfiles[@]}
    # echo  ${numfiles}

    # echo "sleep 1"
    

    # echo "/cvmfs/soft.computecanada.ca/gentoo/2020/usr/bin/time -f'TIMING %e %M' seqc decompose.seq ${i} ${output}/colors_elems/  >${output}/colors_logs/${f}_1 2>${output}/colors_logs/${f}_2"
    # echo "/cvmfs/soft.computecanada.ca/gentoo/2020/usr/bin/time -f'TIMING %e %M' seqc run -release decompose.seq ${i} ${output}/colors_elems/ ${numfiles} >${output}/colors_logs/${f}_1 2>${output}/colors_logs/${f}_2"
    if [[ $i == *.fa ]]; then
        # samtools faidx $i

        echo "/cvmfs/soft.computecanada.ca/gentoo/2020/usr/bin/time -f'TIMING %e %M' ./decompose ${i} ${output}/colors_elems/ >${output}/colors_logs/${f}_1 2>${output}/colors_logs/${f}_2"
    fi
    # /cvmfs/soft.computecanada.ca/gentoo/2020/usr/bin/time -f'TIMING %e %M' seqc decompose.seq ${i} results/hg19_5/colors_elems/ >results/l1
    # break
done | ${TIME} -f'Decomposition time: %E' parallel --will-cite -j 15


# ./decompose2.sh /home/hiseric1/new_sedef/biser/src/results/5_all/same/hg19_hg19/colors_fas/ results/hg19_5_2/
# ./decompose2.sh results/test1 results/hg19_5_2/
