#!/bin/bash
path="$1"
jobs="50"
TIME="/cvmfs/soft.computecanada.ca/gentoo/2020/usr/bin/time"

output="$2"
mkdir ${output}

for i in ${path}/*; do
    # echo ${i}
    f="$(basename -- $i)"
    # echo $f
    echo "spoa -r 1 $i >$output/$f"
done | ${TIME} -f'MSA time: %E' parallel --will-cite -j ${jobs}