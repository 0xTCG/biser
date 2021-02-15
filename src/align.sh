#!/bin/bash


TIME=/usr/bin/time

if [[ ! -f "${TIME}" ]]; then
	TIME="/cvmfs/soft.computecanada.ca/gentoo/2020/usr/bin/time"
fi


path="$1"
fa_path="$2"
fa_path2="$3"
jobs="$4"




# these 2 calls are for mchanging names of hg19_chr1 to chr1
# mkdir "${path}/mergeds2"

# python3 quick_change.py "${path}mergeds/" "${path}mergeds2/"

# exit 1

extension="merged" #merged

input="${path}/${extension}"

rm -rf "${path}/aligned"
rm -rf "${path}/log_align"


mkdir "${path}/aligned"
mkdir "${path}/log_align"


logs="${path}/log_align"

output="${path}/aligned"


files=${input}/*
fa1=$fa_path
fa2=$fa_path2


# output_="test_out/merged" _merged
# exit 1

# echo "HERE ${files}"

for i in $files; do 

    # echo "here:: $i"
    filename1=$(basename -- "$i")
    # echo "${i}"

    # echo "sedef_out/hg19_hg19/aligned2/${filename1}"
    echo "/cvmfs/soft.computecanada.ca/gentoo/2020/usr/bin/time -f'TIMING %e %M' sedef/sedef align generate --k 10 ${fa1} ${i} ${fa2} >${output}/${filename1} 2>${logs}/${filename1}.log"
    
done | ${TIME} -f'Align time: %E' parallel --will-cite -j ${jobs} --joblog "${path}/align.joblog"

# cat ${output}/*.bed_merged >"${path}/final.bed"

sc_time=`grep TIMING ${logs}/*.log | awk '{s+=$2}END{print s}'`
sc_time_h=`echo "${sc_time} / 3600" | bc`
echo "Single-core running time: ${sc_time_h} hours (${sc_time} seconds)"


sc_mem=`grep TIMING ${logs}/*.log | awk '{if($3>m)m=$3}END{print m}'`
sc_mem_k=`echo "${sc_mem} / 1024" | bc`
echo "Memory used: ${sc_mem_k} MB"





