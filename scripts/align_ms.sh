#!/bin/bash

# this is script made for aligning all .bed files in one folder where user specifies input .fa file and folder where .bed files are

TIME=/usr/bin/time

if [[ ! -f "${TIME}" ]]; then
	TIME="/cvmfs/soft.computecanada.ca/gentoo/2020/usr/bin/time"
fi

path="$1"
beds="$2"
jobs="$3"


input="${path}"



files=${input}/*


# output_="test_out/merged" _merged
# exit 1

validchrs=${input}*50.fa
for i in $validchrs; do 

    for j in $validchrs; do 
        if [[ "$i" == "$j" ]]; then #[[ "$i" < "$j" ]] || 
            
            # here it begins
            filename1=$(basename -- "$i")
            filename1="${filename1%_*_*}"
            filename2=$(basename -- "$j")
            filename2="${filename2%_*_*}"
            files=${beds}/${filename1}_${filename2}/
            align_folder=${beds}/${filename1}_${filename2}/aligned/
            mkdir "${beds}/${filename1}_${filename2}/log_align/"
            log_folder=${beds}/${filename1}_${filename2}/log_align/
            # echo "this ${i}    ${j}    ${filename1}    ${filename2} ${files}"

            ${TIME} -f'TIMING %e %M' ./align.sh ${files} ${i} ${j} ${jobs} #>"${beds}/${filename1}_${filename2}/log_align/${filename1}_${filename2}.log"
            cat ${align_folder}* >${beds}/${filename1}_${filename2}/final.bed

        fi
    done
done 
# cat ${output}/*.bed_merged >"${path}/final.bed"

# sc_time=`grep TIMING ${logs}/*.log | awk '{s+=$2}END{print s}'`
# sc_time_h=`echo "${sc_time} / 3600" | bc`
# echo "Single-core running time: ${sc_time_h} hours (${sc_time} seconds)"

# sc_mem=`grep TIMING ${logs}/*.log | awk '{if($3>m)m=$3}END{print m}'`
# sc_mem_k=`echo "${sc_mem} / 1024" | bc`
# echo "Memory used: ${sc_mem_k} MB"


