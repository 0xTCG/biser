#!/bin/bash
beds_path="$1"
fa_path="$2"
output="$3"

jobs=15

TIME=/usr/bin/time

if [[ ! -f "${TIME}" ]]; then
	TIME="/cvmfs/soft.computecanada.ca/gentoo/2020/usr/bin/time"
fi




temp_folder=$(dirname "${beds_path}")
rm -rf ${temp_folder}/colors_dir
rm -rf ${temp_folder}/colors_fas
rm -rf ${output}


mkdir ${temp_folder}/colors_dir
mkdir ${temp_folder}/colors_fas
mkdir ${output}


# first we do big clustering of all SDs - merge 2 SD clades if they overlap or they have pairwise alignment
${TIME} -f 'Big clustering took %E (hh:mm:ss) time and %M KB memmory' ./big_cluster ${beds_path} ${temp_folder}/colors_dir/

# another step is to merge all regions in which SDs are contained and then to extract all sequences that correspond to those regions
${TIME} -f 'Extracting sequences took %E (hh:mm:ss) time and %M KB memmory' python3 chop_regions.py extract_colors ${fa_path} ${temp_folder}/colors_fas/ ${temp_folder}/colors_dir/


${TIME} -f 'Decomposition took %E (hh:mm:ss) time and %M KB memmory' ./decompose2.sh ${temp_folder}/colors_fas/ ${output}


cat ${output}/colors_elems/* >${output}/elementaries_fin.bed
