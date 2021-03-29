#!/bin/bash
#/usr/bin/env bash
# 786

echo "Start: `date`"

GETOPT="getopt"

TIME=/usr/bin/time

if [[ ! -f "${TIME}" ]]; then
	TIME="/cvmfs/soft.computecanada.ca/gentoo/2020/usr/bin/time"
fi


if [[ "$OSTYPE" == "linux-gnu" ]]; then
	:
elif [[ "$OSTYPE" == "darwin"* ]]; then
	GETOPT="$(brew --prefix gnu-getopt)/bin/getopt"
	TIME="gtime"
else
	echo "Unknown environment ${OSTYPE} --- use at your own risk!"
fi

$GETOPT --test > /dev/null
if [[ $? -ne 4 ]]; then
	echo "I’m sorry, `$GETOPT --test` failed in this environment."
	exit 1
fi

PATH="${PATH}:"`pwd`

if ! command -v "samtools" >/dev/null 2>&1 ; then
	echo "Samtools not found in \$PATH (${PATH})"
	# exit 1
	module load samtools
fi

if ! command -v "parallel" >/dev/null 2>&1 ; then
	echo "GNU Parallel not found in \$PATH (${PATH})"
	exit 1
fi

OPTIONS=hj:o:w:fe:t:S:d:g:l:p:r:q:
LONGOPTIONS=help,jobs,output,wgac,force,exclude,translate,stat-params,dynamic,withoutw,filtering,padding,ref,query
PARSED=$($GETOPT --options=$OPTIONS --longoptions=$LONGOPTIONS --name "$0" -- "$@")
if [[ $? -ne 0 ]]; then
	exit 2
fi
eval set -- "$PARSED"

output="sedef_out"
jobs=10
force="n"
wgac=""
translate=""
exclude="^(chr|trans)[0-9XY]+$"
stat_params=""
d="0"
g="0"
l="1"
p="5000"
r="700"
q="100"

while true; do
	case "$1" in
		-h|--help)
			echo "Usage: sedef.sh -o <output directory> -j <max. processes> [-f] <genome.fa>"
			echo "-f removes output directory if it exists; -h shows this message"
			exit 0
			;;
		-f|--force)
			force="y"
			shift
			;;
		-t|--translate)
			translate="$2"
			shift 2
			;;
		-w|--wgac)
			wgac="$2"
			shift 2
			;;
		-d|--dynamic)
			d="$2"
			shift 2
			;;
		-g|--withoutw)
			g="$2"
			shift 2
			;;
		-l|--filtering)
			l="$2"
			shift 2
			;;
		-p|--padding)
			p="$2"
			shift 2
			;;
		-r|--ref)
			r="$2"
			shift 2
			;;
		-q|--query)
			q="$2"
			shift 2
			;;
		-o|--output)
			output="$2"
			shift 2
			;;
		 -j|--jobs)
			jobs="$2"
			shift 2
			;;
		 -e|--exclude)
			exclude="$2"
			shift 2
			;;
		-S|--stat-params)
			stat_params="$2"
			shift 2
			;;
		--)
			shift
			break
			;;
		*)
			echo "Programming error"
			exit 3
			;;
	esac
done

input="$1"

echo "BISER: FASTA=${input}; output=${output}; jobs=${jobs}; force=${force}"


if [ -e "${output}" ]; then
    echo -n "Output file name ${output} exists!"
    if [ "${force}" == "y" ] ; then
    	echo " Removing it."
    	rm -rf "${output}"
    else
    	echo " Please delete ${output} or run with -f/--force if you want to start anew."
    	# exit 1
    fi
fi
mkdir -p "${output}"


mkdir -p "${output}/seeds"
mkdir -p "${output}/merged"
mkdir -p "${output}/log/seeds"


# validchrs=${input}/*50.fa
validchrs="$(cut -f1 "${input}.fai" | awk '$1~/'${exclude}'/')"
export OMP_NUM_THREADS=1

echo "************************************************************************"
if [ ! -f "${output}/seeds.joblog.ok" ] || [ "${force}" == "y" ]; then
	rm -f "${output}/seeds.joblog.ok"
	echo "Running SD seeding..."
	# ARRAY=()
	for i in $validchrs; do 
        
		echo "${TIME} -f'TIMING %e %M' seqc biser_search.seq -p ${p} -d ${d} -g ${g} -f ${l} -k 14 -w 16 -r ${r} -q ${q} -o ${output}/seeds ${input} ${input} $i >${output}/log/seeds/${i}_.log 2>${output}/log/seeds/${i}_2.log"
		
	done | tee "${output}/seeds.comm" | ${TIME} -f'Seeding time: %E' parallel --will-cite -j ${jobs} --bar --joblog "${output}/seeds.joblog"
	
	curent_dic=`pwd`
	echo ${curent_dic}
	time sedef/sedef merge "${curent_dic}/${output}/seeds/" "${curent_dic}/${output}/merged/"

	cat ${output}/merged/*.bed_merged >"${output}/seeds.bed"

	echo "SD seeding done: done running ${proc} jobs!"

	# Get the single-core running time
	sc_time=`grep TIMING ${output}/log/seeds/*.log | awk '{s+=$2}END{print s}'`
	sc_time_h=`echo "${sc_time} / 3600" | bc`
	echo "Single-core search running time: ${sc_time_h} hours (${sc_time} seconds)"

	filtered=`grep FILTERED ${output}/log/seeds/*.log | awk '{s+=$2}END{print s}'`

	echo "FIltered: ${filtered} kmers."

	sc_mem=`grep TIMING ${output}/log/seeds/*.log | awk '{if($3>m)m=$3}END{print m}'`
	sc_mem_k=`echo "${sc_mem} / 1024" | bc`
	echo "Memory used for search: ${sc_mem_k} MB"

	touch "${output}/seeds.joblog.ok"
fi


# now do merging part 
# here path to your sedef file
# ./sedef/sedef merge "${output}/seeds/" "${output}/merged/"


# and align at the end
files=${output}/merged/*
logs=${output}/log/align/

fa1=$input
fa2=$input


mkdir "${output}/aligned/"
mkdir -p "${output}/log/align/"

output2="${output}/aligned/"
echo ${output2}
count_all=0
for i in $files; do 

    filename1=$(basename -- "$i")
	# also put here to sedef exe
    echo "/cvmfs/soft.computecanada.ca/gentoo/2020/usr/bin/time -f'TIMING %e %M' sedef/sedef align generate --k 10 ${fa1} ${i} ${fa2} >${output2}/${filename1} 2>${logs}/${filename1}.log"
	((count_all=count_all+1))
done | ${TIME} -f'Align time: %E' parallel --will-cite -j ${jobs} --joblog "${output}/align.joblog"

# cat ${output}/*.bed_merged >"${path}/final.bed"

sc_time=`grep TIMING ${logs}/*.log | awk '{s+=$2}END{print s}'`
sc_time_h=`echo "${sc_time} / 3600" | bc`
echo "Single-core align running time: ${sc_time_h} hours (${sc_time} seconds)"


sc_mem=`grep TIMING ${logs}/*.log | awk '{if($3>m)m=$3}END{print m}'`
sc_mem_k=`echo "${sc_mem} / 1024" | bc`
echo "Memory used for align: ${sc_mem_k} MB"

#Finished
finished=`grep Finished ${logs}/*.log | awk '{s+=1}END{print s}'`
echo "Finished: ${finished}/${count_all}"

cat ${output2}/* >${output}/final.bed
uf ${output}/final.bed >${output}/elementaries.txt
