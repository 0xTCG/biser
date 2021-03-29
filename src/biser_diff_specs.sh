#!/bin/bash
#/usr/bin/env bash
# 786

echo "Start: `date`"

GETOPT="getopt"
# TIME="/usr/bin/time"

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
r="100"
q="700"

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

# if [[ $# -ne 1 ]]; then
# 	echo "$0: FASTA file is required."
# 	exit 1
# fi
input="$1"
input2="$2"


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

validchrs=${input}*50.fa
# validchrs="$(cut -f1 "${input}.fai" | awk '$1~/'${exclude}'/')"
export OMP_NUM_THREADS=1
mkdir "${output}/log"
mkdir "${output}/log/seeds"
curent_dic=`pwd`

echo "************************************************************************"
if [ ! -f "${output}/seeds.joblog.ok" ] || [ "${force}" == "y" ] ; then # || true ; then
    rm -f "${output}/seeds.joblog.ok"
    echo "Running SD seeding..."
    for i in $validchrs; do


        for j in $validchrs; do 
            filename1=$(basename -- "$i")
            filename1="${filename1%_*_*}"
            filename2=$(basename -- "$j")
            filename2="${filename2%_*_*}"
            # here if one sd regions file has less lines prefear thet one
            lines1=`wc -l ${input2}/SD_regions_${filename1}.fa | cut -f1 -d' '`
            lines2=`wc -l ${input2}/SD_regions_${filename2}.fa | cut -f1 -d' '`
            # echo ${i} ${j} ${lines1}
            if [[ ! "$filename1" == "$filename2" ]]; then
                if [[ "$lines1" < "$lines2" ]] || [[ "$lines1" == "$lines2" ]]; then
                    # here it begins
                    
                    # k=`wc -l "${input2}/SD_regions_${filename1}.fa"`
                    # echo ${filename1} ${filename2} ${i} ${j} ${k} "${input2}/SD_regions_${filename1}.fa" ${lines1} ${lines2}

                    mkdir "${output}/${filename1}_${filename2}"
                    mkdir "${output}/${filename1}_${filename2}/seeds"
                    mkdir "${output}/${filename1}_${filename2}/merged"
                    mkdir "${output}/${filename1}_${filename2}/aligned"
                    samtools faidx ${input2}/SD_regions_${filename1}.fa
                    # final one:
                    echo "/cvmfs/soft.computecanada.ca/gentoo/2020/usr/bin/time -f'TIMING %e %M' seqc biser_search.seq -p ${p} -d ${d} -g ${g} -f ${l} -r ${r} -q ${q} -k 14 -w 16 -o ${output}/${filename1}_${filename2}/seeds ${input2}/SD_regions_${filename1}.fa $j >${output}/log/seeds/${filename1}_${filename2}.log 2>${output}/log/seeds/${filename1}_${filename2}_2.log"
                    
                    # echo "/home/hiseric1/new_sedef/sedef/sedef merge ${curent_dic}/${output}/${filename1}_${filename2}/seeds/ ${curent_dic}/${output}/${filename1}_${filename2}/merged/"
                fi
            fi
        done
    done | tee "${output}/seeds.comm" | ${TIME} -f'Seeding time: %E' parallel --will-cite -j ${jobs} --bar --joblog "${output}/seeds.joblog"

    curent_dic=`pwd`
    echo ${curent_dic}

    # cat ${output}/merged/*.bed_merged >"${output}/seeds.bed"

    # echo "SD seeding done: done running ${proc} jobs!"

    # Get the single-core running time
    sc_time=`grep TIMING ${output}/log/seeds/*.log | awk '{s+=$2}END{print s}'`
    sc_time_h=`echo "${sc_time} / 3600" | bc`
    echo "Single-core running time: ${sc_time_h} hours (${sc_time} seconds)"

    filtered=`grep FILTERED ${output}/log/seeds/*.log | awk '{s+=$2}END{print s}'`

    echo "FIltered: ${filtered} kmers."

    sc_mem=`grep TIMING ${output}/log/seeds/*.log | awk '{if($3>m)m=$3}END{print m}'`
    sc_mem_k=`echo "${sc_mem} / 1024" | bc`
    echo "Memory used: ${sc_mem_k} MB"

    touch "${output}/seeds.joblog.ok"
fi



