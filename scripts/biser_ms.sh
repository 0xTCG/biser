#!/bin/bash

echo "Start: `date`"

GETOPT="getopt"
# TIME="/usr/bin/time"
TIME="time"

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

output="biser_out"
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
r="500"
q="100"

while true; do
    case "$1" in
        -h|--help)
            echo "Usage: biser_ms.sh -o <output directory> -j <max. processes> [-f] <genomes_folder>"
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


# first we find all potential regins within same species and merge them
# ./new_sedef_ms.sh data/genomes/ -o same -l 1 -j 10 
mkdir ${output}/
${TIME} -f'First search took: %E' ./search_ms.sh ${input} -o "${output}/same/" -l ${l} -r ${r} -p ${p} -q ${q} -d ${d} -g ${g} -j ${jobs}

# now we do align of those putative SDs to find alignment and exact locations of those SDs
# ./align_ms.sh data/genomes/ same
${TIME} -f'First align took: %E' ./align_ms.sh ${input} "${output}/same/" ${jobs}

# now we extract all sequences from same species we aligned
# python3 chop_regions.py extract data/genomes sdregions same
mkdir ${output}/sdregions
${TIME} -f'Extracting regions took: %E' python3 chop_regions.py extract ${input} ${output}/sdregions ${output}/same

# now we search SD regions in other whole genomes
#./biser_diff_specs.sh data/genomes/ sdregions8 -o different8 -l 1 -f -j 8
echo 'Searching SDregions in genomes...'
${TIME} -f'Second search took: %E' ./biser_diff_specs.sh ${input} ${output}/sdregions -o ${output}/different -l ${l} -r ${r} -p ${p} -q ${q} -d ${d} -g ${g} -j ${jobs}

# normalize coordinates
python3 chop_regions.py normalize ${output}/different/

# now we align those sequences
# ./align_ms2.sh different8 data/genomes/
${TIME} -f'Second align took: %E' ./align_ms2.sh ${output}/different/ ${input}


# at the end just concatinate everything into one file
python3 chop_regions.py final ${output}/same/ ${output}/different/ ${output}/final.bed

${TIME} -f'Decomposition took: %E' union_find ${output}/final.bed >${output}/elementaries.txt 2>${output}/elementaries_log.txt

tail -1 ${output}/elementaries_log.txt