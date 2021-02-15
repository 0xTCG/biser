#!/bin/bash
jobs=8


TIME=time


input=$1
genomes=$2

mkdir ${input}/log/aligned
for i in ${input}/* ; do
    # echo ${i}
    if [[ -d $i ]]; then
        # echo "$i is a directory"
        for j in ${i}/*.bed ; do
            if [[ -f $j ]]; then
                filename=$(basename -- "$i")
                filename1="${filename%_*}"
                filename2="${filename##*_}"

                echo ${filename1} ${filename2} ${j} ${filename} ${genomes}/${filename1}_hard_50.fa ${genomes}/${filename2}_hard_50.fa ${i}
                if [[ "$filename1" < "$filename2" ]]; then
                    fa1="${genomes}/${filename1}_hard_50.fa"
                    fa2="${genomes}/${filename2}_hard_50.fa"
                else
                    fa1="${genomes}/${filename2}_hard_50.fa"
                    fa2="${genomes}/${filename1}_hard_50.fa"
                fi
                
                # mkdir ${i}/merged

                sedef/sedef merge ${i}/ ${i}/merged/
                # rm -rf ${i}/merged
                # rm -rf ${i}/merged.bed
                # rm  ${i}/merged.bed
                # rm  ${i}/final.bed

                


                # cat ${i}/merged/seeds.bed >${i}/merged.bed
                # echo "/cvmfs/soft.computecanada.ca/gentoo/2020/usr/bin/time -f'TIMING %e %M' /home/hiseric1/new_sedef/sedef/sedef align generate --k 10 ${fa1} ${i}/merged/seeds.bed_merged ${fa2} >${i}/final.bed 2>${input}/log/aligned/${filename1}_${filename2}.log"
                
            
            fi
        done
    fi
done

for i in ${input}/* ; do
    # echo ${i}
    if [[ -d $i ]]; then
        # echo "$i is a directory"
        for j in ${i}/*.bed ; do
            if [[ -f $j ]]; then
                filename=$(basename -- "$i")
                filename1="${filename%_*}"
                filename2="${filename##*_}"

                # echo ${filename1} ${filename2} ${j} ${filename} ${genomes}/${filename1}_hard_50.fa ${genomes}/${filename2}_hard_50.fa ${i}
                if [[ "$filename1" < "$filename2" ]]; then
                    fa1="${genomes}/${filename1}_hard_50.fa"
                    fa2="${genomes}/${filename2}_hard_50.fa"
                else
                    fa1="${genomes}/${filename2}_hard_50.fa"
                    fa2="${genomes}/${filename1}_hard_50.fa"
                fi
                
                # mkdir ${i}/merged

                # echo "/home/hiseric1/new_sedef/sedef/sedef merge ${i}/ ${i}/merged/"
                # rm -rf ${i}/merged
                # rm -rf ${i}/merged.bed
                # rm  ${i}/merged.bed
                # rm  ${i}/final.bed

                


                # cat ${i}/merged/seeds.bed >${i}/merged.bed
                echo "/cvmfs/soft.computecanada.ca/gentoo/2020/usr/bin/time -f'TIMING %e %M' sedef/sedef align generate --k 10 ${fa1} ${i}/merged/seeds.bed_merged ${fa2} >${i}/final.bed 2>${input}/log/aligned/${filename1}_${filename2}.log"
                
            
            fi
        done
    fi
done | tee "align.comm" | ${TIME} -f'All align time time: %E' parallel --will-cite -j ${jobs} --bar --joblog "align.joblog"

sc_time=`grep TIMING ${input}/log/aligned/*.log | awk '{s+=$2}END{print s}'`
sc_time_h=`echo "${sc_time} / 3600" | bc`
echo "Single-core running time: ${sc_time_h} hours (${sc_time} seconds)"

sc_mem=`grep TIMING ${input}/log/aligned/*.log | awk '{if($3>m)m=$3}END{print m}'`
sc_mem_k=`echo "${sc_mem} / 1024" | bc`
echo "Memory used: ${sc_mem_k} MB"
