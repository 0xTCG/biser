
from intervaltree import *
import sys

def how_much_overlaps(i1, i2):
    if i1.overlaps(i2):
        b1 = max(i1.begin, i2.begin)
        e1 = min(i1.end, i2.end)
        return e1 - b1
    return 0

def compare_cores(main_core_file, uf_cores, out_file, sd_file):
    # here build intervsal tree for SDs and check if cores are contained there
    ds_int = IntervalTree()

    for i in open(sd_file, 'r'):
        line = i.split('\t')
        ds_int.add(Interval(int(line[1]), int(line[2])))
        ds_int.add(Interval(int(line[4]), int(line[5])))

    ds_int.merge_overlaps()



    uf_int = IntervalTree()
    d1 = dict()
    cores = []
    for i in open(uf_cores, 'r'):
        line = i.split('\t')
        # print (line)
        if len(line) == 1:
            line = line[0].split(' ')
            cores.append(line[1])
            continue

        s = int(line[2])
        e = int(line[3])
        k = line[0]
        if e>s:
            if uf_int.overlaps(s,e):
                assert False
            uf_int[s:e] = k
        elif s>e:
            # uf_int[e:s] = k
            print (line)
            assert False

        if k in d1:
            d1[k].append(Interval(s,e))
        else:
            d1[k] = [Interval(s,e)]
    print (f'Cores: {cores}')
    d = dict()
    contains_sd = 0
    for i in open(main_core_file, 'r'):
        line = i.split('\t')
        s = int(line[1])
        e = int(line[2])
        k = int(line[3])
        sum_2 = 0
        interval__ = Interval(s,e)
        for ovrlp in ds_int.overlap(s,e):
            sum_2 += how_much_overlaps(ovrlp, interval__)
        if sum_2 > (interval__.end - interval__.begin)*0.9:
            contains_sd += 1


        if k in d:
            d[k].append(Interval(s,e, main_core_file))
        else:
            d[k] = [Interval(s,e, main_core_file)]
    
    main_counter = 0
    main_counter2 = 0
    all_covered = 0
    all_ = 0
    core_containers = 0
    for k in d:
        is_it_ok = True
        overall_sum1 = 0
        overall_sum2 = 0
        alowd_keys = set()
        one_covered = 0
        # print (f'Key: {k}')
        contains_core = False
        for interval in d[k]:
            all_ += 1
            # print (f'Interval: {interval}')
            overlaps = uf_int.overlap(interval.begin, interval.end) # uf_int[interval.begin, interval.end]
            overall_sum1 += interval.end - interval.begin
            sum_ = 0
            
            for i__ in overlaps:
                
                sum_ += how_much_overlaps(i__, interval)
                if i__.data in cores:
                    # print (f'Covered w cores: {i}')
                    contains_core = True
                # else:
                #     print (f'But these: {i}')

                
            if sum_ < (interval.end - interval.begin) * 0.9:
                is_it_ok = False
            else:
                one_covered += 1

            overall_sum2 += sum_
        if overall_sum2 >= overall_sum1 * 0.9: # is_it_ok:
            main_counter+= 1
        if one_covered :
            main_counter2 += 1
        else:
            print ("NOT FOUND")
        if contains_core:
            core_containers += 1
        # else:
        #     print (k, d[k])
        all_covered += one_covered
        # else:
        #     print (d[k])
        #     print( overall_sum2/ overall_sum1 )
    f = open(out_file, 'w')
    f.write(f'{all_covered}\t{all_}\t{main_counter2}\t{main_counter}\t{len(d)}\t{all_covered / all_}\t{main_counter2/len(d)}\t{core_containers}\t{contains_sd}')
    print (main_counter, core_containers, main_counter2, len(d), contains_sd)
    # for i in d:

    #     print(d[i])


print (sys.argv)
if len(sys.argv) >= 4:
    compare_cores(sys.argv[1], sys.argv[2], sys.argv[3])
# else:
#     compare_cores('simulations_fin/simulations/cores#1#20#2391.bed', 'simulations_fin/elementaries/sequence1.simulations#1#20#2391_sequence1.simulations#1#20#2391_n.bed_merged')


    


