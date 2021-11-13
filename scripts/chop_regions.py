from intervaltree import *
import os
import sys

def build_dict(path, d = dict()):
    count_sds = 0
    threshold1 = 5000000
    threshold2 = 10000000
    threshold3 = 20000000

    counters = [0,0,0]
    if os.path.isdir(path):
        for i in os.listdir(path):
            # print (i)
            for j in open(path + '/' + i, 'r'):
                line = j.split('\t')
                chr1 = line[0]
                chr2 = line[3]

                if chr1 == 'chrM' or chr2 == 'chrM':
                    continue
                count_sds += 1
                s1 = int(line[1])
                e1 = int(line[2])

                s2 = int(line[4])
                e2 = int(line[5])

                # if e1 - s1 > 5000000:
                #     print (j)
                if e1 - s1 > threshold1:
                    counters[0] += 1
                if e1 - s1 > threshold2:
                    counters[1] += 1
                if e1 - s1 > threshold3:
                    counters[2] += 1

                if chr1 in d:
                    d[chr1].add(Interval(s1,e1))
                else:
                    d[chr1] = IntervalTree()
                    d[chr1].add(Interval(s1,e1))

                if chr2 in d:
                    d[chr2].add(Interval(s2,e2))
                else:
                    d[chr2] = IntervalTree()
                    d[chr2].add(Interval(s2,e2))
    else:
        for j in open(path, 'r'):
            line = j.split('\t')
            chr1 = line[0]
            if chr1 == '#chr1':
                continue
            chr2 = line[3]

            if chr1 == 'chrM' or chr2 == 'chrM':
                continue
            count_sds += 1
            s1 = int(line[1])
            e1 = int(line[2])

            s2 = int(line[4])
            e2 = int(line[5])

            # if e1 - s1 > 5000000:
            #     print (j)
            if e1 - s1 > threshold1:
                counters[0] += 1
            if e1 - s1 > threshold2:
                counters[1] += 1
            if e1 - s1 > threshold3:
                counters[2] += 1

            if chr1 in d:
                d[chr1].add(Interval(s1,e1))
            else:
                d[chr1] = IntervalTree()
                d[chr1].add(Interval(s1,e1))

            if chr2 in d:
                d[chr2].add(Interval(s2,e2))
            else:
                d[chr2] = IntervalTree()
                d[chr2].add(Interval(s2,e2))
    
    coverage = 0
    # coverage_before = 0
    for i in d:
        # for inerv in d[i]:
        #     coverage_before += inerv.end - inerv.begin
        d[i].merge_overlaps()
        for inerv in d[i]:
            coverage += inerv.end - inerv.begin
    # km = path.split("/")[1]
    # print (coverage_before)
    # print (f'{km}\t{counters[0]}\t{counters[1]}\t{counters[2]}')
    return d, count_sds, coverage

# d, i = build_dict('same8/mm8_mm8/final.bed')
# /home/hiseric1/new_sedef/sedef/hg_2/final.bed same8/mm8_mm8/final.bed
# d__, i__, j__ =  build_dict('/home/hiseric1/new_sedef/sedef/hg_3/final.bed', dict())
# print(i__, j__)

def build_dict_2(path, d = dict()):
    count_sds = 0
    threshold1 = 5000000
    threshold2 = 10000000
    threshold3 = 20000000

    counters = [0,0,0]
    if os.path.isdir(path):
        for i in os.listdir(path):
            # print (i)
            for j in open(path + '/' + i, 'r'):
                line = j.split('\t')
                chr1 = line[0]
                chr2 = line[3]

                if chr1 == 'chrM' or chr2 == 'chrM':
                    continue
                count_sds += 1
                s1 = int(line[1])
                e1 = int(line[2])

                s2 = int(line[4])
                e2 = int(line[5])

                # if e1 - s1 > 5000000:
                #     print (j)
                if e1 - s1 > threshold1:
                    counters[0] += 1
                if e1 - s1 > threshold2:
                    counters[1] += 1
                if e1 - s1 > threshold3:
                    counters[2] += 1

                chr1 += line[8]
                chr2 += line[9] if line[9] == '+' else '_'

                if chr1 in d:
                    d[chr1].add(Interval(s1,e1))
                else:
                    d[chr1] = IntervalTree()
                    d[chr1].add(Interval(s1,e1))

                if chr2 in d:
                    d[chr2].add(Interval(s2,e2))
                else:
                    d[chr2] = IntervalTree()
                    d[chr2].add(Interval(s2,e2))
    else:
        for j in open(path, 'r'):
            line = j.split('\t')
            chr1 = line[0]
            if chr1 == '#chr1':
                continue
            chr2 = line[3]

            if chr1 == 'chrM' or chr2 == 'chrM':
                continue
            count_sds += 1
            s1 = int(line[1])
            e1 = int(line[2])

            s2 = int(line[4])
            e2 = int(line[5])

            # if e1 - s1 > 5000000:
            #     print (j)
            if e1 - s1 > threshold1:
                counters[0] += 1
            if e1 - s1 > threshold2:
                counters[1] += 1
            if e1 - s1 > threshold3:
                counters[2] += 1
            
            chr1 += line[8]
            chr2 += line[9] if line[9] == '+' else '_'

            if chr1 in d:
                d[chr1].add(Interval(s1,e1))
            else:
                d[chr1] = IntervalTree()
                d[chr1].add(Interval(s1,e1))

            if chr2 in d:
                d[chr2].add(Interval(s2,e2))
            else:
                d[chr2] = IntervalTree()
                d[chr2].add(Interval(s2,e2))
    
    coverage = 0
    # coverage_before = 0
    for i in d:
        # for inerv in d[i]:
        #     coverage_before += inerv.end - inerv.begin
        d[i].merge_overlaps()
        for inerv in d[i]:
            coverage += inerv.end - inerv.begin
    km = path.split("/")[1]
    return d, count_sds, coverage


def modify_sequences(main_dict, bed_file):
    from cigar import Cigar
    from Bio.Seq import Seq
    for i in open(bed_file, 'r'):
        line = i.split('\t')
        chr1 = line[0]
        if chr1 == '#chr1':
            continue
        chr2 = line[3]

        if chr1 == 'chrM' or chr2 == 'chrM':
            continue
        s1 = int(line[1])
        e1 = int(line[2])

        s2 = int(line[4])
        e2 = int(line[5])
        cigar = line[12]

        seq1 = main_dict[chr1] [s1:e1]
        seq2 = main_dict[chr2] [s2:e2] if line[9] == '+' else str( Seq(main_dict[chr2] [s2:e2]).reverse_complement() )
        seq_new = ''

        counter_1 = 0
        counter_2 = 0

        for num, let in Cigar(cigar).items():
            # print (let, num)
            if let == 'M':
                seq_new += seq1[counter_1: counter_1+num]
                # if seq1[counter_1: counter_1+num] != seq2[counter_2: counter_2+num]:
                #     print (seq1[counter_1: counter_1+num])
                #     print (seq2[counter_2: counter_2+num])
                counter_1 += num
                counter_2 += num
            elif let == 'D':
                counter_1 += num
            elif let == 'I':
                seq_new += seq2[counter_2: counter_2+num]
                counter_2 += num
        if line[9] == '-':
            seq_new = str ( Seq(seq_new).reverse_complement() )
        main_dict[chr2] = main_dict[chr2][:s2] + seq_new + main_dict[chr2][e2:]
        # print (len(seq_new), len(seq2))
        # print (seq_new)
        # print (seq2)
        assert len(seq_new) == len(seq2)





def chop_this_region(s1,e1,s2,e2, threshold):
    
    while s1 < e1:
        s2_ = s2
        t1 = threshold
        while s2_<e2:
            t2 = threshold
            if s2_ + threshold > e2:
                t2 =e2 - s2_
            if s1 + t1 > e1:
                t1 = e1-s1
            yield ((s1, s1+t1), (s2_, s2_ + t2))
            s2_ += t2
        s1+=t1

# for i in chop_this_region(0,220,300,550, 100):
#     print(i)


def chop_huge_regions():
    path = 'biser_ms'
    threshold = 1000000
    for i in os.listdir(path):
        print (i)
        if os.path.isdir(f'{path}/{i}') and i != 'log' and i != 'aligned' and i.split('_')[0] ==  i.split('_')[1] and i.split('_')[0] == 'mm10':
            
            for files in os.listdir(f'{path}/{i}/merged_ok/'):
                new_sds = []

                for line_ in open(f'{path}/{i}/merged_ok/{files}', 'r'):
                    line = line_.split('\t')
                    chr1 = line[0]
                    chr2 = line[3]

                    if chr1 == 'chrM' or chr2 == 'chrM':
                        continue
                    # count_sds += 1
                    s1 = int(line[1])
                    e1 = int(line[2])

                    s2 = int(line[4])
                    e2 = int(line[5])
                    if e1-s1 > threshold or e2-s2 > threshold:
                        
                        for ((s1, e1),(s2,e2)) in chop_this_region(s1,e1,s2,e2, threshold):
                            new_sds.append(f'{chr1}\t{s1}\t{e1}\t{chr2}\t{s2}\t{e2}\t{line[6]}\t{line[7]}\t{line[8]}\t{line[9]}\t{line[10]}\t{line[11]}\n')

                    else:
                        new_sds.append(f'{chr1}\t{s1}\t{e1}\t{chr2}\t{s2}\t{e2}\t{line[6]}\t{line[7]}\t{line[8]}\t{line[9]}\t{line[10]}\t{line[11]}\n')
                writing_f = open(f'{path}/{i}/merged/{files}', 'w')
                for sd in new_sds:
                    writing_f.write(sd)

# chop_huge_regions()



def chop_refions(path1, path2, path3, out_p):
    #path1 is path to folder where merged seeds are
    # specie1 = path1.split('/')[7].split('_')[0]
    # specie2 = path2.split('/')[7].split('_')[0]
    d = dict()
    threshold = 1# 5000000
    # d[specie1] = build_dict(path1)
    # d[specie2] = build_dict(path2)
    d, count_sds, coverage = build_dict(path1, d)
    d, count_sds1, coverage1 = build_dict(path2, d)
    d3 = dict()

    print ( d.keys(), count_sds, coverage, count_sds1, coverage1 )

    # sys.exit(1)
    new_sds = []
    for i in os.listdir(path3):
        count_of_lines = 0
        # print (i)
        if '_2' in i:
            continue
        new_sds = []
        for j in open(path3 + '/' + i, 'r'):
            # print (j)
            count_of_lines+=1
            line = j.split('\t')
            if len(line) < 4:
                continue
            chr1 = line[0]
            chr2 = line[3]
            if 'chrM' in chr1 or 'chrM' in chr2:
                continue

            

            # specie1 = line[0].split('#')[1]
            # specie2 = line[3].split('#')[1]
            
            s1 = int(line[1])
            e1 = int(line[2])

            s2 = int(line[4])
            e2 = int(line[5])

            # chop only if it is really big
            if e1 - s1 > threshold or e2 - s2 > threshold:
                
                overlaps1 = d[chr1].overlap(s1,e1)

                overlaps2 = d[chr2].overlap(s2,e2)
                # print (f'sum: {e1-s1}, {e2-s2}')
                # sum1 = 0
                # for i in overlaps1:
                #     sum1 += i.end - i.begin
                # print (f'sum1: {sum1}')
            
                # sum2 = 0
                # for i in overlaps2:
                #     sum2 += i.end - i.begin
                # print (f'sum2: {sum2}')
                


                for ov1 in overlaps1:
                    for ov2 in overlaps2:
                        b1 = max(s1, ov1.begin)
                        e1_ = min(e1, ov1.end)

                        b2 = max(s2, ov2.begin)
                        e2_ = min(e2, ov2.end)
                        # print (f'Here now: {b1}, {e1_},,, {b2}, {e2_}')
                        new_sds.append(f'{chr1}\t{b1}\t{e1_}\t{chr2}\t{b2}\t{e2_}\t{line[6]}\t{line[7]}\t{line[8]}\t{line[9]}\t{line[10]}\t{line[11]}\n')
            # else:
            #     overlaps1 = d[chr1].overlap(s1,e1)

            #     overlaps2 = d[chr2].overlap(s2,e2)
            #     if len(overlaps1) == 0 or len(overlaps2) == 0:
            #         continue
            #     new_sds.append(f'{chr1}\t{s1}\t{e1}\t{chr2}\t{s2}\t{e2}\t{line[6]}\t{line[7]}\t{line[8]}\t{line[9]}\t{line[10]}\t{line[11]}\n')

        new_file = open(out_p+f'{i}', 'w')
        if count_of_lines != len(new_sds):
            print (i, count_of_lines, len(new_sds) )
        for new_sd in new_sds:
            new_file.write(new_sd)




# def calc_coverage_rq():
    # path = ''


def calculate_all_regions():
    path = 'biser_ms/'
    main_number = 0
    main_coverage = 0
    for folder in os.listdir(path):
        if os.path.isdir(path + folder):
            # print (folder)
            second_count = 0
            second_coverage = 0
            # if folder == 'calJac3_calJac3':
            print (folder)
            if folder != 'log' and folder.split('_')[0] == folder.split('_')[1]:
                for merged in os.listdir(path + folder):
                    
                    if merged == 'aligned':
                        d, c, cover = build_dict(path + folder +'/' + merged + '/', dict())
                        print (folder, c, cover)
                        second_count += c
                        main_number += c
                        second_coverage += cover
                        main_coverage += second_coverage
                        # for files in os.listdir(path + folder +'/' + merged + '/' ):
                        #     # print (files)
                        #     second_count += len(open(path + folder +'/' + merged + '/' + files,'r').readlines())
                        #     main_number += len(open(path + folder +'/' + merged + '/' + files,'r').readlines())
                        # print (folder, c, cover)
    print (main_number, main_coverage)
                        

# calculate_all_regions()

def quick_calculate():
    l1 = 'data/genomes/hg19_hard_50.fa.fai'
    l2 = 'data/genomes/panTro6_hard_50.fa.fai'
    sum1 = 0
    for i in open(l1, 'r'):
        line = i.split('\t')
        sum1 += int(line[1])
    sum2 = 0
    for i in open(l2, 'r'):
        line = i.split('\t')
        sum2 += int(line[1])
    print (sum1, sum2, sum1+sum2)

# quick_calculate()





# chop_refions( 'biser_ms/calJac3_calJac3/aligned', 'biser_ms/ponAbe3_ponAbe3/aligned','biser_ms/calJac3_ponAbe3/merged', 'biser_ms/calJac3_ponAbe3/chopped2/')

# d__, i__, j__ = build_dict('biser_ms/panTro_hg19/aligned/', dict())
# print(i__, j__)
# d__, i__, j__ = build_dict('biser_ms/calJac3_ponAbe3/chopped2/', dict())
# print(i__, j__)




from fasta_reader import read_fasta


# ok first read whole genome and save ti to dictionary
def extract_sequences(path, out_folder, input_beds):
    # path = 'data/genomes/'
    # out_folder = 'sdregions8'
    # input_beds = 'same8'

    for i in os.listdir(path):
        main_fa_dict = dict()

        path_ = i.split('_')
        # print (path_)
        if len(path_) > 2 and path_[1] == 'hard' and path_[2] == '50.fa':
            specie = path_[0] # 'hg19'
            # print (specie)
            # continue
            new_fa = open(f'{out_folder}/SD_regions_{specie}.fa', 'w')
            for item in read_fasta(f"{path}{specie}_hard_50.fa"):
                main_fa_dict[f'{specie}#{item.defline}'] = item.sequence
                # print(item)
            print (len(main_fa_dict))
            
            d__ = dict()
            d__, i__, j__ = build_dict(f'{input_beds}/{specie}_{specie}/aligned/', dict())

            # print (f'Dictionary built {i__}, {j__}')
            for i in d__:
                for interval in d__[i]:
                    new_fa.write(f'>{i}-{interval.begin}-{interval.end}\n{main_fa_dict[i][interval.begin : interval.end]}\n')
            new_fa.close()


# extract_sequences()
from Bio.Seq import Seq
# ok first read whole genome and save ti to dictionary
def extract_sequences_2(path, out_folder, input_beds):
    # path = 'data/genomes/'
    # out_folder = 'sdregions8'
    # input_beds = 'same8'

    #first we create dict for genome(s):
    main_fa_dict = dict()
    if os.path.isdir(path):
        for i in os.listdir(path):
            path_ = i.split('_')
            # print (path_)
            if len(path_) > 2 and path_[1] == 'hard' and path_[2] == '50.fa':
                specie = path_[0] # 'hg19'
                # print (specie)
                # continue
                for item in read_fasta(f"{path}{specie}_hard_50.fa"):
                    main_fa_dict[f'{specie}#{item.defline}'] = item.sequence
                    # print(item)
                # print (len(main_fa_dict))
    else:
        specie = path.split('/')[-1].split('_')[0]
        # print (specie, path)
        for item in read_fasta(f"{path}"):
            # print (f'{specie}#{item.defline}')
            main_fa_dict[f'{specie}#{item.defline}'] = item.sequence
    
    # Now we iterate trought different colors files and 
    count_array = []
    num_ov_neg = 0
    
    for i in  os.listdir(input_beds):
        # print (i)
        # this is for fixing sequences
        # modify_sequences(main_fa_dict, f'{input_beds}/{i}')
        # sys.exit(1)
        new_fa = open(f'{out_folder}/{i}.fa', 'w')
        d__ = dict()
        d__, i__, j__ = build_dict_2(f'{input_beds}/{i}', dict())
        count_array.append((i__, j__))
        # print (f'Dictionary built {i__}, {j__}')
        for chr_ in d__:
            for interval in d__[chr_]:
                if chr_[-1] == '+':
                    new_fa.write(f'>{chr_}-{interval.begin}-{interval.end}\n{main_fa_dict[chr_[:-1]][interval.begin : interval.end]}\n')
                else:
                    num_ov_neg += 1
                    new_fa.write(f'>{chr_}-{interval.begin}-{interval.end}\n{Seq(main_fa_dict[chr_[:-1]][interval.begin : interval.end]).reverse_complement()}\n')


        new_fa.close()
    count_array.sort(reverse = True)
    print (f'Num of rev comp: {num_ov_neg}')
    print ('10 biggest clusters: (#of sequences, coverage)')
    for i in count_array[:10]:
        print (i)


def change_coordinates(path = 'test_out_2', path2 = 'final2.bed'):
    final = open(path2, 'w')
    print (path)
    for i in os.listdir(path):
        for j in open(path + '/' + i , 'r'):
            line = j.split('\t')
            first_c = int( line[0].split('-')[1] )
            chr1 = line[0].split('-')[0].split('#')
            start = first_c + int(line[1])
            end = first_c + int(line[2])
            s = f'{chr1[1]}#{chr1[2]}\t{start}\t{end}\t{line[3]}\t{line[4]}\t{line[5]}\t{line[6]}\t{line[7]}\t{line[8]}\t{line[9]}'
            final.write(s)


def change_in_folders(out):
    for i in os.listdir(out):
        if i != 'log' and os.path.isdir(out + '/' + i):
            change_coordinates(f'{out}/{i}/seeds/', f'{out}/{i}/seeds.bed')

# change_in_folders('different8')


def statistics(path, path2, write_path):
    # path = 'different8'
    specie = []
    dic = dict()

    w_f = open(write_path, 'w')

    for i in os.listdir(path):
        if os.path.isdir(path + '/' + i) and not 'log' in i:
            # so the first one is SD regions and the bigger one is genome
            path_full = f'{path}/{i}/final.bed'
            d, k, c = build_dict(path_full, d = dict())
            sd_region = i.split('_')[0]
            whole_g = i.split('_')[1]
            specie.append(sd_region)
            specie.append(whole_g)
            dic[i] = (k,c)
            w_f.write( f'{i}\t{k}\t{c}\n')
    for i in os.listdir(path2):
        if os.path.isdir(path2 + '/' + i) and not 'log' in i:
            # so the first one is SD regions and the bigger one is genome
            path_full = f'{path2}/{i}/final.bed'
            d, k, c = build_dict(path_full, d = dict())
            sd_region = i.split('_')[0]
            whole_g = i.split('_')[1]
            specie.append(sd_region)
            specie.append(whole_g)
            dic[i] = (k,c)
            w_f.write( f'{i}\t{k}\t{c}\n')
    w_f.close()
    sys.exit(1)
    
    specie = list(set(specie))
    # specie.sort()
    matrix = [[]]
    for i in specie:
        for j in specie:
            sr = f'{i}_{j}'
            if sr in dic:
                matrix.append(str(dict[sr].first))
            else:
                matrix.append('')
    final_string = ''
    for i in range(0, len(matrix)):
        if i == 0:
            for spec in specie:
                final_string += f'\t{spec}'
            final_string += '\n'
        final_string += f'{specie[i]}\t'
        for j in range(0, len(matrix)):
            final_string += f'{matrix[i][j]}\t'
        final_string += f'{specie[i]}\n'
    print (final_string)
        
# statistics()
        
def statistics2(path):
    # path = 'l1'
    specie = []
    dic = dict()

    for i in open(path, 'r'):
        gen, k, c = i.split('\t')
        sd_region = gen.split('_')[0]
        whole_g = gen.split('_')[1]
        specie.append(sd_region)
        specie.append(whole_g)

        dic[f'{sd_region}_{whole_g}'] = (k,c)
        # print( f'{i}\t{k}\t{c}')
    specie = list(set(specie))
    print (specie)
    print (dic)
    # sys.exit(1)
    specie.sort()
    matrix = []
    for i in specie:
        help_row = []
        for j in specie:
            sr = f'{i}_{j}'
            sr3 = f'{j}_{i}'

            if sr in dic:
                help_row.append(str(dic[sr][0]).strip())
            elif sr3 in dic:
                help_row.append(str(dic[sr3][0]).strip())

            else:
                help_row.append('')
        matrix.append(help_row)
    print (matrix)
    final_string = ''
    for i in range(0, len(matrix)):
        print (specie)
        if i == 0:
            for spec in specie:
                final_string += f'\t{spec}'
            final_string += '\n'
        final_string += f'{specie[i]}\t'
        for j in range(0, len(matrix)):
            final_string += f'{matrix[i][j]}\t'
        final_string += f'\n'
    print (final_string)
        
# statistics2()

def extract_all(out1, out2, out_fin):
    out = out1 # 'same8'
    # 'different8/final_all.bed'
    final_one = open(out_fin, 'w')
    for i in os.listdir(out):
        if i != 'log' and os.path.isdir(out + '/' + i):
            # print (f'Here: {i}')
            count1 = 0
            for line in open(f'{out}/{i}/final.bed'):
                final_one.write(line)
                count1+= 1
            # print (count1)
            

    out = out2 # 'different8'
    for i in os.listdir(out):
        if i != 'log' and os.path.isdir(out + '/' + i):
            # print (f'Here: {i}')
            count1 = 0
            for line in open(f'{out}/{i}/final.bed', 'r').readlines():
                final_one.write(line)
                count1+= 1
            # print (count1)

        

# extract_all()
# statistics2()


def main():
    import sys
    if len(sys.argv) > 1:
        if sys.argv[1] == 'extract':
            # def extract_sequences(path, out_folder, input_beds):
            # path = 'data/genomes/'
            # out_folder = 'sdregions8'
            # input_beds = 'same8'
            # print (sys.argv[1] , sys.argv[2], sys.argv[3], sys.argv[4])
            extract_sequences( sys.argv[2], sys.argv[3], sys.argv[4])
        elif sys.argv[1] == 'normalize':
            change_in_folders(sys.argv[2])
        elif sys.argv[1] == 'final':
            extract_all(sys.argv[2], sys.argv[3], sys.argv[4])
        elif sys.argv[1] == 'statistics':
            # statistics(sys.argv[2], sys.argv[3], sys.argv[4])
            statistics2( sys.argv[4])

        elif sys.argv[1] == 'create':
            d, i, s = build_dict(sys.argv[2])
            print (i, s)
        elif sys.argv[1] == 'extract_colors':
            # path to genome(s), out_folder, input_beds
            extract_sequences_2( sys.argv[2], sys.argv[3], sys.argv[4])

main()