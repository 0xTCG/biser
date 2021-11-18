# from random import *
import random

random.seed(10)

maxSED = 15 #maximum SNP/single indel error (perecnt)
maxLED = 15 #Maximum gap error (percent)

min_gap_size = 300


def get_random_sequence(len_):
    return ''.join([ random.choice("ACTG") for _ in range(len_) ])

def delete_(s, delta):
    how_many_mutateable = int(len(s) * delta)

    index = random.randint(how_many_mutateable, len(s))
    s = s[:index - how_many_mutateable] + s[index:]
    return s

def insert_(s, delta):
    how_many_mutateable = int(len(s) * delta)

    index = random.randint(how_many_mutateable, len(s))
    
    s = s[:index ] + get_random_sequence(how_many_mutateable) + s[index:]
    return s

def insert_chunk(s, insert_seq_len):
    index = random.randint(0, len(s))
    s = s[:index ] + get_random_sequence(insert_seq_len) + s[index:]
    return s, index

def delete_chunk(s, del_seq_len):
    index = random.randint(del_seq_len, len(s))
    s = s[:index - del_seq_len] + s[index:]
    return s, index

def large_mutations(s, delta):
    updates = []
    how_many_mutateable = int(len(s) * delta)

    if how_many_mutateable <= min_gap_size:
        if random.randint(0, 4) % 2 == 0:
            s, index = insert_chunk(s, how_many_mutateable)
            updates.append(('I', index, how_many_mutateable))
            return s, updates
        else:
            s, index = delete_chunk(s, how_many_mutateable)
            updates.append(('D', index, how_many_mutateable))

            return s, updates
    
    while how_many_mutateable > min_gap_size:
        chunk_size = random.randint (min_gap_size, how_many_mutateable)
        how_many_mutateable -= chunk_size
        if random.randint(0, 4) % 2 == 0:
            s, index = insert_chunk(s, chunk_size)
            updates.append(('I', index, chunk_size))

        else:
            s, index = delete_chunk(s, chunk_size)
            updates.append(('D', index, chunk_size))


    
    if how_many_mutateable <= min_gap_size:
        if random.randint(0, 4) % 2 == 0:
            s, index = insert_chunk(s, how_many_mutateable)
            updates.append(('I', index, how_many_mutateable))

            return s, updates
        else:
            s, index = delete_chunk(s, how_many_mutateable)

            updates.append(('D', index, how_many_mutateable))

            return s, updates
    return s, updates

# print(large_mutations('AACCAACCGGTT', 0.3))
# print(large_mutations('AACCAACCGGTT', 0.2))

# print(large_mutations('AACCAACCGGTT', 0.1))



def snps(s, delta):
    s2 = list(s)
    d = {'A': 'CTG', 'C': 'ATG', 'T': 'AGC', 'G': 'TAC'}
    how_many_mutateable = int(len(s) * delta)

    # print (s, how_many_mutateable, delta)
    indexes = random.sample(range(0, len(s) - 1), how_many_mutateable)

    for i in indexes:
        s2[i] = d[s[i]][random.randint(0,2)]
    s2_ret = ''.join(s2)
    return s2_ret

# print(snps('AACCAACCGGTT', 0.3))
# print(snps('AACCAACCGGTT', 0.2))

# print(snps('AACCAACCGGTT', 0.1))



def mutate(s, delta):
    # print (f'Initial:\n{s}')
    # print (delta)
    delta_i = int(delta * 100)
    sED = random.randint(max(0, delta_i - maxLED),min(maxSED, delta_i)) / 100
    s = snps(s, sED)
    # print (sED)
    # print (delta / 2)
    # print ( delta - sED)
    # here we do random insertions and deletions
    s, updates = large_mutations(s, delta - sED)
    # if random.randint(0, 4) % 2 == 0:
    #     s = insert_(s, delta - sED)
    # else:
    #     s = delete_(s, delta - sED)
    # print (f'After:\n{s}')

    return s, updates

# print(mutate('AACCAACCGGTT', 0.3))
# print(mutate('AACCAACCGGTT', 0.2))

# print(mutate('AACCAACCGGTT', 0.1))


def get_random_sequence_from_sequence(s, index):
    len_ = random.randint(1000, min(100000, len(s) - index))
    return s[index: index + len_]
    

def get_complex_elems( min_core_len = 100, max_core_len= 2000, out = 'align/results/test_elems/elems_test.fa', out_elems_path = 'align/results/test_elems/elems_test_rez.bed'):
    list_of_elems = [get_random_sequence(random.randint(min_core_len, max_core_len)) for i in range(0,5)]
    # for i in list_of_elems:
    #     print (i, len(i))

    # chrA
    seq1 = list_of_elems[0] + list_of_elems[1] + list_of_elems[2] + list_of_elems[3] + list_of_elems[2] + list_of_elems[4]
    # chrB
    seq2 = list_of_elems[0] + list_of_elems[1] + list_of_elems[2] + list_of_elems[4]
    # chrC
    seq3 = list_of_elems[1] + list_of_elems[2] + list_of_elems[3] + list_of_elems[2]



    out_w = open(out, 'w')
    out_w.write(f'>chrA+-0-{len(seq1)}\n{seq1}\n>chrB+-0-{len(seq2)}\n{seq2}\n>chrC+-0-{len(seq3)}\n{seq3}')
    out_elems = open(out_elems_path, 'w')
    # 1st elem
    out_elems.write(f'chrA\t0\t{len(list_of_elems[0])}\t0\t{len(list_of_elems[0])}\t+\n')
    out_elems.write(f'chrB\t0\t{len(list_of_elems[0])}\t0\t{len(list_of_elems[0])}\t+\n')
    out_elems.write('\n')

    # 2nd elem
    out_elems.write(f'chrA\t{len(list_of_elems[0])}\t{len(list_of_elems[0]) + len(list_of_elems[1])}\t1\t{len(list_of_elems[1])}\t+\n')
    out_elems.write(f'chrB\t{len(list_of_elems[0])}\t{len(list_of_elems[0]) + len(list_of_elems[1])}\t1\t{len(list_of_elems[1])}\t+\n')
    out_elems.write(f'chrC\t0\t{len(list_of_elems[1])}\t1\t+\n')
    out_elems.write('\n')
    # 3rd elem
    s1 = len(list_of_elems[0] + list_of_elems[1])
    e1 = len(list_of_elems[0] + list_of_elems[1] + list_of_elems[2])
    out_elems.write(f'chrA\t{s1}\t{e1}\t2\t{len(list_of_elems[2])}\t+\n')
    out_elems.write(f'chrB\t{s1}\t{e1}\t2\t{len(list_of_elems[2])}\t+\n')
    s1 += len(list_of_elems[2] + list_of_elems[3])
    e1 += len( list_of_elems[3]+ list_of_elems[2])
    out_elems.write(f'chrA\t{s1}\t{e1}\t2\t{len(list_of_elems[2])}\t+\n')
    out_elems.write(f'chrC\t{len(list_of_elems[1])}\t{len(list_of_elems[1]) + len(list_of_elems[2])}\t2\t{len(list_of_elems[2])}\t+\n')
    s1 = len(list_of_elems[1] + list_of_elems[2] + list_of_elems[3])
    e1 = len(list_of_elems[1] + list_of_elems[2] + list_of_elems[3] + list_of_elems[2])
    out_elems.write(f'chrC\t{s1}\t{e1}\t2\t{len(list_of_elems[2])}\t+\n')
    out_elems.write('\n')


    # 4th elem
    s1 = len(list_of_elems[0] + list_of_elems[1] + list_of_elems[2])
    e1 = len(list_of_elems[0] + list_of_elems[1] + list_of_elems[2] + list_of_elems[3] )
    out_elems.write(f'chrA\t{s1}\t{e1}\t3\t{len(list_of_elems[3])}\t+\n')
    s1 = len(list_of_elems[1] + list_of_elems[2])
    e1 = len(list_of_elems[1] + list_of_elems[2] + list_of_elems[3] )
    out_elems.write(f'chrC\t{s1}\t{e1}\t3\t{len(list_of_elems[3])}\t+\n')
    out_elems.write('\n')

    #5th elem
    s1 = len(list_of_elems[0] + list_of_elems[1] + list_of_elems[2] + list_of_elems[3] + list_of_elems[2])
    e1 = len(list_of_elems[0] + list_of_elems[1] + list_of_elems[2] + list_of_elems[3] + list_of_elems[2] + list_of_elems[4] )
    out_elems.write(f'chrA\t{s1}\t{e1}\t4{len(list_of_elems[4])}\t\t+\n')
    s1 = len(list_of_elems[0] + list_of_elems[1] + list_of_elems[2])
    e1 = len(list_of_elems[0] + list_of_elems[1] + list_of_elems[2] + list_of_elems[4] )
    out_elems.write(f'chrB\t{s1}\t{e1}\t4\t{len(list_of_elems[4])}\t+\n')






get_complex_elems()


def main_funct(path, seq = ''):
    # write_file = open('test/simulations.fa', 'w')
    # path =  r'test/simulations/fa_files2/'

    count = 1
    if seq == '':
        for i in random.sample(range(1000,100000), 1000):
            print (count)
            s = get_random_sequence(i)
            write_file1 = open(f'{path}{count}.fa' , 'w')
            write_file1.write(f'>{count}\n{s}\n')
            write_file1.close()

            write_file2 = open(f'{path}{count}_sds.fa' , 'w')

            count2 = 1
            for delta_i in range(1, 31):
                delta = float(delta_i / 100)

                s2 = mutate(s, delta)[0]
                write_file2.write(f'>{count}c{count2}\n{s2}\n')
                count2 += 1
            count += 1
            write_file2.close()
    else:
        # here we take indexes
        for i in random.sample(range(0,len(seq) - 1000), 1000):
            print (count)
            s = get_random_sequence_from_sequence(seq, i)
            write_file1 = open(f'{path}{count}.fa' , 'w')
            write_file1.write(f'>{count}\n{s}\n')
            write_file1.close()

            write_file2 = open(f'{path}{count}_sds.fa' , 'w')

            count2 = 1
            for delta_i in range(1, 31):
                delta = float(delta_i / 100)

                s2 = mutate(s, delta)[0]
                write_file2.write(f'>{count}c{count2}\n{s2}\n')
                count2 += 1
            count += 1
            write_file2.close()


# this is part of the function used to simulate 30K reads
# import fastaparser
# with open("data/genomes/hg19_hard_50.fa") as fasta_file:
#     print ("data/genomes/hg19_hard_50.fa")
#     parser = fastaparser.Reader(fasta_file, parse_method='quick')
#     for seq in parser:
#         # seq is a FastaSequence object
#         # print (seq.header)

#         if seq.header == '>chr1':
#             print (seq.header)
#             main_funct(seq.sequence)
            # print('Sequence:', seq.sequence_as_string())

def generate_sequences(path, chr_ = ''):
    import fastaparser
    if chr_ != '':
        with open("data/genomes/hg19_hard_50.fa") as fasta_file:
            parser = fastaparser.Reader(fasta_file, parse_method='quick')
            for seq in parser:
                if seq.header == chr_:
                    print (seq.header)
                    main_funct(path, seq.sequence)
    else:
        main_funct(path)



class Core:
    def __init__(self, begin, end):
        self.begin = begin
        self.end = end

delta_main = 0.0

def copy(initial_sequence , s1_begin, s1_end, s1_mate_begin):
    
    # print (delta_main)
    s_hold, updates =  mutate(initial_sequence[s1_begin : s1_end] , delta_main)
    return initial_sequence[:s1_mate_begin] + s_hold + initial_sequence[s1_mate_begin:] + get_random_sequence(10000), len(s_hold), updates
    

def mutate_sequence(initial_sequence):
    return mutate(initial_sequence, 0.1)

def fun_s(i):
    return i[1]

def get_core_coordinates( begin, end, updates):
    updates.sort(key = fun_s)
    if begin == 0:
        for i in updates:
            if i[1] > end:
                break
            if i[0] == 'I':
                end += i[2]
            elif i[1] + i[2] < end:
                end -= i[2]
            # if deleted part goes beyond core len, shift it to starting index of that deletion
            elif i[1] + i[2] >= end:
                end = i[1]
    else:
        for i in updates:
            if i[0] == 'I' and i[1] < begin:
                end += i[2]
                begin += i[2]
            elif i[0] == 'I' and i[1] >= begin:
                end += i[2]
            elif i[0] == 'D' and i[1] + i[2] < begin:
                end -= i[2]
                begin -= i[2]
            # if deleted part goes beyond core len, shift it to starting index of that deletion
            elif i[0] == 'D' and i[1]  >= begin:
                end -= i[2]
            elif i[0] == 'D' and i[1]  < begin and i[1] + i[2] > begin:
                begin = i[1]
                end -= i[2]

    return begin, end




def one_iteration(sd_tree,initial_sequence, sd_coordinates,core_coordinates, main_len, min_core_len, max_sd_len, max_core_len, generation, min_sd_len ):


    core_len = random.randint(min_core_len, max_core_len)
    core_begin = random.randint(max_sd_len, main_len - core_len)
    # print (sd_tree)
    while len(sd_tree.overlap ( core_begin - max_core_len - max_sd_len - min_sd_len, core_begin + max_core_len + max_sd_len + min_sd_len ) ) > 0:
        core_len = random.randint(min_core_len, max_core_len)
        core_begin = random.randint(max_sd_len, main_len - core_len)
        if max_sd_len > 2000:
            max_sd_len -= 500
    core_coordinates[generation] = [(core_begin,core_begin + core_len)]
    assert max_sd_len > min_sd_len
    # core = initial_sequence[ core_begin : core_begin + core_len ]

    # that core has to be created as overlab between 2 sequences; so we first create overlaps and calculate their coordinates
    sd_len = random.randint(max(min_sd_len,core_len + min_core_len), core_len + min_core_len + max_sd_len)
    s1_begin = core_begin
    s1_end = core_begin  + sd_len


    s2_begin =  core_begin + core_len - sd_len
    s2_end = core_begin + core_len

    # now we define coordinates of copies od previously defined SDs.
    s1_mate_begin = main_len # random(0, main_len - max_sd_len)
    
    
    initial_sequence, len_1, updates = copy(initial_sequence , s1_begin, s1_end, s1_mate_begin)
    s1_mate_end = s1_mate_begin + len_1
    main_len = len(initial_sequence)

    b1, e1 = get_core_coordinates( 0, core_len, updates)
    core_coordinates[generation].append((s1_mate_begin + b1, s1_mate_begin + e1))



    s2_mate_begin = main_len # random(0, main_len - max_sd_len)
    

    initial_sequence, len_2, updates = copy(initial_sequence , s2_begin, s2_end, s2_mate_begin)
    s2_mate_end = s2_mate_begin + len_2
    main_len = len(initial_sequence)

    b1, e1 = get_core_coordinates( len_2 - core_len, len_2, updates)
    core_coordinates[generation].append(( s2_mate_end - core_len , s2_mate_end))


    # here we now just update tree and dict of SDs
    if (sd_tree.overlaps(s1_begin, s1_end) or sd_tree.overlaps(s1_mate_begin, s1_mate_end) or sd_tree.overlaps(s2_begin, s2_end) or sd_tree.overlaps(s2_mate_begin, s2_mate_end)):
        print (s1_begin, s1_end)
        print (sd_tree)
        print (sd_tree.overlaps(s1_begin, s1_end), sd_tree.overlaps(s1_mate_begin, s1_mate_end), sd_tree.overlaps(s2_begin, s2_end) , sd_tree.overlaps(s2_mate_begin, s2_mate_end))
        assert False
    sd_tree[s1_begin: s1_end] = True
    sd_tree[s1_mate_begin: s1_mate_end] = True

    sd_tree[s2_begin: s2_end] = True
    sd_tree[s2_mate_begin: s2_mate_end] = True

    sd_coordinates[generation] = [((s1_begin, s1_end), (s1_mate_begin, s1_mate_end)), ((s2_begin, s2_end), (s2_mate_begin, s2_mate_end))]


    # initial_sequence = mutate_sequence(initial_sequence)
    return initial_sequence, main_len







from intervaltree import *

def main_funct_2(main_len, num_of_generations, seed, out_folder = 'test/simulation/', max_error =0.25):
    random.seed(seed)
    global delta_main
    delta_main = 0.0
    print (main_len, num_of_generations, seed, out_folder)
    # max_core_len = 5 # 3000
    # max_sd_len = 20 # 10000
    # main_len = 200 # 1000000
    # min_core_len = 3 # 300
    
    max_core_len = 2000
    max_sd_len = 10000
    min_sd_len = 1000
    # print (max_sd_len)
    # main_len = 100000000 # 100000, 500000, 1000000
    min_core_len = 300
    
    # max_error = 0.3
    # this was first way

    # write_file = open(f'{out_folder}simulations#{int(main_len /100000)}#{num_of_generations}#{seed}.fa', 'w')
    # sds = open(f'{out_folder}sds#{int(main_len /100000)}#{num_of_generations}#{seed}.bed', 'w')
    # cores = open(f'{out_folder}cores#{int(main_len /100000)}#{num_of_generations}#{seed}.bed', 'w')
    max_error_input = int(max_error*100)
    write_file = open(f'{out_folder}simulations#{max_error_input}#{num_of_generations}#{seed}.fa', 'w')
    sds = open(f'{out_folder}sds#{max_error_input}#{num_of_generations}#{seed}.bed', 'w')
    cores = open(f'{out_folder}cores#{max_error_input}#{num_of_generations}#{seed}.bed', 'w')
    print (sds)
    core_coordinates = dict()
    sd_coordinates = dict()

    sd_tree = IntervalTree()

    count = 1
    # num_of_generations = 10

    increase = max_error / num_of_generations

    
    initial_sequence = get_random_sequence(main_len)
    

    for i in range(0, num_of_generations):
        initial_sequence, main_len = one_iteration(sd_tree,initial_sequence, sd_coordinates, core_coordinates, main_len, min_core_len, max_sd_len, max_core_len, i ,min_sd_len)
        delta_main += increase
        
        # print ("111---")

        # print (initial_sequence)
        # print (sd_tree)
        # print (sd_coordinates)
        # print (core_coordinates)
        # print (main_len)
    write_file.write(f'>sequence1\n{initial_sequence}\n')
    for j in sd_coordinates:
        for i in sd_coordinates[j]:

            sds.write(f'sequence1\t{i[0][0]}\t{i[0][1]}\tsequence1\t{i[1][0]}\t{i[1][1]}\n')
    
    for j in core_coordinates:
        for i in core_coordinates[j]:

            cores.write(f'sequence1\t{i[0]}\t{i[1]}\t{j}\n')


# lens = [100000, 500000, 1000000]
# gerenarions = [10, 20, 50]
# seeds = [200, 2391, 139123, 99384, 78767]


# for i in lens:
#     for j in gerenarions:
#         for k in seeds:

#             main_funct_2(i, j, k)




# import fastaparser
from intervaltree import Interval, IntervalTree
import os
def calculate_coverages(path1, path2, threshold, output):
    d = dict() # str: len // name of sequence: its length
    file_name = path1 # 'simulations/main_fa.fa.fai'
    # here we get length of all sequences
    f = open(file_name, 'r').readlines()
    for i in f:
        line = i.split('\t')
        d[line[0]] = int(line[1])
    directory = path2 #'simulations/aligned/'
    main_d = dict() # this is dict that has error as a key and number of successful coverages as value
    count = 0
    
    for filename in os.listdir(directory):
        # get names of 2 seqs from name and their lens with dict d
        if 'y' in filename:
            continue
        name1 = filename.split('_')[0].split('.')[0]
        name2 = filename.split('_')[1].split('.')[0]

        len1 = d[name1]
        len2 = d[name2]


        error = int(name2.split('c')[1])

        # now check the coverages
        l = dict() # dict, name: intervaltree
        for alignment in open(directory + filename ,'r'):
            count += 1
            # print (count)
            # print (alignment)
            line = alignment.split('\t')
            if not line[0] in l:
                # print (line)
                l[ line[0] ] = IntervalTree(  )
                l[line[0]].add(Interval(int(line[1]), int(line[2])))
            else:

                l[ line[0] ].add(Interval(int(line[1]), int(line[2])))
            
            
            if not line[3] in l:
                l[ line[3] ] = IntervalTree (  )
                l[ line[3] ].add(Interval(int(line[4]), int(line[5])))

                
            else:

                l[ line[3] ].add(Interval(int(line[4]), int(line[5])))
            # print (line)
            # print (l)
        mini_d = {line[0]: 0, line[3]: 0}
        for i in l:
            l[i].merge_overlaps()
            for j in l[i]:
                # print (j)
                mini_d[i] += j.end - j.begin
        # print(mini_d)
        # print (mini_d[line[0]] + mini_d[line[3]] , len1 + len2)
        if (mini_d[line[0]] + mini_d[line[3]] ) / (len1 + len2) > threshold:
            if error in main_d:
                main_d[error] += 1
            else:
                main_d[error] = 1
        # else:
        #     print ((mini_d[line[0]] + mini_d[line[3]] ) / (len1 + len2) )
                    
    
    outpur_file = open(output, 'w')
    outpur_file.write('error;hits\n')
    print (main_d, count)
    for i in main_d:
        outpur_file.write(f'{i};{main_d[i]}\n')

            
                
            



# calculate_coverages()