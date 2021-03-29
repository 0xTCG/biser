import os
import sys
from shutil import copyfile


def schedule(path_to_beds, num_of_bins, output):
    bins = 0
    num_of_files = len([name for name in os.listdir('.') if os.path.isfile(name)])
    files_per_bin = num_of_files / num_of_bins
    d=dict()
    for i in range(num_of_bins):
        os.makedirs(f'{output}/{i}')
        d[i] = 0
    bin = 0
    for i in os.listdir( path_to_beds ):
        # os.rename(f'{path_to_beds}/{i}', f"{output}/{bin}/{i}")
        # print (i[-3:] )
        if i[-3:] == '.fa':
            print i
            copyfile(f'{path_to_beds}/{i}', f"{output}/{bin}/{d[bin]}.fa")
            copyfile(f'{path_to_beds}/{i}.fai', f"{output}/{bin}/{d[bin]}.fa.fai")
            d[bin] += 1
            bin += 1
            if bin == num_of_bins:
                bin = 0

schedule(sys.argv[1], int(sys.argv[2]), sys.argv[3])