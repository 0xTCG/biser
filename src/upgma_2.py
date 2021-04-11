

from fastcluster import average
import numpy as np
from scipy.cluster.hierarchy import dendrogram
from matplotlib import pyplot as plt
import scipy
import pandas as pd 
import scipy.spatial
import scipy.cluster
import json
import matplotlib.pyplot as plt
from functools import reduce

# X = [[0,0.4,0.5], [0.4,0,0.1],[0.5,0.1,0]]
# labels = ['1__','2__','3_-____---']

import sys
# print(sys.getrecursionlimit())
# sys.exit(1)
sys.setrecursionlimit(200000)



#Loading data

file_ = open('results/output_distances.txt', 'r').readlines()
matrix_len = file_[0]
b = np.ones((int(matrix_len),int(matrix_len)))
for i in file_[1:]:
    line = i.split('\t')
    c1 = int(line[0])
    c2 = int(line[1])
    num = float(line[2])
    b[c1][c2] = num
    b[c2][c1] = num
X = b
names = []
for i in open('results/output_sds.txt', 'r'):
	line = i.split('\t')
	# if (abs(int(line[1])-int(line[2])) < 2):
	# 	continue
	names.append(f'{line[0]}:{line[1]}-{line[2]}&{line[3]}:{line[4]}-{line[5]}')

labels = names


# calculating UPGMA
x = average (X)
file_1 = open('clustered_data2.txt','w')
for i in x:
    file_1.write(f'{int(i[0])}\t{int(i[1])}\t{i[2]}\t{int(i[3])}\n')


print("Done avg")
fig = plt.figure(figsize=(350,120),  dpi=100)
# figsize=(200, 200)
dn = dendrogram(x, labels=labels)
plt.savefig('results/image2.png', bbox_inches='tight')





# JSON approach:
# T = scipy.cluster.hierarchy.to_tree( x )
# # Create a nested dictionary from the ClusterNode's returned by SciPy
# def add_node(node, parent ):
# 	# First create the new node and append it to its parent's children
# 	newNode = dict( node_id=node.id, children=[] )
# 	parent["children"].append( newNode )

# 	# Recursively add the current node's children
# 	if node.left: add_node( node.left, newNode )
# 	if node.right: add_node( node.right, newNode )

# # Initialize nested dictionary for d3, then recursively iterate through tree
# d3Dendro = dict(children=[], name="Root1")
# add_node( T, d3Dendro )


# id2name = dict(zip(range(len(labels)), labels))

# # Label each node with the names of each leaf in its subtree
# def label_tree( n ):
# 	# If the node is a leaf, then we have its name
# 	if len(n["children"]) == 0:
# 		leafNames = [ id2name[n["node_id"]] ]
	
# 	# If not, flatten all the leaves in the node's subtree
# 	else:
# 		leafNames = reduce(lambda ls, c: ls + label_tree(c), n["children"], [])

# 	# Delete the node id since we don't need it anymore and
# 	# it makes for cleaner JSON
# 	del n["node_id"]

# 	# Labeling convention: "-"-separated leaf names
# 	n["name"] = name = "-".join(sorted(map(str, leafNames)))
	
# 	return leafNames

# label_tree( d3Dendro["children"][0] )

# # Output to JSON
# json.dump(d3Dendro, open("d3-dendrogram2.json", "w"), sort_keys=True, indent=4)
# # plt.savefig('figura_main.png', bbox_inches='tight')