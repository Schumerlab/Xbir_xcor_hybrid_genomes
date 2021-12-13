#! /usr/bin/env python3

import sys
import msprime, pyslim
import random
import csv
import pandas as pd
import pyslim, tskit
import numpy as np
import matplotlib.pyplot as plt

assert len(sys.argv) == 4, 'Incorrect input; please specify "path/to/slim_genetree_to_indiv_ancestries.py input_file.trees number_of_samples output_file.tsv"'

#random.seed(15)

ts = pyslim.load(sys.argv[1]).simplify()
samp_size = int(sys.argv[2])

mixtime = max(ts.individual_times) # Better than specific value for now, but could be generalized more?
was_founder = [x.id for x in ts.nodes() if ((x.population == 1) and (x.time == mixtime))]
today = [x.id for x in ts.nodes() if x.time == 0.0]
today_inds = [ind.id for ind in ts.individuals() if ind.time == 0.0]
samp_inds = random.sample(today_inds, samp_size)
samp_nodes = np.concatenate([ind.nodes for ind in ts.individuals() if ind.id in samp_inds]).tolist()
all_ancestries=[]
for tree in ts.trees():
    leaves_list=[y for z in [tree.leaves(x) for x in was_founder] for y in z] # flattens the nested list of leaves
    ancestries=[1 if x in leaves_list else 0 for x in samp_nodes]
    temp = [sum(ancestries[i:i + 2]) for i in range(0, len(ancestries), 2)]
    outrow = list(tree.interval) + temp
    all_ancestries.append(outrow)

tall_ancestries=list(map(list,zip(*all_ancestries)))
n = len(tall_ancestries) - 2
ancestryDF = pd.DataFrame(all_ancestries, columns=['Start', 'End'] + samp_inds)

ancestryLong = ancestryDF.melt(id_vars=['Start','End'])
long_ancestries = ancestryLong.values.tolist()

clean_ancestries = []
end = long_ancestries[0][1]
rowcount = 0
interval = long_ancestries[0]

while rowcount < len(long_ancestries) - 1:
    nextinterval = long_ancestries[rowcount+1]
    if interval[3] == nextinterval[3] and interval[2] == nextinterval[2]:
        rowcount += 1
    else:
        interval[1] = long_ancestries[rowcount][1]
        clean_ancestries.append(interval)
        interval = nextinterval
        rowcount += 1


interval[1] = long_ancestries[rowcount][1]
clean_ancestries.append(interval)


with open(sys.argv[3], 'w') as file:
    file.writelines('\t'.join(list(map(str,i))) + '\n' for i in clean_ancestries)
