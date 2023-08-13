# -*- coding: utf-8 -*-
"""
Created on Fri Feb 18 11:11:45 2022

@author: Study
"""

import numpy as np
from function import readiter,filtermz

file_name = './result/1000filter.msp'

msp = readiter(file_name)
MZ = []
INT= []
diff = []
for block in msp:
    mz,intensity = filtermz(block[10:])
    MZ.append(mz)
    INT.append(intensity)
    diff.append(np.append(np.diff(mz),0))



mz = np.concatenate(MZ,axis=0)
intensity = np.concatenate(INT,axis=0)
diff = np.concatenate(diff,axis=0)



#%%

munique, mcounts = np.unique(mz, return_counts=True)
m = np.asarray((munique, mcounts)).T
dunique, dcounts = np.unique(diff, return_counts=True)
d = np.asarray((dunique, dcounts)).T
d = d[1:]

diff1000 = np.around(diff,decimals=4)
dunique, dcounts = np.unique(diff1000, return_counts=True)
d = np.asarray((dunique, dcounts)).T[1:]
#%%

import matplotlib.pyplot as plt


counts, edges = np.histogram(diff, bins= np.arange(0.805, 1.205, 0.005))
plt.stairs(counts, edges, fill=True)
d1 = np.asarray((edges, counts)).T

counts, edges = np.histogram(diff, bins= np.arange(0.805, 10.5, 0.005))
plt.stairs(counts, edges, fill=True)

counts, edges = np.histogram(diff, bins= np.arange(10.005, 41, 0.005))
plt.stairs(counts, edges, fill=True)

counts, edges = np.histogram(diff, bins= np.arange(1., 41, 1))
plt.stairs(counts, edges, fill=True)

counts, edges = np.histogram(diff, bins= np.arange(1.,380, 1))
plt.stairs(counts, edges, fill=True)


