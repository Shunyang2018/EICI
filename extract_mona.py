# -*- coding: utf-8 -*-
"""
Created on Sun Sep 25 08:44:28 2022

@author: Study
"""

import pandas as pd 
import re

from function_mgf import rdmsp

mona = './example/MoNA-export-GC-MS_Spectra.msp'

with open('./fullpattern/std.txt','r') as f:
    std = f.read().splitlines()

with open(mona,'r', encoding="utf8") as f:
    msps = f.readlines()


formula = []
mass = []    
inchikeylist = []
RI = []
mzslist = []
ri = False
for msp in rdmsp(msps):
    if len(msp) != 0:
        for line in msp:
            if 'InChIKey' in line:
                inchikey = line.split(':')[-1].strip()
                print(inchikey)
        if inchikey in std:
            mzs = ''
            inchikeylist.append(inchikey)
            num = msp.index([i for i in msp if re.findall(r'Num Peaks*', i)][0])
            for line in msp:
                if 'Formula:' in line:
                    formula.append(line.split(':')[-1].strip())
                if 'ExactMass:' in line:
                    mass.append(line.split(':')[-1].strip())
                if 'Retention_index:' in line:
                    RI.append(line.split(':')[-1].strip())
                    ri = True
            for mz in msp[num+1:]:
                mzs += mz.replace('\n','\t').replace(' ', ':')
            mzslist.append(mzs)
            if not ri:
                RI.append(0)
            ri = False
        
        
df = pd.DataFrame({'InChIKey':inchikeylist, 'formula':formula,'mass':mass, 'spectrum':mzslist, 'RI':RI})


df['RI'] = df['RI'].astype('float')
# merge duplicates

df2 = df.groupby(['InChIKey'])['spectrum'].apply(''.join).reset_index()

df2['formula'] = df['formula']
df2['mass'] = df['mass']
df2['RI'] = df['RI']
df.to_excel('mona_overlap.xlsx')
