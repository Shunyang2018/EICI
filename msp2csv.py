# -*- coding: utf-8 -*-
"""
Created on Mon Sep 26 13:07:30 2022

@author: Study
"""

import pandas as pd 
import re

from function_mgf import rdmsp

mona = './example/MoNA-export-RTX5_Fiehnlib.msp'



with open(mona,'r', encoding="utf8") as f:
    msps = f.readlines()


formula = []
mass = []    
inchikeylist = []

mzslist = []
for msp in rdmsp(msps):
    if len(msp) != 0:

            
        mzs = ''
        num = msp.index([i for i in msp if re.findall(r'Num Peaks*', i)][0])
        for line in msp:
            if 'formula:' in line:
                formula.append(line.split(':')[-1].strip())
            if 'Deriv MW:' in line:
                mass.append(line.split(':')[-1].strip())
            if 'InChiKey' in line:
                inchikey = line.split(':')[-1].strip()
                print(inchikey)
                inchikeylist.append(inchikey)
        for mz in msp[num+1:]:
            mzs += mz.replace('\t', ':').replace('\n','\t')
        mzslist.append(mzs)
        
        
df = pd.DataFrame({'InChIKey':inchikeylist, 'formula':formula,'mass':mass, 'spectrum':mzslist})

df.to_csv(mona.replace('.msp', '.csv'))