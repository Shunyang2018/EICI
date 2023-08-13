# -*- coding: utf-8 -*-
"""
Created on Mon Sep 26 13:07:30 2022

@author: Study
"""

import pandas as pd 
import re
import os 
from function_mgf import rdmsp


path = './fullpattern/msp/'


for mona in os.listdir(path):
    if '.msp' in mona:
        with open(path + mona,'r', encoding="utf8") as f:
            msps = f.readlines()
        
        
        formula = []
        mass = []    
        inchikeylist = []
        i = 0
        mzslist = []
        for msp in rdmsp(msps):
            if len(msp) != 0:
        
                    
                mzs = ''
                num = msp.index([i for i in msp if re.findall(r'Num Peaks*', i)][0])
                for line in msp:
                    if 'RETENTIONTIME:' in line:
                        formula.append(line.split(':')[-1].strip())
                    if 'INTEGRATEDHEIGHT:' in line:
                        mass.append(line.split(':')[-1].strip())
                    if 'InChiKey' in line:
                        inchikey = line.split(':')[-1].strip()
                        print(inchikey)
                        inchikeylist.append(inchikey)
                for mz in msp[num+1:]:
                    mzs += mz.replace('\t', ':').replace('\n','\t')
                mzslist.append(mzs)
                inchikeylist.append(i)
                i += 1
                
                
                
        df = pd.DataFrame({'name':inchikeylist, 'rt':formula,'height':mass, 'spectrum':mzslist})
        
        df.to_csv(path+ mona.replace('.msp', '.csv'))