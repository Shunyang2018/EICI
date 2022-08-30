#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import os 
import sys 
import numpy as np
import pandas as pd
from function_mgf import *


    


if len(sys.argv) <2:
    print('ERROR \n Usage: python msp2mgf.py <directory containing msp> <mgf or ms>')
    raise Exception('Missing input option, stop...')
elif len(sys.argv)==2:
    sys.argv.append('mgf')
    
    
# sys.argv = ['.','./test','mgf']
    
path = sys.argv[1]
print(f'working under {path}...')   
    
for file in os.listdir(path):
    if file.endswith('.msp'):
        if sys.argv[2] == 'mgf':
            name, inchikey,formula,ionization = wmgf(path, file)
        else:
            name, inchikey = wms(path, file)
        
        csv = path + '/' + file.replace('msp','csv')
        
        df = pd.DataFrame({'name':name,'inchikey':inchikey,'formula':formula,'ionization':ionization})
        df.to_csv(csv)
        print(f'writing csv to save the identifiers...\n for data safetyS, all identifiers are removed from mgf files')     
        # with open(csv,'w') as f:
        #     f.writelines('Index\tName\tInChIKey\n')
            
        #     for i in range(len(name)):
        #         f.writelines(str(i+1)+'\t')
        #         f.writelines(name[i].strip() + '\t')S
        #         f.writelines(formula[i].strip() + '\t')
        #         f.writelines(inchikey[i].strip() + '\n')
    
            
                

