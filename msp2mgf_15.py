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

def wms15(path, file):
    msp = path+'/'+file
    ms = path + '/' + file.replace('msp','ms')
    print(f'reading msp: {file} ...' )
    with open(msp,'r') as f:
        tmp = f.readlines()
        f.close()

    print(f'writing ms: {ms} ...')
    with open(ms, 'w') as f:
        name = []
        inchikey = []
        i=0
        
        for block in rdmsp(tmp):
            m15 = True
            MS = False

            for line in block:
                if not MS:
                    if 'Comment' in line:
                        pass
                    else:
                        linea, lineb = line.split(':')
                        if linea == 'PRECURSORMZ':
                            f.writelines('>parentmass'+lineb)
                        elif linea == 'NAME':
                            i +=1
                            name.append(lineb)
                            f.writelines('>compound\t'+str(i)+'\n')
                        elif linea == 'InChiKey':
                            inchikey.append(lineb)
                        elif linea == 'Formula':
                            f.writelines('>formula'+ lineb)
                        elif linea == 'PRECURSORTYPE':

                            f.writelines('>ionization' + lineb)
                            ion = lineb
                            if 'M-' in ion:
                                m = -1
                            elif 'M+' in ion:
                                m = 1
                            elif 'M]' in ion:
                                m = 0
                        elif linea == 'Num Peaks':
                            MS = True
                            f.writelines('\n>ms2\n')
                        elif 'Deriv' in line: # find max intensity of int. deriv mass peak as accurate mass
                            mass0 = line.split(':')[1].strip()# add adducts mass to derive mass
                            #mass0 = str(int(mass0)+m)
                            mass0 = str(int(mass0)-15)# 15 loss
                            ms = pd.DataFrame([sub.split('\t') for sub in block[12:]])
                            ms[1] = ms[1].astype('int')
                            ms['mass']=ms[0].astype('float').round(decimals=0)
                            try:
                                ms0 = ms[ms['mass'] == int(mass0)]
                            
                                idmax = ms0[1].idxmax()
                                accmass = ms0[0][idmax]
                                ms = ms[:][idmax:idmax+3]
    
                                f.writelines('>parentmass\t'+accmass+'\n')
                            except:
                                m15= False
                else:
                    if m15:
                        f.writelines(line)
            if m15:
                f.writelines('\n>ms1\n')
                for index, row in ms.iterrows():
                    f.writelines(row[0]+ '\t' +str( row[1])+ '\n')
    
                f.writelines('\n')


    return name, inchikey


sys.argv = ['.','./test/ms1_15/','ms']

path = sys.argv[1]
print(f'working under {path}...')

for file in os.listdir(path):
    if file.endswith('.msp'):
        if sys.argv[2] == 'mgf':
            name, inchikey, precursormz,formula ,adductformula ,precursor,RT = wmgf_15(path, file)
        else:
            name, inchikey = wms15(path, file)

        csv = path + '/' + file.replace('msp','csv')
        print(f'writing csv to save the identifiers...\n for data safety, all identifiers are removed from mgf files')
        df = pd.DataFrame({'name':name,'inchikey':inchikey,
                            'precursormz':precursormz,'formula':formula ,
                            'adductformula':adductformula , 'precursor':precursor})
        df.to_csv(csv,index=True)
