# -*- coding: utf-8 -*-
"""
Created on Fri Oct 28 23:42:17 2022

@author: Study
"""

import os
import sys
import numpy as np
import pandas as pd
from function_mgf import *


msp = './test/filter/filter_all.msp'

with open(msp,'r') as f:
    tmp = f.readlines()
    f.close()
    

i = 0
name = []
inchikey = []
formula = []
ionization = []
addmass = []
int15 = []
intm = []
intp = []
int29 = []
int41 = []
int0 = []
tms = []
addmass = []
for block in rdmsp(tmp):
    MS = False

    mass0 = '99999'
    for line in block:
        if 'PRECURSORTYPE' in line:
            ion = line.split(':')[1].strip()
            if 'M-' in ion:
                m = -1
            elif 'M+' in ion:
                m = 1
            elif 'M]' in ion:
                m=0
            else:
                print('not found: '+ion )
                exit()
            
    # decide M       
    for line in block:

        if 'Deriv' in line: # deriv is the mass without adducts
            # mass spec:
            
            ms = pd.DataFrame([sub.split('\t') for sub in block[12:]])
            ms['mass']=ms[0].astype('float').round(decimals=0)
            ms[1] = ms[1].astype('int')
            mass0 = line.split(':')[1].strip() # deriv M
            derivm.append(mass0)
            # accmass
            accmass = str(int(mass0)+m)
            accms = ms[ms['mass'] == int(accmass)]
            try:
                idmax = accms[1].idxmax()
                accmass = accms[0][idmax]

            except ValueError:
                accmass = '0'
            addmass.append(accmass)
            
            # M-15
            mass15 = str(int(mass0)-15)
            
            try:
                tmp15 = ms[1][ms['mass'] == int(mass15)].max().item()/100
                int15.append( tmp15 )
            except:
                int15.append('no M-15!')
            try:
                int0.append(ms[1][ms['mass'] == int(mass0)].max().item()/tmp15)
            except AttributeError:
                int0.append(0)
            try:
                intm.append(ms[1][ms['mass'] == int(mass0) - 1].max().item()/tmp15)
            except AttributeError:
                intm.append(0)
            try:
                intp.append(ms[1][ms['mass'] == int(mass0) + 1].max().item()/tmp15)
            except AttributeError:
                intp.append(0)
            try:
                int29.append(ms[1][ms['mass'] == int(mass0) + 29].max().item()/tmp15)
            except AttributeError:
                int29.append(0)
            try:
                int41.append(ms[1][ms['mass'] == int(mass0) + 41].max().item()/tmp15)
            except AttributeError:
                int41.append(0)  

        

    # extract other information
    for line in block:
        if not MS:
            if 'Comment' in line:
                pass
            else:
                tmp = line.split(':')
                linea, lineb = tmp[0], tmp[1]
    
                if linea == 'ionization':
                    f.writelines('Ionization='+lineb)
                    # f.writelines('Ionization=\t[M]+\n')
                    ionization.append(lineb.strip())
                elif linea == 'NAME':
                    name.append(lineb.strip())
                elif linea == 'InChiKey':
                    inchikey.append(lineb.strip())
                elif linea == 'formula':
                    formula.append(lineb.strip())
                elif linea == 'TMS':
                    tms.append(lineb.strip())
                elif linea == 'Num Peaks':
                    MS = True
        else:
            pass


df = pd.DataFrame({'name':name, 'inchikey':inchikey,'TMS':tms, 'precursorM':addmass, 
                   'int_M':int0,'int_M+H':intp,'int_M-H':intm,'int_M-15':int15,'int_M+29':int29,'int_M+41':int41})


df.to_csv('filter_all_pattern_analysis.csv')
