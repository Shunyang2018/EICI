# -*- coding: utf-8 -*-
"""
Created on Tue Sep  6 13:23:22 2022

@author: Study
"""
from function_mgf import *

h0 = 1.007825

msp = './annotate_all.msp'
i=0
t=0
with open(msp,'r') as f:
    tmp = f.readlines()
    f.close()

with open('test/filter_all.msp','w') as f: 
    
    for block in rdmsp(tmp):
        t += 1
        MS = False
        for line in block:
            if 'ionization' in line:
                ion = line.split(':')[1].strip()
                if 'M-' in ion:
                    m = -1
                    h = -h0
                elif 'M+' in ion:
                    m = 1
                    h = h0
                elif 'M]' in ion:
                    m=0
                    h=0
                else:
                    print('not found: '+ion )
                    exit()
                    
            if 'Deriv' in line: # find max intensity of int. deriv mass peak as accurate mass
                mass0 = line.split(':')[1].strip()# add adducts mass to derive mass
                mass0 = str(int(mass0)+m)
                ms = pd.DataFrame([sub.split('\t') for sub in block[12:]])
                
                ms[1] = ms[1].astype('int')
                ms['mass']=ms[0].astype('float').round(decimals=0)
                ms = ms[ms['mass'] == int(mass0)]
                try:
                    idmax = ms[1].idxmax()
                    accmass = ms[0][idmax]
                except ValueError:
                    
                    accmass = 0
            if not MS:    
                linea, lineb = line.split(':')
                if linea == 'formula':
                    formula = lineb.strip()
                    mass_calc = calculate_mass(formula)
                elif linea == 'Num Peaks':
                    MS = True
        try:       
            ppm = abs(mass_calc-float(accmass)+h)/mass_calc * 1e6
            
        except:
            ppm = 9999
        print(ppm)
        if ppm <= 20:
            f.writelines(block)
            f.writelines('\n')
            i += 1
        # else:
        #     break
        