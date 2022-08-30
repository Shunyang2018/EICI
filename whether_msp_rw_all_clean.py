# -*- coding: utf-8 -*-
"""
Created on Fri July 18 15:06:13 2022

@author: Study
use water.csv to generate msp library of all annotated CI spectra
"""
import numpy as np
from function import readiter,filtermzdf,mzdf,ifciMax
import pandas as pd
import re
import pubchempy as pcp
import collections

ms =8 # ms block
rti = 2
mixblock = 6
#  '../result2/false/'+month+'CI.msp'
#%%
#intial
monthlist = ['Sept','Nov','Aug','Feb','Core','July']
with open('../result2/false/annotate.msp','w', encoding='utf-8') as f:
    for month in monthlist: 
        print(month)
        df = pd.read_csv('../result2/'+month+'_water.csv')
        df2 = pd.read_csv('../CI2/'+month+'_all.csv')
        file_name = '../result2/false/2'+month+'Falsefilter.msp'
        msp = readiter(file_name)
        for block in msp:
            
            mix = int(block[mixblock].split('-')[-1])
            # if mix >1:
            #     break
            subdf = df[df.mix == mix]
            rt = float(block[rti].split(':')[-1])

            if rt in subdf.RT.tolist():
                
                
                name = subdf.Name[subdf.RT==rt].iloc[0].replace('water-','(water loss)')
                name0 = name
                f.writelines('NAME:\t'+name+'\n')
                name = name.replace('(water loss)','')
                
                tms = subdf.TMS[subdf.RT==rt].iloc[0]
                meox = int(subdf.MeOx[subdf.RT==rt].iloc[0])
                inchikey = df2.InChIKey[df2['Metabolite Name']==name].iloc[0]
                formula = inchikey2formula(inchikey)
                if formula:
                    if '(water loss)' in name0:
                        formula = water(formula)
                    formula = addTMS(formula,tms)
                    formula = addmeox(formula,meox)
                    f.writelines('formula:\t'+formula+'\n')
                else:
                    f.writelines('formula:\tnan\n')
                f.writelines('Retention time:\t'+str(rt)+'\n')
                f.writelines('TMS:\t'+str(tms)+'\n')
                f.writelines('MeOx:\t'+str(meox)+'\n')
                f.writelines('ionization:\t['+str(subdf.mtype[subdf.RT==rt].iloc[0])+']+\n')
                f.writelines('MW:\t'+str(df2.MW[df2['Metabolite Name']==name].iloc[0])+'\n')
                
                f.writelines('InChiKey:\t'+str(inchikey)+'\n')
                f.writelines('Deriv MW:\t'+str(subdf.M[subdf.RT==rt].iloc[0])+'\n')
                f.writelines(block[6]+'\n')
                f.writelines(block[7]+'\n')
                for l in block[8:]:
                    f.writelines(l+'\n')
                f.writelines('\n')   
        # break
#%%    

def inchikey2formula(inchikey):
    try: 
        result = pcp.get_compounds(inchikey,'inchikey')
        formula = result[0].molecular_formula
        return formula.split('-')[0]
    except:
        return False

    

def formula2dict(formula):
    lst = re.findall(r'[A-Z][a-z]*|\d+', re.sub('[A-Z][a-z]*(?![\da-z])', r'\g<0>1', formula))
    res_dct = {lst[i]: int(lst[i + 1]) for i in range(0, len(lst), 2)}
    return res_dct


def water(formula):
    d1 = collections.Counter(formula2dict(formula))
    d2 = collections.Counter({'H':2,'O':1})
    counter = d1 - d2
    formulanew = ''
    for key,value in counter.items():
        formulanew += key+str(value)
    return formulanew

def addTMS(formula,i):
    dict1 = [formula2dict(formula)]
    tms = formula2dict('C3H8Si1')
    for j in range(i):
        dict1.append(tms)
    counter = collections.Counter()
    for d in dict1: 
        counter.update(d)
    formulanew = ''
    for key,value in counter.items():
        formulanew += key+str(value)
    return formulanew


def addmeox(formula,i):
    dict1 = [formula2dict(formula)]
    tms = formula2dict('CH3N')
    for j in range(i):
        dict1.append(tms)
    counter = collections.Counter()
    for d in dict1: 
        counter.update(d)
    formulanew = ''
    for key,value in counter.items():
        formulanew += key+str(value)
    return formulanew


r = addTMS(formula,0)


    
    
    
    
    
    
    
    
    
    
    
    
    
    
    