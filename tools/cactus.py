#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  5 19:23:54 2022

@author: shunyang
https://stackoverflow.com/questions/54930121/converting-molecule-name-to-smiles
https://cactus.nci.nih.gov/chemical/structure
"""
from urllib.request import urlopen
from urllib.parse import quote
import time
def CIRconvert(ids):
    try:
        url = 'http://cactus.nci.nih.gov/chemical/structure/' + quote(ids) + '/smiles'
        ans = urlopen(url).read().decode('utf8')
        return ans
    except:
        return 'Did not work'

identifiers  = ['3-Methylheptane', 'Aspirin', 'Diethylsulfate', 'Diethyl sulfate', '50-78-2', 'Adamant']
with open('name.txt' ,'r') as f:
    identifiers = f.readlines()
f.close()

out = []
i=0
with open('smi.txt' ,'w') as f:
    for ids in identifiers :
        i +=1
        smi = CIRconvert(ids.strip())
        print(i, smi)
        identifiers = f.writelines(str(i)+ ' ' +smi+'\n')
        out.append(smi)
        time.sleep(0.2)
    

    