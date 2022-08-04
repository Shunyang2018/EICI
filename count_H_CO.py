# -*- coding: utf-8 -*-
"""
Created on Mon May 16 23:43:51 2022

@author: Study
"""
from rdkit import Chem
import rdkit.Chem.Fragments
import rdkit.Chem.Lipinski
import pandas as pd
from rdkit.Chem import Descriptors

with open('smi.txt') as f:
    tmp = f.readlines()
    f.close()
#%%
co = []
TMS = []
smilist = []
m = []
for s in tmp:
    smi = s.strip().split()[-1]
    smilist.append(smi)
    mol = Chem.MolFromSmiles(smi)
    co .append( rdkit.Chem.Fragments.fr_C_O_noCOO(mol))
    TMS.append(rdkit.Chem.Lipinski.NumHDonors(mol))
    m.append(round(Descriptors.MolWt(mol),2))


df = pd.DataFrame({'smi':smilist,'MeOx':co,'TMS':TMS,'mass':m})
df.to_csv('smitms.csv')
