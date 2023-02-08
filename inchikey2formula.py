# -*- coding: utf-8 -*-
"""
Created on Mon Oct 24 22:32:12 2022

@author: Study
"""

import pubchempy
import time
import pandas
file_input_address = 'inchikey.txt'
file_output_address = file_input_address.replace('txt','csv')
inchikey_file=open(file_input_address,'r')
my_dict={'InChIKey':[],'molecular_formula':[]}
#for line in file, call pubchem
counter=0
for line in inchikey_file:
    #if counter%100 == 0:
    print(counter)
    counter+=1
    time.sleep(0.2)
    temp_inchikey=line.strip()
    #print(temp_inchikey)
    temp_compound=pubchempy.get_compounds(temp_inchikey,'inchikey')


    try:
        temp_compound_dict=temp_compound[0].to_dict()
        my_dict['InChIKey'].append(temp_compound_dict['inchikey'])
        my_dict['molecular_formula'].append(temp_compound_dict['molecular_formula'])
    except IndexError:
        
        my_dict['InChIKey'].append(temp_inchikey)
        my_dict['molecular_formula'].append('null')
    #print(my_dict)
my_panda=pandas.DataFrame.from_dict(data=my_dict,orient='columns')
#print(my_panda)
my_panda.to_csv(file_output_address,sep='¬')