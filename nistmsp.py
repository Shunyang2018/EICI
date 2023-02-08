# -*- coding: utf-8 -*-
"""
Created on Sat Sep 10 23:29:48 2022

@author: Study
"""

msp = 'example/mainlib_msms.msp'

with open(msp, 'r') as f:
    tmp = f.readlines()
    

with open('nist_main_ci.msp', 'w') as f:
    for line in tmp:
        if 'CAS' in line:
            pass
        elif 'Synon' in line:
            pass
        elif ';' in line:
            tmpline = line.replace(';\n', '').split(';')
            for l in tmpline: 
                f.writelines(l+ '\n')
        elif 'Comments' in line:
            pass
        else:
            f.writelines(line)
        
            
            
            