#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  1 12:51:09 2022

@author: shunyang.wang
"""
import re
import pandas as pd
atom_mass = {'H': 1.00782503223,
             'C': 12,
             'N': 14.00307400443,
             'O': 15.99491461957,
             'F': 18.99840316273,
             'P': 30.97376199842,
             'S': 31.9720711744,
             'Cl': 34.968852682,
             'Ar': 39.9623831237,
             'K': 38.9637064864,
             'Ca': 39.962590863,
             'Si': 27.976927,
             'Na':22.98976928
             }

def calculate_mass(formula: str):
    #from yuanyue
    all_atom_nums = re.findall('([A-Z][a-z]*)([0-9]*)', formula)
    mol_mass = 0.
    try:
        for atom_num in all_atom_nums:
            n = atom_num[1]
            if n == '':
                mol_mass += atom_mass[atom_num[0]]
            else:
                mol_mass += int(n) * atom_mass[atom_num[0]]
    except KeyError as e:
        print("Atom {} is not known".format(e.args[0]))
    return round(mol_mass,6)


def rdmsp(lines):
    block = []
    for line in lines:
        if len(line)==1:
            yield block
            block = []
        else:
            block.append(line)

def wmgf(path,file):
    msp = path+'/'+file
    mgf = path + '/' + file.replace('msp','mgf')
    print(f'reading msp: {file} ...' )
    with open(msp,'r') as f:
        tmp = f.readlines()
        f.close()

    print(f'writing mgf: {mgf} ...')
    with open(mgf, 'w') as f:
        i = 0
        name = []
        inchikey = []
        formula = []
        ionization = []
        addmass = []
        for block in rdmsp(tmp):
            f.writelines('BEGIN IONS\n')
            MS = False
            mass0 = '99999'
            for line in block:
                if 'ionization' in line:
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
                        accmass = '0'
                    addmass.append(accmass)
                # if mass0 + '.' in line:
                #     accmass = line.split('\t')[0]

            for line in block:
                if not MS:
                    if 'Comment' in line:
                        pass
                    else:
                        linea, lineb = line.split(':')
                        if linea == 'PRECURSORMZ':
                            f.writelines('PEPMASS='+lineb)
                        elif linea == 'ionization':
                            f.writelines('Ionization='+lineb)
                            ionization.append(lineb.strip())
                        elif linea == 'NAME':
                            i +=1
                            name.append(lineb.strip())
                            f.writelines('NAME='+str(i)+'\n')
                        elif linea == 'InChiKey':
                            inchikey.append(lineb.strip())
                        elif linea == 'formula':
                            formula.append(lineb.strip())
                        elif linea == 'Num Peaks':
                            MS = True
                            f.writelines('PEPMASS='+accmass+'\n')
                            f.writelines('CHARGE=1+\n')
                        # else:
                        #     f.writelines(linea+'='+lineb)

                else:
                    f.writelines(line)


            f.writelines('END IONS \n\n')

    return name, inchikey, formula, ionization,addmass

def wms(path, file):
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

            MS = False

            for line in block:
                if not MS:
                    if 'Comment' in line:
                        pass
                    else:
                        linea, lineb = line.split(':')
                        if linea == 'PRECURSORMZ':
                            f.writelines('>parentmass'+lineb)
                        elif linea == 'Name':
                            i +=1
                            name.append(lineb)
                            f.writelines('>compound'+str(i)+'\n')
                        elif linea == 'InChIKey':
                            inchikey.append(lineb)
                        elif linea == 'Formula':
                            f.writelines('>formula'+ lineb)
                        elif linea == 'Precursor_type':
                            f.writelines('>ionization' + lineb)
                        elif linea == 'Num Peaks':
                            MS = True
                            f.writelines('\n>ms2\n')

                else:
                    f.writelines(line)

    return name, inchikey

