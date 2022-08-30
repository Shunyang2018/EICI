# -*- coding: utf-8 -*-
"""
Created on Fri Feb 18 10:48:08 2022

@author: Study
"""
import os
import numpy as np
import pandas as pd

import math
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
def readiter(file):
    #if length =0, then yield block
    with open(file, 'r') as f:
        block = []
        msp = f.readlines()
        for tmp in msp:
            line = tmp.strip()
            if len(line) != 0:
                block.append(line)
            else:

                yield block
                block = []


def filtermz(mzblock):
    mzlist = []
    valuelist = []
    mz = 0
    value = 0
    for line in mzblock:
        tmp = line.split()
        if mz != float(tmp[0]):
            mzlist.append(mz)
            valuelist.append(value)
            mz = float(tmp[0])
            value = float(tmp[1])
        else:
            mz = float(tmp[0])
            value += float(tmp[1])
    df = pd.DataFrame({'mz':mzlist[1:],'Intensity':valuelist[1:]})

    return df['mz'].array, df['Intensity'].array

def filtermzdf(mzblock):
    mzlist = []
    valuelist = []
    for line in mzblock:
        tmp = line.split()      
        mzlist.append(float(tmp[0]))
        valuelist.append(float(tmp[1]))

    df = pd.DataFrame({'mz':mzlist,'Intensity':valuelist})

    return df



def filtermzdfold(mzblock):
    mzlist = []
    valuelist = []
    mz = 0
    value = 0
    for line in mzblock:
        tmp = line.split()
        if mz != float(tmp[0]):
            mzlist.append(mz)
            valuelist.append(value)
            mz = float(tmp[0])
            value = float(tmp[1])
        else:
            mz = float(tmp[0])
            value += float(tmp[1])
    df = pd.DataFrame({'mz':mzlist[1:],'Intensity':valuelist[1:]})

    return df

def mzdf(mzblock):
    mzlist = []
    valuelist = []

    for line in mzblock:
        tmp = line.split()

        mzlist.append(tmp[0])
        valuelist.append(tmp[1])

    df = pd.DataFrame({'mz':mzlist[1:],'Intensity':valuelist[1:]})

    return df




def filtermz85(mzblock):
    mzlist = []
    valuelist = []
    mz = 0
    value = 0
    for line in mzblock:
        tmp = line.split()
        if mz != float(tmp[0]):
            mzlist.append(mz)
            valuelist.append(value)
            mz = float(tmp[0])
            value = float(tmp[1])
        else:
            mz = float(tmp[0])
            value += float(tmp[1])
    df = pd.DataFrame({'mz':mzlist[1:],'Intensity':valuelist[1:]})
    df = df[df['mz']>=85]
    mz_max = df['Intensity'].max()
    df2 =df [df['Intensity']>mz_max/filtersize]

    return zip(df2['mz'].tolist(), df2['Intensity'].tolist()), int(df2['Intensity'].size)

def filtermz85max(mzblock):
    mzlist = []
    valuelist = []
    mz = 0
    value = 0
    for line in mzblock:
        tmp = line.split()
        if mz != float(tmp[0]):
            mzlist.append(mz)
            valuelist.append(value)
            mz = float(tmp[0])
            value = float(tmp[1])
        else:
            mz = float(tmp[0])
            value += float(tmp[1])
    df = pd.DataFrame({'mz':mzlist[1:],'Intensity':valuelist[1:]})
    df = df[df['mz']>=85]
    mz_max = df['Intensity'].max()
    df2 =df [df['Intensity']>mz_max/filtersize]

    return zip(df2['mz'].tolist(), df2['Intensity'].tolist()), int(df2['Intensity'].size),mz_max

def filtermz85size(mzblock,filtersize):
    mzlist = []
    valuelist = []
    mz = 0
    value = 0
    for line in mzblock:
        tmp = line.split()
        if mz != float(tmp[0]):
            mzlist.append(mz)
            valuelist.append(value)
            mz = float(tmp[0])
            value = float(tmp[1])
        else:
            mz = float(tmp[0])
            value += float(tmp[1])
    df = pd.DataFrame({'mz':mzlist[1:],'Intensity':valuelist[1:]})
    df = df[df['mz']>=85]
    mz_max = df['Intensity'].max()
    if not filtersize:
        df2 = df
    else:
        df2 =df [df['Intensity']>mz_max/filtersize]

    return zip(df2['mz'].tolist(), df2['Intensity'].tolist()), int(df2['Intensity'].size),mz_max




#%% CI detect


def getmz(a):
    return a.mz[a.max_index]

def getint(a):
    return a.Intensity[a.max_index]



def pattern(tmp):
    sudoM = -1
    # find three patterns
    mtype = ''
    cilist = []
    # M-H
    tmplist = np.where(tmp == 14)[0]
    CI = 0
    for i in tmplist:
        if tmp[i-1] == 18:
            CI += 1
            if tmp[i-2] == 12:
                CI += 1
                mtype += 'M-H'
                cilist.append(CI)
                sudoM = i
                if tmp[i-3] == 12:
                    CI += 1


    #M+H
    tmplist = np.where(tmp == 16)[0]
    CI = 0
    for i in tmplist:
        if tmp[i-1] == 16:
            CI += 1
            if tmp[i-2] == 12:
                CI += 1
                mtype += 'M+H'
                cilist.append(CI)
                sudoM = i
                if tmp[i-3] == 12:
                    CI += 1

    #M
    tmplist = np.where(tmp == 15)[0]
    CI = 0
    for i in tmplist:
        if tmp[i-1] == 17:
            CI += 1
            if tmp[i-2] == 12:
                CI += 1
                mtype += 'M'
                cilist.append(CI)
                sudoM = i
                if tmp[i-3] == 12:
                    CI += 1
    return cilist, mtype, sudoM

def pattern2(nonoise):
    tmp = np.array(nonoise['diff'])
    sudoM = []
    # find three patterns
    mtypelist = []
    sudoMtype = []
    # M+H
    tmplist = np.where(tmp == 16)[0]

    for i in tmplist:
        mtype = np.zeros(4) # [0,0,0,0] is [M-15,sudoM,M+17,M+29,M=41] relatively
        mass = nonoise['max_mz'].iloc[i]
        if (32 in abs(nonoise['max_mz'][:i]-mass).round(0).values):
            mtype[1] = 1
        if (44 in abs(nonoise['max_mz'][:i]-mass).round(0).values):
            mtype[2] = 1
        if (56 in abs(nonoise['max_mz'][:i]-mass).round(0).values):
            mtype[3] = 1
        if np.sum(mtype) >0:
            mtypelist.append(mtype)
            sudoM.append(nonoise['max_mz'].iloc[i-1])
            sudoMtype.append('M+H')
        # M-H
    tmplist = np.where(tmp == 14)[0]
    # tmplist = [12]
    for i in tmplist:
        mtype = np.zeros(4) # [0,0,0,0] is [M-15,sudoM,M+17,M+29,M=41] relatively
        mass = nonoise['max_mz'].iloc[i]
        if (32 in abs(nonoise['max_mz'][:i]-mass).round(0).values):
            mtype[1] = 1
        if (44 in abs(nonoise['max_mz'][:i]-mass).round(0).values):
            mtype[2] = 1
        if (56 in abs(nonoise['max_mz'][:i]-mass).round(0).values):
            mtype[3] = 1
        if np.sum(mtype) >0:
            mtypelist.append(mtype)
            sudoM.append(nonoise['max_mz'].iloc[i-1])
            sudoMtype.append('M-H')


        # M
    tmplist = np.where(tmp == 15)[0]
    # tmplist = [12]
    for i in tmplist:
        mtype = np.zeros(4) # [0,0,0,0] is [M-15,sudoM,M+17,M+29,M=41] relatively
        mass = nonoise['max_mz'].iloc[i]
        if (32 in abs(nonoise['max_mz'][:i]-mass).round(0).values):
            mtype[1] = 1
        if (44 in abs(nonoise['max_mz'][:i]-mass).round(0).values):
            mtype[2] = 1
        if (56 in abs(nonoise['max_mz'][:i]-mass).round(0).values):
            mtype[3] = 1
        if np.sum(mtype) >0:
            mtypelist.append(mtype)
            sudoM.append(nonoise['max_mz'].iloc[i-1])
            sudoMtype.append('M')



    return mtypelist,sudoMtype,sudoM

def patternM(nonoise):
    tmp = np.array(nonoise['diff'])
    M = []
    # find three patterns
    mtypelist = []
    sudoMtype = []
    # M+H
    tmplist = np.where(tmp == 16)[0]

    for i in tmplist:
        mtype = np.zeros(4) # [0,0,0,0] is [M-15,sudoM,M+17,M+29,M=41] relatively
        mass = nonoise['max_mz'].iloc[i]
        if (32 in abs(nonoise['max_mz'][:i]-mass).round(0).values):
            mtype[1] = 1
        if (44 in abs(nonoise['max_mz'][:i]-mass).round(0).values):
            mtype[2] = 1
        if (56 in abs(nonoise['max_mz'][:i]-mass).round(0).values):
            mtype[3] = 1
        if np.sum(mtype) >0:
            mtypelist.append(mtype)
            M.append(nonoise['max_mz'].iloc[i-1]-1)
            sudoMtype.append('M+H')
        # M-H
    tmplist = np.where(tmp == 14)[0]
    # tmplist = [12]
    for i in tmplist:
        mtype = np.zeros(4) # [0,0,0,0] is [M-15,sudoM,M+17,M+29,M=41] relatively
        mass = nonoise['max_mz'].iloc[i]
        if (32 in abs(nonoise['max_mz'][:i]-mass).round(0).values):
            mtype[1] = 1
        if (44 in abs(nonoise['max_mz'][:i]-mass).round(0).values):
            mtype[2] = 1
        if (56 in abs(nonoise['max_mz'][:i]-mass).round(0).values):
            mtype[3] = 1
        if np.sum(mtype) >0:
            mtypelist.append(mtype)
            M.append(nonoise['max_mz'].iloc[i-1]+1)
            sudoMtype.append('M-H')


        # M
    tmplist = np.where(tmp == 15)[0]
    # tmplist = [12]
    for i in tmplist:
        mtype = np.zeros(4) # [0,0,0,0] is [M-15,sudoM,M+17,M+29,M=41] relatively
        mass = nonoise['max_mz'].iloc[i]
        if (32 in abs(nonoise['max_mz'][:i]-mass).round(0).values):
            mtype[1] = 1
        if (44 in abs(nonoise['max_mz'][:i]-mass).round(0).values):
            mtype[2] = 1
        if (56 in abs(nonoise['max_mz'][:i]-mass).round(0).values):
            mtype[3] = 1
        if np.sum(mtype) >0:
            mtypelist.append(mtype)
            M.append(nonoise['max_mz'].iloc[i-1])
            sudoMtype.append('M')
    if len(M) > 1:
        tmp = list(map(np.sum,mtypelist))
        m = max(tmp)
        maxindex = [i for i, j in enumerate(tmp) if j == m]
        M = [M[i] for i in maxindex]
        sudoMtype = [sudoMtype[i] for i in maxindex]
        mtypelist = [mtypelist[i] for i in maxindex]
    return mtypelist,sudoMtype,M
def ifciM(block,ms,rt):
    #print all molecular ions
    RT = 'no'
    df = filtermzdf(block[ms:])# mz blocks will vary according to different msp files
    mz, intensity = df['mz'].array, df['Intensity'].array
    df['diff'] = np.append(0,np.diff(mz)).round(0)

    # unique, counts = np.unique(difftest, return_counts=True)

    # find out isotopic peaks and fragments
    df.loc[df['diff']< 1,'iso'] = 0
    df.loc[df['diff'] == 1,'iso']= 1
    df.loc[df['diff'] > 1,'iso']= -1

    iso1 = df[df['iso'] == -1]
    isoindex = df['mz'][df['iso'] == -1].index
    df['part'] = 0
    for i in isoindex:
        df['part'].iloc[:i] +=1


    df['value_grp'] = (df['part'].diff(1) != 0).astype('int').cumsum()

    group = df.groupby('part').agg(list)
    group['len'] = group['mz'].apply(len)
    nonoise = group[group['len']>2]   # noise peak remove
    if nonoise.empty:
        return 'no','no','-1',RT
    else:
        nonoise['max_index'] = nonoise['Intensity'].apply(lambda x: x.index(max(x)))
        nonoise['max_mz'] = nonoise.apply(getmz,axis=1)
        nonoise['max_int'] = nonoise.apply(getint,axis=1)
        nonoise['diff'] = np.append(0,-np.diff(nonoise.max_mz)).round(0)
        # cilist, mtype, sudoM = pattern(np.array(nonoise['diff']))
        cilist, mtype, sudoM = patternM(nonoise)

        RT = block[rt].split(':')[-1].strip()
        return cilist,mtype,sudoM,RT
    
    
def ifci(block,ms,rt):
    RT = 'no'
    df = filtermzdf(block[ms:])# mz blocks will vary according to different msp files
    mz, intensity = df['mz'].array, df['Intensity'].array
    df['diff'] = np.append(0,np.diff(mz)).round(0)

    # unique, counts = np.unique(difftest, return_counts=True)

    # find out isotopic peaks and fragments
    df.loc[df['diff']< 1,'iso'] = 0
    df.loc[df['diff'] == 1,'iso']= 1
    df.loc[df['diff'] > 1,'iso']= -1

    iso1 = df[df['iso'] == -1]
    isoindex = df['mz'][df['iso'] == -1].index
    df['part'] = 0
    for i in isoindex:
        df['part'].iloc[:i] +=1


    df['value_grp'] = (df['part'].diff(1) != 0).astype('int').cumsum()

    group = df.groupby('part').agg(list)
    group['len'] = group['mz'].apply(len)
    nonoise = group[group['len']>2]   # noise peak remove
    if nonoise.empty:
        return 'no','no','-1',RT
    else:
        nonoise['max_index'] = nonoise['Intensity'].apply(lambda x: x.index(max(x)))
        nonoise['max_mz'] = nonoise.apply(getmz,axis=1)
        nonoise['max_int'] = nonoise.apply(getint,axis=1)
        nonoise['diff'] = np.append(0,-np.diff(nonoise.max_mz)).round(0)
        # cilist, mtype, sudoM = pattern(np.array(nonoise['diff']))
        cilist, mtype, sudoM = pattern2(nonoise)

        RT = block[rt].split(':')[-1].strip()
        return cilist,mtype,sudoM,RT
    
    
#%%

def patternMax(nonoise):
    tmp = np.array(nonoise['diff'])
    M = []
    # find three patterns
    mtypelist = []
    sudoMtype = []
    # M+H
    tmplist = np.where(tmp == 16)[0]

    for i in tmplist:
        
        mtype = np.zeros(4) # [0,0,0,0] is [M-15,sudoM,M+17,M+29,M=41] relatively
        if i == 1:
            mtype[0]=1
        mass = nonoise['max_mz'].iloc[i]
        if (32 in abs(nonoise['max_mz'][:i]-mass).round(0).values):
            mtype[1] = 1
        if (44 in abs(nonoise['max_mz'][:i]-mass).round(0).values):
            mtype[2] = 1
        if (56 in abs(nonoise['max_mz'][:i]-mass).round(0).values):
            mtype[3] = 1
        if np.sum(mtype) >0:
            mtypelist.append(mtype)
            M.append(nonoise['max_mz'].iloc[i-1]-1)
            sudoMtype.append('M+H')
        # M-H
    tmplist = np.where(tmp == 14)[0]
    # tmplist = [12]
    for i in tmplist:
        mtype = np.zeros(4) # [0,0,0,0] is [M-15,sudoM,M+17,M+29,M=41] relatively
        if i == 1:
            mtype[0]=1
        mass = nonoise['max_mz'].iloc[i]
        if (32 in abs(nonoise['max_mz'][:i]-mass).round(0).values):
            mtype[1] = 1
        if (44 in abs(nonoise['max_mz'][:i]-mass).round(0).values):
            mtype[2] = 1
        if (56 in abs(nonoise['max_mz'][:i]-mass).round(0).values):
            mtype[3] = 1
        if np.sum(mtype) >0:
            mtypelist.append(mtype)
            M.append(nonoise['max_mz'].iloc[i-1]+1)
            sudoMtype.append('M-H')


        # M
    tmplist = np.where(tmp == 15)[0]
    # tmplist = [12]
    for i in tmplist:
        mtype = np.zeros(4) # [0,0,0,0] is [M-15,sudoM,M+17,M+29,M=41] relatively
        if i == 1:
            mtype[0]=1
        mass = nonoise['max_mz'].iloc[i]
        if (32 in abs(nonoise['max_mz'][:i]-mass).round(0).values):
            mtype[1] = 1
        if (44 in abs(nonoise['max_mz'][:i]-mass).round(0).values):
            mtype[2] = 1
        if (56 in abs(nonoise['max_mz'][:i]-mass).round(0).values):
            mtype[3] = 1
        if np.sum(mtype) >0:
            mtypelist.append(mtype)
            M.append(nonoise['max_mz'].iloc[i-1])
            sudoMtype.append('M')
    if len(M) > 0:
        maxindex = np.argmax(M)
        
        return mtypelist[maxindex],sudoMtype[maxindex],M[maxindex]
    else:
        return 'no','no','-1'
def ifciMax(block,ms,rt):
    #only output the largest molecular ion
    RT = 'no'
    df = filtermzdf(block[ms:])# mz blocks will vary according to different msp files
    mz, intensity = df['mz'].array, df['Intensity'].array
    df['diff'] = np.append(0,np.diff(mz)).round(0)

    # unique, counts = np.unique(difftest, return_counts=True)

    # find out isotopic peaks and fragments
    df.loc[df['diff']< 1,'iso'] = 0
    df.loc[df['diff'] == 1,'iso']= 1
    df.loc[df['diff'] > 1,'iso']= -1

    iso1 = df[df['iso'] == -1]
    isoindex = df['mz'][df['iso'] == -1].index
    df['part'] = 0
    for i in isoindex:
        df['part'].iloc[:i] +=1


    df['value_grp'] = (df['part'].diff(1) != 0).astype('int').cumsum()

    group = df.groupby('part').agg(list)
    group['len'] = group['mz'].apply(len)
    nonoise = group[group['len']>2]   # noise peak remove
    if nonoise.empty:
        return 'no','no','-1',RT
    else:
        nonoise['max_index'] = nonoise['Intensity'].apply(lambda x: x.index(max(x)))
        nonoise['max_mz'] = nonoise.apply(getmz,axis=1)
        nonoise['max_int'] = nonoise.apply(getint,axis=1)
        nonoise['diff'] = np.append(0,-np.diff(nonoise.max_mz)).round(0)
        # cilist, mtype, sudoM = pattern(np.array(nonoise['diff']))
        cilist, mtype, sudoM = patternMax(nonoise)

        RT = block[rt].split(':')[-1].strip()
        return cilist,mtype,sudoM,RT
    
    
#%%

def findtype(mz):
    A = np.array(mz).round()
    ilist = []
    jlist = []
    # do for every element in the list
    for diff in [14,15,16]:
        for i in range(len(A)):

            # check if pair with the given difference `(i, i-diff)` exists
            if A[i] + diff in A:
                print(i,A[i])
                ilist.append(i)
                jlist.append(np.where(A == A[i] + diff)[0].item())
    if len(ilist) != 0:
        print(ilist)
        i = top.Intensity.iloc[ilist].idxmax()
        j = top.Intensity.iloc[jlist].idxmax()
        mtype = A[j] - A[i]
        if mtype == 16:
            return i,j, 'M+H'
        elif mtype == 15:
            return i,j, 'M'
        elif mtype == 14:
            return i,j, 'M-H'
    else:
        return False, 'no'

def findPair(mz,diff):
    A = np.array(mz).round()
    ilist = []
    # do for every element in the list
    for i in range(len(A)):

        # check if pair with the given difference `(i, i-diff)` exists
        if A[i] + diff in A:
            return True
    return False


def patterntopold(df):
    # generate top intensity df
    top = int(df.shape[0]/25)
    top = df.nlargest(top,'Intensity').sort_values(by=['mz'])
    top = top.reset_index()

    # initial
    mtypelist = []
    sudoM = []
    sudoMtype = []
    mtype = np.zeros(4) # [0,0,0,0] is [M-15,sudoM,M+17,M+29,M=41] relatively
    ilist = []

    i,j, t = findtype(top.mz)

    if findPair(top.mz[i:],32):
        mtype[1] = 1
    if findPair(top.mz[i:],44):
        mtype[2] = 1
    if findPair(top.mz[i:],56):
        mtype[3] = 1
    if np.sum(mtype) >1:
        mtypelist.append(mtype)
        sudoM.append(top.mz.iloc[j])#top.mz.iloc[i]+16

        sudoMtype.append(t)

    return mtypelist,sudoMtype,sudoM



def ifcitop(block,ms,rt):
    RT = 'no'
    df = filtermzdf(block[ms:])# mz blocks will vary according to different msp files

    cilist, mtype, sudoM = patterntop(df)

    RT = block[rt].split(':')[-1].strip()
    return cilist,mtype,sudoM,RT


def patterntop(df):
    # generate top intensity df
    top = int(df.shape[0]/25)
    top = df.nlargest(top,'Intensity').sort_values(by=['mz'])
    top = top.reset_index()

    # initial
    mtypelist = []
    sudoM = []
    sudoMtype = []
    mtype = np.zeros(4) # [0,0,0,0] is [M-15,sudoM,M+17,M+29,M=41] relatively
    jlist = []
    mz =top.mz
    A = np.array(mz).round()

    typelist = ['M-H','M','M+H']
    # do for every element in the list
    for diff in [14,15,16]:
        ilist = []

        for i in range(len(A)):

            # check if pair with the given difference `(i, i-diff)` exists
            if A[i] + diff in A:

                ilist.append(i)

        if len(ilist) != 0:

            i = top.Intensity.iloc[ilist].idxmax()
            jlist.append(i)
            if findPair(top.mz[i:],32):
                mtype[1] = 1
            if findPair(top.mz[i:],44):
                mtype[2] = 1
            if findPair(top.mz[i:],56):
                mtype[3] = 1
            if np.sum(mtype) >0:
                mtypelist.append(mtype)
                sudoM.append(top.mz.iloc[i])#top.mz.iloc[i]+16
            else:
                mtypelist.append('no')
                sudoM.append(-1)
    if sum(sudoM)> 0:
        j = top.Intensity.iloc[jlist].idxmax()
        final = jlist.index(j)
        sudoMtype = typelist[final]
        sudoM = sudoM[final]
        mtypelist = mtypelist[final]
        return mtypelist, sudoMtype,sudoM
    else:
        return 'no','no','-1'































