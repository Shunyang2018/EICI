# -*- coding: utf-8 -*-
"""
Created on Fri May 13 14:21:10 2022

@author: Study
clean the code for april run and publication

"""

import numpy as np
from function import readiter,filtermzdf,mzdf,ifciMax,filtermz85size
import pandas as pd
import os
# the msp file from ms-dial
path_name = '../../april/' 
filtersize = False
ms =7 # ms block
rti = 2
mixblock = 6
#%% filter
for file_name in os.listdir(path_name+'msdial'):
    name = file_name.split('.')
    if name[-1] == 'msp':
        NAME = []
        SCANNUMBER = []
        RETENTIONTIME = []
        MZ = []
        MODELION = []
        MODELIONHEIGHT = []
   
     
        print(file_name)            
        msp = readiter(path_name+'msdial/'+file_name)
        with open(path_name + 'result_all/'+'nofilter'+name[0]+'.msp','a') as f:
            for block in msp:                              
                # if (float(block[4].split(':')[-1]) > 80000) or (float(block[7].split(':')[-1]) > 10000000):
                # if (float(block[4].split(':')[-1]) < 80000) and (float(block[7].split(':')[-1]) < 10000000):  
                mzblock,size,out = filtermz85size(block[10:],filtersize)
                NAME.append(block[0])#.split(':')[1]
                SCANNUMBER.append(block[1])  
                RETENTIONTIME.append(block[2])
                MODELION.append(block[3])
                MODELIONHEIGHT.append(block[4])
                f.writelines(block[0]+'\n')
                f.writelines(block[1]+'\n')
                f.writelines(block[2]+'\n')
                f.writelines(block[3]+'\n')
                f.writelines(block[4]+'\n')
                f.writelines('Max Intensity:\t'+str(out)+'\n')
                f.writelines('Num Peaks:\t'+str(size)+'\n')

                
                MZ.append(block[10:])
                for l in mzblock:
                    f.writelines(str(l[0])+'\t'+str(int(l[1]))+'\n')
                f.writelines('\n')
#%% CI
for file_name in os.listdir(path_name + 'result_all/'):
    name = file_name.split('.')
    alli = 0   
    cii = 0
    cilist = []
    mtypelist = [] 
    mix = []
    sudoM = []
    rtlist = []
    msp = readiter(path_name + 'result_all/'+file_name)
    with open(path_name + 'CI/'+ 'CI' + name[0] + '.msp','w') as f:
        for block in msp:
            alli += 1
            print(block[0])#just report running progress
            CI,mtype,M,RT = ifciMax(block,ms,rti)
            tmp = block[mixblock].split('-')[-1].strip()
            
            sudoM.append(M)
            rtlist.append(RT)
            if len(CI) == 0:
                cilist.append('no')
                mtypelist.append('no')
              
            else:
                cii +=1
                mtypelist.append(mtype)
                cilist.append(CI)
                
                f.writelines(block[0]+'\n')
                for i in block[1:]:
                    f.writelines(i+'\n')
    
    df = pd.DataFrame({'cilist':cilist,
                       'sudoM':sudoM,
                       'mtype':mtypelist,
                   
                       'RT':rtlist}) 
    df = df[df['RT']!='no']
    df = df.astype({'sudoM':'int','RT':'float','mtype':'str'})
    
    df.to_csv(path_name + 'CI/'+ 'CI' + name[0] +'.csv') 
    df.to_pickle(path_name + 'CI/'+ 'CI' + name[0] +'.pkl')
    
    #%% annotate
    
data = '../../april/agilent.csv'
std = pd.read_csv(data)
std.MW = std.MW.round(0)
std['TMS'] = std['TMS'].fillna(value=0)
std['MeOx'] = std['MeOx'].fillna(value=0)
std['meoxM'] = std['MeOx']*29 + std.MW
std['Estimated RT'] = std['Estimated RT'].fillna(value=0)    
output = np.array([])
for file_name in os.listdir(path_name + 'shunyangmix/'):
    name = file_name.split('.')
    df = pd.read_csv(path_name + 'shunyangmix/'+ file_name) 
    


    subdf = df#[df['mix']== i]
    substd = std#[std['mix'] == i]

    for index, row in substd.iterrows():        
        mass = row.meoxM + 72*np.arange(row.TMS + 1,step=1)
        water = mass -18
        rt = row['Estimated RT']
        rtsubdf = subdf
        # if rt != 0:
        #     rtsubdf = subdf[subdf['RT'].between(rt-3 , rt+3)]
        # else:
        #     rtsubdf = subdf
        for j in mass:
            tmp = rtsubdf[rtsubdf.sudoM == j]
            
            if tmp.size > 0:
                print(j,tmp)
                print('found')
                
                
                tmp['Name'] = row['Metabolite Name']
                tmp['Estimated RT'] = rt
                tmp['TMS'] = np.where(mass == j)[0].item()
                tmp['mix'] = row['mix']
                tmp['shunyang'] = name[0]
                output = np.append(output,tmp)
                print(tmp)

        for j in water:
            tmp = rtsubdf[rtsubdf.sudoM == j]
            if tmp.size > 0:
                
                
                tmp['Name'] = 'water-' +row['Metabolite Name']
                tmp['Estimated RT'] = rt
                tmp['TMS'] = np.where(water == j)[0].item()
                tmp['mix'] = row['mix']
                tmp['shunyang'] = name[0]
                output = np.append(output,tmp)
                    
result = pd.DataFrame(output.reshape(-1,10),columns=['index','cilist','M','mtype','RT','Name','estimate RT','TMS','mix','shunyang'])
result.to_csv(path_name+'annotate/'+'agilent_annotate_all.csv')
                        
