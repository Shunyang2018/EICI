# -*- coding: utf-8 -*-
"""
Created on Fri July 18 15:06:13 2022

@author: Study
can generate annotation.msp and filter.msp
unit debug is testmsp_debug.py
"""
import numpy as np
from function import readiter,filtermzdf,mzdf,ifciMax
import pandas as pd
# file_name = '../result2/false/2JulyFalsefilter.msp'


ms =8 # ms block
rti = 2
mixblock = 6

#%%
#intial
monthlist = ['Sept','Nov','Aug','Feb','Core','July']

for month in monthlist: 
    
    file_name = '../result2/false/2'+month+'Falsefilter.msp'
    alli = 0   
    cii = 0
    cilist = []
    mtypelist = [] 
    mix = []
    sudoM = []
    rtlist = []
    msp = readiter(file_name)
    with open('../result2/false/'+month+'CI.msp','w') as f:
        for block in msp:
            alli += 1
            #print(block[0])#just report running progress
            CI,mtype,M,RT = ifciMax(block,ms,rti)
            tmp = block[mixblock].split('-')[-1].strip()
            mix.append(tmp)# mix information
            sudoM.append(M)
            rtlist.append(RT)
            if len(CI) == 0:
                cilist.append('no')
                mtypelist.append('no')
              
            else:
                cii +=1
                mtypelist.append(mtype)
                cilist.append(CI)
                
                f.writelines(block[0]+tmp+'\n')
                for i in block[1:]:
                    f.writelines(i+'\n')
    
    df = pd.DataFrame({'cilist':cilist,
                       'sudoM':sudoM,
                       'mtype':mtypelist,
                       'mix':mix,
                       'RT':rtlist}) 
    df = df[df['RT']!='no']
    df = df.astype({'mix': 'int','sudoM':'int','RT':'float','mtype':'str'})
    
    df.to_csv('../result2/'+month+'.csv') 
    df.to_pickle('../result2/'+month+'.pkl')
# water csv
    df = pd.read_csv('../result2/'+month+'.csv') 
    data = '../CI2/'+month+'_all.csv'
    std = pd.read_csv(data)
    std.MW = std.MW.round(0)
    std['TMS'] = std['TMS'].fillna(value=0)
    std['MeOx'] = std['MeOx'].fillna(value=0)
    std['meoxM'] = std['MeOx']*29 + std.MW
    std['Estimated RT'] = std['Estimated RT'].fillna(value=0)
    output = np.array([])
    for i in range(18):
        i = i+1
        subdf = df[df['mix']== i]
        substd = std[std['mix'] == i]
        #test
        # i=5
        # output = np.array([])
        # substd = std[std['mix'] == i]
        # # substd = substd.iloc[5:,:]
        for index, row in substd.iterrows():        
            mass = row.meoxM + 72*np.arange(row.TMS + 1,step=1)
            water = mass -18
            rt = row['Estimated RT']
            if rt != 0:
                rtsubdf = subdf[subdf['RT'].between(rt-4 , rt+4)]
            else:
                rtsubdf = subdf
            
            for j in mass:
                tmp = rtsubdf[rtsubdf.sudoM == j]
                
                if tmp.size > 0:
                    #print(j,tmp)
                    #print('found')
                    
                    tmp['Name'] = row['Metabolite Name']
                    tmp['Estimated RT'] = rt
                    tmp['TMS'] = np.where(mass == j)[0].item()
                    tmp['inchikey'] = row['InChIKey']
                    tmp['MeOx'] = row['MeOx']
                    output = np.append(output,tmp)
               
    
            for j in water:
                tmp = rtsubdf[rtsubdf.sudoM == j]
                if tmp.size > 0:
                    
                    tmp['Name'] = 'water-' +row['Metabolite Name']
                    tmp['Estimated RT'] = rt
                    tmp['TMS'] = np.where(water == j)[0].item()
                    tmp['inchikey'] = row['InChIKey']
                    tmp['MeOx'] = row['MeOx']
                    output = np.append(output,tmp)
            
                    
    result = pd.DataFrame(output.reshape(-1,11),columns=['index','cilist','M','mtype','mix','RT','Name','estimate RT','TMS','inchikey','MeOx'])
    result.to_csv('../result2/'+month+'_water.csv')
        




