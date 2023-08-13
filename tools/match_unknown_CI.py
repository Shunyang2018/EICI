# -*- coding: utf-8 -*-
"""
Created on Tue Sep 27 08:45:24 2022

@author: Study
"""

import pandas as pd
import spectral_entropy
import numpy as np

def specarray(t):
    mz = t.split('\t')[:-1]
    mzs = pd.DataFrame([sub.split(':') for sub in mz])
    mzs[0] = mzs[0].astype('float')
    mzs = mzs[mzs[0]>85]
    return spectral_entropy.clean_spectrum(mzs.to_numpy(dtype=np.float32))


ei = './example/GCMS_GITRACT_SignificantUnknowns_ForID.xlsx'
ci = './jake.csv'
ei = pd.read_excel(ei,index_col=0)
ci = pd.read_csv(ci)



#%%
entropy = []
rtlist = []
cirt = []
dot = []
dr = []
ue = []
ciri = []
for index, row in ci.iterrows():
    rt = row.RICI
    rtci = row.cirt
    spec_query  =  specarray(row.CI_spectrum)

    eilist = ei[(ei['Retention Index (GCMS)']< rt + 2000)& (ei['Retention Index (GCMS)'] > rt - 2000)]
    for index,rowei in eilist.iterrows():

            
        spec_reference =  specarray(rowei['MS spectrum (GCMS)'].replace(' ', '\t'))
        entropy.append(spectral_entropy.similarity(spec_query, spec_reference, method='entropy',
                                             ms2_da=1))
        dot.append(spectral_entropy.similarity(spec_query, spec_reference, method='dot_product',
                                             ms2_da=1))
        dr.append(spectral_entropy.similarity(spec_query, spec_reference, method='dot_product_reverse',
                                             ms2_da=1))
        ue.append(spectral_entropy.similarity(spec_query, spec_reference, method='unweighted_entropy',
                                             ms2_da=1))
        rtlist.append(rowei['Retention Index (GCMS)'])
        ciri.append(rt)
        cirt.append(rtci)
    

df = pd.DataFrame({'cirt':cirt,'ciri':ciri,'eirt':rtlist,'entropy':entropy,
                   'dot_product':dot, 'dot_product_reverse': dr, 'unweighted_entropy':ue})

df['mean'] = df.mean(axis=1)

df.sort_values(['cirt','mean'],ascending=False).groupby(by='cirt').head(3)

df.to_csv('jake_unknown.csv')
























