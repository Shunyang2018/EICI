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
        
    return spectral_entropy.clean_spectrum(mzs.to_numpy(dtype=np.float32))


ei = './fullpattern/msp/050322_Jake Pool sample_EI.csv'
ci = './fullpattern/msp/042822_Jake Pool sample.csv'
ei = pd.read_csv(ei,index_col=0)
ci = pd.read_csv(ci,index_col=0)


ei = ei[ei['rt']> 7]
ci = ci[ci['rt']> 7]
#%%
entropy = []
rtlist = []
cirt = []
dot = []
dr = []
ue = []
for index, row in ci.iterrows():
    rt = row.rt
    spec_query  =  specarray(row.spectrum)

    eilist = ei[(ei['rt']< rt + 0.5)& (ei['rt'] > rt - 0.5)]
    for index,rowei in eilist.iterrows():

            
        spec_reference =  specarray(rowei.spectrum)
        entropy.append(spectral_entropy.similarity(spec_query, spec_reference, method='entropy',
                                             ms2_da=1))
        dot.append(spectral_entropy.similarity(spec_query, spec_reference, method='dot_product',
                                             ms2_da=1))
        dr.append(spectral_entropy.similarity(spec_query, spec_reference, method='dot_product_reverse',
                                             ms2_da=1))
        ue.append(spectral_entropy.similarity(spec_query, spec_reference, method='unweighted_entropy',
                                             ms2_da=1))
        rtlist.append(rowei.rt)
        cirt.append(rt)
   

df = pd.DataFrame({'cirt':cirt,'eirt':rtlist,'entropy':entropy,
                   'dot_product':dot, 'dot_product_reverse': dr, 'unweighted_entropy':ue})

df['mean'] = df.mean(axis=1)

df.sort_values(['cirt','mean'],ascending=False).groupby(by='cirt').head(3)

df.to_csv('jack.csv')
























