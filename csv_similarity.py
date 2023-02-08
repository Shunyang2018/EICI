# -*- coding: utf-8 -*-
"""
Created on Mon Sep 26 13:22:35 2022

@author: Study
"""

import numpy as np
import scipy.stats
import spectral_entropy
import pandas as pd

ref = pd.read_csv('mona_overlap.csv')
ci = pd.read_csv('./test/filter_all.csv')


def specarray(t):
    mz = t.split('\t')[:-1]
    mzs = pd.DataFrame([sub.split(':') for sub in mz])
        
    return spectral_entropy.clean_spectrum(mzs.to_numpy(dtype=np.float32))


def specarraymax(t,n=10,ms=False):
    
    mzs = pd.DataFrame(specarray(t))
    mzs[1] = mzs[1].astype(float)
    if ms:
        mzs[0] = mzs[0].astype(float)
        mzs = mzs[mzs[0] < ms]
    mzs = mzs.sort_values(1,ascending=False)
    return mzs[:n+1].to_numpy(dtype=np.float32)
#%%
method = [ 'dot_product_reverse']
for m in method:
    score = []
    for index,row in ci.iterrows():
        inchikey = row.InChIKey
        tmp = ref[ref['InChIKey'] == inchikey]
        if len(tmp) != 0:
            spec_query  =  specarray(row.spectrum)
            spec_reference =  specarray(tmp.spectrum.item())
            score.append(spectral_entropy.similarity(spec_query, spec_reference, method=m,
                                                 ms2_da=0.05))
        else:
            score.append('-1')
        
    
    
    ci[m+'0.05'] = score
    
#%%


method = ['dot_product_reverse']
i = 1
for m in method:
    score = []
    q = []
    r = []
    for index,row in ci.iterrows():
        inchikey = row.InChIKey
        tmp = ref[ref['InChIKey'] == inchikey]
        if len(tmp) != 0:
            q.append(row.spectrum)
            r.append(tmp.spectrum.item())
            spec_query  =  specarraymax(row.spectrum,5,row.mass)
            spec_reference =  specarraymax(tmp.spectrum.item(),5,row.mass)
            score.append(spectral_entropy.similarity(spec_query, spec_reference, method=m,
                                                 ms2_da=1))
        else:
            score.append('-1')
            r.append('none')
    break
    

    #     if i >2 :
    #         break
    # break
    ci[m+'remove_5'] = score
    
#%%
ci['mona_spec'] = r
    
ci.to_csv('ci_ref_score.csv')
