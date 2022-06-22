"""Setup and load datasets for FRB analysis
Jessie Thwaites 8/20/21
GFU, v2p6
Define needed functions in multiple files
"""

import numpy as np
import csky as cy
import pandas as pd

#define sources for analysis (options: single or list, spatial priors
def sources(t_s, mjd, ra, dec, prior=None, deg=True): #t_s is time window (sec)
    if prior !=None:
        if isinstance(mjd,list):
            src = cy.utils.Sources(
                ra=ra,dec=dec, mjd=[m-(t_s/(84600.*2)) for m in mjd],
                t_100=[t_s/84600.]*len(ra), sigma_t=[0.]*len(ra),
                prior=prior, deg=deg)
        else:
            src = cy.utils.Sources(
                ra=ra,dec=dec, mjd=mjd,
                t_100=[t_s/84600.],sigma_t=[0.],
                prior=prior, deg=deg)
    else: 
        if isinstance(mjd,list):
            src = cy.utils.Sources(
                ra=ra,dec=dec, mjd=[m-(t_s/(84600.*2)) for m in mjd],
                t_100=[t_s/84600.]*len(ra), sigma_t=[0.]*len(ra), deg=deg)
        else:
            src = cy.utils.Sources(
                ra=ra,dec=dec, mjd=mjd,
                t_100=[t_s/84600.],sigma_t=[0.], deg=deg)
    return src

#####################################################################        
#load catalog csv file - default is load all if not specified
def load_frbs(**kw):
    if 'spatial_priors' in kw.keys(): 
        frbs = pd.read_csv('/home/jthwaites/FRB/catalog/spatial_priors_frbs.csv')
    elif 'frb121102' in kw.keys():
        frbs = pd.read_csv('/home/jthwaites/FRB/catalog/frb121102_bursts.csv')
    elif 'all' in kw.keys(): 
        frbs = pd.read_csv('/home/jthwaites/FRB/catalog/frbs_all.csv')
    else: 
        print('Loading FRBs from combined catalogs') 
        frbs = pd.read_csv('/home/jthwaites/FRB/catalog/frbs_all.csv')   
    return frbs