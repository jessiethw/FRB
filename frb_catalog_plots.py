"""
Make plots related to the catalog of FRBs used
Jessie Thwaites, 8/23/21
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl

import catalogs_to_csv as cat
mpl.rcParams['font.size'] = 20

frbs=cat.load_frbs()
rep=np.where(frbs['repeater'])[0]
#non_rep=np.where(frbs['repeater'],False,True)[0]
#print(non_rep)

def repeater_bursts():
    unique_rep, counts=np.unique([frbs['src'][i] for i in rep], return_counts=True)
    
    fig, ax = plt.subplots(figsize=(10,6))
    plt.hist(counts, histtype='step', lw=2., bins=100)
    plt.xlabel('number of bursts')
    plt.ylabel(r'$N$')
    plt.title('Total repeating sources: %i'%len(unique_rep))
    plt.savefig('/home/jthwaites/public_html/frb_param_plots/rep_bursts.png')