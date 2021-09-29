"""Setup and load datasets for FRB analysis
Jessie Thwaites 8/20/21
GFU, v2p6
Define needed functions in multiple files
"""

import numpy as np
import csky as cy
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.time import Time

ana_dir = cy.utils.ensure_dir('/home/jthwaites/csky_cache')
ana = cy.get_analysis(cy.selections.repo, 'version-002-p06', cy.selections.GFUDataSpecs.GFU_IC86, 
                dir=ana_dir)

conf = {'extended': True, #use extended LLH due to low time window
        'space': "ps",
        'time': "transient",
        'sig': 'transient',
        'ana':ana,
        'mp_cpus': 5 #some functions fail with >1 (?)
        }
cy.CONF.update(conf)

#######################################################################
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
        frbs = pd.read_csv('./catalog/spatial_priors_frbs.csv')
    elif 'frb121102' in kw.keys():
        frbs = pd.read_csv('./catalog/frb121102_bursts.csv')
    elif 'others' in kw.keys(): 
        frbs = pd.read_csv('./catalog/frbs_excl_sppriors_121102.csv')
    elif 'all' in kw.keys(): frbs = pd.read_csv('./catalog/frbs_all.csv')
    else: 
        frbs = pd.read_csv('./catalog/frbs_all.csv')
        print('Loading all FRBs from catalogs')    
    return frbs

#######################################################################
#make skymap visualization of loaded frbs
def make_frb_skymap(frbs):
    unique_frbs, frb_ind_all, n_frbs = np.unique(frbs['src'].values, 
                                                 return_index=True, return_counts=True)
    frb_ind_all=frbs['src'].index[frb_ind_all]
        
    unique_frbs_coord = SkyCoord(ra=frbs['ra_deg'].reindex(frb_ind_all).values, 
                      dec=frbs['dec_deg'].reindex(frb_ind_all).values, 
                      unit=(u.hourangle, u.deg))

    ra= unique_frbs_coord.ra.rad
    dec=unique_frbs_coord.dec.rad
    src_map = hl.heal.hist(256, dec, ra)

    fig, ax = plt.subplots(figsize=(12,15),subplot_kw=dict(projection='aitoff'))
    sp = cy.plotting.SkyPlotter(pc_kw=dict(cmap='winter',vmin=0))
    mesh, cb = sp.plot_map(ax, np.where(src_map.map>0, src_map.map, np.nan), n_ticks=2)
    kw = dict(color='.5', alpha=.5)
    sp.plot_gp(ax, lw=.5, **kw)
    sp.plot_gc(ax, **kw)
    ax.grid(**kw)
    plt.title('Skymap of FRB locations \n')
    cb.set_label(r'FRBs/bin')
    
    if print_plot==False: plt.show()
    else: plt.savefig('./plots/frb_skymap.png')
        

