#!/usr/bin/env python

""" Generate bg trials with CHIME spatial maps
    all-sky scans, using csky
    Record TS, ns, and best fit location
    Jessie Thwaites, 1/4/22
"""

import numpy as np
import matplotlib as mpl
mpl.use('agg')
import matplotlib.pyplot as plt
import csky as cy
import pandas as pd
from scipy import stats
import argparse

#making healpy maps
import h5py as h5
import healpy as hp
import pickle as pkl

#import from other files
import sys
frb_scripts_path='/home/jthwaites/FRB/scripts'
if frb_scripts_path not in sys.path:
    sys.path.append(frb_scripts_path)

#setup analysis object, load analysis
import setup_analysis 
setup_analysis.reload_ana()

import chime_localizations as loc
import general

######################### configure arguments #############################
parser = argparse.ArgumentParser(description='Run background for FRB with spatial prior')
parser.add_argument("--source", default="FRB20190416A", type=str, 
                    help="FRB source name (tns_name from CHIME) (default=FRB20190416A)")
parser.add_argument("--ntrials", default=1000, type=int,
                    help="Number of trials (default=1000)")
parser.add_argument('--deltaT', type=float, default=86400.,
                    help="Time window in seconds (default=86400s=1d)")
parser.add_argument('--nside', type=int, default=256, 
                    help="nside to use when making healpix maps")
parser.add_argument('--make_ts_hist', type=bool, default=True,
                    help="Make TS distribution and save as png (default=True)")
args = parser.parse_args()
###########################################################################

##Load catalog of FRBs with spatial priors
frbs=general.load_frbs(spatial_priors=True)
mpl.rcParams['font.size'] = 20
ana=cy.CONF['ana']

index=np.where(frbs['src']==args.source)[0]
if len(index)==0:
    print("No match for source with spatial priors found")
    print("Check source name")
    sys.exit()

frb_probs, msk=loc.make_healpix_map(args.source, new_nside=args.nside, max_cl=0.9997)
src=general.sources(args.deltaT, frbs['mjd'].values[index[0]], 
                    frbs['ra_deg'].values[index[0]], frbs['dec_deg'].values[index[0]], 
                    prior=frb_probs)

trials={}
n_trials=args.ntrials

################# set up trial runners ###########
##sp_tr used for getting llh prior, injecting events
sp_tr = cy.get_spatial_prior_trial_runner(src_tr=src, llh_priors=frb_probs, 
            get_pixmask=True, conf={"extra_keep":["energy"]})
##sstr used for each scan
sstr = cy.get_sky_scan_trial_runner(ana=ana, nside=args.nside, 
            src_kw={'mjd':src['mjd'], 't_100':src['t_100'], 'sigma_t':0.}, pixmask=msk)

trials=[]
for seed in range(n_trials):
    trial = sp_tr.get_one_trial(0., seed=seed)
    scan = sstr.get_one_scan_from_trial(trial, seed=seed, mp_cpus=15, logging=False)
    
    ts_with_prior=[scan[1][i]+sp_tr.llh_prior_term[0][i] for i in range(len(sp_tr.llh_prior_term[0]))]
    trials.append(cy.utils.Arrays(init={'mlog10p':scan[0], 'ts':ts_with_prior,
                                        'ns':scan[2], 'gamma':scan[3]}))
                                        
#    if max(scan[1])>0.:
#        ts_with_prior=np.zeros(len(scan[1]))
#        w=np.where(scan[1]>0.)[0]
#        for i in w:
#            pixel_ts=scan[1][i] + sp_tr.llh_prior_term[0][i]
#            ts_with_prior[i]=max([pixel_ts,0])
#        scan_max=max(ts_with_prior)
#        if scan_max>0.: 
#            hottest_pix=np.where(ts_with_prior==scan_max)[0][0]
#            hottest_loc=hp.pix2ang(args.nside,hottest_pix,lonlat=True)
#        else: 
#            hottest_pix=None
#            hottest_loc=None
#        del ts_with_prior #very large array, delete when no longer needed
#    else: 
#        scan_max=max(scan[1])
#        hottest_pix=None
#        hottest_loc=None
#    trials[seed]={'max_TS':scan_max, 'hottest_pix':hottest_pix, 'hottest_loc':hottest_loc} 
print('Done with TS background scans.')    

with open('/home/jthwaites/FRB/background_trials/spatial_prior/%s_bg_%i.pkl'%(args.source,args.deltaT), 'wb') as outfile:
    pkl.dump(trials, outfile)
    
if args.make_ts_hist:
    """Histogram of BG TS values"""
    import histlite as hl
    
    fig, ax = plt.subplots(figsize=(9,6))
    tss=[trials[i]['max_TS'] for i in range(len(trials))]
        
    bg = cy.dists.Chi2TSD(tss)
    h = bg.get_hist(bins=50)
    hl.plot1d(ax, h, crosses=True, label='%i bg trials'%(bg.n_total))

    # chi2 fit
    x = h.centers[0][1:]
    norm = h.integrate().values #normalization for chi-sq
    ax.semilogy(x, norm * bg.pdf(x), lw=1, 
            label=r'$\chi^2$[%.2f dof, $\eta$=%.3f]'%(bg.ndof, bg.eta))

    ax.set_xlabel(r'TS')
    ax.set_ylabel(r'$N$')
    plt.title(f'BG TS distribution, {args.source}, with prior ({args.deltaT}s)')
    ax.legend()

    plt.savefig('/home/jthwaites/FRB/plots/bg_ts_distributions/spatial_prior/%s_bgts_%is_w_prior.png' \
                %(args.source,args.deltaT))

