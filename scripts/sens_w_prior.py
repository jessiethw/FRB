#!/usr/bin/env python

""" Runs trials to calculate sensitivity with CHIME spatial maps
    This script will inject one value of inject flux and record passing
    fraction. Meant to be run in parallel with more jobs at different values
    of injected flux.
    Jessie Thwaites, 1/4/22
"""  

import numpy as np
import matplotlib as mpl
mpl.use('agg')
#import matplotlib.pyplot as plt
import csky as cy
import h5py as h5
import healpy as hp
import pickle as pkl

import pandas as pd
from scipy import stats
import argparse

import sys
frb_scripts_path='/home/jthwaites/FRB/scripts'
if frb_scripts_path not in sys.path:
    sys.path.append(frb_scripts_path)

#setup analysis object, load analysis
import setup_analysis 
setup_analysis.reload_ana()

import chime_localizations as loc
import general

@np.vectorize
def find_n_sig(dec_deg, sig, bg, beta=0.9, nsigma=None):
    # get signal trials, background distribution, and trial runner
    sig_trials = cy.bk.get_best(sig, 'dec', dec_deg, 'nsig')
    b = cy.bk.get_best(bg, dec_deg)
    tr = cy.bk.get_best(trs, dec_deg)
    # determine ts threshold
    if nsigma is not None:
        ts = b.isf_nsigma(nsigma)
    else:
        ts = b.median()
    # include background trials in calculation
    trials = {0: b.trials}
    trials.update(sig_trials)
    # get number of signal events
    # (arguments prevent additional trials from being run)
    result = sptr.find_n_sig(ts, beta, max_batch_size=0, logging=False, trials=trials, n_bootstrap=1)
    # return flux
    return tr.to_E2dNdE(result, E0=100, unit=1e3)

######################### configure arguments #############################
parser = argparse.ArgumentParser(description='Run background for FRB with spatial prior')
parser.add_argument("--source", default="FRB20190416A", type=str, 
                    help="FRB source name (tns_name from CHIME) (default=FRB20190416A)")
parser.add_argument("--ntrials", default=1000, type=int,
                    help="Number of trials (default=1000)")
parser.add_argument('--deltaT', type=float, default=86400.,
                    help="Time window in seconds (default=86400s=1d)")
#commented out for now - seed set later w/in code
#parser.add_argument('--seed', type=int, default=123, help="Random number seed")
parser.add_argument('--nside', type=int, default=256, 
                    help="nside to use when making healpix maps")
parser.add_argument('--ns_max', type=float, default=5.0,
                    help='max n_signal events for running trials (default=5.0)')
parser.add_argument('--step', type=float, default=0.5,
                    help='step size for generating signal trials (default=0.5)')
args = parser.parse_args()
###########################################################################

##Load catalog of FRBs with spatial priors
frbs=general.load_frbs(spatial_priors=True)
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
#n_trials=args.ntrials

################# set up trial runners ###########
##sptr used for getting llh prior, injecting events
sptr = cy.get_spatial_prior_trial_runner(src_tr=src, llh_priors=frb_probs, get_pixmask=True)
##sstr used for each scan
sstr = cy.get_sky_scan_trial_runner(ana=ana, nside=args.nside, 
            src_kw={'mjd':src['mjd'], 't_100':src['t_100'], 'sigma_t':0.}, pixmask=msk)

##run signal trials
print('Running %i signal trials'%(args.ntrials))
for nsig in np.arange(0., args.ns_max+args.step, args.step):
    print('ns = %s'%(str(nsig)))
    for i in range(args.ntrials):
        trial=sptr.get_one_trial(nsig, poisson=True)
        scan = sstr.get_one_scan_from_trial(trial, mp_cpus=15, logging=False)
        trial_obj=cy.utils.Arrays(init={'mlog10p':scan[0], 'ts':scan[1],
                                        'ns':scan[2], 'gamma':scan[3]})
        if i==0:
            trials[nsig]=[trial_obj]
        else: 
            trials[nsig].append(trial_obj)
print('Starting sensitivity calculation')
result = sptr.find_n_sig(np.median(trials[0]), 0.9, max_batch_size=0, 
                        logging=True, trials=trials, n_bootstrap=1)
print(result)
#trial=sptr.get_one_trial(1)
#scan=sstr.get_one_scan_from_trial(trial, mp_cpus=15, logging=False)
#sstr.find_n_sig(0,.9)
"""
for seed in range(n_trials):
    scan = sstr.get_one_scan(seed = seed, mp_cpus=15, logging=False)
    
    if max(scan[1])>0.:
        ts_with_prior=np.zeros(len(scan[1]))
        w=np.where(scan[1]>0.)[0]
        for i in w:
            pixel_ts=scan[1][i] + sp_tr.llh_prior_term[0][i]
            ts_with_prior[i]=max([pixel_ts,0])
        scan_max=max(ts_with_prior)
        if scan_max>0.: 
            hottest_pix=np.where(ts_with_prior==scan_max)[0][0]
            hottest_loc=hp.pix2ang(args.nside,hottest_pix,lonlat=True)
        else: 
            hottest_pix=None
            hottest_loc=None
        
        del ts_with_prior #very large array, delete when no longer needed
    else: 
        scan_max=max(scan[1])
        hottest_pix=None
        hottest_loc=None

    trials[seed]={'max_TS':scan_max, 'hottest_pix':hottest_pix, 'hottest_loc':hottest_loc} 
    del scan #also very large, free memory when unneeded
print('Done with TS background scans.')    

with open('/home/jthwaites/FRB/background_trials/%s_bg_%i.pkl'%(args.source,args.deltaT), 'wb') as outfile:
    pkl.dump(trials, outfile)
"""