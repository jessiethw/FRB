#!/usr/bin/env python

""" Generate bg trials with CHIME spatial maps
    all-sky scans, using csky
    Record TS, ns, and best fit location
    Jessie Thwaites, 1/4/22
"""

%matplotlib inline
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import histlite as hl
import csky as cy
import pandas as pd
from scipy import stats

#making healpy maps
import h5py as h5
import healpy as hp
import pickle as pkl

#import from other files
import setup
import chime_localizations as loc

######################### configure arguments #############################
parser = argparse.ArgumentParser(description='FRB w spatial prior background')
parser.add_argument("--FRB_src", default="FRB20190416A", type=str, 
                    help="FRB source name (tns_name from CHIME) (default=FRB20190416A)")
parser.add_argument("--ntrials", default=1000, type=int,
                    help="Number of trials (default=1000)")
parser.add_argument('--deltaT', type=float, default=86400.,
                    help="Time window in seconds (default=86400s=1d)")
#commented out for now - seed set later w/in code
#parser.add_argument('--seed', type=int, default=123, help="Random number seed")
parser.add_argument('--nside', type=int, default=256, 
                    help="nside to use when making healpix maps")
args = parser.parse_args()
###########################################################################

##Load catalog of FRBs with spatial priors
frbs=setup.load_frbs(spatial_priors=True)
mpl.rcParams['font.size'] = 20
ana=cy.CONF['ana']

index=np.where(frbs['src']==args.FRB_src)
index=index[0][0]
src=setup.sources(args.deltaT, frbs['mjd'].values[index], frbs['ra_deg'].values[index], 
                  frbs['dec_deg'].values[index])

frb_probs, msk=loc.make_healpix_map(args.FRB_src, new_nside=args.nside, max_cl=0.9997)

trials={}
n_trials=args.ntrials
make_plot=True

################# set up trial runners ###########
##sp_tr used for getting llh prior, injecting events
sp_tr = cy.get_spatial_prior_trial_runner(src_tr=src, llh_priors=frb_probs, get_pixmask=True)
##sstr used for each scan
sstr = cy.get_sky_scan_trial_runner(ana=ana, nside=args.nside, 
            src_kw={'mjd':src['mjd'], 't_100':src['t_100'], 'sigma_t':0.}, pixmask=msk)

"""? should this be a random number for the seed, or is this okay?"""
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
            hottest_pix=''
            hottest_loc=''
        
        del ts_with_prior #very large array, delete when no longer needed
    else: 
        scan_max=max(scan[1])
        hottest_pix=''
        hottest_loc=''

    trials[seed]={'max_TS':scan_max, 'hottest_pix':hottest_pix, 'hottest_loc':hottest_loc} 
    del scan #also very large, free memory when unneeded
        
with open('%s_bg_%i.pkl'%(args.FRB_src,args.deltaT), 'wb') as outfile:
    pkl.dump(trials, outfile)
    
"""Histogram of BG TS values"""
fig, ax = plt.subplots(figsize=(9,6))
tss=[]
for i in range(len(trials)):
    tss.append(trials[i]['max_TS'])
        
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
plt.title(r'BG TS distribution, %s, with prior (%is)'%(args.FRB_src,args.deltaT))
ax.legend()

plt.savefig('/home/jthwaites/FRB/plots/bg_ts_distributions/%s_bgts_%is_w_prior.png'%(args.FRB_src,args.deltaT))



