#!/usr/bin/env python

""" Point source like FRB script with csky
    Jessie Thwaites, Sept 2022
"""

import numpy as np
import matplotlib as mpl
mpl.use('agg')
import matplotlib.pyplot as plt
import histlite as hl
import csky as cy
from scipy import stats
import pickle as pkl
import argparse
import sys

import general
frb_scripts_path='/home/jthwaites/FRB/scripts'
if frb_scripts_path not in sys.path:
    sys.path.append(frb_scripts_path)

#setup analysis object, load analysis
import setup_analysis 
setup_analysis.reload_ana()

######################### configure arguments #############################
parser = argparse.ArgumentParser(description='Run background for FRB with spatial prior')
parser.add_argument("--source", default="FRB20190416A", type=str, 
                    help="FRB source name (tns_name from CHIME) (default=FRB20190416A)")
parser.add_argument("--ntrials", default=1000, type=int,
                    help="Number of trials (default=1000)")
parser.add_argument('--deltaT', type=float, default=86400.,
                    help="Time window in seconds (default=86400s=1d)")
#parser.add_argument('--seed', type=int, default=0, help="Random number seed")
args = parser.parse_args()

#######################################################################
#generate backgroud trials, w/out spatial prior
def scan_bg(src, n_trials=10000, frb_name = 'test_frb', print_plot=False):
    tr=cy.get_trial_runner(cy.CONF,ana=cy.CONF['ana'],src=src)
    
    deltaT=round(src['t_100'][0]*84600.)
    #running bg trials
    trials=tr.get_many_fits(n_trials, seed=0)
    bg = cy.dists.Chi2TSD(trials)
    if np.count_nonzero(trials.ts)==0:
        print('Warning: no nonzero TS values')
        print_plot=False
    
    if print_plot==True: 
        fig, ax = plt.subplots(figsize=(9,6))
        h = bg.get_hist(bins=50)
        hl.plot1d(ax, h, crosses=True, label='%i bg trials'%(bg.n_total))

        # compare with the chi2 fit:
        x = h.centers[0][1:] #remove zero TS bin from curve: not fitted here
        norm = h.integrate().values #normalization for chi-sq
        ax.semilogy(x, norm * bg.pdf(x), lw=1, 
                    label=r'$\chi^2$[%.2f dof, $\eta$=%.3f]'%(bg.ndof, bg.eta))

        ax.set_xlabel(r'TS')
        ax.set_ylabel(r'$N$')
        if src['t_100'][0]==1.: plt.title(r'BG TS distribution, %s (1d)'%(frb_name))
        else: plt.title(r'BG TS distribution, %s (%is)'%(frb_name,src['t_100'][0]*84600.))
        ax.legend()
        plt.savefig('/home/jthwaites/public_html/Background_TS/%s_bgts_%is.png'
                    %(frb_name,int(src['t_100'][0]*84600.)))
        
    with open('/home/jthwaites/FRB/background_trials/point_source/'+
             f'{frb_name}_bg_{deltaT}.pkl', 'wb') as outfile:
        pkl.dump(trials, outfile)
    return bg

##################################################################################
#calculate sensitivity and discovery potential w/out spatial prior
# beta is %, nsigma is # of sigma for dp - defaults are 0.9 (avoid flip-flop) and 5sigma DP
def get_sensitivity(src, beta=0.9, nsigma=5, gamma=2., n_trials=10000, logging=False): 
    
    bg=scan_bg(src, n_trials=n_trials, frb_name = args.source, print_plot=True)
    tr=cy.get_trial_runner(cy.CONF,ana=cy.CONF['ana'], src=src, 
                           inj_conf={'flux':cy.hyp.PowerLawFlux(gamma)})

    #90% sensitivity
    sens=tr.find_n_sig(bg.median(),0.9, tol=0.03,n_batches=10,n_sig_step=1, logging=logging)
    #discovery potential
    disc = tr.find_n_sig(bg.isf_nsigma(nsigma), beta, tol=0.03,n_batches=10,n_sig_step=1, logging=logging)
    return sens, disc

def plot_passing_fraction(src, sens, disc, nsigma=5, gamma=2., frb_name=' ', 
            print_plot=False, show_chisq=True):
    fig, ax = plt.subplots(figsize=(9,6))
    mpl.rcParams['font.size'] = 15
    
    xs = np.linspace(0., max(sens['info']['n_sigs']), 500)
    xs_dp = np.linspace(0., max(disc['info']['n_sigs']), 500)
    
    chi2cdf = lambda n: stats.chi2.cdf(n, *sens['info']['params'])
    chi2cdf_dp = lambda n: stats.chi2.cdf(n, *disc['info']['params'])
    
    plt.errorbar(sens['info']['n_sigs'], sens['info']['CLs'], 
             yerr=sens['info']['sigmas'], label = 'Sensitivity')
    plt.errorbar(disc['info']['n_sigs'],disc['info']['CLs'], 
             yerr=disc['info']['sigmas'], label = r'%i$\sigma$ Discovery potential'%nsigma)
    
    if show_chisq==True:
        plt.plot(xs, chi2cdf(xs), label = 'sensitivity Chi2CDF fit')
        plt.axvline(sens['info']['n_sig_chi2cdf'], ls='--')
        plt.plot(xs_dp, chi2cdf_dp(xs_dp), label = 'DP Chi2CDF fit')
        plt.axvline(disc['info']['n_sig_chi2cdf'], ls='--')
        
    plt.legend(loc=4)
    plt.xlabel(r'$n_{\mathrm{inj}}$')
    plt.ylabel(r'Fraction TS $>$ threshold')
    
    if int(src['t_100'][0]) == 1: 
        plt.title(r'passing fraction, %s (1d) $\gamma$=%.1f'%(frb_name,gamma))
    else: 
        plt.title(r'passing fraction, %s (%is) $\gamma$=%.1f'%(frb_name,src['t_100'][0]*84600., gamma))
    
    if print_plot==False: plt.show()
    else: plt.savefig('/home/jthwaites/FRB/plots/passing_frac/%s_pf_%is_%.1f.png'
                      %(frb_name,int(src['t_100'][0]*84600.),gamma))
    
    print('gamma=%.1f'%gamma)
    print(r'Sensitivity: %.3f +/- %.3f | DP: %.3f +/- %.3f'
          %(sens['n_sig'], sens['n_sig']*sens['n_sig_error'],
            disc['n_sig'], disc['n_sig']*disc['n_sig_error']))

#######################
#running functions for sens, dp and plotting
frbs=general.load_frbs(all=True)

index=np.where(frbs['src']==args.source)[0]
if len(index)==0:
    print("No match for source found")
    print("Check source name")
    sys.exit()

src=general.sources(args.deltaT, frbs['mjd'].values[index[0]], 
                    frbs['ra_deg'].values[index[0]], frbs['dec_deg'].values[index[0]])
sens, dp=get_sensitivity(src)
plot_passing_fraction(src, sens, dp)

with open('/home/jthwaites/FRB/sens_trials/point_source/'+
             f'{args.source}_bg_{args.deltaT}.pkl', 'wb') as outfile:
    pkl.dump(sens, dp, outfile)