""" Script to calculate sensitivies
    for FRBs without spatial priors
    and plot passing fractions 
    Jessie Thwaites, 8/23/21
    NOT COMPLETE
"""

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import histlite as hl
import csky as cy

import bg_trials

##################################################################################
#calculate sensitivity and discovery potential w/out spatial prior
# beta is %, nsigma is # of sigma for dp - defaults are 0.9 (avoid flip-flop) and 5sigma DP
def get_sensitivity(src, beta=0.9, nsigma=5, gamma=2., n_trials=10000, logging=False): 
    
    bg=bg_trials.scan_bg(src, n_trials=n_trials)
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
    else: plt.savefig('/home/jthwaites/public_html/passing_frac/%s_pf_%is_%.1f.png'
                      %(frb_name,int(src['t_100'][0]*84600.),gamma))
    
    print('gamma=%.1f'%gamma)
    print(r'Sensitivity: %.3f +/- %.3f | DP: %.3f +/- %.3f'
          %(sens['n_sig'], sens['n_sig']*sens['n_sig_error'],
            disc['n_sig'], disc['n_sig']*disc['n_sig_error']))
