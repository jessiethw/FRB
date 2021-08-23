""" Script to generate bg trials
    for FRBs without spatial priors
    Jessie Thwaites, 8/23/21
    NOT COMPLETE
"""

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import histlite as hl
import csky as cy


#######################################################################
#generate backgroud trials
def scan_bg(src, n_trials=10000, frb_name = ' ', print_plot=False):
    tr=cy.get_trial_runner(cy.CONF,ana=cy.CONF['ana'],src=src)
    
    #running bg trials
    trials=tr.get_many_fits(n_trials)
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
    return bg

