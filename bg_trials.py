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
#generate backgroud trials, w/out spatial prior
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

#######################################################################
#generate bg with spatial prior
def sp_scan_bg(src, ntrials=100, logging=False, save_trials=True):
    tr = cy.get_spatial_prior_trial_runner(src_tr=src, 
                  llh_priors=frb_probs, refine_max=False, get_pixmask=True)

    trials=[]
    for i in range(ntrials): trials.append(tr.get_one_trial())
    if logging==True: print('Done running trials')

    sstr = cy.get_sky_scan_trial_runner(nside=256, src_tr=src)
    tss={'pixel':[],'max_ts':[]}
    count=0
    for trial in trials: 
        scan = sstr.get_one_scan_from_trial(trial, mp_cpus=15, logging=False)
        ts_with_prior = scan[1] + tr.llh_prior_term[0]
    
        new_tss=np.maximum(0, ts_with_prior)
        if max(new_tss)!=0.: 
            w=np.where(new_tss==max(new_tss))[0]
            tss['pixel'].append(w)
            tss['max_ts'].append(new_tss[w])
        else: 
            tss['pixel'].append(None)
            tss['max_ts'].append(0.)
        count+=1
        if count%20==0 and logging==True: print('%i/%i scans done'%(count,ntrials))
    
    if save_trials==True: #save trials as csv
        ts_spprior=pd.DataFrame.from_dict(data=tss)
        ts_spprior.to_csv('./bg_sp_trials.csv',index=False)
        
    bg = cy.dists.Chi2TSD(trials)
    return bg