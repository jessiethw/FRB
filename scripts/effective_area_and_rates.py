import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.style.use('/home/jthwaites/FRB/scripts/plots_style.mplstyle')

def effective_area_gfu():
    exp15 = np.load('/data/ana/analyses/gfu/version-002-p06/IC86_2015_data.npy')
    mc = np.load('/data/ana/analyses/gfu/version-002-p06/IC86_2011_MC.npy')

    #ow units: 1/(GeV cm^2 sr)
    nbins=28
    trueE_bins=np.logspace(2,9,num=nbins+1)
    bin_corr=(9.-2.)/nbins

    #split north/south at -5 deg
    decmsk1=mc['trueDec']>(-5*np.pi/180)
    decmsk2=mc['trueDec']<=(-5*np.pi/180)

    cuts=[np.pi/2,-5*np.pi/180,-np.pi/2]
    dec_msks=[decmsk1,decmsk2]

    A_eff_w=[[],[]]
    for a in range(len(dec_msks)):
        solid_ang=2*np.pi*(np.sin(cuts[a])-np.sin(cuts[a+1]))
        weight=mc['ow'][dec_msks[a]]/(mc['trueE'][dec_msks[a]]*solid_ang*bin_corr*1e4*np.log(10))
        for w in weight: A_eff_w[a].append(w)

    #digitized data from Sam's paper (L2 paper, transient tracks)
    sam_effA_n = pd.read_csv('/home/jthwaites/FRB/catalog/sam_nsky.csv', header=None, names=['Enu','Aeff'])
    sam_effA_s = pd.read_csv('/home/jthwaites/FRB/catalog/sam_ssky.csv', header=None, names=['Enu','Aeff'])

    fig, ax = plt.subplots(figsize=(12,8))
    mpl.rcParams['font.size'] = 20

    plt.plot(sam_effA_n.Enu, sam_effA_n.Aeff, label=r'transient tracks, $\delta>-5$')
    plt.plot(sam_effA_s.Enu, sam_effA_s.Aeff,color='red', label=r'transient tracks, $\delta<-5$')
    plt.hist(mc['trueE'][dec_msks[0]], bins=trueE_bins, weights=A_eff_w[0],
            histtype='step', lw=2., label=r'GFU , $\delta>-5$')
    plt.hist(mc['trueE'][dec_msks[1]], bins=trueE_bins, weights=A_eff_w[1],
            histtype='step', lw=2., label=r'GFU, $\delta<-5$')

    plt.semilogx()
    plt.semilogy()
    plt.legend(loc=0, fontsize='small')
    plt.xlabel(r'$E_\nu$ [GeV]')
    plt.ylabel(r'$A_{eff}$ [$m^2$]') #r means latex only input
    ax.set_title(r'Effective area')

    plt.savefig('/home/jthwaites/FRB/plots/effective_areas_northsouth.png')

    #make finer cuts on dec
    decmsk1=mc['trueDec']>(-5*np.pi/180)
    decmsk2=mc['trueDec']<(-5*np.pi/180)
    decmsk3=mc['trueDec']<(-30*np.pi/180)
    decmsk4=mc['trueDec']>(-30*np.pi/180)
    decmsk5=mc['trueDec']<(30*np.pi/180)
    decmsk6=mc['trueDec']>(30*np.pi/180)

    cuts2=[np.pi/2,30*np.pi/180,-5*np.pi/180,-30*np.pi/180,-np.pi/2]
    dec_msks2=[decmsk3,decmsk4&decmsk2,decmsk1&decmsk5,decmsk6]

    A_eff_w2=[[],[],[],[]]
    for a in range(len(dec_msks2)):
        solid_ang=2*np.pi*(np.sin(cuts2[a])-np.sin(cuts2[a+1]))
        weight=mc['ow'][dec_msks2[a]]/(mc['trueE'][dec_msks2[a]]*solid_ang*bin_corr*1e4*np.log(10))
        for w in weight: A_eff_w2[a].append(w)

    fig, ax = plt.subplots(figsize=(12,8))
    plt.hist(mc['trueE'][dec_msks2[0]], bins=trueE_bins, weights=A_eff_w2[0],
            histtype='step', lw=2., label=r'$-90<\delta<-30$')
    plt.hist(mc['trueE'][dec_msks2[1]], bins=trueE_bins, weights=A_eff_w2[1],
            histtype='step', lw=2., label=r'$-30<\delta<-5$')
    plt.hist(mc['trueE'][dec_msks2[2]], bins=trueE_bins, weights=A_eff_w2[2],
            histtype='step', lw=2., label=r'$-5<\delta<30$')
    plt.hist(mc['trueE'][dec_msks2[3]], bins=trueE_bins, weights=A_eff_w2[3],
            histtype='step', lw=2., label=r'$30<\delta<90$')

    plt.semilogx()
    plt.semilogy()
    plt.legend(loc=2, fontsize='small')
    plt.xlabel(r'$E_\nu$ [GeV]')
    plt.ylabel(r'$A_{eff}$ [$m^2$]') #r means latex only input
    ax.set_title(r'Effective area', fontsize=18)

    plt.savefig('/home/jthwaites/FRB/plots/effective_areas_fine.png')

def effective_area_comparison():
    from matplotlib.gridspec import GridSpec
    #comparison of effective areas with different datasamples
    mc = np.load('/data/ana/analyses/gfu/version-002-p06/IC86_2011_MC.npy')
    mc_ps=np.load('/data/ana/analyses/ps_tracks/version-004-p00/IC86_2016_MC.npy')
    mc_tt=np.load('/data/ana/analyses/transient_tracks/version-001-p02/IC86_2011_North_nugen.npy')

    nbins=28
    trueE_bins=np.logspace(2,9,num=nbins+1) #ow units: 1/(GeV cm^2 sr)
    bin_corr=(9.-2.)/nbins
    decmsk1=mc['trueDec']>(-5*np.pi/180) 
    decmsk2=mc['trueDec']<=(-5*np.pi/180)
    dec_msks=[decmsk1,decmsk2] #GFU
    
    decmsk3=mc_ps['trueDec']>(-5*np.pi/180) 
    decmsk4=mc_ps['trueDec']<=(-5*np.pi/180)
    dec_msksps=[decmsk3,decmsk4]#PS Tracks
    
    decmsk5=mc_tt['trueDec']>(-5*np.pi/180) 
    decmsk6=mc_tt['trueDec']<=(-5*np.pi/180)
    dec_mskstt=[decmsk5,decmsk6] #Transient tracks
    
    cuts=[np.pi/2,-5*np.pi/180,-np.pi/2]
    A_eff_w=[[],[]]
    A_eff_w2=[[],[]]
    A_eff_w3=[[],[]]
    labels=[r'GFU, dec$>-5$',r'GFU, dec$<-5$', 
            r'PS Tracks, dec$>-5$', r'PS Tracks, dec$<-5$', 
            r'Transient Tracks, $\delta>-5$', r'Transient Tracks, $\delta<-5$']
    
    for a in range(len(dec_msks)):
        solid_ang=2*np.pi*(np.sin(cuts[a])-np.sin(cuts[a+1]))
        weight=mc['ow'][dec_msks[a]]/(mc['trueE'][dec_msks[a]]*solid_ang*bin_corr*1e4*np.log(10))
        for w in weight: A_eff_w[a].append(w)
            
    for a in range(len(dec_msksps)):
        solid_ang=2*np.pi*(np.sin(cuts[a])-np.sin(cuts[a+1]))
        weight=mc_ps['ow'][dec_msksps[a]]/(mc_ps['trueE'][dec_msksps[a]]*solid_ang*bin_corr*1e4*np.log(10))
        for w in weight: A_eff_w2[a].append(w)
    
    for a in range(len(dec_mskstt)):
        solid_ang=2*np.pi*(np.sin(cuts[a])-np.sin(cuts[a+1]))
        weight=mc_tt['ow'][dec_mskstt[a]]/(mc_tt['trueE'][dec_mskstt[a]]*solid_ang*bin_corr*1e4*np.log(10))
        for w in weight: A_eff_w3[a].append(w)
    
    fig = plt.figure(dpi=200)
    gs = GridSpec(4,4)
    ax_effa = fig.add_subplot(gs[0:3,0:3])
    ax_ratio = fig.add_subplot(gs[3:4,0:3]) 
        
    ratio_xvals=[(trueE_bins[i+1]+trueE_bins[i])/2 for i in range(len(trueE_bins)-1)]
    for i in range(len(cuts)-1):
        h_ps=ax_effa.hist(mc_ps['trueE'][dec_msksps[i]], bins=trueE_bins, weights=A_eff_w2[i],
                          histtype='step', lw=1.,color='C%i'%(i+2), label=labels[i+2])
        #h_tt=ax_effa.hist(mc_tt['trueE'][dec_mskstt[i]], bins=trueE_bins, weights=A_eff_w3[i],
        #                  histtype='step', lw=1.,color='C%i'%(i+4), label=labels[i+4])
        h_gfu=ax_effa.hist(mc['trueE'][dec_msks[i]], bins=trueE_bins, weights=A_eff_w[i],
                           histtype='step', lw=1., color='C%i'%(i),label=labels[i])
        ax_ratio.plot(ratio_xvals, [h_ps[0][i]/h_gfu[0][i] for i in range(len(h_ps[0]))],
                         '.', color='C%i'%(i+2))
        #ax_ratio.plot(ratio_xvals, [h_tt[0][i]/h_gfu[0][i] for i in range(len(h_tt[0]))],
        #                 '.', color='C%i'%(i+4))
    ax_ratio.plot(ratio_xvals,[1]*len(ratio_xvals),'--',lw=0.5,color='black')
    
    plt.setp(ax_effa.get_xticklabels(), visible=False)
    ax_effa.set_yscale('log')
    ax_effa.set_xscale('log')
    ax_ratio.set_xscale('log')
    ax_effa.set_xlim([50,2e9])
    ax_ratio.set_xlim([50,2e9])

    ax_effa.set_ylabel(r'$A_{eff}$ [$m^2$]')
    ax_ratio.set_ylabel('Ratio:\n Dataset/GFU')
    ax_effa.legend(loc=4, fontsize=7)
    #ax_ratio.legend(loc=0, fontsize=5)
    ax_ratio.set_xlabel(r'$E_\nu$ [GeV]')
    ax_effa.set_title(r'Effective area')
    
    plt.savefig('/home/jthwaites/FRB/plots/effarea_comp.png')

def bg_pdfs():
    exp15 = np.load('/data/ana/analyses/gfu/version-002-p06/IC86_2015_data.npy')
    exp16= np.load('/data/ana/analyses/gfu/version-002-p06/IC86_2016_data.npy')
    exp17= np.load('/data/ana/analyses/gfu/version-002-p06/IC86_2017_data.npy')
    exp18= np.load('/data/ana/analyses/gfu/version-002-p06/IC86_2018_data.npy')
    exp19= np.load('/data/ana/analyses/gfu/version-002-p06/IC86_2019_data.npy')

    decs15=exp15['dec']
    decs16=exp16['dec']
    decs17=exp17['dec']
    decs18=exp18['dec']
    decs19=exp19['dec']

    fig, ax = plt.subplots(figsize=(12,8))

    plt.hist(np.sin(decs15), histtype='step', lw=2., label='2015')
    plt.hist(np.sin(decs16), histtype='step', lw=2., label='2016')
    plt.hist(np.sin(decs17), histtype='step', lw=2., label='2017')
    plt.hist(np.sin(decs18), histtype='step', lw=2., label='2018')
    plt.hist(np.sin(decs19), histtype='step', lw=2., label='2019')

    plt.legend(loc=2)
    plt.xlabel(r'sin(dec)')
    plt.ylabel(r'$N$') #r means latex only input
    ax.set_title(r'Background PDFs', fontsize=18)
    
    plt.savefig('/home/jthwaites/FRB/plots/bg_spatial_pdf.png')

def gfu_rate_plot():
    grl = np.load('/data/ana/analyses/gfu/version-002-p06/GRL/IC86_2015_data.npy')
    
    fig, ax = plt.subplots(figsize=(12,8))
    plt.errorbar(grl['start'], grl['events']/grl['livetime']/86400., 
                yerr=0.002, linestyle='', marker='*')
    plt.xlabel(r'MJD')
    plt.ylabel(r'rate [Hz]') #r means latex only input
    ax.set_title(r'Background rate', fontsize=18)

    plt.savefig('/home/jthwaites/FRB/plots/gfu_rate.png')


####Run functions to make plots
effective_area_gfu()
effective_area_comparison()
bg_pdfs()
gfu_rate_plot()