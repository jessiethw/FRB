import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.style.use('/home/jthwaites/FRB/scripts/plots_style.mplstyle')

def effective_area():
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
effective_area()
bg_pdfs()
gfu_rate_plot()