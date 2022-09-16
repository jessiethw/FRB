"""Plotting functions to make various plots
Jessie Thwaites 
June 2022
"""

import numpy as np
import matplotlib as mpl
mpl.use('agg')

import csky as cy
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy import units as u
import healpy as hp
import histlite as hl
import meander
import matplotlib.pyplot as plt
import os
import pickle as pkl

mpl.style.use('/home/jthwaites/FRB/scripts/plots_style.mplstyle')

#######################################################################
#plot sensitivites over time windows for FRBs
def sens_tw_plot(frb_name, time_windows, path='/home/jthwaites/FRB/sens_trials/point_source'):
    #time_windows is an array of all timewindows to load trials from
    for tw in time_windows:
        tw_str=str(round(tw,1))
        pkl_path=path+f'/{frb_name}_sens_{tw_str}.pkl'

        if not os.path.exists(pkl_path):
            #skip time window if not found
            print(f'pickle not found for tw={tw_str}')
            continue

        with open(pkl_path, 'rb') as f:
            sens = pkl.load(f)
        print(sens)

#######################################################################
#make ts dist by loading pickle file
def make_ts_dist_from_pickle(path):
    import os
    filename=os.path.basename(path)
    source=filename[:12]
    tw=float(filename.split('_')[-1][:-4])
    if tw<2.:
        tw=round(tw,2)
    else:
        tw=int(tw)
    
    trials=np.load(path, allow_pickle=True)

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
    plt.title(r'BG TS distribution, %s, with prior (%is)'%(source,tw))
    ax.legend()

    plt.savefig(f'/home/jthwaites/FRB/plots/bg_ts_distributions/{source}_bgts_{tw}s_w_prior.png')

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
    
    plt.savefig('/home/jthwaites/plots/frb_skymap.png')

########################################################################
#make visualizations for a smaller skymap area than full skymap
##all fxns below needed for this
def plot_zoom(scan, ra, dec, levels=[0.68,0.90], title='', reso=3, 
              var="pVal", ts_range=None,cmap='Blues', contour_scan=None, CL=None, 
              col_label=r"Probability"):
    if ts_range ==None: 
        if max(scan)<= 0.: 
            ts_range=[np.around(min(scan),decimals=3),0.]
            cmap=cmap+'_r'
        else: 
            ts_range=[0.,np.around(max(scan),decimals=3)]
            cmap=cmap
    elif max(scan)<=0.: cmap=cmap+'_r'
   
#    if cmap is None:
#        pdf_palette = sns.color_palette(cmap, 500)
#        cmap = mpl.colors.ListedColormap(pdf_palette)
    hp.gnomview(scan, rot=(np.degrees(ra), np.degrees(dec), 0),
                    cmap="Blues",
                    max=max(scan),
                    min=min(scan),
                    reso=reso,
                    title=title,
                    notext=True,
                    cbar=False,
                    nest=False
                    #unit=r""
                    )
    levels=levels # 68%,90% containment 
    
    if CL is not None:
        #FIX HERE
        theta, phi = plot_contours_with_mask(CL, levels)
    elif contour_scan is None:
        theta, phi = plot_contours(levels,scan)
    elif contour_scan is not False:
        theta, phi = plot_contours(levels,contour_scan)
    else: 
        theta=[0.]
        phi=[0.]
    
    if len(theta) and len(phi) != 0:
        hp.projplot(theta[0],phi[0],linewidth=1,c='k',label='%i%% C.L'%(levels[0]*100))
        for i in np.arange(1,len(theta)): 
            hp.projplot(theta[i],phi[i],c='k',linewidth=1)
        
    plt.plot(4.95/3.*reso*np.radians([-1, 1, 1, -1, -1]), 
             4.95/3.*reso*np.radians([1, 1, -1, -1, 1]), color="k", ls="-", lw=3)

    draw_axes(dec, ra, reso)
        
    plot_color_bar(cmap=cmap, labels=ts_range,col_label=col_label)
    #hp.graticule(verbose=False)

def draw_axes(src_dec, src_ra, reso, axis_labels=True):
    plt.plot(
        4.95/3.*reso*np.radians([-1, 1, 1, -1, -1]), 
        4.95/3.*reso*np.radians([1, 1, -1, -1, 1]), 
        color="k", ls="-", lw=3
        )
    ra_scale_factor = 3 if np.degrees(np.abs(src_dec)) < 30. else 10
    num_ra_lines = ra_scale_factor*reso
    num_dec_lines = 3*reso
    ra_axes = np.linspace(
        np.degrees(src_ra)-360.,
        np.degrees(src_ra)+360.,
        721
        )
    ra_axes = ra_axes[
        (ra_axes > (np.degrees(src_ra) - num_ra_lines)) & \
        (ra_axes < (np.degrees(src_ra) + num_ra_lines))
        ]
    ra_axes = np.where(ra_axes > 360., ra_axes - 360., ra_axes)
    dec_axes = np.linspace(
        np.degrees(src_dec) - 180.,
        np.degrees(src_dec) + 180.,
        361.
        )
    dec_axes = dec_axes[
        (dec_axes > (np.degrees(src_dec) - num_dec_lines)) & \
        (dec_axes < (np.degrees(src_dec) + num_dec_lines))
        ]
    dec_axes = dec_axes[(dec_axes > -90.) & (dec_axes < 90.)]
    for tmp_ra in ra_axes:
        tmp_line_ra = np.radians(np.ones(500)*tmp_ra)
        tmp_line_dec = np.radians(np.linspace(dec_axes[0], dec_axes[-1], 500))
        hp.projplot(np.pi/2. - tmp_line_dec, tmp_line_ra, linewidth=0.5, 
                    color='k', linestyle='dotted', coord='C')
    for tmp_dec in dec_axes:
        tmp_line_dec = np.radians(np.ones(500)*tmp_dec)
        tmp_line_ra = np.radians(np.linspace(ra_axes[0], ra_axes[-1], 500))
        hp.projplot(np.pi/2. - tmp_line_dec, tmp_line_ra, linewidth=0.5, 
                    color='k', linestyle='dotted', coord='C')
    plot_labels(src_dec, src_ra, reso, with_axis_labels=axis_labels)    

def plot_color_bar(labels=[0.,2.,4.,6.], col_label=r"Probability", range=[0,6], cmap=None, offset=-35):
    fig = plt.gcf()
    #ax = fig.add_axes([0.25, -0.03, 0.5, 0.03])
    ax = fig.add_axes([0.95, 0.2, 0.03, 0.6])
    cb = mpl.colorbar.ColorbarBase(ax, cmap='Blues' if cmap is None else cmap,
                        #norm=mpl.colors.Normalize(vmin=range[0], vmax=range[1]), 
                        orientation="vertical")
    #cb.ax.minorticks_on()
    cb.set_label(col_label,labelpad=offset, fontsize=18)
    cb.set_ticks([0., 1.])
    cb.set_ticklabels(labels)
    cb.update_ticks()
    #cb.ax.get_xaxis().set_ticklabels(labels)

def plot_labels(src_dec, src_ra, reso, with_axis_labels=True):
    """Add labels to healpy zoom"""
    fontsize = 20
    
    if np.degrees(src_dec) > 35.:
        ras = [
            np.degrees(src_ra) - reso*3,
            np.degrees(src_ra),
            np.degrees(src_ra) + reso*3
            ]
    else:
        ras = [
            np.degrees(src_ra) - reso,
            np.degrees(src_ra),
            np.degrees(src_ra) + reso
            ]
    decs = [
        np.degrees(src_dec) - reso,
        np.degrees(src_dec),
        np.degrees(src_dec) + reso
        ]
    for ra in ras:
        dec_text = np.pi/2. - src_dec + (np.radians(reso) * 5.5/3.)
        dec_offset = np.abs(np.radians(ra - np.degrees(src_ra)))*np.sin(src_dec)*reso*0.01
        ra_text = r"%.f$^\circ$"%ra if ra < 360. else r"%.f$^\circ$"%ra - 360.
        hp.projtext(dec_text + dec_offset, np.radians(ra),
                    ra_text, lonlat=False,
                    fontsize=20, ha='center')
    for dec in decs:
        scale = np.degrees(hp.rotator.angdist(
            [np.pi/2. - np.radians(dec), src_ra], 
            [np.pi/2. - np.radians(dec), src_ra + np.radians(1.)]
            ))
        ra_text = np.radians(np.degrees(src_ra) + 1.7*reso/scale)
        hp.projtext(np.pi/2. - np.radians(dec), ra_text,
            r"%.2f$^\circ$"%dec, lonlat=False,
            fontsize=20, ha='right')
    
    if with_axis_labels:
        plt.text(-1.05*np.radians(2.4*reso), np.radians(0), r"declination",
                    ha='center', va='center', rotation=90, fontsize=fontsize)
        plt.text(np.radians(0), np.radians(-2*reso), r"right ascension",
                    ha='center', va='center', fontsize=fontsize)
        
def plot_contours(proportions,samples):
    r''' Plot containment contour around desired level.
    E.g 90% containment of a PDF on a healpix map    
    Parameters:
    -----------
    proportions: list
        list of containment level to make contours for.
        E.g [0.68,0.9]
    samples: array
        array of values read in from healpix map
        E.g samples = hp.read_map(file)
    Returns:
    --------
    theta_list: list
        List of arrays containing theta values for desired contours
    phi_list: list
        List of arrays containing phi values for desired contours
    '''    
    levels = []
    sorted_samples = list(reversed(list(sorted(samples))))
    
    nside = hp.pixelfunc.get_nside(samples)
    sample_points = np.array(hp.pix2ang(nside,np.arange(len(samples)))).T
    
    for proportion in proportions:
        level_index = (np.cumsum(sorted_samples) > proportion).tolist().index(True)
        level = (sorted_samples[level_index] + (sorted_samples[level_index+1] 
                                        if level_index+1 < len(samples) else 0)) / 2.0
        levels.append(level)
    contours_by_level = meander.spherical_contours(sample_points, samples, levels)    
    
    theta_list = []; phi_list=[]
    for contours in contours_by_level:
        for contour in contours:
            theta, phi = contour.T
            phi[phi<0] += 2.0*np.pi
            theta_list.append(theta)
            phi_list.append(phi)    
    return theta_list, phi_list     

