""" Functions for computations and visualizations of Chime spatial priors
Jessie Thwaites, 8/20/21
"""

import numpy as np
import h5py as h5
import healpy as hp
import os

#########################################################################
#Globals
chime_loc='/data/user/jthwaites/chime_localization_data/'
chime_loc_files=os.listdir(chime_loc)

#########################################################################
#full extent of 68% contour of map, not just central lobe
def get_new_localizations():
    new_localizations={'tns_name':[],'ra_unc':[],'dec':[],'dec_unc':[]}
    
    for file in chime_loc_files:
        f = h5.File(chime_loc+file,'r')
    
        ra_loc=f['/'].attrs['ra_error']
        dec_loc=f['/'].attrs['dec_error']

        extent=[[0,361],[0,91]]
        for name, contour in f['/contours/68/'].items():
            contour=contour[:]
            max_ra=max(contour[0])
            if max_ra> extent[0][0]: extent[0][0]=max_ra
            min_ra=min(contour[0])
            if min_ra< extent[0][1]: extent[0][1]=min_ra
            max_dec=max(contour[1])
            if max_dec> extent[1][0]: extent[1][0]=max_dec
            min_dec=min(contour[1])
            if min_dec < extent[1][1]: extent[1][1]=min_dec
    
        if (extent[0][0]-extent[0][1]) > ra_loc: ra_loc=extent[0][0]-extent[0][1]
        if (extent[1][0]-extent[1][1]) > dec_loc: dec_loc=extent[1][0]-extent[1][1]
    
        new_localizations['tns_name'].append(f['/'].attrs['tns_name'])
        new_localizations['dec'].append(f['/'].attrs['dec'])
        new_localizations['ra_unc'].append(ra_loc)
        new_localizations['dec_unc'].append(dec_loc)
    return new_localizations

############################################################################
#make a healpix map for one frb
def make_healpix_map(tns_name, new_nside=512, max_cl=0.90):
    new_localizations=get_new_localizations()
    
    index=np.where(tns_name==np.asarray(new_localizations['tns_name']))[0][0]
    frb = chime_loc + tns_name+'_localization.h5'
    f = h5.File(frb, 'r')

    ra = np.radians(f['/'].attrs['ra'])
    dec = np.radians(f['/'].attrs['dec'])

    original_nside = f['healpix'].attrs['nside']
    ipix, CL = f['healpix/ipix'][()], f['healpix/CL'][()]
    #for some reason, it appears there are duplicates of these in all the files. 
    #To remove, use np.unqiue fxn
    unique_CLs, CL_ind=np.unique(ipix, return_index=True)

    probs = (1. - CL)/ (1. - CL[CL_ind]).sum()
    frb_probs = np.zeros(hp.nside2npix(original_nside))
    frb_probs[ipix] = probs
    new_CL = np.ones_like(frb_probs)
    new_CL[ipix] = CL
    
    frb_probs = hp.ud_grade(frb_probs,
                new_nside, order_in='NESTED', order_out='RING', power=-2)
    new_CL = hp.ud_grade(new_CL,
                new_nside, order_in='NESTED', order_out='RING')

    #get mask to use
    msk = new_CL < max_cl
    
    #print('sum of map probabilities: %.5f'%(frb_probs.sum()))
    return frb_probs, msk


        