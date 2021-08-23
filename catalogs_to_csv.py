"""
File which creates .csv files (to catalog folder) of all FRBs.
Outputs files to /home/jthwaites/FRB/catalog/
frbs_all.csv - all FRBs and bursts
frb121102_bursts.csv - all bursts of FRB121102
spatial_priors_frbs.csv - Chime bursts with spatial priors
frbs_excl_sppriors_121102.csv - remaining FRBs
"""


import os
import numpy as np
import pandas as pd
import json
from astropy.time import Time
from astropy.coordinates import SkyCoord
from astropy import units as u


def get_catalogs(rem_subbursts=True, return_dataframes=False):
###############################################################################
#Repeaters data from Chime catalog
    with open('./catalog/repeaters.txt') as json_file:
        repeaters = json.load(json_file)

    chime_repeaters={'source':[],'n_bursts':[], 'flux':[]}
    ra_dec=[[],[]]
    for source in repeaters:
        ra_dec[0].append(repeaters[source]['ra']['value'])
        ra_dec[1].append(repeaters[source]['dec']['value'])
        count=0
        for key in repeaters[source]:
            if key != 'dm' and key !='ymw16' and key!='dec' and key!='gl' and key!='gb' \
            and key!= 'ra' and key!='localized' and key!='last_burst_date' \
            and key!='publication' and key!='ne2001' and key!='previous_name': 
                count+=1
                if repeaters[source][key]['flux']['value'] != {}: 
                    chime_repeaters['flux'].append(repeaters[source][key]['flux']['value'])
                    
        if source=='190907.J08+46': chime_repeaters['source'].append('FRB'+source)
        else:chime_repeaters['source'].append(source)
        chime_repeaters['n_bursts'].append(count)
        
    #source and number data, coordinates for each burst
    chime_repeaters_coord=SkyCoord(ra=ra_dec[0], dec=ra_dec[1], unit=(u.hourangle, u.deg))

    #mjd array for repeater bursts
    rep_utc=[]
    for source in repeaters:
        for event in repeaters[source].keys():
            if 'timestamp' in repeaters[source][event].keys():
                rep_utc.append(repeaters[source][event]['timestamp']['value'])
            
    chime_repeaters_srcs={'src':[],'ra_deg':[],'dec_deg':[],'mjd':[]}
    for i in range(len(chime_repeaters['source'])):
        brst=[chime_repeaters['source'][i]]*chime_repeaters['n_bursts'][i]
        [chime_repeaters_srcs['src'].append(b) for b in brst]
        ra=[chime_repeaters_coord[i].ra.deg]*chime_repeaters['n_bursts'][i]
        [chime_repeaters_srcs['ra_deg'].append(c) for c in ra]
        dec=[chime_repeaters_coord[i].dec.deg]*chime_repeaters['n_bursts'][i]
        [chime_repeaters_srcs['dec_deg'].append(c) for c in dec]
    chime_repeaters_srcs['mjd'] = np.asarray([Time(t, format='iso').mjd for t in rep_utc])

    #all repeaters as a dataframe
    chime_rep=pd.DataFrame.from_dict(data=chime_repeaters_srcs)

###############################################################################
#FRBCat
    frbcat = pd.read_csv('./catalog/frbcat_20210519_all.csv')

    #CHIME repeaters: remove repeaters from FRBCat to not double-count events
    chimes=np.where(frbcat['telescope']=='CHIME/FRB')
    frbcheck=[] #names from frbcat
    for n in chimes[0]:
        frbcheck.append(frbcat['frb_name'][n][0:9])
    
    frbcat_format=''
    rem=[[],[]]
    remove_121102=np.where(frbcat['frb_name']=='FRB121102')[0]
    for source in chime_repeaters['source']:
        if source[0:3]=='FRB':
            frbcat_format=source[0:3]+source[5:-1]
        else:
            frbcat_format='FRB'+source[0:6]
    
        res = any(ele in frbcat_format for ele in frbcheck)
        if bool(res) == 1:
            rem[0].append(frbcheck[frbcheck.index(frbcat_format)])
    
    for ind in range(len(frbcat['frb_name'])):
        res = any(ele in frbcat['frb_name'][ind] for ele in rem[0])
        if bool(res)==1:
            rem[1].append(ind)
        
    frbcat=frbcat.drop(rem[1])
    loc_121102=SkyCoord(ra=frbcat['rop_raj'][remove_121102[0]], 
                    dec=frbcat['rop_decj'][remove_121102[0]], unit=(u.hourangle, u.deg))
    frbcat=frbcat.drop(remove_121102)

    frbs_mjd = np.asarray([Time(t.replace('/', '-'), format='iso').mjd for t in frbcat['utc']]) 
    frbs_coord = SkyCoord(ra=frbcat['rop_raj'], dec=frbcat['rop_decj'],unit=(u.hourangle, u.deg))

###############################################################################
# Chime catalog 1
    chime_cat1 = pd.read_csv('./catalog/chimefrbcat1.csv')

    if rem_subbursts==True: remove_subburst=[]
    remove_121102=[]
    chimecat1_names=[]
    for i in range(len(chime_cat1)):
        #remove subbursts -> has to be first to rem repeater subbursts
        #for now - change to lightcurve in future?
        if rem_subbursts==True and chime_cat1['sub_num'].values[i]!=0: 
            remove_subburst.append(chime_cat1['sub_num'].index[i])
        elif chime_cat1['previous_name'][i]== '190907.J08+46': #not present
            chimecat1_names.append(chime_cat1['previous_name'][i])
        elif chime_cat1['previous_name'][i]== 'FRB121102':
            remove_121102.append(i)
        elif chime_cat1['repeater_name'][i]!='-9999':
            chimecat1_names.append(chime_cat1['repeater_name'][i])
        else: chimecat1_names.append(chime_cat1['tns_name'][i])
        
    
    if rem_subbursts==True: chime_cat1=chime_cat1.drop(remove_subburst+remove_121102)
    else: chime_cat1=chime_cat1.drop(remove_121102)
        
    #remove repeater bursts that appear in both chime catalogs
    chime1_times=np.asarray([Time('2018-07-25 00:00:00', format='iso').mjd,
                         Time('2019-07-01 23:59:59', format='iso').mjd ]) 
    remove_rep_duplicates=[]
    for i in range(len(chime_repeaters_srcs['mjd'])):
        if chime_repeaters_srcs['mjd'][i] > chime1_times[0] and \
        chime_repeaters_srcs['mjd'][i] <chime1_times[1]:
            remove_rep_duplicates.append(i)

    chime_rep=chime_rep.drop(remove_rep_duplicates)
    
###############################################################################
#FRB121102 catalog: 234 bursts
    frb121102_brsts = pd.read_csv('./catalog/frb121102_bursts.csv')
    
###############################################################################
#combined DataFrame for all catalogs
#repeater is a bool array: True=rep, False=single burst
    frbs_all = {'src':[], 'ra_deg':[],'dec_deg':[],'mjd':[],'repeater':[]} 

    [frbs_all['src'].append(name_f) for name_f in frbcat['frb_name']]
    [frbs_all['src'].append(name_c) for name_c in chimecat1_names]
    [frbs_all['src'].append(name_r) for name_r in chime_rep['src']]
    [frbs_all['src'].append('FRB121102') for i in range(len(frb121102_brsts))]

    [frbs_all['ra_deg'].append(ra_f) for ra_f in frbs_coord.ra.deg]
    [frbs_all['ra_deg'].append(ra_c) for ra_c in chime_cat1['ra']]
    [frbs_all['ra_deg'].append(ra_r) for ra_r in chime_rep['ra_deg']]
    [frbs_all['ra_deg'].append(loc_121102.ra.deg) for i in range(len(frb121102_brsts))]

    [frbs_all['dec_deg'].append(dec_f) for dec_f in frbs_coord.dec.deg]
    [frbs_all['dec_deg'].append(dec_c) for dec_c in chime_cat1['dec']]
    [frbs_all['dec_deg'].append(dec_r) for dec_r in chime_rep['dec_deg']]
    [frbs_all['dec_deg'].append(loc_121102.dec.deg) for i in range(len(frb121102_brsts))]

    [frbs_all['mjd'].append(mjd_f) for mjd_f in frbs_mjd]
    [frbs_all['mjd'].append(mjd_c) for mjd_c in chime_cat1['mjd_400']]
    [frbs_all['mjd'].append(mjd_r) for mjd_r in chime_rep['mjd']]
    [frbs_all['mjd'].append(mjd_121102) for mjd_121102 in frb121102_brsts.mjd.values]

    #creating repeater bool array: True=repeater, False=non-repeater
    unique_frbs, ind, n_frbs = np.unique(frbs_all['src'], return_counts=True, return_index=True)
    msk=n_frbs!=1
    rep=np.sort([frbs_all['src'][m] for m in ind[np.where(msk)[0]]])
    
    for name in frbs_all['src']: 
        if name in rep: frbs_all['repeater'].append(True)
        else: frbs_all['repeater'].append(False)

    frbs_all_data=pd.DataFrame.from_dict(data=frbs_all)
    frbs_all_data.to_csv('./catalog/frbs_all.csv',index=False)  
    
###############################################################################
#make a cat of frb121102 bursts only, spatial priors, and rem those from main cat
    frb121102_bursts={'src':[], 'ra_deg':[],'dec_deg':[],'mjd':[]} 
    for mjd in frb121102_brsts.mjd.values:
        frb121102_bursts['mjd'].append(mjd) 
        frb121102_bursts['src'].append('FRB121102')
        frb121102_bursts['ra_deg'].append(loc_121102.ra.deg)
        frb121102_bursts['dec_deg'].append(loc_121102.dec.deg)
    frb121102_bursts=pd.DataFrame.from_dict(data=frb121102_bursts)
    frb121102_bursts.to_csv('./catalog/frb121102_bursts.csv',index=False)

    j=np.arange(len(frbs_all_data['src'].values)-len(frb121102_brsts),
            len(frbs_all_data['src'].values))
    frbs_excl=frbs_all_data.drop(j)
###############################################################################
#Spatial priors from Chime cat1
    dat_loc='/data/user/jthwaites/chime_localization_data/'
    dat_files=os.listdir(dat_loc)
    spatial_priors=[]
    [spatial_priors.append(filename[0:12]) for filename in dat_files]
    spatial_priors_frbs={'src':[], 'ra_deg':[],'dec_deg':[],'mjd':[], 'repeater':[]} 
    msk=[]
    for k in range(len(spatial_priors)):
        index=np.where(spatial_priors[k] == frbs_excl['src'].values)[0]
        for i in index:
            msk.append(i)
            spatial_priors_frbs['src'].append(spatial_priors[k])
            spatial_priors_frbs['ra_deg'].append(frbs_excl['ra_deg'].values[i])
            spatial_priors_frbs['dec_deg'].append(frbs_excl['dec_deg'].values[i])
            spatial_priors_frbs['mjd'].append(frbs_excl['mjd'].values[i])
            spatial_priors_frbs['repeater'].append(frbs_excl['repeater'].values[i])

    spatial_priors_frbs=pd.DataFrame.from_dict(data=spatial_priors_frbs)
    spatial_priors_frbs.to_csv('./catalog/spatial_priors_frbs.csv',index=False)

###############################################################################
#remaining FRBs not in either of FRB121102 or spatial priors catalog
    frbs_excl=frbs_excl.drop(msk)
    frbs_excl.to_csv('./catalog/frbs_excl_sppriors_121102.csv',index=False)
    
    if return_dataframes==True: 
        return frbs_all_data, frb121102_bursts, spatial_priors_frbs, frbs_excl
    else: return
    

#####################################################################        
#load catalog csv file - default is load all if not specified
def load_frbs(**kw):
    if 'spatial_priors' in kw.keys(): 
        frbs = pd.read_csv('./catalog/spatial_priors_frbs.csv')
    elif 'frb121102' in kw.keys():
        frbs = pd.read_csv('./catalog/frb121102_bursts.csv')
    elif 'others' in kw.keys(): 
        frbs = pd.read_csv('./catalog/frbs_excl_sppriors_121102.csv')
    elif 'all' in kw.keys(): frbs = pd.read_csv('./catalog/frbs_all.csv')
    else: 
        frbs = pd.read_csv('./catalog/frbs_all.csv')
        print('Loading all FRBs from catalogs')    
    return frbs
    