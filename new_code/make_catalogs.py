"""
Collects and writes all FRBs to csv.
Outputs files to /home/jthwaites/FRB/catalog/
Files:
 - frbs_all.csv - all FRBs and bursts
 - frb121102_bursts.csv - all bursts of FRB121102
 - spatial_priors_frbs.csv - Chime bursts with spatial priors

Jessie Thwaites, updated June 2022
"""

import numpy as np
import json
import pandas as pd
from astropy.time import Time
from astropy.coordinates import SkyCoord
from astropy import units as u


remove_subbursts=True #flag to remove subbursts from catalogs

##############################################
#Chime repeaters catalog

with open('/home/jthwaites/FRB/catalog/repeaters.txt') as json_file:
    repeaters = json.load(json_file)

repeater_names=[source for source in repeaters]
repeater_names[repeater_names=='190907.J08+46']='FRB190907.J08+46'

repeater_locations=SkyCoord(
    ra=[repeaters[source]['ra']['value'] for source in repeaters],
    dec=[repeaters[source]['dec']['value'] for source in repeaters],
    unit=(u.hourangle, u.deg)
    )

chime_repeaters={
    'frb_name':repeater_names,
    'ra_deg':repeater_locations.ra.deg,
    'dec_deg':repeater_locations.dec.deg,
    'n_bursts':[], 
    'flux':[],
    'mjd':[]
    }

#remove repeater bursts that appear in both chime catalogs
chime1_times=np.asarray([Time('2018-07-25 00:00:00', format='iso').mjd,
                         Time('2019-07-01 23:59:59', format='iso').mjd ]) 

for source in repeaters:
    count=0  #count number of bursts for each repeater
    bursts_mjd=[] #collect mjd of each burst
    fluxes=[] #collect flux value for each burst (if given)
    for key in repeaters[source]:
        if key != 'dm' and key !='ymw16' and key!='dec' and key!='gl' and key!='gb' \
        and key!= 'ra' and key!='localized' and key!='last_burst_date' \
        and key!='publication' and key!='ne2001' and key!='previous_name': 
            count+=1
            if repeaters[source][key]['flux']['value'] != {}: 
                fluxes.append(repeaters[source][key]['flux']['value'])
            else:
                fluxes.append(None)
            burst_mjd=Time(repeaters[source][key]['timestamp']['value'], format='iso').mjd
            if burst_mjd < chime1_times[0] or \
                burst_mjd > chime1_times[1]:
                bursts_mjd.append(burst_mjd)
    
    chime_repeaters['n_bursts'].append(count)
    chime_repeaters['mjd'].append(bursts_mjd)
    chime_repeaters['flux'].append(fluxes)

chime_repeaters=pd.DataFrame(data=chime_repeaters)
#if no bursts are left after dropping those in Chime cat 1, remove from list
n=np.where(np.asarray([len(mjd_list) for mjd_list in chime_repeaters['mjd'].values])==0)[0]
chime_repeaters=chime_repeaters.drop(n)

##############################################
#FRBCat
frbcat_original = pd.read_csv('/home/jthwaites/FRB/catalog/frbcat_20210519_all.csv')

#drop bursts from CHIME/FRB - these are in other catalogs
frbcat=frbcat_original.drop(frbcat_original[frbcat_original['telescope']=='CHIME/FRB'].index)

#drop bursts from FRB121102 - these are in other catalogs
frbcat=frbcat.drop(frbcat[frbcat['frb_name'].values=='FRB121102'].index)
frbcat_mjd = np.asarray([Time(t.replace('/', '-'), format='iso').mjd for t in frbcat['utc']]) 
frbcat_locations = SkyCoord(ra=frbcat['rop_raj'], dec=frbcat['rop_decj'],unit=(u.hourangle, u.deg))

##############################################
#FRB121102 catalog: 234 bursts
frb121102 = pd.read_csv('/home/jthwaites/FRB/catalog/frb121102_bursts.csv')
i=frbcat_original[frbcat_original['frb_name'].values=='FRB121102'].index[0]
location_121102=SkyCoord(ra=frbcat_original['rop_raj'][i], 
                    dec=frbcat_original['rop_decj'][i], unit=(u.hourangle, u.deg))

##############################################
# CHIME Catalog 1
chime_cat1 = pd.read_csv('/home/jthwaites/FRB/catalog/chimefrbcat1.csv')

if remove_subbursts==True: remove_subburst=[]
remove_121102=[]
chimecat1_names=[]
for i in range(len(chime_cat1)):
    if remove_subbursts==True and chime_cat1['sub_num'].values[i]!=0: 
        remove_subburst.append(chime_cat1['sub_num'].index[i])
    elif chime_cat1['previous_name'][i]== '190907.J08+46': #not present
        chimecat1_names.append(chime_cat1['previous_name'][i])
    elif chime_cat1['previous_name'][i]== 'FRB121102':
        remove_121102.append(i)
    elif chime_cat1['repeater_name'][i]!='-9999':
        chimecat1_names.append(chime_cat1['repeater_name'][i])
    else: chimecat1_names.append(chime_cat1['tns_name'][i])

if remove_subbursts==True: chime_cat1=chime_cat1.drop(remove_subburst+remove_121102)
else: chime_cat1=chime_cat1.drop(remove_121102)

##############################################
# saving all of these in DataFrame objects
# repeater is a bool array: True=rep, False=single burst

frbs_all = {'src':[], 'ra_deg':[],'dec_deg':[],'mjd':[],'catalog':[],'repeater':[]} 

[frbs_all['src'].append(name_f) for name_f in frbcat['frb_name']]
[frbs_all['src'].append(name_c) for name_c in chimecat1_names]
[frbs_all['src'].append(name_r) for name_r in chime_repeaters['src']]
[frbs_all['src'].append('FRB121102') for i in range(len(frb121102))]

[frbs_all['ra_deg'].append(ra_f) for ra_f in frbcat_locations.ra.deg]
[frbs_all['ra_deg'].append(ra_c) for ra_c in chime_cat1['ra']]
[frbs_all['ra_deg'].append(ra_r) for ra_r in chime_repeaters['ra_deg']]
[frbs_all['ra_deg'].append(location_121102.ra.deg) for i in range(len(frb121102))]

[frbs_all['dec_deg'].append(dec_f) for dec_f in frbcat_locations.dec.deg]
[frbs_all['dec_deg'].append(dec_c) for dec_c in chime_cat1['dec']]
[frbs_all['dec_deg'].append(dec_r) for dec_r in chime_repeaters['dec_deg']]
[frbs_all['dec_deg'].append(location_121102.dec.deg) for i in range(len(frb121102))]

[frbs_all['mjd'].append(mjd_f) for mjd_f in frbcat_mjd]
[frbs_all['mjd'].append(mjd_c) for mjd_c in chime_cat1['mjd_400']]
[frbs_all['mjd'].append(mjd_r) for mjd_r in chime_repeaters['mjd']]
[frbs_all['mjd'].append(mjd_121102) for mjd_121102 in frb121102.mjd.values]
    
[frbs_all['catalog'].append('FRBCat') for i in range(len(frbcat['frb_name']))]
[frbs_all['catalog'].append('CHIME_1') for i in range(len(chimecat1_names))]
[frbs_all['catalog'].append('CHIME_rep') for i in range(len(chime_repeaters['src']))]
[frbs_all['catalog'].append('FRB121102') for i in range(len(frb121102))]

#creating repeater bool array: True=repeater, False=non-repeater
unique_frbs, ind, n_frbs = np.unique(frbs_all['src'], return_counts=True, return_index=True)
msk=n_frbs!=1
rep=np.sort([frbs_all['src'][m] for m in ind[np.where(msk)[0]]])
    
for name in frbs_all['src']: 
    if name in rep: frbs_all['repeater'].append(True)
    else: frbs_all['repeater'].append(False)

frbs_all_data=pd.DataFrame.from_dict(data=frbs_all)
frbs_all_data.to_csv('/home/jthwaites/FRB/catalog/frbs_all.csv',index=False) 