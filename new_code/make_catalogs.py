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
import os

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

frbs_all = {
    'src':np.concatenate((np.asarray(frbcat['frb_name'].values),np.asarray(chimecat1_names),
                        chime_repeaters['frb_name'].values,np.asarray(['FRB121102']*len(frb121102)))),
    'ra_deg':np.concatenate((frbcat_locations.ra.deg,chime_cat1['ra'],
                        chime_repeaters['ra_deg'],[location_121102.ra.deg]*len(frb121102))),
    'dec_deg':np.concatenate((frbcat_locations.dec.deg,chime_cat1['dec'],
                        chime_repeaters['dec_deg'],[location_121102.dec.deg]*len(frb121102))),
    'mjd':np.concatenate((frbcat_mjd,chime_cat1['mjd_400'], chime_repeaters['mjd'],frb121102.mjd.values)),
    'catalog':np.concatenate((['FRBCat']*len(frbcat['frb_name']),['CHIME_Catalog1']*len(chimecat1_names),
                        ['CHIME_repeaters']*len(chime_repeaters['frb_name']), ['FRB121102']*len(frb121102))),
    'repeater':[]} 

#creating repeater bool array: True=repeater, False=non-repeater
unique_frbs, ind, n_frbs = np.unique(frbs_all['src'], return_counts=True, return_index=True)
repeaters=frbs_all['src'][ind[n_frbs!=1]]
for name in frbs_all['src']: 
    if name in repeaters: frbs_all['repeater'].append(True)
    else: frbs_all['repeater'].append(False)

all_frbs=pd.DataFrame.from_dict(data=frbs_all)
all_frbs.to_csv('/home/jthwaites/FRB/catalog/frbs_all.csv',index=False) 

##############################################
# catalog of all those with spatial priors
localization_files=os.listdir('/data/user/jthwaites/chime_localization_data/')
spatial_priors=[]
[spatial_priors.append(filename[0:12]) for filename in localization_files]

spatial_priors_frbs={'src':[], 'ra_deg':[],'dec_deg':[],'mjd':[], 'repeater':[]} 
for k in range(len(spatial_priors)):
    index=np.where(spatial_priors[k] == all_frbs['src'].values)[0]
    for i in index:
        spatial_priors_frbs['src'].append(spatial_priors[k])
        spatial_priors_frbs['ra_deg'].append(all_frbs['ra_deg'].values[i])
        spatial_priors_frbs['dec_deg'].append(all_frbs['dec_deg'].values[i])
        spatial_priors_frbs['mjd'].append(all_frbs['mjd'].values[i])
        spatial_priors_frbs['repeater'].append(all_frbs['repeater'].values[i])

spatial_priors_frbs=pd.DataFrame.from_dict(data=spatial_priors_frbs)
spatial_priors_frbs.to_csv('/home/jthwaites/FRB/catalog/spatial_priors_frbs.csv',index=False)

##############################################
# catalog of frb121102 bursts
frb121102_bursts=pd.DataFrame.from_dict(data={
    'src':['FRB121102']*len(frb121102), 
    'ra_deg':[location_121102.ra.deg]*len(frb121102),
    'dec_deg':[location_121102.dec.deg]*len(frb121102),
    'mjd':frb121102.mjd.values})

frb121102_bursts.to_csv('/home/jthwaites/FRB/catalog/frb121102_bursts.csv',index=False)