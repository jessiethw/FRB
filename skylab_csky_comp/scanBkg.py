#!/usr/bin/env python

r"""
Run trials for background only, all-sky scans. Record TS and 
number of true and fitted events around best fit location

"""


import numpy  as np
import healpy as hp
import argparse,time

from skylab.priors        import SpatialPrior
from config_GW            import config
from scipy.optimize       import curve_fit
from skylab.ps_injector   import PointSourceInjector
from scipy.stats          import chi2

######################### CONFIGURE ARGUEMENTS #############################
p = argparse.ArgumentParser(description="Calculates Sensitivity and Discovery"
                            " Potential Fluxes for Background Gravitational wave/Neutrino Coincidence study",
                            formatter_class=argparse.RawTextHelpFormatter)
p.add_argument("--ntrials", default=1000, type=int,
                help="Number of trials (default=1000")
p.add_argument("--pid", default=0, type=int,
                help="Process ID to save unique numpy array after running (Default=0)")
p.add_argument("--gw", default=151226, type=int,
                help="Date of GW event (yr-m-d) (default=151226)")
p.add_argument('--tw', default=3., type=float,
                help='Log10(time_window) to perform search in (default=3 which corresponds to 10^3 seconds)')
args = p.parse_args()
###########################################################################

###################### CONFIGURE LLH  ########################
seasons = ['GFU_v002p05','IC86, 2011-2018']
erange  = [0,10]
index = 2.
GW_time_dict = dict({150914:57279.41024306,151226:57382.152001999784,
                     151012:57307.4130028,170104:57757.42498380,
                     170608:57912.08421863,170729:57963.78922801,
                     170809:57974.35303009,170814:57979.43800382,
                     170817:57982.52852350,170818:57983.10079977,
                     170823:57988.55137153,})
GW_time = GW_time_dict[args.gw]
time_window = 10**args.tw/2/86400 #500 seconds in days
time_mask = [time_window,GW_time]

# For 2 week analysis
#time_window = 7.05 #half of 14+0.1 day time window
#time_mask = [time_window,GW_time+6.95] #[-0.1,14] day time window

#############################################################

############# LIGO SKYMAP ###############
fitsFile = '/data/user/rhussain/fitsFiles/GW%s_skymap.fits' % args.gw

### Read map and get probabilities
probs = hp.read_map(fitsFile)
nside = hp.pixelfunc.get_nside(probs)
probs = hp.pixelfunc.ud_grade(probs,nside_out=256,power=-2)
nside = 256


### Set up spatial prior to be used in scan
spatial_prior = SpatialPrior(probs,allow_neg=True,interpolated_ts_norm=True)
##########################################

llh = config(seasons,ncpu=7, days=5,timescramble=True,
             time_mask=time_mask,poisson=True)

### Set range of for loop so that multiple jobs 
### on the cluster will still all have unique 
### seedse
ntrials = args.ntrials
stop = ntrials * (args.pid+1)
start = stop-ntrials

tsList = []
for i in range(start,stop):
    val = llh.scan(0.0,0.0, scramble = True, seed = i,spatial_prior=spatial_prior,
                    time_mask = time_mask, pixel_scan=[nside,3.0])
    try:
        tsList.append(val['TS_spatial_prior_0'].max())
    except ValueError:
        tsList.append(-np.inf)

np.save('/home/rhussain/icecube/dump/gw%s/bkgTrials/maxTS_%s.npy' % (args.gw,args.pid), tsList)
