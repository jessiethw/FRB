""" Writes submit file for bg trials and sensitivities
    for point-source like FRBs

    Jessie Thwaites, Sept 2022
"""

import pycondor
import numpy as np
import pandas as pd
import argparse
import pwd
import os

""" Not needed rn
parser = argparse.ArgumentParser(
    description='Submit script')
parser.add_argument(
    '--output', type=str,
    default='/data/user/jthwaites/2022_GRECO_nova_analysis/csky_trials/systematics/',
    #default='/data/ana/analyses/NuSources/2022_GRECO_nova_analysis/' \
    #    'csky_trials/',
    help="Where to store output"
)
args = parser.parse_args()
"""

username = pwd.getpwuid(os.getuid())[0]
if not os.path.exists(f'/scratch/{username}/'):
    os.mkdir(f'/scratch/{username}/')
if not os.path.exists(f'/scratch/{username}/frb/'):
    os.mkdir(f'/scratch/{username}/frb/')

error = f'/scratch/{username}/frb/error'
output = f'/scratch/{username}/frb/output'
log = f'/scratch/{username}/frb/log'
submit = f'/scratch/{username}/frb/submit'

job = pycondor.Job(
    'sensitivity_frb_point_source',
    '../scripts/point_source_script.py',
    error=error,
    output=output,
    log=log,
    submit=submit,
    getenv=True,
    universe='vanilla',
    verbose=2,
    request_memory=6000,
    extra_lines=[
        'should_transfer_files = YES',
        'when_to_transfer_output = ON_EXIT']
    )

frbs = pd.read_csv('/home/jthwaites/FRB/catalog/frbs_all.csv')
##want to remove those FRBs with priors
localization_files=os.listdir('/data/user/jthwaites/chime_localization_data/')
spatial_priors=[]
[spatial_priors.append(filename[0:12]) for filename in localization_files]

for (index,frb) in frbs.iterrows():
    if frb['src'] in spatial_priors:
        frbs=frbs.drop([index])

repeater_names, rep_index= np.unique(frbs['src'][frbs['repeater']==True], return_index=True)
all_repeats = frbs[frbs['repeater']==True].index
first_repeat= frbs['src'][frbs['repeater']==True].index[rep_index]
for i in all_repeats:
    if i not in first_repeat:
        frbs=frbs.drop([i])

frbs = frbs.reset_index()

for frb_name in frbs['src']:
    for time_window in np.logspace(0,6,num=7):
        if time_window<10.:
            ntrials=1000000
        elif time_window<10000.:
            ntrials=100000
        else:
            ntrials=10000
        job.add_arg( f'--source={frb_name} --deltaT={time_window} --ntrials={ntrials}')

dagman = pycondor.Dagman(
    'FRB_sens_trials_point_source',
    submit=submit, verbose=2)
dagman.add_job(job)
dagman.build_submit()
