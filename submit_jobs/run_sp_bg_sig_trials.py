""" Writes submit file for bg trials and sensitivities
    for (CHIME) spatial prior FRBs

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
if not os.path.exists(f'/scratch/{username}/frb/priors/'):
    os.mkdir(f'/scratch/{username}/frb/priors/')

error = f'/scratch/{username}/frb/priors/error'
output = f'/scratch/{username}/frb/priors/output'
log = f'/scratch/{username}/frb/priors/log'
submit = f'/scratch/{username}/frb/priors/submit'

dagman = pycondor.Dagman(
    'FRB_bg_sens_trials_prior',
    submit=submit, verbose=2)

bg_job = pycondor.Job(
    'bg_trials_frb_prior',
    '../scripts/bg_trials_w_prior.py',
    error=error,
    output=output,
    log=log,
    submit=submit,
    getenv=True,
    universe='vanilla',
    verbose=2,
    request_memory=10000,
    extra_lines=[
        'should_transfer_files = YES',
        'when_to_transfer_output = ON_EXIT'],
    dag=dagman
    )
sens_job = pycondor.Job(
    'sensitivity_trials_frb_prior',
    '../scripts/sens_w_prior.py',
    error=error,
    output=output,
    log=log,
    submit=submit,
    getenv=True,
    universe='vanilla',
    verbose=2,
    request_memory=10000,
    extra_lines=[
        'should_transfer_files = YES',
        'when_to_transfer_output = ON_EXIT'],
    dag=dagman
    )

frbs = pd.read_csv('/home/jthwaites/FRB/catalog/spatial_priors_frbs.csv')

repeater_names, rep_index= np.unique(frbs['src'][frbs['repeater']==True], return_index=True)
all_repeats = frbs[frbs['repeater']==True].index
first_repeat= frbs['src'][frbs['repeater']==True].index[rep_index]
for i in all_repeats:
    if i not in first_repeat:
        frbs=frbs.drop([i])

for frb_name in frbs['src']:
    for time_window in np.logspace(0,6,num=7):
        if time_window<100.:
            ntrials=100000
        else:
            ntrials=10000
        bg_job.add_arg( f'--source={frb_name} --deltaT={time_window} --ntrials={ntrials}')
        sens_job.add_arg( f'--source={frb_name} --deltaT={time_window}')

sens_job.add_parent(bg_job)

dagman.build_submit()
