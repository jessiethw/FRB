""" Writes submit file for bg trials with spatial
    Jessie Thwaites, 1/4/22
"""

import pycondor

error = '/scratch/jthwaites/bg_error.error'
output = '/scratch/jthwaites/bg_output.out' 
log = '/scratch/jthwaites/bg_logs.log' 
submit = '/home/jthwaites/FRB/submit_jobs/submit.sub'

dagman = pycondor.Dagman(
    'csky_sp_frb_bg_jobs', submit=submit, verbose=2
    )

background_job = pycondor.Job( 
    name='background_trials', 
    executable='/home/jthwaites/FRB/bg_trials_w_prior.py', 
    submit=submit, 
    error=error, 
    output=output, 
    log=log, 
    #getenv=True,
    #universe='vanilla',
    verbose=2,
    request_memory=8000,
    request_cpus=15,
    extra_lines=[
        'should_transfer_files = YES',
        'when_to_transfer_output = ON_EXIT'],
    dag=dagman 
    ) 

background_job.add_arg(f'--FRB_src=FRB20190416A --ntrials={10000} '
                        + f'--deltaT={1000}')

dagman.add_job(background_job)
dagman.build_submit()
