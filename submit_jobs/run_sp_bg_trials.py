import pycondor

error = '/scratch/jthwaites/bg_error.error'
output = '/scratch/jthwaites/bg_output.out' 
log = '/scratch/jthwaites/bg_logs.log' 
submit = '/home/jthwaites/FRB/submit_jobs/submit.sub'

# bg_trials.py is a script that takes one argument,  
# the random seed, and runs many BG trials 
background_job = Job( 
    name='background_trials', 
    executable='/path/to/bg_trials_w_prior.py', 
    submit=submit, 
    error=error, 
    output=output, 
    log=log, 
    dag=dagman 
    ) 
"""
# Run many trials with different seeds 
for seed in range(10): 
    background_job.add_arg(f'{seed}')
"""
    
# Submit the DAG 
dagman.build_submit()