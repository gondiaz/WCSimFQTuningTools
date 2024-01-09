## **Simulations**

This directory contains scripts to generate and run the simulations.The simulations use **WCSimFQTuner**, which requires a macro of configuration parameters. This utility creates the *macros* and *tasks* for each desired configuration and the *jobs* to be run.

The simulations are configured by the `config.jl` file, with the parameters described below:

- **fqtunerdir, wcsimdir, g4dir, rootdir**: paths to the different software components needed to run the simulations. In the **task_template** provided with this repository, it is assumed that **ROOT** was instaled using a conda environment named **wcte**.

- **prod_basedir**: directory where macro, task, job and simulation output files are saved.

- **verbose, nevents_per_task, nsubtasks, ntasks_per_job**: names are self-explanatory. **nsubtasks** represents the number of identical tasks for each task configuration. Random seeds for the simulations are unique for each task configuration and subtask. **ntasks_per_job** is the number of parallelized tasks to launch in each job.

- **base_mac, config_mac, task_template, job_template**: look at the provided examples to understand the meaning of each parameter.

- **config_variables**: pairs of variable name and values. A macro will be created for each variable values. The variable names must appear in the **config_mac** with a preceeding $ symbol.

- **queue_command**: the result of this command must return the number of currently running and queued jobs.

- **max_jobs_queue**: the job launcher halts if the number of jobs is higher than this value. It continues launching jobs once the number of jobs gets lower.

The procedure to create the macro files, tasks and jobs is the following:
1) First modify the `config.jl` file parameters as required.
2) Modify the macro, task and job templates as desired. Notice that we define two types of macros, the *base* and the *config* macros. The *base* macro defines the set of fixed parameters for the particular simulation type, and the *config* the *settable* parameters, which are defined with a preceding `$` symbol.
The *_base* macros correspond to WCTE simulations used for the fiTQun tuning.
For your own geometry (HK, SK, other), you might need to define your own macros and change these in the *config.jl* file.
3) Run the scripts in the following order to produce the jobs
    1) `julia generate_macro_files.jl`
    2) `julia generate_task_files.jl`
    3) `julia generate_job_files.jl`
4) Launch the jobs through `julia launch_jobs.jl`

To run the simulations interactively, see **Simulations/local**.

**TODO**: add documentation about the different simulations required for the tuning files