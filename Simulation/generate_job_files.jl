include("config.jl")
using .Parameters

# define paths for jobs, logs and errors
jobs_dir    = joinpath(abspath(prod_basedir), "jobs")
logs_dir    = joinpath(abspath(prod_basedir), "logs")
errs_dir    = joinpath(abspath(prod_basedir), "errs")
mkpath(jobs_dir)
mkpath(logs_dir)
mkpath(errs_dir)

# read task filenames and sort them
# notice that the default path for tasks is used
task_files= readdir(joinpath(abspath(prod_basedir), "tasks"), join=true)
sort!(task_files, by=sorter)

# read and replace ntasks_per_job jobtemplate
jobtemplate = 
let fin = open(abspath(job_template), "r")
    read(fin, String)
end
jobtemplate = replace(jobtemplate, "NTASKS" => ntasks_per_job)

# iterate over taskfile partitions and create jobs
for taskfiles in Iterators.partition(task_files, ntasks_per_job)
    
    # compute taskid of the first partition
    # assumes task files are named task_var1_var2_..._varN.sh
    jobid = basename(taskfiles[1])[6:end-3]

    # add line for each task 
    tasks = ""
    for f in taskfiles
        tasks *= "srun -n 1 --exclusive $f &\n"
    end
#     tasks = chopsuffix(tasks, "&\n")
    tasks *= "\nwait"
    
    # replace variables into job and write it
    job = 
    replace( jobtemplate
           , "JOBNAME"       => "$jobid"
           , "LOGFILENAME"   => joinpath(logs_dir, "$jobid.log")
           , "ERRORFILENAME" => joinpath(errs_dir, "$jobid.err")
           , "jobid"         => jobid
           , "TASKS"         => tasks)
    job_fname = joinpath(jobs_dir, "job_$jobid.sh")
    if verbose println("Writing $job_fname") end
    write(job_fname, job)
end

