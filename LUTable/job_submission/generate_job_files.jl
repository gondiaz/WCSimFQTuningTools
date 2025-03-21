using ArgParse

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "tasksdir"
            help = "Input directory"
            required = true
        "--job_template", "-j"
            help = ""
            arg_type = String
            default  = "job_template.sh"
        "--ntasks_per_job", "-n"
            help = ""
            arg_type = Int
            default  = 10
        "--verbose", "-v"
            help = "verbose"
            action = :store_true
    end
    return parse_args(s)
end

function sorter(fname)
    m    = match(r"_(.+)\.", basename(fname))
    @assert length(m.captures) == 1
    vars = split(m.captures[1], "_")
    return [(tryparse(Int    , v) === nothing) ?
           ((tryparse(Float64, v) === nothing) ? v : parse(Float64, v)) : parse(Int, v) for v in vars]
end

function main()
    
    parsed_args = parse_commandline()

    tasksdir    = joinpath(abspath(parsed_args["tasksdir"]), "")
    jobtemplate = abspath(parsed_args["job_template"])
    ntasks_pj   =         parsed_args["ntasks_per_job"]
    verbose     =         parsed_args["verbose"]

    basedir  = dirname(rstrip(tasksdir, '/'))
    jobs_dir = joinpath(basedir, "jobs")
    logs_dir = joinpath(basedir, "logs")
    errs_dir = joinpath(basedir, "errs")
    mkpath(jobs_dir)
    mkpath(logs_dir)
    mkpath(errs_dir)

    jobtemplate = 
    let fin = open(jobtemplate, "r")
        read(fin, String)
    end

    jobtemplate = replace( jobtemplate, "NTASKS" => ntasks_pj)

    task_files = readdir(tasksdir, join=true)
    task_files = [file for file in task_files if (match(r"task_(.+)_(\d*\.?\d+)_(\d+).sh", basename(file)) !== nothing)]
    sort!(task_files, by=sorter)

    # group tasks by energy
    grouped_tasks = Dict{Number, Vector{String}}()
    for file in task_files
        energy = sorter(file)[2]
        if haskey(grouped_tasks, energy)
            push!(grouped_tasks[energy], file)
        else
            grouped_tasks[energy] = [file]
        end
    end

    # loop in grouped tasks
    for (energy, task_files) in grouped_tasks
        sort!(task_files, by=sorter)

        # loop in ntasks_per_job partitions
        for files in Iterators.partition(task_files, ntasks_pj)
            # compute taskid of the first partition
            # assumes task files are named task_var1_var2_..._varN.sh
            jobid = join(sorter(files[1]), "_")

            # add line for each task 
            tasks = ""
            for f in files
                tasks *= "srun -n 1 --exclusive $f &\n"
            end
            tasks = chopsuffix(tasks, "&\n")
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
            if verbose println("Writting $job_fname") end
            write(job_fname, job)
        end
    end

end


main()

