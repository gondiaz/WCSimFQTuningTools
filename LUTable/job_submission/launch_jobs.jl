using ArgParse

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "jobs_dir"
            help = "Input directory"
            required = true
        "--max_jobs_queue", "-m"
            help = ""
            arg_type = Int
            default  = 100
        "--verbose", "-v"
            help = "verbose"
            action = :store_true
    end
    return parse_args(s)
end

# hardcoded, move to arguments
queue_command = pipeline(`squeue -ah`, `wc -l`) 
get_njobs_in_queue() = parse(Int, replace(readchomp(queue_command), " " => ""))

function sorter(fname)
    m    = match(r"_(.+)\.", basename(fname))
    @assert length(m.captures) == 1
    vars = split(m.captures[1], "_")
    return [(tryparse(Int    , v) === nothing) ?
           ((tryparse(Float64, v) === nothing) ? v : parse(Float64, v)) : parse(Int, v) for v in vars]
end


function main()
    
    parsed_args = parse_commandline()

    jobs_dir       = abspath(parsed_args["jobs_dir"])
    max_jobs_queue =         parsed_args["max_jobs_queue"]
    verbose        =         parsed_args["verbose"]

    job_files = readdir(jobs_dir, join=true)
    sort!(job_files, by=sorter)

    for job in job_files

        free_queue = @task while 
            get_njobs_in_queue() > max_jobs_queue
            sleep(10)
        end

        schedule(free_queue)
        wait(free_queue)

        run(`sbatch $job`)

        if verbose println("$job launched") end
    end
end

main()

