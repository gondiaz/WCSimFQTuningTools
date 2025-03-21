using ArgParse

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "indir"
            help = "Input directory"
            required = true
        "outdir"
            help = "Output directory"
            required = true
        "--task_template", "-t"
            help = ""
            arg_type = String
            default  = "task_template.sh"
        "--files_per_task", "-n"
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

    indir          = joinpath(abspath(parsed_args["indir"]) , "")
    outdir         = joinpath(abspath(parsed_args["outdir"]), "")
    task_template  = abspath(parsed_args["task_template"])
    files_per_task = parsed_args["files_per_task"]
    verbose        = parsed_args["verbose"]

    if !isdir(indir) error("Input directory does not exists") end

    if verbose
        println("Input directory : ", indir)
        println("Output directory: ", outdir)
    end 

    tasktemplate = 
    let fin = open(task_template, "r")
        read(fin, String)
    end

    inputfiles = readdir(indir, join=true)
    inputfiles = [file for file in inputfiles if (match(r"out_(.+)_(\d*\.?\d+)_(\d+).h5", basename(file)) !== nothing)]
    sort!(inputfiles, by=sorter)

    mkpath(joinpath(outdir, "out"))
    mkpath(joinpath(outdir, "tasks"))

    grouped_inputfiles = Iterators.partition(inputfiles, files_per_task)

    for (i, infiles) in enumerate(grouped_inputfiles)

        taskid = join(sorter(infiles[1]), "_")
        
        task = replace( tasktemplate
                      , "input_files" => join(infiles, " ")
                      , "outfilename" => joinpath(outdir, "out/qdistributions_$i.h5"))

        task_fname = joinpath(outdir, "tasks", "task_$taskid.sh")
        if verbose println("Writting $task_fname") end
        write(task_fname, task)
        chmod(task_fname, 0o744)
    end

end


main()