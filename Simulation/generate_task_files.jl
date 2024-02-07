include("config.jl")
using .Parameters

# define tasks dir and create path
tasks_dir = joinpath(abspath(prod_basedir), "tasks")
mkpath(tasks_dir)

# read macro filenames and sort them
macro_files = readdir(joinpath(abspath(prod_basedir), "macros"), join=true)
sort!(macro_files, by=sorter)

# read task template
tasktemplate = 
let fin = open(abspath(task_template), "r")
    read(fin, String)
end

tasktemplate =
replace( tasktemplate
       , "PROD_BASEDIR"       => abspath(prod_basedir)
       , "FQTUNER_INSTALLDIR" => abspath(fqtunerdir)
       , "WCSIM_INSTALLDIR"   => abspath(wcsimdir)
       , "G4_INSTALLDIR"      => abspath(g4dir)
       , "ROOT_INSTALLDIR"    => abspath(rootdir))

for macrofile in macro_files
    
    # taskid assumes macro files are named macro_var1_var2_..._varN.mac
    taskid = basename(macrofile)[7:end-4]

    # create and write task
    task = replace( tasktemplate
                  , "taskid"    => taskid
                  , "macrofile" => macrofile)

    task_fname = joinpath(tasks_dir, "task_$taskid.sh")
    if verbose println("Writing $task_fname") end
    write(task_fname, task)
    chmod(task_fname, 0o744)
end

