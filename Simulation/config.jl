module Parameters

       export fqtunerdir, wcsimdir, g4dir, rootdir, prod_basedir
       export verbose, nevents_per_task, nsubtasks, ntasks_per_job
       export base_mac, config_mac, task_template, job_template
       export config_variables, queue_command, max_jobs_queue

       softdir = "$(ENV["HOME"])/Software/"

       fqtunerdir   = "$(softdir)/HK_Software/WCSimFQTuner/install"
       wcsimdir     = "$(softdir)/HK_Software/WCSim/install-Linux_x86_64-gcc_9-python_3.10.13"
       g4dir        = "$(softdir)/Geant4/install"
       rootdir      = "$(softdir)/ROOT/install"
       prod_basedir = "$(ENV["LUSTRE"])/CProfiles/"

       verbose           = true
#        nevents_per_task  = 5 -> specified in each tuning step
       nsubtasks         = 100
       ntasks_per_job    = 10
       # Cherenkov profile
#        nevents_per_task  = 50
#        base_mac          = abspath("templates/cprofile_base.mac")
#        config_mac        = abspath("templates/cprofile_config.mac")
       # charge PDF
#        nevents_per_task  = 5
#        base_mac          = abspath("templates/charge_base_hk.mac")
#        config_mac        = abspath("templates/charge_config.mac")
       # Scattering tables
       nevents_per_task  = 50
       base_mac          = abspath("templates/stable_base_hk.mac")
       config_mac        = abspath("templates/stable_config.mac")

       task_template     = abspath("templates/task_template.sh")
       job_template      = abspath("templates/job_template.sh")

       # For Cherenkov profiles
#        config_variables  = Dict( "momentum" => range(50, 1000, step=10)
#                                , "particle" => ["e-"])
       # For scattering tables
       config_variables = Dict( "particle" => ["e-"])
       # For Charge PDF
#        config_variables = Dict( "true_charge" => vcat( range(0.1, 1.9, step=0.1)
#                                                      , range(2.0, 4.5, step=0.5)
#                                                      , range(5.0, 9.0, step=1)
#                                                      , range(  10, 19, step=2)
#                                                      , range(  20, 50, step=5)))
       # config_variables = Dict("particle" => ["e-"], "energy" => range(100, 1000, step=100))

       queue_command  = pipeline(`squeue -ah`, `wc -l`)
       max_jobs_queue = 100
end

# dont need to modify, it is needed by the generate_job_files.jl and launch_jobs.jl scripts
function sorter(fname)
       # returns the values of the config_variables for a given fname
       m    = match(r"_(.+)\.", basename(fname))
       @assert length(m.captures) == 1
       vars = split(m.captures[1], "_")
       return [(tryparse(Float64, v) === nothing) ? v : parse(Float64, v) for v in vars]
end
