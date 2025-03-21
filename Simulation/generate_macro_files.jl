include("config.jl")
using .Parameters

# function that replaces fullpaths by basenames in dictionary
function full_to_base(dict::Dict)
    return Dict(k => (v isa String && ispath(v) ? basename(v) : v) for (k, v) in dict)
end

# create macro and out directories
macro_dir = joinpath(abspath(prod_basedir), "macros")
out_dir   = joinpath(abspath(prod_basedir), "out")
mkpath(macro_dir)
mkpath(out_dir)

# read config macro template
mactemplate = 
let f = open(abspath(config_mac), "r")
    read(f, String)
end
# get the variables that the macro file expects to be replaced (those starting with $)
macvars = [replace(m.match, "\$" =>"") for m in eachmatch(r"\$\S+", mactemplate)]

# get names of config variables
ks = collect(keys(config_variables))

# check that all ks are in macvars
@assert issubset(ks, macvars)

# define macro and root filenames using the config_variables
# such that macro_var1_var2_var3..._subtask.mac
macro_fname = joinpath(macro_dir, string("macro_", join(ks, "_"), "_subtask.mac"))
out_fname   = joinpath(  out_dir, string(  "out_", join(ks, "_"), "_subtask.root"))

global rseed = 1
# loop over all possible combinations of config_variables values
for vals in Iterators.product(values(config_variables)...)
    # define dict with (config-variable, value)
    d = Dict{String, Any}(zip(ks, vals))

    if verbose println("Writting files for $(join(["$key = $value" for (key, value) in d], ", "))") end

    # subtask for each set of values
    for subtask in range(1, nsubtasks)
        
        # add values to dic
        d["rseed"]          = rseed
        d["subtask"]        = subtask
        d["out_root_fname"] = replace(replace(out_fname, collect(full_to_base(d))...), " " => "_")

        # add macvars to dictionary
        for var in setdiff(macvars, collect(keys(d)), ["out_root_fname"])
            d[var] = eval(Symbol(var))
        end
        
        # replace variable values into macro and write it
        mac = replace(mactemplate, [Pair("\$" * k, v) for (k,v) in d]...)

        # write macro file
        write(replace(replace(macro_fname, collect(full_to_base(d))...), " " => "_"), mac)

        global rseed += 1
    end
end
