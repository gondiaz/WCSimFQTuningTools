"""
    merge_cprofiles_parallel.py
"""
import re
import glob
import argparse
import itertools
import multiprocessing
import concurrent.futures

from itertools import groupby
from os        import cpu_count
from os.path   import expandvars, join, basename

import uproot
import ROOT
import numpy as np


def process_momentum(pid, momentum, files, htype, lock, verbose=False):

    if verbose: 
        print(f"Processing {momentum} MeV/c files for {pid}...".ljust(50))

    nevents = 0
    # start with array to accumulate histogram entries
    with uproot.open(files[0]) as fin:
        nevents += fin["nevents"].counts()[0]
        h, thbins, sbins = fin[f"{htype}g"].to_numpy()
        # transform distance units from mm to cm
        sbins = np.divide(sbins, 10)
    # loop through the rest of files
    for file in files[1:]:
        with uproot.open(file) as fin:
            nevents += fin["nevents"].counts()[0]
            h_, _, _ = fin[f"{htype}g"].to_numpy()
            h += h_ # accumulate entries

    th2d = ROOT.TH2D( f"g_{momentum}", f"g_{momentum}"
                    , len(thbins)-1, thbins
                    , len(sbins) -1, sbins)
    for ix, iy in itertools.product(range(1, len(thbins)), range(1, len(sbins))):
        th2d.SetBinContent(ix, iy, h[ix-1, iy-1])

    with lock:
        fout = ROOT.TFile(f"cprofiles_{pid}_merged.root", "UPDATE")
        fout.WriteObject(th2d, f"g_{momentum}")
        fout.Close()

    return nevents


def main():

    ############ Program arguments ############
    parser = argparse.ArgumentParser( prog        = f"{basename(__file__)}"
                                    , description = "description"
                                    , epilog      = "Text at the bottom of help")

    parser.add_argument("-v", "--verbose", action="store_true")

    parser.add_argument( "indir", type=str, nargs=1, help = "directory containing produced files")
    parser.add_argument(   "pid", type=str,          help = "particle type (eg e-, mu+...", default="all")
    parser.add_argument(  "type", type=str, nargs=1, help = "type of cherenkov profile, true (tr) or weighted (wt)") # modify this, it can be tr/wt

    args = parser.parse_args()
    ##########################################

    # check pid and define output file
    if args.pid == "all":
        if args.verbose:
            print("Merging all files in %s", args.indir[0])
        output_file = "cprofiles_all_merged.root"
    elif args.pid in ["e-", "e+", "mu-", "mu+", "pi-", "pi+"]:
        if args.verbose: 
            print(f"Merging {args.pid} files")
        output_file = f"cprofiles_{args.pid}_merged.root"
    else:
        raise Exception(f"Unknown parameter: {args.pid}, should be e-, e+, mu-, mu+, pi-, pi+")

    # check if type is tr or wt
    if args.type[0] not in ["tr", "wt"]:
        raise Exception("Type must be one of tr or wt")

    # sorter function
    get_momentum_and_index = lambda fname: list(map(float, re.findall(r'\d+(?:\.\d+)?', basename(fname))[:2]))

    # get all files in input dir, filter out the '_flat.root' and sort them based on (momentum, index)
    infiles = glob.glob(join(expandvars(args.indir[0]), "*"))
    if args.pid == "all":
        infiles = [f for f in infiles if re.match("out_([^_]+)_\d+(?:\.\d+)?_\d+.root", basename(f))]
    else:
        infiles = [f for f in infiles if re.match(f"out_{args.pid}_\d+(?:\.\d+)?_\d+.root", basename(f))]
    infiles = sorted(infiles, key=get_momentum_and_index)

    # split input files in momentum groups
    momenta = np.unique([get_momentum_and_index(f)[0] for f in infiles])
    groups = [list(group) for key, group in groupby(infiles, key=lambda x: get_momentum_and_index(x)[0])]

    # create output file
    fout = ROOT.TFile(output_file, "RECREATE")
    fout.Close()

    # paralellize loop over momenta to fill histograms
    results = []
    with multiprocessing.Manager() as manager:
        lock = manager.Lock()
        with concurrent.futures.ProcessPoolExecutor(max_workers=cpu_count()) as executor:
            future_to_momentum = {executor.submit(process_momentum, args.pid, momentum, files, args.type[0], lock, args.verbose): momentum for momentum, files in zip(momenta, groups)}

            # get results once all the tasks are finished
            for future in concurrent.futures.as_completed(future_to_momentum):
                nevents = future.result()
                momentum = future_to_momentum[future]
                results.append((momentum, nevents))

    # sort by momentum values
    results = np.array(results, dtype=[("momentum", float), ("nevents", float)])
    results.sort(order="momentum")

    # get values
    momenta = results["momentum"].copy()
    nevents = results["nevents"] .copy()

    fout = ROOT.TFile(output_file, "UPDATE")

    # save number of events per momentum
    g = ROOT.TGraph(len(momenta), momenta, nevents)
    g.SetTitle("Number of events per momentum")
    fout.WriteObject(g, "g_nevents")

    # close output file
    fout.Close()
    return


if __name__ == "__main__":
    main()
