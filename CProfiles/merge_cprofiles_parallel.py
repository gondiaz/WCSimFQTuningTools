import os
import re
import glob
import argparse
import itertools
import uproot
import ROOT
import numpy  as np

import multiprocessing
import concurrent.futures
from itertools import groupby
from os.path   import expandvars, realpath, join, basename


def process_momentum(pid, momentum, files, htype, lock, verbose=False):

    if verbose: print(f"Processing {momentum} MeV/c files for {pid}...".ljust(50))

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
            
    # mean nphotons per momentum
    nphotons = h.sum()/nevents

    # create, fill and save 2D histogram for this momentum
    thbinw = thbins[1] - thbins[0]
    sbinw =   sbins[1] -  sbins[0]
    h *= (1/h.sum())*(1/(thbinw*sbinw)) # normalize histogram

    th2d = ROOT.TH2D( f"g_{momentum}", f"g_{momentum}"
                    , len(thbins)-1, thbins
                    , len(sbins) -1, sbins)
    for ix, iy in itertools.product(range(1, len(thbins)), range(1, len(sbins))): th2d.SetBinContent(ix, iy, h[ix-1, iy-1])

    with lock:
        fout = ROOT.TFile(f"cprofiles_{pid}_merged.root", "UPDATE")
        fout.WriteObject(th2d, f"g_{momentum}")
        fout.Close()

    # compute s distance at which 90% of photons were emitted
    # (needed in the parabolic approximation for the angular response)
    sproj = h.sum(axis=0)*thbinw # projected pdf in s
    assert len(sproj) == len(sbins) - 1

    scum = np.cumsum(sproj)*sbinw
    sindex = np.argwhere(scum>0.9).flatten()[0]
    sthr = (sbins[sindex] + sbins[sindex+1])/2.

    return nphotons, sthr


def main():

    ############ Program arguments ############
    parser = argparse.ArgumentParser( prog        = "merge_cprofiles"
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
        if args.verbose: print("Merging all files")
        output_file = "cprofiles_merged.root"
    elif args.pid in ["e-", "e+", "mu-", "mu+", "pi-", "pi+"]:
        if args.verbose: print(f"Merging {args.pid} files")
        output_file = f"cprofiles_{args.pid}_merged.root"
    else:
        print(f"Unknown parameter: {args.pid}, should be e-, e+, mu-, mu+, pi-, pi+")
        exit(1)

    # check if type is tr or wt
    if args.type[0] not in ["tr", "wt"]: raise Exception("Type must be one of tr or wt")

    # sorter function
    get_momentum_and_index = lambda fname: list(map(float, re.findall(r'\d+(?:\.\d+)?', basename(fname))[:2]))

    # get all files in input dir, filter out the '_flat.root' and sort them based on (momentum, index)
    infiles = glob.glob(join(args.indir[0], "*"))
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
        with concurrent.futures.ProcessPoolExecutor(max_workers=os.cpu_count()) as executor:
            future_to_momentum = {executor.submit(process_momentum, args.pid, momentum, files, args.type[0], lock, args.verbose): momentum for momentum, files in zip(momenta, groups)}

            # get results once all the tasks are finished
            for future in concurrent.futures.as_completed(future_to_momentum):
                nphotons, sthr = future.result()
                momentum = future_to_momentum[future]
                results.append((momentum, nphotons, sthr))

    # sort by momentum values
    results = np.array(results, dtype=[("momentum", float), ("nphotons", float), ("sthrs", float)])
    results.sort(order="momentum")

    # get values
    momenta  = results["momentum"].copy()
    nphotons = results["nphotons"].copy()
    sthrs    = results["sthrs"]   .copy()

    fout = ROOT.TFile(output_file, "UPDATE")

    # save mean number of photons per momentum
    g = ROOT.TGraph(len(momenta), momenta, nphotons)
    g.SetTitle("Mean number of photons per event")
    fout.WriteObject(g, "gNphot")

    # save s thresholds
    g = ROOT.TGraph(len(momenta), momenta, sthrs)
    g.SetTitle("Length of track")
    fout.WriteObject(g, "gsthr")

    # close output file
    fout.Close()
    return


if __name__ == "__main__":
    main()