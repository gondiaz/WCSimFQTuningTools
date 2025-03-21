import re
import sys
import glob
import argparse
import warnings
import concurrent.futures

from os.path   import expandvars, basename, join
from itertools import groupby

import ROOT
import tables as tb
import pandas as pd
import numpy  as np

def get_total_charge(filename):
    rootf = ROOT.TFile(filename, "read")
    tree  = rootf.GetKey("wcsimT").ReadObj()
    tree.GetBranch("wcsimrootevent").SetAutoDelete(True)
    nevents = tree.GetEntries()

    qs = np.zeros(nevents)
    for event in range(nevents):
        tree.GetEvent(event)
        trigger  = tree.wcsimrootevent.GetTrigger(0)
        qs[event] = trigger.GetSumQ()
    rootf.Close()
    return qs


def process_energy(energy, filenames, verbose):
    """TODO: add description"""
    
    if verbose: print(f"--> Processing energy = {energy}".ljust(50))

    qs = [None]*len(filenames)
    for i, filename in enumerate(filenames):
        qs[i] = get_total_charge(filename)
    qs = np.concatenate(qs)

    if verbose: print(f"<-- Done for energy = {energy}".ljust(50))
    return qs


def main():

    ############ Program arguments ############
    parser = argparse.ArgumentParser( prog        = "get_conversion_distributions"
                                    , description = ""
                                    , epilog      = """""")
    
    parser.add_argument("-v", "--verbose", action="store_true")
    parser.add_argument( "indir", type=str, nargs=1, help = "directory containing input files")
    parser.add_argument( "--wcsimlib",   type=str, nargs="?", help = "WCSim lib path", default="$HOME/Software/WCSim/install/lib")
    
    args = parser.parse_args()
    ##########################################

    # Load WCSimRoot.so
    ROOT.gSystem.AddDynamicPath(expandvars(args.wcsimlib))
    ROOT.gSystem.Load          ("libWCSimRoot.dylib" if sys.platform == "darwin" else "libWCSimRoot.so")

    get_energy_and_index = lambda fname: list(map(float, re.findall(r'\d+(?:\.\d+)?', basename(fname))[:2]))

    # get all files in input dir, filter out the '_flat.root' and sort them based on (p, index)
    infiles = glob.glob(join(args.indir[0], "*"))
    infiles = [f for f in infiles if re.match("out_([^_]+)_\d+(?:\.\d+)?_\d+.root", basename(f))]
    infiles = sorted(infiles, key=get_energy_and_index)

    # split input files in mu groups
    energies = np.unique([get_energy_and_index(f)[0] for f in infiles])
    groups = [list(group) for key, group in groupby(infiles, key=lambda x: get_energy_and_index(x)[0])]

    # loop over energy and get data
    results = []
    with concurrent.futures.ProcessPoolExecutor() as executor:
        # parallelize
        future_to_E = {executor.submit(process_energy, E, files, args.verbose): E for E, files in zip(energies, groups)}

        # get results once all the tasks are finished
        for future in concurrent.futures.as_completed(future_to_E):
            qs = future.result()
            E  = future_to_E[future]
            results.append((E, qs))

    # save results
    warnings.filterwarnings("ignore", category=tb.NaturalNameWarning)

    if args.verbose: print("Saving distributions...")
    with tb.open_file("charge_distributions.h5", "w") as f:
        for E, qs in results:
            f.create_array(f.root, f"E_{E}", qs, f"charges for E = {E} MeV")
    return



if __name__ == "__main__":
    main()