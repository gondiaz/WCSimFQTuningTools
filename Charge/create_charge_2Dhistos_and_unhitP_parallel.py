"""
    create_charge_2Dhistos_and_unhitP_parallel.py
"""
import sys
import re
import glob
import argparse
import itertools
import concurrent.futures

from os      import cpu_count
from os.path import expandvars, join, basename, exists

import ROOT
import numpy  as np
from scipy.optimize import curve_fit

def get_DigiHitQs_and_nHits(filename):
    """Get DigiHit charges and number of hits for each PMT and event in the root file"""
    rootf = ROOT.TFile(filename, "read")
    tree  = rootf.GetKey("wcsimT").ReadObj()
    tree.GetBranch("wcsimrootevent").SetAutoDelete(True) # required to dealloc
    nevents = tree.GetEntries()

    Qs   = []
    nHit = []
    for event in range(nevents):
        tree.GetEvent(event)
        trigger  = tree.wcsimrootevent.GetTrigger(0)
        # fill DigiHit charges
        DigiHits = trigger.GetCherenkovDigiHits()
        for dhit in DigiHits: Qs.append(dhit.GetQ())
        # fill nHit
        nHit.append(trigger.GetNcherenkovdigihits())

    rootf.Close()
    return np.array(Qs), np.array(nHit)


def get_nPMTs(filename):
    """Get number of PMTs in the detector"""
    rootf = ROOT.TFile(filename, "read")
    tree  = rootf.GetKey("wcsimGeoT").ReadObj()
    tree.GetBranch("wcsimrootgeom").SetAutoDelete(True) # required to dealloc
    tree.GetEntry(0)
    geom  = tree.wcsimrootgeom
    return geom.GetWCNumPMT()


def process_mu(mu, files, qbins, nPMTs, verbose=False):
    """Get q distribution and hit prob. for list of files (which are supposed to have same mu value)"""
    if verbose:
        print(f"--> Processing mu = {mu}".ljust(50))

    Qs = []
    nHits = []
    for filename in files:
        Qs_, nHits_ = get_DigiHitQs_and_nHits(filename)
        Qs.extend(Qs_)
        nHits.extend(nHits_)

    Qs = np.array(Qs)
    Phit = np.mean(nHits) / nPMTs
    hq, _ = np.histogram(Qs, bins=qbins, density=True)

    if verbose:
        print(f"<-- Done for mu = {mu}".ljust(50))
    return mu, Phit, hq


def main():

    ############ Program arguments ############
    parser = argparse.ArgumentParser( prog        = f"{basename(__file__)}"
                                    , description = "Compute 2D histogram (predicted vs measured charge) and predicted Q vs unhitP"
                                    , epilog      = """""")

    parser.add_argument("-v", "--verbose", action="store_true")
    parser.add_argument(   "indir", type=str, nargs=1, help = "directory containing input files")
    parser.add_argument( "--qbins", type=str, nargs=1, help = "qbins file", default="qbins_wcte.txt")
    parser.add_argument( "--wcsimlib",   type=str, nargs="?", help = "WCSim lib path")

    args = parser.parse_args()
    ##########################################

    # Load WCSimRoot.so
    ROOT.gSystem.AddDynamicPath(expandvars(args.wcsimlib))
    ROOT.gSystem.Load          ("libWCSimRoot.dylib" if sys.platform == "darwin" else "libWCSimRoot.so")

    # in what follows, "mu" stands for predicted and "q" for measured charge
    # sorter function
    get_mu_and_index = lambda fname: list(map(float, re.findall(r'\d+(?:\.\d+)?', basename(fname))[:2]))

    # get all files in input dir, filter out the '_flat.root' and sort them based on (mu, index)
    infiles = glob.glob(join(args.indir[0], "*"))
    infiles = [f for f in infiles if re.match("out_\d+(?:\.\d+)?_\d+.root", basename(f))]
    infiles = sorted(infiles, key=get_mu_and_index)

    # get number of PMTs in the detector (used to compute the hit probability)
    nPMTs = get_nPMTs(infiles[0])

    # split input files in mu groups
    mus = np.unique([get_mu_and_index(f)[0] for f in infiles])
    groups = [list(group) for key, group in itertools.groupby(infiles, key=lambda x: get_mu_and_index(x)[0])]

    # read qbins from file
    exists(args.qbins[0])
    qbins = np.loadtxt(args.qbins[0])
    # loop over mu and fill histograms for each one
    results = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=cpu_count()) as executor:
        # parallelize
        future_to_mu = {executor.submit(process_mu, mu, files, qbins, nPMTs, args.verbose): mu for mu, files in zip(mus, groups)}

        # get results once all the tasks are finished
        for future in concurrent.futures.as_completed(future_to_mu):
            mu, Phit, hq = future.result()
            mu = future_to_mu[future]
            results.append((mu, Phit, hq))

    results = np.array(results, dtype=[("mu", float), ("Phit", float), ("hq", np.ndarray)])
    results.sort(order="mu")

    # get 2D histogram and hit prob.
    Hq   = np.vstack(results["hq"])
    Phit = results["Phit"]

    # Create output file
    fout = ROOT.TFile("charge2D_and_unhit.root", "RECREATE")

    # Save 2D charge histogram
    # compute mu bins such that mus array lie at the center
    mubins = (mus[:-1] + mus[1:])/2.
    mubins = np.insert(mubins,           0, mus[0]  - (mubins[0]  - mus[0]))
    mubins = np.insert(mubins, len(mubins), mus[-1] - (mubins[-1] - mus[-1]))

    # define and fill histogram
    th2d = ROOT.TH2D( "charge2D", "charge2D"
                    , len(mubins)-1, mubins
                    , len(qbins) -1, qbins)
    for ix, iy in itertools.product(range(1, len(mubins)), range(1, len(qbins))):
        th2d.SetBinContent(ix, iy, Hq[ix-1, iy-1])

    # save to file
    fout.WriteObject(th2d, "charge2D")

    # Save UnHit probability, values and fit result
    tgraph = ROOT.TGraph(len(mus), mus, 1.-Phit)
    fout.WriteObject(tgraph, "PUnhit")

    # Fit UnHit probability
    # define function to fit
    def PUnhit_func(x, a1, a2, a3):
        return (1. + a1*x + a2*x**2 + a3*x**3)*np.exp(-x)

    # perform fit
    popt, _ = curve_fit(PUnhit_func, mus, 1.-Phit, bounds=(0, [np.inf, np.inf, np.inf]))

    # Save fit parameters
    th1d = ROOT.TH1D("hPunhitPar", "c_n for P(unhit|#mu)", 10, 0.5, 10.5)
    for i, par in enumerate(popt, 1):
        th1d.SetBinContent(i, par)
    fout.WriteObject(th1d, "hPunhitPar")

    return


if __name__ == "__main__":
    main()
