import sys
import argparse
import uproot
import ROOT
import numpy  as np

from scipy.optimize import curve_fit

from os.path   import expandvars


def get_fitting_function(mulow, muup):
    """Returns function of mu to be fitted with mu bounds (mulow, muup)"""
    def func(mu, *pars):
        npars = len(pars)
        if (mulow<=mu) & (mu<=muup):
            return np.poly1d(np.flip(pars))(mu)
        elif (mu<mulow):
            c0 = sum([        pars[i]*mulow**i        for i in range(0, npars)])
            c1 = sum([i*      pars[i]*mulow**(i-1)    for i in range(1, npars)])
            c2 = sum([i*(i-1)*pars[i]*mulow**(i-2)/2. for i in range(2, npars)])
            if (c2>0.): c2 = 0
            dmu = (mu-mulow)
            return c0 + dmu*(c1 + c2*dmu)
        elif (muup<mu):
            c0 = sum([        pars[i]*muup**i        for i in range(0, npars)])
            c1 = sum([i*      pars[i]*muup**(i-1)    for i in range(1, npars)])
            c2 = sum([i*(i-1)*pars[i]*muup**(i-2)/2. for i in range(2, npars)])
            if (c2>0.): c2 = 0
            dmu = (mu-muup)
            return c0 + dmu*(c1 + c2*dmu)
    return np.vectorize(func)


def main():

    ############ Program arguments ############
    parser = argparse.ArgumentParser( prog        = "fit_charge_pdfs"
                                    , description = ""
                                    , epilog      = """""")

    parser.add_argument("-v", "--verbose", action="store_true")
    parser.add_argument(  "--infile", type=str, nargs="?", help = "charge 2D file", default="charge2D_and_unhit.root")
    parser.add_argument(  "--npars", type=int  , nargs="+", help = "number of parameters per range")
    parser.add_argument("--qranges", type=float, nargs="+", help = "range limits for q")

    args = parser.parse_args()
    ##########################################

    # Config variables
    infilename = args.infile
    npars      = args.npars
    qranges    = args.qranges

    nqranges = len(npars)
    assert nqranges == len(qranges) - 1

    if args.verbose:
        print("-------------------------")
        print(f"Using {len(npars)} q-ranges")
        print("-------------------------")
        for i, n in enumerate(npars):
            print(f"Range {i} between [{qranges[i]}, {qranges[i+1]}] with {n} parameters")

    # in what follows, "mu" stands for predicted and "q" for measured charge

    # read 2D histogram
    with uproot.open(infilename) as file:
        Hq, mubins, qbins = file["charge2D"].to_numpy()

    # compute bin centers
    mus = (mubins[1:] + mubins[:-1])/2.
    qs  = ( qbins[1:] +  qbins[:-1])/2.

    # output file
    outfilename = "fitted_cpdf.root"
    fout = ROOT.TFile(outfilename, "RECREATE")

    # define q fit ranges and save to output
    hCPDFrange = ROOT.TH1D("hCPDFrange", "", nqranges, np.array(qranges))
    for i, n in enumerate(npars, 1):
        hCPDFrange.SetBinContent(i, n)
    fout.WriteObject(hCPDFrange, "hCPDFrange")

    # loop in q ranges
    # - gParams: for each range and parameter in this range: q vs parameter
    # - gmuthrs: for each range and parameter in this range: q vs low and q vs high thresholds
    globalpi = 1 # global index of the qbin (unlike in range index): represents the qbin index
    for rang in range(nqranges):
        if args.verbose:
            print()
            print(f"Processing range {rang}")
            print("-----------------------")

        # select qs in range
        qmin, qmax = qranges[rang], qranges[rang+1]
        selqbins = np.argwhere((qmin <= qs) & (qs < qmax)).flatten()

        # graphs to fill for this range
        npar = npars[rang]
        gParams = [ROOT.TGraph() for par in range(npar)] # one graph per parameter
        gmuthrs = [ROOT.TGraph() for i in range(2)] # lower and upper bounds

        # fit for each q in range
        for pi, qindex in enumerate(selqbins, 0): # pi stands for point-index
            # get the qvalue, define and fill thresholds
            q = qs[qindex]
            if args.verbose: print(f"Processing bin {globalpi}, q = {q}")
            # get pdf values for this q
            pdf = Hq[:, qindex]
            # remove zeros (to avoid inf when doing the log)
            sel = ~(pdf == 0)
            pdf  = pdf[sel]
            mus_ = mus[sel]
            # log(pdf) instead of pdf is fitted
            logpdf = np.log(pdf)

            # save graph for logpdf vs mu to output file
            if len(mus_) != 0: glogpdf = ROOT.TGraph(len(mus_), mus_, logpdf)
            else:              glogpdf = ROOT.TGraph()            
            fout.WriteObject(glogpdf, f"glogPDF_Rang{rang}_{globalpi}")

            # define mu thresholds
            mulow = q - 4.*np.sqrt(q)
            muup  = q + 4.*np.sqrt(q)
            gmuthrs[0].SetPoint(pi, q, mulow)
            gmuthrs[1].SetPoint(pi, q, muup)

            # Fit
            fitfunc = get_fitting_function(mulow, muup)
            if len(mus_) > 0: pars, _ = curve_fit(fitfunc, mus_, logpdf, p0=np.zeros(npar), method="trf")
            else            : pars = np.zeros(npar)

            # save fitted parameters
            for par, gPar in zip(pars, gParams): gPar.SetPoint(pi, q, par)
            globalpi += 1

        # save to output file
        for  pari, gPar   in enumerate(gParams): fout.WriteObject(  gPar, f"gParam_Rang{rang}_{pari}")
        for bound, gmuthr in enumerate(gmuthrs): fout.WriteObject(gmuthr, f"gmuthr_Rang{rang}_{bound}")

    # due to a bug in fiTQun, also save Rang = nqrange with empty values
    gParams = [ROOT.TGraph() for par in range(npars[rang])]
    gmuthrs = [ROOT.TGraph() for i in range(2)]
    rang = nqranges
    for  pari, gPar   in enumerate(gParams): fout.WriteObject(  gPar, f"gParam_Rang{rang}_{pari}")
    for bound, gmuthr in enumerate(gmuthrs): fout.WriteObject(gmuthr, f"gmuthr_Rang{rang}_{bound}")


    # copy PUnhit to output file (copying a TGraph is not implemented in uproot)
    fin = ROOT.TFile(infilename)
    tgraph = fin.Get("PUnhit")
    fout.WriteObject(tgraph, "PUnhit")
    fout.Close()

    # copy 2D charge histogram and PUnhit parameters to output file
    with (uproot.open(infilename) as fin, uproot.update(outfilename) as fout):
        fout["charge2D"]   = fin["charge2D"]
        fout["hPunhitPar"] = fin["hPunhitPar"]

    return 


if __name__ == "__main__":
    main()