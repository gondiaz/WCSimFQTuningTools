import re
import argparse
import uproot
import ROOT
import numpy  as np

from os.path import basename

def histogram2d_to_func(hist, xbins, ybins):
    '''
    Converts a numpy 2D histogram into a function
    Example: todo
    '''
    def getter(x, y):
        xbin = np.digitize(x, xbins) - 1
        ybin = np.digitize(y, ybins) - 1
        return hist[xbin, ybin]
    return getter

def main():

    ############ Program arguments ############
    parser = argparse.ArgumentParser( prog        = f"{basename(__file__)}"
                                    , description = "description"
                                    , epilog      = "Text at the bottom of help")
    
    parser.add_argument("-v", "--verbose", action="store_true")
    parser.add_argument( "-i", "--infile", type=str, nargs="?", help = "", default="cprofiles_merged.root")

    parser.add_argument(    "r0max", type=float, nargs=1, help = "in cm")
    parser.add_argument(  "nr0bins", type=int  , nargs=1, help = "")
    parser.add_argument( "nth0bins", type=int  , nargs=1, help = "")
    parser.add_argument( "--max_relerr", type=float, nargs=1, help = "", default=0.1)
    
    args = parser.parse_args()
    ##########################################

    # collect the arguments
    r0max      = args.r0max   [0]
    nr0bins    = args.nr0bins [0]
    nth0bins   = args.nth0bins[0]
    max_relerr = args.max_relerr

    # define r0, th0 bining
    r0bins  = np.linspace(0, r0max*nr0bins/(nr0bins-1), nr0bins+1)
    th0bins = np.linspace(-1, 1, nth0bins+1)

    # open input file and get momenta
    fin = uproot.open(args.infile)
    momenta = np.sort([float(re.findall(r'\d+(?:\.\d+)?', key)[0]) for key in fin.keys() if re.match("g_\d+", key)])
    # momentum bins such that momenta are the bin centers
    mbins = (momenta[1:] + momenta[:-1])/2.
    mbins = np.insert(mbins,          0, momenta [0] - (mbins [0] - momenta [0]))
    mbins = np.insert(mbins, len(mbins), momenta[-1] - (mbins[-1] - momenta[-1]))

    # define 3D histograms to save the integrals
    for n in range(3):
        locals()[f"th3d_{n}"] = ROOT.TH3D(     f"I_{n}", f"I_{n}"
                                         ,      nr0bins, r0bins
                                         ,     nth0bins, th0bins
                                         , len(momenta), mbins)
        if n == 0: continue
        locals()[f"th1d_{n}"] = ROOT.TH1D( f"I_iso_{n}", f"I_iso_{n}"
                                         , len(momenta), mbins)

    # compute integrals as a function of (momentum, r0, th0)
    if args.verbose: print(f"Performing Cherenkov integrals...")
    for mbin, momentum in enumerate(momenta, 1):
        if args.verbose: print(f"Integrating {momentum} MeV/c...")

        # get cherenkov profile 2D histogram
        th2d = fin[f"g_{momentum}"]
        g, thbins, sbins = th2d.to_numpy()
        g, thbins, sbins = g.copy(), thbins.copy(), sbins.copy()
        thbinw = thbins[1] - thbins[0]
        sbinw  =  sbins[1] -  sbins[0]

        with np.errstate(divide='ignore', invalid='ignore'):
            err = 1./np.sqrt(g)
            g[err>max_relerr] = 0
            if g.sum() == 0:
                if args.verbose: print(f"Warning: Lacking statistics for {momentum} MeV/c...")
                continue
            # normalize the 2D cherenkov profile
            g = histogram2d_to_func(g*(1./g.sum())*(1./(thbinw*sbinw)), thbins, sbins)

        # values of s for evaluation
        s  = (sbins[1:] + sbins[:-1])/2.
        ds = sbins[1] - sbins[0]

        # compute bin centers (avoids errors as th lies alway -1<th<1)
        th0s = (th0bins[1:] + th0bins[:-1])/2.
        r0s  = ( r0bins[1:] +  r0bins[:-1])/2.

        # loop over th0 and r0 values
        for th0bin, th0 in enumerate(th0s, 1):
            for r0bin, r0 in enumerate(r0s, 1):
                # from simple trigonometry
                r2 = r0**2 + s**2 - 2*r0*s*th0
                th = (r0*th0 - s) / np.sqrt(r2)

                # compute integrals
                i = g(th, s)
                I = [sum(i*s**n)*ds for n in range(3)]

                # save computed values
                for n in range(3): locals()[f"th3d_{n}"].SetBinContent(r0bin, th0bin, mbin, I[n])


    # compute isotropic integrals as a function of momentum
    gsthr    = []
    nphotons = []
    momenta_, nevents  = fin[f"g_nevents"].values()
    if args.verbose: print(f"Performing Isotropic integrals...")
    for mbin, momentum in enumerate(momenta, 1):
        if args.verbose: print(f"Integrating {momentum} MeV/c...")

        # get cherenkov profile 2D histogram 
        th2d = fin[f"g_{momentum}"]
        g, thbins, sbins = th2d.to_numpy()
        g, thbins, sbins = g.copy(), thbins.copy(), sbins.copy()
        thbinw = thbins[1] - thbins[0]
        sbinw =   sbins[1] -  sbins[0]

        # mean number of photons
        nphotons.append(g.sum() / nevents[(momentum == momenta_)][0])

        # normalize projected distribution
        g = g*(1./g.sum())*(1./(thbinw*sbinw))

        # values of s for evaluation
        s  = (sbins[1:] + sbins[:-1])/2.
        ds = sbins[1] - sbins[0]

        # get normalized 1D isotropic projection (s)
        iso = np.sum(g, axis=0) * (thbins[1] - thbins[0])

        # compute integrals
        I = [sum(iso*s**n)*ds for n in range(1, 3)]

        # save computed values
        for n in range(1, 3): locals()[f"th1d_{n}"].SetBinContent(mbin, I[n-1])

        # compute gsthr
        gsthr.append(s[np.cumsum(iso)*ds>0.9][0])

    # open output file to save data
    fout = ROOT.TFile("cprofiles_integrals.root", "RECREATE")

    # write 3D integral histograms to output file
    for n in range(3): fout.WriteObject(locals()[f"th3d_{n}"], f"I_{n}")

    # write 1D integral histograms to output file
    for n in range(1, 3): fout.WriteObject(locals()[f"th1d_{n}"], f"I_iso_{n}")

    # write gsthr
    fout.WriteObject(ROOT.TGraph(len(momenta), momenta, np.array(gsthr)), "gsthr")

    # save mean number of photons per momentum
    g = ROOT.TGraph(len(momenta), momenta, np.array(nphotons))
    g.SetTitle("Mean number of photons per event")
    fout.WriteObject(g, "gNphot")

    fout.Close()
    fin .close()
    return


if __name__ == "__main__":
    main()
    