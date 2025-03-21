import argparse
import uproot
import ROOT
import numpy as np

from os.path import basename


def histogram2d_to_func(hist, xbins, ybins):
    """
    Converts a numpy 2D histogram into a function
    Example: todo
    """
    def getter(x, y):
        xbin = np.digitize(x, xbins) - 1
        ybin = np.digitize(y, ybins) - 1
        return hist[xbin, ybin]
    return getter


def split_array(array, N, nmin=1):

    """Given an array, it split it in N partitions
    with nmin minimun number of elements.
    It returns L%N partitions with (L//N + 1) elements
    and N-L%N with L//N.
    
    The last element of each partition coincides with
    the first element of the next.

    For example:
    >>> array = [1, 2, 3, 4, 5]
    >>> split_array(array, 2)
    [[1, 2, 3, 4], [4, 5]]

    Notice that in the result the first partition contains 
    2 more values because of the first element of the fist partiton
    is also added.
    """
    splitted = []
    n = len(array)//N
    r = len(array)% N
    if n<nmin:
        n = nmin
        N = int(np.ceil(len(array)/n))
        r = len(array)% N
    iu = 0
    while (iu<len(array)-1):
        il = iu
        iu = il + n
        if r: iu += 1
        splitted.append(array[il: iu+1])
        r-=1
    if len(splitted[-1]) < nmin:
        splitted[-2] = np.concatenate((splitted[-2], splitted[-1][1:]))
        splitted.pop()
    return splitted


def compute_goodness_of_fit(xs, ys, parss):
    """
    This function computes the coefficient of determination
    https://en.wikipedia.org/wiki/Coefficient_of_determination
    for a polynomial fit.

    It computes the coefficient for a sectioned polynolial fit
    xs: list with the x values for each partition
    ys: list with the y values for each partition
    parss: polynomial parameters for each partition
    """
    Sres = 0
    Stot = 0
    n    = 0
    ysum = 0
    for s, (x, y, pars) in enumerate(zip(xs, ys, parss)):
        pol = np.poly1d(pars)
        Sres += np.sum((y - pol(x))**2)
        ysum += np.sum(y)
        n    += len(x)
    ymean = ysum/n
    for y in ys: Stot += np.sum((y - ymean)**2)
    if   (Stot == 0.): return 1. # i.e avoid 0/0
    else             : return 1. - Sres/Stot


def main():

    ############ Program arguments ############
    parser = argparse.ArgumentParser( prog        = f"{basename(__file__)}"
                                    , description = "description"
                                    , epilog      = "Text at the bottom of help")
    
    parser.add_argument("-v", "--verbose", action="store_true")
    parser.add_argument( "-i", "--infile", type=str, nargs="?", help = "", default="cprofiles_integrals.root")

    parser.add_argument(      "npars", type=int,   nargs=1, help = "")
    parser.add_argument( "--nsectmax", type=int, nargs="?", help = "", default=10)
    parser.add_argument(     "--Rmin", type=float, nargs="?", help = "", default=0.99)
    
    args = parser.parse_args()
    ##########################################

    npars    = args.npars[0]
    nsectmax = args.nsectmax
    Rmin     = args.Rmin    

    # open input and output files
    fin  = uproot.open(args.infile)
    fout = ROOT.TFile("cprofiles_fits.root", "RECREATE")

    # fit parameters are saved in 3D histo, thus define the bins
    # add bins to save the fit bounds (as required by fiTQun)
    nparbins = npars*nsectmax + nsectmax + 1
    parsbins = np.arange(-0.5, nparbins+0.5, 1)

    # loop over s exponent in the integrals
    for n in range(3):

        if args.verbose: print(f"Fitting n={n}...")

        h, r0bins, th0bins, mbins = fin[f"I_{n}"].to_numpy()

        th3d = ROOT.TH3D( f"Js_{n}", f"Js_{n}"
                        ,  len(r0bins)-1,  r0bins
                        , len(th0bins)-1, th0bins
                        ,       nparbins, parsbins)
        
        th2d_chi2 = ROOT.TH2D( f"chi2_{n}", f"chi2_{n}"
                             ,  len(r0bins)-1,  r0bins
                             , len(th0bins)-1, th0bins)
        
        th2d_nsect = ROOT.TH2D( f"nsect_{n}", f"nsect_{n}"
                              ,  len(r0bins)-1,  r0bins
                              , len(th0bins)-1, th0bins)

        # fit I_n vs log(p)
        log_momenta = np.log((mbins[1:] + mbins[:-1])/2.)

        # select momenta (only those with values of I, i.e. different from 0)
        sel_momenta = (h.sum(axis=(0, 1)) > 0.)
        log_momenta = log_momenta[sel_momenta]
        h           = h[:,:, sel_momenta]

        # loop over r0 and th0 bins
        for th0bin, _ in enumerate(th0bins[:-1], 0):
            for r0bin, _ in enumerate(r0bins[:-1], 0):
                # get integrals for (r0, th0) bins
                I = h[r0bin, th0bin]
                
                # polynomial fit, full region
                js = np.polyfit(log_momenta, I, npars-1)
                J = [(0, log_momenta[0], log_momenta[-1], js)]
                R = compute_goodness_of_fit([log_momenta], [I], [js])

                # perform sectioned fits, if needed
                nsect = 1
                # stop if number of sections saturates if number of points lower than parameters
                stop  = False
                while ((R<Rmin) & (nsect<nsectmax) & (~stop)):
                    nsect += 1
                    J = []
                    # create sections
                    sect_momenta = split_array(log_momenta, nsect, npars)
                    sect_I       = split_array(          I, nsect, npars)
                    # saturation of number of sections
                    if len(sect_momenta)<nsect:
                        nsect = len(sect_momenta)
                        stop  = True
                    # polynomial fit for each section
                    for s, (x, y) in enumerate(zip(sect_momenta, sect_I)):
                        j = np.polyfit(x, y, npars-1)
                        J.append((s, x[0], x[-1], j))
                    # compute coefficient of determination
                    R = compute_goodness_of_fit(sect_momenta, sect_I, [section[-1] for section in J])

                # fill 3D histogram with parameters and bounds
                for section in J:
                    s  = section[0]
                    js = section[3]
                    for i, ji in enumerate(np.flip(js), s*npars+1): th3d.SetBinContent(r0bin+1, th0bin+1, i, ji)
                    th3d.SetBinContent(r0bin+1, th0bin+1, nsect*npars+s+1, section[1])
                th3d.SetBinContent(r0bin+1, th0bin+1, nsect*(npars+1)+1, section[2])
                for i in range((npars+1)*nsect+2, nparbins+1): th3d.SetBinContent(r0bin+1, th0bin+1, i, -9999.)

                # save 1 in nsect histo (to be fiTQun compatible)
                th2d_nsect.SetBinContent(r0bin+1, th0bin+1, nsect)

                # TODO: compute chi2
                th2d_chi2.SetBinContent(r0bin+1, th0bin+1, 1)

        fout.WriteObject(      th3d, f"hI3d_par_{n}")
        fout.WriteObject( th2d_chi2, f"hI3d_par_{n}_chi2")
        fout.WriteObject(th2d_nsect, f"hI3d_nsect_{n}")

    # Add hprofinf
    hprofinf = ROOT.TH1F("hprofinf", "hprofinf", 6, np.arange(-0.5, 6.5, 0.5))
    hprofinf.SetBinContent(1, npars)
    hprofinf.SetBinContent(2, nsectmax)          # not really used in fiTQun
    hprofinf.SetBinContent(3, 0)                 # momentum offset
    hprofinf.SetBinContent(4, np.min(np.diff(np.log(mbins)))) # momentum step (?)
    hprofinf.SetBinContent(5, np.log(mbins[0]))  # min momentum
    hprofinf.SetBinContent(6, np.log(mbins[-1])) # max momentum
    fout.WriteObject(hprofinf, "hprofinf")
    
    # fit and save isotropic integrals, gNphot and gsthr
    # add +2 bins to save the fit bounds (to be fiTQun compatible)
    parsbins = np.arange(-0.5, npars+0.5+2, 1)
    # fit and save isotopric integrals
    for n in range(1, 3):
        h, mbins = fin[f"I_iso_{n}"].to_numpy()
        momenta = (mbins[1:] + mbins[:-1])/2.
        pars = np.polyfit(np.log(momenta), np.log(h), npars-1)
        th1d = ROOT.TH1D(f"gI_iso_{n}_pars", f"gI_iso_{n}_pars",  len(parsbins)-1, parsbins)
        for i, pi in enumerate(np.flip(pars), 1): th1d.SetBinContent(i, pi)
        th1d.SetBinContent(npars+1, np.log(momenta[0]))
        th1d.SetBinContent(npars+2, np.log(momenta[-1]))
        fout.WriteObject(th1d, f"gI_iso_{n}_pars")

    # copy and fit mean number of photons per momentum (gNphot)
    x, y = fin["gNphot"].values()
    g = ROOT.TGraph(len(x), x, y)
    g.SetTitle("Mean number of photons per event")
    fout.WriteObject(g, "gNphot")
    # fit, fill histogram and save
    pars = np.polyfit(np.log(x), np.log(y), npars-1)
    th1d = ROOT.TH1D("gNphot_pars", "gNphot_pars",  len(parsbins)-1, parsbins)
    for i, pi in enumerate(np.flip(pars), 1): th1d.SetBinContent(i, pi)
    th1d.SetBinContent(npars+1, np.log(x[0]))
    th1d.SetBinContent(npars+2, np.log(x[-1]))
    fout.WriteObject(th1d, "gNphot_pars")

    # copy and fit gsthr
    x, y = fin["gsthr"].values()
    g = ROOT.TGraph(len(x), x, y)
    fout.WriteObject(g, "gsthr")
    # fit, fill histogram and save
    pars = np.polyfit(np.log(x), np.log(y), npars-1)
    th1d = ROOT.TH1D("gsthr_pars", "gsthr_pars",  len(parsbins)-1, parsbins)
    for i, pi in enumerate(np.flip(pars), 1): th1d.SetBinContent(i, pi)
    th1d.SetBinContent(npars+1, np.log(x[0]))
    th1d.SetBinContent(npars+2, np.log(x[-1]))
    fout.WriteObject(th1d, "gsthr_pars")
    
    # copy the original integral histograms to allow easier checking
    for n in range(0, 3): fout.WriteObject(fin[f"I_{n}"]    .to_pyroot(), f"I_{n}")
    for n in range(1, 3): fout.WriteObject(fin[f"I_iso_{n}"].to_pyroot(), f"hI_iso_{n}")
    
    fout.Close()
    fin .close()
    return

if __name__ == "__main__":
    main()