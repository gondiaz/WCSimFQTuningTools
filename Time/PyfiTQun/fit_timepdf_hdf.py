import re
import argparse
import warnings
import uproot
import ROOT
import numpy  as np
import tables as tb
from scipy.optimize import curve_fit

from os.path import basename

warnings.filterwarnings("ignore", category=tb.NaturalNameWarning)

# define gaussian pdf
def gaussian(x, mu, sig):
    return (1/(sig*np.sqrt(2.*np.pi)))*np.exp(-(x-mu)**2/(2*sig**2))

# define indirect light timepdf
def indirect_timepdf(t, delta, sig, gamma):
    out = np.zeros(len(t))
    dt = t - delta
    A = 1./(np.sqrt(np.pi/2.)*sig + 2*gamma)
    sel = dt<0
    out[sel]  = A*np.exp(-dt[sel]**2/(2.*sig**2))
    out[~sel] = A*(dt[~sel]/gamma + 1.)*np.exp(-dt[~sel]/gamma)
    return out


def main():

    ############ Program arguments ############
    parser = argparse.ArgumentParser( prog        = f"{basename(__file__)}"
                                    , description = ""
                                    , epilog      = """""")
    
    parser.add_argument("-v", "--verbose", action="store_true")
    parser.add_argument(  "--infile",   type=str, nargs="?", help = "time 2D file", default="tpdf_histograms.h5")

    parser.add_argument( "--min_entries", type=int, nargs="?", help = "min number of entries", default=100)
    parser.add_argument("npars_gauss", type=int, nargs=1, help = "gaussian fits")
    parser.add_argument("npars_pars" , type=int, nargs=1, help = "parameter fits")
    
    args = parser.parse_args()
    ##########################################

    # Config variables
    infilename  = args.infile
    min_entries = args.min_entries
    npars_gauss = args.npars_gauss[0]
    npars_pars  = args.npars_pars [0]

    # read momentum values and define momentum bins
    with tb.open_file(infilename) as f:
        momenta_keys = list(f.root.direct._v_children)
        for key in list(f.root.indirect._v_children): assert key in momenta_keys 
        momenta = np.sort([float(re.findall(r"p_(\d+\.\d+)", key)[0]) for key in momenta_keys])

    pbins = (momenta[:-1] + momenta[1:])/2.
    pbins = np.insert(pbins,            0, momenta[0]  - (pbins[0]  - momenta[0]))
    pbins = np.insert(pbins, len(momenta), momenta[-1] - (pbins[-1] - momenta[-1]))

    # read 2D histogram bins (tres vs μ)
    with tb.open_file(infilename) as f:
        tresbins = f.root.bins.tres.read()
        μbins    = f.root.bins.μ   .read()
    tress = (tresbins[1:] + tresbins[:-1])/2.
    μs    = (μbins   [1:] + μbins   [:-1])/2.

    # define outfile
    outfilename = "fitted_timepdf.h5"
    with tb.open_file(outfilename, "w") as f:
        f.create_group("/", "gaussian_fits")
        f.create_group("/", "polynomial_fits")
        f.create_group("/", "bins")
        f.create_group("/", "parameter_fits")
        f.create_array("/bins", "momentum", pbins)
        f.create_array("/bins",     "tres", tresbins)
        f.create_array("/bins",        "μ", μbins)

    # define fit parameters 2D histograms to fill
    hmeans  = np.zeros((len(pbins)-1, len(μbins)-1))
    hsigmas = np.zeros((len(pbins)-1, len(μbins)-1))

    ################# Direct light ##############
    # First loop in momenta: fill Gaussian fit parameters for each (p, μ)
    if args.verbose: print("1) Performing gaussian fits...")
    for pi, p in enumerate(momenta, 0):
        # read 2D histogram for this momentum
        with tb.open_file(infilename) as f:
            H = getattr(f.root.direct, f"p_{p}").read()

        # loop in true charge
        for μi, μ in enumerate(μs, 0):
            # get projection on residual time direction
            proj = H[:, μi]
            norm = np.sum(proj * (tresbins[1:] - tresbins[:-1]))

            # do the fit (only if number of entries is enough)
            if (proj.sum()<min_entries): pars = [np.nan, np.nan]
            else:
                try:
                    pars, cov = curve_fit(gaussian, tress, proj/norm, p0=[0, 0.1], method="trf")
                    if ((pars[0]<tresbins[0])&(tresbins[-1]<pars[0])): pars = [np.nan, np.nan]
                except RuntimeError: pars = [np.nan, np.nan] # if RuntimeError: minimization fails
                
            # save result
            hmeans [pi, μi] = pars[0]
            hsigmas[pi, μi] = pars[1]

    # save gaussian fit results
    with tb.open_file(outfilename, "a") as f:
        f.create_array("/gaussian_fits",  "means", hmeans)
        f.create_array("/gaussian_fits", "sigmas", hsigmas)
    if args.verbose: print("==> 1) gaussian fits finished")

    # Second loop in momenta: Polynomial fits, for each p
    if args.verbose: print("2) Performing polynomial fits for each momentum...")

    hgauss = [hmeans, hsigmas]
    hpars  = [np.zeros((len(pbins)-1, npars_gauss)) for i in range(2)]
    for i, h in enumerate(hgauss):
        for pi, p in enumerate(momenta):
            # get parameter (mean or sigma)
            gausspar = h[pi, :]
            # polynomial fit
            sel = ~np.isnan(gausspar)
            pars = np.flip(np.polyfit(μs[sel], gausspar[sel], npars_gauss-1))
            # save values
            hpars[i][pi] = pars

    # save parameters to output file
    with tb.open_file(outfilename, "a") as f:
        f.create_array("/polynomial_fits",  "means", hpars[0])
        f.create_array("/polynomial_fits", "sigmas", hpars[1])
    if args.verbose:
        print("==> 2) polynomial fits for each momentum finished")


    # Polynomial fits p vs fit parameters
    if args.verbose:
        print("3) Performing polynomial fits for each parameter...")

    for i, parname in enumerate(["means", "sigmas"]):
        parameters = np.zeros((npars_gauss, npars_pars))
        for pari in range(npars_gauss):
            npars = npars_pars + 1
            while npars>2:
                npars -= 1
                if npars < 2:
                    raise Exception("Cannot perform 3) polynomial fits for each momentum")
                try:
                    pars = np.flip(np.polyfit(momenta, hpars[i][:, pari], npars - 1))
                    break
                except np.RankWarning:
                    continue
            parameters[pari, :npars] = pars

        # save parameters to output file
        with tb.open_file(outfilename, "a") as f:
            f.create_array("/parameter_fits", parname, parameters)

    if args.verbose:
        print("==> 3) polynomial fits for each parameter finished")


    # save input file direct light 2D histograms as 3D histogram
    # TODO: normalize
    Hdirect   = np.zeros((len(pbins)-1, len(tresbins)-1, len(μbins)-1))
    Hindirect = np.zeros((len(pbins)-1, len(tresbins)-1, len(μbins)-1))
    for pi, p in enumerate(momenta, 0):
        # read 2D histogram for this momentum
        with tb.open_file(infilename) as f:
            H = getattr(f.root.direct  , f"p_{p}").read()
            Hdirect[pi]   = H
            H = getattr(f.root.indirect, f"p_{p}").read()
            Hindirect[pi] = H
    with tb.open_file(outfilename, "a") as f:
        f.create_array("/",   "direct_timepdfs", Hdirect)
        f.create_array("/", "indirect_timepdfs", Hindirect)
    ################# Direct light ##############


    ################# Indirect light ##############
    # Fit timepdf for indirect light
    # get normalized distributions for each momentum
    timepdfs = np.zeros((len(momenta_keys), len(tresbins)-1))
    for pi, p_key in enumerate(momenta_keys, 0):
        # read 2D histogram for this momentum
        with tb.open_file(infilename) as f:
            H = getattr(f.root.indirect, p_key).read()
            
        # project and normalize
        h    = H.sum(axis=1)
        norm = np.sum(h*(tresbins[1:]-tresbins[:-1]))
        timepdfs[pi] = h/norm

    # merge momentum distributions and fit
    pdf = np.zeros(len(tresbins)-1)
    for tpdf in timepdfs: pdf += tpdf
    norm = np.sum(pdf*(tresbins[1:]-tresbins[:-1]))

    try:
        pdf  = pdf/norm
        pars, cov = curve_fit(indirect_timepdf, tress, pdf)
    except:
        warnings.warn("Indirect pdf fit didn't succed, setting parameters to zero value")
        pars = np.zeros(3)
    
    # save indirect pdf parameters
    with tb.open_file(outfilename, "a") as f:
        f.create_array("/", "indirect_fit", pars)

    return

if __name__ == "__main__":
    main()