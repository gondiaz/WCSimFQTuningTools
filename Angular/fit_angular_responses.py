import sys
import argparse
import uproot
import ROOT
import numpy  as np
import matplotlib.pyplot as plt

from os.path import basename

from scipy.optimize import curve_fit


def polynomial(x, *coefficients):
    return np.polyval(coefficients, x)


def main():

    ############ Program arguments ############
    parser = argparse.ArgumentParser( prog        = f"{basename(__file__)}"
                                    , description = "description"
                                    , epilog      = """""")
    
    parser.add_argument("-v", "--verbose", action="store_true")
    parser.add_argument("-i",  "--infile", type=str, nargs="?", help = "", default="angular.root")
    parser.add_argument("-l",    "--lims", type=float, nargs=2, help="")

    parser.add_argument("polydeg", type=int, nargs=1, help = "")

    args = parser.parse_args()
    ##########################################

    lims    = args.lims
    polydeg = args.polydeg[0]

    # read angular responses file
    with uproot.open(args.infile) as infile:
        # read angular response histogram computed from all the PMTs (side and caps)
        h, rbins, etabins = infile["AngularResponse_all"].to_numpy()
        # select 
        il = np.digitize(lims[0], rbins)
        iu = np.digitize(lims[1], rbins)-1
        e = np.mean(h[il:iu, :], axis=0)

    # compute bin centers and fit angular response
    etas = (etabins[1:] + etabins[:-1])/2.
    pars, _ = curve_fit(polynomial, etas, e, p0=np.ones(polydeg + 1)/2., bounds=(0, np.ones(polydeg + 1)))

    # create polynomial TF1 and set the computed parameters
    tf1 = ROOT.TF1("angResp", f"pol{polydeg}", etabins[0], etabins[-1])
    tf1.SetParameters(*np.flip(pars))

    if args.verbose:
        # plot results
        plt.figure(figsize=[5, 3])
        plt.scatter(etas, e, color="k")
        x = np.linspace(etabins[0], etabins[-1], 100)
        fitted = [tf1.Eval(eta) for eta in x]
        plt.plot(x, fitted, color="r")
        plt.xlabel(r"$\cos \eta$")
        plt.ylabel(r"$\epsilon(\cos \eta)$")
        plt.tight_layout()
        plt.show()

    # write output
    fout = ROOT.TFile("fitted_angular.root", "RECREATE")
    fout.WriteObject(tf1, "angResp")
    fout.Close()

    return


if __name__ == "__main__":
    main()