import sys
import glob
import itertools
import argparse
import uproot
import ROOT
import numpy  as np

from os.path   import expandvars, join, basename

from STable.STable_tools import read_wcsim_geometry


def main():

    ############ Program arguments ############
    parser = argparse.ArgumentParser( prog        = f"{basename(__file__)}"
                                    , description = "description"
                                    , epilog      = """""")
    
    parser.add_argument("-v", "--verbose", action="store_true")

    parser.add_argument(      "indir",   type=str,   nargs=1, help = "directory containing produced files")
    parser.add_argument(    "--nbins",   type=int,   nargs=2, help = "histogram bins for the r and eta variables", default=[10, 10])
    parser.add_argument(    "--vaxis",   type=int, nargs="?", help = "detector vertical axis (0)", default=2)
    parser.add_argument(    "--zedge", type=float, nargs="?", help = "histogram edge for vertical (z) dimensions in cm", default=136.95-20.) # default wcte
    parser.add_argument(    "--redge", type=float, nargs="?", help = "histogram edge for radial (r) dimensions in cm"  , default=172.05-20.) # default wcte
    parser.add_argument( "--wcsimlib",   type=str, nargs="?", help = "WCSim lib path"         , default="$HOME/Software/WCSim/install/lib")
    parser.add_argument("--maxnfiles", type=int, nargs="?", help = "max files to proccess", default=None)
    
    args = parser.parse_args()
    ##########################################

    # Load WCSimRoot.so
    ROOT.gSystem.AddDynamicPath(expandvars(args.wcsimlib))
    ROOT.gSystem.Load          ("libWCSimRoot.dylib" if sys.platform == "darwin" else "libWCSimRoot.so")

    # get simulation filenames, removing those ending in '_flat.root'
    infiles = glob.glob(join(args.indir[0], "*"))
    infiles = sorted([f for f in infiles if "_flat.root" not in basename(f)])
    if args.maxnfiles: infiles = infiles[:args.maxnfiles]

    # Select rotation matrix based on vertical axis
    vaxis = args.vaxis
    if   vaxis == 0:
        R = np.array([[0, 0, -1], [0, 1, 0], [1, 0, 0]])
    elif vaxis == 1:
        R = np.array([[1, 0, 0], [0, 0, -1], [0, 1, 0]])
    elif vaxis == 2: pass
    else: raise Exception("Invalid value for vertical axis (vaxis)")

    # read PMT information (ids, positions, orientations)
    _, pmts_df = read_wcsim_geometry(infiles[0])
    pmts_df = pmts_df.set_index("TubeNo")

    # # get ids of the tubes normal to the vessel (assumed to be the central PMTs)
    # inormal = pmts_df[pmts_df.mPMT_PMTNo == 1].index.values.astype(int)

    # these lines are commented because these parameter definitions ere misleading in WCSim
    # tube_ztop = df.loc["WCCylLength", "WC"]/2.  # detector cylinder half-length
    # tube_rad  = df.loc["WCPMTRadius", "WC"]     # PMT module radius
    # cyl_rad   = df.loc["WCCylRadius", "WC"]     # detector cylinder radius

    # define histogram bins
    detLength = args.zedge*2.
    detRad    = args.redge
    rbins     = np.linspace(0, min(detRad, detLength/2.), args.nbins[0]+1)
    etabins   = np.linspace(0,                         1, args.nbins[1]+1)

    # empty histograms to be filled
    Hside = np.zeros((len(rbins)-1, len(etabins)-1))
    Hcaps = np.zeros((len(rbins)-1, len(etabins)-1))
    Hall  = np.zeros((len(rbins)-1, len(etabins)-1))
    # Loop over simulation files and fill the histograms
    for n, filename in enumerate(infiles, 1):
        if args.verbose: 
            print(f"Processing file {n}/{len(infiles)}...", basename(filename))
        # read data
        with uproot.open(filename) as f:
            # contains the following data if the optical photon reaches a PMT
            # either directly or indirectly
            # photon emission position
            oppos   = f["sttree/oppos"]  .array().to_numpy()
            # PMT number/pos
            ihPMT   = f["sttree/ihPMT"]  .array().to_numpy()
            tubepos = f["sttree/tubepos"].array().to_numpy()
            # ISTORY=(# of refl)*1000 + (# of scat); 0 if direct hit
            isct = f["sttree/isct"].array().to_numpy()

        # select direct photons
        sel_direct = (isct == 0)
        oppos   = oppos  [sel_direct]
        ihPMT   = ihPMT  [sel_direct]
        tubepos = tubepos[sel_direct]

        # select only PMTs normal to the vessel (central PMTs)
        sel_normal = np.isin(ihPMT, inormal)
        oppos   = oppos  [sel_normal]
        ihPMT   = ihPMT  [sel_normal]
        tubepos = tubepos[sel_normal]

        # get PMT orientation
        tubedir = pmts_df.loc[ihPMT, ["Orientation_x0", "Orientation_x1", "Orientation_x2"]].values.astype(float)

        # transform from mm to cm
        oppos   *= 0.1
        tubepos *= 0.1

        # compute variables of interest:
        #    - distance between PMT and emission point
        #    - cosine of the angle of incidence w.r.t. the PMT normal
        # don't need to rotate since module and angles are invariant
        # compute distance
        r    = oppos - tubepos
        modr = np.sqrt(np.sum(r**2, axis=1))
        # select only photons within maximum possible distance
        sel_r = modr < rbins[-1]  # the < instead <= is used to avoid out of bound at np.digitize below
        r       =       r[sel_r]
        modr    =    modr[sel_r]
        tubepos = tubepos[sel_r]
        tubedir = tubedir[sel_r]
        # compute angles
        eta  = np.array([np.dot(v1, v2) for (v1, v2) in zip(r, tubedir)])/modr

        # rotate tube positions to perform fiducial selections based on vertical and radial directions below
        if vaxis != 2:
            tubepos = np.matmul(R, tubepos.T).T

        # select photons in valid shells
        rlim = rbins[np.digitize(modr, rbins)]
        sel_fiducial_shell_side = (rlim + abs(tubepos[:, 2]))                             <= detLength/2.
        sel_fiducial_shell_caps = (rlim + np.sqrt(np.sum(tubepos[:, (0, 1)]**2, axis=1))) <= detRad

        # 2D histogram in distance and angle variables
        hside, _, _ = np.histogram2d(modr[sel_fiducial_shell_side], eta[sel_fiducial_shell_side], bins=[rbins, etabins])
        hcaps, _, _ = np.histogram2d(modr[sel_fiducial_shell_caps], eta[sel_fiducial_shell_caps], bins=[rbins, etabins])
        hall = hside + hcaps

        # add histograms
        Hside += hside
        Hcaps += hcaps
        Hall  += hall

    # normalize histogram as follow: for each value of r, normalize the angular distribution
    # w.r.t normal incidence, i.e. the entry corresponding to cos(eta) = 1
    np.seterr(invalid="ignore") # ignore invalid value when dividing by zero
    Hside = Hside/Hside[:, -1][:, np.newaxis]
    Hcaps = Hcaps/Hcaps[:, -1][:, np.newaxis]
    Hall  = Hall /Hall [:, -1][:, np.newaxis]

    # save photon tables
    if args.verbose:
        print("Writting histograms...")
    
    # define root 2D histograms
    th2d_side = ROOT.TH2D("AngularResponse_side", "AngularResponse_side", len(rbins)-1, rbins, len(etabins)-1, etabins)
    th2d_caps = ROOT.TH2D("AngularResponse_caps", "AngularResponse_caps", len(rbins)-1, rbins, len(etabins)-1, etabins)
    th2d_all  = ROOT.TH2D("AngularResponse_all" , "AngularResponse_all" , len(rbins)-1, rbins, len(etabins)-1, etabins)

    # fill histograms
    for ix, iy in itertools.product(range(1, len(rbins)), range(1, len(etabins))):
        th2d_side.SetBinContent(ix, iy, Hside[ix-1, iy-1])
        th2d_caps.SetBinContent(ix, iy, Hcaps[ix-1, iy-1])
        th2d_all .SetBinContent(ix, iy, Hall [ix-1, iy-1])

    # write output
    fout = ROOT.TFile("angular.root", "RECREATE")
    fout.WriteObject(th2d_side, "AngularResponse_side")
    fout.WriteObject(th2d_caps, "AngularResponse_caps")
    fout.WriteObject(th2d_all , "AngularResponse_all")
    fout.Close()

    return


if __name__ == "__main__":
    main()
    