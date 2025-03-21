import sys
import glob
import argparse
import uproot
import ROOT
import numpy as np

from os.path import expandvars, join, basename

from STable.STable_tools import split_tubeids, clockwise_azimuth_angle


def main():

    ############ Program arguments ############
    parser = argparse.ArgumentParser( prog        = f"{basename(__file__)}"
                                    , description = "description"
                                    , epilog      = """ The table variables are the following:
                                                    zs: source z position
                                                    Rs: source R position
                                                    (bottom and top table) R_PMT: R position of the PMT
                                                    (side table) z_PMT: zposition of the PMT
                                                    phi: angle between source and PMT
                                                    zd : source z direction
                                                    theta: source direction angle w.r.t PMT""")
    
    parser.add_argument("-v", "--verbose", action="store_true")
    parser.add_argument("indir", type=str, nargs=1, help = "directory containing produced files")
    parser.add_argument( "--wcsimlib", type=str, nargs=1, help="WCSim lib path")
    parser.add_argument("--fitqunlib", type=str, nargs=1, help="fiTQun lib path")
    parser.add_argument("--vaxis", type=int, nargs="?", help = "detector vertical axis (0)", default=2)
    parser.add_argument("--nbins", type=int, nargs=7, help = "histogram bins in the following order: zs,Rs,R_PMT,z_PMT,phi,zd,theta")
    parser.add_argument("--zedge", type=float, nargs=1, help = "histogram edge for vertical (z) dimensions in cm")
    parser.add_argument("--redge", type=float, nargs=1, help = "histogram edge for radial (r) dimensions in cm")
    
    args = parser.parse_args()
    ##########################################

    # Load WCSimRoot.so
    ROOT.gSystem.AddDynamicPath(expandvars(args.wcsimlib[0]))
    ROOT.gSystem.Load          ("libWCSimRoot.dylib" if sys.platform == "darwin" else "libWCSimRoot.so")

    # Load libfiTQunLib.so to access TScatTable class
    ROOT.gSystem.AddDynamicPath(expandvars(args.fitqunlib[0]))
    ROOT.gSystem.Load          ("libfiTQunLib.dylib" if sys.platform == "darwin" else "libfiTQunLib.so")

    # get simulation filenames, removing those ending in '_flat.root'
    infiles = glob.glob(join(args.indir[0], "*"))
    infiles = sorted([f for f in infiles if "_flat.root" not in basename(f)])

    # Select rotation matrix based on vertical axis
    vaxis = args.vaxis
    if   vaxis == 0:
        R = np.array([[0, 0, -1], [0, 1, 0], [1, 0, 0]])
    elif vaxis == 1:
        R = np.array([[1, 0, 0], [0, 0, -1], [0, 1, 0]])
    elif vaxis == 2: pass
    else: raise Exception("Invalid value for vertical axis (vaxis)")

    # get bottom, top and side tube-ids from the first file
    tubeid_bottom, tubeid_top, tubeid_side = split_tubeids(infiles[0], vaxis=vaxis)

    # Define scatering table binning
    # 6 variables
    # zs: source z position
    # Rs: source R position
    # (bottom and top table) R_PMT
    # (side table) z_PMT: zposition of the PMT
    # phi: angle between source and PMT
    # zd : source z direction
    # theta: source direction angle w.r.t PMT
    nzs    = args.nbins[0]
    nRs    = args.nbins[1]
    nR_PMT = args.nbins[2]
    nz_PMT = args.nbins[3]
    nphi   = args.nbins[4]
    nzd    = args.nbins[5]
    ntheta = args.nbins[6]

    zedge = args.zedge[0]
    redge = args.redge[0]

    zs    = [-zedge, zedge]
    Rs    = [0., redge]
    R_PMT = [0., redge]
    z_PMT = [-zedge, zedge]
    phi   = [-np.pi, +np.pi]
    zd    = [-1., 1.]
    theta = [-np.pi, +np.pi]

    # Define scattering tables, notice the only change beside name and title is the use of R_PMT (bottom and top) and z_PMT (side)
    # bottom
    Tbottom_indirect = ROOT.TScatTableF( "bottomnphot", "Number of Photons for Tubes on the Bottom of the Tank"
                                       , nzs, zs[0], zs[1], nRs, Rs[0], Rs[1], nR_PMT, R_PMT[0], R_PMT[1]
                                       , nphi, phi[0], phi[1], nzd, zd[0], zd[1], ntheta, theta[0], theta[1])
    
    Tbottom_direct   = ROOT.TScatTableF( "bottomnphot", "Number of Photons for Tubes on the Bottom of the Tank"
                                       , nzs, zs[0], zs[1], nRs, Rs[0], Rs[1], nR_PMT, R_PMT[0], R_PMT[1], nphi, phi[0], phi[1])
    
    # top
    Ttop_indirect    = ROOT.TScatTableF( "topnphot", "Number of Photons for Tubes on the Top of the Tank"
                                       , nzs, zs[0], zs[1], nRs, Rs[0], Rs[1], nR_PMT, R_PMT[0], R_PMT[1]
                                       , nphi, phi[0], phi[1], nzd, zd[0], zd[1], ntheta, theta[0], theta[1])
    
    Ttop_direct      = ROOT.TScatTableF( "topnphot", "Number of Photons for Tubes on the Top of the Tank"
                                       , nzs, zs[0], zs[1], nRs, Rs[0], Rs[1], nR_PMT, R_PMT[0], R_PMT[1], nphi, phi[0], phi[1])

    # side
    Tside_indirect = ROOT.TScatTableF( "sidenphot", "Number of Photons for Tubes on the Side of the Tank"
                                     , nzs, zs[0], zs[1], nRs, Rs[0], Rs[1], nz_PMT, z_PMT[0], z_PMT[1]
                                     , nphi, phi[0], phi[1], nzd, zd[0], zd[1], ntheta, theta[0], theta[1])
    
    Tside_direct   = ROOT.TScatTableF( "sidenphot", "Number of Photons for Tubes on the Side of the Tank"
                                     , nzs, zs[0], zs[1], nRs, Rs[0], Rs[1], nz_PMT, z_PMT[0], z_PMT[1], nphi, phi[0], phi[1])

    # Loop over simulation files in order to fill the table histograms
    for n, filename in enumerate(infiles, 1):
        if args.verbose: print(f"Processing file {n}/{len(infiles)}...", basename(filename))
        # read data
        with uproot.open(filename) as f:
            try:
                # contains the following data if the optical photon reaches a PMT
                # either directly or indirectly
                # parent pos and dir
                srcpos = f["sttree/srcpos"].array().to_numpy()
                srcdir = f["sttree/srcdir"].array().to_numpy()
                # photon pos and dir
                #oppos  = f["sttree/oppos"] .array().to_numpy()
                #opdir  = f["sttree/opdir"] .array().to_numpy()
                # PMT number/pos
                ihPMT   = f["sttree/ihPMT"]  .array().to_numpy()
                tubepos = f["sttree/tubepos"].array().to_numpy()
                # distance covered? (not needed)
                deltaL = f["sttree/deltaL"].array().to_numpy()
                # ISTORY=(# of refl)*1000 + (# of scat); 0 if direct hit
                isct = f["sttree/isct"].array().to_numpy()
            except uproot.exceptions.KeyInFileError: continue


        # Rotate with vertical direction given by vaxis
        if vaxis != 2:
            srcpos  = np.matmul(R,  srcpos.T).T
            srcdir  = np.matmul(R,  srcdir.T).T
            #oppos   = np.matmul(R,   oppos.T).T
            #opdir   = np.matmul(R,   opdir.T).T
            tubepos = np.matmul(R, tubepos.T).T

        # Compute scattering table variables
        zs    =  srcpos[:, 2]
        Rs    = np.sum(srcpos[:, (0, 1)]**2, axis=1)**0.5
        # z_PMT used for side table, R_PMT used for bottom and top
        z_PMT = tubepos[:, 2]
        R_PMT = np.sum(tubepos[:, (0, 1)]**2, axis=1)**0.5
        phi   = clockwise_azimuth_angle(tubepos, srcpos)
        zd    = srcdir[:, 2]
        theta = clockwise_azimuth_angle(srcdir, tubepos-srcpos)

        # transform from mm to cm
        zs    *= 0.1
        Rs    *= 0.1
        z_PMT *= 0.1
        R_PMT *= 0.1

        # transform angles from [0, 2pi] to [-pi, pi]
        phi  [  phi>np.pi] -= 2.*np.pi
        theta[theta>np.pi] -= 2.*np.pi

        # Fill indirect tables
        sel_indirect = (isct != 0) # select indirect photons
        # bottom
        sel = sel_indirect & np.isin(ihPMT, tubeid_bottom)
        for vals in zip(zs[ sel], Rs[ sel], R_PMT[ sel], phi[ sel], zd[ sel], theta[ sel]):Tbottom_indirect.Fill(*vals)
        # top
        sel = sel_indirect & np.isin(ihPMT, tubeid_top)
        for vals in zip(zs[ sel], Rs[ sel], R_PMT[ sel], phi[ sel], zd[ sel], theta[ sel]):Ttop_indirect   .Fill(*vals)
        # side
        sel = sel_indirect & np.isin(ihPMT, tubeid_side)
        for vals in zip(zs[ sel], Rs[ sel], z_PMT[ sel], phi[ sel], zd[ sel], theta[ sel]):Tside_indirect  .Fill(*vals)

        # Fill direct tables
        # bottom
        sel = ~sel_indirect & np.isin(ihPMT, tubeid_bottom)
        for vals in zip(zs[ sel], Rs[ sel], R_PMT[ sel], phi[ sel]):Tbottom_direct.Fill(*vals)
        # top
        sel = ~sel_indirect & np.isin(ihPMT, tubeid_top)
        for vals in zip(zs[ sel], Rs[ sel], R_PMT[ sel], phi[ sel]):Ttop_direct   .Fill(*vals)
        # side
        sel = ~sel_indirect & np.isin(ihPMT, tubeid_side)
        for vals in zip(zs[ sel], Rs[ sel], z_PMT[ sel], phi[ sel]):Tside_direct  .Fill(*vals)

    # save photon tables
    if args.verbose: print("Writting tables...")

    fout = ROOT.TFile("scattables.root", "RECREATE")
    fout.WriteObject(Tbottom_indirect, "bottom_indirect")
    fout.WriteObject(Tbottom_direct  , "bottom_direct")
    fout.WriteObject(Ttop_indirect   , "top_indirect")
    fout.WriteObject(Ttop_direct     , "top_direct")
    fout.WriteObject(Tside_indirect  , "side_indirect")
    fout.WriteObject(Tside_direct    , "side_direct")

    Tbottom_indirect.DivideUnnormalized4D(Tbottom_direct)
    Ttop_indirect   .DivideUnnormalized4D(Ttop_direct)
    Tside_indirect  .DivideUnnormalized4D(Tside_direct)

    fout.WriteObject(Tbottom_indirect, "botscattable")
    fout.WriteObject(Ttop_indirect   , "topscattable")
    fout.WriteObject(Tside_indirect  , "sidescattable")
    fout.Close()

    return


if __name__ == "__main__":
    main()