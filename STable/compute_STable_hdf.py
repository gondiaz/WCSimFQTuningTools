import os
import sys
import glob
import argparse
import uproot
import ROOT
import concurrent.futures
import numpy  as np
import tables as tb

from os.path   import expandvars, join, basename

from STable.STable_tools import split_tubeids, clockwise_azimuth_angle, read_wcsim_geometry


def process_table(infiles, tabname, bins, tubeids, vaxis, R, verbose):

    if verbose: print(f">> Creating {tabname:6} table...")

    # Define scattering tables, notice the only change beside name and title is the use of R_PMT (bottom and top) and z_PMT (side)
    h_indirect = np.zeros(shape=tuple(  [len(bins_)-1 for bins_ in bins]))
    h_direct   = np.zeros(shape=tuple([*[len(bins_)-1 for bins_ in bins[:3]], len(bins[-1])-1]))

    # Loop over simulation files in order to fill the table histograms
    for n, filename in enumerate(infiles, 1):
        if verbose: print(f"  {(tabname):7}: Processing file {n}/{len(infiles)}...", basename(filename))

        # read data
        with uproot.open(filename) as f:
            try:
                # photons reaching a PMT either directly or indirectly
                # parent pos and dir
                srcpos = f["sttree/srcpos"].array().to_numpy()
                srcdir = f["sttree/srcdir"].array().to_numpy()
                # PMT number/pos
                ihPMT   = f["sttree/ihPMT"]  .array().to_numpy()
                tubepos = f["sttree/tubepos"].array().to_numpy()
                # ISTORY=(# of refl)*1000 + (# of scat); 0 if direct hit
                isct = f["sttree/isct"].array().to_numpy()
            except uproot.exceptions.KeyInFileError: continue

        # Rotate with vertical direction given by vaxis
        if vaxis != 2:
            srcpos  = np.matmul(R,  srcpos.T).T
            srcdir  = np.matmul(R,  srcdir.T).T
            tubepos = np.matmul(R, tubepos.T).T

        # Compute scattering table variables (transform from mm to cm)
        zs    = 0.1*srcpos[:, 2]
        Rs    = 0.1*np.sqrt(np.sum(srcpos[:, (0, 1)]**2, axis=1))
        phi   = clockwise_azimuth_angle(tubepos, srcpos)
        zd    = srcdir[:, 2]
        theta = clockwise_azimuth_angle(srcdir, tubepos-srcpos)
        # z_PMT used for side table, R_PMT used for bottom and top
        z_PMT = 0.1*tubepos[:, 2]
        R_PMT = 0.1*np.sqrt(np.sum(tubepos[:, (0, 1)]**2, axis=1))

        # transform angles from [0, 2pi] to [-pi, pi]
        phi  [  phi>np.pi] -= 2.*np.pi
        theta[theta>np.pi] -= 2.*np.pi

        # histograms
        if tabname == "side": pmtpos = z_PMT
        else                : pmtpos = R_PMT
        sel_indirect = (isct != 0)
        # indirect
        sel  = sel_indirect & np.isin(ihPMT, tubeids)
        h, _ = np.histogramdd((zs[sel], Rs[sel], phi[sel], zd[sel], theta[sel], pmtpos[sel]), bins=bins)
        h_indirect += h
        # direct
        sel  = ~sel_indirect & np.isin(ihPMT, tubeids)
        h, _ = np.histogramdd((zs[sel], Rs[sel], phi[sel], pmtpos[sel]), bins=(*bins[:3], bins[-1]))
        h_direct += h

    # Add direction axis to direct light histogram
    h_direct = np.expand_dims(h_direct, axis=(3, 4))
    return h_indirect, h_direct


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
    parser.add_argument("--vaxis", type=int, nargs="?", help = "detector vertical axis (0)", default=2)
    parser.add_argument("--nbins", type=int, nargs=5, help = "histogram bins in the following order: zs,Rs,phi,zd,theta")
    parser.add_argument("--maxnfiles", type=int, nargs="?", help = "max files to proccess", default=None)
    
    args = parser.parse_args()
    ##########################################

    # Load WCSimRoot.so (needed by STable_tools.py functions)
    ROOT.gSystem.AddDynamicPath(expandvars(args.wcsimlib[0]))
    ROOT.gSystem.Load          ("libWCSimRoot.dylib" if sys.platform == "darwin" else "libWCSimRoot.so")

    # get simulation filenames, removing those ending in '_flat.root'
    infiles = glob.glob(join(args.indir[0], "*"))
    infiles = sorted([f for f in infiles if "_flat.root" not in basename(f)])
    if args.maxnfiles: infiles = infiles[:args.maxnfiles]

    # get bottom, top and side tube-ids from the first file
    tabnames = ["bottom", "top", "side"]
    tubeids = split_tubeids(infiles[0], vaxis=args.vaxis)
    tubeids = dict(zip(tabnames, tubeids))

    # Select rotation matrix based on vertical axis
    if   args.vaxis == 0:
        R = np.array([[0, 0, -1], [0, 1, 0], [1, 0, 0]])
    elif args.vaxis == 1:
        R = np.array([[1, 0, 0], [0, 0, -1], [0, 1, 0]])
    elif args.vaxis == 2: pass
    else: raise Exception("Invalid value for vertical axis (vaxis)")

    # get detector length and radius
    df, pmts = read_wcsim_geometry(infiles[0])
    length = df.loc["WCDetHeight", "WC"]
    radius = df.loc["WCDetRadius", "WC"]

    # rotate pmt positions and orientations
    poscolumns = [f"Position_x{i}"    for i in range(3)]
    orscolumns = [f"Orientation_x{i}" for i in range(3)]
    if args.vaxis != 2:
        pmts.loc[:, poscolumns] = np.matmul(R, pmts.loc[:, poscolumns].values.T).T
        pmts.loc[:, orscolumns] = np.matmul(R, pmts.loc[:, orscolumns].values.T).T

    # Define scatering table binning
    # 6 variables
    # zs: source z position
    # Rs: source R position
    # phi: angle between source and PMT
    # zd : source z direction
    # theta: source direction angle w.r.t PMT
    # (bottom and top table) R_PMT
    # (side table) z_PMT: zposition of the PMT
    # define bounds
    zs    = [-length/2., length/2.]
    Rs    = [0., radius]
    phi   = [-np.pi, +np.pi]
    zd    = [-1., 1.]
    theta = [-np.pi, +np.pi]
    bounds = [zs, Rs, phi, zd, theta]
    # define bins
    bins = 5*[None]
    for dim, nbins in enumerate(args.nbins):
        bins[dim] = np.linspace(bounds[dim][0], bounds[dim][1], nbins+1)

    # define PMT bins
    pmts_sel = pmts.set_index("TubeNo").loc[tubeids["side"]]
    zPMTbinc = np.sort(pmts_sel.Position_x2.unique().astype(float))
    zPMTbins = (zPMTbinc[1:] + zPMTbinc[:-1])/2.
    zPMTbins = np.insert(zPMTbins,             0, zPMTbinc [0] - (zPMTbins [0]-zPMTbinc [0]))
    zPMTbins = np.insert(zPMTbins, len(zPMTbins), zPMTbinc[-1] - (zPMTbins[-1]-zPMTbinc[-1]))

    pmts_sel = pmts.set_index("TubeNo").loc[tubeids["top"]] # assuming top and bottom positions are the same
    RPMTbinc = np.sort(np.unique(np.sum(pmts_sel.loc[:, ("Position_x0", "Position_x1")]**2, axis=1)**0.5).astype(float))
    RPMTbins = (RPMTbinc[1:] + RPMTbinc[:-1])/2.
    RPMTbins = np.insert(RPMTbins,             0, RPMTbinc [0] - (RPMTbins [0]-RPMTbinc [0]))
    RPMTbins = np.insert(RPMTbins, len(RPMTbins), RPMTbinc[-1] - (RPMTbins[-1]-RPMTbinc[-1]))

    # compute scale factor for direct light 4D->6D (assumes uniform bin size) 
    # given by dx/Dx == nbins for the direction dimensions
    scale_factor = 1.
    for i in range(3, 4+1): scale_factor *= (bins[i][1]-bins[i][0])/(bins[i][-1] - bins[i][0])

    # Open output file
    filters = tb.Filters(complevel=5, complib="zlib")
    atom = tb.Float64Atom()
    fout = tb.open_file("scattables.h5", mode="w", title="STable")
    # Write bins to output file
    bins_group = fout.create_group("/", "bins", "binning")
    for dim, bins_ in enumerate(bins):
        fout.create_carray(bins_group, f"bins_{dim}", atom, bins_.shape, obj=bins_, filters=filters)
    fout.create_carray(bins_group, "bins_zPMT", atom, zPMTbins.shape, obj=zPMTbins, filters=filters)
    fout.create_carray(bins_group, "bins_RPMT", atom, RPMTbins.shape, obj=RPMTbins, filters=filters)
    fout.flush()
    # define tables group tables
    tables_group = fout.create_group("/", "tables", "tables")

    # paralellize loop over tables to fill histograms
    results         = dict()
    future_to_table = dict()

    # with concurrent.futures.ProcessPoolExecutor(max_workers=os.cpu_count()) as executor:
    #     for tabname in tabnames:
    #         bins_ = bins.copy()
    #         if tabname == "side": bins_.append(zPMTbins)
    #         else                : bins_.append(RPMTbins)
    #         future_to_table[executor.submit(process_table, infiles, tabname, bins_, tubeids[tabname], args.vaxis, R, args.verbose)] = tabname

    #     # get results once all the tasks are finished
    #     for future in concurrent.futures.as_completed(future_to_table):
    #         h_indirect, h_direct = future.result()
    #         tabname = future_to_table[future]
    #         results[tabname] = (h_indirect, h_direct)

    for tabname in tabnames:
        bins_ = bins.copy()
        if tabname == "side": bins_.append(zPMTbins)
        else                : bins_.append(RPMTbins)
        results[tabname] = process_table(infiles, tabname, bins_, tubeids[tabname], args.vaxis, R, args.verbose)

        # Scattering tables
        if args.verbose: print(f">> Writting {tabname} table...")

        h_indirect = results[tabname][0]
        h_direct   = results[tabname][1]

        g = fout.create_group(tables_group, tabname, "tables")
        fout.create_carray(g, "indirect", atom, h_indirect.shape, obj=h_indirect, filters=filters)
        fout.create_carray(g,   "direct", atom,   h_direct.shape, obj=  h_direct, filters=filters)
        stable = np.divide( h_indirect, h_direct*scale_factor
                          , out  = np.zeros_like(h_indirect)
                          , where= (h_direct != 0))
        fout.create_carray(g, "stable", atom, stable.shape, obj=stable, filters=filters)
        fout.flush()

    fout.close()
    return


if __name__ == "__main__":
    main()