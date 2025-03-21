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


def read_photon_information(filename):

    rootf = ROOT.TFile(filename)

    tree = rootf.Get("sttree")
    try                  : n = tree.GetEntries()
    except AttributeError: return None

    srcpos = np.zeros((n, 3), dtype=np.float32)
    srcdir = np.zeros((n, 3), dtype=np.float32)
    ihPMT  = np.zeros(n, dtype=int)
    isct   = np.zeros(n, dtype=int)

    for i in range(n):
        tree.GetEvent(i)
        srcpos[i] = np.frombuffer(tree.srcpos, dtype=np.float32)
        srcdir[i] = np.frombuffer(tree.srcdir, dtype=np.float32)
        ihPMT [i] = tree.ihPMT
        isct  [i] = tree.isct
    
    rootf.Close()
    return srcpos, srcdir, ihPMT, isct


def process_table(infiles, tabname, bins, tubeids, R, mPMTmap, mPMTpositions, mPMTbins, verbose):

    if verbose: print(f">> Creating {tabname:6} table...")

    # Define scattering tables, notice the only change beside name and title is the use of R_PMT (bottom and top) and z_PMT (side)
    h_indirect = np.zeros(shape=tuple(  [len(bins_)-1 for bins_ in bins]))
    h_direct   = np.zeros(shape=tuple([*[len(bins_)-1 for bins_ in bins[:3]], *[len(bins_)-1 for bins_ in bins[-2:]]]))

    # Loop over simulation files in order to fill the table histograms
    for n, filename in enumerate(infiles, 1):
        if verbose: print(f"  {(tabname):7}: Processing file {n}/{len(infiles)}...", basename(filename))

        # read data
        photon_information = read_photon_information(filename)
        if photon_information is None: continue
        srcpos, srcdir, ihPMT, isct = photon_information

        # select photons in tubeids
        sel =  np.isin(ihPMT, tubeids)
        srcpos = srcpos[sel]
        srcdir = srcdir[sel]
        ihPMT  = ihPMT [sel]
        isct   = isct  [sel]

        # get mPMT variables
        mPMTinfo    = mPMTmap.loc[ihPMT, ("mPMTNo", "mPMT_PMTNo")].values.astype(int)
        mPMTNos     = mPMTinfo[:, 0]
        mPMT_PMTNos = mPMTinfo[:, 1]
        mPMTpos     = mPMTpositions.loc[mPMTNos].values
        mPMTzR      = mPMTbins     .loc[mPMTNos].values

        # Rotate with vertical direction given by vaxis
        if R is not None:
            srcpos  = np.matmul(R, srcpos.T).T
            srcdir  = np.matmul(R, srcdir.T).T

        # Compute scattering table variables (transform from mm to cm)
        zs      = 0.1*srcpos[:, 2]
        Rs      = 0.1*np.sqrt(np.sum(srcpos[:, (0, 1)]**2, axis=1))
        phi     = clockwise_azimuth_angle(mPMTpos, srcpos)
        phi_dir = clockwise_azimuth_angle(srcdir, mPMTpos-srcpos)
        z_dir   = srcdir[:, 2]
        
        sel_indirect = (isct != 0)
        # indirect
        sel  = sel_indirect
        h, _ = np.histogramdd((zs[sel], Rs[sel], phi[sel], phi_dir[sel], z_dir[sel], mPMTzR[sel], mPMT_PMTNos[sel]), bins=bins)
        h_indirect += h

        # direct
        sel  = ~sel_indirect
        h, _ = np.histogramdd((zs[sel], Rs[sel], phi[sel], mPMTzR[sel], mPMT_PMTNos[sel]), bins=(*bins[:3], *bins[-2:]))
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
    R = None
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
    zs      = [-length/2., length/2.]
    Rs      = [0., radius]
    phi     = [0., 2.*np.pi]
    phi_dir = [0., 2.*np.pi]
    z_dir   = [-1, +1]
    bounds  = [zs, Rs, phi, phi_dir, z_dir]
    # define bins
    bins = 5*[None]
    for dim, nbins in enumerate(args.nbins):
        bins[dim] = np.linspace(bounds[dim][0], bounds[dim][1], nbins+1)

    # define PMT bins
    pmts_sel = pmts.set_index("TubeNo").loc[tubeids["side"]]
    side_mPMT_zs = pmts_sel.set_index("mPMT_PMTNo").loc[1].set_index("mPMTNo").loc[:, "Position_x2"].astype(float)
    zPMTbinc = np.sort(side_mPMT_zs.unique())
    zPMTbins = (zPMTbinc[1:] + zPMTbinc[:-1])/2.
    zPMTbins = np.insert(zPMTbins,             0, zPMTbinc [0] - (zPMTbins [0]-zPMTbinc [0]))
    zPMTbins = np.insert(zPMTbins, len(zPMTbins), zPMTbinc[-1] - (zPMTbins[-1]-zPMTbinc[-1]))

    pmts_sel = pmts.set_index("TubeNo").loc[tubeids["top"]]
    top_mPMT_Rs = pmts_sel.groupby("mPMTNo").apply(lambda df: np.linalg.norm(df.set_index("mPMT_PMTNo").loc[1].to_frame().loc[["Position_x0", "Position_x1"]]))
    pmts_sel = pmts.set_index("TubeNo").loc[tubeids["bottom"]]
    bottom_mPMT_Rs = pmts_sel.groupby("mPMTNo").apply(lambda df: np.linalg.norm(df.set_index("mPMT_PMTNo").loc[1].to_frame().loc[["Position_x0", "Position_x1"]]))
    # assuming top and bottom positions are the same
    RPMTbinc = np.sort(bottom_mPMT_Rs.unique())
    RPMTbins = (RPMTbinc[1:] + RPMTbinc[:-1])/2.
    RPMTbins = np.insert(RPMTbins,             0, RPMTbinc [0] - (RPMTbins [0]-RPMTbinc [0]))
    RPMTbins = np.insert(RPMTbins, len(RPMTbins), RPMTbinc[-1] - (RPMTbins[-1]-RPMTbinc[-1]))

    # define mPMT map (TubeNo->mPMTNo)
    mPMTmap = pmts.set_index("TubeNo").loc[:, ("mPMTNo", "mPMT_PMTNo")]

    # define mPMT_PMTNobins (implicitly assumes all modules contain the same number of PMTs)
    mPMT_PMTNobinc = mPMTmap.set_index("mPMTNo").loc[1].mPMT_PMTNo.values.astype(int)
    mPMT_PMTNobins = (mPMT_PMTNobinc[1:] + mPMT_PMTNobinc[:-1])/2.
    mPMT_PMTNobins = np.insert(mPMT_PMTNobins, 0, mPMT_PMTNobinc[0] - (mPMT_PMTNobins[0]-mPMT_PMTNobinc[0]))
    mPMT_PMTNobins = np.insert(mPMT_PMTNobins, len(mPMT_PMTNobins), mPMT_PMTNobinc[-1] - (mPMT_PMTNobins[-1]-mPMT_PMTNobinc[-1]))

    # group the mPMT bins into a dictionary
    mPMTbins = dict()
    mPMTbins["side"]   = side_mPMT_zs
    mPMTbins["top"]    = top_mPMT_Rs
    mPMTbins["bottom"] = bottom_mPMT_Rs

    # define mPMT positions as the central PMT position
    mPMTpositions = pmts.set_index("mPMT_PMTNo").loc[1, ("mPMTNo", "Position_x0", "Position_x1", "Position_x2")].astype(float).set_index("mPMTNo")

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
    fout.create_carray(bins_group, "bins_mPMT_PMTNo", atom, mPMT_PMTNobins.shape, obj=mPMT_PMTNobins, filters=filters)
    fout.flush()
    # define tables group tables
    tables_group = fout.create_group("/", "tables", "tables")

    # paralellize loop over tables to fill histograms
    results         = dict()
    future_to_table = dict()

    with concurrent.futures.ProcessPoolExecutor(max_workers=os.cpu_count()) as executor:
        for tabname in tabnames:
            bins_ = bins.copy()
            if tabname == "side": bins_.append(zPMTbins)
            else                : bins_.append(RPMTbins)
            bins_.append(mPMT_PMTNobins)
            future_to_table[executor.submit(process_table, infiles, tabname, bins_, tubeids[tabname], R, mPMTmap, mPMTpositions, mPMTbins[tabname], args.verbose)] = tabname

        # get results once all the tasks are finished
        for future in concurrent.futures.as_completed(future_to_table):
            h_indirect, h_direct = future.result()
            tabname = future_to_table[future]
            results[tabname] = (h_indirect, h_direct)

    for tabname in tabnames:
        # bins_ = bins.copy()
        # # add zRPMTbins
        # if tabname == "side": bins_.append(zPMTbins)
        # else                : bins_.append(RPMTbins)
        # # add mPMT_PMTNo bins
        # bins_.append(mPMT_PMTNobins)
        # results[tabname] = process_table(infiles, tabname, bins_, tubeids[tabname], R, mPMTmap, mPMTpositions, mPMTbins[tabname], args.verbose)

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