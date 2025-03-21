import os
import re
import sys
import glob
import argparse
import warnings
import uproot
import ROOT
import configparser
import numpy  as np
import pandas as pd
import tables as tb
import concurrent.futures

from itertools import groupby
from os.path   import expandvars, join, basename, exists

from STable.STable_tools import split_tubeids, clockwise_azimuth_angle, read_wcsim_geometry

from fitqun.logger                   import get_logger, c0
from fitqun_fit.predicted_charges    import compute_predicted_charges, is_contained
from fitqun_io.tuning_charge_readers import CProfile, AngularResponse, ScatteringTables

logger = get_logger(__name__)

warnings.filterwarnings("ignore", category=tb.NaturalNameWarning)

def process_momentum( filenames, outfilename
                    , pmts, R, radius, length, pmtradius, trigger_offset
                    , attenuation_length, refraction_index, QE
                    , cprof, ang, scat
                    , tresbins, μbins
                    , max_contained_statistics):

    # Create tuning file instances
    ang   = AngularResponse (ang, pmtradius, attenuation_length)
    scat  = ScatteringTables(scat)
    cprof = CProfile        (cprof)

    # Define empty variables (useful later)
    x = np.zeros(7)
    vertex    = np.zeros((1, 3))
    direction = np.zeros((1, 3))

    # read momentum from first event (instead of computing it from the energy)
    rootf = ROOT.TFile(filenames[0])
    tree  = rootf.GetKey("wcsimT").ReadObj()
    tree.GetBranch("wcsimrootevent").SetAutoDelete(True)
    tree.GetEvent(0)
    trigger = tree.wcsimrootevent.GetTrigger(0)
    tracks = trigger.GetTracks()

    # find parent track
    for i in range(tracks.GetEntries()):
        track = tracks[i]
        if (track.GetParenttype()==0)&(track.GetId()==1):
            break
    if ((track.GetParenttype()!=0) or (track.GetId()!=1)):
        logger.error("No parent track for first event")
    momentum = track.GetP()
    rootf.Close()

    if (np.log(momentum)<=cprof.logpbins[0]) | (cprof.logpbins[-1]<=np.log(momentum)):
        logger.warning(f"Momentum {momentum} out of tuning limits. Skipping.")
        return

    # get max travel distance
    max_distance = cprof.get_maxdistance(np.log(momentum))

    logger.info(f"Processing momentum {round(momentum, 2)} MeV/c...")

    # define histograms to be filled
    hdirect   = np.zeros((len(tresbins)-1, len(μbins)-1))
    hindirect = np.zeros((len(tresbins)-1, len(μbins)-1))
    contained_counter = 0
    total_counter     = 0
    for filename in filenames:

        logger.info("Processing %s", basename(filename))

        rootf = ROOT.TFile(filename)
        tree  = rootf.GetKey("wcsimT").ReadObj()
        tree.GetBranch("wcsimrootevent").SetAutoDelete(True)
        nevents = tree.GetEntries()

        for event_i in range(nevents):

            tree.GetEvent(event_i)
            ntriggers = tree.wcsimrootevent.GetNumberOfEvents()

            trigger = tree.wcsimrootevent.GetTrigger(0) # only trigger 0 is considered
            tracks = trigger.GetTracks()

            # find parent track
            for i in range(tracks.GetEntries()):
                track = tracks[i]
                if ((track.GetParenttype()==0)&(track.GetId()==1)):
                    break

            if ((track.GetParenttype()!=0) or (track.GetId()!=1)):
                logger.error(f"No parent track for event index {event_i}")
                continue

            # get vertex and direction
            for i in range(3):
                vertex   [0, i] = track.GetStart(i)
                direction[0, i] = track.GetDir  (i)

            if R is not None:
                vertex    = np.matmul(R,    vertex.T).T
                direction = np.matmul(R, direction.T).T

            x[1:4] = vertex
            x[4]   = np.arccos (direction[:, 2])[0]
            x[5]   = np.arctan2(direction[:, 1], direction[:, 0])[0]
            x[6]   = track.GetP()

            # if is contained, process the event
            if is_contained(vertex[0], direction[0], max_distance, radius, length):
                contained_counter += 1

                # get digi hits
                digihits = trigger.GetCherenkovDigiHits()

                # loop over hits and compute residual times
                # tresidual = np.zeros((trigger.GetNcherenkovdigihits(), 2)) # saves residual time and tubeid
                tresidual = []
                midpos = vertex + direction * (max_distance/2.)
                t      = trigger_offset - trigger.GetHeader().GetDate()
                for hit_i in range(trigger.GetNcherenkovdigihits()):
                    
                    digihit = digihits[hit_i]
                    tubeid  = digihit.GetTubeId()

                    if tubeid not in pmts.index.values:
                        continue

                    pmtpos = pmts.loc[tubeid].loc[["x", "y", "z"]].values

                    # save residual time and tubeid
                    midtrack_pmt_distance = np.linalg.norm(pmtpos - midpos)
                    # tresidual[hit_i, 0] = tubeid
                    # tresidual[hit_i, 1] = digihit.GetT() - t - (max_distance/2.)/c0 - midtrack_pmt_distance*refraction_index/c0
                    tresidual.append([ tubeid
                                     , digihit.GetT() - t - (max_distance/2.)/c0 - midtrack_pmt_distance*refraction_index/c0])

                # sort by tubeid
                if len(tresidual) == 0:
                    continue
                tresidual = np.array(tresidual)
                sorted_indices = np.argsort(tresidual[:, 0])
                tresidual = tresidual[sorted_indices]

                # compute predicted charges
                μ_direct, μ_indirect = compute_predicted_charges( x, pmts
                                                                , ang, cprof, scat
                                                                , attenuation_length
                                                                , refraction_index, QE)
                
                # select only hitted pmts
                sel = np.isin(pmts.index, tresidual[:, 0])
                μ_direct   = μ_direct  [sel]
                μ_indirect = μ_indirect[sel]

                # fill histograms (add 1e-20 to silence warning)
                h, _, _ = np.histogram2d(tresidual[:, 1], np.log10(μ_direct   + 1.e-20), bins=[tresbins, μbins])
                hdirect += h
                h, _, _ = np.histogram2d(tresidual[:, 1], np.log10(μ_indirect + 1.e-20), bins=[tresbins, μbins])
                hindirect += h

            total_counter += 1

            # check enough statistics
            if (contained_counter >= max_contained_statistics): break

        # check enough statistics
        rootf.Close()
        if (contained_counter >= max_contained_statistics): break

    # Save results
    logger.info(f"Saving {round(momentum, 2)} MeV/c histograms")
    with tb.open_file(outfilename, "a") as f:
        f.create_array(f.root.direct  , f"p_{round(momentum, 3)}", hdirect)
        f.create_array(f.root.indirect, f"p_{round(momentum, 3)}", hindirect)

        f.root.events.append([(momentum, contained_counter, total_counter)])
        f.flush()

    return



def main():

    ############ Program arguments ############
    parser = argparse.ArgumentParser( prog        = f"{basename(__file__)}"
                                    , description = "description"
                                    , epilog      = """ """)
    
    parser.add_argument("-v", "--verbose", action="store_true")
    parser.add_argument("--wcsimlib"        , type=str,  nargs=1  , help = "WCSim lib path")
    parser.add_argument("-i", "--infiles"   , type=str,  nargs="+", help = "Defines the input files"     , required=True)
    parser.add_argument("-o", "--outdir"    , type=str            , help = "Defines the output directory", required=True)
    parser.add_argument("-c", "--configfile", type=str, nargs=None, help = "configuration parameters"    , required=True)
    parser.add_argument("-n", "--max-contained-statistics", type=int, default=float("inf"), help = "Number of events to process")

    args = parser.parse_args()
    ##########################################

    # Load libWCSimRoot library
    ROOT.gSystem.AddDynamicPath(expandvars(args.wcsimlib[0]))
    ROOT.gSystem.Load          ("libWCSimRoot.dylib" if sys.platform == "darwin" else "libWCSimRoot.so")

    # filenames sorter function
    get_energy_and_index = lambda fname: list(map(float, re.findall(r'\d+(?:\.\d+)?', basename(fname))[:2]))

    # read configuration file and get parameters
    if not exists(args.configfile):
        logger.critical("Missing config file")
    logger.info("Reading parameters from configuration file")
    config = configparser.ConfigParser()
    config.read(args.configfile)

    # define output filename
    outfilename = join(expandvars(args.outdir), "tpdf_histograms.h5")
    # maximum number of contained events used for each histogram
    max_contained_statistics = args.max_contained_statistics

    # group files in energies
    infiles = [f for f in args.infiles if re.match("out_.+_\d+(?:\.\d+)?_\d+.root", basename(f))]
    infiles = sorted(infiles, key=get_energy_and_index)
    grouped_filenames = [list(group) for key, group in groupby(infiles, key=lambda x: get_energy_and_index(x)[0])]
    energies = [get_energy_and_index(filenames[0])[0] for filenames in grouped_filenames]

    # get particle type, assert all files contain same particle
    particle = basename(infiles[0]).split("_")[1]
    for filename in infiles: 
        assert basename(filename).split("_")[1] == particle

    # get needed data from configuration file
    vaxis              =   int(config["detector"].get("vaxis"))
    attenuation_length = float(config["physics"].get("attenuation_length"))
    refraction_index   = float(config["physics"].get("refraction_index"))
    QE                 = float(config["physics"].get("quantum_efficiency"))

    # Select rotation matrix based on vertical axis
    R = None
    if   vaxis == 0:
        R = np.array([[0, 0, -1], [0, 1, 0], [1, 0, 0]])
    elif vaxis == 1:
        R = np.array([[1, 0, 0], [0, 0, -1], [0, 1, 0]])
    elif vaxis == 2: pass
    else: raise Exception("Invalid value for vertical axis (vaxis)")

    # get a random filename and extract needed pmt information
    logger.info("Reading PMTs information...")
    df, pmts = read_wcsim_geometry(infiles[0])
    pmts.drop("CylLoc", axis=1, inplace=True)
    pmts.rename({"Position_x0": "x", "Position_x1": "y", "Position_x2": "z"}, axis=1, inplace=True)
    pmts.rename({"Orientation_x0": "ox", "Orientation_x1": "oy", "Orientation_x2": "oz"}, axis=1, inplace=True)
    pmts.rename({"TubeNo": "tubeid", "mPMTNo": "mPMTid", "mPMT_PMTNo": "mPMT_tubeid"}, axis=1, inplace=True)
    pmts.sort_values(by="tubeid", inplace=True)
    # select only central PMTs
    pmts = pmts.loc[pmts.mPMT_tubeid == 1]

    pmtradius = float(df.loc["WCPMTRadius", "WC"])
    radius    = float(df.loc["WCDetRadius", "WC"])
    length    = float(df.loc["WCDetHeight", "WC"])

    # rotate positions and orientations
    if R is not None:
        logger.info("Rotating detector")
        pos = pmts.loc[:, ( "x",  "y",  "z")].values
        ors = pmts.loc[:, ("ox", "oy", "oz")].values
        pmts.loc[:, ( "x",  "y",  "z")] = np.matmul(R, pos.T).T
        pmts.loc[:, ("ox", "oy", "oz")] = np.matmul(R, ors.T).T

    # add tube labels to pmts pd.DataFrame
    tubeid_bottom, tubeid_top, tubeid_side = split_tubeids(infiles[0], vaxis=vaxis)
    pmts["label"] = 0
    pmts.loc[np.isin(pmts.tubeid, tubeid_side  ), "label"] = 0
    pmts.loc[np.isin(pmts.tubeid, tubeid_bottom), "label"] = 1
    pmts.loc[np.isin(pmts.tubeid, tubeid_top   ), "label"] = 2
        
    pmts.set_index("tubeid", inplace=True)

    # create tuning instances
    logger.info("Reading charge tuning file names...")
    ang   = expandvars(config["tuning_files"]["Angular"])
    scat  = expandvars(config["tuning_files"]["STable"])
    cprof = expandvars(config["tuning_files.cprofiles"][particle])

    # define 2D histogram bins
    ntresbins = int  (config["histogram_bins"]["ntresbins"])
    nμbins    = int  (config["histogram_bins"]["nμbins"])
    treslow   = float(config["histogram_bins"]["treslow"])
    tresup    = float(config["histogram_bins"]["tresup"])
    μlow      = float(config["histogram_bins"]["μlow"])
    μup       = float(config["histogram_bins"]["μup"])
    tresbins = np.linspace(treslow, tresup, ntresbins+1)
    μbins    = np.linspace(   μlow,    μup, nμbins   +1)

    # get trigger offset
    rootf = ROOT.TFile(infiles[0])
    tree = rootf.GetKey("wcsimRootOptionsT").ReadObj()
    tree.GetEvent(0)
    options = tree.wcsimrootoptions
    trigger_offset = options.GetTriggerOffset()
    rootf.Close()

    logger.info(f"Got trigger offset {trigger_offset} ns")

    # create groups for direct and indirect pdfs
    logger.info(f"Creating output file in running directory: {basename(outfilename)}")
    with tb.open_file(outfilename, "w") as f:
        f.create_group("/",   "direct",   "direct light histogram")
        f.create_group("/", "indirect", "indirect light histogram")
        f.create_group("/",     "bins", "histogram bins")
        f.create_array("/bins", "tres", tresbins)
        f.create_array("/bins",    "μ", μbins)

        f.create_earray("/",  "events", shape=(0, 3), atom=tb.Int64Atom(), expectedrows=len(energies))

    ##########################################
    # loop on energies (parallelized)
    ##########################################
    logger.info("Launching parallel processes")
    with concurrent.futures.ProcessPoolExecutor() as executor:
        for _, filenames in zip(energies, grouped_filenames):
            executor.submit(process_momentum, filenames, outfilename
                                            , pmts, R, radius, length, pmtradius, trigger_offset
                                            , attenuation_length, refraction_index, QE
                                            , cprof, ang, scat
                                            , tresbins, μbins
                                            , max_contained_statistics)

    logger.info("Processing completed succesfully")
    return


if __name__ == "__main__":
    main()
