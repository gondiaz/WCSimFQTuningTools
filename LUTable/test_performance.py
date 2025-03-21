"""
    tool for temporal studies yet, notice it depends on Pyfitqun
"""
import time
import re
import glob
import argparse

from os.path   import expandvars, basename, join
from itertools import groupby, product
from functools import partial

import tables as tb
import pandas as pd
import numpy  as np

from fitqun_io.input_readers         import WCSimReader
from fitqun_io.tuning_charge_readers import CProfile
from fitqun_fit.predicted_charges    import is_contained
from fitqun.logger                   import get_logger

logger = get_logger(__name__)


def get_charges(filename, cprofile):

    # by convention, the vertical direction of the cylinder is z
    vdir = np.array([0., 0., 1.])

    # initialize reader
    reader = WCSimReader(filename)
    radius = reader.radius
    length = reader.length
    # select only side pmts since for top and bottom the z axis corresponds with the normal
    pmts   = reader.pmts.set_index("tubeid").loc[reader.tubeid_side].reset_index()

    # read mc info, assert equal momentum and get travel distance
    tracks = reader.read_mc_truth()
    p = tracks.loc[0, "P"]
    assert all(tracks.P.values == tracks.loc[0, "P"])
    dmax = cprofile.get_maxdistance(np.log(p))/0.9

    rs, dphis, thpos, thdir, qs = [], [], [], [], []

    # loop on events
    while True:
        next_event = reader.read_next_event_hits()
        if next_event is None:
            break
        event, hits = next_event

        # check if contained
        pos = tracks.loc[tracks.EvtNum == event, ("Start_x0", "Start_x1", "Start_x2")].values[0]
        dir = tracks.loc[tracks.EvtNum == event, (  "Dir_x0",   "Dir_x1",   "Dir_x2")].values[0]

        if not is_contained(pos, dir, dmax, radius, length):
            logger.warning("Event %s not contained", event)
            continue

        # create allhits dataframe such that unhited pmts contain zero charge
        hits = hits.set_index("tubeid").loc[np.intersect1d(hits.tubeid, pmts.index)].reset_index()
        allhits = pd.DataFrame(dict({col: pd.Series(dtype=typ)
                                     for col, typ in zip(hits.columns, hits.dtypes)}))
        allhits.tubeid = pmts.tubeid
        allhits.mPMTid = pmts.mPMTid
        allhits.q      = 0.
        allhits.t      = 0.
        allhits.set_index("tubeid", inplace=True)
        allhits.loc[:, ("x", "y", "z")]    = pmts.loc[:, ("x", "y", "z")]
        allhits.loc[:, ("ox", "oy", "oz")] = pmts.loc[:, ("ox", "oy", "oz")]
        allhits.loc[hits.tubeid, "q"] = hits.set_index("tubeid").q
        allhits.loc[hits.tubeid, "t"] = hits.set_index("tubeid").t

        # compute distance and angular variables
        ors   = allhits.loc[:, ("ox", "oy", "oz")].values
        xdirs = np.cross(ors, vdir)

        r   = allhits.loc[:, ("x", "y", "z")].values - pos
        rs_ =  np.sqrt(np.sum(r**2, axis=1))

        # azimuth angles (phi): [-pi, pi] --> [0, 2pi]
        phipos_ = np.arctan2(  np.matmul(r, vdir), np.einsum( "ij,ij->i", r, xdirs))
        phidir_ = np.arctan2(np.matmul(dir, vdir), np.einsum( "j,ij->i", dir, xdirs))
        phipos_[phipos_ < 0.] += 2.*np.pi
        phidir_[phidir_ < 0.] += 2.*np.pi
        dphis_  = abs(phipos_ - phidir_)

        # polar angles (theta) (using the cosine)
        thpos_ = np.matmul(  r, vdir)/rs_
        thdir_ = np.matmul(dir, vdir)/rs_

        rs   .extend(rs_)
        dphis.extend(dphis_)
        thpos.extend(thpos_)
        thdir.extend(thdir_)
        qs   .extend(allhits.q.values)

    rs    = np.array(rs)
    dphis = np.array(dphis) / (2.*np.pi)
    thpos = np.array(thpos)
    thdir = np.array(thdir)
    qs    = np.array(qs)
    return rs, dphis, thpos, thdir, qs


def main():

    ############ Program arguments ############
    parser = argparse.ArgumentParser( prog        = "get_conversion_distributions"
                                    , description = ""
                                    , epilog      = """""")

    parser.add_argument("-v", "--verbose", action="store_true")
    parser.add_argument( "indir", type=str, help = "directory containing input files")

    parser.add_argument(  "--cprofile", type=str, help = "")

    args = parser.parse_args()
    ##########################################

    # define input filename structure
    filename_pattern     = "out_([^_]+)_\d+(?:\.\d+)?_\d+.h5"
    get_energy_and_index = lambda fname: list(map(float, re.findall(r'\d+(?:\.\d+)?', basename(fname))[:2]))

    # get all files in input dir, filter out the '_flat.root' and sort them based on (p, index)
    infiles = glob.glob(join(args.indir, "*"))
    infiles = [f for f in infiles if re.match(filename_pattern, basename(f))]
    infiles = sorted(infiles, key=get_energy_and_index)

    # split input files in mu groups
    energies = np.unique([get_energy_and_index(f)[0] for f in infiles])
    groups = [list(group) for key, group in groupby(infiles, key=lambda x: get_energy_and_index(x)[0])]

    # instantiate CProfile
    cprofile = CProfile(args.cprofile)

    # loop on energies
    for (energy, filenames) in zip(energies, groups):
        if args.verbose:
            logger.info("--> Processing energy %s MeV", energy)

        # loop on files
        try:
            # the actual loop on files
            for filename in filenames:
                t0 = time.time()
                if args.verbose:
                    logger.info("--> Processing filename %s", basename(filename))

                # get charges
                rs, dphis, thpos, thdirs, qs = get_charges(filename, cprofile)

                if args.verbose:
                    logger.info("Processing time %f s", (time.time() - t0))

        except KeyboardInterrupt:
            logger.warning("Proccesing interrupted by user")
    return



if __name__ == "__main__":
    main()