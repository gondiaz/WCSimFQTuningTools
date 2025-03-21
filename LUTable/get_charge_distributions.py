"""
    tool for temporal studies yet, notice it depends on Pyfitqun
"""
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
        hits = hits.set_index("tubeid").loc[np.intersect1d(hits.tubeid, pmts.tubeid)].reset_index()
        allhits = pd.DataFrame(dict({col: pd.Series(dtype=typ)
                                     for col, typ in zip(hits.columns, hits.dtypes)}))
        allhits.tubeid = pmts.tubeid
        allhits.mPMTid = pmts.mPMTid
        allhits.q      = 0.
        allhits.t      = 0.
        allhits.set_index("tubeid", inplace=True)
        allhits.loc[:, ( "x",  "y",  "z")] = pmts.set_index("tubeid").loc[:, ( "x",  "y",  "z")]
        allhits.loc[:, ("ox", "oy", "oz")] = pmts.set_index("tubeid").loc[:, ("ox", "oy", "oz")]
        allhits.loc[hits.tubeid, "q"]      = hits.set_index("tubeid").q
        allhits.loc[hits.tubeid, "t"]      = hits.set_index("tubeid").t

        # compute distance and angular variables
        ors   = allhits.loc[:, ("ox", "oy", "oz")].values
        xdirs = np.cross(ors, vdir)

        r   = pos - allhits.loc[:, ("x", "y", "z")].values
        rs_ =  np.sqrt(np.sum(r**2, axis=1))

        # azimuth angles (phi), transform from [-pi, pi] --> [0, 2pi] for direction
        phipos_ = np.arctan2(np.einsum("ij,ij->i",   r, ors), np.einsum("ij,ij->i",   r, xdirs))
        phidir_ = np.arctan2(np.einsum( "j,ij->i", dir, ors), np.einsum( "j,ij->i", dir, xdirs))
        phidir_[phidir_ < 0.] += 2.*np.pi
        assert np.all((0 <= phipos_) & (phipos_ <= np.pi))
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
    parser = argparse.ArgumentParser( prog        = f"{basename(__file__)}"
                                    , description = ""
                                    , epilog      = """""")

    parser.add_argument("-v", "--verbose", action="store_true")
    parser.add_argument("-i", "--infiles", type=str, nargs="+"
                       , help="Defines the input files", required=True)
    parser.add_argument("--outfile", type=str
                        , help = "name of output file", required=True)
    
    parser.add_argument(  "--cprofile", type=str, help = "", required=True)
    parser.add_argument(      "--rmax", type=int, help = "", default=100)
    parser.add_argument(    "--nrbins", type=int, help = "source-pmt distance"           , default=10)
    parser.add_argument( "--ndphibins", type=int, help = "source angle wrt pmt normal"   , default=10)
    parser.add_argument("--nthposbins", type=int, help = "source angle wrt vertical axis", default=10)
    parser.add_argument("--nthdirbins", type=int, help = ""                              , default=10)
    parser.add_argument(      "--Nmax", type=int
                       , help = "maximum number of charge values to save per bin", default=100)

    args = parser.parse_args()
    ##########################################

    # define input filename structure
    filename_pattern     = "out_([^_]+)_\d+(?:\.\d+)?_\d+.h5"
    get_energy_and_index = lambda fname: list(map(float, re.findall(r'\d+(?:\.\d+)?', basename(fname))[:2]))

    # sort input files in (energy, index)
    infiles = [expandvars(f) for f in args.infiles]
    infiles = [f for f in infiles if re.match(filename_pattern, basename(f))]
    infiles = sorted(infiles, key=get_energy_and_index)

    # split input files in mu groups. Assumes only one energy value.
    energies = np.unique([get_energy_and_index(f)[0] for f in infiles])
    assert len(energies) == 1

    # instantiate CProfile
    cprofile = CProfile(args.cprofile)

    # create outfile
    rbins     = np.linspace( 0., args.rmax, args.nrbins    +1)
    dphibins  = np.linspace( 0.,        1., args.ndphibins +1)
    thposbins = np.linspace(-1.,        1., args.nthposbins+1)
    thdirbins = np.linspace(-1.,        1., args.nthdirbins+1)

    # loop on files
    with tb.open_file(expandvars(args.outfile), mode='w') as f:
        f.create_array(f.root,     "rbins", rbins)
        f.create_array(f.root,  "dphibins", dphibins)
        f.create_array(f.root, "thposbins", thposbins)
        f.create_array(f.root, "thdirbins", thdirbins)

        # 2D array storing predicted charges
        charges = f.create_array( f.root, 'charges', atom = tb.Float32Atom(args.Nmax)
                                , shape=( args.nrbins, args.ndphibins
                                        , args.nthposbins, args.nthdirbins))
        # 2D array storing statistics
        stats = f.create_array( f.root, 'stats', atom = tb.IntAtom()
                                , shape=( args.nrbins, args.ndphibins
                                        , args.nthposbins, args.nthdirbins))
        # fill arrays
        for ri, dphii, thposi, thdiri in product( range(args.nrbins)
                                                , range(args.ndphibins)
                                                , range(args.nthposbins)
                                                , range(args.nthdirbins)):
            charges[ri, dphii, thposi, thdiri] = -1
            stats  [ri, dphii, thposi, thdiri] = 0

        try:
            # the actual loop on files
            for filename in infiles:
                if args.verbose:
                    logger.info("--> Processing filename %s", basename(filename))

                # get charges
                rs, dphis, thpos, thdirs, qs = get_charges(filename, cprofile)

                # select values within bin limits
                sel    = rs <= args.rmax
                rs     =     rs[sel]
                dphis  =  dphis[sel]
                thpos  =  thpos[sel]
                thdirs = thdirs[sel]
                qs     =     qs[sel]

                # get bin indexes
                rid     = np.digitize(    rs, rbins)     - 1
                dphiid  = np.digitize( dphis, dphibins)  - 1
                thposid = np.digitize( thpos, thposbins) - 1
                thdirid = np.digitize(thdirs, thdirbins) - 1

                # fill charges into array
                for (ri, dphii, thposi, thdiri, q) in zip(rid, dphiid, thposid, thdirid, qs):
                    # ignore filled index
                    if stats[ri, dphii, thposi, thdiri] == args.Nmax:
                        continue
                    charges[ri, dphii, thposi, thdiri, stats[ri, dphii, thposi, thdiri]] = q
                    # update qid
                    stats  [ri, dphii, thposi, thdiri] += 1

        except KeyboardInterrupt:
            logger.warning("Proccesing interrupted by user")

        # write the charges
        charges.flush()
        stats  .flush()
    return



if __name__ == "__main__":
    main()
