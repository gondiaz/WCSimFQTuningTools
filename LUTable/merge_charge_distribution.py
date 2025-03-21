"""
    tool for temporal studies yet, notice it depends on Pyfitqun
"""
import sys
import argparse

import multiprocessing
import concurrent.futures

from os        import cpu_count
from os.path   import expandvars, basename, join
from itertools import groupby, product
from functools import partial

import tables as tb
import pandas as pd
import numpy  as np

from fitqun.logger import get_logger

logger = get_logger(__name__)


def proccess_index(flatindex, indices, infiles, outfilename, lock):

    logger.info(f"Processing bin {flatindex+1}")

    # loop on files to get the charges
    qs = []
    for filename in infiles:
        with tb.open_file(filename) as f:
            n   = f.root.stats  [indices]
            if n == 0: continue
            qs_ = f.root.charges[indices]
            qs_ = qs_[:n]
            qs.extend(qs_)

    # save to outfile
    with lock:
        with tb.open_file(outfilename, mode='a') as f:

            logger.info("Saving")

            charges_vlarray = f.root.charges
            indices_earray  = f.root.indices

            charges_vlarray.append(qs)
            indices_earray .append([flatindex])

            charges_vlarray.flush()
            indices_earray .flush()
    return


def main():
    """
    TODO: describe
    """

    ############ Program arguments ############
    parser = argparse.ArgumentParser( prog        = f"{basename(__file__)}"
                                    , description = ""
                                    , epilog      = """""")

    parser.add_argument("-v", "--verbose", action="store_true")
    parser.add_argument("-i", "--infiles", type=str, nargs="+"
                       , help="Defines the input files", required=True)
    parser.add_argument("-o", "--outfile", type=str
                        , help = "name of output file", required=True)

    args = parser.parse_args()
    ##########################################

    # define input variables
    infiles     = sorted([expandvars(f) for f in args.infiles])
    outfilename = expandvars(args.outfile)

    # get bins from first file
    with tb.open_file(infiles[0]) as f:
        rbins     = f.root.rbins    .read()
        dphibins  = f.root.dphibins .read()
        thposbins = f.root.thposbins.read()
        thdirbins = f.root.thdirbins.read()

    nrows = (len(rbins) - 1) * (len(dphibins) - 1) * \
            (len(thposbins) - 1) * (len(thdirbins) - 1)

    # create output file
    with tb.open_file(outfilename, mode='w') as f:
        f.create_array  (f.root,     "rbins", rbins)
        f.create_array  (f.root,  "dphibins", dphibins)
        f.create_array  (f.root, "thposbins", thposbins)
        f.create_array  (f.root, "thdirbins", thdirbins)
        f.create_earray (f.root,   "indices", tb.Int32Atom(), shape=(0,), expectedrows=nrows)
        f.create_vlarray(f.root,   "charges", tb.Float32Atom(), " ", tb.Filters(1), nrows)


    # paralellize loop over bins
    with multiprocessing.Manager() as manager:
        lock = manager.Lock()

        logger.info(f"Number of bins to process {nrows}")

        with concurrent.futures.ProcessPoolExecutor(max_workers=cpu_count()) as executor:

            for flatindex, indices in enumerate(product(range(len(rbins)    -1)
                                                      , range(len(dphibins) -1)
                                                      , range(len(thposbins)-1)
                                                      , range(len(thdirbins)-1)), 0):
                executor.submit(proccess_index, flatindex, indices, infiles, outfilename, lock)

                # proccess_index(flatindex, indices, infiles, outfilename, lock)

                # if flatindex == 0: break

    return



if __name__ == "__main__":
    main()
