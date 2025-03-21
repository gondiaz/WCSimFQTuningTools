"""
    tool for temporal studies yet, notice it depends on Pyfitqun
"""
import sys
import glob
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


def checker(filename):

    try:
        with tb.open(filename) as f:
            pass
    except:
        logger.critical(f"CORRUPT {basename(filename)}")

   
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
    parser.add_argument("-i", "--indir", type=str, nargs=1
                       , help="Defines the input directory", required=True)

    args = parser.parse_args()
    #########################################

    infiles = sorted(glob.glob(join(args.indir[0], "*")))

    # paralellize
    with multiprocessing.Manager() as manager:
        lock = manager.Lock()

        with concurrent.futures.ProcessPoolExecutor(max_workers=cpu_count()) as executor:

            for filename in infiles:
                executor.submit(checker, filename)

    return



if __name__ == "__main__":
    main()
