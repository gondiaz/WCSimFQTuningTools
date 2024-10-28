import ROOT
import uproot
import itertools
import numpy  as np
import pandas as pd

from os.path import expandvars

from numba import jit, vectorize, float64

@vectorize([float64(float64, float64)])
def azimuth_angle(x, y):
    """
    Computes the azimuth angle [0, 2pi) for a vector of components x, y
    """
    mod = np.sqrt(x**2 + y**2)
    if (mod == 0.): return 0.
    phi = np.arccos(x/mod)
    if (y != 0): phi  = np.sign(y)*phi
    if (phi< 0): phi += 2.*np.pi
    return phi

@jit(cache=True)
def clockwise_azimuth_angle(v1, v2):
    """Returns clockwise angle between the two vectors
    TODO: improve description"""
    angles = np.zeros(len(v1))
    a1 = azimuth_angle(v1[:, 0], v1[:, 1])
    a2 = azimuth_angle(v2[:, 0], v2[:, 1])
    
    da = a2 - a1
    sel = da>0
    angles[sel] = 2.*np.pi - da[sel]
    sel = da<0
    angles[sel] = -da[sel]
    return angles


def read_wcsim_geometry(filename):

    """
    This function reads the geometrical information (in cm)
    from **filename** and returns it in two dataframes.
    The first dataframe contains general information about the detector,
    and the second returns the information for each PMT

    Requires WCSimRoot library loaded
    """
    filename = expandvars(filename)
    rootf = ROOT.TFile(filename, "read")

    tree  = rootf.GetKey("wcsimGeoT").ReadObj()
    tree.GetEvent(0)
    geom  = tree.wcsimrootgeom

    # general info
    df = pd.DataFrame()
    df = pd.DataFrame(columns=["WC"])
    df.loc["WCCylRadius"] = geom.GetWCCylRadius()
    df.loc["WCCylLength"] = geom.GetWCCylLength()
    df.loc["Geo_Type"]    = geom.GetGeo_Type   ()
    df.loc["WCNumPMT"]    = geom.GetWCNumPMT   ()
    df.loc["WCPMTRadius"] = geom.GetWCPMTRadius()
    for i in range(3): df.loc[f"WCOffset{i}"] = geom.GetWCOffset(i)

    with uproot.open(filename) as f:
        for key in ["WCDetRadius", "WCDetHeight"]:
            df.loc[key] = 0.1*f[f"Settings/{key}"].array()[0]

    # mPMT info
    columns = ["TubeNo", "mPMTNo", "mPMT_PMTNo", "CylLoc"]
    for i in range(3): columns.append(f"Orientation_x{i}")
    for i in range(3): columns.append(f"Position_x{i}")

    pmts_empty = np.zeros(geom.GetWCNumPMT(), dtype=list(zip(columns, [*4*[int], *6*[float]])))
    pmts_df = pd.DataFrame(pmts_empty)

    for i in range(geom.GetWCNumPMT()):
        pmt = geom.GetPMT(i)
        pmts_df.loc[i, "TubeNo"]     = int(pmt.GetTubeNo())
        pmts_df.loc[i, "mPMTNo"]     = int(pmt.GetmPMTNo())
        pmts_df.loc[i, "mPMT_PMTNo"] = int(pmt.GetmPMT_PMTNo())
        pmts_df.loc[i, "CylLoc"]     = int(pmt.GetCylLoc())

        for j in range(3):
            pmts_df.loc[i, f"Orientation_x{j}"] = float(pmt.GetOrientation(j))
            pmts_df.loc[i, f"Position_x{j}"]    = float(pmt.GetPosition(j))
    
    rootf.Close()
    return (df, pmts_df)


def split_tubeids(filename, vaxis=2):

    """ 
    This function reads the geometrical information from the file
    and returns the tubeids for bottom, top and side tubes in this order.

    CAUTION:
    This function assumes that the central PMT orientation is perpendicular to the vessel
    and its internal number in the mPMT module is 1 (ie mPMT_PMTNo == 1).

    This assumption might not be correct in future software versions or detector geometries.

    The ideal implementation would be to perform the splitting using a location flag for each PMT.
    """
    _, pmts_df = read_wcsim_geometry(expandvars(filename))

    mPMTids       = pmts_df.mPMTNo.unique()
    mPMTid_bottom = pmts_df.loc[(pmts_df.mPMT_PMTNo == 1) & (pmts_df.loc[:, f"Orientation_x{vaxis}"] == +1)].mPMTNo.values
    mPMTid_top    = pmts_df.loc[(pmts_df.mPMT_PMTNo == 1) & (pmts_df.loc[:, f"Orientation_x{vaxis}"] == -1)].mPMTNo.values
    mPMTid_side   = mPMTids[~np.isin(mPMTids, np.concatenate([mPMTid_bottom, mPMTid_top]))]

    tubeid_bottom = pmts_df.set_index("mPMTNo").loc[mPMTid_bottom].TubeNo.values
    tubeid_top    = pmts_df.set_index("mPMTNo").loc[mPMTid_top   ].TubeNo.values
    tubeid_side   = pmts_df.set_index("mPMTNo").loc[mPMTid_side  ].TubeNo.values

    assert len(tubeid_side) + len(tubeid_bottom) + len(tubeid_top) == len(pmts_df)

    return tubeid_bottom, tubeid_top, tubeid_side


def read_stable(filename, names=["botscattable", "topscattable", "sidescattable"]):
    """
    Read STables from file as numpy arrays
    the output is a dictionary with the bins and the tables 
    for the bottom, top and side PMTs

    Requires fiTQun library loaded
    """
    out = dict()

    file = ROOT.TFile(expandvars(filename))
    for tabname in names:
        tab = file.Get(tabname)

        nbins  = [tab.GetNBins(dim) for dim in range(6)]
        bounds = [[tab.GetAxisBound(dim, i) for i in range(2)] for dim in range(6)]

        histo6d = np.ndarray(shape=nbins)
        for index in itertools.product(*[range(n) for n in nbins]): histo6d[index] = tab.GetElement(*index)
        
        out[tabname + "_bins"] = (nbins, bounds)
        out[tabname]           = histo6d
    return out