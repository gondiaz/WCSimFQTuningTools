import ROOT
import itertools
import numpy  as np
import pandas as pd

from os.path   import expandvars

def read_wcsim_geometry(filename):

    """
    This function reads the geometrical information (in mm)
    from **filename** and returns it in two dataframes.
    The first dataframe contains general information about the detector,
    and the second returns the information for each PMT
    """
    rootf = ROOT.TFile(expandvars(filename), "read")

    tree  = rootf.GetKey("wcsimGeoT").ReadObj()
    tree.GetEvent(0)
    geom  = tree.wcsimrootgeom

    # general info
    df = pd.DataFrame()
    df = pd.DataFrame(columns=["WC"])
    df.loc["WCCylRadius"] = geom.GetWCCylRadius() * 10
    df.loc["WCCylLength"] = geom.GetWCCylLength() * 10
    df.loc["Geo_Type"]    = geom.GetGeo_Type   ()
    df.loc["WCNumPMT"]    = geom.GetWCNumPMT   ()
    df.loc["WCPMTRadius"] = geom.GetWCPMTRadius()
    for i in range(3): df.loc[f"WCOffset{i}"] = geom.GetWCOffset(i)

    # mPMT info
    columns = ["TubeNo", "mPMTNo", "mPMT_PMTNo", "CylLoc"]
    for i in range(3): columns.append(f"Orientation_x{i}")
    for i in range(3): columns.append(f"Position_x{i}")
    pmts_df = pd.DataFrame(columns=columns)

    for i in range(geom.GetWCNumPMT()):
        pmt = geom.GetPMT(i)
        pmts_df.loc[i, "TubeNo"]     = pmt.GetTubeNo()
        pmts_df.loc[i, "mPMTNo"]     = pmt.GetmPMTNo()
        pmts_df.loc[i, "mPMT_PMTNo"] = pmt.GetmPMT_PMTNo()
        pmts_df.loc[i, "CylLoc"]     = pmt.GetCylLoc()

        for j in range(3):
            pmts_df.loc[i, f"Orientation_x{j}"] = pmt.GetOrientation(j)
            pmts_df.loc[i, f"Position_x{j}"]    = pmt.GetPosition(j)*10
    
    rootf.Close()
    return (df, pmts_df)


def split_tubeids(filename, vaxis=2):

    """ 
    This function reads the geometrical information from the file
    and returns the tubeids for bottom, top and side tubes in this order
    """
    _, pmts_df = read_wcsim_geometry(expandvars(filename))

    # select bottom, top and side pmtids
    zm = getattr(pmts_df.groupby("mPMTNo"), f"Position_x{vaxis}").mean()
    mPMT_bottom = zm[zm == zm.min()].index.values
    mPMT_top    = zm[zm == zm.max()].index.values

    tubeid_bottom = pmts_df[np.isin(pmts_df.mPMTNo, mPMT_bottom)].TubeNo.values
    tubeid_top    = pmts_df[np.isin(pmts_df.mPMTNo, mPMT_top)]   .TubeNo.values

    # assuming same number of top and bottom mPMTs and same number of PMTs in each module
    assert len(tubeid_bottom) == len(tubeid_top)

    # side pmts are those not in top or bottom
    tubeid_side = pmts_df.TubeNo.values
    tubeid_side = tubeid_side[~(np.isin(tubeid_side, tubeid_bottom) | np.isin(tubeid_side, tubeid_top))]

    assert len(tubeid_side) + len(tubeid_bottom) + len(tubeid_top) == len(pmts_df)

    return tubeid_bottom, tubeid_top, tubeid_side


def azimuth_angle(v):
    """
    Computes the azimuth angle [0, 2pi) in radians for the vector **v**
    Assumes that first and second entries are the X and Y projections
    **v** can be an array of vectors, ie an array of shape (# of vectors, vector dimension)
    """
    # between [-pi, pi]
    m   = np.linalg.norm(v[:, (0, 1)], axis=1)
    px  = v[:, 0] / m
    phi = np.arccos(px)
    sel = (v[:, 1] != 0)
    phi[sel] = np.sign(v[sel, 1]) * phi[sel]
    # transform to [0, 2pi)
    phi[phi<0] += 2.*np.pi
    return phi


def clockwise_angle(v1, v2):
    """Returns clockwise angle between the two vectors"""
    angles = np.zeros(len(v1))
    a1 = azimuth_angle(v1)
    a2 = azimuth_angle(v2)
    
    da = a2 - a1
    sel = da>0
    angles[ sel] = 2.*np.pi - da[sel]
    sel = da<0
    angles[sel] = -da[sel]
    return angles


def read_stable(filename, names=["botscattable", "topscattable", "sidescattable"]):
    """
    Read STables from file as numpy arrays
    the output is a dictionary with the bins and the tables 
    for the bottom, top and side PMTs
    """
    # TO-DO: Load TScatTable_c.so library inside function
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