import pandas as pd
import xarray as xr
import numpy as np
from xarray_einstats import linalg


import pyspedas
import pytplot
from pytplot import (
    tplot,
    get_data,
    get_timespan,
    data_exists,
    tplot_copy,
    options,
    split_vec,
    join_vec,
    time_clip
)

from pyspedas.cotrans.minvar_matrix_make import minvar_matrix_make
from pyspedas import tvector_rotate

import logging
from typing import List, Optional

def data_exists_tr(tvar: str, tstart: str, tstop: str) -> bool:
    """
    Checks if a tplot variable exists and if the time range is within the specified time range.
    """
    if tvar in pytplot.data_quants.keys():
        _tstart = pytplot.data_quants[tvar].attrs['plot_options']['trange'][0]
        _tstop = pytplot.data_quants[tvar].attrs['plot_options']['trange'][1]

        tstart_timestamp = pd.Timestamp(tstart)
        tstop_timestamp = pd.Timestamp(tstop)
        _tstart_timestamp = pd.Timestamp(_tstart, unit='s')
        _tstop_timestamp = pd.Timestamp(_tstop, unit='s')

        if tstart_timestamp > _tstart_timestamp and tstop_timestamp < _tstop_timestamp:
            return True
    return False

def fetch_fgm_data(probe: str, trange: List[str], datatype: str) -> None:
    thx_fgm = f"th{probe}_{datatype}_gsm"
    thx_fgm_btotal = f"th{probe}_{datatype}_btotal"
    
    if not data_exists_tr(thx_fgm_btotal, *trange):
        pyspedas.themis.fgm(probe=probe, trange=trange, varnames=[thx_fgm, thx_fgm_btotal])

THX_MAG_DATA_TYPES = ["fgh", "fgl", "fgs"]

def thx_mag_prod(probe: str, tstart: str, tstop: str, datatype: Optional[str] = None, 
                 coord: str = "mva", xarray: bool = True):
    """
    This function loads the magnetic field data for a given probe and time range, and converts it to MVA coordinates.
    If the datatype is not specified, it will try to load the highest resolution data available.
    """
    trange = [tstart, tstop]
    
    if datatype:
        fetch_fgm_data(probe, trange, datatype)
        if not data_exists_tr(f"th{probe}_{datatype}_btotal", *trange):
            logging.warning(f"Data for {probe} and {datatype} not available for the specified time range.")
            return None
    else:
        for dt in THX_MAG_DATA_TYPES:
            fetch_fgm_data(probe, trange, dt)
            if data_exists_tr(f"th{probe}_{dt}_btotal", *trange):
                datatype = dt
                break
    
    # Clip magnetic field data
    tplot_copy(f"th{probe}_{datatype}_gsm", "thx_fgm")
    tplot_copy(f"th{probe}_{datatype}_btotal", "thx_fgm_btotal")
    time_clip("thx_fgm", tstart, tstop, suffix="")
    time_clip("thx_fgm_btotal", tstart, tstop, suffix="")

    # Convert magnetic field from GSM to MVA coordinates
    if coord == "mva":
        minvar_matrix_make("thx_fgm")
        tvector_rotate("thx_fgm_mva_mat", "thx_fgm")
        split_vec("thx_fgm_rot")
        join_vec(
            [
                "thx_fgm_rot_x",
                "thx_fgm_rot_y",
                "thx_fgm_rot_z",
                "thx_fgm_btotal",
            ],
            new_tvar="thx_fgm_all",
        )
        options("thx_fgm_all", "ytitle", f"TH{probe.upper()} {datatype.upper()}")
        options("thx_fgm_all", "color", ["blue", "green", "red", "black"])
        options("thx_fgm_all", "legend_names", [r"$B_l$", r"$B_m$", r"$B_n$", r"$B_{total}$"])
        options("thx_fgm_all", "ysubtitle", "[nT LMN]")

    if xarray:
        return get_data("thx_fgm_all", xarray=True)
    else:
        return get_data("thx_fgm_all")



def thx_part_prod(probe: str, tstart: str, tstop: str, datatype: str = "peir", coord: str = "mva", xarray: bool = True):
    """
    This function loads the particle data for a given probe and time range.
    """
    trange = [tstart, tstop]
    
    # ion velocity
    thx_part = f"th{probe}_{datatype}_velocity_gsm"
    
    if not data_exists_tr(thx_part, *trange):
        pyspedas.themis.esa(probe=probe, trange=trange, datatype=datatype, varnames=[thx_part])
        if not data_exists_tr(thx_part, *trange):
            logging.warning(f"Data for {probe} and {datatype} not available for the specified time range.")
            return None
    
    # Clip particle data
    tplot_copy(thx_part, "thx_part")
    time_clip("thx_part", tstart, tstop, suffix="")
    
    if coord == "mva":
        tvector_rotate("thx_fgm_mva_mat", "thx_part")

    if xarray:
        return get_data("thx_part", xarray=True)
    else:
        return get_data("thx_part")

#%%

def PVI_map(vec, tau_range, resample_frequency=None):
    """_summary_

    Args:
        vec (_type_): _description_
    """
    if resample_frequency = None:
        PVI_series = xr.concat([calculate_PVI_xr(vec, tau) for tau in tau_range], dim='tau')
    else:
        PVI_series = xr.concat([calculate_PVI_xr(vec, tau, resample_frequency) for tau in tau_range], dim='tau')
    return PVI_series
        
        
    

def calculate_mag_change_xr(vec, tau):
    
    # Calculate the magnitudes of the vectors
    mag = linalg.norm(vec, dims='v_dim')
    # Group the vector at the given time lag (tau)
    vec_groups = mag.resample(time=tau)

    mag_change = (vec_groups.max() - vec_groups.min()) / vec_groups.mean()
    
    if 'units' in mag_change.attrs:
        del mag_change.attrs['units']

    return mag_change.rename('Mag change')


def calculate_PVI_series(vec, tau, interval_of_averaging=None):
    """
    This function calculates the Partial Variance of Increments (PVI) series for a given time series.

    Parameters:
    vec (np.array): The input time series, represented as a numpy array.
    tau (int): The time lag, typically selected to lie in the inertial range of the fluctuations.
    interval_of_averaging (int): The number of samples over which to compute the trailing average. It's often chosen to be comparable to, or greater than, a correlation length (or time) of the signal.

    Returns:
    PVI_series (np.array): The resulting PVI series.
    """
    
    # Calculate the increments in the vector for the given time lag (tau)
    increments = vec[tau:] - vec[:-tau]
    
    # Calculate the magnitudes of these increments
    mag_increments = np.linalg.norm(increments, axis=1)
    
    # Square the magnitudes of the increments, and compute a moving average over the specified interval
    if interval_of_averaging is None:
        normalized_factor = np.sqrt(np.mean(np.square(mag_increments)))
        PVI_series = mag_increments / normalized_factor
    else:
        avg_square_increments = np.convolve(np.square(mag_increments[1:]), np.ones(interval_of_averaging), 'valid') / interval_of_averaging

        # Divide the magnitude of the increments by the square root of the averaged squared increments
        # This gives the PVI series
        PVI_series = mag_increments[:-interval_of_averaging] / np.sqrt(avg_square_increments)

    return PVI_series


def calculate_mag_change(vec, tau, interval_of_averaging):
    
    # Calculate the magnitude of vec
    mag = np.linalg.norm(vec, axis=1)
    mag_s_tau = mag[tau:] - mag[:-tau]

    # Moving average
    mag_avg = np.convolve(np.square(mag[1:-tau]), np.ones(interval_of_averaging), 'valid') / interval_of_averaging

    return mag_s_tau[:-interval_of_averaging] / np.sqrt(mag_avg)

# vec = np.random.rand(1000, 3)
# tau = 10
# interval_of_averaging = 100
# PVI_series = calculate_PVI_series(vec, tau, interval_of_averaging)
# calculate_mag_change(vec, tau, interval_of_averaging)
# np.mean(np.square(PVI_series))

#%%
def themis_image_url(date, probe, type):
    """
    Returns the URL for a THEMIS image for a given date, probe, and type.

    Examples:
    http://themis.ssl.berkeley.edu/themisdata/overplots/2019/01/01/thc_l2_moms_20190101_0024.png
    http://themis.ssl.berkeley.edu/themisdata/overplots/2020/02/20/thc_l2_overview_20200220_0024.png
    http://themis.ssl.berkeley.edu/themisdata/thg/l0/asi/2020/02/orbit_moon_multi_mission_2020-02-20_0024.gif
    http://themis.ssl.berkeley.edu/themisdata/thg/l0/asi/2020/02/orbit_multi_mission_2020-02-20_0024.gif
    """
    
    base_url = f"http://themis.ssl.berkeley.edu/themisdata"
    match type:
        case 'overview_moms':
            return f"{base_url}/overplots/{date.strftime('%Y/%m/%d')}/th{probe}_l2_moms_{date.strftime('%Y%m%d')}_0024.png"
        case 'orbit_moon':
            return f"{base_url}/thg/l0/asi/{date.strftime('%Y/%m')}/orbit_moon_multi_mission_{date.strftime('%Y-%m-%d')}_0024.gif"