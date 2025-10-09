from typing import Tuple

import numpy as np
import sg2  # type: ignore
from scipy import interpolate  # type: ignore


def night_filter(
    ssrd: np.ndarray,
    t2m: np.ndarray,
    time_: np.ndarray,
    dt_orig: int,
    dt_downscale: int,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Filter nighttime values except first before/after sunrise/sunset.

    Parameters
    ----------
    ssrd : np.ndarray
        Surface solar radiation downwelling.
    t2m : np.ndarray
        Air temperature at 2-m height.
    time_ : np.ndarray
        Timeseries timestamps.
    dt_orig : int
        Data's original time resolution, in minutes.
    dt_downscale : int
        Time resolution to which data is downscaled, in minutes.

    Returns
    -------
    ix_day : np.ndarray
        index for only daytime values
        (and value right before/after sunrise/sunset).
    ssrd_day : np.ndarray
        ssrd variable with only daytime values
        (and value right before/after sunrise/sunset).
    t2m_day : np.ndarray
        t2m variable with only daytime values
        (and value right before/after sunrise/sunset).
    time_day : np.ndarray
        time_ variable with only daytime values
        (and value right before/after sunrise/sunset).
    time_day_ds : np.ndarray
        Downscaled daytime timestamps for later data downscaling.

    """
    # index for left-side zeros each day
    l_ix = np.where((ssrd[1:, 0] > 0) & (ssrd[:-1, 0] == 0))[0]
    # index for right-side zeros each day
    r_ix = np.where((ssrd[:-1, 0] > 0) & (ssrd[1:, 0] == 0))[0] + 1

    if ssrd[0, 0] > 0:
        l_ix = np.insert(l_ix, 0, 0)
    if ssrd[-1, 0] > 0:
        r_ix = np.append(r_ix, ssrd.shape[0] - 1)

    # indices for only daytime values
    n_days = range(len(l_ix))
    ix_day = np.concatenate([list(range(l_ix[i], r_ix[i] + 1)) for i in n_days])

    # time for downscaled sub-hourly data
    if len(ix_day) == ssrd.shape[0]:  # if there is full daylight month
        time_day_ds = np.arange(
            time_[0] + np.timedelta64(dt_downscale - dt_orig, "m"),
            time_[-1] + np.timedelta64(dt_downscale, "m"),
            np.timedelta64(dt_downscale, "m"),
        )
    else:
        delta_t = dt_downscale - dt_orig
        time_day_ds = np.concatenate(
            [
                np.arange(
                    time_[l_ix[i]] + np.timedelta64(delta_t, "m"),
                    time_[r_ix[i]] + np.timedelta64(dt_downscale, "m"),
                    np.timedelta64(dt_downscale, "m"),
                )
                for i in range(len(l_ix))
            ]
        )

    # night filtering
    ssrd_day = ssrd[ix_day, :]
    t2m_day = t2m[ix_day, :]
    time_day = time_[ix_day]

    return (
        ix_day,
        ssrd_day,
        t2m_day,
        time_day,
        time_day_ds,
    )


def astro_calc(
    vLon: np.ndarray,
    vLat: np.ndarray,
    time_: np.ndarray,
) -> dict:
    """
    Calculate solar astronomy variables.

    Parameters
    ----------
    vLon : np.ndarray
        Latitude for which variables are to be calculated.
    vLat : np.ndarray
        Longitude for which variables are to be calculated.
    time_: np.ndarray
        Timestamps for which variables are to be calculated.

    Returns
    -------
    astro_out : dict
        Set of solar astronomy variables, with keys "SAZ" (solar azimuth),
        "SEL" (solar elevation), and "SZA" (solar zenith),
        "TOA_h" (top of atmosphere irradiance for horizontal surface)

    """
    altitude = 0
    gp = np.array([[vLon, vLat, altitude]])  # N by 3 (in this case, N=1)

    out_list = ["topoc.alpha_S", "topoc.gamma_S0", "topoc.toa_hi"]

    p = sg2.sun_position(gp, time_, out_list)

    astro_out = {}
    # Solar Azimuth
    astro_out["SAZ"] = np.degrees(p.topoc.alpha_S)
    # Solar Elevation
    astro_out["SEL"] = np.degrees(p.topoc.gamma_S0)
    # Solar Zenith
    astro_out["SZA"] = 90 - astro_out["SEL"]

    # Top of Atmosphere horizontal
    astro_out["TOA_h"] = p.topoc.toa_hi
    # clip to non-negative values
    astro_out["TOA_h"] = np.where(
        astro_out["TOA_h"] < 0,
        0,
        astro_out["TOA_h"],
    )

    for key in astro_out.keys():
        astro_out[key] = astro_out[key].reshape(-1, 1)

    return astro_out


def downscale_t2m(
    t2m: np.ndarray,
    dt_orig: int,
    dt_downscale: int,
) -> np.ndarray:
    """
    Downscale air temperature linear interpolation.

    Parameters
    ----------
    t2m: np.ndarray
        Air temperature at 2-m height.
    dt_orig : int
        Original time resolution, in minutes.
    dt_downscale : int
        Time resolution to which data is to be downscaled, in minutes.

    Returns
    -------
    t2m_ds : np.ndarray
        Downscaled air temperature at 2-m height.

    """
    # downscaling factor
    ds_f = dt_orig // dt_downscale

    # downscaling by linear interpolation
    # t=0 needs to be extrapoled to the previous four 15-min steps.
    # for those, repeat T2M(t=0)
    f = interpolate.interp1d(
        np.arange(0, t2m.shape[0]),
        t2m,
        axis=0,
        fill_value=(t2m[0, :], np.nan),
        bounds_error=False,
    )

    t2m_ds = f(np.arange(-1 + 1 / ds_f, t2m.shape[0] - 1 + 1 / ds_f, 1 / ds_f))

    return t2m_ds
