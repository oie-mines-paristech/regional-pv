from typing import Optional, Union

import numpy as np
import pvlib  # type: ignore

from regional_pv.utils.pv_functions import (
    account_electric_losses,
    account_optical_losses,
    account_stc_efficiency,
    account_thermal_losses,
    compute_Tmodule,
)
from regional_pv.utils.radiation_functions import (
    irrad_transpose,
    ssrd_decompose,
    ssrd_downscale,
)
from regional_pv.utils.support_functions import (
    astro_calc,
    downscale_t2m,
    night_filter,
)


def spv_workflow(
    pv_type: str,
    ssrd: np.ndarray,
    t2m: np.ndarray,
    meta: list,
    azim: Optional[np.ndarray],
    tilt: Optional[np.ndarray],
    w_orient: Optional[np.ndarray],
    k: float,
    dt_orig: Union[int, dict] = 60,
    dt_downscale: int = 15,
) -> np.ndarray:
    """
    Modelling chain converting ssrd and t2m variables to PV capacity factor.

    Parameters
    ----------
    pv_type: str
        IDs PV typology being computed.
    ssrd : np.ndarray
        Surface solar radiation downwelling for a given location.
    t2m : np.ndarray
        2-m air temperature for a given location.
    meta : list
        Three elements: timestamp, longitude, and latitude np.ndarrays.
    azim : Optional[np.ndarray]
        Module azimuth to be considered (degrees, 180º: South, 90º East).
        None if tracking typology.
    tilt : Optional[np.ndarray]
        Module tilt to be considered (degrees, 0º horizontal).
        None if tracking typology.
    w_orient: Optional[np.ndarray]
        Weights for combining multiple module orientations.
        None if single orientation or tracking typology.
    k: float
        PV Ross parameter.
    dt_orig : int
        Original time resolution of weather data, in minutes. Can be dictionary
        if ssrd and t2m differ.
        Default value assumes 1 hour (60).
    dt_downscale : int
        Finer resolution to which the data is downscaled to better capture
        intra-day losses (mainly optical and thermal).
        Default value assumes downscaling to 15-min data.

    Returns
    -------
    nPac_full : np.ndarray
        Normalized PV capacity factor, accounting for inverter losses. Assumes
        no overplanting, meaning that it can be used as MW/MWp upon which
        clipping losses are calculated.

    """
    # time and coordinates for the location of interest
    time_ = meta[0]
    lon = meta[1]
    lat = meta[2]

    if isinstance(dt_orig, int):
        dt_orig_ssrd = dt_orig
        dt_orig_t2m = dt_orig
    elif isinstance(dt_orig, dict):
        # checks if dict has proper keys
        msg = "dict keys must contain ssrd and t2m."
        assert {"ssrd", "t2m"}.issubset(dt_orig), msg

        dt_orig_ssrd = dt_orig["ssrd"]
        dt_orig_t2m = dt_orig["t2m"]
    else:
        raise TypeError("dt_orig must be either int or dict.")

    # to shape final output variable
    n_datapoints = ssrd.shape[0]

    if dt_orig_ssrd < 24 * 60:  # for sub-daily data
        # filters night values except last/first sunrise/sunset moment of each day
        # also prepares downscaled night-ignoring time vector
        ix_day, ssrd_daytime, t2m_daytime, time__, time_day_ds = night_filter(
            ssrd, t2m, time_, dt_orig_ssrd, dt_downscale
        )

        # calculates solar position and top of atmosphere (TOA) radiation
        # for downscaled temporal resolution
        astro_out = astro_calc(lon, lat, time_day_ds)

    else:
        time_ds = np.arange(
            time_[0] + np.timedelta64(-24, "h"),
            time_[-1],
            np.timedelta64(15, "m"),
        )

        # calculates solar position and top of atmosphere (TOA) radiation
        # for downscaled temporal resolution
        astro_out = astro_calc(lon, lat, time_ds)

    # downscales ssrd in time by linear interpolation of clearness index (Kt)
    # also outputs downscaled Kt to use later
    ssrd_daytime_ds, tKT = ssrd_downscale(
        ssrd_daytime, astro_out["TOA_h"], dt_orig_ssrd, dt_downscale
    )

    # T2M temporal downscale
    t2m_daytime_ds = downscale_t2m(t2m_daytime, dt_orig_t2m, dt_downscale)

    # once downscaling is done, we can fully ignore nighttime
    ix_day_ds = ssrd_daytime_ds.flatten() > 0
    for key in astro_out.keys():
        astro_out[key] = astro_out[key][ix_day_ds, :]
    ssrd_daytime_ds = ssrd_daytime_ds[ix_day_ds, :]
    tKT = tKT[ix_day_ds, :]
    t2m_daytime_ds = t2m_daytime_ds[ix_day_ds, :]

    # TODO: define in global, read it in main, adapt for 1D, 2D, 3D cases
    if pv_type == "utility_track_1axis":
        out = pvlib.tracking.singleaxis(
            astro_out["SZA"].flatten(),
            astro_out["SAZ"].flatten(),
            axis_tilt=0,
            axis_azimuth=180,
            backtrack=True,
            gcr=0.35,
        )

        # in non-tracking this is global (so not explicit input of function)
        tilt_ = out["surface_tilt"].reshape(-1, 1)
        azim_ = out["surface_azimuth"].reshape(-1, 1)
    else:
        tilt_ = tilt
        azim_ = azim
        # aoi = out['aoi'].reshape(-1,1)

    # prep code for/if when 2-axis tracking is included
    # elif pv_type == 'utility_track_2axis':
    # TODO: check if reshape is needed
    # TODO: check if 2-axis often uses backtracking and operational constraints
    # pv_tilt = astro_out['SZA'].reshape(-1,1)
    # pv_azim = astro_out['SAZ'].reshape(-1,1)
    # aoi = out['aoi'].reshape(-1,1)

    # infers diffuse and direct components from global irradiance
    decomp_out = ssrd_decompose(
        ssrd_daytime_ds,
        tKT,
        astro_out["SEL"],
    )

    # transposes horizontal irradiance components to given tilted plane
    # (or plane-of-array,POA). it also models reflected irradiance
    POA_aoi, POA_dir, POA_dif, POA_ref = irrad_transpose(
        tilt_, azim_, ssrd_daytime_ds, decomp_out, astro_out
    )
    del ssrd, ssrd_daytime, ssrd_daytime_ds, decomp_out, astro_out

    # estimate PV module temperature, T2M must be in ºC
    tmod = compute_Tmodule(t2m_daytime_ds, POA_dir + POA_dif + POA_ref, k)
    del t2m, t2m_daytime, t2m_daytime_ds

    # Calculates effective GTI, after accounting for optical losses
    GTI_eff = account_optical_losses(tilt_, POA_aoi, POA_dif, POA_dir, POA_ref)
    del POA_aoi, POA_dir, POA_dif, POA_ref

    # normalized PV capacity factor, considering STC
    # (standard test conditions)
    nPdc25 = account_stc_efficiency(GTI_eff)
    del GTI_eff

    # normalized effective DC capacity, after accounting for thermal losses
    nPdc = account_thermal_losses(nPdc25, tmod)
    del nPdc25, tmod

    # effective PV capacity factor, accounting for electric losses.
    # still dismisses overplanting and possible clipping
    nPac = account_electric_losses(nPdc)
    del nPdc

    # aggregate different geometries
    if nPac.shape[1] > 1:
        nPac = (nPac * w_orient).sum(axis=1)

    # rebuild timeseries for proper reconversion to hourly
    temp = np.zeros(([ix_day_ds.shape[0], 1]))
    temp[ix_day_ds, 0] = nPac.flatten()
    nPac = temp.reshape(-1, 4).mean(axis=1)

    # full night-including shape
    nPac_full = np.zeros((n_datapoints,))
    nPac_full[ix_day] = nPac

    return nPac_full
