from typing import Tuple

import numpy as np
import pvlib as pvlib  # type: ignore


def ssrd_downscale(
    ssrd: np.ndarray,
    TOA_h: np.ndarray,
    dt_orig: int,
    dt_downscale: int,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Downscale surface solar radiation time-wise.

    Parameters
    ----------
    ssrd : np.ndarray
        Surface solar radiation downwelling.
    TOA_h : np.ndarray
        Top of Atmosphere horizontal irradiance, before impact of atmosphere.
    dt_orig : int
        Original time resolution, in minutes.
    dt_downscale : int
        Finer resolution to which data is downscaled.

    Returns
    -------
    ssrd_ds : np.ndarray
        Downscaled surface solar radiation downwelling.
    Kt : np.ndarray
        Clearness index, which is used when processing ssrd later on.

    """
    # downscaling factor
    n_reps = dt_orig // dt_downscale

    # reshape to easily calculate hourly averages
    TOA_h_hourly = TOA_h.reshape(-1, n_reps).mean(axis=1).reshape(-1, 1)

    TOA_h_hourly[TOA_h_hourly == 0] = np.nan
    Kt = ssrd / TOA_h_hourly
    Kt = np.where(Kt > 1, 1, Kt)
    Kt[np.isinf(Kt)] = np.nan

    # Kt is assumed to be constant in each hour
    # TOA_h captures the daily cycle at higher resolution
    Kt = np.repeat(Kt, 4, axis=0)
    ssrd_ds = Kt * TOA_h

    ssrd_ds[np.isnan(ssrd_ds)] = 0

    # makes later calculations + efficient, nan will correspond to nightime
    Kt[ssrd_ds == 0] = np.nan

    return ssrd_ds, Kt


def ssrd_decompose(ssrd: np.ndarray, Kt: np.ndarray, SEL: np.ndarray) -> dict:
    """
    Estimate diffuse and direct irradiance components from global irradiance.

    Parameters
    ----------
    ssrd : np.ndarray
        Surface solar radiation downwelling.
    Kt : np.ndarray
        Clearness index.
    SEL : np.ndarray
        Solar elevation angle.

    Returns
    -------
    decomp_out : dict
        Irradiance componentes, with keys "DNI" (direct normal irradiance) and
        "DHI" (diffuse horizontal irradiance).

    """
    # diffuse fraction
    Kd = np.zeros(ssrd.shape)

    # Skarveit Olseth (1987), https://doi.org/10.1016/0038-092X(87)90049-1
    k_c = 0.87 - 0.56 * np.exp(-0.06 * SEL)
    d_c = 0.15 + 0.43 * np.exp(-0.06 * SEL)
    k_0 = 0.20
    a_1 = 1.09
    a_2 = 0.27
    K = 0.5 * (1 + np.sin(np.pi * ((Kt - k_0) / (k_c - k_0) - 0.5)))

    # three regimes
    ix1 = Kt <= k_0
    ix2 = (Kt > k_0) & (Kt <= a_1 * k_c)
    ix3 = Kt > a_1 * k_c

    Kd[ix1] = 1

    temp1 = 1 - d_c[ix2]
    temp2 = a_2 * np.sqrt(K[ix2])
    temp3 = (1 - a_2) * K[ix2] ** 2

    Kd[ix2] = 1 - temp1 * (temp2 + temp3)
    del temp1, temp2, temp3

    temp1 = a_1 * k_c[ix3] - k_0
    temp2 = k_c[ix3] - k_0  # estava aqui o erro
    K_3 = 0.5 * (1 + np.sin(np.pi * (temp1 / temp2 - 0.5)))

    F_K = 1 - (1 - d_c[ix3]) * (a_2 * np.sqrt(K_3) + (1 - a_2) * K_3**2)
    Kd[ix3] = 1 - a_1 * k_c[ix3] * (1 - F_K) / (0.0001 + Kt[ix3])

    # clips diffuse fraction between 0 and 1
    Kd = np.clip(Kd, 0, 1)

    # diffuse horizontal irradiance
    dhi = ssrd * Kd
    dhi = np.minimum(dhi, ssrd)

    # direct normal irradiance
    dni = np.zeros(ssrd.shape)
    ix = SEL > 1  # daytime
    dni[ix] = (ssrd[ix] - dhi[ix]) / np.sin(np.radians(SEL[ix]))
    dni = np.where(dni > 1370, 1370, dni)

    decomp_out = {"DHI": dhi, "DNI": dni}

    return decomp_out


def irrad_transpose(
    pv_tilt: np.ndarray,
    pv_azim: np.ndarray,
    ssrd: np.ndarray,
    decomp_out: dict,
    astro_out: dict,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Transpose irradiance components to given tilted plane.

    Parameters
    ----------
    pv_tilt : np.ndarray
        PV module tilt, in degrees.
    pv_azim : np.ndarray
        PV module azimuth, in degrees.
    ssrd : np.ndarray
        Surface solar radiation downwelling.
    decomp_out : dict
        Solar radiation components.
    astro_out : dict
        Solar astronomy variables, such as solar elevation and azimuth.

    Returns
    -------
    POA_aoi : np.ndarray
        Angle of incidence between module surface and sun.
    POA_dir : np.ndarray
        Direct irradiance at module plane-of-array.
    POA_dif : np.ndarray
        Diffuse irradiance at module plane-of-array.
    POA_ref : np.ndarray
        Reflected irradiance at module plane-of-array.

    """
    # AOI: angle of incidence between sun and given surface
    # POA: plane-of array
    POA_aoi = pvlib.irradiance.aoi(
        pv_tilt,
        pv_azim,
        astro_out["SZA"],
        astro_out["SAZ"],
    )

    # Direct tilted irradiance
    POA_dir = pvlib.irradiance.beam_component(
        pv_tilt, pv_azim, astro_out["SZA"], astro_out["SAZ"], decomp_out["DNI"]
    )

    # Diffuse Tilted Irradiance
    # Klucher (1979), https://doi.org/10.1016/0038-092X(79)90110-5
    POA_dif = pvlib.irradiance.klucher(
        pv_tilt,
        pv_azim,
        decomp_out["DHI"],
        ssrd,
        astro_out["SZA"],
        astro_out["SAZ"],
    )

    # ground albedo
    ALB = 0.2
    # Reflected Tilted Irradiance
    POA_ref = pvlib.irradiance.get_ground_diffuse(pv_tilt, ssrd, albedo=ALB)

    return POA_aoi, POA_dir, POA_dif, POA_ref
