import numpy as np


def compute_Tmodule(
    T2M_degC: np.ndarray,
    GTI: np.ndarray,
    k: float,
) -> np.ndarray:
    """
    Estimate PV module temperature.

    Parameters
    ----------
    T2M_degC : np.ndarray
        2-m high air temperature, must be in degrees Celsius.
    GTI : np.ndarray
        Global tilted irradiance, in W.m-2.
    k: float
        PV thermal (Ross) coefficient.

    Returns
    -------
    Tmod : np.ndarray
          PV module temperature, in degrees Celsius.

    """
    # Calculation of the module temperature
    Tmod = T2M_degC + k / 1000 * GTI

    return Tmod


def account_optical_losses(
    pv_tilt: np.ndarray,
    POA_aoi: np.ndarray,
    POA_dif: np.ndarray,
    POA_dir: np.ndarray,
    POA_ref: np.ndarray,
) -> np.ndarray:
    """
    Estimate effective incident radiation, accounting for optical losses.

    Parameters
    ----------
    pv_tilt : np.ndarray
        Module tilt to be considered (degrees, 0º horizontal).
    POA_aoi : np.ndarray
        Angle of incidence between module and sun (degrees, 0º normal).
    POA_dif : np.ndarray
        Diffuse irradiance over the module (tilted) plane.
    POA_dir : np.ndarray
        Direct irradiance over the module (tilted) plane.
    POA_ref : np.ndarray
        Reflected irradiance over the module (tilted) plane.

    Returns
    -------
    Geff : np.ndarray
        Effective global tilted irradiance.

    """
    # Martin N, Ruiz JM (2001), https://doi.org/10.1016/S0927-0248(00)00408-6
    # Corrigendum (2012), https://doi.org/10.1016/j.solmat.2012.11.002
    # values extracted from Table 3 p.32
    arOpt = 0.17
    c1Opt = 4 / (3 * np.pi)
    c2Opt = -0.069
    exp_arOpt = np.exp(-1 / arOpt)

    # Module tilt angle in radian
    B = np.deg2rad(pv_tilt)
    B = np.where(B < 0.00001, 0.00001, B)

    cos_B = np.cos(B)
    sin_B = np.sin(B)

    aoi_rad = POA_aoi * np.pi / 180

    Fb = 1 - (np.exp(-np.cos(aoi_rad) / arOpt) - exp_arOpt) / (1 - exp_arOpt)

    Fr = 1 - np.exp(
        -c1Opt / arOpt * (sin_B + (B - sin_B) / (1 - cos_B))
        - c2Opt / arOpt * (sin_B + (B - sin_B) / (1 - cos_B)) ** 2
    )

    Fd = 1 - np.exp(
        -c1Opt / arOpt * (sin_B + (np.pi - B - sin_B) / (1 + cos_B))
        - c2Opt / arOpt * (sin_B + (np.pi - B - sin_B) / (1 + cos_B)) ** 2
    )

    # avoid negative coefficients
    Fb = np.where(Fb < 0, 0, Fb)
    Fd = np.where(Fd < 0, 0, Fd)
    Fr = np.where(Fr < 0, 0, Fr)

    # effective incident radiation
    Geff = POA_dir * Fb + POA_dif * Fd + POA_ref * Fr

    return Geff


def account_stc_efficiency(Geff: np.ndarray) -> np.ndarray:
    """
    Estimate PV capacity factor at standard test conditions (STC).

    Parameters
    ----------
    Geff : np.ndarray
        Effective incident radiation.

    Returns
    -------
    nPdc_stc : np.ndarray
        PV capacity factor, considering standard test conditions (STC).

    """
    # source ?
    ai = [1.4306, -1.0084, 1.0121, -0.4401, 0.1979]

    # for the Pdc25 calculation
    Geff = Geff / 1000

    nPdc_stc = np.maximum(
        0,
        ai[0] * Geff
        + ai[1] * Geff**2
        + ai[2] * Geff**3
        + ai[3] * Geff**4
        + ai[4] * Geff * np.log(Geff),
    )

    return nPdc_stc


def account_thermal_losses(
    nPdc_stc: np.ndarray,
    Tmod: np.ndarray,
) -> np.ndarray:
    """
    Account for thermal PV losses on top of optics+STC capacity factor.

    Parameters
    ----------
    nPdc_stc : np.ndarray
        PV capacity factor, considering optical losses and
        standard test conditions (STC).
    Tmod : np.ndarray
        PV module temperature, in degrees Celsius.

    Returns
    -------
    nPdc : np.ndarray
        PV capacity factor, considering optical and thermal losses.

    """
    # generic thermal coefficient (source?)
    Mu_T = -0.45 / 100

    # thermal relative losses
    loss_T = Mu_T * (Tmod - 25)

    # Consideration of the effect of the module temperature on the DC power
    nPdc = np.maximum(0, nPdc_stc * (1 + loss_T))

    return nPdc


# normalized effective AC capacity, after accounting for electric losses
def account_electric_losses(nPdc: np.ndarray) -> np.ndarray:
    """
    Account for inverter losses on top of PV capacity factor.

    Parameters
    ----------
    nPdc : np.ndarray
        PV capacity factor considering optical and thermal losses.

    Returns
    -------
    nPac : np.ndarray
        PV capacity factor, considering optical, thermal, & electric losses.
        It still does not account for overplanting, nor for clipping losses.

    """
    # accounts for inverter efficiency and an additional derating factor
    # (for other losses).

    bi = [-0.0005, 1.0269, -0.00217]
    kDerating = 0.8

    nPac = kDerating * (bi[0] + bi[1] * nPdc + bi[2] * nPdc**2)
    nPac = np.where(nPac < 0, 0, nPac)

    return nPac
