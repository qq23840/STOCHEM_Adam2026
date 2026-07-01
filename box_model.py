"""
Box model for atmospheric chemistry kinetics and Global Warming Potential calculations.

This module implements a 1-box model for simulating the oxidation of organic species
by atmospheric OH and O3, tracking their degradation to CF3CHO and HFC-23. It includes
functions for kinetic integration and GWP calculations over specified time horizons.

Author: StoChem Project
License: (specify your license)
"""

import json
from pathlib import Path
from functools import partial
from typing import Tuple, List, Dict, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.integrate import odeint

# ============================================================================
# Constants
# ============================================================================

# Time conversion factors
SECONDS_PER_YEAR = 60 * 60 * 24 * 365
SECONDS_PER_HOUR = 3600

# Atmospheric concentrations (ppb and per cm^3 conversions)
O3_CONCENTRATION_PPB = 50  # typical tropospheric ozone concentration
PPBV_TO_PER_CM3 = 2.5e10
PPTV_TO_PER_CM3 = 2.5e7
MOLAR_CONCENTRATION_CM3 = 2.5e19  # for conversion factors

# Atmospheric hydroxyl radical concentration
OH_CONCENTRATION_PER_CM3 = 1e6

# Lifetimes (years)
HFC23_LIFETIME_YEARS = 243
CF3CHO_LIFETIME_HOURS = 72

# Molar masses (g/mol)
MOLAR_MASS_HFC143A = 84.4
MOLAR_MASS_CF3CHO = 90.0
MOLAR_MASS_HFC23 = 70.3

# Radiative efficiency and CO2 forcing (WMO 2022)
HFC23_RADIATIVE_EFFICIENCY = 0.192e-3  # W m^-2 ppt^-1
CO2_GWP_FORCING_20_YR = 2.434e-14  # W m^-2 yr kg^-1
CO2_GWP_FORCING_100_YR = 8.947e-14  # W m^-2 yr kg^-1
CO2_GWP_FORCING_500_YR = 3.138e-13  # W m^-2 yr kg^-1

# Time horizons and corresponding parameters
TIME_HORIZONS = {
    20: {'end_date': '1720-01-01', 'gwp_forcing': CO2_GWP_FORCING_20_YR},
    100: {'end_date': '1800-01-01', 'gwp_forcing': CO2_GWP_FORCING_100_YR},
    500: {'end_date': '2200-01-01', 'gwp_forcing': CO2_GWP_FORCING_500_YR},
}

# File paths (relative to module directory)
SPECIES_DATA_FILE = Path(__file__).parent / 'species_data.json'


# ============================================================================
# Utility functions
# ============================================================================

def load_species_data(filepath: Path = SPECIES_DATA_FILE) -> Dict:
    """
    Load species kinetic data from JSON file.
    
    Parameters
    ----------
    filepath : Path, optional
        Path to species data JSON file. Defaults to SPECIES_DATA_FILE.
    
    Returns
    -------
    dict
        Dictionary containing kinetic parameters for each species.
    
    Raises
    ------
    FileNotFoundError
        If the species data file is not found.
    """
    if not filepath.exists():
        raise FileNotFoundError(
            f"Species data file not found: {filepath}\n"
            f"Please ensure species_data.json is in the module directory."
        )
    
    with open(filepath, 'r') as f:
        return json.load(f)


# ============================================================================
# Box model kinetics
# ============================================================================

def kin(
    species: str,
    S: np.ndarray,
    hfc23_yield: float,
    cf3cho_lifetime: float,
    use_lifetime: bool,
    temperature: float,
    y: np.ndarray,
    t: float,
    species_data: Dict = None,
) -> Tuple[float, float, float]:
    """
    Kinetics function for a 1-box atmospheric chemistry model.
    
    Simulates the reaction of a source gas species with atmospheric OH and O3,
    producing CF3CHO as an intermediate, which subsequently photolyzes to HFC-23.
    
    Parameters
    ----------
    species : str
        Name of the source gas species (e.g., 'HFC-143a').
    S : ndarray
        Source terms (emissions) for each species [S_parent, S_CF3CHO, S_HFC23] (Gg/yr).
    hfc23_yield : float
        Fraction of CF3CHO photolysis that produces HFC-23 (unitless).
    cf3cho_lifetime : float
        Atmospheric lifetime of CF3CHO (hours).
    use_lifetime : bool
        If True, use species lifetime instead of explicit OH/O3 kinetics.
    temperature : float
        Temperature (K) for temperature-dependent rate constant calculation.
    y : ndarray
        Current concentrations [parent, CF3CHO, HFC-23] (molecules cm^-3).
    t : float
        Current time (seconds).
    species_data : dict, optional
        Preloaded species kinetic data. If None, will be loaded from file.
    
    Returns
    -------
    tuple of float
        Time derivatives [dA/dt, dB/dt, dC/dt] (molecules cm^-3 s^-1).
        A: parent species concentration
        B: CF3CHO concentration
        C: HFC-23 concentration
    
    Notes
    -----
    The model uses:
    - Temperature-dependent OH rate constant calculation (Arrhenius-type)
    - Fixed O3 concentration in troposphere (~50 ppb)
    - Fixed OH radical concentration (~1e6 cm^-3)
    - Two pathways: OH oxidation and O3 addition/ozonolysis
    """
    if species_data is None:
        species_data = load_species_data()
    
    T = temperature  # K
    
    # Extract kinetic parameters for this species
    k_A = species_data[species]['k_A']  # cm^3 molecule^-1 s^-1
    T_frac = species_data[species]['T_frac']  # K
    T_power = species_data[species]['T_power']  # unitless
    exp_E = species_data[species]['exp_E']  # K
    ozonolysis_yield = species_data[species]['ozonolysis_yield']  # unitless
    cf3cho_yield = species_data[species]['cf3cho_yield']  # unitless
    
    # Calculate temperature-dependent OH rate constant
    kOH = k_A * 1e-12 * (298 / T_frac) ** T_power * np.exp(exp_E / T)
    kO3 = species_data[species]['kO3']  # cm^3 molecule^-1 s^-1
    
    # Atmospheric concentrations
    o3_per_cm3 = O3_CONCENTRATION_PPB * PPBV_TO_PER_CM3
    oh_per_cm3 = OH_CONCENTRATION_PER_CM3
    
    # Loss rate constants
    kCF3CHO = 1.0 / (cf3cho_lifetime * SECONDS_PER_HOUR)  # s^-1
    kHFC23 = 1.0 / (HFC23_LIFETIME_YEARS * SECONDS_PER_YEAR)  # s^-1
    
    # Convert emissions from Gg/yr to molecules cm^-3 s^-1
    S_converted = S * SECONDS_PER_YEAR
    
    # Apply source term only during first second
    if t < 1:
        SA, SB, SC = S_converted
    else:
        SA = SB = SC = 0.0
    
    A, B, C = y
    
    # Calculate derivatives using either lifetime or explicit kinetics
    if use_lifetime:
        kA = 1.0 / (species_data[species]['lifetime'] * SECONDS_PER_YEAR)  # s^-1
        dA = SA - kA * A - kO3 * o3_per_cm3 * A
        dB = SB + kA * A * cf3cho_yield - kCF3CHO * B
        dC = SC + kO3 * o3_per_cm3 * A * ozonolysis_yield + hfc23_yield * kCF3CHO * B - kHFC23 * C
    else:
        # Use explicit OH and O3 kinetics
        dA = SA - kOH * oh_per_cm3 * A - kO3 * o3_per_cm3 * A
        dB = SB + kOH * oh_per_cm3 * A * cf3cho_yield - kCF3CHO * B
        dC = SC + kO3 * o3_per_cm3 * A * ozonolysis_yield + hfc23_yield * kCF3CHO * B - kHFC23 * C
    
    return dA, dB, dC



def integrate_chem(
    species: str,
    S: np.ndarray,
    initial_concs: np.ndarray,
    hfc23_yield: float,
    cf3cho_lifetime: float,
    use_lifetime: bool,
    temperature: float,
    starttime: str,
    endtime: str,
    horizon: int = 100,
    species_data: Dict = None,
) -> Tuple[pd.Series, pd.Series, pd.Series]:
    """
    Integrate a box model chemical system over time.
    
    Parameters
    ----------
    species : str
        Name of the species being modeled (e.g., 'HFC-143a').
    S : ndarray
        Source terms (emissions) for each species (Tg/yr): [S_parent, S_CF3CHO, S_HFC23].
    initial_concs : ndarray
        Initial concentrations for each species (pptv): [A_0, B_0, C_0].
    hfc23_yield : float
        Fraction of CF3CHO photolysis yielding HFC-23 (0-1).
    cf3cho_lifetime : float
        Atmospheric lifetime of CF3CHO (hours).
    use_lifetime : bool
        If True, use species lifetime; if False, use explicit OH/O3 kinetics.
    temperature : float
        Temperature (K) for rate calculations.
    starttime : str
        Start date as string (e.g., '2000-01-01').
    endtime : str
        End date as string.
    horizon : int, optional
        Integration period in years (default: 100).
    species_data : dict, optional
        Preloaded species kinetic data. If None, will be loaded from file.
    
    Returns
    -------
    tuple of pd.Series
        Time series of concentrations (pptv):
        - A: parent species burden
        - B: CF3CHO burden
        - C: HFC-23 burden
    """
    if species_data is None:
        species_data = load_species_data()
    
    lifetime = species_data[species]['lifetime']  # years
    molar_mass = species_data[species]['molar_mass']  # g/mol
    
    # Source terms: convert from Tg/yr to molecules cm^-3 s^-1
    S1_in_Tg, S2_in_Tg, S3_in_Tg = S
    S1 = S1_in_Tg * 1e12 / 5e21 * 28.8 / molar_mass * MOLAR_CONCENTRATION_CM3 / SECONDS_PER_YEAR
    S2 = S2_in_Tg * 1e12 / 5e21 * 28.8 / MOLAR_MASS_CF3CHO * MOLAR_CONCENTRATION_CM3 / SECONDS_PER_YEAR
    S3 = S3_in_Tg * 1e12 / 5e21 * 28.8 / MOLAR_MASS_HFC23 * MOLAR_CONCENTRATION_CM3 / SECONDS_PER_YEAR
    S_converted = np.array([S1, S2, S3])
    
    # Initial concentrations: convert from pptv to molecules cm^-3
    Ai, Bi, Ci = initial_concs
    y_init = np.array([
        Ai * PPTV_TO_PER_CM3,
        Bi * PPTV_TO_PER_CM3,
        Ci * PPTV_TO_PER_CM3
    ])
    
    # Create time array
    times = pd.date_range(start=starttime, end=endtime, freq='D')
    time_array = np.linspace(0.0, horizon * SECONDS_PER_YEAR, len(times))
    
    # Integrate the ODE system
    kinetics_func = partial(kin, species, S_converted, hfc23_yield,
                           cf3cho_lifetime, use_lifetime, temperature,
                           species_data=species_data)
    y_integrated = odeint(kinetics_func, y_init, time_array)
    
    # Convert back to pptv and create pandas Series
    A = pd.Series(data=y_integrated[:, 0] / PPTV_TO_PER_CM3, index=times)
    B = pd.Series(data=y_integrated[:, 1] / PPTV_TO_PER_CM3, index=times)
    C = pd.Series(data=y_integrated[:, 2] / PPTV_TO_PER_CM3, index=times)
    
    return A, B, C



def GWP_species(
    species: str,
    hfc23_yield: float,
    cf3cho_lifetime: float,
    use_lifetime: bool,
    temperature: float,
    horizon: int = 100,
    species_data: Dict = None,
) -> float:
    """
    Calculate the additional Global Warming Potential (GWP) of a species.
    
    Computes the cumulative radiative forcing from HFC-23 produced via CF3CHO photolysis
    from a pulse emission of the parent species, relative to CO2 over a specified horizon.
    
    Parameters
    ----------
    species : str
        Name of the species (e.g., 'HFC-143a').
    hfc23_yield : float
        Fraction of CF3CHO photolysis that produces HFC-23 (0-1).
    cf3cho_lifetime : float
        Atmospheric lifetime of CF3CHO (hours).
    use_lifetime : bool
        If True, use species lifetime; if False, use explicit OH/O3 kinetics.
    temperature : float
        Temperature (K) for rate calculations.
    horizon : int, optional
        Time horizon in years. Must be 20, 100, or 500 (default: 100).
    species_data : dict, optional
        Preloaded species kinetic data. If None, will be loaded from file.
    
    Returns
    -------
    float
        GWP value relative to CO2 over the specified horizon.
    
    Raises
    ------
    ValueError
        If horizon is not in {20, 100, 500}.
    
    References
    ----------
    Radiative efficiency and CO2 forcing values from:
    WMO (World Meteorological Organization) 2022 Scientific Assessment of Ozone Depletion.
    """
    if species_data is None:
        species_data = load_species_data()
    
    if horizon not in TIME_HORIZONS:
        raise ValueError(f"Horizon must be 20, 100, or 500 years. Got {horizon}")
    
    config = TIME_HORIZONS[horizon]
    end_date = config['end_date']
    co2_gwp_forcing = config['gwp_forcing']
    
    # Integrate a pulse emission (1e-9 Tg) to get HFC-23 burden
    A, B, C = integrate_chem(
        species=species,
        S=np.array([1e-9, 0, 0]),
        initial_concs=np.array([0, 0, 0]),
        hfc23_yield=hfc23_yield,
        cf3cho_lifetime=cf3cho_lifetime,
        use_lifetime=use_lifetime,
        temperature=temperature,
        starttime='1700-01-01',
        endtime=end_date,
        horizon=horizon,
        species_data=species_data,
    )
    
    # Integrate HFC-23 burden using trapezoidal rule
    nyears = horizon
    end_idx = min(nyears * 365, len(C) - 1)  # limit to horizon years
    C_horizon = C.iloc[:end_idx]
    
    dt_years = np.diff(C_horizon.index) / np.timedelta64(1, 'D') / 365.0
    C_mid = (C_horizon.values[:-1] + C_horizon.values[1:]) / 2.0
    integrated_C_ppt_years = np.sum(C_mid * dt_years)  # ppt-years per kg
    
    # Calculate GWP
    gwp = HFC23_RADIATIVE_EFFICIENCY * integrated_C_ppt_years / co2_gwp_forcing
    
    return gwp



def plot_species(
    species: str,
    hfc23_yield: float,
    cf3cho_lifetime: float,
    use_lifetime: bool,
    temperature: float,
    horizon: int = 100,
    species_data: Dict = None,
) -> Tuple[plt.Figure, Tuple]:
    """
    Plot time evolution of a species and its degradation products.
    
    Creates three subplots showing the atmospheric burden over time of:
    1. The parent species
    2. CF3CHO intermediate
    3. HFC-23 product
    
    Parameters
    ----------
    species : str
        Name of the species to plot (e.g., 'HFC-143a').
    hfc23_yield : float
        Fraction of CF3CHO photolysis yielding HFC-23 (0-1).
    cf3cho_lifetime : float
        Atmospheric lifetime of CF3CHO (hours).
    use_lifetime : bool
        If True, use species lifetime; if False, use explicit OH/O3 kinetics.
    temperature : float
        Temperature (K) for rate calculations.
    horizon : int, optional
        Integration period in years (default: 100).
    species_data : dict, optional
        Preloaded species kinetic data. If None, will be loaded from file.
    
    Returns
    -------
    fig : plt.Figure
        The matplotlib figure object.
    axes : tuple of Axes
        The three subplots (ax1, ax2, ax3).
    
    Notes
    -----
    A 1 kg pulse emission is used for all species for comparison purposes.
    The plots use 300 dpi for publication quality.
    """
    if species_data is None:
        species_data = load_species_data()
    
    if horizon not in TIME_HORIZONS:
        raise ValueError(f"Horizon must be 20, 100, or 500 years. Got {horizon}")
    
    end_date = TIME_HORIZONS[horizon]['end_date']
    
    # Integrate a pulse emission
    A, B, C = integrate_chem(
        species=species,
        S=np.array([1e-9, 0, 0]),
        initial_concs=np.array([0, 0, 0]),
        hfc23_yield=hfc23_yield,
        cf3cho_lifetime=cf3cho_lifetime,
        use_lifetime=use_lifetime,
        temperature=temperature,
        starttime='1700-01-01',
        endtime=end_date,
        horizon=horizon,
        species_data=species_data,
    )
    
    # Create figure with three subplots
    fig, (ax1, ax2, ax3) = plt.subplots(ncols=3, figsize=(14, 4), dpi=300)
    
    # Plot parent species
    ax1.plot(A, linewidth=1.5)
    ax1.set_title(f'{species} Burden', fontsize=12, fontweight='bold')
    ax1.set_ylabel(f'[{species}] / pptv', fontsize=11)
    ax1.set_xlabel('Time / years', fontsize=11)
    ax1.grid(True, alpha=0.3)
    
    # Plot CF3CHO
    ax2.plot(B, linewidth=1.5, color='orange')
    ax2.set_title('CF3CHO Burden from Degradation', fontsize=12, fontweight='bold')
    ax2.set_ylabel('[CF3CHO] / pptv', fontsize=11)
    ax2.set_xlabel('Time / years', fontsize=11)
    ax2.grid(True, alpha=0.3)
    
    # Plot HFC-23
    ax3.plot(C, linewidth=1.5, color='green')
    ax3.set_title('HFC-23 from CF3CHO Photolysis', fontsize=12, fontweight='bold')
    ax3.set_ylabel('[HFC-23] / pptv', fontsize=11)
    ax3.set_xlabel('Time / years', fontsize=11)
    ax3.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    return fig, (ax1, ax2, ax3)
