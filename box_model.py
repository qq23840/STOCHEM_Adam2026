# Functions to plot box model results over time for combined kinetics (OH and O3)
# Also GWP calculations. 

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint
import pandas as pd
import json

def kin(species,
        S,
        hfc23_yield,
        cf3cho_lifetime,
        use_lifetime,
        temperature,
        y, 
        t):
    """
    Kinetics function for a 1-box model. Reads in a single source gas species, 
    and outputs kinetics for its reaction with Ozone and OH to produce CF3CHO
    
    args:
        species: string, name of species
        S: input for species in Gg
        y: list of initial concentrations [A,B,C]
        t: time array in seconds
    returns:
        dA, dB, dC: derivatives of concentrations with respect to time"""
    
    with open('/user/work/qq23840/stochem_output/ben_box_model/species_data.json') as f:
        species_data = json.load(f)

    T = temperature # K

    k_A = species_data[species]['k_A'] # in cm3/molecule/s
    T_frac = species_data[species]['T_frac'] # in K
    T_power = species_data[species]['T_power'] # unitless
    exp_E = species_data[species]['exp_E'] # in K
    kOH = k_A * 1e-12* (298/T_frac)**T_power * np.exp(exp_E/T) # in cm3/molecule/s
    kO3 = species_data[species]['kO3'] # in cm3/molecule/s
    ozonolysis_yield = species_data[species]['ozonolysis_yield'] # unitless fraction
    cf3cho_yield = species_data[species]['cf3cho_yield'] # unitless fraction

    o3_conc = 50 # ppb
    ppbv_to_per_cm3 = 2.5e10
    o3_per_cm3 = o3_conc * ppbv_to_per_cm3

    oh_per_cm3 = 1e6 # 1e6 per cm3

    kCF3CHO = 1/(cf3cho_lifetime*3600) # CHECK THIS, current 72 hours

    kHFC23 = 1/(243*365*24*3600) # 243 years

    yr_in_s = 60*60*24*365
    S = S * yr_in_s

    if t<1:
        SA, SB, SC = S
    else:
        SA = 0.
        SB = 0.
        SC = 0.
    
    
    A,B,C = y

    # using OH concentrations
    # dA = SA - kOH*oh_per_cm3*A - kO3*o3_per_cm3*A # loss of main species A due to OH and O3
    # dB = SB + kOH*oh_per_cm3*A*cf3cho_yield - kCF3CHO*B  # B is CF3CHO
    # dC = SC + kO3*o3_per_cm3*A*ozonolysis_yield + hfc23_yield*kCF3CHO*B- kHFC23*C # HFC-23. 

    if use_lifetime:

        kA = 1/(species_data[species]['lifetime']*yr_in_s) # in s-1

        dA = SA - kA*A - kO3*o3_per_cm3*A # loss of main species A due to OH and O3
        dB = SB + kA*A*cf3cho_yield - kCF3CHO*B  # B is CF3CHO
        dC = SC + kO3*o3_per_cm3*A*ozonolysis_yield + hfc23_yield*kCF3CHO*B- kHFC23*C # HFC-23. 
    else: 
        dA = SA - kOH*oh_per_cm3*A - kO3*o3_per_cm3*A # loss of main species A due to OH and O3
        dB = SB + kOH*oh_per_cm3*A*cf3cho_yield - kCF3CHO*B  # B is CF3CHO
        dC = SC + kO3*o3_per_cm3*A*ozonolysis_yield + hfc23_yield*kCF3CHO*B- kHFC23*C # HFC-23. 
    return dA,dB,dC

def integrate_chem(species, S, initial_concs, fin, starttime, endtime, hfc23_yield, cf3cho_lifetime, use_lifetime, temperature, horizon=100):
    """
    Integrate a box model chemical system over time.
    
    Parameters:
        species: str, name of the species being modeled
        S: array-like, source terms (e.g., emissions) for each species (in Tg/yr)
        initial_concs: array-like, initial concentrations for each species (in pptv)
        fin: function, kinetic function returning derivatives
        starttime: str, start date for time series
        endtime: str, end date for time series
        hfc23yield: float, yield of HFC-23 from CF3CHO breakdown
        horizon: float, number of years to integrate (20, 100 or 500)
    Returns:
        A, B, C: pandas.Series, time series of concentrations (in pptv) for each species
    """
    from functools import partial

    with open('/user/work/qq23840/stochem_output/ben_box_model/species_data.json') as f:
        species_data = json.load(f)

    lifetime = species_data[species]['lifetime'] # in years
    molar_mass = species_data[species]['molar_mass'] # in g/mol

    # conversion factors
    pptv_to_per_cm3 = 2.5e7
    yr_in_s = 60*60*24*365
    # source terms (if required)
    S1_in_Tg,S2_in_Tg,S3_in_Tg = S
    
    S1 = S1_in_Tg*1e12/5e21*28.8/molar_mass*2.5e19/yr_in_s # molar mass of HFC-143a is 84.4
    S2 = S2_in_Tg*1e12/5e21*28.8/90*2.5e19/yr_in_s # molar mass of CF3CHO is 90.0
    S3 = S3_in_Tg*1e12/5e21*28.8/70.3*2.5e19/yr_in_s # molar mass of HFC-23 is 70.3
    Sinit = np.array([S1,S2,S3])
    # initial concentrations in pptv need converting to chemical scale
    Ai, Bi, Ci = initial_concs
    # convert to per cm3 from pptv
    yinit  = np.array([Ai*pptv_to_per_cm3 , Bi*pptv_to_per_cm3 , Ci*pptv_to_per_cm3])

    times = pd.date_range(start=starttime, end=endtime, freq='D')

    # set up time integration
    time = np.linspace(start=0.0, stop=horizon*yr_in_s, num=len(times))
    # solve the chemistry
    y2 =  odeint(partial(fin, species, Sinit, hfc23_yield, cf3cho_lifetime, use_lifetime, temperature), yinit, time)
    # generate results as pandas timeseries of mole fraction in pptv from the per cm3 output
    times = pd.date_range(start=starttime, end=endtime, freq='D')
    A = pd.Series(data=y2[:,0]/pptv_to_per_cm3, index=times)
    B = pd.Series(data=y2[:,1]/pptv_to_per_cm3, index=times)
    C = pd.Series(data=y2[:,2]/pptv_to_per_cm3, index=times)
    
    return A,B,C

def GWP_species(species, hfc23_yield, cf3cho_lifetime, use_lifetime, temperature,horizon=100):
    """
    Calculate the additional Global Warming Potential (GWP) of a species over a given time horizon
    due to its degradation to HFC-23.
    
    Parameters:
        species: str
            Name of the species to calculate GWP for.
        horizon: int, optional
            Time horizon in years (default: 100). Supported: 20, 100, 500.
    Returns:
        float
            GWP value relative to CO2 over the specified horizon.
    """

    if horizon == 20:
        endtime='1720-01-01'
        CO2_GWP = 2.434e-14 # 20 year horizon, WMO2022, in W m-2 yr kg-1
    elif horizon == 100:
        endtime='1800-01-01'
        CO2_GWP = 8.947e-14 # 100 year horizon, WMO2022, in W m-2 yr kg-1
    elif horizon == 500:
        endtime='2200-01-01'
        CO2_GWP = 3.138e-13 # 500 year horizon, WMO2022, in W m-2 yr kg-1

    A, B, C = integrate_chem(species, 
                             S=[1e-9,0,0], 
                             
                             initial_concs=[0,0,0], 
                             fin=kin, starttime='1700-01-01',
                             endtime=endtime, 
                             hfc23_yield=hfc23_yield, 
                             cf3cho_lifetime=cf3cho_lifetime,
                             use_lifetime=use_lifetime,
                             horizon=horizon,
                             temperature=temperature)

    r_eff = 0.192*1e-3 # W per m2 per ppt for HFC-23 (from WMO 2022 report)

    # Integrate HFC-23 mole fraction (ppt) over 100 years to get ppt-years using daily time index
    nyears = 100
    end_date = C.index[0] + pd.Timedelta(days=nyears*365)
    C_horizon = C[C.index <= end_date]
    dt_years = np.diff(C_horizon.index) / np.timedelta64(1, 'D') / 365.0  # time steps in years (daily index)
    C_mid = (C_horizon.values[:-1] + C_horizon.values[1:]) / 2  # trapezoidal rule
    integrated_C_ppt_years = np.sum(C_mid * dt_years) # in ppt years per kg pulse

    GWP_extra = r_eff * integrated_C_ppt_years / CO2_GWP

    return GWP_extra

def plot_species(species, hfc23_yield, cf3cho_lifetime, use_lifetime, temperature, horizon=100):
    """
    Plot the time evolution of a species and its degradation products in the box model.
    
    Parameters:
        species: str
            Name of the species to plot.
        horizon: int, optional
            Time horizon in years (default: 100). Supported: 20, 100, 500.
    Returns:
        None
            Displays matplotlib plots of the species, CF3CHO, and HFC-23 burdens over time.
    """

    if horizon == 20:
        endtime='1720-01-01'
    elif horizon == 100:
        endtime='1800-01-01'
    elif horizon == 500:
        endtime='2200-01-01'

    A, B, C = integrate_chem(species, 
                             S=[1e-9,0,0], 
                             initial_concs=[0,0,0], 
                             fin=kin, 
                             starttime='1700-01-01', 
                             endtime=endtime, 
                             hfc23_yield=hfc23_yield,
                             cf3cho_lifetime=cf3cho_lifetime,
                             use_lifetime=use_lifetime,
                             horizon=horizon,
                             temperature=temperature)

    fig, (ax1, ax2, ax3) = plt.subplots(ncols=3, figsize=(12,5), dpi=300)
    ax1.plot(A)
    ax1.set_title(f'{species} burden')
    ax1.set_ylabel(f'[{species}] / pptv')
    ax2.plot(B)
    ax2.set_title(f'CF3CHO burden from {species} degradation')
    ax2.set_ylabel('[CF3CHO] / pptv')
    ax3.plot(C)
    ax3.set_title('HFC-23 from CF3CHO photolysis')
    ax3.set_ylabel('[HFC-23] / pptv')
    plt.tight_layout()