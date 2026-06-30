"""
Module for reading and analyzing STOCHEM model output files.

This module provides functions to extract and process data from STOCHEM output files:
- stoch3d.bud: Budget/flux data
- stoch3d.dat: Species concentration distributions
- stoch3d.sta: Station timeseries data
- stoch3d.tot: Inventory totals
"""

import numpy as np
import pandas as pd
from pathlib import Path
from typing import Union, List, Optional
import warnings

# Grid dimensions constants
NLAT = 36  # Number of latitude points
NLON = 72  # Number of longitude points
NHEIGHT = 9  # Number of vertical levels
NTIME = 12  # Number of time steps (months)


def flux_dist(index: int, folder: Union[str, Path]) -> np.ndarray:
    """
    Extract 4D flux distribution array for a specific reaction from stoch3d.bud.
    
    Args:
        index: Reaction index corresponding to the reaction in stoch3d.dat and flux10_OH.dat
        folder: Path to the directory containing stoch3d.bud
    
    Returns:
        np.ndarray: Array of shape (NLAT, NLON, NHEIGHT, NTIME) containing flux values
                   in each grid cell at each time step for the specified reaction.
                   Units depend on the reaction type.
    
    Raises:
        FileNotFoundError: If stoch3d.bud file doesn't exist
        ValueError: If no data found for the specified index
    """
    folder = Path(folder)
    fname = folder / 'stoch3d.bud'
    
    if not fname.exists():
        raise FileNotFoundError(f"File not found: {fname}")
    
    # Find line indices for the reaction
    lines = []
    try:
        with open(fname, 'r') as f:
            for j, line in enumerate(f):
                # Match reaction index with proper spacing
                if (line.startswith(f' {index} ') or line.startswith(f'  {index} ')) and '->' in line:
                    lines.append(j + 2)
    except Exception as e:
        raise IOError(f"Error reading {fname}: {e}")
    
    if not lines:
        raise ValueError(f"No data found for reaction index {index} in {fname}")
    
    # Initialize output array
    a_tot = np.zeros((NLAT, NLON, NHEIGHT, NTIME))
    
    # Constants for file structure
    HEIGHT_OFFSET = 446  # Lines between height levels
    LON_OFFSETS = [225, 299, 373, 3, 77, 151]  # Offsets for longitude chunks
    LON_CHUNKS = [(0, 12), (12, 24), (24, 36), (36, 48), (48, 60), (60, 72)]
    
    try:
        for t, s0 in enumerate(lines):
            a_time = np.zeros((NLAT, NLON, NHEIGHT))
            
            for i in range(NHEIGHT):
                s = s0 + HEIGHT_OFFSET * i
                a = np.zeros((NLAT, NLON))
                
                # Read longitude chunks
                for offset, (lon_start, lon_end) in zip(LON_OFFSETS, LON_CHUNKS):
                    data = np.loadtxt(fname, skiprows=s + offset, max_rows=NLON)
                    a[:, lon_start:lon_end] = data.reshape((NLAT, 12))
                
                a_time[:, :, i] = a
            
            a_tot[:, :, :, t] = a_time
    except Exception as e:
        raise IOError(f"Error parsing flux data from {fname}: {e}")
    
    return a_tot


def spec_dist(species: str, folder: Union[str, Path]) -> np.ndarray:
    """
    Extract 4D concentration distribution array for a species from stoch3d.dat.
    
    Args:
        species: Species name as it appears in stoch3d.dat (e.g., 'OH', 'O3')
                Must match the species name in chem10_soa.dat
        folder: Path to the directory containing stoch3d.dat
    
    Returns:
        np.ndarray: Array of shape (NLAT, NLON, NHEIGHT, NTIME) containing concentration
                   values (typically in molecules/cm³) in each grid cell at each time step.
    
    Raises:
        FileNotFoundError: If stoch3d.dat file doesn't exist
        ValueError: If no data found for the specified species
    """
    folder = Path(folder)
    fname = folder / 'stoch3d.dat'
    
    if not fname.exists():
        raise FileNotFoundError(f"File not found: {fname}")
    
    # Find line indices for the species
    lines = []
    try:
        with open(fname, 'r') as f:
            for j, line in enumerate(f):
                if species in line:
                    lines.append(j + 2)
    except Exception as e:
        raise IOError(f"Error reading {fname}: {e}")
    
    if not lines:
        raise ValueError(f"No data found for species '{species}' in {fname}")
    
    # Initialize output array
    a_tot = np.zeros((NLAT, NLON, NHEIGHT, NTIME))
    
    # Constants for file structure
    HEIGHT_OFFSET = 230  # Lines between height levels
    LON_OFFSETS = [117, 155, 193, 3, 41, 79]  # Offsets for longitude chunks (order matters!)
    LON_CHUNKS = [(0, 12), (12, 24), (24, 36), (36, 48), (48, 60), (60, 72)]
    
    try:
        for t, s0 in enumerate(lines):
            a_time = np.zeros((NLAT, NLON, NHEIGHT))
            
            for i in range(NHEIGHT):
                s = s0 + HEIGHT_OFFSET * i
                a = np.zeros((NLAT, NLON))
                
                # Read longitude chunks in specific order
                for offset, (lon_start, lon_end) in zip(LON_OFFSETS, LON_CHUNKS):
                    data = np.loadtxt(fname, skiprows=s + offset, max_rows=NLAT)
                    a[:, lon_start:lon_end] = data
                
                a_time[:, :, i] = a
            
            a_tot[:, :, :, t] = a_time
    except Exception as e:
        raise IOError(f"Error parsing concentration data from {fname}: {e}")
    
    return a_tot


def spec_inventory(species: str, folder: Union[str, Path]) -> List[float]:
    """
    Extract time-averaged inventory values for a species from stoch3d.tot.
    
    Args:
        species: Species name as it appears in stoch3d.tot
        folder: Path to the directory containing stoch3d.tot
    
    Returns:
        List[float]: Time-averaged inventory values (typically in kg or Tg) for each time step.
                    Returns the last 12 values if more are found (from spin-up periods).
    
    Raises:
        FileNotFoundError: If stoch3d.tot file doesn't exist
        ValueError: If no data found for the specified species
    
    Notes:
        Returns every third value to extract averaged inventory, skipping TOTM0 and TOTMASS.
    """
    folder = Path(folder)
    fname = folder / 'stoch3d.tot'
    
    if not fname.exists():
        raise FileNotFoundError(f"File not found: {fname}")
    
    totals = []
    try:
        with open(fname, 'r') as f:
            for line in f:
                if line.strip().startswith(species):
                    try:
                        # Extract value from position 48 to end (excluding newline)
                        value = float(line[48:].strip())
                        totals.append(value)
                    except (ValueError, IndexError) as e:
                        warnings.warn(f"Could not parse value from line: {line.strip()}")
    except Exception as e:
        raise IOError(f"Error reading {fname}: {e}")
    
    if not totals:
        raise ValueError(f"No data found for species '{species}' in {fname}")
    
    # Extract averaged inventory (every third value)
    output = totals[2::3]
    
    if len(output) > NTIME:
        warnings.warn(
            f"Found {len(output)} values for '{species}', expected {NTIME}. "
            f"Returning last {NTIME} values (likely excluding spin-up)."
        )
        return output[-NTIME:]
    
    return output


def flux_inventory(reaction: str, folder: Union[str, Path]) -> List[float]:
    """
    Extract total flux inventory values for a reaction from stoch3d.tot.
    
    Args:
        reaction: Reaction name as it appears in stoch3d.tot (starting at column 5)
        folder: Path to the directory containing stoch3d.tot
    
    Returns:
        List[float]: Total flux values for each time step. Returns the last 12 values
                    if more are found (from spin-up periods).
    
    Raises:
        FileNotFoundError: If stoch3d.tot file doesn't exist
        ValueError: If no data found for the specified reaction
    """
    folder = Path(folder)
    fname = folder / 'stoch3d.tot'
    
    if not fname.exists():
        raise FileNotFoundError(f"File not found: {fname}")
    
    totals = []
    try:
        with open(fname, 'r') as f:
            for line in f:
                if line[5:].startswith(reaction):
                    try:
                        # Extract value from position 53 to end (excluding newline)
                        value = float(line[53:].strip())
                        totals.append(value)
                    except (ValueError, IndexError) as e:
                        warnings.warn(f"Could not parse value from line: {line.strip()}")
    except Exception as e:
        raise IOError(f"Error reading {fname}: {e}")
    
    if not totals:
        raise ValueError(f"No data found for reaction '{reaction}' in {fname}")
    
    if len(totals) > NTIME:
        warnings.warn(
            f"Found {len(totals)} values for '{reaction}', expected {NTIME}. "
            f"Returning last {NTIME} values (likely excluding spin-up)."
        )
        return totals[-NTIME:]
    
    return totals


def mod_timeseries(
    site: str, 
    folder: Union[str, Path], 
    year: int = 2021, 
    average: bool = True
) -> pd.DataFrame:
    """
    Extract modelled timeseries data for all species at a specific observation site.
    
    Args:
        site: Site name as it appears in stoch3d.sta (e.g., 'MHD', 'CGO')
        folder: Path to the directory containing stoch3d.sta
        year: Year to use for datetime construction (default: 2021)
        average: If True, average multiple values per timestep. If False, retain all values.
    
    Returns:
        pd.DataFrame: Time-indexed DataFrame with columns for each species.
                     Index is pandas DatetimeIndex. Values are typically mole fractions (ppt).
    
    Raises:
        FileNotFoundError: If stoch3d.sta file doesn't exist
        ValueError: If site not found in file or file format is unexpected
    
    Notes:
        Multiple model values may exist for each timestep. By default these are averaged
        to produce a single value per timestamp, but set average=False to keep all values.
    """
    folder = Path(folder)
    fname = folder / 'stoch3d.sta'
    
    if not fname.exists():
        raise FileNotFoundError(f"File not found: {fname}")
    
    # Parse header to find site index and species list
    species = []
    site_index = None
    
    try:
        with open(fname, 'r') as f:
            n_sites = int(f.readline().strip())
            
            # Find the site index
            for i in range(n_sites):
                site_name = f.readline().strip()
                if site_name == site:
                    site_index = str(i + 1)
            
            if site_index is None:
                raise ValueError(f"Site '{site}' not found in {fname}")
            
            n_species = int(f.readline().strip())
            for i in range(n_species):
                species_name = f.readline().strip()
                species.append(species_name)
    except (ValueError, StopIteration) as e:
        raise ValueError(f"Error parsing header from {fname}: {e}")
    
    # Extract timeseries data for the site
    times = []
    data = {sp: [] for sp in species}
    
    try:
        with open(fname, 'r') as f:
            for line in f:
                parts = line.split()
                # Check if line is a site header with 7 columns
                if len(parts) == 7 and parts[0] == site_index:
                    # Parse date/time information
                    month = int(parts[2])
                    day_with_hour = float(parts[3])
                    hour_frac, day = np.modf(day_with_hour)
                    hour = int(hour_frac * 24)
                    
                    # Create datetime
                    try:
                        datetime = pd.to_datetime(f"{year}-{month:02d}-{int(day):02d} {hour:02d}:00:00")
                        times.append(datetime)
                    except ValueError:
                        warnings.warn(f"Invalid date: year={year}, month={month}, day={int(day)}, hour={hour}")
                        continue
                    
                    # Read species data (one species per line following the header)
                    for sp in species:
                        species_line = f.readline()
                        try:
                            data[sp].append(float(species_line.strip()))
                        except ValueError:
                            warnings.warn(f"Could not parse value for {sp}: {species_line.strip()}")
                            data[sp].append(np.nan)
    except Exception as e:
        raise IOError(f"Error parsing timeseries data from {fname}: {e}")
    
    if not times:
        raise ValueError(f"No timeseries data found for site '{site}' in {fname}")
    
    # Create DataFrame
    df = pd.DataFrame(index=times, data=data)
    
    # Average duplicate timestamps if requested
    if average:
        df = df.groupby(df.index).mean()
    
    return df

