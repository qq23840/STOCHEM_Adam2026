"""
Module for reading and analyzing STOCHEM model output files.

This module provides functions to extract and process data from STOCHEM model output files,
including:
- stoch3d.bud: Budget/flux data for reactions
- stoch3d.dat: Species concentration distributions
- stoch3d.sta: Station timeseries observations
- stoch3d.tot: Inventory and flux totals

Functions support extraction of 4D distributions (lat, lon, height, time) and inventory data.

Author: StoChem Project
License: (specify your license)
"""

import numpy as np
import pandas as pd
from pathlib import Path
from typing import Union, List, Optional, Tuple
import warnings

# ============================================================================
# Grid Dimensions
# ============================================================================
NLAT = 36  # Number of latitude points
NLON = 72  # Number of longitude points
NHEIGHT = 9  # Number of vertical levels
NTIME = 12  # Number of time steps (months)

# ============================================================================
# File Structure Constants
# ============================================================================
# stoch3d.bud file parameters
BUD_HEIGHT_STEP = 446  # Number of lines between height levels in stoch3d.bud
BUD_LON_OFFSETS = [225, 299, 373, 3, 77, 151]  # Line offsets for longitude chunks
BUD_LON_CHUNKS = [(0, 12), (12, 24), (24, 36), (36, 48), (48, 60), (60, 72)]  # Longitude ranges

# stoch3d.dat file parameters
DAT_HEIGHT_STEP = 230  # Number of lines between height levels in stoch3d.dat
DAT_LON_OFFSETS = [117, 155, 193, 3, 41, 79]  # Line offsets for longitude chunks
DAT_LON_CHUNKS = [(0, 12), (12, 24), (24, 36), (36, 48), (48, 60), (60, 72)]  # Longitude ranges

# stoch3d.tot file parameters (column positions for parsing)
TOT_SPECIES_VALUE_COLUMN = 48  # Starting column for species inventory values
TOT_FLUX_VALUE_COLUMN = 53  # Starting column for flux inventory values
TOT_REACTION_START_COLUMN = 5  # Starting column for reaction name



def flux_dist(reaction_index: int, folder: Union[str, Path]) -> np.ndarray:
    """
    Extract 4D flux distribution array for a specific reaction from stoch3d.bud.
    
    Parameters
    ----------
    reaction_index : int
        Reaction index corresponding to the reaction in stoch3d.dat and flux10_OH.dat.
        Used to locate the reaction header in stoch3d.bud.
    folder : str or Path
        Path to the directory containing stoch3d.bud file.
    
    Returns
    -------
    np.ndarray
        Array of shape (NLAT, NLON, NHEIGHT, NTIME) containing flux values
        in each grid cell at each time step for the specified reaction. Units
        of molecules per grid cell per 25 days. 
    
    Raises
    ------
    FileNotFoundError
        If stoch3d.bud file doesn't exist in the specified folder.
    ValueError
        If no data found for the specified reaction index.
    IOError
        If there are errors reading or parsing the file.
    
    Notes
    -----
    The stoch3d.bud file has a complex structure with reaction headers followed by
    concentration data arranged by latitude, longitude, and height levels.
    """
    folder = Path(folder)
    file_path = folder / 'stoch3d.bud'
    
    if not file_path.exists():
        raise FileNotFoundError(f"File not found: {file_path}")
    
    # Find line indices where the reaction data begins
    reaction_line_indices = []
    try:
        with open(file_path, 'r') as f:
            for line_num, line in enumerate(f):
                # Match reaction index with proper spacing and reaction arrow marker
                if (line.startswith(f' {reaction_index} ') or 
                    line.startswith(f'  {reaction_index} ')) and '->' in line:
                    # Data starts 2 lines after the header
                    reaction_line_indices.append(line_num + 2)
    except Exception as e:
        raise IOError(f"Error reading {file_path}: {e}")
    
    if not reaction_line_indices:
        raise ValueError(
            f"No data found for reaction index {reaction_index} in {file_path}. "
            f"Check that the index is correct."
        )
    
    # Initialize output array for all times
    flux_data = np.zeros((NLAT, NLON, NHEIGHT, NTIME))
    
    try:
        for time_idx, data_start_line in enumerate(reaction_line_indices):
            time_step_flux = np.zeros((NLAT, NLON, NHEIGHT))
            
            # Read data for each height level
            for height_idx in range(NHEIGHT):
                # Calculate starting line for this height level
                start_line = data_start_line + BUD_HEIGHT_STEP * height_idx
                
                # Initialize grid for this height level
                vertical_level_data = np.zeros((NLAT, NLON))
                
                # Read longitude chunks (they are scattered throughout the file)
                for offset, (lon_start, lon_end) in zip(BUD_LON_OFFSETS, BUD_LON_CHUNKS):
                    grid_data = np.loadtxt(file_path, skiprows=start_line + offset, max_rows=NLAT)
                    vertical_level_data[:, lon_start:lon_end] = grid_data.reshape((NLAT, 12))
                
                time_step_flux[:, :, height_idx] = vertical_level_data
            
            flux_data[:, :, :, time_idx] = time_step_flux
    except Exception as e:
        raise IOError(f"Error parsing flux data from {file_path}: {e}")
    
    return flux_data




def spec_dist(species: str, folder: Union[str, Path]) -> np.ndarray:
    """
    Extract 4D concentration distribution array for a species from stoch3d.dat.
    
    Parameters
    ----------
    species : str
        Species name as it appears in stoch3d.dat (e.g., 'OH', 'O3', 'NO2').
        Must match exactly the species name in chem10_soa.dat.
    folder : str or Path
        Path to the directory containing stoch3d.dat file.
    
    Returns
    -------
    np.ndarray
        Array of shape (NLAT, NLON, NHEIGHT, NTIME) containing concentration values
        in each grid cell at each time step. Units of mol/mol.
    
    Raises
    ------
    FileNotFoundError
        If stoch3d.dat file doesn't exist in the specified folder.
    ValueError
        If no data found for the specified species.
    IOError
        If there are errors reading or parsing the file.
    
    Notes
    -----
    The stoch3d.dat file structure is similar to stoch3d.bud but with different
    spacing between height levels and longitude offsets. Concentration data is
    arranged by height levels with longitude chunks scattered throughout.
    """
    folder = Path(folder)
    file_path = folder / 'stoch3d.dat'
    
    if not file_path.exists():
        raise FileNotFoundError(f"File not found: {file_path}")
    
    # Find line indices where the species data begins
    species_line_indices = []
    try:
        with open(file_path, 'r') as f:
            for line_num, line in enumerate(f):
                if species in line:
                    # Data starts 2 lines after the header
                    species_line_indices.append(line_num + 2)
    except Exception as e:
        raise IOError(f"Error reading {file_path}: {e}")
    
    if not species_line_indices:
        raise ValueError(
            f"No data found for species '{species}' in {file_path}. "
            f"Check that the species name matches exactly."
        )
    
    # Initialize output array for all times
    concentration_data = np.zeros((NLAT, NLON, NHEIGHT, NTIME))
    
    try:
        for time_idx, data_start_line in enumerate(species_line_indices):
            time_step_conc = np.zeros((NLAT, NLON, NHEIGHT))
            
            # Read data for each height level
            for height_idx in range(NHEIGHT):
                # Calculate starting line for this height level
                start_line = data_start_line + DAT_HEIGHT_STEP * height_idx
                
                # Initialize grid for this height level
                vertical_level_data = np.zeros((NLAT, NLON))
                
                # Read longitude chunks in specific order (not sequential in file)
                for offset, (lon_start, lon_end) in zip(DAT_LON_OFFSETS, DAT_LON_CHUNKS):
                    grid_data = np.loadtxt(file_path, skiprows=start_line + offset, max_rows=NLAT)
                    vertical_level_data[:, lon_start:lon_end] = grid_data
                
                time_step_conc[:, :, height_idx] = vertical_level_data
            
            concentration_data[:, :, :, time_idx] = time_step_conc
    except Exception as e:
        raise IOError(f"Error parsing concentration data from {file_path}: {e}")
    
    return concentration_data



def spec_inventory(species: str, folder: Union[str, Path]) -> List[float]:
    """
    Extract time-averaged inventory values for a species from stoch3d.tot.
    
    Parameters
    ----------
    species : str
        Species name as it appears in stoch3d.tot (must match column 1-5).
    folder : str or Path
        Path to the directory containing stoch3d.tot file.
    
    Returns
    -------
    List[float]
        Time-averaged inventory values (in molecules) for the entire model domain,
        one per time step. If multiple sets are found (from spin-up), returns the
        last 12 values (final month).
    
    Raises
    ------
    FileNotFoundError
        If stoch3d.tot file doesn't exist in the specified folder.
    ValueError
        If no data found for the specified species.
    IOError
        If there are errors reading or parsing the file.
    
    Notes
    -----
    The stoch3d.tot file contains three entries per species per timestep:
    TOTM0 (initial), TOTMASS (current), and TOTAVER (average). This function
    extracts every third value (TOTAVER) which represents the time-averaged inventory.
    """
    folder = Path(folder)
    file_path = folder / 'stoch3d.tot'
    
    if not file_path.exists():
        raise FileNotFoundError(f"File not found: {file_path}")
    
    inventory_values = []
    try:
        with open(file_path, 'r') as f:
            for line in f:
                if line.strip().startswith(species):
                    try:
                        # Extract value from column 48 onwards (as per file format)
                        value = float(line[TOT_SPECIES_VALUE_COLUMN:].strip())
                        inventory_values.append(value)
                    except (ValueError, IndexError) as e:
                        warnings.warn(
                            f"Could not parse inventory value for '{species}' from line: "
                            f"{line.strip()}"
                        )
    except Exception as e:
        raise IOError(f"Error reading {file_path}: {e}")
    
    if not inventory_values:
        raise ValueError(f"No data found for species '{species}' in {file_path}")
    
    # Extract averaged inventory: every third value (TOTAVER entries)
    averaged_inventory = inventory_values[2::3]
    
    # Check if we have extra data (likely from spin-up periods)
    if len(averaged_inventory) > NTIME:
        warnings.warn(
            f"Found {len(averaged_inventory)} values for '{species}', expected {NTIME}. "
            f"Returning last {NTIME} values (excluding spin-up)."
        )
        return averaged_inventory[-NTIME:]
    
    return averaged_inventory



def flux_inventory(reaction: str, folder: Union[str, Path]) -> List[float]:
    """
    Extract total flux inventory values for a reaction from stoch3d.tot.
    
    Parameters
    ----------
    reaction : str
        Reaction name as it appears in stoch3d.tot (starting at column 5).
        Format: "A + B => C + D" or similar.
    folder : str or Path
        Path to the directory containing stoch3d.tot file.
    
    Returns
    -------
    List[float]
        Total integrated flux values for each time step. If multiple sets are found
        (from spin-up), returns the last 12 values (final month).
    
    Raises
    ------
    FileNotFoundError
        If stoch3d.tot file doesn't exist in the specified folder.
    ValueError
        If no data found for the specified reaction.
    IOError
        If there are errors reading or parsing the file.
    
    Notes
    -----
    Flux values are integrated over the entire model domain and represent
    the total rate of reaction per time step.
    """
    folder = Path(folder)
    file_path = folder / 'stoch3d.tot'
    
    if not file_path.exists():
        raise FileNotFoundError(f"File not found: {file_path}")
    
    flux_values = []
    try:
        with open(file_path, 'r') as f:
            for line in f:
                # Check if line starts with reaction at column 5
                if line[TOT_REACTION_START_COLUMN:].startswith(reaction):
                    try:
                        # Extract value from column 53 onwards (as per file format)
                        value = float(line[TOT_FLUX_VALUE_COLUMN:].strip())
                        flux_values.append(value)
                    except (ValueError, IndexError) as e:
                        warnings.warn(
                            f"Could not parse flux value for '{reaction}' from line: "
                            f"{line.strip()}"
                        )
    except Exception as e:
        raise IOError(f"Error reading {file_path}: {e}")
    
    if not flux_values:
        raise ValueError(f"No data found for reaction '{reaction}' in {file_path}")
    
    # Check if we have extra data (likely from spin-up periods)
    if len(flux_values) > NTIME:
        warnings.warn(
            f"Found {len(flux_values)} values for '{reaction}', expected {NTIME}. "
            f"Returning last {NTIME} values (excluding spin-up)."
        )
        return flux_values[-NTIME:]
    
    return flux_values



def mod_timeseries(
    site: str,
    folder: Union[str, Path],
    year: int = 2021,
    average: bool = True
) -> pd.DataFrame:
    """
    Extract modelled timeseries data for all species at a specific observation site.
    
    Parameters
    ----------
    site : str
        Site name as it appears in stoch3d.sta (e.g., 'MHD', 'CGO', 'SPO').
    folder : str or Path
        Path to the directory containing stoch3d.sta file.
    year : int, optional
        Year to use for datetime construction (default: 2021). Used because the
        STOCHEM output files contain day-of-year and hour information without year.
    average : bool, optional
        If True (default), average multiple values per timestep to produce one value
        per timestamp. If False, retain all values (may have multiple per timestamp).
    
    Returns
    -------
    pd.DataFrame
        Time-indexed DataFrame with datetime index and columns for each species.
        Values are typically mole fractions in mol/mol.
        Index is pandas DatetimeIndex sorted chronologically.
    
    Raises
    ------
    FileNotFoundError
        If stoch3d.sta file doesn't exist in the specified folder.
    ValueError
        If site not found in file or file format is unexpected.
    IOError
        If there are errors reading or parsing the file.
    
    Notes
    -----
    The stoch3d.sta file structure includes a header with:
    - Number of observation sites
    - Site names
    - Number of species
    - Species names
    
    Data lines follow with site index and species values. Multiple model
    sample times may produce multiple entries per calendar timestamp;
    setting average=True will combine these into one value per time.
    """
    folder = Path(folder)
    file_path = folder / 'stoch3d.sta'
    
    if not file_path.exists():
        raise FileNotFoundError(f"File not found: {file_path}")
    
    # Parse file header to find site and species information
    species_list = []
    site_index_str = None
    
    try:
        with open(file_path, 'r') as f:
            n_sites = int(f.readline().strip())
            
            # Find the site index in the header
            for site_idx in range(n_sites):
                site_name = f.readline().strip()
                if site_name == site:
                    site_index_str = str(site_idx + 1)
            
            if site_index_str is None:
                raise ValueError(f"Site '{site}' not found in {file_path}")
            
            # Read species names
            n_species = int(f.readline().strip())
            for _ in range(n_species):
                species_name = f.readline().strip()
                species_list.append(species_name)
    except (ValueError, StopIteration) as e:
        raise ValueError(f"Error parsing header from {file_path}: {e}")
    
    # Extract timeseries data for the site
    timestamps = []
    data_dict = {sp: [] for sp in species_list}
    
    try:
        with open(file_path, 'r') as f:
            for line in f:
                parts = line.split()
                # Check if line is a site data header (7 columns)
                if len(parts) == 7 and parts[0] == site_index_str:
                    # Parse date/time information
                    month = int(parts[2])
                    day_with_hour_fraction = float(parts[3])
                    
                    # Extract day and hour from day.hour_fraction format
                    hour_fraction, day = np.modf(day_with_hour_fraction)
                    hour = int(hour_fraction * 24)
                    
                    # Create datetime for this entry
                    try:
                        timestamp = pd.to_datetime(
                            f"{year}-{month:02d}-{int(day):02d} {hour:02d}:00:00"
                        )
                        timestamps.append(timestamp)
                    except ValueError:
                        warnings.warn(
                            f"Invalid date parsing: year={year}, month={month}, "
                            f"day={int(day)}, hour={hour}. Skipping entry."
                        )
                        # Skip this entry and the following species lines
                        for _ in species_list:
                            f.readline()
                        continue
                    
                    # Read species data (one species per line following the header)
                    for sp in species_list:
                        species_line = f.readline()
                        try:
                            value = float(species_line.strip())
                            data_dict[sp].append(value)
                        except ValueError:
                            warnings.warn(
                                f"Could not parse value for {sp}: {species_line.strip()}"
                            )
                            data_dict[sp].append(np.nan)
    except Exception as e:
        raise IOError(f"Error parsing timeseries data from {file_path}: {e}")
    
    if not timestamps:
        raise ValueError(f"No timeseries data found for site '{site}' in {file_path}")
    
    # Create DataFrame with datetime index
    timeseries_df = pd.DataFrame(index=timestamps, data=data_dict)
    
    # Average duplicate timestamps if requested
    if average:
        timeseries_df = timeseries_df.groupby(timeseries_df.index).mean()
    
    return timeseries_df

