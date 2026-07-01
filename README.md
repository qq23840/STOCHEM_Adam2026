# STOCHEM Code

### Description

This repository contains code to reproduce the results from "Further constraining the role of in-atmosphere production on the global HFC-23 budget" (Adam et al., 2026) published in Atmospheric Chemistry and Physics. 

The results from that study are generated using STOCHEM, a 3D model simulating chemistry and transport in the troposphere. For more information about STOCHEM, see references below. 

This repository contains

- Model run scripts for all twenty model runs described in the manuscript, including the base simulation `sensitivity_base.for`. These are written in Fortran and can be found in the `runs` directory. 
- Auxiliary model files called by the run scripts, including emissions and chemical fields. These can be found in the `model_files` directory.
- Config files for each model run script, setting (for example) output directories. These files are found in the `configs` directory. 
- Empty output directories for each of the sensitivity runs, contained in `outputs`. 
- Python code used to read in model outputs, such as global budgets and 4D flux and concentration fields. This is found in `analyse.py`.
- The box model used to estimate indirect global warming potentials in the study, written in Python. This is found in `box_model.py`, and reads in model configuration data from `species_data.json`. 

In order to run the model, the meteorological data used to drive the model needs to be downloaded and unzipped from its Zenodo repository (doi: 10.5281/zenodo.21100197), as it's too big for GitHub. This should be extracted into a directory called `stochem_met`. 

### Warning

This model was originally developed in the early 1990s and has been developed in a linear fashion, as laid out in the first few lines. As a result, it is not particularly modular and the documentation is patchy. There is work ongoing at the University of Bristol to improve these aspects, but if you are interested in using the model then **please first contact Ben Adam (benjamin.adam@bristol.ac.uk) or Rayne Holland (rayne.holland@bristol.ac.uk) to discuss**. 

### Compiling and Running

Before compiling or running, ensure the output directory stated in the relevant `INOPERVM10_*` file 

The model was compiled using the Intel-OneAPI ifx compiler, using:

```
ifx runs/sensitivity_base -o sensitivity_base.o
```

The spin-up year was run using

```
# the model asks for the year (97), month (1), day (1) and run length in days (365)
printf '97\n1\n1\n365\n' | ./sensitivity_base.o
```
Once the model run has completed, the output files `stoch3d.bud`, `stoch3d.dat`, `stoch3d.out`, `stoch3d.sta` and `stoch3d.tot` were deleted from the output directory, and the `DUMP.BIN` file was renamed to `RESTART.BIN`. Then, the model was re-run for the analysis year using

```
printf '98\n1\n1\n365\n' | ./sensitivity_base.o
```

### References

1. Collins, W. J., Stevenson, D. S., Johnson, C. E., & Derwent, R. G. (1997). Tropospheric Ozone in a Global-Scale Three-Dimensional Lagrangian Model and Its Response to NOx Emission Controls. Journal of Atmospheric Chemistry, 26(3), 223–274. https://doi.org/10.1023/A:1005836531979
2. Derwent, R. G., Jenkin, M. E., Stevenson, D. S., Utembe, S. R., Khan, A. H., & Shallcross, D. E. (2025). Influence of the oxidation of non-methane volatile organic compounds on tropospheric hydrogen: A STOCHEM-CRI global Lagrangian model study. Atmospheric Environment, 352, 121214. https://doi.org/10.1016/j.atmosenv.2025.121214
3. Derwent, R. G., Stevenson, D. S., Doherty, R. M., Collins, W. J., & Sanderson, M. G. (2008). How is surface ozone in Europe linked to Asian and North American NOx emissions? Atmospheric Environment, 42(32), 7412–7422. https://doi.org/10.1016/j.atmosenv.2008.06.037
4. Holland, R., Khan, M. A. H., Chhantyal-Pun, R., Orr-Ewing, A. J., Percival, C. J., Taatjes, C. A., & Shallcross, D. E. (2020). Investigating the Atmospheric Sources and Sinks of Perfluorooctanoic Acid Using a Global Chemistry Transport Model. Atmosphere, 11(4), 407. https://doi.org/10.3390/atmos11040407
5. Holland, R., Khan, M. A. H., Driscoll, I., Chhantyal-Pun, R., Derwent, R. G., Taatjes, C. A., Orr-Ewing, A. J., Percival, C. J., & Shallcross, D. E. (2021). Investigation of the Production of Trifluoroacetic Acid from Two Halocarbons, HFC-134a and HFO-1234yf and Its Fates Using a Global Three-Dimensional Chemical Transport Model. ACS Earth and Space Chemistry, 5(4), 849–857. https://doi.org/10.1021/acsearthspacechem.0c00355
6. Khan, M. a. H., Jenkin, M. E., Foulds, A., Derwent, R. G., Percival, C. J., & Shallcross, D. E. (2017). A modeling study of secondary organic aerosol formation from sesquiterpenes using the STOCHEM global chemistry and transport model. Journal of Geophysical Research: Atmospheres, 122(8), 4426–4439. https://doi.org/10.1002/2016JD026415
7. Khan, M. A. H., Lyons, K., Chhantyal-Pun, R., McGillen, M. R., Caravan, R. L., Taatjes, C. A., Orr-Ewing, A. J., Percival, C. J., & Shallcross, D. E. (2018). Investigating the Tropospheric Chemistry of Acetic Acid Using the Global 3-D Chemistry Transport Model, STOCHEM-CRI. Journal of Geophysical Research: Atmospheres, 123(11), 6267–6281. https://doi.org/10.1029/2018JD028529

