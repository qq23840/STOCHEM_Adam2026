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

In order to run the model, the meteorological data used to drive the model needs to be downloaded and unzipped from its Zenodo repository, as it's too big for GitHub. This should be extracted into a directory called `stochem_met`. 

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
printf '97\n1\n1\n365\n' | ./sensitiviy_base.o
```
Once the model run has completed, the output files `stoch3d.bud`, `stoch3d.dat`, `stoch3d.out`, `stoch3d.sta` and `stoch3d.tot` were deleted from the output directory, and the `DUMP.BIN` file was renamed to `RESTART.BIN`. Then, the model was re-run for the analysis year using

```
printf '98\n1\n1\n365\n' | ./sensitiviy_base.o
```

### References
