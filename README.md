# STOCHEM Code

### Description

This repository contains code to reproduce the results from "Further constraining the role of in-atmosphere production on the global HFC-23 budget" (Adam et al., 2026) published in Atmospheric Chemistry and Physics. 

The results from that study are generated using STOCHEM, a 3D model simulating chemistry and transport in the troposphere. For more information about STOCHEM, see references below. 

This repository contains

- Model run scripts for all twenty model runs described in the manuscript, including the base simulation `sensitivity_base.for`. These are written in Fortran and can be found in the `runs` directory. 
- Auxiliary model files called by the run scripts, including emissions and chemical fields. These can be found in the `model_files` directory.
- Meteorological data used to drive the model. This is located in `stochem_met.tar`, and needs to be unzipped into a directory called `XXXXXX` before use.
- Python code used to read in model outputs, such as global budgets and 4D flux and concentration fields. This is found in `analyse.py`
- The box model used to estimate indirect global warming potentials in the study, written in Python. This is found in `box_model.py`. 


### Warning

### Compiling

### References
