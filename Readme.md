
VIC-Res: prepresentation and operation of reservoirs in the Variable Infiltration Capacity model

-------------- Requirements ------------

The routing component of VIC-Res is developed in FORTRAN. It requires the G77 compiler.
The automatic calibration for model parameters is developed in Python. It requires the following Python packages: (i) Numpy, (ii) os, (iii) multiprocessing, (iv) Matplotlib (optional for plotting), (v) scipy, (vi) numba, and (vii) plotly (optional for parallel coordinate plots). 
It also requires Platypus - a platform for Multi-Objective Evolutionary Algorithms (MOEAs) in Python.
VIC-Res has been tested on Linux Ubuntu 16.04 (Xenial Xerus) operating systems.

-------------- How to run -------------


-------------- Citation --------------

If you use VIC-Res, please cite the following papers:

T.D.Dang, A.K.Chowdhury, and  S.Galelli. On  the  representation  of  waterreservoir storage and operations in large-scale hydrological models:  implicationson model parameterization and climate change impact assessments.Hydrology andEarth System Sciences Discussions, 2019:1–34, 2019. doi:  10.5194/hess-2019-334. 

-------------- Acknowledgement --------------

VIC-Res development is supported by Singapore's Ministry of Education (MoE) via the Tier 2 project "Lining water availability to hydropower supply - an engineering system approach" (Award No. MOE2017-T2-1-143).