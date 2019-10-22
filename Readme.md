 ![Build](https://img.shields.io/badge/VIC--Res-developing-orange) ![Version](https://img.shields.io/badge/version-1.0-blue)

# VIC-Res: Representation and operation of reservoirs in the Variable Infiltration Capacity (VIC) model

Variable Infiltration Capacity (VIC) is a popular macro-scale hydrological model developed and maintained by the University of Washington, US. VIC has been used extensively to model river basins around the world. However, the operation of artificial reservoirs cannot be modelled in the original version of VIC. Since the environmental impact of reservoirs comes under ever-increasing scrutiny as the global demand for water and energy increases, we develop a variant of VIC which allows modelling explicitly reservoirs. More details about how we incorporate reservoirs into the model are presented in [Dang et al. (2019)](https://www.hydrol-earth-syst-sci-discuss.net/hess-2019-334/).

### Prerequisites

We choose VIC 4.2 which is a prominent version in terms of modelling speed and the number of active users to develop VIC-Res. In this version, the rainfall-runoff component is written in C, and the routing component is developed in FORTRAN. Both VIC and VIC-Res have been developed for uses on Linux and Unix platforms, and they require the GNU and G77 compilers. VIC-Res has been tested on Linux Ubuntu 16.04 (Xenial Xerus) operating systems.

### How to run

Running VIC and VIC-Res is considered to be a "formidable task", and this task can be implemented more easily in steps:
1. Downloading/compiling the source code (by using GNU and G77)
2. Preparing input files:
* Climate forcings (e.g. precipitation, maximum/minimum temperature, wind, etc.)
* Land use
* Soil parameters
* Flow direction
* Location of basin outlet
* Location of reservoirs (new)
* Reservoir parameters (new)
3. Configuring output files

More information on how to run the original version of VIC can be found here:
```
https://vic.readthedocs.io/en/master/
```
The sequence of syntax to run VIC-Res is similar to VIC (these two versions are only different on how they are customized).

### How to customize

Reservoir representation:
The location of reservoirs is stored in *reservoirlocation.txt* which has the same number of rows and columms as the flow direction matrix. This file contains integer numbers:
* '1' - '9998' represents for the ID of reservoirs
* '9999' represents open water surface (reservoir cell)
* '0' is other types of land use

Reservoir operation:
For each reservoir ID (e.g. 1), there is a reservoir configuration file (e.g. res1.txt), containing parameters such as maximum water level (*m*), minimum water level (*m*), storage capacity (*1000 m3*), water head (*m*; optional for hydropower production estimation), design discharge (*m3/s*), commision year, initial water volume (*1000 m3*), name of reservoir, seepage rate (*m3/s*), infiltration rate (*m3/s*), and characteristics of rule curves. 

### Example 

We use a simplified catchment with a dimesion of 5 x 5 (cell) to demonstrate the capacity of VIC-Res (see *Rainfall-runoffSetup* and *RoutingSetup* folders). 

Reservoir representation: the 2D matrix in *reservoirlocation.txt* shows that there are two cascade reservoirs (Reservoir 1 and Reservoir 2, IDs 1 and 2). The number of reservoir cells are 1 and 2, respectively.

Reservoir operation: The parameter files of these two reservoirs are named *res1.txt* and *res2.txt*.

### Citation

If you use VIC-Res, please cite the following paper:

*T.D. Dang, A.K. Chowdhury, and  S. Galelli. On  the  representation  of  waterreservoir storage and operations in large-scale hydrological models:  implications on model parameterization and climate change impact assessments. Hydrology andEarth System Sciences Discussions, 2019:1–34, 2019.* ![DOI](https://img.shields.io/badge/DOI-doi.org%2F10.5194%2Fhess--2019--334-lightgrey)

## Acknowledgments

VIC-Res development is supported by Singapore's Ministry of Education (MoE) via the Tier 2 project "Linking water availability to hydropower supply - an engineering system approach" (Award No. MOE2017-T2-1-143).

### Contact

Questions and feedback related to VIC-Res (not VIC) and requests to fix possible bugs, please send to thanhiwer@gmail.com (Thanh Dang), Resilient Water Systems Group, Singapore University of Technology and Design, Singapore.

### Possible future works

* Groundwater modelling
* New Graphic User Interface (GUI) in Java 
