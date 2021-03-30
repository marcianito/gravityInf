Repository containing scripts related to _Reich et al. (2021)_ publication
======================================================================
**Article:**  
Marvin Reich<sup>1</sup>, Michal Mikolaj<sup>1</sup>, Theresa Blume<sup>1</sup>, Andreas Güntner<sup>1,2</sup>: **2021**. **Field-scale subsurface flow processes inferred from continuous gravity monitoring during a sprinkling experiment**, (in review)

**Affiliation**  
<sup>1</sup>Helmholtz Centre Potsdam GFZ German Research Centre for Geosciences, Section Hydrology, 14473 Potsdam, Germany  
<sup>2</sup>University of Potsdam, Institute of Earth and Environmental Science, 14476 Potsdam, Germany
> _Please cite this article when using the here provided scripts_

## Description

This is a small R-package which aims at supporting the anaylsis of a sprinkling (infiltration) experiment in
combination with simultaneous and continious gravity measurements, 
presented in the above mentioned paper.
With this package you can easily walk through the necessary steps in order to set up an infiltration scenario,
maybe based on your own sprinkling / irrigation experiment and carry out simple hydrological modelling of water distribution
in 3D in the subsurface.
An observed gravity time series is needed for the model in order to fit and thus identify
the dominant infiltration process for your research area.

A model functionality and limitations can be found in the above mentioned publication.


For bug fixes, comments or further development please contact: mreich@posteo.de.

## Installation

1. Start R

2. Install all dependencies

3. Install package via devtools: 
`devtools::install_github("marcianito/gravityInf")`

4. load packages: 
`library(HyGra)`; 
`library(gravityInf)`

## Dependencies

### Computationally
* r-base version 3.3.1
* gravity R-package [HyGra](http://github.com/marcianito/HyGra) (for standard routines): `devtools::install_github("marcianito/HyGra")`
* [ppso](http://github.com/TillF/ppso): `devtools::install_github("TillF/ppso")`
* following R-packages: devtools, dplyr, ggplot2, reshape2, viridis, xts, zoo, doParralel, foreach, ptinpoly, spacetime, data.table, automap
* system libraries for devtools

in debian install using: 
`apt-get update && apt-get install libssl-dev`

**Warning**: depending on your model discretization, in both space and time, it might be
necessary to run this analysis on a cluster (or have at least a high performance machine).

### Data-wise
It is necessary to have a time series of observed gravity data (could be synthetic).

## Processing

1. Start R
2. Load infiltration_example.r script
3. Modify according to description below
4. Run script and look at outputted results

## Infiltration modeling procedure (computations)
#### For more details, please look at the vignette or the corresponding help-files (within R).

All changes should be done in a new file following (or a copy of) the infiltration_example.r file.

General setup:

* Directory and configs (input / output, file extentions, enable plotting)

Gravity grid setup: 

* Gravimeter coordinates (x, y, z + height of sensor)
* Model domain (x and y extensions) and discretization
* Input file to load for DEM input
* Input file settings (data dimensions, order of data columns from input files)

Inverse sprinkling modeling setup:

* Input file of observed gravity signal
* Input file water intensity distribution
* Interpolation method
* Number of model runs
* Objective function
* Senario types for different layers
* Duration of sprinkling
* Select combined processes and/or lateral flow
* Allowed mass balance error
* Input water volume of experiment
* Model parameter values (min, max)
* Model parameter starting values
* Plotting options

