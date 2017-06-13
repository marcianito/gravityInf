Repository containing scripts related to _Reich et al., (2017, ???)_ publication
======================================================================
**Article:**  
Marvin Reich<sup>1</sup>, Michal Mikolaj<sup>1</sup>, Theresa Blume<sup>1</sup>, Andreas Güntner<sup>1,4</sup>: **_TITLE_**, published in: [???](http://link.com), **2017**

**Affiliation**  
<sup>1</sup>Helmholtz Centre Potsdam GFZ German Research Centre for Geosciences, Section 5.4 Hydrology, 14473 Potsdam, Germany  
<sup>2</sup>University of Potsdam, Institute of Earth and Environmental Science, 14476 Potsdam, Germany
> _Please cite this article when using here provided scripts_

## Description

This is a small R-package which aims at supporting the anaylsis of a sprinkling (infiltration) experiment in
combination with simultaneous and continious gravity measurements, 
presented in the above mentioned paper.
With this package you can easily walk through the necessary steps in order to set up an infiltration scenario,
maybe based on your own sprinkling / irrigation experiment and carry out simple hydrological modelling of water distribution
in 3D in the subsurface.
An observed gravity time series is needed for the model in order to fit and thus identify
the dominant infiltration process for your research area.

For bug fixes, comments or further development please contact: mreich@posteo.de.

## Installation

1. Start R
2. Install package-dependence via devtools: 
`devtools::install_github("marcianito/UmbrellaEffect")`

3. Install package via devtools: 
`devtools::install_github("marcianito/gravityInf")`

4. load packages: 
`library(UmbrellaEffect)`; 
`library(gravityInf)`

## Dependencies

### Computationally
* r-base version 3.3.1
* other gravity R-package, containing basic functions for gravity grid setups: UmbrellaEffect
* for further dependencies of UmbrellaEffect package, please visit [here](http://github.com/marcianito/UmbrellaEffect)
* following R-packages: devtools, dplyr, ggplot2, gstat, ..., reshape2, viridis, xts, zoo, doParralel, foreach, ...
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

All changes should be done in a new file following (or a copy of) the reduction_example.r file.

(Step 2 is only explanatory for what the script does; no modifications necessary.)

1. Setup: 
	* _Directory and configs (input / output, file extentions, enable plotting)_
	* _Gravimeter coordinates (x, y, z + height of sensor)_
	* _Model domain (x and y extensions)_
	* _Discretization and vertical model extent_
	* __
	* _Set correct file to load for DEM input_
	* _Input file settings (data dimensions, order of data columns from input files)_
	* _Input file names (DEM, observed gravity signal, .. etc !?)_
2. Calculations: 
	* _Construct surface grid_
	* _Create gravity component grid_
	* __
	* _Plot all time series_
5. Run the entire script and look at output files

