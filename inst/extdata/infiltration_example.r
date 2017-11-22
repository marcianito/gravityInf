#!/usr/bin/Rscript

#############################################################################
### example script: hydrological modeling of infiltration into subsurface ###
#############################################################################

####################
## Instructions
####################
# Please modify or (better) copy this file to start a new reduction routine
# All lines in the SETUP section have to be filled, according to SG setup
# On completion, the script can be sourced
# or independently run within an R console
# Once it is finished, output files are stored in the set folders
##
# Optional: If a gravity data observation time series was supplied,
# data will be reduced by output results (directly by this script)
##
# For questions, comments or bugs,
# please visit https://github.com/marcianito/gravityInf
# or write to mreich@posteo.de
####################

## initial message
print("Welcome to gravityInf!")
print("This package will now run an infiltration scenario model, in order to find site dominant infiltration dynamics, based on your input.")

## developing
# library(HyGra)
library(devtools)
setwd("/home/mreich/Dokumente/written/ResearchRepos/")
load_all("UmbrellaEffect")
load_all("gravityInf")
# create package
# devtools::create("gravityInf")
# create docu
# library(roxygen2)
# setwd("/home/mreich/Dokumente/written/ResearchRepos/gravityInf")
# document()
# install package locally
# not working !? -> use load_all() above
# install()

####################
## load libraries
message("Loading necessary libraries..")
library(zoo); Sys.setenv(TZ = "GMT")
library(xts)
library(dplyr)
library(raster)
# library(UmbrellaEffect)
library(reshape2)
library(ggplot2)
library(viridis)
library(gstat)
library(ptinpoly)
#
library(plot3D)

library(grid)
library(gridExtra)
library(scales)
library(hydroGOF)
library(data.table)
library(ppso)
# for kriging
library(spacetime)
library(sp)
library(automap)
message("done.")
####################

#########################################
## SETUP
####################
message("Initializing setup..checking input data..")

## Output and input settings
# Directory
# path should be absolute
# (if not, it will be relative to the packages library path)
# use "test-data" for dir_input to use supplied example files within the package
# dir_input = "test-data"
dir_input = "~/temp/GI/"
dir_output = "~/temp/GI/"
# Output file type
# set to "csv", if output should also be saved as .csv (besides .rData)
output_type = ""
# Plotting option: should a plot of all time series be shown (and saved) in the end?
plot_data = TRUE

## Gravimeter location
# in [m]
# relativ, local coordinate sytem
SG_x = 0
SG_y = 0
SG_Z = 0
SG_SensorHeight = 1.05 
# UTM coordinate system
# SG_x = 4564082.00
# SG_y = 5445669.70
# SG_Z = 609.755
# SG_SensorHeight = 1.5 

## Model domain
# in [m]
# local grid or UTM, depending on the coordinates of the SG !
sprinklingArea_x = c(-7.5, 7.5) # min, max
sprinklingArea_y = c(-7.5, 7.5) # min, max
# grid3d_depth = c(-3, 0) # min, max
# UTM
# sprinklingArea_x = c(SG_x - 7.5, SG_x + 7.5) # min, max
# sprinklingArea_y = c(SG_y - 7.5, SG_y + 7.5) # min, max
# grid3d_depth = c(SG_Z, SG_Z - 3) # min, max

## Model discretization
# in [m]
grid3d_discr = data.frame(x = .5, y = .5, z = .1)
grid3d_depth = c(-3, 0) # min, max

## Boundaries of SG pillar
# please use same units as in DEM and model domain
# if pillar has the structure of a rectangular pillar
# local grid
Building_SGpillar_x = c(-1, 1) # min, max
Building_SGpillar_y = c(-1, 1) # min, max
Building_SGpillar_z = c(-1.2, 0) # min, max
# UTM
Building_SGpillar_x = c(SG_x - 1, SG_x + 1) # min, max
Building_SGpillar_y = c(SG_y - 1, SG_y + 1) # min, max
Building_SGpillar_z = c(SG_Z - 1.2, SG_Z) # min, max
# if pillar has the structure of a cylinder
# this is independent of local grid or UTM coordinates
# in [m]
thres_radius = 0.5
thres_depth = -1.2

## Input files
## general settings
# in case using .csv data, the special information has to supplied, in which columns the spatial information is stored
# the settings below are valid for 2d data files
# in the vector, the order is: x, y, z
spatial_col = c(1, 2, NA)
# in all cases, a column has to be specified, containing the observation data
# columns of observation data (in .csv and .rData files)
data_col = 3
# columns of observation data (in .tsf files)
data_tsf = 7 
# if the .csv has special characters for separating columns
# or decimal places, etc.
# the have to be EXPLICITLY specified in the read_data-function
# using sep = "??"
# using dec = "??"
# for further usage see ?read.csv

## DEM input file
# file name including its path
# should be absolute
# if left empty, a flat topographie will be assumed
DEM_input_file = ""
# DEM_input_file = "WE_UP_TO_300m_05m.asc"

## Water intensity distribution file
# IntensityDistribution_file = "waterIntensity_measured.csv"
IntensityDistribution_file = "waterIntensity_measured.rData"
# which interpolation algorithm should be used:
# inverse distance weight (IDW) or kriging (krige)
interpolation_method = "IDW"
# set percentage of how many zeros should be added at border / side of grid
# if not desired, set to 0
Zeros_border_density = .2

## Observed gravity data time series
# this is optional and can be left empty if no automatized reduction is desired
gravityObservations_file = "iGrav006_obs_60sec.tsf"
# gravityObservations_input_file = "iGrav006_obs_60sec.rData"
## Information for optimization algorithm
# use inverse or conversion mode
# if set to FALSE, a single conversion run of the infiltration model is exectuted
# in this case, infiltration parameters supplied in 'starting values' are used as model input
inverse = FALSE
# number of iterations of the algorithm
model_runs = 10

## set infiltration model parameters

# set scenarios to include in optimization routine
# use macro pore flow on top?
use_macro = TRUE
# 2 macro pore layers area only used if two_macro is set TRUE
two_macro = TRUE
# further infiltration dynamics
# possible options are: wetting front advancement (wfa), by-pass flow and perched water table
# in the optimization, input is treated and dealt with in terms of real numbers
# consequentially, discrete values as infiltration process descriptions have to be translated accordingly
# including or exlucing a process is therefore realized via adjusting the boundary values of the following vector
# Wfa = 1
# perched water table = 2
# by-pass flow = 3
# the example setting therefore includes both Wfa and perched water table scenarios
inf_dynamics_min = 1
inf_dynamics_max = 3

## modeling time (duration of sprinkling experiment)
# [min]
# precip_time = 360 
precip_time = 10 
## water input per timestep
# in [m³/min]
water_vol_min = 0.035

# set permitted error for mass balance
mb_permitted_error = 0.05

## Defines soil parameter boundaries
# min and max values, defining the search boundaries for the optimization algorithm
# Saturation deficit (dtheta)
# macropore layer
dtheta_macro_min = 0.05 #[VWC]
dtheta_macro_max = 0.20 #[VWC]
dtheta_macro2_min = 0.02 #[VWC]
dtheta_macro2_max = 0.10 #[VWC]
# Depth of processes
# in the case of no macro pore flow layers, the parameter 'mdepth' will be used for
# the depth of the single infiltration processes
mdepth_min = -0.5 #[m]
mdepth_max = -0.1 #[m]
mdepth2_min = -1.5 #[m]
mdepth2_max = -0.2 #[m]
# secondary infiltration process
# vertical bounaries
# if use_macro is set FALSE, this will be the only process used
dtheta_other_min = 0.02 #[VWC]
dtheta_other_max = 0.25 #[VWC]
# other infiltration processes ("pipe") are now spatially DIRECTLY connected below the macro pore layer
# if other is desired (e.g. gap between macro and pipe), this has to be activated again !
# with the following lines uncommented, a gap between macro and pipe layer is allowed
# break up criteria when pipedepth < mdepth is implemented in objective function (inf_model_3d)
# pipedepth_min = 0.2 #[m]
# pipedepth_max = 4.5 #[m]
pdepth_min = -2
pdepth_max = -1.0

# min max values for dividing water into horizontal / lateral parts (factor)
# in [%]
latflow_fac_min = 0.0 #[1]
latflow_fac_max = 1.0 #[1]

# Starting values of above set infiltration model parameters
dtheta_macro_start = 0.1
dtheta_macro2_start = 0.05
mdepth_start = -0.3
mdepth2_start = -1.0
dtheta_other_start = 0.05
# pipedepth_start = 0.4
latflow_fac_start = 0.4
# infiltration process
inf_dynamics_start = 2
# other process starting depth
pdepth_start = -1.5

## plotting options
plot_interval = precip_time / 10
plot_transect_loc = SG_y 

message("done.")
## end SETUP
#########################################

#########################################
## CALCULATIONS
####################
## nothing has to be changed from here on !!
message("Starting with calculation routine..")

# set working directory
if(dir_input == "test-data"){
  dir_wd = system.file("data", package="gravityInf")
  setwd(dir_wd)
}else{
  setwd(dir_input)
}

#########################################
## Gravimeter location
#########################################
SG_z = SG_Z + SG_SensorHeight
SGloc = data.frame(x=SG_x, y=SG_y, z=SG_z)

#########################################
## Generate cropped DEM and surface grid
#########################################
message("Generate cropped DEM and surface grid..")

surface_grid = surface_grid(
            DEM = DEM_input_file,
            grid_domain_x = sprinklingArea_x,
            grid_domain_y = sprinklingArea_y,
            grid_discretization = grid3d_discr,
            input_dir = dir_input,
            output_dir = dir_output
            # , sep = "a", etc.
)

if(!is.null(surface_grid)){
  if(plot_data){
    ggplot(surface_grid, 
           aes(x=x, y=y)) + 
           geom_tile(aes(fill = z))
  }
}

message("done.")
#########################################
# Generate 3d gravity component grid 
#########################################
message("Generate 3d gravity component grid..")

gravity_component_grid3d = gravity_comp_grid(
            surface = surface_grid,
            SG_coordinates = SGloc,
            grid_discretization = grid3d_discr,
            grid_depth = grid3d_depth,
            range_coords_x = sprinklingArea_x,
            range_coords_y = sprinklingArea_y
)

if(plot_data){
  ggplot(gravity_component_grid3d, 
         aes(x=x, y=y)) + 
         geom_tile(aes(fill = z))
         # geom_point(aes(color = z))
}

message("done.")
#########################################
## Correct gravity component grid for SG pillar 
#########################################
message("Removing SG pillar from gravity component grid..")

gravity_component_grid3d = correct_SGpillar(
            gravity_comp3d = gravity_component_grid3d,
            # Pillar_x = NA,
            # Pillar_y = NA,
            # Pillar_z = NA,
            correct_radius = thres_radius,
            correct_depth = thres_depth,
            SG_X = SG_x,
            SG_Y = SG_y #,
            # grid_discretization = NA
)

# save grid
save(gravity_component_grid3d, file = paste0(dir_output, "gravity_component_grid3d.rData"))

if(plot_data){
  message("Plotting transect of gravity component grid and saving plot to output directory..")
  plot_gcomp_grid(
                  grid_input = gravity_component_grid3d,
                  yloc = SG_y,
                  output_dir = dir_output,
                  grid_discretization = grid3d_discr
)
}

message("done.")
#########################################
## Create water intensity distribution grid
#########################################
message("Creating grid for water intensity distribution..")

Intensity_distribution_interpolated = create_WaterIntensityGrid(
            input_file = IntensityDistribution_file,
            intp_method = interpolation_method,
            surface = gravity_component_grid3d,
            zerosBorder = Zeros_border_density,
            spat_col = spatial_col,
            dat_col = data_col,
            UTM_gridcenter_x = SG_x,
            UTM_gridcenter_y = SG_y,
            input_dir = dir_input, 
            output_dir = dir_output
)

if(plot_data){
  ggplot(Intensity_distribution_interpolated, 
         aes(x=x, y=y)) + 
         geom_tile(aes(fill = intensity))
         # geom_point(aes(color = intensity))
}

# save grid
save(Intensity_distribution_interpolated, file = paste0(dir_output, "Intensity_distribution_interpolated.rData"))

message("done.")
#########################################
## Run infiltration model
#########################################

# combine data for config file
configfile = data.frame(dir_input,
                        dir_output,
                        precip_time,
                        water_vol_min,
                        IntensityDistribution_file = "Intensity_distribution_interpolated.rData",
                        gcomp_file = "gravity_component_grid3d.rData",
                        gravityObservations_file,
                        data_tsf,
                        spatial_col_x = spatial_col[1],
                        spatial_col_y = spatial_col[2],
                        spatial_col_z = spatial_col[3],
                        data_col,
                        discr_x = grid3d_discr$x,
                        discr_y = grid3d_discr$y,
                        discr_z = grid3d_discr$z,
                        mb_permitted_error,
                        use_macro,
                        two_macro,
                        model_runs,
                        plot_interval,
                        plot_transect_loc,
                        stringsAsFactors=FALSE)
save(configfile, file=paste0(dir_output, "configfile.rdata"))

## run model in conversion or inversion mode
print("Running infiltration model..")

if(!inverse){
  print("Model is run in conversion mode.")
  model_result = run_model_conversion(
              dtheta_macro = dtheta_macro_start,
              dtheta_macro2 = dtheta_macro2_start,
              mdepth = mdepth_start,
              mdepth2 = mdepth2_start,
              dtheta_other = dtheta_other_start,
              latflow_fac = latflow_fac_start,
              inf_dynamics = inf_dynamics_start,
              pdepth = pdepth_start,
              output_dir = dir_output
  )

}else{
## Run optimization algorithm
  print("Model is run in inversion mode.")
  print("Run optimization algorithm..")
  print("..this will take some time..")
  
  ## run model within optimization algorithm
  # for changing optimization function additional parameters, please see ?optim_dds
  model_result = run_model_inversion(
              dtheta_macro_min = dtheta_macro_min,
              dtheta_macro_max = dtheta_macro_max,
              dtheta_macro2_min = dtheta_macro2_min,
              dtheta_macro2_max = dtheta_macro2_max,
              dtheta_other_min = dtheta_other_min,
              dtheta_other_max = dtheta_other_max,
              mdepth_min = mdepth_min,
              mdepth_max = mdepth_max,
              mdepth2_min = mdepth2_min,
              mdepth2_max = mdepth2_max,
              latflow_fac_min = latflow_fac_min,
              latflow_fac_max = latflow_fac_max,
              inf_dynamics_min = inf_dynamics_min,
              inf_dynamics_max = inf_dynamics_max,
              pdepth_min = pdepth_min,
              pdepth_max = pdepth_max,
              dtheta_macro_start = dtheta_macro_start,
              dtheta_macro2_start = dtheta_macro2_start,
              mdepth_start = mdepth_start,
              mdepth2_start = mdepth2_start,
              dtheta_other_start = dtheta_other_start,
              latflow_fac_start = latflow_fac_start,
              inf_dynamics_start = inf_dynamics_start,
              pdepth_start = pdepth_start,
              input_dir = dir_input,
              output_dir = dir_output,
              inner_inum = 1
  )

  print("Finished optimization.")
# end of inversion / conversion mode infiltration model runs
}

# save model data (parameters in- and output)
save(model_result, file=paste0(dir_output, "Model_stats.rdata"))
write.table(model_result, file=paste0(dir_output, "Model_stats.csv"), sep="\t", dec=".", row.names = F, col.names = T, append = F)

print("done.")
#########################################
## Plot: Gravity reponse (observed and modeled)
#########################################

if(plot_data){
  message("Plotting modeled and observed gravity signal..")

if(!inverse){
  plot_gravity_responses(
              gravity_obs = gravityObservations_file,
              gravity_mod = paste0("model_output/GravityResponse_Infiltration_model_1.rData"),
              input_dir = dir_input,
              output_dir = dir_output
  )
}else{
# plot LAST (optimized) model run scenario
  plot_gravity_responses(
              gravity_obs = gravityObservations_file,
              gravity_mod = paste0("model_output/GravityResponse_Infiltration_model_", (n_param - 1), ".rData"),
              # gravity_mod = paste0("model_output/GravityResponse_Infiltration_model_9.rData"),
              input_dir = dir_input,
              output_dir = dir_output
  )
}
  message("done.")
}else{
  message("No plotting desired.")
}

#########################################
## Plot: model output along a 2d transect
#########################################

if(plot_data){
  message("Plotting 2d transect of modeled soil moisture data..")

if(!inverse){
  plot_transects_2d(
              soilmoisture_mod = paste0("model_output/Infiltration_model_output_1.rData"),
              plot_int = plot_interval,
              y_pos = SG_y,
              vert_limit = NA,
              output_dir = dir_output
  )
}else{
  plot_transects_2d(
              soilmoisture_mod = paste0("model_output/Infiltration_model_output_", (n_param - 1), ".rData"),
              # soilmoisture_mod = paste0("model_output/Infiltration_model_output_9.rData"),
              plot_int = plot_interval,
              y_pos = SG_y,
              vert_limit = NA,
              # input_dir = dir_input,
              output_dir = dir_output
  )
}
  message("done.")
}else{
  message("No plotting desired.")
}

# remove iteration parameter for inner optimization function calls
rm(n_param)

#########################################
## end CALCULATIONS
#########################################

message("All calculations have finished.")
message(paste0("Please have a look at the output file, located at: ", dir_output))

message("If you use this software in your publication, please cite this package.
        Information can be obtained using citation()")


##########
### prepare input TS data

# tt = load(file=paste0(dir_input, "igrav_raw_data_exp3.rdata"))
# igrav_exp3
# write.table(intensities_measured, file=paste0(dir_output, IntensityDistribution_file), row.names=F)

# sm = read_data(soilMoisture_input_file, dir_input)
# sm_ts = unique(sm$datetime)
# gg = read_data(gravityObservations_input_file, dir_input)
# gg_ts = unique(gg$datatime)
# sm = sm %>% dplyr::filter(datetime != 0)
# change_dates = data.frame(
#                           datetime = sm_ts[2:745],
#                           datatime = gg_ts)
# sm_mod = left_join(change_dates, sm)
## !! check if z has NEGATIVE COORDINATES DOWNWARDS !!
# SoilMoisture_input_1d = data.frame(
#                                    datetime = sm_mod$datatime,
#                                    z = sm_mod$z,
#                                    data = sm_mod$value)
# write.table(SoilMoisture_input_1d, file="SMdata_TS_1d.csv", row.names=F)
# save(SoilMoisture_input_1d, file="SMdata_TS_1d.rdata")
# getwd()

# tt = read_data(data_in = gravityObservations_input_file,
#                data_dir = dir_input,
#                dat_tsf = data_tsf
# )
