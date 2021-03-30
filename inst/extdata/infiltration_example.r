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
setwd("path/to/local/R/package/locations/")
load_all("HyGra")
load_all("gravityInf")
# alternatively, if installed via github, run the following 2 lines instead:
library(HyGra)
library(gravityInf)

# create package
# devtools::create("gravityInf")
# create docu
# library(roxygen2)
# setwd("/home/mreich/Documents/written/ResearchRepos/gravityInf")
# document()
# install package locally
# not working !? -> use load_all() above

####################
## load libraries
message("Loading necessary libraries..")
library(zoo); Sys.setenv(TZ = "GMT")
library(xts)
library(dplyr)
library(raster)
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
#
library(spacetime)
library(sp)
library(automap)
library(akima)
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
dir_input = "test-data"
# dir_input = ""
dir_output = "/path/to/output/storage/"
# Output file type
# set to "csv", if output should also be saved as .csv (besides .rData)
output_type = ""
# Plotting option: should a plot of all time series be shown (and saved) in the end?
plot_data = TRUE
#
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
# SG_SensorHeight = 1.05
#
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
#
## Model discretization
# in [m]
# grid3d_discr = data.frame(x = .5, y = .5, z = .1) # 10 runs; 1.68 mins
# grid3d_discr = data.frame(x = .25, y = .25, z = .1) # 10 runs; 5.30 mins
grid3d_discr = data.frame(x = .1, y = .1, z = .1) # 10 runs;  mins
# in REAL modeling so far: allDir = 0.1; depth up to 5 (?10) m;
grid3d_depth = c(-3, 0) # min, max
#
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
#
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
#
## DEM input file
# file name including its path
# should be absolute
# if left empty, a flat topographie will be assumed
DEM_input_file = ""
# DEM_input_file = "WE_UP_TO_300m_05m.asc"
#
## Water intensity distribution file
# IntensityDistribution_file = "waterIntensity_measured.csv"
IntensityDistribution_file = "waterIntensity_measured.rData"
# which interpolation algorithm should be used:
# inverse distance weight (IDW) or kriging (krige)
interpolation_method = "IDW"
# set percentage of how many zeros should be added at border / side of grid
# if not desired, set to 0
Zeros_border_density = .2
#
## Observed gravity data time series
# this is optional and can be left empty if no automatized reduction is desired
gravityObservations_file = "iGrav006_obs_60sec.tsf"
# gravityObservations_input_file = "iGrav006_obs_60sec.rData"
## Information for optimization algorithm
# use inverse or conversion mode
# if set to FALSE, a single conversion run of the infiltration model is exectuted
# in this case, infiltration parameters supplied in 'starting values' are used as model input
inverse = TRUE
# number of iterations of the algorithm
model_runs = 10
## objective function
# options are: KGE or gini or mNSE
# standard should be "mNSE"
objFunction = "mNSE"
# delete results and log-file of previous inverse model runs?
# recommended to turn this on, experienced problems with the algoritum to 
# deal with previous input
# but it should be possible to even continue and extend previous runs
delete_prevRuns = TRUE
#
## set infiltration model parameters
#
### model type (complexity)
# use lateral flow?
use_lateral_flow = T
# combine two infiltration processes ?
combine_processes = T
#
# set scenarios to include in optimization routine
# options are: wfa, bypass, perched, macropores
use_scenario = "macropores"
use_scenario2 = "bypass"
#
## modeling time (duration of sprinkling experiment)
# [min]
# precip_time = 360 
precip_time = 10
# precip_time = 3 
## water input per timestep
# in [m³/min]
# water_vol_min = 0.035
# water vol from experiment 3
# dir_Infinput = "/home/hydro/mreich/Irrigation/input/"
# exp3_meta = read.table(file=paste0(dir_Infinput, "Irrigation_precondition_dry"), skip= 5, nrows=8, dec=".", colClasses = character(), stringsAsFactors=F)
exp3_meta = read.table(file=paste0(dir_input, "Irrigation_precondition_dry"), skip= 5, nrows=8, dec=".", colClasses = character(), stringsAsFactors=F)
water_total_experiment = as.numeric(exp3_meta[6,2]) + as.numeric(exp3_meta[7,2]) # [m³]
water_vol_min = water_total_experiment / precip_time #[m³/min]
#
# set permitted error for mass balance
mb_permitted_error = 0.05
#
## Defines soil parameter boundaries
# min and max values, defining the search boundaries for the optimization algorithm
# Saturation deficit (dtheta)
# for selected process 1
dtheta_min = 0.05 #[VWC]
dtheta_max = 0.25 #[VWC]
# for selected process 2
dtheta2_min = 0.05 #[VWC]
dtheta2_max = 0.25 #[VWC]
# Depth of processes
# in the case of no macro pore flow layers, the parameter 'mdepth' will be used for
# the depth of the single infiltration processes
# for selected process 1
pdepth_min = -1.5 #[m]
pdepth_max = -0.1 #[m]
# for selected process 2
pdepth2_min = -4.9 #[m]
pdepth2_max = -1.5 #[m]
# min max values for dividing water into horizontal / lateral parts (factor)
# in [%]
latflow_fac_min = 0.0 #[1]
latflow_fac_max = 1.0 #[1]
#
# Starting values of above set infiltration model parameters
dtheta_start = 0.15
dtheta2_start = 0.15
pdepth_start = -0.5
pdepth2_start = -1.5
latflow_fac_start = 0.4
#
## plotting options
plot_interval = precip_time / 6
plot_transect_loc = SG_y 
#
#
message("done.")
## end SETUP
#########################################

#########################################
## CALCULATIONS
####################
## nothing has to be changed from here on !!
message("Starting with calculation routine..")
#
# set working directory
if(dir_input == "test-data"){
  dir_wd = system.file("data", package="gravityInf")
  setwd(dir_wd)
}else{
  setwd(dir_input)
}
#
#########################################
## Gravimeter location
#########################################
SG_z = SG_Z + SG_SensorHeight
SGloc = data.frame(x=SG_x, y=SG_y, z=SG_z)
#
#########################################
## Create Gravity Grid and correct for SG pillar 
#########################################
message("Generate 3d gravity component grid..")
gravity_component_grid3d = create_singleGravityGrid(
            DEM_input_file = DEM_input_file,
            DEM_dir = dir_input,
            SGloc = SGloc,
            grid_discretizations = grid3d_discr,
            grid_vertical = grid3d_depth,
            range_coords_x = sprinklingArea_x,
            range_coords_y = sprinklingArea_y,
            correct_SGpillar = c(thres_radius, thres_depth)
            )
#
# save grid
save(gravity_component_grid3d, file = paste0(dir_output, "gravity_component_grid3d.rData"))
#
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
#
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
#
if(plot_data){
  ggplot(Intensity_distribution_interpolated, 
         aes(x=x, y=y)) + 
         geom_tile(aes(fill = intensity))
         # geom_point(aes(color = intensity))
}
#
# save grid
save(Intensity_distribution_interpolated, file = paste0(dir_output, "Intensity_distribution_interpolated.rData"))

message("done.")
#########################################
## Run infiltration model
#########################################
#
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
                        use_scenario,
                        use_scenario2,
                        model_runs,
                        plot_interval,
                        plot_transect_loc,
                        objFunction,
                        stringsAsFactors=FALSE)
save(configfile, file=paste0(dir_output, "configfile.rdata"))
#
## run model in conversion or inversion mode
print("Running infiltration model..")
#
## use lateral flow module in model?
# if yes:
if(use_lateral_flow){
  print("Module for lateral flow was selected")
  if(combine_processes){
  print("Module for combined infiltration processes was selected")
  if(!inverse){
    print("Model is run in conversion mode.")
    model_result = run_model_conversion_combinedProcPP(
                dtheta = dtheta_start,
                latflow_fac = latflow_fac_start,
                pdepth = pdepth_start,
                output_dir = dir_output
                )
  #
  }else{
  ## Run optimization algorithm
    print("Model is run in inversion mode.")
    print("Run optimization algorithm..")
    print("..this will take some time..")
   # 
    ## run model within optimization algorithm
    # for changing optimization function additional parameters, please see ?optim_dds
    model_result = run_model_inversion_combinedProcPP(
                dtheta_min = dtheta_min,
                dtheta_max = dtheta_max,
                dtheta2_min = dtheta2_min,
                dtheta2_max = dtheta2_max,
                latflow_fac_min = latflow_fac_min,
                latflow_fac_max = latflow_fac_max,
                pdepth_min = pdepth_min,
                pdepth_max = pdepth_max,
                pdepth2_min = pdepth2_min,
                pdepth2_max = pdepth2_max,
                dtheta_start = dtheta_start,
                dtheta2_start = dtheta2_start,
                latflow_fac_start = latflow_fac_start,
                pdepth_start = pdepth_start,
                pdepth2_start = pdepth2_start,
                input_dir = dir_input,
                output_dir = dir_output,
                inner_inum = 1,
                del_prev = delete_prevRuns
    )
  ## run FINAL model, with optimized parameter values
    model_final = inf_model_3d_combinedProcLatFlow_final(
                                              param_1 = as.numeric(model_result$dtheta),
                                              param_2 = as.numeric(model_result$dtheta2),
                                              param_3 = as.numeric(model_result$pdepth),
                                              param_4 = as.numeric(model_result$pdepth2),
                                              param_5 = as.numeric(model_result$latflow)
                                                )
    print(paste0("final performance measure: ", model_final))
    message(paste0("final performance measure: ", model_final))
  }
  }else{
  print("Module for single infiltration process was selected")
  if(!inverse){
    print("Model is run in conversion mode.")
    model_result = run_model_conversion_singleProcPP(
                dtheta = dtheta_start,
                latflow_fac = latflow_fac_start,
                pdepth = pdepth_start,
                output_dir = dir_output
                )
  #
  }else{
  ## Run optimization algorithm
    print("Model is run in inversion mode.")
    print("Run optimization algorithm..")
    print("..this will take some time..")
   # 
    ## run model within optimization algorithm
    # for changing optimization function additional parameters, please see ?optim_dds
    model_result = run_model_inversion_singleProcPP(
                dtheta_min = dtheta_min,
                dtheta_max = dtheta_max,
                latflow_fac_min = latflow_fac_min,
                latflow_fac_max = latflow_fac_max,
                pdepth_min = pdepth_min,
                pdepth_max = pdepth_max,
                dtheta_start = dtheta_start,
                latflow_fac_start = latflow_fac_start,
                pdepth_start = pdepth_start,
                input_dir = dir_input,
                output_dir = dir_output,
                inner_inum = 1,
                del_prev = delete_prevRuns
    )
  ## run FINAL model, with optimized parameter values
    model_final = inf_model_3d_singleProcLatFlow_final(
                                              param_1 = as.numeric(model_result$dtheta),
                                              param_2 = as.numeric(model_result$pdepth),
                                              param_3 = as.numeric(model_result$latflow)
                                                )
    print(paste0("final performance measure: ", model_final))
    message(paste0("final performance measure: ", model_final))
  }
  }
}else{ # if not using lateral flow:
  print("Lateral flow is NOT considered")
  if(!inverse){
    print("Model is run in conversion mode.")
    model_result = run_model_conversion_singleProcPP(
                dtheta = dtheta_start,
                pdepth = pdepth_start,
                output_dir = dir_output
                )
  #
  }else{
  ## Run optimization algorithm
    print("Model is run in inversion mode.")
    print("Run optimization algorithm..")
    print("..this will take some time..")
   # 
    ## run model within optimization algorithm
    # for changing optimization function additional parameters, please see ?optim_dds
    # time keeping
    st = Sys.time()
    model_result = run_model_inversion_singleProcPP(
                dtheta_min = dtheta_min,
                dtheta_max = dtheta_max,
                pdepth_min = pdepth_min,
                pdepth_max = pdepth_max,
                dtheta_start = dtheta_start,
                pdepth_start = pdepth_start,
                input_dir = dir_input,
                output_dir = dir_output,
                inner_inum = 1,
                del_prev = delete_prevRuns
    )
  #
  en = Sys.time()
  time_elapsed = en - st
    print("Finished optimization.")
    print(time_elapsed)
  # end of inversion / conversion mode infiltration model runs
  ## run FINAL model, with optimized parameter values
    model_final = inf_model_3d_singleProc_final(
                                              param_1 = as.numeric(model_result$dtheta),
                                              param_2 = as.numeric(model_result$pdepth)
                                                )
    print(paste0("final performance measure: ", model_final))
    message(paste0("final performance measure: ", model_final))
  }
# end of decision: lateral flow module: yes or no
}
#
# save model data (parameters in- and output)
save(model_result, file=paste0(dir_output, "Model_stats.rdata"))
write.table(model_result, file=paste0(dir_output, "Model_stats.csv"), sep="\t", dec=".", row.names = F, col.names = T, append = F)

print("done.")
#########################################
## Plot: Gravity reponse (observed and modeled)
#########################################

if(plot_data){
  message("Plotting modeled and observed gravity signal..")
#
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
              gravity_mod = paste0("GravityResponse_Infiltration_model_final.rData"),
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
              soilmoisture_mod = paste0("Infiltration_model_output_final.rData"),
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
# rm(n_param)

#########################################
## end CALCULATIONS
#########################################

message("All calculations have finished.")
message(paste0("Please have a look at the output file, located at: ", dir_output))

message("If you use this software in your publication, please cite this package.
        Information can be obtained using citation()")


## close R
print("Closing R")
mpi.quit()
