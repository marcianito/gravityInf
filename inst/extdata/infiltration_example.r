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

## developing
# library(HyGra)
library(devtools)
setwd("/home/mreich/Dokumente/written/ResearchRepos/")
load_all("UmbrellaEffect")
load_all("gravityInf")
# create package
# devtools::create("gravityInf")
# create docu
library(roxygen2)
setwd("/home/mreich/Dokumente/written/ResearchRepos/gravityInf")
document()
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
SG_x = 3
SG_y = 3
SG_Z = 0
# UTM coordinate system
SG_x = 4564041.87 
SG_y = 5445662.88 
SG_Z = 606.471
SG_SensorHeight = 1.5 

## Model domain
# in [m]
# local grid or UTM, depending on the coordinates of the SG !
spinklingArea_x = c(0, 6) # min, max
spinklingArea_y = c(0, 6) # min, max
# grid3d_depth = c(-3, 0) # min, max
# UTM
spinklingArea_x = c(SG_x - 3, SG_x + 3) # min, max
spinklingArea_y = c(SG_y - 3, SG_y + 3) # min, max
# grid3d_depth = c(SG_Z, SG_Z - 3) # min, max

## Model discretization
# in [m]
grid3d_discr = data.frame(x = .5, y = .5, z = .5)
grid3d_depth = c(-3, 0) # min, max

## Boundaries of SG pillar
# please use same units as in DEM and model domain
# if pillar has the structure of a rectangular pillar
# local grid
Building_SGpillar_x = c(2, 4) # min, max
Building_SGpillar_y = c(2, 4) # min, max
Building_SGpillar_z = c(-1, 0) # min, max
# UTM
Building_SGpillar_x = c(SG_x - 1, SG_x + 1) # min, max
Building_SGpillar_y = c(SG_y - 1, SG_y + 1) # min, max
Building_SGpillar_z = c(SG_Z - 1, SG_Z) # min, max
# if pillar has the structure of a cylinder
# this is independent of local grid or UTM coordinates
# in [m]
thres_radius = 0.5
thres_depth = 1.2

## Input files
## general settings
# in case using .csv data, the special information has to supplied, in which columns the spatial information is stored
# the settings below are valid for 1d data files
# in the vector, the order is: x, y, z
spatial_col = c(NA, NA, 2)
# in all cases, a column has to be specified, containing the observation data
# columns of observation data (in .csv and .rData files)
data_col = 3
# columns of observation data (in .tsf files)
data_tsf = 13 
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
DEM_input_file = "WE_UP_TO_300m_05m.asc"

## Water intensity distribution file
IntensityDistribution_file = "FILENAME"
# which interpolation algorithm should be used:
# inverse distance weight (IDW) or kriging (krige)
interpolation_method = "IDW"
# set percentage of how many zeros should be added at border / side of grid
# if not desired, set to 0
Zeros_border_density = 6

## Observed gravity data time series
# this is optional and can be left empty if no automatized reduction is desired
gravityObservations_input_file = "SG030_TS_1month.tsf"

##########################################
## set infiltration model parameters
##########################################

# ## single
# # scen = "macro"
# # scen = "wfa"
# # scen = "bypass"
# # scen = "perched"
# ## combined
# scen = "macro_bypass"
# # scen = "macro_wfa"
# # scen = "macro_perched"

## modeling time (duration of sprinkling experiment)
# [min]
precip_time = 360 
## water input per timestep
# in [m³/min]
water_vol_min = 0.035

# set permitted error for mass balance
mb_permitted_error = 0.05

## Defines soil parameter boundaries
# min and max values, defining the search boundaries for the optimization algorithm
# Saturation deficit (theta)
# use macro pore flow on top?
use_macro = TRUE
# 2 macro pore layers area only used if two_macro is set TRUE
two_macro = TRUE
# macropore layer
dtheta_macro_min = 0.05 #[VWC]
dtheta_macro_max = 0.20 #[VWC]
dtheta_macro2_min = 0.02 #[VWC]
dtheta_macro2_max = 0.10 #[VWC]
mdepth_min = 0.1 #[m]
mdepth_max = 0.5 #[m]
mdepth2_min = 0.2 #[m]
mdepth2_max = 1.5 #[m]
# secondary infiltration process
# vertical bounaries
# if use_macro is set FALSE, this will be the only process used
dtheta_pipe_min = 0.2 #[VWC]
dtheta_pipe_max = 0.25 #[VWC]
# other infiltration processes ("pipe") are now spatially DIRECTLY connected below the macro pore layer
# if other is desired (e.g. gap between macro and pipe), this has to be activated again !
# with the following lines uncommented, a gap between macro and pipe layer is allowed
# break up criteria when pipedepth < mdepth is implemented in objective function (inf_model_3d)
# pipedepth_min = 0.2 #[m]
# pipedepth_max = 4.5 #[m]

# min max values for dividing water into horizontal / lateral parts (factor)
# in [%]
latflow_fac_min = 0.0 #[1]
latflow_fac_max = 1.0 #[1]

# Starting values of above set infiltration model parameters
dtheta_macro_start = 0.1
dtheta_macro2_start = 0.05
mdepth_start = 0.3
mdepth2_start = 1.0
dtheta_pipe_start = 0.05
# pipedepth_start = 0.4
latflow_fac_start = 0.4

message("done.")
## end SETUP
#########################################

#########################################
## CALCULATIONS
####################
## nothing has to be changed from here on !!
message("Starting with calculation routine..")

# set working directory
dir_wd = system.file("data", package="gravityInf")
setwd(dir_wd)

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

message("done.")
#########################################
# Generate 3d gravity component grid 
#########################################
message("Generate 3d gravity component grid..")

gravity_component_grid3d = gravity_comp_grid(
            surface = surface_grid,
            SG_coordinates = SGloc,
            grid_discretization = grid3d_discr,
            grid_depth = grid3d_depth
)

message("done.")
#########################################
## Correct gravity component grid for SG pillar 
#########################################
message("Removing SG pillar from gravity component grid..")

gravity_component_grid3d = correct_SGpillar(
            gcompgrid = gravity_component_grid3d,

)

if(circular){
  gcomp_irrigation_domain = dplyr::mutate(gcomp_irrigation_domain, gcomp = ifelse((x - igrav_x)^2 + (y - igrav_y)^2 < thres_radius & zgrid <= thres_depth,0,gcomp))
}else{
}


save(gravity_component_grid3d, file = paste0(dir_output, "gravity_component_grid3d.rData"))

message("done.")
#########################################
## Create water intensity distribution grid
#########################################
message("Creating grid for water intensity distribution..")

create_WaterIntensityGrid(
            input_file = IntensityDistribution_file,
            intp_method = interpolation_method,
            surface_grid = surface_grid,
            zerosBorder = Zeros_border_density,
            output_dir = dir_output,
            ...
)



# load file: IntensityDistribution_file as Iintensities


Iintensities_addZeros = select(Iintensities, xrel, yrel, total_weight_dif)
# number of "zeros" per border-side
nzeros_border = 6
Zeros = data.frame(xrel = c(seq(min(grid_surface$xrel),max(grid_surface$xrel),length.out=nzeros_border),seq(min(grid_surface$xrel),max(grid_surface$xrel),length.out=nzeros_border), rep(min(grid_surface$xrel),nzeros_border),rep(max(grid_surface$xrel),nzeros_border)),
                   yrel = c(rep(min(grid_surface$yrel),nzeros_border),rep(max(grid_surface$yrel),nzeros_border),seq(min(grid_surface$yrel),max(grid_surface$yrel),length.out=nzeros_border),seq(min(grid_surface$yrel),max(grid_surface$yrel),length.out=nzeros_border)),
                   total_weight_dif = 0
                   )
Iintensities_Zeros = rbind(Iintensities_addZeros, Zeros)

## IWD
# # without extra added Zeros
# idw.gstat = gstat(formula = total_weight_dif ~ 1, locations = ~ xrel + yrel, data = Iintensities, nmax = 10, set = list(idp = 2))
# with extra Zeros (as stützstellen at the borders of the irrigation area)
idw.gstat = gstat(formula = total_weight_dif ~ 1, locations = ~ xrel + yrel, data = Iintensities_Zeros, nmax = 10, set = list(idp = 2))
#idw.gstat = gstat(formula = dataraw ~ 1, locations = ~ x + y + Depth, data = data_in, nmax = 10, nmin=3, maxdist=1.1, set = list(idp = 2))
data_interpolated <- predict(idw.gstat, grid_surface)[,-4]
colnames(data_interpolated)[3] = "weight_intp"

## calculate intensities per min out of total_weight_dif
# Becher_radius = 0.035 #[m]
num_cell = length(grid_surface$x)
TotalWaterOnBecher = sum(data_interpolated$weight_intp, na.rm=T)
weight_avg = TotalWaterOnBecher / num_cell
# intensity as relatin of:
# cell weight / average cell weight
Iintensities_interpolated_IDW = dplyr::mutate(data_interpolated, intensity = weight_intp / weight_avg)

# check sum of weights after interpolation
# should be equal to number of cells
sum(Iintensities_interpolated_IDW$intensity, na.rm=T)

## kriging

# load irrigation intensities
load(file="../../output/irrigation/Iintesities_spatialPointDistribution.rdata")

grid_surface_krige = cbind(grid_surface, dummy=2, fummy=5)
gridded(grid_surface_krige) = ~xrel + yrel
coordinates(Iintensities)= ~xrel + yrel
#hydrus domain grid
vario = autofitVariogram(total_weight_dif ~ 1, Iintensities,
	      model = c("Sph", "Exp", "Gau", "Ste", "Mat"),
              #model = c("Mat"),
	      verbose=T)
# load irrigation intensities WITH Zeros added
# one has to be used for variogram (normal raw data)
# one has to be used for modeling (normal data + Zero data)
load(file="../../output/irrigation/Iintesities_spatialPointDistribution_withZeros.rdata")
coordinates(Iintensities_Zeros)= ~xrel + yrel
#interpolate data to new grid
data_pred = krige(total_weight_dif ~ 1, Iintensities_Zeros, grid_surface_krige, vario$var_model)

#####
var <- variogram(total_weight_dif~1,data=Iintensities,assumeRegular=F,na.omit=T) 
vario_krige = fit.variogram(var, vario$var_model)
vario_krige = vgm(340677.1, "Mat", 12.29943, 1.7)

data_pred = krige(total_weight_dif ~ 1, Iintensities_Zeros, grid_surface_krige, vario_krige)
data_pred = krige(total_weight_dif ~ 1, Iintensities, grid_surface_krige, vario_krige)
plot(data_pred)
####

# construct data.frame
data_interpolated = data.frame(xrel = coordinates(data_pred)[,1], yrel = coordinates(data_pred)[,2], weight_intp = data_pred$var1.pred)

## calculate intensities per min out of total_weight_dif
Becher_radius = 0.035 #[m]
TotalWaterOnBecher = sum(data_interpolated$weight_intp, na.rm=T)

# set negative values to 0 (zero)
data_interpolated$weight_intp[which(data_interpolated$weight_intp < 0)] = 0
# check sum of weights after interpolation
sum(data_interpolated$weight_intp, na.rm=T)

# intensity in mm/min
Iintensities_interpolated_krige = mutate(data_interpolated, intensity = weight_intp/1000 / (pi*Becher_radius^2)) %>%
		mutate(intensity_percent = weight_intp / TotalWaterOnBecher)


# save data
save(Intensity_distribution_interpolated, file = paste0(dir_output, "Intensity_distribution_interpolated.rData"))

message("done.")
#########################################
## 
#########################################

#########################################
## 
#########################################

##########################################



##########################################
## Run optimization algorithm
##########################################
print("Run optimization algorithm..")
print("..this will take some time..")

# combine data for config file
configfile = data.frame(dir_input,
                        dir_output,
                        precip_time,
                        water_vol_min,
                        IntensityDistribution_file,
                        "gravity_component_grid3d.rData",
                        gravity_observations_file,
                        mb_permitted_error,
                        stringsAsFactors=FALSE)
save(configfile, file=paste0(dir_output, "configfile.rdata"))

## run model within optimization algorithm
# for changing optimization function additional parameters, please see ?optim_dds
opt_result = run_model(
            dtheta_macro_min = dtheta_macro_min,
            dtheta_macro_max = dtheta_macro_max,
            dtheta_macro2_min = dtheta_macro2_min,
            dtheta_macro2_max = dtheta_macro2_max,
            dtheta_pipe_min = dtheta_pipe_min,
            dtheta_pipe_max = dtheta_pipe_max,
            mdepth_min = mdepth_min,
            mdepth_max = mdepth_max,
            mdepth2_min = mdepth2_min,
            mdepth2_max = mdepth2_max,
            latflow_fac_min = latflow_fac_min,
            latflow_fac_max = latflow_fac_max,
            dtheta_macro_start = dtheta_macro_start,
            dtheta_macro2_start = dtheta_macro2_start,
            mdepth_start = mdepth_start,
            mdepth2_start = mdepth2_start,
            dtheta_pipe_start = dtheta_pipe_start,
            latflow_fac_start = latflow_fac_start,
            macopores = use_macro,
            macro2 = two_macro,
            ...
)

# prepare model input data and parameter boundaries
if(macropores){
  if(macro2){
    # combine input parameters
    param_bounds = data.frame(minimum = c(dtheta_macro_min, dtheta_macro2_min, dtheta_pipe_min, mdepth_min, mdepth2_min, latflow_fac_min),
                          maximum = c(dtheta_macro_max, dtheta_macro2_max, dtheta_pipe_max, mdepth_max, mdepth2_max, latflow_fac_max))
    # combine start parameter values
    param_startvalue = c(dtheta_macro_start, dtheta_macro2_start, dtheta_pipe_start, mdepth_start, mdepth2_start, latflow_fac_start)
    # set name of model to use
    model = "inf_model_3d_2macro"
  }else{
    # combine input parameters
    param_bounds = data.frame(minimum = c(dtheta_macro_min, dtheta_pipe_min, mdepth_min, latflow_fac_min),
                          maximum = c(dtheta_macro_max, dtheta_pipe_max, mdepth_max, latflow_fac_max))
    # combine start parameter values
    param_startvalue = c(dtheta_macro_start, dtheta_pipe_start, mdepth_start, latflow_fac_start)
    # set name of model to use
    model = "inf_model_3d_macro"
  }
}else{
  # combine input parameters
  param_bounds = data.frame(minimum = c(dtheta_pipe_min, mdepth_min, latflow_fac_min),
                        maximum = c(dtheta_pipe_max, mdepth_max, latflow_fac_max))
  # combine start parameter values
  param_startvalue = c(dtheta_pipe_start, mdepth_start, latflow_fac_start)
  # set name of model to use
  model = "inf_model_3d_single"
}

# set counting parameter
# this is used for naming figures and files for individiual model runs within the optimization routine
n_param = 1 
# n_param <<- n_param
# param_bounds <<- param_bounds
# param_startvalue <<- param_startvalue

# set working directory
# setwd(dir_output)

## run optimization
opt_result = optim_dds(
    # objective_function = inf_model_3d_2macro, #set model to use, according to input parameters supplied
    objective_function = get(model), #set model to use, according to input parameters supplied
    number_of_parameters = length(param_bounds[,1]),
    number_of_particles =  1,
    max_number_function_calls= 200,
    r=0.2,
    abstol = -Inf,
    reltol = -Inf,
    max_wait_iterations=50,
    parameter_bounds = param_bounds,
    initial_estimates = param_startvalue,
    lhc_init=FALSE,
    do_plot = NULL,
    wait_for_keystroke = FALSE,
    # logfile=  ppso_log,
    # projectfile = "projectfile_to_resume",
    load_projectfile = "try",
    break_file=NULL,
    plot_progress=FALSE,
    tryCall=FALSE)

return(opt_result)



##########################################
## save results
##########################################
print("combine and save results and run information")
# scenario information
model_info = list(dir_input = dir_input,
                     dir_output = dir_output,
                     duration = precip_time,
                     total_water_volume = water_vol_min,
                     water_distribution_file = IntensityDistribution_file,
                     gravity_component_grid_file = "gravity_component_grid3d.rData",
                     gravity_observations_file = gravity_observations_file,
                     permitted_massbalance_error = mb_permitted_error)

# combine optimization with scenario information
stats = cbind(model_info, opt_result)

save(stats, file=paste0(dir_output, "OptModel_stats.rdata"))
write.table(stats, file=paste0(dir_output, "OptModel_stats.csv"), sep="\t", dec=".", row.names = F, col.names = T, append = T)

print("Finished optimization.")
#########################################
## Plot: Gravity reponse (observed and modeled)
#########################################

if(plot_data){
  message("Plotting modeled and observed gravity signal..")


load(file=paste0(dir_scen, "modOptimized_2macro/gmod_Irrigation_macropiping_27.rdata"))

# load igrav time series in same period
load(file=paste0(dir_output, "gravity/igrav_raw_data_exp3.rdata"))
# set same column name for joining datasets
colnames(igrav_exp3)[2] = "gmod"
igrav_exp3_cor = mutate(igrav_exp3, gmod = gmod - min(gmod))

igrav_timesteps = data.frame(Timestep = 1:length(igrav_exp3_cor$gmod[-1]), gmod = igrav_exp3_cor$gmod[-1])

## create measurement uncertainty of +- 10 %
# this will be due to unprecise estimated total water mass used in experiment
gsignal_irrigation_macropiping_over = gsignal_irrigation_macropiping
gsignal_irrigation_macropiping_over$gmod = gsignal_irrigation_macropiping_over$gmod * 1.1
gsignal_irrigation_macropiping_under = gsignal_irrigation_macropiping
gsignal_irrigation_macropiping_under$gmod = gsignal_irrigation_macropiping_under$gmod * 0.9 

## combine datasets
gmod = rbind(
         # cbind(igrav_exp3_cor, Scenario="iGrav (observed)"),
         cbind(igrav_timesteps, Scenario="Observed gravity response"),
	     cbind(gsignal_irrigation_macropiping, Scenario="Model: Macropores & By-pass flow"),
	     cbind(gsignal_irrigation_macropiping_over, Scenario="uncertainty"),
	     cbind(gsignal_irrigation_macropiping_under, Scenario="uncertainty")
         )
gmod$Scenario = factor(gmod$Scenario, levels=c("uncertainty","Observed gravity response", "Model: Macropores & By-pass flow"))

# png(file=paste0(dir_output, "gravity/plots/Irrigation_complexScenarios_ItensityIdeal_dtheta005_mitiGrav.png"), width=1500, height=1000, res=250)
# png(file=paste0(dir_output, "gravity/plots/Irrigation_combinedScenarios_macropipe_pipedepth150cm_macrodepth50cm_DynamiDdelays_dtehta010_IReal.png"), width=1500, height=1000, res=250)
png(file=paste0(dir_praesi, "Irrigation_combinedExtScenarios.png"), width=1500, height=1800, res=250, bg="transparent")
ggplot(gmod, aes(x=Timestep, y=gmod, colour=Scenario)) + geom_line(size=1.5) + 
	ylab("Gravity [nm/s²]") + xlab("Time since start of experiment [min]") +
    scale_color_manual(values = c("lightgrey","red","blue"), breaks=c("Observed gravity response", "Model: Macropores & By-pass flow")) + 
    figure_style + theme(
     legend.position = "bottom",
     legend.title = element_blank(),
	 legend.text=element_text(size=18),
	 panel.background = element_rect(fill="transparent"),
     panel.grid.major = element_line(colour = "black", linetype = "dotted")
     ) +
    guides(colour = guide_legend(nrow = 2))
dev.off()

  message("done.")
}else{
  message("No plotting desired.")
}

#########################################
## Plot: model output along a 2d transect
#########################################

if(plot_data){
  message("Plotting 2d transect of modeled soil moisture data..")


filename = "Irrigation_combiExt_macropiping_161"
modeloutput = load(file=paste0(dir_input, filename,".rdata"))

# ERT time steps
ERT_exp3 = load(file=paste0(dir_output, "ERT/ERT_profB_exp3_spatiallyAdjusted.rdata"))
ERT_exp3 = load(file=paste0(dir_output, "ERT/ERT_profB_exp3_spatiallyAdjusted_neg.rdata"))
ERT_exp3 = get(ERT_exp3)
ERT_exp3$datetime = as.POSIXct(trunc(ERT_exp3$datetime, units = "mins"))

data_obs_ts = data.frame(datetime = unique(ERT_exp3$datetime), Timestep = unique(ERT_exp3$Timestep))
# limit to first XX timesteps
# 1 timestep is approx. 1 hour
# ts = 0 is reference state BEFORE experiment
tslimit = 6
data_obs_ts = dplyr::filter(data_obs_ts, Timestep <= tslimit) %>%
              dplyr::mutate(Timestep = Timestep * 60) #multiply with 60 minutes

## load sensor coordinates
# in Gauss-Krueger
sensors_coords = read.table(file=paste0(dir_input_std, "DEM/sensor_coords_GK.csv"), sep=",", dec =".", header=T)
ERT_coords = dplyr::filter(sensors_coords, name =="ERT_B") %>%
             dplyr::filter(distance >= (max(distance) - 8)) %>%
             dplyr::select(-org_fid, -id) %>%
             dplyr::rename(x = distance)

SM_coords = dplyr::filter(sensors_coords, name =="SM_interpolated") %>%
            dplyr::filter(distance <= 8) %>%
            dplyr::mutate(x = rev(distance)) %>%
            dplyr::select(-org_fid, -id, -distance)
# save
save(ERT_coords, file=paste0(dir_output,"irrigation/modeloutput_profiles/ERT_exp3_withGKcoords.rdata"))
save(SM_coords, file=paste0(dir_output,"irrigation/modeloutput_profiles/SM_exp3_withGKcoords.rdata"))

#########################################
## limit data to ERT and SM observation time steps
#########################################
# modeloutput = dplyr::inner_join(data_obs_ts, get(modeloutput))
modeloutput = dplyr::inner_join(data_obs_ts, modeloutput)

#########################################
## round X and Y data (coordinates) in order to join datasets
#########################################
modeloutput_round = dplyr::mutate(modeloutput, x = round(x,1)) %>%
         dplyr::mutate(y = round(y,1))

ERT_coords_round = dplyr::rename(ERT_coords, xrel = x) %>%
        dplyr::mutate(x = round(X,1)) %>%
        dplyr::mutate(y = round(Y,1)) %>%
        dplyr::select(-X,-Y)

SM_coords_round = dplyr::rename(SM_coords, xrel = x) %>%
        dplyr::mutate(x = round(X,1)) %>%
        dplyr::mutate(y = round(Y,1)) %>%
        dplyr::select(-X,-Y)

#########################################
## limit data to SM profile
#########################################

modeloutput_SMprofile = dplyr::left_join(SM_coords_round, modeloutput_round)

#########################################
## limit data to ERT profile
#########################################

modeloutput_ERTprofile = dplyr::left_join(ERT_coords_round, modeloutput_round)

# ERT: first 300 cm
model_exp3_300cm = dplyr::filter(modeloutput_ERTprofile, zgrid <= 3.0)
# SM: first 40 cm
    scale_x_continuous(breaks=c(10,12,14,16,18)) +
    scale_fill_gradientn(breaks = c(0.01,0.05,0.1,0.2), colours=rev(viridis(7)), na.value="red") +
	# scale_fill_gradientn(colours=rev(viridis(7)), na.value="red") +
    facet_grid(Timestep ~ ., scale="free") +
	#ylim(3.1,0) + 
    # ylab("Depth [m]") + xlab("Profile length [m]") + 
    ylab("Depth [m]") + xlab("Profile x [m]") + 
    # labs(fill = expression(Delta * "Resistivity [µS/cm]")) +
    labs(fill = expression(Delta * "Soil moisture [%VWC]")) +
    ggtitle("Model output: Soil moisture") +
    figure_style + 
    theme(legend.position ="bottom",
	      legend.text=element_text(size=17),
	      legend.title=element_text(size=19),
	      panel.grid.major = element_line(colour = "black", linetype = "dotted"),
	      panel.grid.minor = element_line(colour = "black", linetype = "dotted"))

model_profiles_SM.gg = ggplot(model_exp3_40cm, aes(x=xrel, y=zgrid)) +
    geom_tile(aes(fill = value, width = .25)) + 
    scale_y_continuous(breaks=c(.4,.2,0), trans="reverse") + 
    scale_x_continuous(breaks=c(0,2,4,6,8)) +
    scale_fill_gradientn(breaks = c(0.01,0.05,0.1,0.2), colours=rev(viridis(7)), na.value="red") +
    facet_grid(Timestep ~ ., scale="free") +
	#ylim(3.1,0) + 
    # ylab("Depth [m]") + xlab("Profile length [m]") + 
    ylab("Depth [m]") + xlab("Profile x [m]") + 
    labs(fill = expression(Delta * "Soil moisture [%VWC]")) +
    ggtitle("Model output: Soil moisture") +
    figure_style + 
    theme(legend.position ="bottom",
	      legend.text=element_text(size=17),
	      legend.title=element_text(size=19),
	      panel.grid.major = element_line(colour = "black", linetype = "dotted"),
	      panel.grid.minor = element_line(colour = "black", linetype = "dotted"))


  message("done.")
}else{
  message("No plotting desired.")
}

#########################################
## end CALCULATIONS
#########################################

message("All calculations have finished.")
message("Please have a look at the output file, located at: ")
message(dir_output)

message("")


##########
### prepare input TS data

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
