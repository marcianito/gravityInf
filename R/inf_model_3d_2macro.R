#' @title Hydro-gravimetric model (incl. lateral flow), for 2 different macro pore layers
#'
#' @description 3D distribution algorithm, modeling the soil infiltration of water from a sprinkling experiment
#'
#' @param paramvec Numeric vector, proving input in this EXACT ORDER of the parameter values: saturation deficit of 1st process (macro pores), saturation deficit of 2nd process, total vertical thickness (extent) of 1st process (macropores), scaling factor (0 to 1) defining the relatinship between lateral and vertical flow after saturation of a cell.
#' @param test
#' @param test
#' ...
#' @details The function was written to be easily usable together with optimization tools.
#' It thus  needs a configfile.rdata-file with the following columns, each consisting of one single value:
#' dir_input, dir_output, precip_time, IntensityDistributionution, water_vol_min, gcompfile, gravityObs, mb_permitted_error.

#' @references Marvin Reich (2017), mreich@@posteo.de
#' @examples missing

inf_model_3d_2macro = function(param_vec){

# ####################
# ## debugging
# param_vec = c(dtheta_macro_start, dtheta_macro2_start, dtheta_other_start, mdepth_start, mdepth2_start, latflow_fac_start, inf_dynamics_start, pdepth_start)
# ####################

# ##########
# ## logging to file
# logfile = paste0("irrigation zgrid: ",unique(Irrigation_in$zgrid))
# write.table(logfile, file=paste0(dir_output, "logfile"), append=T)
print(paste0("n_param: ", n_param))
# ##########

###################
## load model configuration
load(file="configfile.rdata")
input_dir = configfile$dir_input
output_dir = configfile$dir_output
precip_time = configfile$precip_time
IntensityDistribution_file = configfile$IntensityDistribution_file
water_vol_min = configfile$water_vol_min
gcomp_file = configfile$gcomp_file
gravity_obs = configfile$gravityObservations_file
dat_tsf = configfile$data_tsf
spat_col_x = configfile$spatial_col_x
spat_col_y = configfile$spatial_col_y
spat_col_z = configfile$spatial_col_z
dat_col = configfile$data_col
x = configfile$discr_x
y = configfile$discr_y
z = configfile$discr_z
mb_permitted_error = configfile$mb_permitted_error
macropore_layer = configfile$use_macro
macropore_layer2 = configfile$two_macro
n_iterations = configfile$model_runs
plot_interval = configfile$plot_interval
plot_transect_2d = configfile$plot_transect_loc

# construct grid discretization
grid_discretization = data.frame(x, y, z)
# construct spatial input data columns
spat_col = c(spat_col_x, spat_col_y, spat_col_z)
# set plot directory
dir_plots = paste0(output_dir,"model_output/plots/")

# volume of grid cell 
vol_cell = grid_discretization$x * grid_discretization$y * grid_discretization$z # [m³]

## load gravity calculation data
# load gravity effect grid
gcomp_grid = load(file=paste0(input_dir, gcomp_file))
gcomp_grid = get(gcomp_grid)

# load igrav time series in same period
gravity_obs_data = read_data(gravity_obs, input_dir)
# set same column name for joining datasets
colnames(gravity_obs_data)[2] = "gmod"
# reduce time series by first value (in order to start at zero mass change)
# gravity_obs_data = mutate(gravity_obs_data, gmod = gmod - min(gmod))
gravity_obs_data = mutate(gravity_obs_data, gmod = gmod - gmod[1])
# gravity_timesteps = data.frame(datetime = 1:length(gravity_obs_data$gmod[-1]), gmod = gravity_obs_data$gmod[-1])
gravity_timesteps = data.frame(datetime = 1:precip_time, gmod = gravity_obs_data$gmod[2:(precip_time+1)])

Irrigation_grid = dplyr::select(gcomp_grid, x, y, z, Depth)

# vertical layers
zlayers = round(unique(Irrigation_grid$Depth),1)
# number of cells per layer
num_cell = length(unique(Irrigation_grid$x)) *  length(unique(Irrigation_grid$y))

# define factors for scaling between "side" and "corner" cell neighbours
# following (Quinn 1997), where cardinal (side) and diagnol (corner) factors are summed
# all 8 cells are delivered thus leading to the following factors
# for side cells
# 4 * 0.5/SUM(all factors for all cells) = 4 * 0.15
sfactor = 0.6
# for corner cells
# 4 * 0.35/SUM(all factors for all cells) = 4 * 0.1
cfactor = 0.4

# ## logging to file
# logfile = paste0("irrigation Depth: ",unique(Irrigation_in$Depth))
# write.table(logfile, file=paste0(dir_output, "logfile"), append=T)
# ##########

# load water intensity distribution
Intensity_distribution = load(file=paste0(output_dir, IntensityDistribution_file))
Intensity_distribution = get(Intensity_distribution)

## prepare statistis outout
stats = data.frame()

print("Finished loading setup parameters.")

##########################################
## pass parameters from DDS-function param_vec space
##########################################
dtheta_macro = param_vec[1]
dtheta_macro2 = param_vec[2]
dtheta_other = param_vec[3]
# parameter values in meters have to be rounded
# otherwise values will NOT match grid discretization !!
mdepth = round(param_vec[4],1)
mdepth2 = round(param_vec[5],1)
# somewhere define that pipedepth > mdepth !!
# now pipedepth is directly below the macro pore layer
# pipedepth = round(mdepth + 0.1, 1)
# pipedepth is set individually
# this allows a space between macro pore layer and 2nd process
# pipedepth = round(param_vec[6],1)
# lateral flow factor (seperate into lateral and vertical flow)
latflow_fac = param_vec[6]
vertflow_fac = 1 - latflow_fac
# infiltration dynamics
# round to integer values !!
# in order to have discrete values for the individual process
inf_dyn = round(param_vec[7])
# depth of other infiltration process
# accounts / valid for: by-pass flow and perched water table
pdepth = round(param_vec[8],1)

####################
## validity of some assumptions

if(mdepth2 >= mdepth | pdepth >= mdepth2){
  print(paste0("Problem with vertical layer distribution in paramterset: ",n_param))
  n_param <<- n_param + 1
  kge_fit = 1
  return(kge_fit) 
# if everything is okay, run normally
}else{ 

##########################################
## build up model domain space

# depth of macropore layer thickness
# mdepth = .3 # [m]
# determine vertical start
mdepth_layer = which(zlayers == mdepth)
mdepth2_layer = which(zlayers == mdepth2)
# print(mdepth_layer)
# print(zlayers)
# print(dtheta_macro)
# print(dtheta_infDynProc)

# # depth of pipe routing
# otherdepth = mdepth2 - grid_discretization$z # [m]
# depth of underlying process (wfa, by-pass, perched)
pdepth_layer = which(zlayers == pdepth)
# determine vertical start
# other_layer = which(zlayers == otherdepth)
## how many layers are between macro and other process?
# layer_between = other_layer - mdepth_layer - 1
macro_layer_between = mdepth2_layer - mdepth_layer
layer_between_processes = pdepth_layer -  mdepth2_layer
# ####################
# 
tsx = dplyr::mutate(Irrigation_grid,
      # column = rep(gridcells$ncell, length(zlayers))) %>%
      column = rep(1:num_cell, length(zlayers))) %>%
      dplyr::mutate(cnt = 1) %>%
      dplyr::mutate(datetime = 0) %>%
      dplyr::mutate(value = 0) %>%
      dplyr::mutate(prevalue = 0)

## tag vertical layers with infiltration process

## check which infiltration dynamic is beeing applied right now
# and respectively adjust the filling of cells
## availlable are:
# these are (and have to be) hard-coded !!
# Wfa = 1
# perched water table = 2
# by-pass flow = 3

## Wfa
if(inf_dyn == 1){
layer_params = data.frame(Depth = zlayers,
                         infProcess = c(rep("macro", mdepth_layer),rep("macro2", macro_layer_between), rep("other",(length(zlayers) - mdepth2_layer))),
                         nlayer = c(rep(1,mdepth_layer), rep(1, macro_layer_between), seq(1,length.out=(length(zlayers) - mdepth2_layer))),
                         dtheta = c(rep(dtheta_macro,mdepth_layer),rep(dtheta_macro2,macro_layer_between),rep(dtheta_other,(length(zlayers) - mdepth2_layer)))
              )
}
## perched water table
if(inf_dyn == 2){
layer_params = data.frame(Depth = zlayers,
                         infProcess = c(rep("macro", mdepth_layer),rep("macro2", macro_layer_between), rep("other", layer_between_processes), rep("none",(length(zlayers) - pdepth_layer))),
                         nlayer = c(rep(1,mdepth_layer), rep(1, macro_layer_between), rev(seq(1,length.out=layer_between_processes)), rep(10001, (length(zlayers) - pdepth_layer))),
                         dtheta = c(rep(dtheta_macro,mdepth_layer), rep(dtheta_macro2,macro_layer_between), rep(dtheta_other, layer_between_processes), rep(0,(length(zlayers) - pdepth_layer)))
              )
}
## by-pass flow
if(inf_dyn == 3){
layer_params = data.frame(Depth = zlayers,
                         infProcess = c(rep("macro", mdepth_layer),rep("macro2", macro_layer_between), rep("none", layer_between_processes), rep("other",(length(zlayers) - pdepth_layer))),
                         nlayer = c(rep(1,mdepth_layer), rep(1, macro_layer_between), rep(10001, layer_between_processes), seq(1,length.out=(length(zlayers) - pdepth_layer))),
                         dtheta = c(rep(dtheta_macro,mdepth_layer), rep(dtheta_macro2,macro_layer_between), rep(0, layer_between_processes), rep(dtheta_other,(length(zlayers) - pdepth_layer)))
              )
}
## DISTINGUISH order of nlayers for process "other", DEPTH (!?) for perched water table and by-pass flow !!

# old version
# layer_params = data.frame(Depth = zlayers,
#                          infProcess = c(rep("macro", mdepth_layer),rep("macro2", macro_layer_between), rep("other",(length(zlayers) - mdepth2_layer))),
#                          nlayer = c(rep(1,mdepth_layer), rep(1, macro_layer_between), seq(1,length.out=(length(zlayers) - mdepth2_layer))),
#                          # nlayer = c(rep(1,mdepth_layer), seq(1,length.out=(length(zlayers) - mdepth_layer))),
#                          dtheta = c(rep(dtheta_macro,mdepth_layer),rep(dtheta_macro2,macro_layer_between),rep(dtheta_other,(length(zlayers) - mdepth2_layer)))
#               )

## tag cells with layer parameters
tsx = inner_join(tsx, layer_params) # %>%
      # inner_join(Irrigation_grid)

## calculate number of "same" cells per layers in a column
cellnums = group_by(tsx, column, nlayer) %>%
# cellnums = group_by(tsx, column, nlayer, infProcess) %>%
           dplyr::summarize(cellsLayerColumn = sum(cnt, na.rm=T)) %>%
           ungroup()

## join with intensity distribution coefficents
tsx = inner_join(tsx, Intensity_distribution) %>%
      mutate(distrWater = water_vol_min * intensity / num_cell / vol_cell) %>%
      mutate(layerfill = 1) %>%
      inner_join(cellnums) %>%
      dplyr::select(-intensity) %>%
      mutate(sat = F) %>%
      mutate(aboveSat = F) %>%
      mutate(value_lat = NA) %>%
      mutate(lat_water = distrWater * latflow_fac) %>%
      mutate(cell_id = paste0(nlayer,"_",column,"_",infProcess))

## construct database with cell neighbours
celllayer = dplyr::filter(tsx, Depth == 0) %>%
                 dplyr::select(x,y,column)

# find side cells
cellneig_sides =  as.data.frame(adjacent(rasterFromXYZ(as.data.frame(celllayer)), cell = seq(1,num_cell), directions=4, sorted=TRUE))
cellneig_sides$nums = sequence(rle(cellneig_sides$from)$lengths)
cellneig_sides$columns = paste0("lat_s",cellneig_sides$nums)
# find corner cells
cellneig_corners =  as.data.frame(adjacent(rasterFromXYZ(as.data.frame(celllayer)), cell = seq(1,num_cell), directions="bishop", sorted=TRUE))
cellneig_corners$nums = sequence(rle(cellneig_corners$from)$lengths)
cellneig_corners$columns = paste0("lat_c",cellneig_corners$nums)

# cellneig_sides = spread(cellneig_sides, columns, to, fill=NA)
cellneig_sides = dcast(cellneig_sides, from ~ columns ,fill = NA, value.var = "to")
colnames(cellneig_sides)[1] = "column"
# cellneig_corners = spread(cellneig_corners, columns, to, fill=NA)
cellneig_corners = dcast(cellneig_corners, from ~ columns ,fill = NA, value.var = "to")
colnames(cellneig_corners)[1] = "column"

## join adjacent cells with cell-grid
tsx = left_join(tsx, cellneig_sides) %>%
      dplyr::mutate(lat_s1 = ifelse(!is.na(lat_s1),paste0(nlayer,"_",lat_s1), lat_s1))%>%
      dplyr::mutate(lat_s2 = ifelse(!is.na(lat_s2),paste0(nlayer,"_",lat_s2), lat_s2))%>%
      dplyr::mutate(lat_s3 = ifelse(!is.na(lat_s3),paste0(nlayer,"_",lat_s3), lat_s3))%>%
      dplyr::mutate(lat_s4 = ifelse(!is.na(lat_s4),paste0(nlayer,"_",lat_s4), lat_s4))%>%
      left_join(cellneig_corners) %>%
      dplyr::mutate(lat_c1 = ifelse(!is.na(lat_c1),paste0(nlayer,"_",lat_c1), lat_c1))%>%
      dplyr::mutate(lat_c2 = ifelse(!is.na(lat_c2),paste0(nlayer,"_",lat_c2), lat_c2))%>%
      dplyr::mutate(lat_c3 = ifelse(!is.na(lat_c3),paste0(nlayer,"_",lat_c3), lat_c3))%>%
      dplyr::mutate(lat_c4 = ifelse(!is.na(lat_c4),paste0(nlayer,"_",lat_c4), lat_c4))

## create mass balance error file
mb_error = data.frame(datetime = seq(1:precip_time), error = NA, corrected = NA)

####################
## start with for loop and TS1
for(i in 1:precip_time){ 
# for(i in 1:95){ 
  # i=2
## pass previous values to new time step data.frame
tsx$prevalue = tsx$value
# tsx$value = 0
tsx$datetime = i

## fill cells with water of time step x
####################
## vertical infiltration
## filling is done seperately for macro pore cells and other ("other") cells
tsx = dplyr::mutate(tsx, value_macro = ifelse(infProcess == "macro" & nlayer == layerfill, prevalue + (distrWater / cellsLayerColumn), prevalue)) %>%
      dplyr::mutate(value_macro = ifelse(infProcess == "macro2" & nlayer == layerfill, prevalue + (distrWater / cellsLayerColumn), value_macro)) %>%
      dplyr::mutate(value_macro = ifelse(infProcess == "macro" | infProcess == "macro2", value_macro, 0)) %>%
      ## normal filing
      dplyr::mutate(value_other = ifelse(infProcess == "other" & nlayer == layerfill, prevalue + (distrWater / cellsLayerColumn) * vertflow_fac, prevalue)) %>%
      ## adjust water amount for "vertically first filled cell", where all column water goes in !!
      dplyr::mutate(value_other = ifelse(infProcess == "other" & nlayer == layerfill & aboveSat == F, value_other + (distrWater / cellsLayerColumn) * (1 - vertflow_fac), value_other)) %>%
      dplyr::mutate(value_other = ifelse(infProcess == "other", value_other, 0)) %>%
      dplyr::mutate(value = value_macro + value_other)
####################

####################
## lateral water flow
## distribution of water into adjecent neighbouring cells
## all 8 cells in its surounding are considered, distinguishing
## 1) side cells, bordering with a complete cell_width
## 2) corner cells, connecting with the "origin cell" only at its corners
## procedure is:
## 1) find saturated cells
## 2) exclude saturated neighbour-cells
## 3) distribute later water (part of water_vol_min, depending on the factor lat_facflow chosen) into cells
## 4) recombine later water with already existing water in each cell
# sides
lateral_flow_sides = dplyr::filter(tsx, sat == T) %>%
# lateral_flow_sides = dplyr::filter(tt, sat == T) %>%
                    dplyr::select(cell_id, nlayer, lat_water, lat_s1, lat_s2, lat_s3, lat_s4) %>%
                    melt(id=c("cell_id", "lat_water", "nlayer")) #%>%
cnt_NOTna_cells = group_by(lateral_flow_sides, cell_id) %>%
                  dplyr::summarize(num_latcells = length(na.omit(value)))
lateral_flow_sides = left_join(lateral_flow_sides, cnt_NOTna_cells) %>%
                    dplyr::mutate(value_lat = lat_water * sfactor / num_latcells) # %>%
                    # dplyr::mutate(lat_cell = paste0(nlayer,"_",value)) # %>%
                    # dplyr::mutate(org_cell = cell_id)
# corners
lateral_flow_corners = dplyr::filter(tsx, sat == T) %>%
                    dplyr::select(cell_id, nlayer, lat_water, lat_c1, lat_c2, lat_c3, lat_c4) %>%
                    melt(id=c("cell_id", "lat_water", "nlayer")) #%>%
cnt_NOTna_cells = group_by(lateral_flow_corners, cell_id) %>%
                  dplyr::summarize(num_latcells = length(na.omit(value)))
lateral_flow_corners = left_join(lateral_flow_corners, cnt_NOTna_cells) %>%
                    dplyr::mutate(value_lat = lat_water * cfactor / num_latcells) # %>%
                    # dplyr::mutate(lat_cell = paste0(nlayer,"_",value)) # %>%
                    # dplyr::mutate(org_cell = cell_id)
# combine both
lateral_flow_sc = rbind(lateral_flow_sides, lateral_flow_corners)
                    # es muss auch noch geguckt werden, ob alle 4 cells je befüllt werden,
                    # sonst muss noch skaliert werden !!!!
## nach zusammenführen von corner und sides, noch nur 1 value_lat pro cellid schaffen
# lateral_flows = group_by(lateral_flow_sc, lat_cell) %>%
lateral_flows = group_by(lateral_flow_sc, value) %>%
                dplyr::summarize(value_lat = sum(value_lat, na.rm=T)) %>%
                dplyr::mutate(lat_cell = value) %>%
                dplyr::select(lat_cell, value_lat)

tsx = left_join(dplyr::select(tsx, -value_lat), as.data.frame(lateral_flows), by=c("cell_id" = "lat_cell")) %>%
      dplyr::mutate(value = ifelse(is.na(value_lat) | !is.finite(value_lat), value, value + value_lat))
####################

####################
## check mass (water) balance, before saving !!
mb_ts = sum(tsx$value * vol_cell, na.rm = T) - sum(tsx$prevalue * vol_cell, na.rm = T)
mb_error$error[i] = abs(1 - mb_ts / water_vol_min)
# mb_error$corrected[i] = abs(1 - mb_ts / water_vol_min)
# print(error)
## if mass balance error is too big
## adjust distributed water mass in this timestep
## this is done via scaling all distributed water (only of this timestep)
## to match with the available water (water_vol_min)
## this could lead to influenced and "wrong" distribution in cells !?
# while(abs(1 - mb_ts / water_vol_min) < mb_permitted_error){
if(abs(1 - mb_ts / water_vol_min) > mb_permitted_error){
    error_scaling = 1 / (mb_ts / water_vol_min)
    tsx = dplyr::mutate(tsx, value_ts_scaled = (value - prevalue) * error_scaling) %>%
          dplyr::mutate(value = prevalue + value_ts_scaled) %>%
          dplyr::select(-value_ts_scaled)
    ## calculate new mass balance, after re-scaling
    mb_ts_cor = sum(tsx$value * vol_cell, na.rm = T) - sum(tsx$prevalue * vol_cell, na.rm = T)
    mb_error$corrected[i] = abs(1 - mb_ts_cor / water_vol_min)
}

## combine with data from previous time steps
# IN FINAL VERSION: save only selected columns !!! -> memory issue after concating and reading back into R !!
# if(i == 1){
# write.table(tsx, file=paste0(dir_output, "raw/tsx_TS", formatC(i, width=3, flag="0"), ".txt"), sep="\t", row.names=F, col.names=T, append=F)
# }else{
# write.table(tsx, file=paste0(dir_output, "raw/tsx_TS", formatC(i, width=3, flag="0"), ".txt"), sep="\t", row.names=F, col.names=F, append=F)
# }
if(i == 1){
write.table(dplyr::select(tsx, x, y, z, Depth, value, datetime), 
            file=paste0(output_dir, "model_output/raw/tsx_TS", formatC(i, width=3, flag="0"), ".txt"), sep="\t", row.names=F, col.names=T, append=F)
}else{
write.table(dplyr::select(tsx, x, y, z, Depth, value, datetime),
            file=paste0(output_dir, "model_output/raw/tsx_TS", formatC(i, width=3, flag="0"), ".txt"), sep="\t", row.names=F, col.names=F, append=F)
}

## transmitting information to adjacent cells about saturation state
## also setting own state sat = T
# cells_sat = dplyr::mutate(tsx, saturated = ifelse(value >= dtheta, T, F)) %>%
            # dplyr::select(tsx, cell_id, nlayer, saturated, lat_s1, lat_s2, lat_s3, lat_s4, lat_c1, lat_c2, lat_c3, lat_c4)
cells_sat = dplyr::filter(tsx, value >= dtheta) %>%
            dplyr::mutate(belowSatnLayer = paste(cell_id,"_", (nlayer + 1))) %>%
            dplyr::select(cell_id, belowSatnLayer, nlayer, column)

# set saturation state of cell
tsx$sat[which(tsx$cell_id %in% cells_sat$cell_id)] = T
# give cell the information if the direct above cell is already saturated
# except for the macro pores
tsx$belowSatnLayer = paste(tsx$cell_id,"_", (tsx$nlayer))
tsx$aboveSat[which(tsx$belowSatnLayer %in% cells_sat$belowSatnLayer & tsx$infProcess != "macro" & tsx$infProcess != "macro2")] = T
# set information NA, if neighbour cell is already saturated
tsx$lat_s1[which(tsx$lat_s1 %in% cells_sat$cell_id)] = NA
tsx$lat_s2[which(tsx$lat_s2 %in% cells_sat$cell_id)] = NA
tsx$lat_s3[which(tsx$lat_s3 %in% cells_sat$cell_id)] = NA
tsx$lat_s4[which(tsx$lat_s4 %in% cells_sat$cell_id)] = NA
tsx$lat_c1[which(tsx$lat_c1 %in% cells_sat$cell_id)] = NA
tsx$lat_c2[which(tsx$lat_c2 %in% cells_sat$cell_id)] = NA
tsx$lat_c3[which(tsx$lat_c3 %in% cells_sat$cell_id)] = NA
tsx$lat_c4[which(tsx$lat_c4 %in% cells_sat$cell_id)] = NA

## determine which vertical cell of a column gets filled in the next timestep
## depends on saturation state of cell
layerfilling = dplyr::mutate(tsx, unsaturated = ifelse(sat, 10000, nlayer)) %>%
               group_by(column, infProcess) %>%
               dplyr::summarize(layerfill = min(unsaturated, na.rm=T))
### update which layers to fill
tsx = dplyr::select(tsx, - layerfill) %>%
      inner_join(layerfilling)

## determine the number of cells in each column, which are to be filled in the next timestep
## this value is used to divide the avaivable water per timestep
## and thus directly influences errors in mass balance
cellnums_dynamic = dplyr::mutate(layerfilling, ncells = ifelse(infProcess == "macro" & layerfill < 2, mdepth_layer, 0)) %>%
                   dplyr::mutate(ncells = ifelse(infProcess == "macro2" & layerfill < 2, (mdepth2_layer - mdepth_layer), ncells)) %>%
                   dplyr::mutate(ncells = ifelse(infProcess == "other", 1, ncells)) %>%
                   dplyr::group_by(column) %>%
                   dplyr::summarize(cellsLayerColumn = sum(ncells, na.rm=T))
tsx = dplyr::select(tsx, - cellsLayerColumn) %>%
      inner_join(cellnums_dynamic)

} # end for loop

## save mass balance error log
save(mb_error, file=paste0(output_dir, "model_output/MassBalanceError_", n_param, ".rdata"))

# # ####################
# # ## debugging
# # i = 360
# # load(file=paste0(dir_output, "/scenarios/combiExt/tsx_TS", i, ".rdata"))
# # sum(tsx$value) * vol_cell
# # watbal = data.frame(Timestep = 1:360, val=0)
# # for(i in 1:360){
# # for(i in 1:precip_time){
# # load(file=paste0(dir_output, "/scenarios/combiExt/tsx_TS", i, ".rdata"))
# # # water irrigated this timestep
# # val = sum(tsx$value) / i * vol_cell
# # watbal$Timestep[i] = i
# # watbal$val[i] = val
# # }
# # watbal$cumval = cumsum(watbal$val)
# # # plot vertical water distribution
# # cutrel = 8
# # cutposition = min(tsx$y) + cutrel
# # data_plot = dplyr::filter(tsx, y == cutposition)
# # maxy = max(data_plot$zgrid, na.rm=T)
# # ggplot(data_plot, aes(x = x, y = zgrid)) + geom_raster(aes(fill = value)) + ylim(maxy,0) + 
# # theme(legend.position="right") +
# # ylab("Depth [m]") + xlab("x [m]") + 
# # scale_fill_gradientn(colours = rev(viridis(7)))
# # ####################
# #               
# # print("finished modeling")
# print("starting data loading and stichting")
# # ####################
# # ## load and stich all datasets
# # ## define empty result data.frame
# # Irrigation_macropiping = data.frame()
# # 
# # for(i in 1:precip_time){
# # load(file=paste0(dir_output, "tsx_paramset_", n_param, "_TS", i, ".rdata"))
# # Irrigation_macropiping = rbind(Irrigation_macropiping, tsx)
# # }
# # print("finished data loading and stichting")
# # ## remove indexing, not needed columns (to save space)
# # # Irrigation_macropiping = dplyr::select(Irrigation_macropiping, - cnt, - layerfill)
# 
####################
## stich ascii files using bash
print("stiching output files")
systemcall_stich = paste0("cat ", output_dir, "model_output/raw/*.txt > ", output_dir, "model_output/raw/rawdata")
system(systemcall_stich)
####################
## load stitched file
print("read stiched output files")
# regular read.table
# data.tables fread
Infiltration_model_results = fread(file=paste0(output_dir, "model_output/raw/rawdata"), header = T, sep="\t", dec=".", 
                               ## read only columns needed
                               select = c("x","y","z","Depth","value","datetime"))
## convert data.table to data.frame
setDF(Infiltration_model_results)

####################
## check if all water was actually distributed
# water_bal = group_by(Infiltration_macropiping, Timestep) %>%
#                 dplyr::summarize(water_distributed = sum(value / Timestep, na.rm = T) * vol_cell)
# save(water_bal, file=paste0(dir_output,"water_bal_", n_param, ".rdata"))

####################
## reasonability checks ("sub"-unit testing on scenario basis)
# check total water of system
print("Water balance of every time step:")
# # check water of each timestep
plots=F
print(checkWaterVolumes(
        sm_mod_data = Infiltration_model_results,
        whichCheck = "waterpertimestep",
        inf_time = precip_time,
        cell_volume = vol_cell,
        infWater_timestep = water_vol_min)
)

# plot maximum soil moisture of each depth over experiment time
png(file=paste0(dir_plots, "Distribution_soilmoisture_max_", n_param, ".png"), width=1000, height=1000, res=250)
plot(plot_tempWaterDistribution(
        sm_mod_data = Infiltration_model_results,
        plotvar = "max")
)
dev.off()

# plot one transect for many (chosen) time steps
# define time steps to plot
plot_ts = seq(min(unique(Infiltration_model_results$datetime)), max(unique(Infiltration_model_results$datetime)), by = plot_interval)
for(ts_i in plot_ts){
png(file=paste0(dir_plots, "Transect_soilmoisture_TS_", plot_interval, "_", n_param, ".png"), width=1000, height=1000, res=250)
plot(plot_transect_2d(
        sm_mod_data = Infiltration_model_results,
        tstep = ts_i,
        y_pos = plot_transect_loc)
)
dev.off()
}

print(paste0("Saving Infiltration data for n_param: ", n_param))
## save data
save(Infiltration_model_results, file=paste0(output_dir, "model_output/Infiltration_model_output_", n_param, ".rData"))

# 
##debug
#load data
# load(file=paste0(dir_output, "Infiltration_combiExt_macropiping_", n_param, ".rdata"))
# print(str(Infiltration_model_results))
# print(str(gcomp_irrigation_domain))
# tt =dplyr::select(Infiltration_model_results, x,y,z,value,Timestep)
# print(str(tt))
# print(range(tt$Timestep))
# ttt = left_join(gcomp_irrigation_domain, tt)
# print(str(ttt))
# print(range(ttt$Timestep))
# print(range(ttt$value))
# print("other test irrigation dataset:")
# load(file="/home/hydro/mreich/Infiltration/output/combi/Infiltration_simple_macropiping.rdata")
# print("-zgrid")
# tt =dplyr::select(Infiltration_model_results, -zgrid)
# print(str(tt))
# ttt = left_join(gcomp_irrigation_domain, tt)
# print(str(ttt))
# print("select positiv")
# tt =dplyr::select(Infiltration_model_results, x,y,z,value,Timestep)
# print(str(tt))
# ttt = left_join(gcomp_irrigation_domain, tt)
# print(str(ttt))

########################################
## gravity response in nm/s²
########################################
## !!
## due to saving into ascii, reading and stiching (within the processing procedure)
## joining columns (x,y,z) have to be rounded! in order to match with gravity component grid
## !!
round_x = decimalplaces(grid_discretization$x)
round_y = decimalplaces(grid_discretization$y)
round_z = decimalplaces(grid_discretization$z)
Infiltration_model_results$x = round(Infiltration_model_results$x, round_x)
Infiltration_model_results$y = round(Infiltration_model_results$y, round_y)
Infiltration_model_results$z = round(Infiltration_model_results$z, round_z)
gcomp_grid$x = round(gcomp_grid$x, round_x)
gcomp_grid$y = round(gcomp_grid$y, round_y)
gcomp_grid$z = round(gcomp_grid$z, round_z)

# calculating gravity reponse
print("Calculating gravity response..")
infiltration_gmod = calculate_gravity_response(gcomp_grid, Infiltration_model_results)
#save gsignal
save(infiltration_gmod, file=paste0(output_dir, "model_output/GravityResponse_Infiltration_model_", n_param, ".rData"))
# print(str(infiltration_gmod))
# change column name for plotting
colnames(infiltration_gmod)[2] = "gmod"

## combine datasets and plot
gmod = rbind(
         cbind(gravity_timesteps, Scenario="Gravity signal observed"),
	     cbind(infiltration_gmod, Scenario="Gravity signal modeled")
         )

png(file=paste0(dir_plots, "GravityResponses_", n_param, ".png"), width=1500, height=1000, res=250)
plot(ggplot(gmod, aes(x=datetime, y=gmod, colour=Scenario)) + geom_line() + 
	ylab("Gravity [nm/s²]") + xlab("Time since sprinkling start"))
dev.off()

####################
## calculate fit to gobs

## regular KGE
# kge_value = KGE(infiltration_gmod$gmod, igrav_timesteps$gmod)
## changing scaling factor of component BIAS
kge_value = KGE(infiltration_gmod$gmod, gravity_timesteps$gmod, s=c(2.5/6,2.5/6,1/6))
kge_fit = 1 - kge_value

####################
## clean up memory
rm(Infiltration_model_results, infiltration_gmod)
gc()
## move to next n_param value for plot indexing
print(paste0("Finished parameterset: ",n_param))
n_param <<- n_param + 1

## returning quality criteria:
## KGE
return(kge_fit) 

## else statement, runs if otherdepth > mdepth
}

print("Finished model run.")
####################
} # end of function

