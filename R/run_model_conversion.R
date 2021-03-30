#' @title Run infiltration model in conversion mode
#'
#' @description The infiltration model is run one in conversion mode.
#' 
#' @param test test
#' 
#' @return A list (?) of the model parameters and model input is returned.
#' 
#' @details missing
#'
#' @references Marvin Reich (2017), mreich@@posteo.de
#' @examples missing

run_model_conversion = function(
            dtheta_macro,
            dtheta_macro2,
            mdepth,
            mdepth2,
            dtheta_other,
            latflow_fac,
            inf_dynamics,
            pdepth,
            output_dir, 
            ...
){
   # dtheta_macro = dtheta_macro_start
   # dtheta_macro2 = dtheta_macro2_start
   # mdepth = mdepth_start
   # mdepth2 = mdepth2_start
   # dtheta_other = dtheta_other_start
   # latflow_fac = latflow_fac_start
   # output_dir = dir_output

  ## create necessary output directories
  if(!file.exists(paste0(output_dir, "model_output"))){
      dir.create(file.path(output_dir, "model_output"))
      dir.create(file.path(output_dir, "model_output", "plots"))
      dir.create(file.path(output_dir, "model_output", "raw"))
  }

  # load config file
  load(file=paste0(output_dir, "configfile.rdata"))
  macropores = configfile$use_macro
  macro2 = configfile$two_macro

  # set counting parameter
  # this is used for naming figures and files for individiual model runs within the optimization routine
  n_param <<- 1
  # store actual working directory
  wd_now = getwd()
  # set working directory
  setwd(output_dir)

  # prepare model input data and parameter boundaries
  if(macropores){
    if(macro2){
      # combine input parameters
      parameter_input = c(dtheta_macro, dtheta_macro2, dtheta_other, mdepth, mdepth2, latflow_fac, inf_dynamics, pdepth)
      model_info_colnames = c("dtheta_macro", "dtheta_macro2", "dtheta_other", "mdepth", "mdepth2", "latflow_fac", "inf_dynamics", "pdepth")
      # run model 
      model_result = inf_model_3d_2macro(
                        param_vec = parameter_input
                     )
    }else{
      # combine input parameters
      parameter_input = c(dtheta_macro, dtheta_other, mdepth, latflow_fac, inf_dynamics, pdepth)
      model_info_colnames = c("dtheta_macro", "dtheta_other", "mdepth", "latflow_fac", "inf_dynamics", "pdepth")
      # run model 
      model_result = inf_model_3d_macro(
                        param_vec = parameter_input
                     )
    }
  }else{
    # combine input parameters
    parameter_input = c(dtheta_other, mdepth, latflow_fac, inf_dynamics, pdepth)
    model_info_colnames = c("dtheta_other", "mdepth", "latflow_fac", "inf_dynamics", "pdepth")
      # run model 
      model_result = inf_model_3d_single(
                        param_vec = parameter_input
                     )
  }
  # output to user
  print(paste0("Model results and figures were saved to: ", output_dir, "model_output/."))

  # set previous working directory
  setwd(wd_now)

  ## combine and save model results and model input
  # find lower infiltration process investigated and convert from number to word for better result recognition
  switch(inf_dynamics,
         "1" = {inf_process = "Wfa"},
         "2" = {inf_process = "Perched water table"},
         "3" = {inf_process = "By-pass flow"}
         )
  # scenario information
  model_info = data.frame(dir_input = configfile$dir_input,
                       dir_output = configfile$dir_output,
                       duration = configfile$precip_time,
                       total_water_volume = configfile$water_vol_min,
                       water_distribution_file = configfile$IntensityDistribution_file,
                       gravity_component_grid_file = configfile$gcomp_file,
                       gravity_observations_file = configfile$gravityObservations_file,
                       permitted_massbalance_error = configfile$mb_permitted_error,
                       macropore_layer = configfile$use_macro,
                       macropore_layer2 = configfile$two_macro,
                       model_kge = model_result,
                       inf_process_tested = inf_process,
                       plot_interval = configfile$plot_interval,
                       plot_transect_2d = configfile$plot_transect_loc,
                       stringsAsFactors = F
                       )

  # add estimated optimal parameter values
  info_columns_num = length(model_info)
  for(param_num in 1:length(parameter_input)){
    model_info[,(param_num + info_columns_num)] = parameter_input[param_num]
    colnames(model_info)[(param_num + info_columns_num)] = model_info_colnames[param_num]
  }
  # set column names
  # colnames(model_info)[(info_columns_num + 1):(info_columns_num + length(model_info_colnames))] = model_info_colnames
  
  # return model statistics and parameters
  return(model_info)
}

