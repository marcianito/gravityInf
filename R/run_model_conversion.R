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
            dtheta_macro_start = dtheta_macro_start,
            dtheta_macro2_start = dtheta_macro2_start,
            mdepth_start = mdepth_start,
            mdepth2_start = mdepth2_start,
            dtheta_other_start = dtheta_other_start,
            latflow_fac_start = latflow_fac_start,
            output_dir, 
            ...
){
  ## create necessary output directories
  if(!file.exists(paste0(output_dir, "model_output"))){
      dir.create(file.path(output_dir, "model_output"))
      dir.create(file.path(output_dir, "model_output", "plots"))
      # dir.create(file.path(output_dir, "model_output", "raw"))
  }

  # load config file
  load(file=paste0(output_dir, "configfile.rdata"))
  macropores = configfile$use_macro
  macro2 = configfile$two_macro

  # prepare model input data and parameter boundaries
  if(macropores){
    if(macro2){
      # combine input parameters
      param_startvalue = c(dtheta_macro_start, dtheta_macro2_start, dtheta_other_start, mdepth_start, mdepth2_start, latflow_fac_start)
      # run model 
      model_result = inf_model_3d_2macro_conv(
                        param_vec = param_startvalue
                     )
    }else{
      # combine input parameters
      param_startvalue = c(dtheta_macro_start, dtheta_other_start, mdepth_start, latflow_fac_start)
      # run model 
      model_result = inf_model_3d_macro_conv(
                        param_vec = param_startvalue
                     )
    }
  }else{
    # combine input parameters
    param_startvalue = c(dtheta_other_start, mdepth_start, latflow_fac_start)
      # run model 
      model_result = inf_model_3d_single_conv(
                        param_vec = param_startvalue
                     )
  }
    
  ## combine and save model results and model input
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
                       infiltration_dynamics = configfile$inf_dynamics,
                       n_iterations = configfile$model_runs,
                       plot_interval = configfile$plot_interval,
                       plot_transect_2d = configfile$plot_transect_loc,
                       stringsAsFactors = F
                       )

  # add estimated optimal parameter values
  info_columns_num = length(model_info)
  for(param_num in 1:length(param_startvalue)){
    model_info[,(param_num + info_columns_num)] = param_startvalue[param_num]
  }

  # return model statistics and parameters
  return(model_info)
}

