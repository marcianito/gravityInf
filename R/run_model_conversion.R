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
            output_dir, 
            ...
){
  ## create necessary output directories
  if(!file.exists(paste0(output_dir, "model_output"))){
      dir.create(file.path(output_dir, "model_output"))
      dir.create(file.path(output_dir, "model_output", "plots"))
      dir.create(file.path(output_dir, "model_output", "raw"))
  }

  # load config file for model input
  load(file=paste0(output_dir, "configfile.rdata"))
  dir_input = configfile$dir_input
  dir_output = configfile$dir_output
  precip_time = configfile$precip_time
  IntensityDistribution = configfile$IntensityDistribution
  water_vol_min = configfile$water_vol_min
  gcompfile = configfile$gcompfile
  gravityObs = configfile$gravityObs
  mb_permitted_error = configfile$mb_permitted_error

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
    
  ## combine and save model results and model input
  conv_result = list(
  # infiltration parameters

  )
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
  stats = cbind(model_info, conv_result)
  
  # return model statistics and parameters
  return(stats)
}

