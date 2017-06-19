#' @title Run infiltration model in inversion mode
#'
#' @description The infiltration model is run wrapped within a
#'  particle-swarm optimization algorithm, based on R-package 'ppso'.
#'
#' @param test test
#' 
#' @return A list (?) of the optimized model parameters is returned.
#' 
#' @details It can be decided if one, two or none marco pore layer is desired.
#'  For modifications of the optimization algorithm, variables have to be supplied to the function
#'  as described in ?optim_dds.
#'
#' @references Marvin Reich (2017), mreich@@posteo.de
#' @examples missing

run_model_inversion = function(
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
            n_iterations = 100,
            output_dir, 
            ...
){
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
    
    # set working directory
    # setwd(dir_output)
    
    ## run optimization
    opt_result = optim_dds(
        # objective_function = inf_model_3d_2macro, #set model to use, according to input parameters supplied
        objective_function = get(model), #set model to use, according to input parameters supplied
        number_of_parameters = length(param_bounds[,1]),
        number_of_particles =  1,
        max_number_function_calls= n_iterations,
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

  ## combine and save model results and model input
  # load config file
  load(file=paste0(output_dir, "configfile.rdata"))
  dir_input = configfile$dir_input
  dir_output = configfile$dir_output
  precip_time = configfile$precip_time
  IntensityDistribution = configfile$IntensityDistribution
  water_vol_min = configfile$water_vol_min
  gcompfile = configfile$gcompfile
  gravityObs = configfile$gravityObs
  mb_permitted_error = configfile$mb_permitted_error

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
  
  # return model statistics and parameters
  return(stats)
}

