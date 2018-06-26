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

run_model_inversion_singleProcPP = function(
            dtheta_min,
            dtheta_max,
            # dtheta_macro2_min,
            # dtheta_macro2_max,
            # dtheta_other_min,
            # dtheta_other_max,
            # mdepth_min,
            # mdepth_max,
            # mdepth2_min,
            # mdepth2_max,
            # latflow_fac_min,
            # latflow_fac_max,
            # inf_dynamics_min,
            # inf_dynamics_max,
            pdepth_min,
            pdepth_max,
            dtheta_start,
            # dtheta_macro2_start,
            # mdepth_start,
            # mdepth2_start,
            # dtheta_other_start,
            # latflow_fac_start,
            # inf_dynamics_start,
            pdepth_start,
            output_dir, 
            input_dir,
            inner_inum,
            del_prev,
            ...
){
  ## library: necessary for parallel backend to work
  library(devtools)
  # setwd("/home/hydro/mreich/R/")
  load_all("/home/hydro/mreich/R/UmbrellaEffect")
  load_all("/home/hydro/mreich/R/gravityInf")
  ## create necessary output directories
  if(!file.exists(paste0(output_dir, "model_output"))){
      dir.create(file.path(output_dir, "model_output"))
      dir.create(file.path(output_dir, "model_output", "plots"))
      dir.create(file.path(output_dir, "model_output", "raw"))
  }
  ## delete results of previous inversion runs (and its log-files)
  # if not, ppso would try to extent these runs
  if(del_prev & file.exists(paste0(output_dir, "ppso.pro"))){
    file.remove(paste0(output_dir, "ppso.log"))
    file.remove(paste0(output_dir, "ppso.pro"))
  }

  # load config file
  load(file=paste0(output_dir, "configfile.rdata"))
  scenario = configfile$use_scenario
  # macropores = configfile$use_macro
  # macro2 = configfile$two_macro

  # prepare model input data and parameter boundaries
      # combine input parameters
      param_bounds = data.frame(minimum = c(dtheta_min, pdepth_min),
                            maximum = c(dtheta_max, pdepth_max))
      # combine start parameter values
      param_startvalue = c(dtheta_start, pdepth_start)
      # column names for model output
      model_info_colnames = c("dtheta", "pdepth")
      # set name of model to use
      model = "inf_model_3d_singleProcPP"
  
  # set counting parameter
  # this is used for naming figures and files for individiual model runs within the optimization routine
  n_param <<- inner_inum
  
  # store actual working directory
  wd_now = getwd()
  # set working directory
  setwd(output_dir)

  # get number of desired model runs
  n_iterations = configfile$model_runs
  
  ## run optimization
  opt_result = optim_ppso_robust(
      # objective_function = inf_model_3d_2macro, #set model to use, according to input parameters supplied
      objective_function = get(model), #set model to use, according to input parameters supplied
      number_of_parameters = length(param_bounds[,1]),
      number_of_particles =  10,
      max_number_function_calls= n_iterations,
      nslaves = 500,
      # r=0.2,
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
      verbose=TRUE,
      tryCall=FALSE)

  # set previous working directory
  setwd(wd_now)

  ## combine and save model results and model input
  # find lower infiltration process investigated and convert from number to word for better result recognition
  # switch(inf_dynamics,
  #        "1" = {inf_process = "Wfa"},
  #        "2" = {inf_process = "Perched water table"},
  #        "3" = {inf_process = "By-pass flow"}
  #        )
  # scenario information
  model_info = data.frame(dir_input = configfile$dir_input,
                       dir_output = configfile$dir_output,
                       duration = configfile$precip_time,
                       total_water_volume = configfile$water_vol_min,
                       water_distribution_file = configfile$IntensityDistribution_file,
                       gravity_component_grid_file = configfile$gcomp_file,
                       gravity_observations_file = configfile$gravityObservations_file,
                       permitted_massbalance_error = configfile$mb_permitted_error,
                       use_scenario = configfile$use_scenario,
                       # macropore_layer = configfile$use_macro,
                       # macropore_layer2 = configfile$two_macro,
                       n_iterations = configfile$model_runs,
                       plot_interval = configfile$plot_interval,
                       plot_transect_2d = configfile$plot_transect_loc,
                       opt_model_kge = opt_result$value,
                       opt_function_calls = opt_result$function_calls,
                       opt_break_flag = opt_result$break_flag,
                       stringsAsFactors = F
                       )

  # add estimated optimal parameter values
  info_columns_num = length(model_info)
  for(param_num in 1:length(param_bounds[,1])){
    model_info[,(param_num + info_columns_num)] = opt_result$par[param_num]
    colnames(model_info)[(param_num + info_columns_num)] = model_info_colnames[param_num]
  }

  ## close cluster
  print("closing R slaves")
  mpi.close.Rslaves()
  ##
  # return model statistics and parameters
  return(model_info)
}

