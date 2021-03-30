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

run_model_inversion_combinedProcPP_dthetafree = function(
            dtheta1_min,
            dtheta1_max,
            dtheta2_min,
            dtheta2_max,
            dtheta3_min,
            dtheta3_max,
            dtheta4_min,
            dtheta4_max,
            dtheta5_min,
            dtheta5_max,
            dtheta6_min,
            dtheta6_max,
            dtheta7_min,
            dtheta7_max,
            dtheta8_min,
            dtheta8_max,
            dtheta9_min,
            dtheta9_max,
            dtheta10_min,
            dtheta10_max,
            dtheta11_min,
            dtheta11_max,
            dtheta12_min,
            dtheta12_max,
            dtheta13_min,
            dtheta13_max,
            dtheta14_min,
            dtheta14_max,
            dtheta15_min,
            dtheta15_max,
            dtheta16_min,
            dtheta16_max,
            dtheta17_min,
            dtheta17_max,
            dtheta18_min,
            dtheta18_max,
            dtheta19_min,
            dtheta19_max,
            dtheta20_min,
            dtheta20_max,
            dtheta21_min,
            dtheta21_max,
            dtheta22_min,
            dtheta22_max,
            dtheta23_min,
            dtheta23_max,
            dtheta24_min,
            dtheta24_max,
            dtheta25_min,
            dtheta25_max,
            dtheta26_min,
            dtheta26_max,
            dtheta27_min,
            dtheta27_max,
            dtheta28_min,
            dtheta28_max,
            dtheta29_min,
            dtheta29_max,
            dtheta30_min,
            dtheta30_max,
            dtheta31_min,
            dtheta31_max,
            dtheta32_min,
            dtheta32_max,
            dtheta33_min,
            dtheta33_max,
            dtheta34_min,
            dtheta34_max,
            dtheta35_min,
            dtheta35_max,
            dtheta36_min,
            dtheta36_max,
            dtheta37_min,
            dtheta37_max,
            dtheta38_min,
            dtheta38_max,
            dtheta39_min,
            dtheta39_max,
            dtheta40_min,
            dtheta40_max,
            dtheta41_min,
            dtheta41_max,
            dtheta42_min,
            dtheta42_max,
            dtheta43_min,
            dtheta43_max,
            dtheta44_min,
            dtheta44_max,
            dtheta45_min,
            dtheta45_max,
            dtheta46_min,
            dtheta46_max,
            dtheta47_min,
            dtheta47_max,
            dtheta48_min,
            dtheta48_max,
            dtheta49_min,
            dtheta49_max,
            dtheta50_min,
            dtheta50_max,
            latflow_fac_min = NA,
            latflow_fac_max = NA,
            pdepth_min,
            pdepth_max,
            pdepth2_min,
            pdepth2_max,
            dtheta1_start,
            dtheta2_start,
            dtheta3_start,
            dtheta4_start,
            dtheta5_start,
            dtheta6_start,
            dtheta7_start,
            dtheta8_start,
            dtheta9_start,
            dtheta10_start,
            dtheta11_start,
            dtheta12_start,
            dtheta13_start,
            dtheta14_start,
            dtheta15_start,
            dtheta16_start,
            dtheta17_start,
            dtheta18_start,
            dtheta19_start,
            dtheta20_start,
            dtheta21_start,
            dtheta22_start,
            dtheta23_start,
            dtheta24_start,
            dtheta25_start,
            dtheta26_start,
            dtheta27_start,
            dtheta28_start,
            dtheta29_start,
            dtheta30_start,
            dtheta31_start,
            dtheta32_start,
            dtheta33_start,
            dtheta34_start,
            dtheta35_start,
            dtheta36_start,
            dtheta37_start,
            dtheta38_start,
            dtheta39_start,
            dtheta40_start,
            dtheta41_start,
            dtheta42_start,
            dtheta43_start,
            dtheta44_start,
            dtheta45_start,
            dtheta46_start,
            dtheta47_start,
            dtheta48_start,
            dtheta49_start,
            dtheta50_start,
            latflow_fac_start = NA,
            pdepth_start,
            pdepth2_start,
            output_dir, 
            input_dir,
            inner_inum,
            del_prev,
            ...
){
  ## library: necessary for parallel backend to work
  library(devtools)
  # setwd("/home/hydro/mreich/R/")
  load_all("/home/mreich/R/HyGra")
  # load_all("/home/hydro/mreich/R/UmbrellaEffect")
  load_all("/home/mreich/R/gravityInf")
  ## create necessary output directories
  if(!file.exists(paste0(output_dir, "model_output"))){
      dir.create(file.path(output_dir, "model_output"))
      dir.create(file.path(output_dir, "model_output", "plots"))
      dir.create(file.path(output_dir, "model_output", "raw"))
  }
  if(!file.exists(paste0(output_dir, "model_output/plots"))){
      dir.create(file.path(output_dir, "model_output", "plots"))
  }
  if(!file.exists(paste0(output_dir, "model_output/raw"))){
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
  scenario2 = configfile$use_scenario2

  # prepare model input data and parameter boundaries
  # no lateral flow module
  if(is.na(latflow_fac_start)){
      # combine input parameters
      print("this is still not implemented")
      # param_bounds = data.frame(minimum = c(dtheta_min, pdepth_min),
      #                       maximum = c(dtheta_max, pdepth_max))
      # # combine start parameter values
      # param_startvalue = c(dtheta_start, pdepth_start)
      # # column names for model output
      # model_info_colnames = c("dtheta", "pdepth")
      # # set name of model to use
      # model = "inf_model_3d_combinedProcPP"
  # DO use lateral flow module
  }else{
      switch(scenario2,
      bypass = {
        # combine input parameters
        param_bounds = data.frame(minimum = c(dtheta1_min, dtheta2_min, dtheta3_min, dtheta4_min, dtheta5_min, dtheta6_min, dtheta7_min, dtheta8_min, dtheta9_min, dtheta10_min, dtheta11_min, dtheta12_min, dtheta13_min, dtheta14_min, dtheta15_min, dtheta16_min, dtheta17_min, dtheta18_min, dtheta19_min, dtheta20_min, dtheta21_min, dtheta22_min, dtheta23_min, dtheta24_min, dtheta25_min, dtheta26_min, dtheta27_min, dtheta28_min, dtheta29_min, dtheta30_min, dtheta31_min, dtheta32_min, dtheta33_min, dtheta34_min, dtheta35_min, dtheta36_min, dtheta37_min, dtheta38_min, dtheta39_min, dtheta40_min, dtheta41_min, dtheta42_min, dtheta43_min, dtheta44_min, dtheta45_min, dtheta46_min, dtheta47_min, dtheta48_min, dtheta49_min, dtheta50_min, pdepth_min, pdepth2_min, latflow_fac_min),
                              maximum = c(dtheta1_max, dtheta2_max, dtheta3_max, dtheta4_max, dtheta5_max, dtheta6_max, dtheta7_max, dtheta8_max, dtheta9_max, dtheta10_max, dtheta11_max, dtheta12_max, dtheta13_max, dtheta14_max, dtheta15_max, dtheta16_max, dtheta17_max, dtheta18_max, dtheta19_max, dtheta20_max, dtheta21_max, dtheta22_max, dtheta23_max, dtheta24_max, dtheta25_max, dtheta26_max, dtheta27_max, dtheta28_max, dtheta29_max, dtheta30_max, dtheta31_max, dtheta32_max, dtheta33_max, dtheta34_max, dtheta35_max, dtheta36_max, dtheta37_max, dtheta38_max, dtheta39_max, dtheta40_max, dtheta41_max, dtheta42_max, dtheta43_max, dtheta44_max, dtheta45_max, dtheta46_max, dtheta47_max, dtheta48_max, dtheta49_max, dtheta50_max, pdepth_max, pdepth2_max, latflow_fac_max))
        # combine start parameter values
        param_startvalue = c(dtheta1_start, dtheta2_start, dtheta3_start, dtheta4_start, dtheta5_start, dtheta6_start, dtheta7_start, dtheta8_start, dtheta9_start, dtheta10_start, dtheta11_start, dtheta12_start, dtheta13_start, dtheta14_start, dtheta15_start, dtheta16_start, dtheta17_start, dtheta18_start, dtheta19_start, dtheta20_start, dtheta21_start, dtheta22_start, dtheta23_start, dtheta24_start, dtheta25_start, dtheta26_start, dtheta27_start, dtheta28_start, dtheta29_start, dtheta30_start, dtheta31_start, dtheta32_start, dtheta33_start, dtheta34_start, dtheta35_start, dtheta36_start, dtheta37_start, dtheta38_start, dtheta39_start, dtheta40_start, dtheta41_start, dtheta42_start, dtheta43_start, dtheta44_start, dtheta45_start, dtheta46_start, dtheta47_start, dtheta48_start, dtheta49_start, dtheta50_start, pdepth_start, pdepth2_start, latflow_fac_start)
        # column names for model output
        model_info_colnames = c("dtheta1", "dtheta2", "dtheta3", "dtheta4", "dtheta5", "dtheta6", "dtheta7", "dtheta8", "dtheta9", "dtheta10", "dtheta11", "dtheta12", "dtheta13", "dtheta14", "dtheta15", "dtheta16", "dtheta17", "dtheta18", "dtheta19", "dtheta20", "dtheta21", "dtheta22", "dtheta23", "dtheta24", "dtheta25", "dtheta26", "dtheta27", "dtheta28", "dtheta29", "dtheta30", "dtheta31", "dtheta32", "dtheta33", "dtheta34", "dtheta35", "dtheta36", "dtheta37", "dtheta38", "dtheta39", "dtheta40", "dtheta41", "dtheta42", "dtheta43", "dtheta44", "dtheta45", "dtheta46", "dtheta47", "dtheta48", "dtheta49", "dtheta50", "pdepth", "pdepth2", "latflow")
        # set name of model to use
        model = "inf_model_3d_combinedProcLatFlowPP_dthetafree"
             },
      perched = {
        # combine input parameters
        param_bounds = data.frame(minimum = c(dtheta1_min, dtheta2_min, dtheta3_min, dtheta4_min, dtheta5_min, dtheta6_min, dtheta7_min, dtheta8_min, dtheta9_min, dtheta10_min, dtheta11_min, dtheta12_min, dtheta13_min, dtheta14_min, dtheta15_min, dtheta16_min, dtheta17_min, dtheta18_min, dtheta19_min, dtheta20_min, dtheta21_min, dtheta22_min, dtheta23_min, dtheta24_min, dtheta25_min, dtheta26_min, dtheta27_min, dtheta28_min, dtheta29_min, dtheta30_min, dtheta31_min, dtheta32_min, dtheta33_min, dtheta34_min, dtheta35_min, dtheta36_min, dtheta37_min, dtheta38_min, dtheta39_min, dtheta40_min, dtheta41_min, dtheta42_min, dtheta43_min, dtheta44_min, dtheta45_min, dtheta46_min, dtheta47_min, dtheta48_min, dtheta49_min, dtheta50_min, pdepth_min, pdepth2_min, latflow_fac_min),
                              maximum = c(dtheta1_max, dtheta2_max, dtheta3_max, dtheta4_max, dtheta5_max, dtheta6_max, dtheta7_max, dtheta8_max, dtheta9_max, dtheta10_max, dtheta11_max, dtheta12_max, dtheta13_max, dtheta14_max, dtheta15_max, dtheta16_max, dtheta17_max, dtheta18_max, dtheta19_max, dtheta20_max, dtheta21_max, dtheta22_max, dtheta23_max, dtheta24_max, dtheta25_max, dtheta26_max, dtheta27_max, dtheta28_max, dtheta29_max, dtheta30_max, dtheta31_max, dtheta32_max, dtheta33_max, dtheta34_max, dtheta35_max, dtheta36_max, dtheta37_max, dtheta38_max, dtheta39_max, dtheta40_max, dtheta41_max, dtheta42_max, dtheta43_max, dtheta44_max, dtheta45_max, dtheta46_max, dtheta47_max, dtheta48_max, dtheta49_max, dtheta50_max, pdepth_max, pdepth2_max, latflow_fac_max))
        # combine start parameter values
        param_startvalue = c(dtheta1_start, dtheta2_start, dtheta3_start, dtheta4_start, dtheta5_start, dtheta6_start, dtheta7_start, dtheta8_start, dtheta9_start, dtheta10_start, dtheta11_start, dtheta12_start, dtheta13_start, dtheta14_start, dtheta15_start, dtheta16_start, dtheta17_start, dtheta18_start, dtheta19_start, dtheta20_start, dtheta21_start, dtheta22_start, dtheta23_start, dtheta24_start, dtheta25_start, dtheta26_start, dtheta27_start, dtheta28_start, dtheta29_start, dtheta30_start, dtheta31_start, dtheta32_start, dtheta33_start, dtheta34_start, dtheta35_start, dtheta36_start, dtheta37_start, dtheta38_start, dtheta39_start, dtheta40_start, dtheta41_start, dtheta42_start, dtheta43_start, dtheta44_start, dtheta45_start, dtheta46_start, dtheta47_start, dtheta48_start, dtheta49_start, dtheta50_start, pdepth_start, pdepth2_start, latflow_fac_start)
        # column names for model output
        model_info_colnames = c("dtheta1", "dtheta2", "dtheta3", "dtheta4", "dtheta5", "dtheta6", "dtheta7", "dtheta8", "dtheta9", "dtheta10", "dtheta11", "dtheta12", "dtheta13", "dtheta14", "dtheta15", "dtheta16", "dtheta17", "dtheta18", "dtheta19", "dtheta20", "dtheta21", "dtheta22", "dtheta23", "dtheta24", "dtheta25", "dtheta26", "dtheta27", "dtheta28", "dtheta29", "dtheta30", "dtheta31", "dtheta32", "dtheta33", "dtheta34", "dtheta35", "dtheta36", "dtheta37", "dtheta38", "dtheta39", "dtheta40", "dtheta41", "dtheta42", "dtheta43", "dtheta44", "dtheta45", "dtheta46", "dtheta47", "dtheta48", "dtheta49", "dtheta50", "pdepth", "pdepth2", "latflow")
        # set name of model to use
        model = "inf_model_3d_combinedProcLatFlowPP_dthetafree"
             },
      wfa = {
        # combine input parameters
        param_bounds = data.frame(minimum = c(dtheta1_min, dtheta2_min, dtheta3_min, dtheta4_min, dtheta5_min, dtheta6_min, dtheta7_min, dtheta8_min, dtheta9_min, dtheta10_min, dtheta11_min, dtheta12_min, dtheta13_min, dtheta14_min, dtheta15_min, dtheta16_min, dtheta17_min, dtheta18_min, dtheta19_min, dtheta20_min, dtheta21_min, dtheta22_min, dtheta23_min, dtheta24_min, dtheta25_min, dtheta26_min, dtheta27_min, dtheta28_min, dtheta29_min, dtheta30_min, dtheta31_min, dtheta32_min, dtheta33_min, dtheta34_min, dtheta35_min, dtheta36_min, dtheta37_min, dtheta38_min, dtheta39_min, dtheta40_min, dtheta41_min, dtheta42_min, dtheta43_min, dtheta44_min, dtheta45_min, dtheta46_min, dtheta47_min, dtheta48_min, dtheta49_min, dtheta50_min, pdepth_min, latflow_fac_min),
                              maximum = c(dtheta1_max, dtheta2_max, dtheta3_max, dtheta4_max, dtheta5_max, dtheta6_max, dtheta7_max, dtheta8_max, dtheta9_max, dtheta10_max, dtheta11_max, dtheta12_max, dtheta13_max, dtheta14_max, dtheta15_max, dtheta16_max, dtheta17_max, dtheta18_max, dtheta19_max, dtheta20_max, dtheta21_max, dtheta22_max, dtheta23_max, dtheta24_max, dtheta25_max, dtheta26_max, dtheta27_max, dtheta28_max, dtheta29_max, dtheta30_max, dtheta31_max, dtheta32_max, dtheta33_max, dtheta34_max, dtheta35_max, dtheta36_max, dtheta37_max, dtheta38_max, dtheta39_max, dtheta40_max, dtheta41_max, dtheta42_max, dtheta43_max, dtheta44_max, dtheta45_max, dtheta46_max, dtheta47_max, dtheta48_max, dtheta49_max, dtheta50_max, pdepth_max, latflow_fac_max))
        # combine start parameter values
        param_startvalue = c(dtheta1_start, dtheta2_start, dtheta3_start, dtheta4_start, dtheta5_start, dtheta6_start, dtheta7_start, dtheta8_start, dtheta9_start, dtheta10_start, dtheta11_start, dtheta12_start, dtheta13_start, dtheta14_start, dtheta15_start, dtheta16_start, dtheta17_start, dtheta18_start, dtheta19_start, dtheta20_start, dtheta21_start, dtheta22_start, dtheta23_start, dtheta24_start, dtheta25_start, dtheta26_start, dtheta27_start, dtheta28_start, dtheta29_start, dtheta30_start, dtheta31_start, dtheta32_start, dtheta33_start, dtheta34_start, dtheta35_start, dtheta36_start, dtheta37_start, dtheta38_start, dtheta39_start, dtheta40_start, dtheta41_start, dtheta42_start, dtheta43_start, dtheta44_start, dtheta45_start, dtheta46_start, dtheta47_start, dtheta48_start, dtheta49_start, dtheta50_start, pdepth_start, latflow_fac_start)
        # column names for model output
        model_info_colnames = c("dtheta1", "dtheta2", "dtheta3", "dtheta4", "dtheta5", "dtheta6", "dtheta7", "dtheta8", "dtheta9", "dtheta10", "dtheta11", "dtheta12", "dtheta13", "dtheta14", "dtheta15", "dtheta16", "dtheta17", "dtheta18", "dtheta19", "dtheta20", "dtheta21", "dtheta22", "dtheta23", "dtheta24", "dtheta25", "dtheta26", "dtheta27", "dtheta28", "dtheta29", "dtheta30", "dtheta31", "dtheta32", "dtheta33", "dtheta34", "dtheta35", "dtheta36", "dtheta37", "dtheta38", "dtheta39", "dtheta40", "dtheta41", "dtheta42", "dtheta43", "dtheta44", "dtheta45", "dtheta46", "dtheta47", "dtheta48", "dtheta49", "dtheta50", "pdepth", "pdepth2", "latflow")
        # set name of model to use
        model = "inf_model_3d_combinedProcLatFlow_wfaPP_dthetafree"
      }
      )
  }
  
  # set counting parameter
  # this is used for naming figures and files for individiual model runs within the optimization routine
  # n_param <<- inner_inum
  
  # store actual working directory
  wd_now = getwd()
  # set working directory
  setwd(output_dir)

  # get number of desired model runs
  n_iterations = configfile$model_runs

  # debug
  print("before ppso model call..")
  library(Rmpi)
  print("number of available slaves (from mpi system):")
  print(Rmpi::mpi.universe.size())
  print("param bounds: ")
  print(length(param_bounds[,1]))
  #
  
  ## run optimization
  opt_result = optim_ppso_robust(
      # objective_function = inf_model_3d_2macro, #set model to use, according to input parameters supplied
      objective_function = get(model), #set model to use, according to input parameters supplied
      number_of_parameters = length(param_bounds[,1]),
      number_of_particles =  10,
      max_number_function_calls= n_iterations,
      max_number_of_iterations= n_iterations,
      # nslaves = 10,
      # r=0.2,
      abstol = -Inf,
      reltol = -Inf,
      # max_wait_iterations=50,
      max_wait_iterations=n_iterations,
      wait_complete_iteration = FALSE,
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
      tryCall=TRUE)

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
                       use_scenario2 = configfile$use_scenario2,
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
  # delete raw directory
  # this has the reason for unique file names, 
  # thus adding new files for every run and
  # a maximum number of files per folder on cluster of 51 million
  # unlink(paste0(output_dir, "model_output/raw"), recursive = TRUE)

  # return model statistics and parameters
  return(model_info)
}

