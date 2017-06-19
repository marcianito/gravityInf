#' @title Plot gravity response of observed and modeled gravity signal
#'
#' @description Plot gravity time series obtained through measurements during the sprinkling
#' experiment and modeling of infiltration of sprinkled water.
#'
#' @param gravity_obs Character string with the filename of the gravity observation data to be used.
#' @param gravity_mod Character string with the filename of the modeled gravity response data.
#' @param input_dir Character string, containing the path to input file "gravity_obs".
#' @param output_dir Character string, containing the path to the directory where to save the plot.
#' 
#' @return Returns NULL as object. It creates a plot of all input time series.
#' This plot is shown on screen as well as saved to the specified output folder.
#' 
#' @details missing
#' @references Marvin Reich (2017), mreich@@posteo.de
#' @examples missing

plot_gravity_responses = function(
            gravity_obs,
            gravity_mod = "gravity_response_modeled.rData",
            input_dir = dir_input,
            output_dir = dir_output,
            dat_tsf = 7,
            ...
){
    # gravity_obs = gravityObservations_input_file
    # gravity_outside = gravity_response_outside_building
    # gravity_below = gravity_response_below_building
    # gravity_reduced = gravity_data_reduced
    # input_dir = dir_input

    # subtract starting value (of first timestep)
    # subtract starting value (of first timestep)
    # this is done to create relative mass change responses and not absolut
    # gravity reponse outside & below
    # gravity_outside$value = gravity_outside$value - gravity_outside$value[1]
    # gravity_below$value = gravity_below$value - gravity_below$value[1]
    # subtract mean value 
    
    # gravity_outside$value = gravity_outside$value - mean(gravity_outside$value, na.rm = T)
    # gravity_below$value = gravity_below$value - mean(gravity_below$value, na.rm = T)

    ## read in routine for gravity obs data
    # same as in function "reduce_gravity"
    # this time resulting in columns $datetime, $value

    # load data observed gravity data
    gravity_obs_data = read_data(gravity_obs, input_dir, dat_tsf = dat_tsf)
    # set same column name for joining datasets
    colnames(gravity_obs_data)[2] = "gmod"
    # subtract starting value 
    gravity_obs_data = mutate(gravity_obs_data, gmod = gmod - min(gmod))
    # subtract mean value 
    # gravity_obs_data$value = gravity_obs_data$value - mean(gravity_obs_data$value, na.rm = T)

    # time steps of gravity observations
    SG_timesteps = data.frame(Timestep = 1:length(gravity_obs_data$gmod[-1]), gmod = gravity_obs_data$gmod[-1])

    ## load modeled gravity response
    g_mod_data = load(file = paste0(output_dir, gravity_mod))
    gravity_mod_data = get(g_mod_data)

    ## create measurement uncertainty of +- 10 %
    # this will be due to unprecise estimated total water mass used in experiment
    gravity_mod_data_over = gravity_mod_data
    gravity_mod_data_over$gmod = gravity_mod_data_over$gmod * 1.1
    gravity_mod_data_under = gravity_mod_data
    gravity_mod_data_under$gmod = gravity_mod_data_under$gmod * 0.9 
    
    ## combine datasets
    gmod = rbind(
             cbind(SG_timesteps, Scenario="Observed gravity response"),
    	     cbind(gravity_mod_data, Scenario="Infiltration model"),
    	     cbind(gravity_mod_data_over, Scenario="uncertainty"),
    	     cbind(gravity_mod_data_under, Scenario="uncertainty")
             )
    # adjust factor levels for plotting order
    gmod$Scenario = factor(gmod$Scenario, levels=c("uncertainty","Observed gravity response", "Infiltration model"))

    # plot
    gravity_ts. = ggplot(gmod, aes(x=Timestep, y=gmod, colour=Scenario)) + geom_line(size=1.5) + 
    	ylab("Gravity [nm/s²]") + xlab("Time since start of experiment [min]") +
        scale_color_manual(values = c("lightgrey","red","blue"), breaks=c("Observed gravity response", "Infiltration model")) + 
        theme(
         legend.position = "bottom",
         legend.title = element_blank(),
    	 legend.text=element_text(size=18),
    	 panel.background = element_rect(fill="transparent"),
         panel.grid.major = element_line(colour = "black", linetype = "dotted")
         ) +
        guides(colour = guide_legend(nrow = 2))
    
    # save plot
    png(filename = paste0(output_dir, "Gravity_responses.png"),
                      width = 600,
                      height = 400,
                      res = 150)
    print(gravity_ts.gg)
    dev.off()
    # print on screen
    # plot(gravity_ts.gg)
    # return NULL
    return(gravity_ts.gg)
}

