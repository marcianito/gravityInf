#' @title Plot 2d transect of modeled infiltration waster distribution
#'
#' @description Plot a transect in 2d of the 3d modeled infiltration water distribution.
#' One or more specific time steps of the total sprinkling experiment duration can be chosen.
#'
#' @param gravity_mod Character string with the filename of the modeled gravity response data.
#' @param plot_ts Vector of one or more time stamps at which plotting is desired.
#'  They have to be supplied in the unit "XX".
#' @param input_dir Character string, containing the path to input file "gravity_obs".
#' @param output_dir Character string, containing the path to the directory where to save the plot.
#' 
#' @return Returns NULL as object. It creates a plot of all input time series.
#' This plot is shown on screen as well as saved to the specified output folder.
#' 
#' @details missing
#' @references Marvin Reich (2017), mreich@@posteo.de
#' @examples missing

plot_transect_2d = function(
            gravity_mod = "gravity_response_modeled.rData",
            plot_ts,
            input_dir = dir_input,
            output_dir = dir_output,
            ...
){
    # input_dir = dir_input

    ## load modeled gravity response
    g_mod_data = load(file = paste0(output_dir, gravity_mod))
    gravity_mod_data = get(g_mod_data)

    # limit data to desired plotting time steps
    data_plot = dplyr::inner_join(plot_ts, gravity_mod_data)
    
    # limit data vertically
    if(limit supplied){
    data_plot = dplyr::filter(data_plot, zgrid <= ??)
    }
    
    transect_2d.gg = ggplot(data_plot, aes(x=xrel, y=zgrid)) +
        geom_tile(aes(fill = value, width = .25)) + 
        # scale_y_continuous(breaks=c(.4,.2,0), trans="reverse") + 
        # scale_x_continuous(breaks=c(0,2,4,6,8)) +
        # scale_fill_gradientn(breaks = c(0.01,0.05,0.1,0.2), colours=rev(viridis(7)), na.value="red") +
        facet_grid(Timestep ~ ., scale="free") +
        ylab("Depth [m]") + xlab("Profile x [m]") + 
        labs(fill = expression(Delta * "Soil moisture [%VWC]")) +
        theme(legend.position ="bottom",
    	      legend.text=element_text(size=17),
    	      legend.title=element_text(size=19),
    	      panel.grid.major = element_line(colour = "black", linetype = "dotted"),
    	      panel.grid.minor = element_line(colour = "black", linetype = "dotted"))
        
    # save plot
    png(filename = paste0(output_dir, "Modeled_2d_transect.png"),
                      width = 600,
                      height = 400,
                      res = 150)
    print(transect_2d.gg)
    dev.off()
    # print on screen
    # plot(gravity_ts.gg)
    # return NULL
    return(gravity_ts.gg)
}

