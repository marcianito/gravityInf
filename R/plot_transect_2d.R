#' @title Plot 2d transect of modeled infiltration waster distribution
#'
#' @description Plot a transect in 2d of the 3d modeled infiltration water distribution.
#' One or more specific time steps of the total sprinkling experiment duration can be chosen.
#'
#' @param soilmoisture_mod Character string with the filename of the modeled soil moisture data.
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
            sm_mod_data,
            tstep,
            y_pos,
            vert_limit = NA,
            ...
){
    #     soilmoisture_mod = "model_output/Infiltration_model_output_9.rData"
    #     output_dir = dir_output
    #     tstep = 1
    #     y_pos = SG_y

    ## load modeled gravity response
    #     sm_mod_data = load(file = paste0(output_dir, soilmoisture_mod))
    #     sm_mod_data = get(sm_mod_data)

    # limit data to one time step
    data_plot = dplyr::filter(sm_mod_data, datetime == tstep)
    
    # limit data vertically
    if(!is.na(vert_limit)){
    data_plot = dplyr::filter(data_plot, Depth >= vert_limit)
    }

    # limit data to one y profile
    data_plot = dplyr::filter(data_plot, y == y_pos)

    # get data value ranges
    # data_min = min(data_plot$value, na.rm = T)
    # data_max = max(data_plot$value, na.rm = T)
    # data_mean = mean(data_plot$value, na.rm = T)
    data_qt = round(stats::quantile(data_plot$value, na.rm = T, names = F), 2)
    
    transect_2d.gg = ggplot(data_plot, aes(x=x, y=Depth)) +
        geom_tile(aes(fill = value)) + #, width = .25) + 
        scale_fill_gradientn(breaks = data_qt, colours=rev(viridis(7)), na.value="red") +
        # facet_grid(datetime ~ ., scale="free") +
        ylab("Depth [m]") + xlab("Profile x [m]") + 
        labs(fill = expression(Delta * "Soil moisture [%VWC]")) +
        theme(legend.position ="bottom",
              # legend.text=element_text(size=17),
              # legend.title=element_text(size=19),
    	      panel.grid.major = element_line(colour = "black", linetype = "dotted"),
    	      panel.grid.minor = element_line(colour = "black", linetype = "dotted"))
        
    # # save plot
    # png(filename = paste0(output_dir, "Modeled_2d_transect.png"),
    #                   width = 600,
    #                   height = 100 * length(plot_ts$datetime),
    #                   res = 150)
    # print(transect_2d.gg)
    # dev.off()
    # 
    return(transect_2d.gg)
}

