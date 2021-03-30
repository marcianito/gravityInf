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

plot_transects_2d = function(
            soilmoisture_mod,
            plot_int,
            y_pos,
            n_profiles = NA,
            vert_limit = NA,
            # input_dir = dir_input,
            output_dir = dir_output,
            ...
){
    # soilmoisture_mod = "model_output/Infiltration_model_output_9.rData"
    # output_dir = dir_output
    # plot_int = plot_interval
    # y_pos = SG_y
    # output_dir = "/home/mreich/Dokumente/written/beregnungsPaper/data/output/irrigation/plots/model_profiles/"
    # plot_int = 60
    # y_pos = -7
    # n_profiles = 10
    # soilmoisture_mod = "Infiltration_model_output_node161_1_2018-08-08_1502.rData"

    ## load modeled gravity response
    sm_mod_data = load(file = paste0(output_dir, soilmoisture_mod))
    sm_mod_data = get(sm_mod_data)

    # limit data to desired plotting time steps
    plot_ts = data.frame(
                # datetime = seq(sm_mod_data$datetime[1], sm_mod_data$datetime[length(sm_mod_data$datetime)], by = plot_int)
                datetime = seq(plot_int, sm_mod_data$datetime[length(sm_mod_data$datetime)], by = plot_int)
    )
    data_plot = dplyr::inner_join(plot_ts, sm_mod_data)
    
    # limit data vertically
    if(!is.na(vert_limit)){
    data_plot = dplyr::filter(data_plot, Depth >= vert_limit)
    }

if(!is.na(n_profiles)){
    num_tot_profiles = length(unique(data_plot$y))
    num_dis_profiles = trunc(num_tot_profiles / n_profiles)
    select_profiles = seq(1,num_tot_profiles, by = num_dis_profiles)
    select_ypos = unique(data_plot$y)[select_profiles]
    
    for(y_pos in select_ypos){
        # limit data to one y profile
        print(y_pos)
        data_plot_pos = dplyr::filter(data_plot, y == y_pos)
    
        # get data value ranges
        # data_min = min(data_plot$value, na.rm = T)
        # data_max = max(data_plot$value, na.rm = T)
        # data_mean = mean(data_plot$value, na.rm = T)
        # data_qt = round(stats::quantile(data_plot$value * 100, na.rm = T, names = F), 2)
        data_qt = trunc(c(min(data_plot$value * 100), max(data_plot$value * 100)))
        
        transect_2d.gg = ggplot(data_plot_pos, aes(x=x, y=Depth)) +
            geom_tile(aes(fill = value *100)) + #, width = .25) + 
            scale_fill_gradientn(breaks = data_qt, colours=rev(viridis(7)), na.value="red",
                                 limits = c(min(data_plot$value * 100), max(data_plot$value * 100))) +
            facet_grid(datetime ~ ., scale="free") +
            ylab("Depth [m]") + xlab("Profile x [m]") + 
            labs(fill = expression(Delta * "Soil moisture [%VWC]")) +
            theme(legend.position ="bottom",
        	      legend.text=element_text(size=15),
        	      legend.title=element_text(size=15),
        	      panel.grid.major = element_line(colour = "black", linetype = "dotted"),
        	      panel.grid.minor = element_line(colour = "black", linetype = "dotted"))
            
        # save plot
        png(filename = paste0(output_dir, "Modeled_2d_transect_ypos", y_pos, ".png"),
                          width = 600,
                          height = 100 * length(plot_ts$datetime),
                          res = 150)
        print(transect_2d.gg)
        dev.off()
}

}else{
    data_plot = dplyr::filter(data_plot, y == y_pos)

    # get data value ranges
    # data_min = min(data_plot$value, na.rm = T)
    # data_max = max(data_plot$value, na.rm = T)
    # data_mean = mean(data_plot$value, na.rm = T)
    data_qt = round(stats::quantile(data_plot$value, na.rm = T, names = F), 2)
    
    transect_2d.gg = ggplot(data_plot, aes(x=x, y=Depth)) +
        geom_tile(aes(fill = value)) + #, width = .25) + 
        scale_fill_gradientn(breaks = data_qt, colours=rev(viridis(7)), na.value="red") +
        facet_grid(datetime ~ ., scale="free") +
        ylab("Depth [m]") + xlab("Profile x [m]") + 
        labs(fill = expression(Delta * "Soil moisture [%VWC]")) +
        theme(legend.position ="bottom",
    	      legend.text=element_text(size=17),
    	      legend.title=element_text(size=19),
    	      panel.grid.major = element_line(colour = "black", linetype = "dotted"),
    	      panel.grid.minor = element_line(colour = "black", linetype = "dotted"))
        
    # save plot
    png(filename = paste0(output_dir, "Modeled_2d_transect_ypos", y_pos, ".png"),
                      width = 600,
                      height = 100 * length(plot_ts$datetime),
                      res = 150)
    print(transect_2d.gg)
    dev.off()
}

    return(transect_2d.gg)
}

