#' @title test
#'
#' @description test
#'
#' @param test test
#' 
#' @return test
#' 
#' @details missing
#' @references Marvin Reich (2017), mreich@@posteo.de
#' @examples missing

plot_tempWaterDistribution <- function(sm_mod_data, plotvar){
  switch(plotvar,
        valuesum = {
                   data_plot = dplyr::group_by(sm_mod_data, datetime, Depth) %>%
                   dplyr::summarize(value = sum(value, na.rm = T))
                   },
        mean = {
                   data_plot = dplyr::group_by(sm_mod_data, datetime, Depth) %>%
                   dplyr::summarize(value = mean(value, na.rm = T))
                   },
        max = {
                   data_plot = dplyr::group_by(sm_mod_data, datetime, Depth) %>%
                   dplyr::summarize(value = max(value, na.rm = T))
                   }
  )
  # get value ranges
  data_qt = round(stats::quantile(data_plot$value, na.rm = T, names = F), 2)
  # maximal vertical extent
  # min has to be used as maximum extents are given in negative distances (below surface)
  maxy = min(data_plot$Depth, na.rm=T)
  # miny = max(data_plot$Depth, na.rm=T)
  # actual plotting
  data.gg = ggplot(data_plot, aes(x = datetime, y = Depth)) +
               geom_tile(aes(fill = value)) +
               ylim(maxy,0) + 
               # theme(legend.position=leg) +
               ylab("Depth [m]") + xlab("Time since start of experiment") + 
               # scale_fill_gradientn(colours = rev(viridis(7))) + 
               scale_fill_gradientn(breaks = data_qt, colours=rev(viridis(7)), na.value="red") +
               labs(fill = expression(Delta * "Soil moisture [%VWC]")) +
               theme(legend.position ="bottom",
           	      legend.text=element_text(size=17),
           	      legend.title=element_text(size=19),
           	      panel.grid.major = element_line(colour = "black", linetype = "dotted"),
           	      panel.grid.minor = element_line(colour = "black", linetype = "dotted"))
  return(data.gg) 
}
