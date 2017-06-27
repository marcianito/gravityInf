#' @title Check infiltrated water volumen
#'
#' @description test
#'
#' @param test
#' @param test
#' @param test
#' ...
#' @details missing
#' @references Marvin Reich (2017), mreich@@posteo.de
#' @examples missing

checkWaterVolumes = function(sm_mod_data,
                              whichCheck,
                              inf_time,
                              cell_volume,
                              infWater_timestep,
                              ts_check = 1,
                              plotting = F
){
  # set chosen check-method
  switch(whichCheck,
         totalwater = {
                        modelwater = dplyr::summarize(sm_mod_data, modelsum = sum(value, na.rm=T)) 
                        comparison = data.frame(water_model = as.integer(modelwater) * cell_volume,
                                                water_real = sum(seq(1,inf_time)) * infWater_timestep
                                               )
                      },
         waterpertimestep = {
                        modelwater = dplyr::group_by(sm_mod_data, datetime) %>%
                                     dplyr::summarize(modelsum = sum(value/datetime, na.rm=T)) 
                        comparison = data.frame(timestep = modelwater$datetime,
                                                water_model = modelwater$modelsum * cell_volume,
                                                water_real = infWater_timestep
                                               )
                        if(plotting){
                            plot(ggplot(modelwater, aes(x=datetime, y=modelsum)) + geom_line())
                        }
                      },
         tswater = {
                        modelwater = dplyr::filter(sm_mod_data, datetime == ts_check)
                        comparison = data.frame(water_model = sum(modelwater$value, na.rm = T) * cell_volume,
                                                water_real = ts_check * infWater_timestep
                                               )
                      }
  )
  return(comparison)
}
