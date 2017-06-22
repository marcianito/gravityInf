#' @title Create an interpolated grid of the superficial water intensity distribution
#'
#' @description Based on the supplied model domain and measured (or constructed) water intensities,
#'  a surface grid is interpolated for further usage as infiltration model input.
#'
#' @param test test
#' 
#' @return A data.frame, containing spatial information (x, y) and a scale factor of the water intensity for each grid cell.
#' 
#' @details missing
#' @references Marvin Reich (2017), mreich@@posteo.de
#' @examples missing

create_WaterIntensityGrid = function(
            input_file,
            intp_method,
            surface,
            zerosBorder,
            spat_col,
            dat_col,
            UTM_gridcenter_x = NA,
            UTM_gridcenter_y = NA,
            input_dir = dir_input, 
            output_dir = dir_output,
            ...
){

    # input_file = "waterIntensity_measured.csv"
    # input_file = "waterIntensity_measured.rData"
    # intp_method = "IDW"
    # surface = gravity_component_grid3d
    # zerosBorder = .5
    # spat_col = spatial_col
    # dat_col = data_col
    # UTM_gridcenter_x = SG_x
    # UTM_gridcenter_y = SG_y
    # input_dir = dir_input
    # output_dir = dir_output

    # load measured (or set up) water intensity distribution data
    intensities_raw = read_data(input_file, input_dir, spat_col, dat_col)
    # if UTM coordinates are used,
    # the UTM center coordinates of the grid have to be ADDED
    # to the relative intensity grid
    # !!
    # !! add option to ALSO supply intensity grid in UTM
    # !!
    if(!is.na(UTM_gridcenter_x) & !is.na(UTM_gridcenter_y)){
    intensities_raw$x = intensities_raw$x + UTM_gridcenter_x
    intensities_raw$y = intensities_raw$y + UTM_gridcenter_y
    }

    # get surface grid out of gravity component grid
    grid_surface = surface %>%
                dplyr::filter(Depth == max(surface$Depth)) %>%
                dplyr::select(x,y)
    
    # number of "zeros" per border-side
    numPoints_side = min(c(length(unique(grid_surface$x)), length(unique(grid_surface$y))))
    # number of points per side
    nzeros_border = round(numPoints_side * zerosBorder,0)
    # nzeros_border = 7 
    
    Zeros = data.frame(x = c(seq(min(grid_surface$x),max(grid_surface$x),length.out=nzeros_border),seq(min(grid_surface$x),max(grid_surface$x),length.out=nzeros_border), rep(min(grid_surface$x),nzeros_border),rep(max(grid_surface$x),nzeros_border)),
                       y = c(rep(min(grid_surface$y),nzeros_border),rep(max(grid_surface$y),nzeros_border),seq(min(grid_surface$y),max(grid_surface$y),length.out=nzeros_border),seq(min(grid_surface$y),max(grid_surface$y),length.out=nzeros_border)),
                       value = 0)
    # subset intensity raw data
    intensities_addZeros = dplyr::select(intensities_raw, x, y, value)
    # combine subsetted intensity data and artificial added zeros
    intensities_Zeros = rbind(intensities_addZeros, Zeros)
    # plot
    # ggplot(intensities_Zeros, aes(x=x,y=y, color=value)) + geom_point()

    # decide interpolation method
    switch(intp_method,
           # Inverse distance weight
           IDW = {
            ## without extra added Zeros
            # idw.gstat = gstat(formula = total_weight_dif ~ 1, locations = ~ xrel + yrel, data = Iintensities, nmax = 10, set = list(idp = 2))
            # with extra Zeros (as stützstellen at the borders of the irrigation area)
            idw.gstat = gstat(formula = value ~ 1, locations = ~ x + y, data = intensities_Zeros, nmax = 10, set = list(idp = 2))
            #idw.gstat = gstat(formula = dataraw ~ 1, locations = ~ x + y + Depth, data = data_in, nmax = 10, nmin=3, maxdist=1.1, set = list(idp = 2))
            data_interpolated <- predict(idw.gstat, grid_surface)[,-4]
            colnames(data_interpolated)[3] = "value"
            ## alternative way:
            # coordinates(grid_surface) = ~ x + y
            # gridded(grid_surface) = TRUE
            # coordinates(intensities_Zeros) = ~ x + y
            # # IDW
            # idw_interp = idw(formula=value ~ 1, locations=intensities_Zeros, newdata=grid_surface)
            # # structure output data
            # data_interpolated = as.data.frame(idw_interp)
            # names(data_interpolated)[1:3] = c("x","y","value")
            # data_interpolated = data_interpolated[,-4]
            
            # plot
            # ggplot(data_interpolated, aes(x=x,y=y, color=value)) + geom_point()

            ## calculate intensities per min out of total_weight_dif
            # Becher_radius = 0.035 #[m]
            num_cell = length(grid_surface$x)
            totalWater = sum(data_interpolated$value, na.rm=T)
            water_avg = totalWater / num_cell
            # intensity as relatin of:
            # cell weight / average cell weight
            intensities_grid_interpolated = dplyr::mutate(data_interpolated, intensity = value / water_avg)
            
            # check sum of weights after interpolation
            # should be equal to number of cells
            # sum(intensities_grid_interpolated$intensity, na.rm=T)
           },
           linear = {
data_interpolated = approx(intensities_Zeros$x, intensities_Zeros$y, grid_surface, method="linear")
install.packages("akima")
library(akima)
data_interpolated = akima::interp(intensities_Zeros, z = intensities_Zeros$value, surface_grid, linear = T)
dups = duplicated(intensities_Zeros[,1:2])
intensities_Zeros = intensities_Zeros[!dups,]
data_interpolated = akima::interp(x = intensities_Zeros$x, y = intensities_Zeros$y, z = intensities_Zeros$value, xo = grid_surface$x, yo = grid_surface$y, linear = T)
dataplot = data.frame(x = data_interpolated[[1]], y = data_interpolated[[2]], value = data_interpolated[[3]])
tt = data_interpolated[[1]]
tt = data_interpolated[[2]]
tt = data_interpolated[[3]]
            # plot
            ggplot(as.data.frame(data_interpolated), aes(x=x,y=y, color=value)) + geom_point()

           },
           # kriging
           krige = {
            grid_surface_krige = cbind(grid_surface, dummy=2, fummy=5)
            gridded(grid_surface_krige) = ~x + y
            coordinates(intensities_Zeros)= ~x + y
            # remove duplicate locations
            intensities_Zeros <- intensities_Zeros[-zerodist(intensities_Zeros)[,1],]
            # variogram
            vario = autofitVariogram(value ~ 1, intensities_Zeros,
                      model = c("Sph", "Exp", "Gau", "Ste", "Mat"),
            	      verbose=T)
            # interpolate data to new grid
            data_pred = krige(value ~ 1, intensities_Zeros, grid_surface_krige, vario$var_model)
            # construct data.frame
            data_interpolated = data.frame(x = coordinates(data_pred)[,1], y = coordinates(data_pred)[,2], value = data_pred$var1.pred)

            ## calculate intensities per min out of total_weight_dif
            # Becher_radius = 0.035 #[m]
            num_cell = length(grid_surface$x)
            totalWater = sum(data_interpolated$value, na.rm=T)
            water_avg = totalWater / num_cell
            # intensity as relatin of:
            # cell weight / average cell weight
            intensities_grid_interpolated = dplyr::mutate(data_interpolated, intensity = value / water_avg)
            
            # check sum of weights after interpolation
            # should be equal to number of cells
            # sum(intensities_grid_interpolated$intensity, na.rm=T)
            
            # plot
            # ggplot(data_interpolated, aes(x=x,y=y, color=value)) + geom_point()

            # # kriging
            # NOT WORKING !!
            # # create semivariogram
            # # from raw data, WITHOUT zeros
            # # coordinates(intensities_raw)= ~ x + y
            # # intensities_raw <- intensities_raw[-zerodist(intensities_raw)[,1],]
            # semivariog = variogram(value~1, locations=intensities_Zeros, data=intensities_Zeros)
            # # semivariog = variogram(value~1, locations=intensities_raw, data=intensities_raw)
            # # read and estimate range, sill and nugget and use to
            # # create variogram model
            # model.variog<-vgm(psill=60000, model="Exp", nugget=0, range=5)
            # model.variog<-vgm(psill=60000, model="Mat", nugget=0, range=5)
            # model.variog<-vgm(psill=60000, model="Pow", nugget=0, range=5)
            # # fit semivariogram to model
            # fit.variog<-fit.variogram(semivariog, model.variog)
            # plot(semivariog, fit.variog)
            # # interpolate using this model
            # data_interpolated = krige(formula=value ~ 1, locations=intensities_Zeros, newdata=grid_surface, model=model.variog)
            # # structure output data
            # data_interpolated = as.data.frame(idw_interp)
            # names(data_interpolated)[1:3] = c("x","y","value")
            # data_interpolated = data_interpolated[,-4]
           }
    )

    return(intensities_grid_interpolated)
}

