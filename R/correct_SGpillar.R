#' @title Correct gravity component grid for the SG pillar
#'
#' @description Removes gravity components at location of the SG pillar.
#'  It can correct for both rectangular or circular shaped pillars, 
#'  depending on what data is supplied.
#'
#' @param gravity_gcomp3d Data.frame, of the gravity component grid containing columns x, y, z, zgrid and gcomp.
#' @param Pillar_x,y,z Numeric, coordinate of the SG pillar (x and y).
#' @param correct_radius Numeric, radius of the SG pillar.
#' @param correct_depth Numeric, vertical extent of the SG pillar.
#' @param SG_x,y Numeric, coordinate of the SG (x and y).
#' @param grid_discretization test
#' 
#' @return A data.frame of the corrected gravity component grid is return.
#' 
#' @details missing
#' @references Marvin Reich (2017), mreich@@posteo.de
#' @examples missing

correct_SGpillar = function(
            gravity_comp3d,
            Pillar_x = NA,
            Pillar_y = NA,
            Pillar_z = NA,
            correct_radius = NA,
            correct_depth = NA,
            SG_X = SG_x,
            SG_Y = SG_y,
            grid_discretization = NA
){
    # gravity_gcomp3d = gravity_component_grid3d
    # Pillar_x = Building_SGpillar_x
    # Pillar_y = Building_SGpillar_y
    # Pillar_z = Building_SGpillar_z
    # correct_radius = NA,
    # correct_depth = NA,
    # SG_X = SG_x,
    # SG_Y = SG_y,
    # grid_discretization = grid3d_discr

    # check if circular correction data is supplied
    if(is.na(correct_radius) & is.na(correct_depth)){
        circular = FALSE
    }else{
        circular = TRUE
    }
    
    # start correction

    if(circular){
      # circular shape
      gcomp_grid3d = dplyr::mutate(gravity_comp3d,
            gcomp = ifelse((x - SG_X)^2 + (y - SG_Y)^2 < correct_radius & Depth >= correct_depth, 0, gcomp))
    }else{
      # rectangular shape
      # construct faces
      faces = as.matrix(data.frame(
              x = c(1, 7, 1, 6, 1, 4, 2, 8, 3, 8, 5, 8),
              y = c(5, 3, 2, 5, 3, 2, 4, 6, 7, 4, 6, 7),
              z = c(3, 5, 5, 2, 2, 3, 6, 4, 4, 7, 7, 6)
      ))
    
      # round z values in case of UTM coordinates
      round_x = decimalplaces(grid_discretization$y)
      round_y = decimalplaces(grid_discretization$x)
      round_z = decimalplaces(grid_discretization$z)
    
      # construct vertices
      # for SG pillar
      vert_SGpillar = construct_vertices(
                          x_cords = Pillar_x,
                          y_cords = Pillar_y,
                          z_cords = Pillar_z
      )
    
      # correct rectangular SG pillar
      gcomp_grid3d = gravity_comp3d %>%
          dplyr::mutate(SGpillar = pip3d(vert_SGpillar, faces, as.matrix(gravity_gcomp3d[,1:3]))) %>%
          dplyr::mutate(gcomp = ifelse(SGpillar >= 0, 0, gcomp)) %>%
          dplyr::select(-SGpillar)
    }
    
    return(gcomp_grid3d)
}

