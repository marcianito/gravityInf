% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_model_inversion.R
\name{run_model_inversion}
\alias{run_model_inversion}
\title{Run infiltration model in inversion mode}
\usage{
run_model_inversion(dtheta_macro_min = dtheta_macro_min,
  dtheta_macro_max = dtheta_macro_max,
  dtheta_macro2_min = dtheta_macro2_min,
  dtheta_macro2_max = dtheta_macro2_max,
  dtheta_other_min = dtheta_other_min, dtheta_other_max = dtheta_other_max,
  mdepth_min = mdepth_min, mdepth_max = mdepth_max,
  mdepth2_min = mdepth2_min, mdepth2_max = mdepth2_max,
  latflow_fac_min = latflow_fac_min, latflow_fac_max = latflow_fac_max,
  dtheta_macro_start = dtheta_macro_start,
  dtheta_macro2_start = dtheta_macro2_start, mdepth_start = mdepth_start,
  mdepth2_start = mdepth2_start, dtheta_other_start = dtheta_other_start,
  latflow_fac_start = latflow_fac_start, output_dir, input_dir, inner_inum,
  ...)
}
\arguments{
\item{test}{test}
}
\value{
A list (?) of the optimized model parameters is returned.
}
\description{
The infiltration model is run wrapped within a
 particle-swarm optimization algorithm, based on R-package 'ppso'.
}
\details{
It can be decided if one, two or none marco pore layer is desired.
 For modifications of the optimization algorithm, variables have to be supplied to the function
 as described in ?optim_dds.
}
\examples{
missing
}
\references{
Marvin Reich (2017), mreich@posteo.de
}
