###Function for Calculating the Density from Mass Strength and Temperature###
#' Calculate the density from % alcohol by mass
#'
#' @param mass_abm Alcohol % by mass (numeric)
#' @param temperature Temperature in degrees C (numeric)
#' @return a density
#' @export
density_from_mass_abm <- function(mass_abm, temperature) {
  A_comp <- 0
  A_const <- c(9.9820123e2,
               -1.929769495e2,
               3.891238958e2,
               -1.668103923e3,
               1.352215441e4,
               -8.829278388e4,
               3.062874042e5,
               -6.138381234e5,
               7.470172998e5,
               -5.478461354e5,
               2.234460334e5,
               -3.903285426e4)
  for(ka in 2:12){
    A_comp <- A_comp + A_const[ka] * (mass_abm/100)^(ka-1)
  }
  B_comp <- 0
  B_comp <- 0
  B_const <- c(-2.0618513e-1,
               -5.2682542e-3,
               3.6130013e-5,
               -3.8957702e-7,
               7.1693540e-9,
               -9.9739231e-11)
  for(kb in 1:6){
    B_comp <- B_comp + B_const[kb] * (temperature - 20)^kb
  }
  C_comp <- 0
  C_const <- list(
    c(1.693443461530087e-1,
      -1.046914743455169e1,
      7.196353469546523e1,
      -7.047478054272792e2,
      3.924090430035045e3,
      -1.210164659068747e4,
      2.248646550400788e4,
      -2.605562982188164e4,
      1.852373922069467e4,
      -7.420201433430137e3,
      1.285617841998974e3),
    c(-1.193013005057010e-2,
      2.517399633803461e-1,
      -2.170575700536993,
      1.353034988843029e1,
      -5.029988758547014e1,
      1.096355666577570e2,
      -1.422753946421155e2,
      1.080435942856230e2,
      -4.414153236817392e1,
      7.442971530188783),
    c(-6.802995733503803e-4,
      1.876837790289664e-2,
      -2.002561813734156e-1,
      1.022992966719220,
      -2.895696483903638,
      4.810060584300675,
      -4.672147440794683,
      2.458043105903461,
      -5.411227621436812e-1),
    c(4.075376675622027e-6,
      -8.763058573471110e-6,
      6.515031360099368e-6,
      -1.515784836987210e-6),
    c(-2.788074354782409e-8,
      1.345612883493354e-8)
  )
  for (n in seq_along(C_const)){
    for (m in seq_along(C_const[[n]])) {
      C_comp <- C_comp + C_const[[n]][m] * (mass_abm/100)^m * (temperature - 20)^n
    }
  }
  density <- A_const[1] + A_comp + B_comp + C_comp
  return(density)
}

###Function for calculating the density from volume strength and temperature - NB uses mass strength function###
#' Calculate density and alcohol mass fraction from alcohol % by volume
#'
#' @param vol_abv Alcohol % by volume (numeric)
#' @param temperature temperature in degrees C (numeric)
#' @return a list with density at 20 degrees and mass fraction
#' @export
density_from_vol_abv <- function(vol_abv, temperature){
  LL <- 0
  UL <- 1
  mass_frac <- 0.5
  while((UL - LL) > 0.0000001) {
    mass_frac <- (LL + UL)/2
    density <- density_from_mass_abm(mass_frac*100, 20)
    abv_vol <- (density * mass_frac)/789.24*100
    if(abv_vol < vol_abv){LL <- mass_frac} else {UL <- mass_frac}
  }
  density <- density_from_mass_abm(mass_frac*100, temperature)
  return(list(Density_20degrees = density,
              Mass_percent = round(mass_frac*100,2)))
}
