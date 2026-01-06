# packages that might not be installed - install any not already installed
packages <- c("devtools", "rfishbase", "gsw")
install.packages(setdiff(packages, rownames(installed.packages())))
library(dplyr)
library(devtools)
library(gsw)
library(readxl)
library(stringr)
library(lubridate)
library(sf)

calc_po2_sat <- function(salinity, temp, depth, oxygen, lat, long, umol_m3, ml_L) {
  # Input:       S = Salinity (pss-78)
  #              T = Temp (deg C) ! use potential temp
  #depth is in meters

  #Pena et al. ROMS and GLORYS are in mmol per m^3 (needs to be converted to umol per kg) (o2 from trawl data was in mL/L, so had to do extra conversions)
  gas_const = 8.31
  partial_molar_vol = 0.000032
  kelvin = 273.15
  boltz = 0.000086173324

  #convert mmol to umol
   umol_m3 <- oxygen*1000
  #convert m3 to l
   umol_l <- umol_m3/1000
  #convert from molality (moles per volume) to molarity (moles per mass)
   #1 L of water = 1 kg of water, so no equation needed?? Right??
   o2_umolkg <- umol_l * 1/1

  SA = gsw_SA_from_SP(salinity,depth,long,lat) #absolute salinity for pot T calc
  pt = gsw_pt_from_t(SA,temp,depth) #potential temp at a particular depth
  #this is for if using data that has oxygen in ml/L
  #CT = gsw_CT_from_t(SA,temp,depth) #conservative temp
  #sigma0 = gsw_sigma0(SA,CT)
 # o2_umolkg = oxygen*44660/(sigma0+1000)

  O2_Sat0 = gsw_O2sol_SP_pt(salinity,pt)

  #= o2satv2a(sal,pt) #uses practical salinity and potential temp - solubity at p =1 atm
  press = exp(depth*10000*partial_molar_vol/gas_const/(temp+kelvin))
  O2_satdepth = O2_Sat0*press

  #solubility at p=0
  sol0 = O2_Sat0/0.209
  sol_Dep = sol0*press
  po2 = o2_umolkg/sol_Dep
  po2 <- po2 * 101.325 # convert to kPa
  return(po2)

}

# calc o2 solubility, relies on o2 in umol/kg
gsw_O2sol_SP_pt <- function(sal,pt) {
  x = sal
  pt68 = pt*1.00024
  y = log((298.15 - pt68)/(273.15 + pt68))

  a0 =  5.80871
  a1 =  3.20291
  a2 =  4.17887
  a3 =  5.10006
  a4 = -9.86643e-2
  a5 =  3.80369
  b0 = -7.01577e-3
  b1 = -7.70028e-3
  b2 = -1.13864e-2
  b3 = -9.51519e-3
  c0 = -2.75915e-7

  O2sol = exp(a0 + y*(a1 + y*(a2 + y*(a3 + y*(a4 + a5*y)))) + x*(b0 + y*(b1 + y*(b2 + b3*y)) + c0*x))
  return(O2sol)
}

# calc o2 solubility, relies on o2 in umol/kg
calc_o2_sol <- function(sal, temp, depth, long, lat) {
  SA = gsw_SA_from_SP(sal,depth,long,lat) #absolute salinity for pot T calc
  pt = gsw_pt_from_t(SA,temp,depth) #potential temp at a particular depth
  pt68 = pt*1.00024
  y = log((298.15 - pt68)/(273.15 + pt68))

  a0 =  5.80871
  a1 =  3.20291
  a2 =  4.17887
  a3 =  5.10006
  a4 = -9.86643e-2
  a5 =  3.80369
  b0 = -7.01577e-3
  b1 = -7.70028e-3
  b2 = -1.13864e-2
  b3 = -9.51519e-3
  c0 = -2.75915e-7

  O2sol = exp(a0 + y*(a1 + y*(a2 + y*(a3 + y*(a4 + a5*y)))) + sal*(b0 + y*(b1 + y*(b2 + b3*y)) + c0*sal))
  return(O2sol)
}

calc_sigma <- function(s, t, p) {
  constants <- list(B0 = 8.24493e-1, B1 = -4.0899e-3, B2 = 7.6438e-5, B3 = -8.2467e-7, B4 = 5.3875e-9, C0 = -5.72466e-3, C1 = 1.0227e-4, C2 = -1.6546e-6, D0 = 4.8314e-4, A0 = 999.842594, A1 = 6.793952e-2, A2 = -9.095290e-3, A3 = 1.001685e-4, A4 = -1.120083e-6, A5 = 6.536332e-9, FQ0 = 54.6746, FQ1 = -0.603459, FQ2 = 1.09987e-2, FQ3 = -6.1670e-5, G0 = 7.944e-2, G1 = 1.6483e-2, G2 = -5.3009e-4, i0 = 2.2838e-3, i1 = -1.0981e-5, i2 = -1.6078e-6, J0 =1.91075e-4, M0 = -9.9348e-7, M1 = 2.0816e-8, M2 = 9.1697e-10, E0 = 19652.21, E1 = 148.4206, E2 = -2.327105, E3 = 1.360477e-2, E4 = -5.155288e-5, H0 = 3.239908, H1 = 1.43713e-3, H2 = 1.16092e-4, H3 = -5.77905e-7, K0 = 8.50935e-5, K1 =-6.12293e-6, K2 = 5.2787e-8)

  list2env(constants, environment())
  constant_names <- names(constants)
  t2 = t*t
  t3 = t*t2
  t4 = t*t3
  t5 = t*t4
  #  if (s <= 0.0) s = 0.000001
  s32 = s^ 1.5
  p = p / 10.0 # convert decibars to bars */
  sigma = A0 + A1*t + A2*t2 + A3*t3 + A4*t4 + A5*t5 + (B0 + B1*t + B2*t2 + B3*t3 + B4*t4)*s + (C0 + C1*t + C2*t2)*s32 + D0*s*s
  kw = E0 + E1*t + E2*t2 + E3*t3 + E4*t4
  aw = H0 + H1*t + H2*t2 + H3*t3
  bw = K0 + K1*t + K2*t2
  k = kw + (FQ0 + FQ1*t + FQ2*t2 + FQ3*t3)*s + (G0 + G1*t + G2*t2)*s32 + (aw + (i0 + i1*t + i2*t2)*s + (J0*s32))*p + (bw + (M0 + M1*t + M2*t2)*s)*p*p
  val = 1 - p / k
  sigma = sigma / val - 1000.0
  rm(constant_names)
  return(sigma)
}

# convert ml /l to umol / kg
convert_o2 <- function(o2ml_l, sigma){
  (o2ml_l * 44.660)/((1000 + sigma)/1000)
}

rsq <- function(x, y) cor(x,y) * cor(x,y)
