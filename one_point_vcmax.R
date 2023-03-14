# one_point_vcmax.R
## code to calculate one point vcmax following equation 3 from De Kauwe et al. (2015)

### calculate atmospheric pressure (Pa) from elevation (z; m)
calc_patm = function(z) {
  
  kPo = 101325   # standard atmosphere, Pa (Allen, 1973)
  kTo = 298.15   # base temperature, K (Prentice, unpublished)
  kL = 0.0065    # temperature lapse rate, K/m (Allen, 1973)
  kG = 9.80665   # gravitational acceleration, m/s**2 (Allen, 1973)
  kR = 8.3143    # universal gas constant, J/mol/K (Allen, 1973)
  kMa = 0.028963 # molecular weight of dry air, kg/mol (Tsilingiris, 2008)
  
  patm = kPo*(1.0 - kL*z/kTo)**(kG*kMa/(kR*kL))
  
  patm
}

### calculate gammastar (Pa)
calc_gammastar_pa = function(temp, z) {
  
  patm = calc_patm(z)
  rat = calc_patm(z) / calc_patm(0)
  
  #gammastar25 = 42.75  # ppm
  gammastar25 = 4.332 * rat  # Pa
  Hgm=37830 # J mol-1
  R = 8.314        # J K-1 mol-1
  O2 = 2.09476e5 # ppm
  O2_0 = O2 * 1e-6 * calc_patm(0)
  O2_z = O2 * 1e-6 * calc_patm(z)
  
  temp_k = 273.15+ temp
  
  gStar_pa = gammastar25*exp((Hgm/R)*(1/298.15-1/temp_k))
  
  gStar_pa
  
}

### calcualte the Michaelis-Menton coefficient (Pa) for Rubisco (Km) from temperature
calc_km_pa = function(temp, z) {
  
  patm = calc_patm(z) 
  rat = patm / calc_patm(0)
  
  R = 8.314        
  O2 = 2.09476e5      
  Kc25 = 41.03 * rat 
  Ko25 = 28210 * rat 
  Hkc = 79430  
  Hko = 36380 
  
  temp_k = 273.15 + temp
  
  Kc_pa =Kc25 * exp(Hkc * ((temp_k - 298.15) / (298.15 * R * temp_k)))
  Ko_pa =Ko25* exp(Hko * ((temp_k - 298.15) / (298.15 * R * temp_k)))
  
  O2_pa = O2 * (1e-6) * patm 
  
  Km_pa = Kc_pa * (1 + O2_pa/Ko_pa)
  
  Km_pa 
  
}

### calculate one-point vcmax (µmol m-2 s-1) from 
### light-saturated photosynthesis (asat; µmol m-2 s-1),  ci (µmol mol-1), 
### temperature (temp; °C), and elevation (z; m)
calc_onepoint_vcmax <- function(asat, ci, temp, z){
  
  patm = calc_patm(z)
  gammastar_pa = calc_gammastar_pa(temp, z)
  km_pa = calc_km_pa(temp, z)
  ci_pa = ci * 1e-6 * patm
  
  mm_term = (ci_pa + km_pa) / (ci_pa - gammastar_pa)
  
  vcmax = asat * (mm_term - 0.015)
  
  vcmax
  
}





