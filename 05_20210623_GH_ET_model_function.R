# MLFischer 2020-08
# Fischer et al., 2021
# A Phytolith supported Biosphere-Hydrosphere Predictive Model for southern Ethiopia: 
# insights into paleoenvironmental changes and human landscape preferences since the Last Glacial Maximum
###
# ET Function based on Blodgett et al., 1997; Lenters and Cook, 1999; Fischer et al., 2020
###

ETa <- function(t_a_ld    = 291.5, 
                emis_ld   = 0.96, 
                albedo_ld = 0.136,
                rh        = 0.58,
                f_ld      = 0.25,
                cc        = 0.54,
                ws        = 1.42,
                a         = 0.39,
                b         = 0.38,
                a_2       = 0.22,
                b_2       = 2.00,
                cds       = 0.0076,
                p         = 81735,
                r_swc     = 415){
  #t_a_ld    = 291.5     #abaya  # Air Temperature in [K] modern conditions mean
  # Average 24.3 deg C from gridded values (1)
  #r_sw     = CALC      # Net short wave radiation down  
  #r_lup    = CALC      # Net long wave radiation up
  #r_ldw    = CALC      # Net long wave radiation down
  #emis_ld   = 0.96;     #abaya  # Surface emissivity (2) *OK*
  #h        = CALC      # Sensible heat rate
  l         = 2.45E6;   #FIX Latent heat of vaporization in [J kg-1] (4) (FIX)
  #r_swc     = 415;      #  # Cloud-free r_sw in [W/m2] (6,7)
  #albedo_ld = 0.136;    #abaya  # Albedo over land in [%] aus MODIS mean of SW, VIS, NIR
  sigma     = 5.67E-8;  #FIX # Stefan Boltzmann constant in [kg s-3 K-4] (FIX)
  #rh        = 0.58      #abaya 0.57 # Relative humidity in [#/100] (1) *OK*
  #f_ld      = 0.25;     #set schätzwert  # Moisture avail function over land in [#/100] (8)
  #es       = CALC      # Saturation vapour pressure
  #cc        = 0.54;     #abaya  # Cloud cover in [#/100] (2) *OK* http://www.earthenv.org/cloud
  #a         = 0.39;     # Short wave cloud parameters
  #b         = 0.38;     # Short wave cloud parameters
  #a_2       = 0.22;     # Long wave cloud parameters
  #b_2       = 2.00;     # Long wave cloud parameters
  #ws        = 1.42;     #abaya  # Wind speed in 10 m above ground (3.2+/-0.1 m/s) (1)
  #cds       = 0.0076;   #set  # Surface drag coefficient (9)
  cp        = 1005.0;   #FIX # Specific heat capacity of dry air [kJ/kg K] (FIX)
  #p         = 81735     #abaya  # Air pressure in [Pa] (5) 
  r         = 287.0;    #FIX # Gas constant for dry air in [J/K kg] (FIX)
  mmyr      = 3.1536E7  #FIX # Converts evaporation rate from [kg m-^2 s^-1] to
  
  f_rsw_out<-function (){
    
    #global r_swc      % Cloud-free r_sw
    #global albedo_ld  % Albedo over land
    #global cc         % Cloud cover
    #global a          % Short wave cloud parameters
    #global b          % Short wave cloud parameters
    
    r_sw_out = (r_swc * (1 - albedo_ld) * (1 - (a+b*cc)*cc));
    
    return(r_sw_out)
    
  }
  
  f_rlu_out <- function(t){
    
    #global emis_ld    % Surface emissivity
    #global sigma      % Stefan Boltzmann constant
    
    rlu_out=emis_ld*sigma*(t^4);
    
    return(rlu_out)
  }
  
  f_rld_out <-function(t){
    
    #global sigma      % Stefan Boltzmann constant
    #global rh         % Relative humidity
    #global cc         % Cloud cover
    #global a_2        % Long wave cloud parameters
    #global b_2        % Long wave cloud parameters
    
    rld_out=1.24*f_nroot_out(((rh*f_es_out(t))/(100.0*t)),7.0)* sigma*(1+a_2*(cc^b_2))*(t^4.0)
    
    
    return(rld_out)
  }
  
  f_nroot_out <-function(radiant, n){
    
    # Routine calculates nth root of an radiant
    
    nroot_out = radiant^(n^ (-1))
    
    return(nroot_out)
  }
  
  f_es_out<- function(t){
    
    # A few constants
    
    a0 = 6984.505294;
    a1 = -188.903931;
    a2 = 2.133357675;
    a3 = -1.288580973E-2;
    a4 = (4.393587233E-5);
    a5 = (-(8.023923082E-8));
    a6 = (6.136820929E-11);
    
    # The equation
    
    es_out=(100.0*(a0+t*(a1+t*(a2+t*(a3+t*(a4+t*(a5+t*a6)))))));
    
    return(es_out)
  }
  
  
  
  f_rld_out <-function(t){
    
    #global sigma      % Stefan Boltzmann constant
    #global rh         % Relative humidity
    #global cc         % Cloud cover
    #global a_2        % Long wave cloud parameters
    #global b_2        % Long wave cloud parameters
    
    rld_out=1.24*f_nroot_out(((rh*f_es_out(t))/(100.0*t)),7.0)* sigma*(1+a_2*(cc^b_2))*(t^4.0)
    
    
    return(rld_out)
  }
  
  
  
  f_heat_out<- function(t){
    
    #global t_a_ld     % Air Temperature in [K] modern conditions mean
    #global ws         % Wind speed 
    #global cds        % Surface drag coefficient
    #global cp         % Specific heat of dry air 
    #global p          % Air pressure 
    #global r          % Gas constant for dry air 
    
    heat_out=(((p*cds*ws*cp)*(t-t_a_ld))/(r*t_a_ld));
    
    return(heat_out)
  }
  
  
  f_eva_out <- function(t){
    
    #global t_a_ld     % Air Temperature in [K] modern conditions mean
    #global rh         % Relative humidity
    #global f_ld       % Moisture availability function over land
    #global ws         % Wind speed 
    #global cds        % Surface drag coefficient
    #global r          % Gas constant for dry air 
    
    eva_out=(((0.622*cds*ws*f_ld)*(f_es_out(t)-rh*f_es_out(t_a_ld)))/(r*t_a_ld));
    
    return(eva_out)
  }
  
  
  
  
  f_net_out <-function(t){
    
    #global t_a_ld     % Air Temperature in [K] modern conditions mean
    #global emis_ld    % Surface emissivity
    #global l          % Latent heat of vaporization
    
    net_out=(f_rsw_out()-f_rlu_out(t)+emis_ld*f_rld_out(t_a_ld)-f_heat_out(t)-l*f_eva_out(t));
    
    return(net_out)
  }
  
  t_s_ld = as.numeric(uniroot(f_net_out,c(270,330))[1])
  
  #% Best value for e_ld
  
  e_ld = mmyr * f_eva_out(t_s_ld)
  
  
  return(e_ld)
  
  
}

