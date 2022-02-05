########
# LIBS #
########

# Libraries import
import math as mt                                                                                                       # Math lib

########
# DEFS #
########

# Data-structure definition to contain heat-exchanger details
class He:                                                                                                               # Heat-exchanger details data-structure class
  eff_len_m: float                                                                                                      # Heat-exchanger effective length [m]
  glass_pipe_id_m: float                                                                                                # Heat-exchanger glass pipe internal-diameter [m]
  steel_pipes_id_m: float                                                                                               # Heat-exchanger steel pipes internal-diameter [m]
  steel_pipes_ed_m: float                                                                                               # Heat-exchanger steel pipes external-diameter [m]
  steel_pipes_num: float                                                                                                # Number of heat-exchanger steel pipes
  steel_pipes_thick_m: float                                                                                            # Heat-exchanger steel pipes thickness [m]
  steel_pipes_is_m2:float                                                                                               # Heat-exchanger steel pipes internal-surface [m^2]
  steel_pipes_es_m2:float                                                                                               # Heat-exchanger steel pipes external-surface [m^2]

##########
# FUNCTS #
##########

# Function definition to calculate cylindrical lateral surface
def cyl_lat_surf(num_cyl, diam, len):                                                                                   # cyl_lat_surf(Number of cylinders, Cylinder diameter, Cylinder length)
  return num_cyl*mt.pi*diam*len                                                                                         # Return cylindrical lateral surface

# Function definition to convert volume flow rate [l/h] into mass flow rate [kg/s]
def vol_flow_rate_to_mass_flow_rate(vol_flow_rate_l_h, density):                                                        # vol_flow_rate_to_mass_flow_rate(Volume flow rate value [l/h], density value)
  conv_factor = 0.001/3600                                                                                              # Conversion factor from [l/h] to [m^3/s]
  return vol_flow_rate_l_h*conv_factor*density                                                                          # Return mass flow rate [kg/s] --> [m^3/s]*[kg/m^3]=[kg/s]

# Function definition to calculate average between two values
def avg(val1, val2):                                                                                                    # avg(Value1, Value2)
  return (val1+val2)/2                                                                                                  # Return average value btwn Value1 and Value2

# Function to calculate log-mean temperature difference (LMTD)
def lmtd(delta1, delta2):                                                                                               # lmtd(DeltaTemperature1, DeltaTemperature2)
  return (delta1-delta2)/mt.log((delta1/delta2), mt.e)                                                                  # Return log-mean temperature difference (LMTD)

# Function definition to calculate thermal power exchange
def therm_pow(mass_flow_rate_kg_s, cp_kj_kg_c, delta_temp):                                                             # therm_pow(Mass flow rate [kg/s], Specific heat at const-pressure [kJ/(kg*K)], Delta temperature [°C] or [K])
  return mass_flow_rate_kg_s*cp_kj_kg_c*delta_temp                                                                      # Return thermal power [kJ/s]=[kW]  -->  Q-L'=delta(H) with L'=0  =>  Q=delta(H)=m.*cp*delta(T)

# Function to calculate global heat transfer coefficient (global HTC)
def glob_htc_coeff(thermal_pow_kw, surf_m2, lmtd):                                                                      # glob_htc_coeff(Thermal power [kW], Surface [m^2], Log-mean temperature difference [°C] or [K])
  return thermal_pow_kw/(surf_m2*lmtd)                                                                                  # Return global heat transfer coefficient (global HTC) [kW/(m^2*K)]
