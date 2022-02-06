########
# LIBS #
########

# Libraries import
import math as mt                                                                                                       # Math lib
# Project personal libraries import
import libs.data_analysis_lib as da                                                                                     # Data analysis lib

########
# VARS #
########

# Heat-exchanger physical properties called by heat-exchanger class constructor
_eff_len_m = 0.680                                                                                                      # Heat-exchamher effective length [m]
_glass_pipe_id_m = 0.05                                                                                                 # Heat-exchamher glass pipe internal-diameter [m]
_steel_pipes_id_m = 0.008                                                                                               # Heat-exchamher steel pipes internal-diameter [m]
_steel_pipes_ed_m = 0.01                                                                                                # Heat-exchamher steel pipes external-diameter [m]
_steel_pipes_num = 5                                                                                                    # Number of heat-exchamher steel pipes

########
# DEFS #
########

# Heat-exchanger class
class He:                                                                                                               # Heat-exchanger class (attributes, constructor, methods)
  eff_len_m = 0.0                                                                                                       # Heat-exchanger effective length [m]
  glass_pipe_id_m = 0.0                                                                                                 # Heat-exchanger glass pipe internal-diameter [m]
  steel_pipes_id_m = 0.0                                                                                                # Heat-exchanger steel pipes internal-diameter [m]
  steel_pipes_ed_m = 0.0                                                                                                # Heat-exchanger steel pipes external-diameter [m]
  steel_pipes_num = 0.0                                                                                                 # Number of heat-exchanger steel pipes
  steel_pipes_thick_m = 0.0                                                                                             # Heat-exchanger steel pipes thickness [m]
  steel_pipes_is_m2 = 0.0                                                                                               # Heat-exchanger steel pipes internal-surface [m^2]
  steel_pipes_es_m2 = 0.0                                                                                               # Heat-exchanger steel pipes external-surface [m^2]
  def __init__(self):                                                                                                   # Constructor
    self.eff_len_m = _eff_len_m                                                                                         # Heat-exchamher effective length [m] init
    self.glass_pipe_id_m = _glass_pipe_id_m                                                                             # Heat-exchamher glass pipe internal-diameter [m] init
    self.steel_pipes_id_m = _steel_pipes_id_m                                                                           # Heat-exchamher steel pipes internal-diameter [m] init
    self.steel_pipes_ed_m = _steel_pipes_ed_m                                                                           # Heat-exchamher steel pipes external-diameter [m] init
    self.steel_pipes_num = _steel_pipes_num                                                                             # Number of heat-exchamher steel pipes init
    self.steel_pipes_thick_m = self.steel_pipes_ed_m-self.steel_pipes_id_m                                              # Heat-exchamher steel pipes thickness [m] init
    self.steel_pipes_is_m2 = cyl_lat_surf(self.steel_pipes_num, self.steel_pipes_id_m, self.eff_len_m)                  # Heat-exchamher steel pipes internal-surface [m^2] init
    self.steel_pipes_es_m2 = cyl_lat_surf(self.steel_pipes_num, self.steel_pipes_ed_m, self.eff_len_m)                  # Heat-exchamher steel pipes external-surface [m^2] init
    return                                                                                                              # Return nothing
  def print_info(self, dbg_flg):                                                                                        # Info printing method with debug flag
    if (dbg_flg):                                                                                                       # If dbg flg is ena
      print("--> Heat-exchanger effective length: "+str(self.eff_len_m)+" [m]")                                         # Print dbg fbk
      print("--> Heat-exchanger glass pipe internal-diameter: "+str(self.glass_pipe_id_m)+" [m]")                       # Print dbg fbk
      print("--> Heat-exchanger steel pipes internal-diameter: "+str(self.steel_pipes_id_m)+" [m]")                     # Print dbg fbk
      print("--> Heat-exchanger steel pipes external-diameter: "+str(self.steel_pipes_ed_m)+" [m]")                     # Print dbg fbk
      print("--> Number of heat-exchanger steel pipes: "+str(self.steel_pipes_num))                                     # Print dbg fbk
      print("--> Heat-exchanger steel pipes thickness: "+str(self.steel_pipes_thick_m)+" [m]")                          # Print dbg fbk
      print("--> Heat-exchanger steel pipes internal-surface: "+str(self.steel_pipes_is_m2)+" [m^2]")                   # Print dbg fbk
      print("--> Heat-exchanger steel pipes external-surface: "+str(self.steel_pipes_es_m2)+" [m^2]")                   # Print dbg fbk
    return                                                                                                              # Return nothing

##########
# FUNCTS #
##########

# Function definition to convert temperatures from [K] to [°C]
def conv_temp_k_c(temp_array):                                                                                          # conv_temp_k_c(Temperatures array [K])
  conv_factor = -273.15                                                                                                 # Conversion factor from [K] to [°C]
  temp_array = temp_array+conv_factor                                                                                   # Temperature-vals conversion
  return temp_array                                                                                                     # Return temperatures array [°C]

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

# Function definition to calculate log-mean temperature difference (LMTD)
def lmtd(delta1, delta2):                                                                                               # lmtd(DeltaTemperature1, DeltaTemperature2)
  return (delta1-delta2)/mt.log((delta1/delta2), mt.e)                                                                  # Return log-mean temperature difference (LMTD)

# Function definition to calculate thermal power exchange
def therm_pow(mass_flow_rate_kg_s, cp_kj_kg_c, delta_temp):                                                             # therm_pow(Mass flow rate [kg/s], Specific heat at const-pressure [kJ/(kg*K)], Delta temperature [°C] or [K])
  return mass_flow_rate_kg_s*cp_kj_kg_c*delta_temp                                                                      # Return thermal power [kJ/s]=[kW]  -->  Q-L'=delta(H) with L'=0  =>  Q=delta(H)=m.*cp*delta(T)

# Function definition to calculate global heat transfer coefficient (global HTC)
def glob_htc_coeff(thermal_pow_kw, surf_m2, lmtd):                                                                      # glob_htc_coeff(Thermal power [kW], Surface [m^2], Log-mean temperature difference [°C] or [K])
  return thermal_pow_kw/(surf_m2*lmtd)                                                                                  # Return global heat transfer coefficient (global HTC) [kW/(m^2*K)]

## LAST
# Function definition to print measures calcs results
def print_measures_calcs_res(measures, dbg_flg):                                                                        # print_measures_calcs_res(Measures list, Debug flag)
  idx = 0                                                                                                               # Measure index
  for meas in measures:                                                                                                 # Measures scrollin' cycle
    print("\n--> "+da.meas_names[idx]+" measure calculations results:")                                                 # Print dbg fbk
    meas.print_info(dbg_flg)                                                                                            # Print measures calcs results (if debug flag is enabled) by callin' the print-info method of the class
    idx += 1                                                                                                            # Measure index upd
  return                                                                                                                # Return nothing
