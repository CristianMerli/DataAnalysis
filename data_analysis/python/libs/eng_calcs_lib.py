########
# LIBS #
########

# Libraries import
import math as mt                                                                                                       # Math lib
import enum as en                                                                                                       # Enum lib
# Project personal libraries import
import libs.data_analysis_lib as da                                                                                     # Data analysis lib
import libs.output_lib as out                                                                                           # Output lib

########
# VARS #
########

# Heat-exchanger physical properties called by heat-exchanger class constructor
_eff_len_m = 0.7                                                                                                        # Heat-exchamher effective length [m]
_glass_pipe_id_m = 0.05                                                                                                 # Heat-exchamher glass pipe internal-diameter [m]
_glass_pipe_diaph_num  = 13                                                                                             # Number of heat-exchanger diaphragms in glass pipe
_glass_pipe_diaph_pace_m = 0.05                                                                                         # Heat-exchanger glass pipe diaphragms pace [m]
_steel_pipes_id_m = 0.008                                                                                               # Heat-exchamher steel pipes internal-diameter [m]
_steel_pipes_ed_m = 0.01                                                                                                # Heat-exchamher steel pipes external-diameter [m]
_steel_pipes_num = 5                                                                                                    # Number of heat-exchamher steel pipes
_steel_pipes_pace_m = 0.015                                                                                             # Heat-exchanger steel pipes pace [m]
# Rounding vars
rnd_dec_places = 4                                                                                                      # Number of decimal places in rounding operations to avoid floating-point errors
# Engineering calcs vars
re_trans = 2300                                                                                                         # Reynolds number transition value between fluid laminar and turbulent flow [adimensional]

########
# DEFS #
########

# Heat-exchanger side enum definition
class He_side(en.Enum):                                                                                                 # Heat-exchanger side enum class
  tube = 1                                                                                                              # Tube side (Hot fluid in cylindrical steel pipes)
  shell = 2                                                                                                             # Shell side (Cold fluid in glass pipe)
# Heat-exchanger class
class He:                                                                                                               # Heat-exchanger class (attributes, constructor, methods)
  eff_len_m = 0.0                                                                                                       # Heat-exchanger effective length [m]
  glass_pipe_id_m = 0.0                                                                                                 # Heat-exchanger glass pipe internal-diameter [m]
  glass_pipe_diaph_num = 0                                                                                              # Number of heat-exchanger diaphragms in glass pipe
  glass_pipe_diaph_pace_m = 0.0                                                                                         # Heat-exchanger glass pipe diaphragms pace [m]
  steel_pipes_id_m = 0.0                                                                                                # Heat-exchanger steel pipes internal-diameter [m]
  steel_pipes_ed_m = 0.0                                                                                                # Heat-exchanger steel pipes external-diameter [m]
  steel_pipes_num = 0                                                                                                   # Number of heat-exchanger steel pipes
  steel_pipes_pace_m = 0.0                                                                                              # Heat-exchanger steel pipes pace [m]
  steel_pipes_thick_m = 0.0                                                                                             # Heat-exchanger steel pipes thickness [m]
  steel_pipes_interspce = 0.0                                                                                           # Heat-exchanger interspace between steel pipes [m]
  steel_pipes_is_m2 = 0.0                                                                                               # Heat-exchanger steel pipes internal-surface [m^2]
  steel_pipes_es_m2 = 0.0                                                                                               # Heat-exchanger steel pipes external-surface [m^2]
  def __init__(self):                                                                                                   # Constructor
    self.eff_len_m = _eff_len_m                                                                                         # Heat-exchamher effective length [m] init
    self.glass_pipe_id_m = _glass_pipe_id_m                                                                             # Heat-exchamher glass pipe internal-diameter [m] init
    self.glass_pipe_diaph_num = _glass_pipe_diaph_num                                                                   # Number of heat-exchanger diaphragms in glass pipe init
    self.glass_pipe_diaph_pace_m = _glass_pipe_diaph_pace_m                                                             # Heat-exchanger glass pipe diaphragms pace [m] init
    self.steel_pipes_id_m = _steel_pipes_id_m                                                                           # Heat-exchamher steel pipes internal-diameter [m] init
    self.steel_pipes_ed_m = _steel_pipes_ed_m                                                                           # Heat-exchamher steel pipes external-diameter [m] init
    self.steel_pipes_num = _steel_pipes_num                                                                             # Number of heat-exchamher steel pipes init
    self.steel_pipes_pace_m = _steel_pipes_pace_m                                                                       # Heat-exchanger steel pipes pace [m] init
    self.steel_pipes_thick_m = round(self.steel_pipes_ed_m-self.steel_pipes_id_m, rnd_dec_places)                       # Heat-exchamher steel pipes thickness [m] init and round value to avoid floating-point errors
    self.steel_pipes_interspce = round(self.steel_pipes_pace_m-self.steel_pipes_ed_m, rnd_dec_places)                   # Heat-exchanger interspace between steel pipes [m] init and value to avoid floating-point errors
    self.steel_pipes_is_m2 = cyl_lat_surf(self.steel_pipes_num, self.steel_pipes_id_m, self.eff_len_m)                  # Heat-exchamher steel pipes internal-surface [m^2] init
    self.steel_pipes_es_m2 = cyl_lat_surf(self.steel_pipes_num, self.steel_pipes_ed_m, self.eff_len_m)                  # Heat-exchamher steel pipes external-surface [m^2] init
    return                                                                                                              # Return nothing
  def print_info_save_out(self, dbg_flg):                                                                               # Info printing and output saving method with debug flag
    if (dbg_flg):                                                                                                       # If dbg flg is ena
      dbg_str = ("\n--> Heat-exchanger effective length: "+str(self.eff_len_m)+" [m]"\
      +"\n--> Heat-exchanger glass pipe internal-diameter: "+str(self.glass_pipe_id_m)+" [m]"\
      +"\n--> Number of heat-exchanger diaphragms in glass pipe: "+str(self.glass_pipe_diaph_num)\
      +"\n--> Heat-exchanger glass pipe diaphragms pace: "+str(self.glass_pipe_diaph_pace_m)+" [m]"\
      +"\n--> Heat-exchanger steel pipes internal-diameter: "+str(self.steel_pipes_id_m)+" [m]"\
      +"\n--> Heat-exchanger steel pipes external-diameter: "+str(self.steel_pipes_ed_m)+" [m]"\
      +"\n--> Number of heat-exchanger steel pipes: "+str(self.steel_pipes_num)\
      +"\n--> Heat-exchanger steel pipes pace: "+str(self.steel_pipes_pace_m)+" [m]"\
      +"\n--> Heat-exchanger steel pipes thickness: "+str(self.steel_pipes_thick_m)+" [m]"\
      +"\n--> Heat-exchanger interspace between steel pipes: "+str(self.steel_pipes_interspce)+" [m]"\
      +"\n--> Heat-exchanger steel pipes internal-surface: "+str(self.steel_pipes_is_m2)+" [m^2]"\
      +"\n--> Heat-exchanger steel pipes external-surface: "+str(self.steel_pipes_es_m2)+" [m^2]\n")                    # Dbg fbk
      print(dbg_str)                                                                                                    # Print dbg fbk
      out.save_output(out.Output_typ.he, dbg_str)                                                                       # Save dbg output
    return                                                                                                              # Return nothing

##########
# FUNCTS #
##########

# Function definition to convert temperatures from [K] to [°C]
def conv_temp_k_c(temp_array):                                                                                          # conv_temp_k_c(Temperatures array [K])
  conv_factor = -273.15                                                                                                 # Conversion factor from [K] to [°C]
  temp_array = temp_array+conv_factor                                                                                   # Temperature-vals conversion
  return temp_array                                                                                                     # Return temperatures array [°C]

# Function definition to calculate cylindrical lateral surface - EQUATIONS 5.142a and 5.142b 'RelazioniScambiatore'
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

# Function definition to calculate thermal power exchange - EQUATIONS 5.140 and 5.141 'RelazioniScambiatore'
def therm_pow_kw(mass_flow_rate_kg_s, cp_kj_kg_c, delta_temp):                                                          # therm_pow_kw(Mass flow rate [kg/s], Specific heat at const-pressure [kJ/(kg*K)], Delta temperature [°C] or [K])
  return mass_flow_rate_kg_s*cp_kj_kg_c*delta_temp                                                                      # Return thermal power [kJ/s]=[kW]  -->  Q-L'=delta(H) with L'=0  =>  Q=delta(H)=m.*cp*delta(T)

# Function definition to calculate global heat transfer coefficient (global HTC)
def glob_htc_coeff_kw_m2_k(thermal_pow_kw, surf_m2, lmtd):                                                              # glob_htc_coeff_kw_m2_k(Thermal power [kW], Surface [m^2], Log-mean temperature difference [°C] or [K])
  return abs(thermal_pow_kw)/(surf_m2*lmtd)                                                                             # Return positive global heat transfer coefficient (global HTC) [kW/(m^2*K)]

# Function to calculate C-point-min (min[mass-flow-rate*Cp]) and C-point-max (max[mass-flow-rate*Cp])
def cptmin_cptmax(mass_flow_rate_cold, cp_cold, mass_flow_rate_hot, cp_hot):                                            # cptmin_cptmax(Mass flow rate cold fluid [kg/s], Specific heat at const-pressure cold fluid [kJ/(kg*k)], Mass flow rate hot fluid [kg/s], Specific heat at const-pressure hot fluid [kJ/(kg*k)])
  cpt_cold = mass_flow_rate_cold*cp_cold                                                                                # Calculate cold fluid C-point [kJ/(K*s)] --> [kg/s]*[kJ/(kg*K)]=[kJ/(K*s)]
  cpt_hot = mass_flow_rate_hot*cp_hot                                                                                   # Calculate hot fluid C-point [kJ/(K*s)] --> [kg/s]*[kJ/(kg*K)]=[kJ/(K*s)]
  if (cpt_cold < cpt_hot):                                                                                              # If cold fluid C-point is less than hot fluid C-point
    return cpt_cold, cpt_hot                                                                                            # Return (cold C-point, hot C-point) --> (Cptmin, Cptmax)
  else:                                                                                                                 # Else if cold fluid C-point is greater than hot fluid C-point
    return cpt_hot, cpt_cold                                                                                            # Return (hot C-point, cold C-point) --> (Cptmin, Cptmax)

# Function definition to calculate the number of transfer units (NTU)
def ntu(glob_htc_coeff, surf_m2, cpt_min):                                                                              # ntu(Global heat transfer coefficient [kW/(m^2*K)], Surface [m^2], C-point-min [kJ/(K*s)])
  return (glob_htc_coeff*surf_m2)/cpt_min                                                                               # Return calculated NTU

# Function definition to calculate effectiveness
def effectiveness(meas_typ, ntu, cpts_ratio):                                                                           # effectiveness(Measure type: Cocurrent/countercurrent/undefined, Number of tranfer units, C-points-ratio --> Cptmin/Cptmax)
  if (meas_typ == da.Meas_typ.ccurr):                                                                                   # In case of cocurrent measure type
    return (1-mt.exp(-ntu*(1+cpts_ratio)))/(1+cpts_ratio)                                                               # Return calculated effectiveness (epsilon)
  elif (meas_typ == da.Meas_typ.cntcurr):                                                                               # Else in case of countercurrent measure type
    return (1-mt.exp(-ntu*(1-cpts_ratio)))/(1-cpts_ratio*mt.exp(-ntu*(1-cpts_ratio)))                                   # Return calculated effectiveness (epsilon)
  else:                                                                                                                 # Else in case of undefined measure type
    return -1                                                                                                           # Return calculated effectiveness (epsilon)

# Function definition to print and save measures calcs results
def print_save_measures_calcs_res(measures, dbg_flg, out_typ):                                                          # print_save_measures_calcs_res(Measures list, Debug flag, Output save type)
  if (dbg_flg):                                                                                                         # If dbg flg is ena
    dbg_str = str()                                                                                                     # Dbg fbk
    for meas in measures:                                                                                               # Measures scrollin' cycle
      dbg_str += meas.get_info()                                                                                        # Get measures calcs results by callin' get-info method of measure class
    print(dbg_str)                                                                                                      # Print dbg fbk
    out.save_output(out_typ, dbg_str)                                                                                   # Save dbg output
    return                                                                                                              # Return nothing

# Function definition to calculate dynamic viscosity from kinematic viscosity and density vs temp
def dyn_vis(f_ni, f_rho, temp):                                                                                         # dyn_vis(Ni vs temp: kinematic viscosity [m^2/s], Rho vs temp: density [kg/m^3], Temperature [°C])
  return f_ni(temp)*f_rho(temp)                                                                                         # Return dynamic viscosity [kg/(m*s)] --> [m^2/s]*][kg/m^3]=[kg/(m*s)]=[Pa*S]=[Pl]

# Function definition to calculate Reynolds number for fluid in cylindrical pipe
# tube side - EQUATIONS 5.144 and 5.144a 'RelazioniScambiatore'
def re_cyl_pipe(pipe_mass_flow_rate, pipe_int_diam, fluid_dyn_vis):                                                     # re_cyl_pipe(Pipe mass flow rate [kg/s], Pipe internal diameter [m], Mu: fluid dynamic viscosity [kg/(m*S)])
  return (4*pipe_mass_flow_rate)/(mt.pi*pipe_int_diam*fluid_dyn_vis)                                                    # Return calculated Reynolds number for fluid in cylindrical pipe [adimensional] --> [kg/s]/([m]*[kg/(m*s)])=[adimensional]

# Function definition to calculate Reynolds number for fluid in glass pipe shell side
def re_shell():                                                                                                         # re_shell()
  return                                                                                                                # Return ---

# Function definition to calculate Reynolds number
def re(side, pipe_mass_flow_rate, pipe_int_diam, fluid_dyn_vis):                                                        # ---re(Heat-exchanger side, Pipe mass flow rate [kg/s], Pipe internal diameter [m], Mu: fluid dynamic viscosity [kg/(m*S)])
  if (side == He_side.tube):                                                                                            # In case of tube side
    return re_cyl_pipe(pipe_mass_flow_rate, pipe_int_diam, fluid_dyn_vis)                                               # Return tube side calculated Reynolds number
  else:                                                                                                                 # Else in case of shell side
    return re_shell()                                                                                                   # Return shell side calculated Reynolds number

# Function definition to calculate Nusselt number in laminar flow
# tube side - EQUATION 5.145 'RelazioniScambiatore'
def nu_lam_flow_cyl_pipe():                                                                                             # nu_lam_flow_cyl_pipe()
  return                                                                                                                # Return ---

# Function definition to calculate friction factor in turbulent flow
# tube side - EQUATION 5.147 'RelazioniScambiatore'
def frict_fact_turb_flow_cyl_pipe():                                                                                    # frict_fact_turb_flow_cyl_pipe()
  return                                                                                                                # Return ---

# Function definition to calculate Nusselt number in turbulent flow
# tube side - EQUATION 5.146 'RelazioniScambiatore'
def nu_turb_flow_cyl_pipe():                                                                                            # nu_turb_flow_cyl_pipe()
  f = frict_fact_turb_flow_cyl_pipe()                                                                                   # -
  return                                                                                                                # Return ---

# Function definition to calculate Nusselt number shell side - EQUATION 5.148
# (last-one) 'RelazioniScambiatore'
def nu_shell():                                                                                                         # nu_shell()
  return                                                                                                                # Return ---

# Function definition to calculate Nusselt number
def nu(side, re, pr):                                                                                                   # ---nu(Heat-exchanger side, Reynolds number, Prandtl number)
  if (side == He_side.tube):                                                                                            # In case of tube side
    if (re < re_trans):                                                                                                 # And in case of laminar fluid flow
      return nu_lam_flow_cyl_pipe()                                                                                     # Return tube side laminar fluid flow calculated Nusselt number
    else:                                                                                                               # Else in case of tube side and turbulent fluid flow
      return nu_turb_flow_cyl_pipe()                                                                                    # Return tube side turbulent fluid flow calculated Nusselt number
  else:                                                                                                                 # Else in case of shell side
    return nu_shell()                                                                                                   # Return shell side calculated Nusselt number

# Function definition to calculate heat transfer coefficient (h) from Nusselt number Nu=(h*eq_hyd_diam)/therm_cond
# -->  h=(Nu*therm_cond)/eq_hyd_diam [W/(m^2*K)] --> ([adim]*[W/(m*K)])/[m]=[W/(m^2*K)]
def h_from_nu_w_m2_k(nu, eq_hyd_diam, therm_cond):                                                                      # h_from_nu_w_m2_k(Nusselt number, Equivalent hydraulic diameter [m], Thermal conductivity lambda [W/(m*K)])
  return (nu*therm_cond)/eq_hyd_diam                                                                                    # Return heat transfer coefficient (h) [W/(m^2*K)] calculated using Nusselt number

# Function definition to calculate overall heat transfer coefficient - EQUATION 5.143
def overall_htc():                                                                                                      # overall_htc()
  return                                                                                                                # Return ---
