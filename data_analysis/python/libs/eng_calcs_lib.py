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
_eff_len_m = 0.68                                                                                                       # Heat-exchamher effective length [m]
_glass_pipe_id_m = 5*1e-2                                                                                               # Heat-exchamher glass pipe internal-diameter [m]
_supp_glass_pipe_thick_m = 1*1e-3                                                                                       # Supposed heat-exchanger glass pipe thickness [m]
_glass_pipe_diaph_num  = 13                                                                                             # Number of heat-exchanger diaphragms in glass pipe
_glass_pipe_diaph_pace_m = 5*1e-2                                                                                       # Heat-exchanger glass pipe diaphragms pace [m]
_steel_pipes_id_m = 8*1e-3                                                                                              # Heat-exchamher steel pipes internal-diameter [m]
_steel_pipes_ed_m = 1*1e-2                                                                                              # Heat-exchamher steel pipes external-diameter [m]
_steel_pipes_num = 5                                                                                                    # Number of heat-exchamher steel pipes
_steel_pipes_pace_m = 15*1e-3                                                                                           # Heat-exchanger steel pipes pace [m]
_env_temp = 21.3                                                                                                        # Environment temperature [°C]
# Rounding vars
rnd_dec_places = 4                                                                                                      # Number of decimal places in rounding operations to avoid floating-point errors
# Engineering calcs vars
re_trans = 2.3*1e3                                                                                                      # Reynolds number transition value between fluid laminar and turbulent flow [adimensional]
zero_k = -273.15                                                                                                        # Zero Kelvin absolute temperature value [°C]

########
# DEFS #
########

# Heat-exchanger side enum definition
class He_side(en.Enum):                                                                                                 # Heat-exchanger side enum class
  tube = 1                                                                                                              # Tube side (Hot fluid in cylindrical steel pipes)
  shell = 2                                                                                                             # Shell side (Cold fluid in glass pipe)
# Heat-exchanger section enum definition
class He_sect(en.Enum):                                                                                                 # Heat-exchanger section enum class
  one_pipe = 1                                                                                                          # 1-pipe heat-exchanger section (Hot fluid in cylindrical steel pipes)
  three_pipes = 2                                                                                                       # 3-pipes heat-exchanger section (Cold fluid in glass pipe)
# Heat-exchanger class
class He:                                                                                                               # Heat-exchanger class (attributes, constructor, methods)
  eff_len_m = 0.0                                                                                                       # Heat-exchanger effective length [m]
  glass_pipe_id_m = 0.0                                                                                                 # Heat-exchanger glass pipe internal-diameter [m]
  glass_pipe_ibs_m2 = 0.0                                                                                               # Heat-exchanger glass pipe internal-basal surface [m^2]
  glass_pipe_ils_m2 = 0.0                                                                                               # Heat-exchanger glass pipe internal-lateral surface [m^2]
  supp_glass_pipe_thick_m = 0.0                                                                                         # Supposed heat-exchanger glass pipe thickness [m]
  supp_glass_pipe_ed_m = 0.0                                                                                            # Supposed heat-exchanger glass pipe external-diameter [m]
  supp_glass_pipe_els_m2 = 0.0                                                                                          # Supposed heat-exchanger glass pipe external-lateral surface [m^2]
  glass_pipe_diaph_num = 0                                                                                              # Number of heat-exchanger diaphragms in glass pipe
  glass_pipe_diaph_pace_m = 0.0                                                                                         # Heat-exchanger glass pipe diaphragms pace [m]
  steel_pipes_id_m = 0.0                                                                                                # Heat-exchanger steel pipes internal-diameter [m]
  steel_pipes_ed_m = 0.0                                                                                                # Heat-exchanger steel pipes external-diameter [m]
  steel_pipes_num = 0                                                                                                   # Number of heat-exchanger steel pipes
  steel_pipes_pace_m = 0.0                                                                                              # Heat-exchanger steel pipes pace [m]
  steel_pipes_thick_m = 0.0                                                                                             # Heat-exchanger steel pipes thickness [m]
  steel_pipes_interspce = 0.0                                                                                           # Heat-exchanger interspace between steel pipes [m]
  steel_pipes_ibs_m2 = 0.0                                                                                              # Heat-exchanger steel pipes internal-basal surface [m^2]
  steel_pipes_ils_m2 = 0.0                                                                                              # Heat-exchanger steel pipes internal-lateral surface [m^2]
  steel_pipes_els_m2 = 0.0                                                                                              # Heat-exchanger steel pipes external-lateral surface [m^2]
  glass_pipe_1pipe_sect_fs_m2 = 0.0                                                                                     # Heat-exchanger glass pipe 1-pipe section flow surface [m^2]
  glass_pipe_1pipe_sect_whp_m = 0.0                                                                                     # Heat-exchanger glass pipe 1-pipe section wetted hydraulic perimeter [m]
  glass_pipe_1pipe_sect_ehd_m = 0.0                                                                                     # Heat-exchanger glass pipe 1-pipe section equivalent hydraulic diameter [m]
  glass_pipe_1pipe_sect_wtp_m = 0.0                                                                                     # Heat-exchanger glass pipe 1-pipe section wetted thermal perimeter [m]
  glass_pipe_1pipe_sect_etd_m = 0.0                                                                                     # Heat-exchanger glass pipe 1-pipe section equivalent thermal diameter [m]
  glass_pipe_3pipes_sect_fs_m2 = 0.0                                                                                    # Heat-exchanger glass pipe 3-pipes section flow surface [m^2]
  glass_pipe_3pipes_sect_whp_m = 0.0                                                                                    # Heat-exchanger glass pipe 3-pipes section wetted hydraulic perimeter [m]
  glass_pipe_3pipes_sect_ehd_m = 0.0                                                                                    # Heat-exchanger glass pipe 3-pipes section equivalent hydraulic diameter [m]
  glass_pipe_3pipes_sect_wtp_m = 0.0                                                                                    # Heat-exchanger glass pipe 3-pipes section wetted thermal perimeter [m]
  glass_pipe_3pipes_sect_etd_m = 0.0                                                                                    # Heat-exchanger glass pipe 3-pipes section equivalent thermal diameter [m]
  env_temp = 0.0                                                                                                        # Environment temperature [°C]
  def __init__(self):                                                                                                   # Constructor
    self.eff_len_m = _eff_len_m                                                                                         # Heat-exchamher effective length [m] init
    self.glass_pipe_id_m = _glass_pipe_id_m                                                                             # Heat-exchamher glass pipe internal-diameter [m] init
    self.glass_pipe_ibs_m2 = cyl_bas_surf(self.glass_pipe_id_m)                                                         # Heat-exchanger glass pipe internal-basal surface [m^2] init
    self.glass_pipe_ils_m2 = cyl_lat_surf(1, self.glass_pipe_id_m, self.eff_len_m)                                      # Heat-exchanger glass pipe internal-lateral surface [m^2] init
    self.supp_glass_pipe_thick_m = _supp_glass_pipe_thick_m                                                             # Supposed heat-exchanger glass pipe thickness [m] init
    self.supp_glass_pipe_ed_m = round(self.glass_pipe_id_m+2*self.supp_glass_pipe_thick_m, rnd_dec_places)              # Supposed heat-exchanger glass pipe external-diameter [m] init and round value to avoid floating-point errors
    self.supp_glass_pipe_els_m2 = cyl_lat_surf(1, self.supp_glass_pipe_ed_m, self.eff_len_m)                            # Supposed heat-exchanger glass pipe external-lateral surface [m^2] init
    self.glass_pipe_diaph_num = _glass_pipe_diaph_num                                                                   # Number of heat-exchanger diaphragms in glass pipe init
    self.glass_pipe_diaph_pace_m = _glass_pipe_diaph_pace_m                                                             # Heat-exchanger glass pipe diaphragms pace [m] init
    self.steel_pipes_id_m = _steel_pipes_id_m                                                                           # Heat-exchamher steel pipes internal-diameter [m] init
    self.steel_pipes_ed_m = _steel_pipes_ed_m                                                                           # Heat-exchamher steel pipes external-diameter [m] init
    self.steel_pipes_num = _steel_pipes_num                                                                             # Number of heat-exchamher steel pipes init
    self.steel_pipes_pace_m = _steel_pipes_pace_m                                                                       # Heat-exchanger steel pipes pace [m] init
    self.steel_pipes_thick_m = round((self.steel_pipes_ed_m-self.steel_pipes_id_m)/2, rnd_dec_places)                   # Heat-exchamher steel pipes thickness [m] init and round value to avoid floating-point errors
    self.steel_pipes_interspce = round(self.steel_pipes_pace_m-self.steel_pipes_ed_m, rnd_dec_places)                   # Heat-exchanger interspace between steel pipes [m] init and value to avoid floating-point errors
    self.steel_pipes_ibs_m2 = cyl_bas_surf(self.steel_pipes_id_m)                                                       # Heat-exchanger steel pipes internal-basal surface [m^2] init
    self.steel_pipes_ils_m2 = cyl_lat_surf(self.steel_pipes_num, self.steel_pipes_id_m, self.eff_len_m)                 # Heat-exchanger steel pipes internal-lateral surface [m^2] init
    self.steel_pipes_els_m2 = cyl_lat_surf(self.steel_pipes_num, self.steel_pipes_ed_m, self.eff_len_m)                 # Heat-exchanger steel pipes external-lateral surface [m^2] init
    self.glass_pipe_1pipe_sect_fs_m2 = 2*1.5*self.steel_pipes_interspce*self.glass_pipe_diaph_pace_m                    # Heat-exchanger glass pipe 1-pipe section flow surface [m^2] init
    self.glass_pipe_1pipe_sect_whp_m = 4*(1.5*(self.steel_pipes_interspce)+self.glass_pipe_diaph_pace_m)                # Heat-exchanger glass pipe 1-pipe section wetted hydraulic perimeter [m] init
    self.glass_pipe_1pipe_sect_ehd_m = (4*self.glass_pipe_1pipe_sect_fs_m2)/self.glass_pipe_1pipe_sect_whp_m            # Heat-exchanger glass pipe 1-pipe section equivalent hydraulic diameter [m] init
    self.glass_pipe_1pipe_sect_wtp_m = 2*self.glass_pipe_diaph_pace_m                                                   # Heat-exchanger glass pipe 1-pipe section wetted thermal perimeter [m] init
    self.glass_pipe_1pipe_sect_etd_m = (4*self.glass_pipe_1pipe_sect_fs_m2)/self.glass_pipe_1pipe_sect_wtp_m            # Heat-exchanger glass pipe 1-pipe section equivalent thermal diameter [m] init
    self.glass_pipe_3pipes_sect_fs_m2 = 4*self.steel_pipes_interspce*self.glass_pipe_diaph_pace_m                       # Heat-exchanger glass pipe 3-pipes section flow surface [m^2] init
    self.glass_pipe_3pipes_sect_whp_m = 8*(self.steel_pipes_interspce+self.glass_pipe_diaph_pace_m)                     # Heat-exchanger glass pipe 3-pipes section wetted hydraulic perimeter [m] init
    self.glass_pipe_3pipes_sect_ehd_m = (4*self.glass_pipe_3pipes_sect_fs_m2)/self.glass_pipe_3pipes_sect_whp_m         # Heat-exchanger glass pipe 3-pipes section equivalent hydraulic diameter [m] init
    self.glass_pipe_3pipes_sect_wtp_m = 6*self.glass_pipe_diaph_pace_m                                                  # Heat-exchanger glass pipe 3-pipes section wetted thermal perimeter [m] init
    self.glass_pipe_3pipes_sect_etd_m = (4*self.glass_pipe_3pipes_sect_fs_m2)/self.glass_pipe_3pipes_sect_wtp_m         # Heat-exchanger glass pipe 3-pipes section equivalent thermal diameter [m] init
    self.env_temp = _env_temp                                                                                           # Environment temperature [°C] init
    return                                                                                                              # Return nothing
  def print_info_save_out(self, dbg_flg):                                                                               # Info printing and output saving method with debug flag
    if (dbg_flg):                                                                                                       # If dbg flg is ena
      dbg_str = ("\n--> Heat-exchanger effective length: "+str(self.eff_len_m)+" [m]"\
      +"\n--> Heat-exchanger glass pipe internal-diameter: "+str(self.glass_pipe_id_m)+" [m]"\
      +"\n--> Heat-exchanger glass pipe internal-basal surface: "+str(self.glass_pipe_ibs_m2)+" [m^2]"
      +"\n--> Heat-exchanger glass pipe internal-lateral surface: "+str(self.glass_pipe_ils_m2)+" [m^2]"
      +"\n--> Supposed heat-exchanger glass pipe thickness: "+str(self.supp_glass_pipe_thick_m)+" [m]"\
      +"\n--> Supposed heat-exchanger glass pipe external-diameter: "+str(self.supp_glass_pipe_ed_m)+" [m]"\
      +"\n--> Supposed heat-exchanger glass pipe external-lateral surface: "\
        +str(self.supp_glass_pipe_els_m2)+" [m^2]"
      +"\n--> Number of heat-exchanger diaphragms in glass pipe: "+str(self.glass_pipe_diaph_num)\
      +"\n--> Heat-exchanger glass pipe diaphragms pace: "+str(self.glass_pipe_diaph_pace_m)+" [m]"\
      +"\n--> Heat-exchanger steel pipes internal-diameter: "+str(self.steel_pipes_id_m)+" [m]"\
      +"\n--> Heat-exchanger steel pipes external-diameter: "+str(self.steel_pipes_ed_m)+" [m]"\
      +"\n--> Number of heat-exchanger steel pipes: "+str(self.steel_pipes_num)\
      +"\n--> Heat-exchanger steel pipes pace: "+str(self.steel_pipes_pace_m)+" [m]"\
      +"\n--> Heat-exchanger steel pipes thickness: "+str(self.steel_pipes_thick_m)+" [m]"\
      +"\n--> Heat-exchanger interspace between steel pipes: "+str(self.steel_pipes_interspce)+" [m]"\
      +"\n--> Heat-exchanger steel pipes internal-basal surface: "+str(self.steel_pipes_ibs_m2)+" [m^2]"\
      +"\n--> Heat-exchanger steel pipes internal-lateral surface: "+str(self.steel_pipes_ils_m2)+" [m^2]"\
      +"\n--> Heat-exchanger steel pipes external-lateral surface: "+str(self.steel_pipes_els_m2)+" [m^2]"\
      +"\n--> Heat-exchanger glass pipe 1-pipe section flow surface: "\
        +str(self.glass_pipe_1pipe_sect_fs_m2)+" [m^2]"\
      +"\n--> Heat-exchanger glass pipe 1-pipe section wetted hydraulic perimeter: "\
        +str(self.glass_pipe_1pipe_sect_whp_m)+" [m]"
      +"\n--> Heat-exchanger glass pipe 1-pipe section equivalent hydraulic diameter: "\
        +str(self.glass_pipe_1pipe_sect_ehd_m)+" [m]"\
      +"\n--> Heat-exchanger glass pipe 1-pipe section wetted thermal perimeter: "\
        +str(self.glass_pipe_1pipe_sect_wtp_m)+" [m]"
      +"\n--> Heat-exchanger glass pipe 1-pipe section equivalent thermal diameter: "\
        +str(self.glass_pipe_1pipe_sect_etd_m)+" [m]"\
      +"\n--> Heat-exchanger glass pipe 3-pipes section flow surface: "\
        +str(self.glass_pipe_3pipes_sect_fs_m2)+" [m^2]"\
      +"\n--> Heat-exchanger glass pipe 3-pipes section wetted hydraulic perimeter: "\
        +str(self.glass_pipe_3pipes_sect_whp_m)+" [m]"
      +"\n--> Heat-exchanger glass pipe 3-pipes section equivalent hydraulic diameter: "\
        +str(self.glass_pipe_3pipes_sect_ehd_m)+" [m]"\
      +"\n--> Heat-exchanger glass pipe 3-pipes section wetted thermal perimeter: "\
        +str(self.glass_pipe_3pipes_sect_wtp_m)+" [m]"
      +"\n--> Heat-exchanger glass pipe 3-pipes section equivalent thermal diameter: "\
        +str(self.glass_pipe_3pipes_sect_etd_m)+" [m]"\
      +"\n--> Environment temperature: "+str(self.env_temp)+" [°C]\n")                                                  # Dbg fbk
      print(dbg_str)                                                                                                    # Print dbg fbk
      out.save_output(out.Output_typ.he, dbg_str)                                                                       # Save dbg output
    return                                                                                                              # Return nothing

##########
# FUNCTS #
##########

# Function definition to calculate average between two values
def avg(val1, val2):                                                                                                    # avg(Value1, Value2)
  return (val1+val2)/2                                                                                                  # Return average value btwn Value1 and Value2

# Function definition to calculate percentage between partial and total values specifying rounding
def perc(partial, total, dec_cyph):                                                                                     # perc(Partial value, Total value, Number of decimal cyphers)
  if (partial != None and total != 0 and dec_cyph >= 0):                                                                # Check in-vals consistency, if ok
    return round(float(100.0)*float(abs(partial))/float(abs(total)), dec_cyph)                                          # Return calculated percentage between partial and total values with specified rounding
  else:                                                                                                                 # Else if in-vals consistency ain't ok
    return -1                                                                                                           # Retur err val

# Function definition to convert temperatures from [K] to [°C]
def conv_temp_k_c(temp_array):                                                                                          # conv_temp_k_c(Temperatures array [K])
  conv_factor = zero_k                                                                                                  # Conversion factor from [K] to [°C]
  return temp_array+conv_factor                                                                                         # Return converted temperatures array [°C]

# Function definition to convert thermal conductivities (lamba) from [W/(cm*K)] to [W/(m*K)]
def conv_lambda_w_cm_k_w_m_k(lambda_array):                                                                             # conv_lambda_w_cm_k_w_m_k(Lambda: thermal conductivities array [W/(cm*k)])
  conv_factor = 1*1e-2                                                                                                  # Conversion factor from [W/(cm*K)] to [W/(m*K)]
  return lambda_array/conv_factor                                                                                       # Return converted thermal conductivities array [W/(m*k)]

# Function definition to calculate cylindrical lateral surface - EQUATIONS 5.142a and 5.142b 'RelazioniScambiatore'
def cyl_lat_surf(num_cyl, diam, len):                                                                                   # cyl_lat_surf(Number of cylinders, Cylinder diameter, Cylinder length)
  if (num_cyl > 0 and diam > 0 and len > 0):                                                                            # Check in-vals consistency, if ok
    return num_cyl*mt.pi*diam*len                                                                                       # Return cylindrical lateral surface
  else:                                                                                                                 # Else if in-vals consistency ain't ok
    return -1                                                                                                           # Retur err val

# Function definition to calculate cylindrical basal surface
def cyl_bas_surf(diam):                                                                                                 # cyl_bas_surf(Cylinder diameter)
  if (diam > 0):                                                                                                        # Check in-vals consistency, if ok
    return mt.pi*mt.pow(diam/2, 2)                                                                                      # Return cylindrical basal surface
  else:                                                                                                                 # Else if in-vals consistency ain't ok
    return -1                                                                                                           # Retur err val

# Function definition to convert volume flow rate [l/h] into [m^3/s] and mass flow rate [kg/s]
def vol_flow_rate_mass_flow_rate(vol_flow_rate_l_h, density):                                                           # vol_flow_rate_mass_flow_rate(Volume flow rate value [l/h], density value)
  if (vol_flow_rate_l_h > 0 and density > 0):                                                                           # Check in-vals consistency, if ok
    conv_factor = (1*1e-3)/(3.6*1e3)                                                                                    # Conversion factor from [l/h] to [m^3/s]
    vol_flow_rate_m3_s = vol_flow_rate_l_h*conv_factor                                                                  # Volume flow rate conversion from [l/h] to [m^3/s]
    return vol_flow_rate_m3_s, vol_flow_rate_m3_s*density                                                               # Return volume flow rate [m^3/s] and mass flow rate [kg/s] --> [m^3/s]*[kg/m^3]=[kg/s]
  else:                                                                                                                 # Else if in-vals consistency ain't ok
    return -1, -1                                                                                                       # Return err vals

# Function definition to calculate log-mean temperature difference (LMTD)
def lmtd(delta1, delta2):                                                                                               # lmtd(delta(Temperature1), delta(Temperature2))
  if (delta1/delta2 > 0):                                                                                               # Check in-vals consistency, if ok
    return (delta1-delta2)/mt.log((delta1/delta2), mt.e)                                                                # Return log-mean temperature difference (LMTD)
  else:                                                                                                                 # Else if in-vals consistency ain't ok
    return -1                                                                                                           # Return err val

# Function definition to calculate thermal power exchange - EQUATIONS 5.140 and 5.141 'RelazioniScambiatore'
def therm_pow_kw(mass_flow_rate_kg_s, cp_kj_kg_c, delta_temp):                                                          # therm_pow_kw(Mass flow rate [kg/s], Specific heat at const-pressure [kJ/(kg*K)], Delta temperature [°C] or [K])
  if (mass_flow_rate_kg_s > 0 and cp_kj_kg_c > 0 and delta_temp != 0):                                                  # Check in-vals consistency, if ok
    return mass_flow_rate_kg_s*cp_kj_kg_c*delta_temp                                                                    # Return thermal power [kJ/s]=[kW]  -->  Q-L'=delta(H) with L'=0  =>  Q=delta(H)=m.*cp*delta(T)
  else:                                                                                                                 # Else if in-vals consistency ain't ok
    return -1                                                                                                           # Return err val

# Function definition to calculate global heat transfer coefficient (global HTC) [kW/(m^2*K)]
def glob_htc_kw_m2_k(thermal_pow_kw, surf_m2, lmtd):                                                                    # glob_htc_kw_m2_k(Thermal power [kW], Surface [m^2], Log-mean temperature difference [°C] or [K])
  if (thermal_pow_kw != 0 and surf_m2 > 0 and lmtd != 0):                                                               # Check in-vals consistency, if ok
    return abs(thermal_pow_kw)/(surf_m2*lmtd)                                                                           # Return positive global heat transfer coefficient (global HTC) [kW/(m^2*K)]
  else:                                                                                                                 # Else if in-vals consistency ain't ok
    return -1                                                                                                           # Return err val

# Function to calculate C-point-min (min[mass-flow-rate*Cp]) and C-point-max (max[mass-flow-rate*Cp])
def cptmin_cptmax(mass_flow_rate_cold, cp_cold, mass_flow_rate_hot, cp_hot):                                            # cptmin_cptmax(Mass flow rate cold fluid [kg/s], Specific heat at const-pressure cold fluid [kJ/(kg*k)], Mass flow rate hot fluid [kg/s], Specific heat at const-pressure hot fluid [kJ/(kg*k)])
  cpt_cold = mass_flow_rate_cold*cp_cold                                                                                # Calculate cold fluid C-point [kJ/(K*s)] --> [kg/s]*[kJ/(kg*K)]=[kJ/(K*s)]
  cpt_hot = mass_flow_rate_hot*cp_hot                                                                                   # Calculate hot fluid C-point [kJ/(K*s)] --> [kg/s]*[kJ/(kg*K)]=[kJ/(K*s)]
  if (mass_flow_rate_cold > 0 and cp_cold > 0 and mass_flow_rate_hot > 0 and cp_hot > 0):                               # Check in-vals consistency, if ok
    if (cpt_cold < cpt_hot):                                                                                            # If cold fluid C-point is less than hot fluid C-point
      return cpt_cold, cpt_hot                                                                                          # Return (cold C-point, hot C-point) --> (Cptmin, Cptmax)
    else:                                                                                                               # Else if cold fluid C-point is greater than hot fluid C-point
      return cpt_hot, cpt_cold                                                                                          # Return (hot C-point, cold C-point) --> (Cptmin, Cptmax)
  else:                                                                                                                 # Else if in-vals consistency ain't ok
    return -1                                                                                                           # Return err val

# Function definition to calculate the number of transfer units (NTU)
def ntu(glob_htc_coeff, surf_m2, cpt_min):                                                                              # ntu(Global heat transfer coefficient [kW/(m^2*K)], Surface [m^2], C-point-min [kJ/(K*s)])
  if (glob_htc_coeff > 0 and surf_m2 > 0 and cpt_min > 0):                                                              # Check in-vals consistency, if ok
    return (glob_htc_coeff*surf_m2)/cpt_min                                                                             # Return calculated NTU
  else:                                                                                                                 # Else if in-vals consistency ain't ok
    return -1                                                                                                           # Return err val

# Function definition to calculate effectiveness
def effectiveness(meas_typ, ntu, cpts_ratio):                                                                           # effectiveness(Measure type: Cocurrent/countercurrent/undefined, Number of tranfer units, C-points-ratio --> Cptmin/Cptmax)
  if (meas_typ != None and ntu > 0 and cpts_ratio > 0):                                                                 # Check in-vals consistency, if ok
    if (meas_typ == da.Meas_typ.ccurr):                                                                                 # In case of cocurrent measure type
      epsilon = (1-mt.exp(-ntu*(1+cpts_ratio)))/(1+cpts_ratio)                                                          # Calc effectiveness (epsilon)
    elif (meas_typ == da.Meas_typ.cntcurr):                                                                             # Else in case of countercurrent measure type
      epsilon = (1-mt.exp(-ntu*(1-cpts_ratio)))/(1-cpts_ratio*mt.exp(-ntu*(1-cpts_ratio)))                              # Calc effectiveness (epsilon)
    else:                                                                                                               # Else in case of undefined measure type
      return -1                                                                                                         # Return err val
    if (epsilon < 1 and epsilon > 0):                                                                                   # Check effectiveness val consistency, if ok
      return epsilon                                                                                                    # Return calc effectiveness (epsilon)
    else:                                                                                                               # Else if effectiveness val consistency ain't ok
      return -1                                                                                                         # Return err val
  else:                                                                                                                 # Else if in-vals consistency ain't ok
    return -1                                                                                                           # Return err val

# Function definition to print and save measures calcs results
def print_save_measures_calcs_res(measures, dbg_flg, out_typ):                                                          # print_save_measures_calcs_res(Measures list, Debug flag, Output save type)
  if (dbg_flg):                                                                                                         # If dbg flg is ena
    dbg_str = str()                                                                                                     # Dbg fbk
    for meas in measures:                                                                                               # Measures scrollin' cycle
      dbg_str += meas.get_info()                                                                                        # Get measures calcs results by callin' get-info method of measure class
    print(dbg_str)                                                                                                      # Print dbg fbk
    out.save_output(out_typ, dbg_str)                                                                                   # Save dbg output
    return                                                                                                              # Return nothing

# Function definition to approximate measure heat-exchange surface temperature [°C]
# and define shell inlet/outlet sections fluid temperatures [°C]
def define_temp(meas):                                                                                                  # define_temp(Measure)
  if (meas != None or meas.typ != da.Meas_typ.undef):                                                                   # If measure is defined and in-vals consistency is ok
    approx_surf_temp = avg(meas.avg_cold_fl_temp, meas.avg_hot_fl_temp)                                                 # Def heat-exchange average surface approximated temperature [°C] for steel thermal conductivity calculation
    if (meas.typ == da.Meas_typ.ccurr):                                                                                 # In case of measure taken while heat-exchanger was in cocurrent config
      inlet_sect_temp = meas.t1                                                                                         # Def shell inlet section cold fluid temperature (bottom section)
      outlet_sect_temp = meas.t3                                                                                        # Def shell outlet section cold fluid temperature (top section)
    else:                                                                                                               # Else in case of measure taken while heat-exchanger was in countercurrent config
      inlet_sect_temp = meas.t3                                                                                         # Def shell inlet section cold fluid temperature (bottom section)
      outlet_sect_temp = meas.t1                                                                                        # Def shell outlet section cold fluid temperature (top section)
    return approx_surf_temp, inlet_sect_temp, outlet_sect_temp                                                          # Return approximated-heat-exch-surface-temp, shell-inlet-sect-fluid-temp, shell-outlet-sect-fluid-temp
  else:                                                                                                                 # Else if measure is undefined or in-vals consistency ain't ok
    return -1, -1, -1                                                                                                   # Return err vals

# Function definition to calculate dynamic viscosity from kinematic viscosity and density vs temp
def dyn_vis(f_ni, f_rho, temp):                                                                                         # dyn_vis(Ni vs temp: kinematic viscosity [m^2/s], Rho vs temp: density [kg/m^3], Temperature [°C])
  if (f_ni != None and f_rho != None and temp >= zero_k):                                                               # Check in-vals consistency, if ok
    return f_ni(temp)*f_rho(temp)                                                                                       # Return dynamic viscosity [kg/(m*s)] --> [m^2/s]*][kg/m^3]=[kg/(m*s)]=[Pa*S]=[Pl]
  else:                                                                                                                 # Else if in-vals consistency ain't ok
    return -1                                                                                                           # Return err val

# Function definition to calculate Reynolds number for fluid in cylindrical pipe
# tube side - EQUATIONS 5.144 and 5.144a 'RelazioniScambiatore'
def re_cyl_pipe(pipe_mass_flow_rate, pipe_int_diam, fluid_dyn_vis):                                                     # re_cyl_pipe(Pipe mass flow rate [kg/s], Pipe internal diameter [m], Mu: fluid dynamic viscosity [kg/(m*s)])
  if (pipe_mass_flow_rate > 0 and pipe_int_diam > 0 and fluid_dyn_vis > 0):                                             # Check in-vals consistency, if ok
    return (4*pipe_mass_flow_rate)/(mt.pi*pipe_int_diam*fluid_dyn_vis)                                                  # Return calculated Reynolds number for fluid in cylindrical pipe [adimensional] --> [kg/s]/([m]*[kg/(m*s)])=[adimensional]
  else:                                                                                                                 # Else if in-vals consistency ain't ok
    return -1                                                                                                           # Return err val

# Function definition to calculate Reynolds number from EQUATIONS 5.144 and 5.144a 'RelazioniScambiatore' using
# volume flow rate to calculate fluid velocity, pipe equivalent diameter and fluid kinematic viscosity (alternative)
def re_aternative(pipe_vol_flow_rate, fluid_in_surf_m2, pipe_equiv_diam, fluid_kin_vis):                                # re_aternative(Pipe volume flow rate [m^3/s], fluid inlet surface [m^2], Pipe equivalent diameter [m], kinematic viscosity [m^2/s])
  if (pipe_vol_flow_rate > 0 and pipe_equiv_diam > 0 and fluid_kin_vis > 0):                                            # Check in-vals consistency, if ok
    return ((pipe_vol_flow_rate/fluid_in_surf_m2)*pipe_equiv_diam)/fluid_kin_vis                                        # Return calculated Reynolds number [adimensional] --> ([m^3/s]*[m])/([m^2]*[m^2/s])=[adimensional]
  else:                                                                                                                 # Else if in-vals consistency ain't ok
    return -1                                                                                                           # Return err val

# ---Function definition to calculate Reynolds number for fluid inside
# external pipe shell side section
def re_shell_sect(sect, pipe_mass_flow_rate, fluid_dyn_vis, he):                                                        # re_shell_sect(Heat-exchanger section, Pipe mass flow rate [kg/s], Mu: fluid dynamic viscosity [kg/(m*s)], Heat-exchanger)
  if ((sect == He_sect.one_pipe or sect == He_sect.three_pipes) and\
    pipe_mass_flow_rate > 0 and fluid_dyn_vis > 0 and he != None):                                                      # Check in-vals consistency, if ok
    if (sect == He_sect.one_pipe):                                                                                      # In case of 1-pipe section
      shell_sect_eq_hyd_diam = he.glass_pipe_1pipe_sect_ehd_m                                                           # Def shell section equivalent diameter
      shell_sect_flow_surf = he.glass_pipe_1pipe_sect_fs_m2                                                             # Def shell section flow surface
    else:                                                                                                               # Else in case of 3-pipes section
      shell_sect_eq_hyd_diam = he.glass_pipe_3pipes_sect_ehd_m                                                          # Def shell section equivalent diameter
      shell_sect_flow_surf = he.glass_pipe_3pipes_sect_fs_m2                                                            # Def shell section flow surface
    return (shell_sect_eq_hyd_diam*pipe_mass_flow_rate)/(shell_sect_flow_surf*fluid_dyn_vis)                            # Return calculated Reynolds number for fluid in external pipe shell side section
  else:                                                                                                                 # Else if in-vals consistency ain't ok
    return -1                                                                                                           # Return err val

# Function definition to calculate Reynolds number
def re(side, sect, pipe_mass_flow_rate, pipe_int_diam, fluid_dyn_vis, he):                                              # re(Heat-exchanger side, Heat-exchanger section, Pipe mass flow rate [kg/s], Pipe internal diameter [m], Mu: fluid dynamic viscosity [kg/(m*s)], Heat-exchanger)
  if (side == He_side.tube):                                                                                            # In case of tube side
    return re_cyl_pipe(pipe_mass_flow_rate, pipe_int_diam, fluid_dyn_vis)                                               # Return tube side calculated Reynolds number
  elif (side == He_side.shell):                                                                                         # Else in case of shell side
    return re_shell_sect(sect, pipe_mass_flow_rate, fluid_dyn_vis, he)                                                  # Return shell side section calculated Reynolds number
  else:                                                                                                                 # Else in case of unknown side enum val
    return -1                                                                                                           # Return err val

# Function definition to calculate Nusselt number in laminar flow
# tube side - EQUATION 5.145 'RelazioniScambiatore'
def nu_lam_flow_cyl_pipe(re, pr, fluid_dyn_vis, fluid_dyn_vis_s, pipe_int_diam, pipe_len_m):                            # nu_lam_flow_cyl_pipe(Reynolds number, Prandtl number, Mu: fluid dynamic viscosity [kg/(m*s)], Mu-s: fluid dynamic viscosity at surface temperature [kg/(m*s)], Pipe internal diameter [m], Pipe length [m])
  if (re > 0 and re < re_trans and pr > 0.48 and pr < 1.67*1e4 and\
      (fluid_dyn_vis/fluid_dyn_vis_s) > 4.4*1e-3 and (fluid_dyn_vis/fluid_dyn_vis_s) < 9.75 and pipe_int_diam > 0):     # If Reynolds number, Prandtl number, Fluid-dyn-viscosity/Fluid-dyn-viscosity-surf-temp ratio are in formula range and in-vals consistency is ok
    nu = 1.86*mt.pow((pipe_int_diam*re*pr)/pipe_len_m, 1/3)*mt.pow(fluid_dyn_vis/fluid_dyn_vis_s, 0.14)                 # Return calculated Nusselt number for fluid in laminar flow inside cylindrical pipe
    if (nu < 3.66):                                                                                                     # If calculated Nusselt number is less than the formula law change point
      nu = 3.66                                                                                                         # Apply formula law change
    return nu                                                                                                           # Return Nusselt number for fluid in laminar flow inside cylindrical pipe
  else:                                                                                                                 # Else if Reynolds number, Prandtl number, Fluid-dyn-viscosity/Fluid-dyn-viscosity-surf-temp ratio ain't in formula range or in-vals consistency not ok
    return -1                                                                                                           # Return err val

# Function definition to calculate friction factor in turbulent flow
# tube side - EQUATION 5.147 'RelazioniScambiatore'
def frict_fact_turb_flow_cyl_pipe(re):                                                                                  # frict_fact_turb_flow_cyl_pipe(Reynolds number)
  if (re > 0):                                                                                                          # Check in-vals consistency, if ok
    return mt.pow((1.58*mt.log(re, mt.e)-3.28), -2)                                                                     # Return calculated friction factor for fluid in turbulent flow inside cylindrical pipe
  else:                                                                                                                 # Else if in-vals consistency ain't ok
    return -1                                                                                                           # Return err val

# Function definition to calculate Nusselt number in turbulent flow
# tube side - EQUATION 5.146 'RelazioniScambiatore'
def nu_turb_flow_cyl_pipe(re, pr):                                                                                      # nu_turb_flow_cyl_pipe(Reynolds number, Prandtl number)
  f = frict_fact_turb_flow_cyl_pipe(re)                                                                                 # Friction factor for fluid in turbulent flow inside cylindrical pipe
  if (re > 300 and re < 5*1e6 and pr >= 0.5 and pr <= 2*1e3):                                                           # If Reynolds and Prandtl numbers are in formula range
    return ((f/2)*(re-1000)*pr)/(1+12.7*mt.pow(f/2, 1/2)*(mt.pow(pr, 2/3)-1))                                           # Return calculated Nusselt number for fluid in turbulent flow inside cylindrical pipe
  else:                                                                                                                 # Else if Reynolds and Prandtl numbers ain't in formula range
    return -1                                                                                                           # Return err val

# Function definition to calculate Nusselt number shell side - EQUATION 5.148
# (last-one) 'RelazioniScambiatore', Removed formula range: re > 2*1e3, turbulent flow for diaphragms  
def nu_shell(re, pr, fluid_dyn_vis, fluid_dyn_vis_s):                                                                   # nu_shell(Reynolds number, Nusselt number, Mu: fluid dynamic viscosity [kg/(m*s)], Mu-s: fluid dynamic viscosity at surface temperature [kg/(m*s)])
  if (re < 1*1e6 and pr > 0 and fluid_dyn_vis > 0 and fluid_dyn_vis_s > 0):                                             # If Reynolds number is in formula range and in-vals consistency is ok
    return 0.36*mt.pow(re, 0.55)*mt.pow(pr, 1/3)*mt.pow((fluid_dyn_vis/fluid_dyn_vis_s), 0.14)                          # Return calculated Nusselt number for fluid in heat-exchanger shell
  else:                                                                                                                 # Else if Reynolds number ain't in formula range or in-vals consistency not ok
    return -1                                                                                                           # Return err val

# Function definition to calculate Nusselt number
def nu(side, re, pr, fluid_dyn_vis, fluid_dyn_vis_s, pipe_int_diam, pipe_len_m):                                        # nu(Heat-exchanger side, Reynolds number, Prandtl number, Mu: fluid dynamic viscosity [kg/(m*s)], Mu-s: fluid dynamic viscosity at surface temperature [kg/(m*s)], Pipe internal diameter [m], Pipe length [m])
  if (side == He_side.tube):                                                                                            # In case of tube side
    if (re < re_trans):                                                                                                 # And in case of laminar fluid flow
      return nu_lam_flow_cyl_pipe(re, pr, fluid_dyn_vis, fluid_dyn_vis_s, pipe_int_diam, pipe_len_m)                    # Return tube side laminar fluid flow calculated Nusselt number
    else:                                                                                                               # Else in case of tube side and turbulent fluid flow
      return nu_turb_flow_cyl_pipe(re, pr)                                                                              # Return tube side turbulent fluid flow calculated Nusselt number
  elif (side == He_side.shell):                                                                                         # Else in case of shell side
    return nu_shell(re, pr, fluid_dyn_vis, fluid_dyn_vis_s)                                                             # Return shell side calculated Nusselt number
  else:                                                                                                                 # Else in case of unknown side enum val
    return -1                                                                                                           # Return err val

# ---Function definition to calculate heat transfer coefficient (h) from Nusselt number Nu=(h*eq_hyd_diam)/therm_cond
# -->  h=(Nu*therm_cond)/eq_hyd_diam [W/(m^2*K)] --> ([adim]*[W/(m*K)])/[m]=[W/(m^2*K)]
def h_from_nu_w_m2_k(nu, eq_therm_diam, therm_cond):                                                                    # h_from_nu_w_m2_k(Nusselt number, Equivalent thermal diameter [m], Thermal conductivity lambda [W/(m*K)])
  if (nu > 0 and eq_therm_diam > 0 and therm_cond > 0):                                                                 # Check in-vals consistency, if ok
    return (nu*therm_cond)/eq_therm_diam                                                                                # Return heat transfer coefficient (h) [W/(m^2*K)] calculated using Nusselt number
  else:                                                                                                                 # Else if in-vals consistency ain't ok
    return -1                                                                                                           # Return err val

# Function definition to calculate absolute convective resistance [K/W]
def abs_conv_res(h_w_m2_k, surf_m2):                                                                                    # abs_conv_res(Heat transfer coefficient [W/(m^2*K)], Surface [m^2])
  if (h_w_m2_k > 0 and surf_m2 > 0):                                                                                    # Check in-vals consistency, if ok
    return 1/(h_w_m2_k*surf_m2)                                                                                         # Return calculated absolute convective resistance [K/W]
  else:                                                                                                                 # Else if in-vals consistency ain't ok
    return -1                                                                                                           # Return err val

# Function definition to calculate absolute conductive resistance [K/W]
def abs_cond_res(pipe_int_diam, pipe_ext_diam, pipe_len, pipes_num, therm_cond):                                        # abs_cond_res(Pipe internal diameter [m], Pipe external diameter [m], Pipe length [m], Number of pipes, Thermal conductivity lambda [W/(m*K)])
  if (pipe_int_diam > 0 and pipe_ext_diam > 0 and pipe_len > 0 and pipes_num > 0 and therm_cond > 0):                   # Check in-vals consistency, if ok
    return  (mt.log(pipe_ext_diam/pipe_int_diam, mt.e))/(2*mt.pi*therm_cond*pipe_len*pipes_num)                         # Return calculated absolute conductive resistance [K/W]
  else:                                                                                                                 # Else if in-vals consistency ain't ok
    return -1                                                                                                           # Return err val

# Function definition to calculate pipe internal and external surface temperatures [°C]
# using q=delta(T)/R <=> q=(Tm-Ts)/R  =>  Ts=Tm-R*q [°C]
def calc_pipe_int_ext_surf_temp(therm_pow_kw, avg_int_fl_temp, int_fl_conv_r, pipe_cond_r):                             # calc_pipe_int_ext_surf_temp(Exchanged thermal power [kW], Average interal fluid temperature [°C], Internal fluid absolute convective resistance [K/W], Pipe absolute conductive resistance [K/W])
  if (therm_pow_kw != 0 and avg_int_fl_temp >= zero_k and int_fl_conv_r > 0 and pipe_cond_r > 0):                       # Check in-vals consistency, if ok
    therm_pow_w = abs(therm_pow_kw)*1*1e3                                                                               # Thermal power conversion from [kW] into [W]
    int_surf_temp = avg_int_fl_temp-(int_fl_conv_r*therm_pow_w)                                                         # Calculate pipe internal surface temperature [°C]  -->  con [K/W]=[°C/W]  -->  [°C]-([°C/W]*[W])=[°C]
    ext_surf_temp = int_surf_temp-(pipe_cond_r*therm_pow_w)                                                             # Calculate pipe external surface temperature [°C]  -->  con [K/W]=[°C/W]  -->  [°C]-([°C/W]*[W])=[°C]
    return int_surf_temp, ext_surf_temp                                                                                 # Return calculated pipe internal and external surface temperatures [°C]
  else:                                                                                                                 # Else if in-vals consistency ain't ok
    return -1, -1                                                                                                       # Return err vals

# Function definition to calculate overall heat transfer coefficient [kW/(m^2*K)] - EQUATION 5.143
def overall_htc_kw_m2_k(surf_m2, abs_int_conv_res_k_w, abs_cond_res_k_w, abs_ext_conv_res_k_w):                         # overall_htc_kw_m2_k(Surface [m^2], Absolute internal convective resistance [K/W], Absolute conductive resistance [K/W], Absolute external convective resistance [K/W])
  if (surf_m2 > 0 and abs_int_conv_res_k_w > 0 and abs_cond_res_k_w > 0 and abs_ext_conv_res_k_w > 0):                  # Check in-vals consistency, if ok
    return ((1/surf_m2)/(abs_int_conv_res_k_w+abs_cond_res_k_w+abs_ext_conv_res_k_w))*1*1e-3                            # Return calculated overall heat transfer coefficient [kW/(m^2*K)]
  else:                                                                                                                 # Else if in-vals consistency ain't ok
    return -1                                                                                                           # Return err val

# Function definition to calculate transferred heat (thermal power) [kW]
# using overall heat transfer coefficient
def therm_pow_kw_overall_htc(ova_htc_kw_m2_k, surf_m2, lmtd):                                                           # therm_pow_kw_overall_htc(Overall heat transfer coefficient [kW/(m^2*K)], Surface [m^2], Log-mean temperature difference [°C] or [K])
  if (ova_htc_kw_m2_k > 0 and surf_m2 > 0 and lmtd != 0):                                                               # Check in-vals consistency, if ok
    return ova_htc_kw_m2_k*surf_m2*lmtd                                                                                 # Return thermal power [kW]
  else:                                                                                                                 # Else if in-vals consistency ain't ok
    return -1                                                                                                           # Return err val

# Function definition to calculate Grashof adimensional number
def shell_gr(pipe_ext_diam, beta, surf_temp, env_temp, fluid_kin_vis):                                                  # gr(Pipe external diameter [m], Beta: thermal expansion coefficient [1/K], Surface temperature [K], Environment temperature [K], Fluid kinematic viscosity [[m^2/s]])
  if (pipe_ext_diam > 0 and beta != None and surf_temp >= zero_k and env_temp >= zero_k and fluid_kin_vis > 0):         # Check in-vals consistency, if ok
    g = 9.81                                                                                                            # Gravity acceleration value [m/s^2]
    return (g*mt.pow(pipe_ext_diam, 3)*beta*abs(surf_temp-env_temp))/mt.pow(fluid_kin_vis, 2)                           # Return calculated Grashof number [adimensional]  -->  ([m/s^2]*[m^3]*[1/K]*[K])/[m^2/s]=[adimensional]
  else:                                                                                                                 # Else if in-vals consistency ain't ok
    return -1                                                                                                           # Return err val

# Function definition to calculate Rayleigh adimensional number
def shell_ra(gr, pr):                                                                                                   # shell_ra(Grashof number, Prandtl number)
  if (gr != 0 and pr != 0):                                                                                             # Check in-vals consistency, if ok
    return gr*pr                                                                                                        # Return calculated Rayleigh number [adimensional]  -->  [adimensional]*[adimensional]=[adimensional]
  else:                                                                                                                 # Else if in-vals consistency ain't ok
    return -1                                                                                                           # Return err val

# Function definition to calculate Nusselt adimensional number from Rayleigh adimensional number
def shell_nu_from_ra(ra):                                                                                               # shell_nu_from_ra(Rayleigh number)
  if (ra > 0 and ra < 1*1e9):                                                                                           # If Rayleigh number is in formula range
    return 0.59*mt.pow(ra, 1/4)                                                                                         # Return calculated Rayleigh number [adimensional]
  else:                                                                                                                 # Else if Rayleigh number ain't in formula range
    return -1                                                                                                           # Return err val
