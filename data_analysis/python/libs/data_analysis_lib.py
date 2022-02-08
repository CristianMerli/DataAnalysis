########
# LIBS #
########

# Libraries import
import pandas as pd                                                                                                     # Data-analysis panda lib
import enum as en                                                                                                       # Enum lib
# Project personal libraries import
import libs.plotting_lib as pl                                                                                          # Plotting lib
import libs.output_lib as out                                                                                           # Output lib

########
# VARS #
########

# Dataset file-vars: original (.dat) file and formatted (.csv) file
dat_data_filepath = "../dataset/scambiatore26112021_5.dat"                                                              # Original dataset filepath (.dat)
csv_data_filepath = "../dataset/scambiatore26112021_5.csv"                                                              # Formatted dataset filepath (generate new .csv file)
# Modification-vars to perform while converting and formatting (.dat) file into (.csv) file
old_sep_chr = '\t'                                                                                                      # Sep chr in (.dat) file
csv_sep_chr = ';'                                                                                                       # Sep chr in (.csv) file
comments_str = "% "                                                                                                     # Delete comment-strings in (.dat) file when creatin' (.csv) file
oth_op_old_lbl = "Altra operazione"                                                                                     # Replace 'other oper' lbl
# Heat-exchanger operation label-vars
oth_op_lbl = "Altra-operazione"                                                                                         # Other oper lbl
cocurrent_flow_lbl = "Equicorrente"                                                                                     # Cocurrent flow lbl
countercurrent_flow_lbl = "Controcorrente"                                                                              # Countercurrent flow lbl
# Dataset columns vars
time_col = "Time(s)"                                                                                                    # Time col in dataset (time ellapsed vals [s])
f1_col = "F1(l/h)"                                                                                                      # Volume flow rate 1 col in dataset (cold fluid volume flow rate vals [l/h])
f2_col = "F2(l/h)"                                                                                                      # Volume flow rate 2 col in dataset (hot fluid volume flow rate vals [l/h])
t1_col = "T1(degC)"                                                                                                     # Temp 1 col in dataset (cold-in fluid temp vals [°C])
t2_col = "T2(degC)"                                                                                                     # Temp 2 col in dataset (hot-in fluid temp vals [°C])
t3_col = "T3(degC)"                                                                                                     # Temp 3 col in dataset (cold-out fluid temp vals [°C])
t4_col = "T4(degC)"                                                                                                     # Temp 4 col in dataset (hot-out fluid temp vals [°C])
conf_col = "Configurazione"                                                                                             # Config col in dataset
# Measures strings definition
meas_str = "Measure "                                                                                                   # Measure string def (ONLY FOR PRINTS/PLOTS B4 MEASURES DEF)
meas_typ_str = ["Undefined measure", "Cocurrent measure", "Countercurrent measure"]                                     # Measure-types strings array

########
# DEFS #
########

# Measure type enum definition
class Meas_typ(en.Enum):                                                                                                # Measure type enum class
  undef = 0                                                                                                             # Undefined measure type
  ccurr = 1                                                                                                             # Cocurrent measure type
  cntcurr = 2                                                                                                           # Countercurrent measure type
# Measure-variables class
class Meas_vars:                                                                                                        # Measure-vars class (attributes, constructor, methods)
  typ = Meas_typ.undef                                                                                                  # Measure type
  typ_counter = 0                                                                                                       # Measure type counter
  name = ""                                                                                                             # Measure name
  f1 = 0.0                                                                                                              # Mean F1 var value --> Cold fluid volume flow rate [l/h] --> mass flow rate [kg/s]
  f2 = 0.0                                                                                                              # Mean F2 var value --> Hot fluid volume flow rate [l/h] --> mass flow rate [kg/s]
  t1 = 0.0                                                                                                              # Mean T1 var value --> Cold fluid inlet temperature [°C]
  t2 = 0.0                                                                                                              # Mean T2 var value --> Hot fluid inlet temperature [°C]
  t3 = 0.0                                                                                                              # Mean T3 var value --> Cold fluid outlet temperature [°C]
  t4 = 0.0                                                                                                              # Mean T4 var value --> Hot fluid outlet temperature [°C]
  avg_cold_fl_temp = 0.0                                                                                                # Average cold fluid temperature [°C]
  avg_hot_fl_temp = 0.0                                                                                                 # Average hot fluid temperature [°C]
  cold_fl_delta_temp = 0.0                                                                                              # Cold fluid delta temperature [°C] or [K]
  hot_fl_delta_temp = 0.0                                                                                               # Hot fluid delta temperature [°C] or [K]
  lmtd = 0.0                                                                                                            # Log-mean temperature difference (LMTD) [°C] or [K] (approximative)
  cold_fl_tr_heat = 0.0                                                                                                 # Cold fluid transferred heat (thermal power) [kW]
  hot_fl_tr_heat = 0.0                                                                                                  # Hot fluid transferred heat (thermal power) [kW]
  heat_loss = 0.0                                                                                                       # Heat loss (thermal power) [kW]
  glob_htc = 0.0                                                                                                        # Global heat transfer coefficient (global HTC) [kW/(m^2*K)]
  cpt_min = 0.0                                                                                                         # C-point-min min[mass-flow-rate*Cp] [kJ/(K*s)] (approximative)
  cpt_max = 0.0                                                                                                         # C-point-max max[mass-flow-rate*Cp] [kJ/(K*s)] (approximative)
  ntu = 0.0                                                                                                             # Number of transfer units (NTU) (approximative)
  epsilon = 0.0                                                                                                         # Effectiveness (epsilon) (approximative)
  re_int = 0.0                                                                                                          # Reynolds number inside steel pipes
  nu_int = 0.0                                                                                                          # Nusselt number inside steel pipes
  pr_int = 0.0                                                                                                          # Prandtl number inside steel pipes
  h_int = 0.0                                                                                                           # Heat transfer coefficient (h) inside steel pipes
  re_ext = 0.0                                                                                                          # Reynolds number inside glass pipe
  nu_ext = 0.0                                                                                                          # Nusselt number inside glass pipe
  pr_ext = 0.0                                                                                                          # Prandtl number inside glass pipe
  h_ext = 0.0                                                                                                           # Heat transfer coefficient (h) inside glass pipe
  def __init__(self):                                                                                                   # Constructor
    return                                                                                                              # Return nothing
  def get_info(self):                                                                                                   # Measure class method to get measure info
    dbg_str = ("\n--> "+self.name+" calculations results:"\
    +"\n- Cold fluid mass flow rate: "+str(self.f1)+" [kg/s]"\
    +"\n- Hot fluid mass flow rate: "+str(self.f2)+" [kg/s]"\
    +"\n- Cold fluid inlet temperature: "+str(self.t1)+" [°C]"\
    +"\n- Hot fluid inlet temperature: "+str(self.t2)+" [°C]"\
    +"\n- Cold fluid outlet temperature: "+str(self.t3)+" [°C]"\
    +"\n- Hot fluid outlet temperature: "+str(self.t4)+" [°C]"\
    +"\n- Average cold fluid temperature: "+str(self.avg_cold_fl_temp)+" [°C]"\
    +"\n- Average hot fluid temperature: "+str(self.avg_hot_fl_temp)+" [°C]"\
    +"\n- Cold fluid delta temperature: "+str(self.cold_fl_delta_temp)+" [°C]"\
    +"\n- Hot fluid delta temperature: "+str(self.hot_fl_delta_temp)+" [°C]"\
    +"\n- Log-mean temperature difference (LMTD): "+str(self.lmtd)+" [°C] (approximative)"\
    +"\n- Cold fluid transferred heat (thermal power): "+str(self.cold_fl_tr_heat)+" [kW]"\
    +"\n- Hot fluid transferred heat (thermal power): "+str(self.hot_fl_tr_heat)+" [kW]"\
    +"\n- Heat losses (thermal power): "+str(self.heat_loss)+" [kW]"\
    +"\n- Global heat transfer coefficient (global HTC): "+str(self.glob_htc)+" [kW/(m^2*K)] (approximative)"\
    +"\n- C-point-min min[mass-flow-rate*Cp]: "+str(self.cpt_min)+" [kJ/(K*s)] (approximative)"\
    +"\n- C-point-max max[mass-flow-rate*Cp]: "+str(self.cpt_max)+" [kJ/(K*s)] (approximative)"\
    +"\n- Number of transfer units (NTU): "+str(self.ntu)+" (approximative)"\
    +"\n- Effectiveness (epsilon): "+str(self.epsilon)+" (approximative)"\
    +"\n- Reynolds number inside steel pipes: "+str(self.re_int)\
    +"\n- Nusselt number inside steel pipes: "+str(self.nu_int)\
    +"\n- Prandtl number inside steel pipes: "+str(self.pr_int)\
    +"\n- Heat transfer coefficient (h) inside steel pipes: "+str(self.h_int)+" [W/(m^2*K)]"\
    +"\n- Reynolds number inside glass pipe: "+str(self.re_ext)\
    +"\n- Nusselt number inside glass pipe: "+str(self.nu_ext)\
    +"\n- Prandtl number inside glass pipe: "+str(self.pr_ext)\
    +"\n- Heat transfer coefficient (h) inside glass pipe: "+str(self.h_ext)+" [W/(m^2*K)]\n")                          # Dbg fbk
    return dbg_str                                                                                                      # Return dbg fbk

##########
# FUNCTS #
##########

# Function definition to convert and format dataset: open and manipulate data inside
# (.dat) file, write formatted data inside (.csv) file and return DataFrame var
def load_dataset_data():                                                                                                # load_dataset_data()
  with open(dat_data_filepath,'r') as src_fl:                                                                           # Open src file in read mode: (.dat) file
    with open(csv_data_filepath,'w') as dest_fl:                                                                        # Open dest file in write mode: (.csv) file
        next(src_fl)                                                                                                    # Skip header line in src file
        for line in src_fl:                                                                                             # Read src file line-by-line
            line = line.replace(oth_op_old_lbl, oth_op_lbl)                                                             # Replace target string (other operation label)
            line = line.replace(old_sep_chr, csv_sep_chr)                                                               # Replace data separator-chars
            line = line.replace(comments_str, '')                                                                       # Delete target string (comment string in src file)
            dest_fl.write(line)                                                                                         # Write each formatted data line inside dest file (.csv)
  return pd.read_csv(csv_data_filepath, sep=csv_sep_chr, encoding="utf8")                                               # Return DataFrame var containing measure dataset

# Data intervals detection: extract data in different operating conditions
# and save delimiting idxs
def find_plt_save_measures(he_data, plt_flg, dbg_flg):                                                                  # find_plt_save_measures(Heat-exchanger DataFrame, Plotting flag, Debug flag)
  he_measures_data = list()                                                                                             # New heat-exchanger measures data list declaration and following definition
  delim_idxs = list()                                                                                                   # New delimiter-indexes list (list containing dataframe row idxs: configurations start/end idxs)
  old_oper_str = str()                                                                                                  # Operation string in previous dataframe row
  oper_str = str()                                                                                                      # Operation string in current dataframe row
  rows_scroll_index = 0                                                                                                 # Index to trace current row in dataframe rows scrolling cycle
  for rows_scroll_index, row in he_data.iterrows():                                                                     # Cycle to scroll rows in dataframe, tracing row index
    old_oper_str = oper_str                                                                                             # Update previous dataframe row string
    oper_str = row[conf_col]                                                                                            # Update current dataframe row string
    if (rows_scroll_index > 1):                                                                                         # Skip first row (previous and current dataframe row strings are initialized with the same val)
      if ((old_oper_str == oth_op_lbl) & (oper_str == cocurrent_flow_lbl)):                                             # Cocurrent flow configuration start detecting condition
        delim_idxs.append(rows_scroll_index)                                                                            # Insert current row index in delimiter-indexes list
      if ((old_oper_str == cocurrent_flow_lbl) & (oper_str == oth_op_lbl)):                                             # Cocurrent flow configuration end detecting condition
        delim_idxs.append(rows_scroll_index)                                                                            # Insert current row index (+1 for data split) in delimiter-indexes list
      if ((old_oper_str == oth_op_lbl) & (oper_str == countercurrent_flow_lbl)):                                        # Countercurrent flow configuration start detecting condition
        delim_idxs.append(rows_scroll_index)                                                                            # Insert current row index in delimiter-indexes list
      if ((old_oper_str == countercurrent_flow_lbl) & (oper_str == oth_op_lbl)):                                        # Countercurrent flow configuration end detecting condition
        delim_idxs.append(rows_scroll_index)                                                                            # Insert current row index (+1 for data split) in delimiter-indexes list
  if (len(delim_idxs)%2 == 1):                                                                                          # In case last configuration lasts 'till dataframe tail row
      delim_idxs.append(rows_scroll_index+1)                                                                            # Set last configuration end idx as tail row idx inside delimiter-idxs list (+1 for data split)
  for i in range(1, len(delim_idxs), 2):                                                                                # Heat-exchanger measures data extraction cycle
    he_measures_data.append(he_data[delim_idxs[i-1]:delim_idxs[i]])                                                     # Extract and add measures (btwn indexs in delimiter-indexes list)
  for he_data_measure in he_measures_data:                                                                              # Heat-exchanger measures data scrollin' cycle
    he_data_measure = he_data_measure[~he_data_measure[conf_col].isin([oth_op_lbl])]                                    # Remove eventual data corresponding to other-operation in each heat-exchanger data measure in list
  if (plt_flg):                                                                                                         # If plotting flag is ena
    pl.plot_data_flt(he_data, pl.plt_title, None, delim_idxs, None, pl.Plt_mode.complete)                               # Function call to graphically plot data filtering operations (complete plottin' mode)
  if (dbg_flg):                                                                                                         # If dbg flg is ena
    dbg_str = ("\n--> Measures-delimiting indexes: "+str(delim_idxs)+'\n')                                              # Dbg fbk
    print(dbg_str)                                                                                                      # Print dbg fbk
    out.save_output(out.Output_typ.delim_idxs, dbg_str)                                                                 # Save dbg output
  return he_measures_data, delim_idxs                                                                                   # Return heat-exchanger measures data list and measures delimiters indexes list

# Function definition to find min value index in vals list (non-zero idx)
def find_min_idx(vals):                                                                                                 # find_min_idx(Vals list)
  min_idx = 0                                                                                                           # Min val index init
  for i in range(len(vals)):                                                                                            # Vals scrollin' cycle
    if (vals[i] < vals[min_idx]):                                                                                       # Min val idx upd cond
      min_idx = i                                                                                                       # Min val idx upd oper
  return min_idx                                                                                                        # Return min val idx in vals list (non-zero idx)

# Function definition to find and plot optimal steady-conditions data windows, by calulating mean
# standard-deviations in each window, avoiding the first transitory data-window and saving std-devs output
def find_plt_save_stdy_cond_win(dbs, win_span, plt_flg, dbg_flg):                                                       # find_plt_save_stdy_cond_win(Measures-datablocks to split, Data-windows span [samples], Plotting flag, Debug flag)
  sc_windows = list()                                                                                                   # New steady-conditions data-windows list declaration and following definition
  idx = 0                                                                                                               # Measure index
  if (dbg_flg):                                                                                                         # If dbg flg is ena
    dbg_str = str()                                                                                                     # Dbg fbk
  for db in dbs:                                                                                                        # Measures-data datablocks scrollin' cycle
    if (dbg_flg):                                                                                                       # If dbg flg is ena
      dbg_str += ("\n-----------------------------------------------------------------------"\
      +"\n--> NEW 'find_stdy_cond()' FUNCTION CALL FOR "+str(meas_str+str(idx+1))\
      +"\n-----------------------------------------------------------------------\n")                                   # Dbg fbk
    win_list_size = len(db) // win_span                                                                                 # Data-windows list size (zero-idx)
    max_win_idx = db.index[-1]-db.index[0]                                                                              # Max data-windows dataframe-idx
    windows = list()                                                                                                    # Data-windows list
    win_stddevs = list()                                                                                                # Data-windows mean standard deviations list to detect the smallest one
    win_start_idxs = list()                                                                                             # Data-windows starting-idxs (by data-windows span [samples])
    win_end_idxs = list()                                                                                               # Data-windows ending-idxs (by data-windows span [samples])
    for i in range(win_list_size, 0, -1):                                                                               # Cycle to scroll data-windows backwards, excluding the first incomplete set (tail-transitory)
      win_start_idxs.append(max_win_idx-((i-1)*win_span)-win_span+1)                                                    # Data-windows starting-idx (by data-windows span [samples])
      win_end_idxs.append(max_win_idx-((i-1)*win_span)+1)                                                               # Data-windows ending-idx (by data-windows span [samples])
      windows.append(db[win_start_idxs[-1]:win_end_idxs[-1]])                                                           # Populate data-windows list
      f1_stddev = windows[-1][f1_col].std()                                                                             # Calc f1 data standard-deviation
      f2_stddev = windows[-1][f2_col].std()                                                                             # Calc f2 data standard-deviation
      t1_stddev = windows[-1][t1_col].std()                                                                             # Calc t1 data standard-deviation
      t2_stddev = windows[-1][t2_col].std()                                                                             # Calc t2 data standard-deviation
      t3_stddev = windows[-1][t3_col].std()                                                                             # Calc t3 data standard-deviation
      t4_stddev = windows[-1][t4_col].std()                                                                             # Calc t4 data standard-deviation
      win_stddevs.append((f1_stddev+f2_stddev+t1_stddev+t2_stddev+t3_stddev+t4_stddev) / 6)                             # Calc the data mean standard deviations, to look for the best data window (steady cond)
      if (dbg_flg):                                                                                                     # If dbg flg is ena
        dbg_str += ("\nf1_stddev: "+str(f1_stddev)+"\nf2_stddev: "+str(f2_stddev)+"\nt1_stddev: "+str(t1_stddev)\
        +"\nt2_stddev: "+str(t2_stddev)+"\nt3_stddev: "+str(t3_stddev)+"\nt4_stddev: "+str(t4_stddev)\
        +"\nmean_stddev: "+str(win_stddevs[-1])+'\n')                                                                   # Dbg fbk
    min_stddevs_win_idx = find_min_idx(win_stddevs)                                                                     # Function call to find data-window with min standard deviations avg in data-windows list by-index (non-zero idx)
    sc_windows.append(db[win_start_idxs[min_stddevs_win_idx]:win_end_idxs[min_stddevs_win_idx]])                        # Add optimal steady-conditions data window for each heat-exchanger data measure
    if (plt_flg):                                                                                                       # If plotting flag is ena
      pl.plot_data_flt(db, meas_str+str(idx+1), win_start_idxs,
                       win_end_idxs, min_stddevs_win_idx, pl.Plt_mode.detailed)                                         # Function call to graphically plot data filtering operations (detailed plottin' mode)
    if (dbg_flg):                                                                                                       # If dbg flg is ena
      dbg_str += ("\nmin_stddevs_datablocks_idx: "+str(min_stddevs_win_idx)+'\n')                                       # Dbg fbk
    idx += 1                                                                                                            # Measure index upd
  if (dbg_flg):                                                                                                         # If dbg flg is ena
    print(dbg_str)                                                                                                      # Print dbg fbk
    out.save_output(out.Output_typ.std_devs, dbg_str)                                                                   # Save dbg output
  return sc_windows                                                                                                     # Return steady-conditions data-windows

# Function definition to determine and save measured variables values
# (temperatures and volume flow rates) and define measure-vars data-structures list
def def_save_meas_vars(sc_windows, dbg_flg):                                                                            # def_save_meas_vars(Steady-conditions data-windows, Debug flag)
  measures = list()                                                                                                     # New measures list declaration and following definition
  ccurr_meas_counter = 0                                                                                                # Cocurrent-measures counter
  cntcurr_meas_counter = 0                                                                                              # Countercurrent-measures counter
  undef_meas_counter = 0                                                                                                # Undefined-measures counter
  if (dbg_flg):                                                                                                         # If dbg flg is ena
    dbg_str = str()                                                                                                     # Dbg fbk
  for sc_win in sc_windows:                                                                                             # Steady-conditions data-windows scrollin' cycle
    measure = Meas_vars()                                                                                               # Define new measure-values data-structure
    if (sc_win[conf_col].iloc[0] == cocurrent_flow_lbl):                                                                # Cocurrent measure detectin' cond
      ccurr_meas_counter += 1                                                                                           # Cocurrent-measures counter upd
      measure.typ = Meas_typ.ccurr                                                                                      # Measure type definition
      measure.typ_counter = ccurr_meas_counter                                                                          # Measure type counter definition
      measure.name = meas_typ_str[measure.typ.value]+" "+str(measure.typ_counter)                                       # Measure name definition
    elif (sc_win[conf_col].iloc[0] == countercurrent_flow_lbl):                                                         # Countercurrent measure detectin' cond
      cntcurr_meas_counter += 1                                                                                         # Countercurrent-measures counter upd
      measure.typ = Meas_typ.cntcurr                                                                                    # Measure type definition
      measure.typ_counter = cntcurr_meas_counter                                                                        # Measure type counter definition
      measure.name = meas_typ_str[measure.typ.value]+" "+str(measure.typ_counter)                                       # Measure name definition
    else:                                                                                                               # Undefined measure detectin' cond
      undef_meas_counter += 1                                                                                           # Undefined-measures counter upd
      measure.typ = Meas_typ.undef                                                                                      # Measure type definition
      measure.typ_counter = undef_meas_counter                                                                          # Measure type counter definition
      measure.name = meas_typ_str[measure.typ.value]+" "+str(measure.typ_counter)                                       # Measure name definition
    measure.f1 = sc_win[f1_col].mean()                                                                                  # Calc mean F1 var value in optimal steady-conditions data-window and populate measure data-structure
    measure.f2 = sc_win[f2_col].mean()                                                                                  # Calc mean F2 var value in optimal steady-conditions data-window and populate measure data-structure
    measure.t1 = sc_win[t1_col].mean()                                                                                  # Calc mean T1 var value in optimal steady-conditions data-window and populate measure data-structure
    measure.t2 = sc_win[t2_col].mean()                                                                                  # Calc mean T2 var value in optimal steady-conditions data-window and populate measure data-structure
    measure.t3 = sc_win[t3_col].mean()                                                                                  # Calc mean T3 var value in optimal steady-conditions data-window and populate measure data-structure
    measure.t4 = sc_win[t4_col].mean()                                                                                  # Calc mean T4 var value in optimal steady-conditions data-window and populate measure data-structure
    measures.append(measure)                                                                                            # Add measure data-structure in measures list
    if (dbg_flg):                                                                                                       # If dbg flg is ena
      dbg_str += ("\n--> "+str(measure.name)+" mean vals:"+"\nF1: "+str(measure.f1)+" [l/h]"\
      +"\nF2: "+str(measure.f2)+" [l/h]"+"\nT1: "+str(measure.t1)+" [°C]"+"\nT2: "+str(measure.t2)+" [°C]"\
      +"\nT3: "+str(measure.t3)+" [°C]"+"\nT4: "+str(measure.t4)+" [°C]\n")                                             # Dbg fbk  
  if (dbg_flg):                                                                                                         # If dbg flg is ena
    print(dbg_str)                                                                                                      # Print dbg fbk
    out.save_output(out.Output_typ.measures, dbg_str)                                                                   # Save dbg output
  return measures                                                                                                       # Return measure-vars data-structures list
