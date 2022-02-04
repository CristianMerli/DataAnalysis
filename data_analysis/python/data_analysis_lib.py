########
# LIBS #
########

# Libraries import
import pandas as pd                                                                                     # Data-analysis panda lib
# Project personal libraries import
import plotting_lib as pl                                                                               # Plotting lib

########
# VARS #
########

# Dataset file-vars: original (.dat) file and formatted (.csv) file
dat_data_filepath = "../dataset/scambiatore26112021_5.dat"                                              # Original dataset filepath (.dat)
csv_data_filepath = "../dataset/scambiatore26112021_5.csv"                                              # Formatted dataset filepath (generate new .csv file)
# Modification-vars to perform while converting and formatting (.dat) file into (.csv) file
old_sep_chr = '\t'                                                                                      # Sep chr in (.dat) file
csv_sep_chr = ';'                                                                                       # Sep chr in (.csv) file
comments_str = "% "                                                                                     # Delete comment-strings in (.dat) file when creatin' (.csv) file
oth_op_old_lbl = "Altra operazione"                                                                     # Replace 'other oper' lbl
# Heat-exchanger operation label-vars
oth_op_lbl = "Altra-operazione"                                                                         # Other oper lbl
cocurrent_flow_lbl = "Equicorrente"                                                                     # Cocurrent flow lbl
countercurrent_flow_lbl = "Controcorrente"                                                              # Countercurrent flow lbl
# Dataset columns vars
time_col = "Time(s)"                                                                                    # Time col in dataset (time ellapsed vals [s])
f1_col = "F1(l/h)"                                                                                      # Volume flow rate 1 col in dataset (cold fluid volume flow rate vals [l/h])
f2_col = "F2(l/h)"                                                                                      # Volume flow rate 2 col in dataset (hot fluid volume flow rate vals [l/h])
t1_col = "T1(degC)"                                                                                     # Temp 1 col in dataset (cold-in fluid temp vals [°C])
t2_col = "T2(degC)"                                                                                     # Temp 2 col in dataset (hot-in fluid temp vals [°C])
t3_col = "T3(degC)"                                                                                     # Temp 3 col in dataset (cold-out fluid temp vals [°C])
t4_col = "T4(degC)"                                                                                     # Temp 4 col in dataset (hot-out fluid temp vals [°C])
conf_col = "Configurazione"                                                                             # Config col in dataset
# Measure-names strings array definition
meas_names = ["Cocurrent measure 1", "Countercurrent measure 1",
              "Countercurrent measure 2", "Cocurrent measure 2"]                                        # Measure-names strings array def

########
# DEFS #
########

# Data-structure definition to contain measures variables values
class Meas_vars:                                                                                        # Measure vars data-structure class
  f1: float                                                                                             # Mean F1 var value --> Cold fluid volume flow rate [l/h] --> Cold fluid mass flow rate [kg/s]
  f2: float                                                                                             # Mean F2 var value --> Hot fluid volume flow rate [l/h] --> Hot fluid mass flow rate [kg/s]
  t1: float                                                                                             # Mean T1 var value --> Cold fluid inlet temperature [°C]
  t2: float                                                                                             # Mean T2 var value --> Hot fluid inlet temperature [°C]
  t3: float                                                                                             # Mean T3 var value --> Cold fluid outlet temperature [°C]
  t4: float                                                                                             # Mean T4 var value --> Hot fluid outlet temperature [°C]
  avg_cold_fl_temp:float                                                                                # Average cold fluid temperature [°C]
  avg_hot_fl_temp:float                                                                                 # Average hot fluid temperature [°C]
  cold_fl_delta_temp:float                                                                              # Cold fluid delta temperature [°C] or [K]
  hot_fl_delta_temp:float                                                                               # Hot fluid delta temperature [°C] or [K]
  lmtd:float                                                                                            # Log-mean temperature difference (LMTD) [°C] or [K]
  cold_fl_tr_heat:float                                                                                 # Cold fluid transferred heat (thermal power) [kW]
  hot_fl_tr_heat:float                                                                                  # Hot fluid transferred heat (thermal power) [kW]
  heat_loss:float                                                                                       # Heat loss (thermal power) [kW]
  glob_htc:float                                                                                        # Global heat transfer coefficient (global HTC) [kW/(m^2*K)]

##########
# FUNCTS #
##########

# Function definition to convert and format dataset: open and manipulate data inside
# (.dat) file, write formatted data inside (.csv) file and return DataFrame var
def load_dataset_data():                                                                                # load_dataset_data()
  with open(dat_data_filepath,'r') as src_fl:                                                           # Open src file in read mode: (.dat) file
    with open(csv_data_filepath,'w') as dest_fl:                                                        # Open dest file in write mode: (.csv) file
        next(src_fl)                                                                                    # Skip header line in src file
        for line in src_fl:                                                                             # Read src file line-by-line
            line = line.replace(oth_op_old_lbl, oth_op_lbl)                                             # Replace target string (other operation label)
            line = line.replace(old_sep_chr, csv_sep_chr)                                               # Replace data separator-chars
            line = line.replace(comments_str, '')                                                       # Delete target string (comment string in src file)
            dest_fl.write(line)                                                                         # Write each formatted data line inside dest file (.csv)
  return pd.read_csv(csv_data_filepath, sep=csv_sep_chr, encoding="utf8")                               # Return 

# Function definition to find min value index in vals list (non-zero idx)
def find_min_idx(vals):                                                                                 # find_min_idx(Vals list)
  min_idx = 0                                                                                           # Min val index init
  for i in range(len(vals)):                                                                            # Vals scrollin' cycle
    if (vals[i] < vals[min_idx]):                                                                       # Min val idx upd cond
      min_idx = i                                                                                       # Min val idx upd oper
  return min_idx                                                                                        # Rerturn min val idx in vals list (non-zero idx)

# Function definition to find and plotoptimal steady conditions data windows, by calulating
# mean standard-deviations in datablocks and avoiding first transitory datablock
def find_stdy_cond_plt(db, dbs_span, dbg_flg, call_str):                                                # find_stdy_cond(Datablock to split, Datablocks span [samples], Debug flag, Function call string)
  if (dbg_flg):                                                                                         # If dbg flg is ena
    print("\n-----------------------------------------------------------------------")                  # Print dbg fbk
    print("--> NEW 'find_stdy_cond()' FUNCTION CALL FOR", call_str)                                     # Print dbg fbk
    print("-----------------------------------------------------------------------\n")                  # Print dbg fbk
  dbs = list()                                                                                          # Datablocks list
  stddevs_dbs = list()                                                                                  # Datablocks mean standard deviations list to detect the smallest one
  start_idxs_dbs = list()                                                                               # Datablocks starting-idxs (by datablocks span [samples])
  end_idxs_dbs = list()                                                                                 # Datablocks ending-idxs (by datablocks span [samples])
  num_datablocks = len(db) // dbs_span                                                                  # Datablocks list size (zero-idx)
  max_idx = db.index[-1]-db.index[0]                                                                    # Max datablock dataframe idx
  for i in range(num_datablocks, 0, -1):                                                                # Cycle to scroll datablocks backwards, excluding the first incomplete set (tail-transitory)
    start_idxs_dbs.append(max_idx-((i-1)*dbs_span)-dbs_span+1)                                          # Datablocks starting-idx (by datablocks span [samples])
    end_idxs_dbs.append(max_idx-((i-1)*dbs_span)+1)                                                     # Datablocks ending-idx (by datablocks span [samples])
    dbs.append(db[start_idxs_dbs[-1]:end_idxs_dbs[-1]])                                                 # Populate datablocks list
    f1_stddev = dbs[-1][f1_col].std()                                                                   # Calc f1 data standard-deviation
    f2_stddev = dbs[-1][f2_col].std()                                                                   # Calc f2 data standard-deviation
    t1_stddev = dbs[-1][t1_col].std()                                                                   # Calc t1 data standard-deviation
    t2_stddev = dbs[-1][t2_col].std()                                                                   # Calc t2 data standard-deviation
    t3_stddev = dbs[-1][t3_col].std()                                                                   # Calc t3 data standard-deviation
    t4_stddev = dbs[-1][t4_col].std()                                                                   # Calc t4 data standard-deviation
    stddevs_dbs.append((f1_stddev+f2_stddev+t1_stddev+t2_stddev+t3_stddev+t4_stddev) / 6)               # Calc the data mean standard deviations, to look for the best data window (steady cond)
    if (dbg_flg):                                                                                       # If dbg flg is ena
      print("f1_stddev: ", f1_stddev, "\nf2_stddev: ", f2_stddev, "\nt1_stddev: ", t1_stddev)           # Print dbg fbk
      print("t2_stddev: ", t2_stddev, "\nt3_stddev: ", t3_stddev, "\nt4_stddev: ", t4_stddev)           # Print dbg fbk
      print("mean_stddev: ", stddevs_dbs[-1], '\n')                                                     # Print dbg fbk
  min_stddevs_dbs_idx = find_min_idx(stddevs_dbs)                                                       # Function call to find min mean standard deviations in datablocks (non-zero idx)
  if (dbg_flg):                                                                                         # If dbg flg is ena
    print("min_stddevs_datablocks_idx: ", min_stddevs_dbs_idx)                                          # Print dbg fbk
  pl.plot_data_flt(db, call_str, start_idxs_dbs, end_idxs_dbs,
                   min_stddevs_dbs_idx, pl.Plt_mode.detailed)                                           # Function call to graphically plot data filtering operations (detailed plottin' mode)
  return db[start_idxs_dbs[min_stddevs_dbs_idx]:end_idxs_dbs[min_stddevs_dbs_idx]]                      # Return datablock with best steady conditions

# Function definition to determine measured variables values (temperatures and volume flow rates)
def calc_meas_vars_vals(sc_db_win_interval):                                                            # calc_meas_vars_vals(Steady-conditions datablock window interval)
  m_vars_vals = Meas_vars()                                                                             # Define new measure values data-structure
  m_vars_vals.f1 = sc_db_win_interval[f1_col].mean()                                                    # Calc mean F1 var value in passed datablock and populate data-structure
  m_vars_vals.f2 = sc_db_win_interval[f2_col].mean()                                                    # Calc mean F2 var value in passed datablock and populate data-structure
  m_vars_vals.t1 = sc_db_win_interval[t1_col].mean()                                                    # Calc mean T1 var value in passed datablock and populate data-structure
  m_vars_vals.t2 = sc_db_win_interval[t2_col].mean()                                                    # Calc mean T2 var value in passed datablock and populate data-structure
  m_vars_vals.t3 = sc_db_win_interval[t3_col].mean()                                                    # Calc mean T3 var value in passed datablock and populate data-structure
  m_vars_vals.t4 = sc_db_win_interval[t4_col].mean()                                                    # Calc mean T4 var value in passed datablock and populate data-structure
  return m_vars_vals                                                                                    # Return mean vars values data-structure
