########
# LIBS #
########

# Libraries import
import enum as en                                                                                                       # Enum lib
import numpy as np                                                                                                      # Data-analysis numpy lib
import matplotlib.pyplot as plt                                                                                         # Math plottin' lib
import seaborn as sns                                                                                                   # Seaborn plottin' lib
# Project personal libraries import
import libs.data_analysis_lib as da                                                                                     # Data analysis lib
import libs.eng_calcs_lib as ec                                                                                         # Engineering calcs lib
import libs.poly_approx_lib as pa                                                                                       # Poly-approximation lib

########
# VARS #
########

# Plotting settings-vars
plt_style = "seaborn-dark"                                                                                              # Plotting style
bkg_col = "#00071c"                                                                                                     # Background color (hex)
grid_col = "#1c338a"                                                                                                    # Grid color (hex)
text_col = "#e1e4ed"                                                                                                    # Text color (hex)
leg_txt_col = "#ffffff"                                                                                                 # Legend text color (hex)
plt_f1_col = "#00fff7"                                                                                                  # F1 data trend color (hex)
plt_f2_col = "#ff8000"                                                                                                  # F2 data trend color (hex)
plt_t1_col = "#ff1100"                                                                                                  # T1 data trend color (hex)
plt_t2_col = "#0ac700"                                                                                                  # T2 data trend color (hex)
plt_t3_col = "#f2ff00"                                                                                                  # T3 data trend color (hex)
plt_t4_col = "#ff00f2"                                                                                                  # T4 data trend color (hex)
tp_var_col = "#ff00e8"                                                                                                  # Thermophysical variable trend color (hex)
tp_intp_var_col = "#ff7c00"                                                                                             # Thermophysical variable interpolated/fitted trend color (hex)
dbs_c = "#0045c4"                                                                                                       # Datablocks lines color (hex)
dbs_sel_c = "#4fbd37"                                                                                                   # Selected datablocks lines color (hex)
dbs_sel_col = "#0cf700"                                                                                                 # Selected datablocks fill color (hex)
interval_box_col = "#ae00ff"                                                                                            # Interval text-box color (hex)
poly_approx_res_col = ["#ff0000", "#00a013"]                                                                            # Poly-approximation result text-box colors: Discarded/Accepted (hex)
plt_size_x = 26                                                                                                         # Plotting figure X size
plt_size_y = 14                                                                                                         # Plotting figure Y size
dbs_w = 3                                                                                                               # Datablocks lines width
dbs_sel_w = 9                                                                                                           # Selected datablocks lines width
dbs_sel_alpha = 0.15                                                                                                    # Selected datablocks fill alpha
interval_txt_y_offs_p = 70                                                                                              # Interval txt Y pos offset (percentage)
sel_interval_txt_y_offs_p = -10                                                                                         # Selected interval txt Y pos offset (percentage)
interval_txt_size = 20                                                                                                  # Interval text size
interval_box_alpha = 0.8                                                                                                # Interval text-box alpha
poly_approx_txt_size = 17                                                                                               # Poly-approximation text size
ploy_approx_box_alpha = 0.7                                                                                             # Poly-approximation text-box alpha
plt_line = "--"                                                                                                         # Plotting line type
plt_marker = 'o'                                                                                                        # Plotting marker type
# Data-acquisition plotting label-vars
plt_da_title = "Experimental data from LabView data-acquisition on heat-exchanger with measure-windows"                 # Title lbl
plt_time_lbl = "Time [s]"                                                                                               # X-axis lbl
plt_temp_flow_lbl = "Temperatures [°C]   /   Volume flow rates [l/h]"                                                   # Y-axis lbl
plt_f1_lbl = "F1 - Cold fluid volume flow rate [l/h]"                                                                   # Cold fluid vol flow rate (F1) lbl
plt_f2_lbl = "F2 - Hot fluid volume flow rate [l/h]"                                                                    # Hot fluid vol flow rate (F2) lbl
plt_t1_lbl = "T1 - Cold-in fluid temperature [°C]"                                                                      # Cold-in fluid temp (T1) lbl
plt_t2_lbl = "T2 - Hot-in fluid temperature [°C]"                                                                       # Hot-in fluid temp (T2) lbl
plt_t3_lbl = "T3 - Cold-out fluid temperature [°C]"                                                                     # Cold-out fluid temp (T3) lbl
plt_t4_lbl = "T4 - Hot-out fluid temperature [°C]"                                                                      # Hot-out fluid temp (T4) lbl
sel_interval_lbl = "SELECTED INTERVAL"                                                                                  # Selected interval lbl
# Thermophysical-variables plotting label-vars
plt_temp_lbl = "Temperature - T [°C]"                                                                                   # X-axis lbl: temperature
plt_rho_lbl = "density - ρ [$\mathregular{kg/m^3}$]"                                                                    # Y-axis lbl: thermophysical variable
plt_cp_lbl = "specific heat at constant pressure - Cp [kJ/(kg*K)]"                                                      # Y-axis lbl: thermophysical variable
plt_lambda_lbl = "thermal conductivity - λ [W/(m*K)]"                                                                   # Y-axis lbl: thermophysical variable
plt_ni_lbl = "kinematic viscosity - ν [$\mathregular{m^2/s}$]"                                                          # Y-axis lbl: thermophysical variable
plt_beta_lbl = "thermal expansion coefficient - β [1/K]"                                                                # Y-axis lbl: thermophysical variable
plt_pr_lbl = "prandtl number - Pr"                                                                                      # Y-axis lbl: thermophysical variable
plt_title_sep = "  /  "                                                                                                 # Title separator in str concat
plt_intp_rho_lbl = "interpolated density - ρ [$\mathregular{kg/m^3}$]"                                                  # Y-axis lbl: interpolated thermophysical variable
plt_intp_cp_lbl = "interpolated specific heat at constant pressure - Cp [kJ/(kg*K)]"                                    # Y-axis lbl: interpolated thermophysical variable
plt_intp_lambda_lbl = "interpolated thermal conductivity - λ [W/(m*K)]"                                                 # Y-axis lbl: interpolated thermophysical variable
plt_intp_ni_lbl = "interpolated kinematic viscosity - ν [$\mathregular{m^2/s}$]"                                        # Y-axis lbl: interpolated thermophysical variable
plt_intp_beta_lbl = "interpolated thermal expansion coefficient - β [1/K]"                                              # Y-axis lbl: interpolated thermophysical variable
plt_intp_pr_lbl = "interpolated Prandtl number - Pr"                                                                    # Y-axis lbl: interpolated thermophysical variable
plt_fit_rho_lbl = "fitted density - ρ [$\mathregular{kg/m^3}$]"                                                         # Y-axis lbl: poly-fitted thermophysical variable
plt_fit_cp_lbl = "fitted specific heat at constant pressure - Cp [kJ/(kg*K)]"                                           # Y-axis lbl: poly-fitted thermophysical variable
plt_fit_lambda_lbl = "fitted thermal conductivity - λ [W/(m*K)]"                                                        # Y-axis lbl: poly-fitted thermophysical variable
plt_fit_ni_lbl = "fitted kinematic viscosity - ν [$\mathregular{m^2/s}$]"                                               # Y-axis lbl: poly-fitted thermophysical variable
plt_fit_beta_lbl = "fitted thermal expansion coefficient - β [1/K]"                                                     # Y-axis lbl: poly-fitted thermophysical variable
plt_fit_pr_lbl = "fitted Prandtl number - Pr"                                                                           # Y-axis lbl: poly-fitted thermophysical variable
materials = ["Atm pressure air ", "Water ", "AISI-316 ", "Pyrex glass "]                                                # Materials-lbls (Air-atmp/Water/AISI-316-stainless-steel/Borosilicate (pyrex) glass)
poly_approx_res = ["BAD APPROXIMATION - DISCARDED", "GOOD APPROXIMATION - ACCEPTED"]                                    # Poly-approximation result (Discarded/Accepted)
# Approssimative temperature trend in heat-exchanger section plotting label-vars
plt_he_temp_trend_titles = [" approssimative temperature trend in heat-exchanger inlet (bottom) section",
                            " approssimative temperature trend in heat-exchanger central section",
                            " approssimative temperature trend in heat-exchanger outlet (top) section"]                 # Title lbls
plt_he_dist_lbl = "Distance - d [m]"                                                                                    # X-axis lbl
plt_he_temp_lbl = "Temperature - T [°C]"                                                                                # Y-axis lbl
he_temp_trend_lbls = ["approssimative temperature trend in heat-exchanger bottom section",
                      "approssimative temperature trend in heat-exchanger central section",
                      "approssimative temperature trend in heat-exchanger top section"]                                 # Heat-exchanger temperature trend lbls
he_steel_pipe_lbl = "Steel pipe profile"                                                                                # Heat-exchanger steel pipe profile lbl
he_glass_pipe_lbl = "Glass pipe profile"                                                                                # Heat-exchanger glass pipe profile lbl
he_temp_trend_col = "#13ff00"                                                                                           # Heat-exchanger temperature trend color
he_steel_pipe_sect_c = "#ff00d8"                                                                                        # Heat-exchanger steel pipe section lines color (hex)
he_glass_pipe_sect_c = "#8700ff"                                                                                        # Heat-exchanger glass pipe section lines color (hex)
he_sect_w = 2                                                                                                           # Heat-exchanger section lines width
he_sect_box_ypos_top_offs = 5                                                                                           # Heat-exchanger pipe section text-box Y-position offset from the top
he_sect_box_txt_size = 25                                                                                               # Heat-exchanger pipe section text-box text size
he_sect_box_alpha = 0.9                                                                                                 # Heat-exchanger pipe section text-box alpha

########
# DEFS #
########

# Data-acquisition plotting-mode enum definition
class Plt_mode(en.Enum):                                                                                                # Plottin' mode enum class
  complete = 1                                                                                                          # Complete plotting mode
  detailed = 2                                                                                                          # Detailed plotting mode
# Temperature trend enum definition
class Temp_trend(en.Enum):                                                                                              # Temperature trend enum class
  inlet = 0                                                                                                             # Inlet (bottom) temperature trend
  central = 1                                                                                                           # Central temperature trend
  outlet = 2                                                                                                            # Outlet (top) temperature trend

##########
# FUNCTS #
##########

# Function definition to customize plotting style
def set_plt_style(init_flg):                                                                                            # set_plt_style(Initialization flag)
  plt.style.use(plt_style);                                                                                             # Set defined plotting style
  sns.set(rc={"figure.figsize":(plt_size_x, plt_size_y)})                                                               # Plottin' size
  for param in ["figure.facecolor", "axes.facecolor", "savefig.facecolor"]:                                             # Chg bkg col
      plt.rcParams[param] = bkg_col;                                                                                    # Set col
  for param in ["text.color", "axes.labelcolor", "xtick.color", "ytick.color"]:                                         # Chg bkg txt
      plt.rcParams[param] = text_col;                                                                                   # Set col
  plt.grid(color=grid_col);                                                                                             # Set grid col
  if (not init_flg):                                                                                                    # If initialization flag ain't set
    leg = plt.legend();                                                                                                 # Mod legend
    plt.setp(leg.get_texts(), color=leg_txt_col);                                                                       # Set legend txt col
  return                                                                                                                # Return nothing

# Function to initialize personalized plotting style
def init_plt_style():                                                                                                   # init_plt_style()
  set_plt_style(True);                                                                                                  # Function call to set personalized plotting style with init flg
  plt.cla();                                                                                                            # Clear graph plotted axes
  plt.clf();                                                                                                            # Clear graph plotted to set personalized plotting style
  plt.close();                                                                                                          # Close plotted graph
  return                                                                                                                # Return nothing

# Function definition to graphically plot data filtering operations
def plot_data_flt(db, call_str, start_idxs_dbs, end_idxs_dbs, min_stddevs_dbs_idx, mode):                               # plot_data_flt(Datablock to split, Datablocks start idxs, Datablocks end idxs, Function call string, Plotting mode: complete/detailed)
  plt.title(call_str)                                                                                                   # Plot title
  plt.xlabel(plt_time_lbl)                                                                                              # X-axis lbl
  plt.ylabel(plt_temp_flow_lbl)                                                                                         # Y-axis lbl
  sns.lineplot(x=da.time_col, y=da.f1_col, data=db, label=plt_f1_lbl, color=plt_f1_col)                                 # Plot cold fluid vol flow rate (F1)
  sns.lineplot(x=da.time_col, y=da.f2_col, data=db, label=plt_f2_lbl, color=plt_f2_col)                                 # Plot hot fluid vol flow rate (F2)
  sns.lineplot(x=da.time_col, y=da.t1_col, data=db, label=plt_t1_lbl, color=plt_t1_col)                                 # Plot cold-in fluid temp (T1)
  sns.lineplot(x=da.time_col, y=da.t2_col, data=db, label=plt_t2_lbl, color=plt_t2_col)                                 # Plot hot-in fluid temp (T2)
  sns.lineplot(x=da.time_col, y=da.t3_col, data=db, label=plt_t3_lbl, color=plt_t3_col)                                 # Plot cold-out fluid temp (T3)
  sns.lineplot(x=da.time_col, y=da.t4_col, data=db, label=plt_t4_lbl, color=plt_t4_col)                                 # Plot hot-out fluid temp (T4)
  if (mode == Plt_mode.detailed):                                                                                       # In case of detailed plottin' mode selected
    plt.axvline(db[da.time_col].values[start_idxs_dbs[0]]-1, linewidth=dbs_w, color=dbs_c)                              # Plot first datablock-extra line
    for s_idx in start_idxs_dbs:                                                                                        # Datablocks staring-lines plottin' cycle
      if (s_idx == start_idxs_dbs[min_stddevs_dbs_idx]):                                                                # In case of selected datablock
        plt.axvline(db[da.time_col].values[s_idx], linewidth=dbs_sel_w, color=dbs_sel_c)                                # Plot selected datablock-line
      else:                                                                                                             # Else in case of unselected datablocks
        plt.axvline(db[da.time_col].values[s_idx], linewidth=dbs_w, color=dbs_c)                                        # Plot normal datablocks-lines
  idx = 0                                                                                                               # Detailed plotting mode list scrolling index
  tgt_idx = 0                                                                                                           # Detailed plotting mode list scrolling target index
  for e_idx in end_idxs_dbs:                                                                                            # Datablocks ending-lines plottin' cycle
    if (mode == Plt_mode.detailed and e_idx == end_idxs_dbs[min_stddevs_dbs_idx]):                                      # In case of selected datablock (and detailed plottin' mode selected)
      plt.axvline(db[da.time_col].values[e_idx-1], linewidth=dbs_sel_w, color=dbs_sel_c)                                # Plot selected datablock-line
      start_col = db[da.time_col].values[start_idxs_dbs[min_stddevs_dbs_idx]]                                           # Define fill-color starting point
      end_col = db[da.time_col].values[end_idxs_dbs[min_stddevs_dbs_idx]-1]                                             # Define fill-color ending point
      bottom, top = plt.ylim()                                                                                          # Get Y-axis limits to calc interval text-box pos
      plt.axvspan(start_col, end_col, alpha=dbs_sel_alpha, color=dbs_sel_col)                                           # Plot fill-color
      plt.text((start_col+end_col)/2, (bottom+top)/2+(((bottom+top)/2)/100)*sel_interval_txt_y_offs_p,
               sel_interval_lbl, fontsize=interval_txt_size, ha="center", va="center",
               bbox=dict(facecolor=interval_box_col, alpha=interval_box_alpha))                                         # Plot selected interval txt and box
    else:                                                                                                               # Else in case of unselected datablocks and/or complete plottin' mode selected
      plt.axvline(db[da.time_col].values[e_idx-1], linewidth=dbs_w, color=dbs_c)                                        # Plot normal datablocks-lines
      if (mode == Plt_mode.complete and idx == tgt_idx and idx < len(end_idxs_dbs)-1):                                  # Check fill-color plottin' cond
        start_col = db[da.time_col].values[end_idxs_dbs[idx]]                                                           # Define fill-color starting point
        end_col = db[da.time_col].values[end_idxs_dbs[idx+1]-1]                                                         # Define fill-color ending point
        bottom, top = plt.ylim()                                                                                        # Get Y-axis limits to calc interval text-box pos
        plt.axvspan(start_col, end_col, alpha=dbs_sel_alpha, color=dbs_sel_col)                                         # Plot fill-color
        plt.text((start_col+end_col)/2, (bottom+top)/2+(((bottom+top)/2)/100)*interval_txt_y_offs_p,
                 da.meas_str+str(int(idx/2)+1), fontsize=interval_txt_size, ha="center", va="center",
                 bbox=dict(facecolor=interval_box_col, alpha=interval_box_alpha))                                       # Plot interval txt and box
        tgt_idx += 2                                                                                                    # Upd detailed plotting mode list scrolling target index
    idx += 1                                                                                                            # Upd detailed plotting mode list scrolling index
  if (mode == Plt_mode.detailed):                                                                                       # In case of detailed plottin' mode selected
    plt.axvline(db[da.time_col].values[e_idx-1]+1, linewidth=dbs_w, color=dbs_c)                                        # Plot last datablock-extra line
  set_plt_style(False)                                                                                                  # Function call to set personalized plotting style without init flg
  plt.figure()                                                                                                          # Plot figure
  return                                                                                                                # Return nothing
  
# Function definition to graphically plot thermophysical variables interpolation/fitting
def plot_tp_vars(x_arr, y_arr, f_intp_fit, intp_fit_pts, mat_typ, x_lbl, y_lbl, y_intp_fit_lbl, res, tb_y_pos_offs):    # plot_tp_vars(X-array, Y-array, Interpolation/fitting function, Number of interpolation/fitting points to plot, Material type: Air at atm-pressure/water/AISI-316-stainless-steel/Borosilicate (pyrex) glass, X-label, Y-label, Y_interp_fitting-label, Poly-approximation result: discarded/accepted, Textbox Y-pos-offset)
  plt.title(materials[mat_typ.value]+y_lbl+plt_title_sep+materials[mat_typ.value]+y_intp_fit_lbl)                       # Plot title
  plt.xlabel(x_lbl)                                                                                                     # X-axis lbl
  plt.ylabel(materials[mat_typ.value]+y_lbl)                                                                            # Y-axis lbl
  sns.lineplot(x=x_arr, y=y_arr, label=materials[mat_typ.value]+y_lbl,
              marker=plt_marker, linestyle=plt_line, color=tp_var_col)                                                  # Plot thermophysical variable data
  x_inpt_fit_arr = np.linspace(x_arr[0], x_arr[-1], intp_fit_pts)                                                       # Generate interpolation/fitting X-array
  y_intp_fit_arr = f_intp_fit(x_inpt_fit_arr)                                                                           # Define interpolation/fitting Y-array
  sns.lineplot(x=x_inpt_fit_arr, y=y_intp_fit_arr,
               label=materials[mat_typ.value]+y_intp_fit_lbl, color=tp_intp_var_col)                                    # Plot thermophysical variable interpolated/fitted data
  left, right = plt.xlim()                                                                                              # Get X-axis limits to calc interval text-box pos
  bottom, top = plt.ylim()                                                                                              # Get Y-axis limits to calc interval text-box pos
  plt.text((left+right)/2, (bottom+top)/2+tb_y_pos_offs, poly_approx_res[res.value],
          fontsize=poly_approx_txt_size, ha="center", va="center",
          bbox=dict(facecolor=poly_approx_res_col[res.value], alpha=ploy_approx_box_alpha))                             # Plot poly-approximation result txt and box
  set_plt_style(False)                                                                                                  # Function call to set personalized plotting style without init flg
  plt.figure()                                                                                                          # Plot figure
  return                                                                                                                # Return nothing

# Function definition to graphically plot approssimative temperature trend in heat-exchanger section
def plot_he_sect_approx_temp_trend(meas_arr, he, section):                                                              # plot_he_sect_approx_temp_trend(Measures array, Heat exchanger, Heat-exchanger section: inlet(bottom)/outlet(top))
  for meas in meas_arr:                                                                                                 # Measures scrollin' cycle
    plt.title(meas.name+plt_he_temp_trend_titles[section.value])                                                        # Plot title
    plt.xlabel(plt_he_dist_lbl)                                                                                         # X-axis lbl
    plt.ylabel(plt_he_temp_lbl)                                                                                         # Y-axis lbl
    x_arr = [0.0, he.steel_pipes_id_m/2, he.steel_pipes_ed_m/2,
             ec.avg(he.steel_pipes_ed_m/2, (he.steel_pipes_ed_m/2)+he.steel_pipes_interspce),
             (he.steel_pipes_ed_m/2)+he.steel_pipes_interspce,
             (he.steel_pipes_ed_m/2)+he.steel_pipes_interspce+he.supp_glass_pipe_thick_m,
             (he.steel_pipes_ed_m/2)+(2*he.steel_pipes_interspce)+he.supp_glass_pipe_thick_m]                           # Distances array def [m]
    if (section == Temp_trend.inlet):                                                                                   # In case of inlet (bottom) section plotting req
      y_arr = [meas.t2, meas.int_pipes_int_surf_temp, meas.int_pipes_ext_surf_temp,
               meas.ext_in_sect_fl_temp, meas.ext_pipe_int_surf_temp, meas.ext_pipe_ext_surf_temp, he.env_temp]         # Inlet section temperatures array def [°C]
    elif (section == Temp_trend.central):                                                                               # Else in case of central section plotting req
      y_arr = [meas.avg_hot_fl_temp, meas.int_pipes_int_surf_temp, meas.int_pipes_ext_surf_temp,
               meas.avg_cold_fl_temp, meas.ext_pipe_int_surf_temp, meas.ext_pipe_ext_surf_temp, he.env_temp]            # Central section temperatures array def [°C]
    else:                                                                                                               # Else in case of outlet (top) section plotting req
      y_arr = [meas.t4, meas.int_pipes_int_surf_temp, meas.int_pipes_ext_surf_temp,
               meas.ext_out_sect_fl_temp, meas.ext_pipe_int_surf_temp, meas.ext_pipe_ext_surf_temp, he.env_temp]        # Outlet section temperatures array def [°C]
    plt.xlim(x_arr[0], x_arr[-1])                                                                                       # X-axis plotting limits mod
    plt.axvline(x_arr[1], linewidth=he_sect_w, color=he_steel_pipe_sect_c)                                              # Plot heat-exchanger section vertical-line (Steel pipe internal surface)
    plt.axvline(x_arr[2], linewidth=he_sect_w, color=he_steel_pipe_sect_c)                                              # Plot heat-exchanger section vertical-line (Steel pipe external surface)
    plt.axvline(x_arr[4], linewidth=he_sect_w, color=he_glass_pipe_sect_c)                                              # Plot heat-exchanger section vertical-line (Glass pipe internal surface)
    plt.axvline(x_arr[5], linewidth=he_sect_w, color=he_glass_pipe_sect_c)                                              # Plot heat-exchanger section vertical-line (Glass pipe supposed external surface)
    sns.lineplot(x=x_arr, y=y_arr, label=he_temp_trend_lbls[section.value],
                 color=he_temp_trend_col, marker=plt_marker, linestyle=plt_line)                                        # Plot heat-exchanger section temperature trend interpolated data
    _, top = plt.ylim()                                                                                                 # Get Y-axis limits to calc text-boxes pos
    plt.text(ec.avg(x_arr[1], x_arr[2]), top-he_sect_box_ypos_top_offs,
             he_steel_pipe_lbl, fontsize=he_sect_box_txt_size, ha="center", va="center",
             bbox = dict(facecolor=he_steel_pipe_sect_c, alpha=he_sect_box_alpha))                                      # Plot steel pipe profile txt and box
    plt.text(ec.avg(x_arr[4], x_arr[5]), top-he_sect_box_ypos_top_offs,
             he_glass_pipe_lbl, fontsize=he_sect_box_txt_size, ha="center", va="center",
             bbox = dict(facecolor=he_glass_pipe_sect_c, alpha=he_sect_box_alpha))                                      # Plot glass pipe profile txt and box
    set_plt_style(False)                                                                                                # Function call to set personalized plotting style without init flg
    plt.figure()                                                                                                        # Plot figure
  return                                                                                                                # Return nothing
