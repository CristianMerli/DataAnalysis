########
# LIBS #
########

# Libraries import
import pandas as pd                                                                                                     # Data-analysis panda lib
import numpy as np                                                                                                      # Data-analysis numpy lib
import enum as en                                                                                                       # Enum lib
# Project personal libraries import
import libs.poly_approx_lib as pa                                                                                       # Poly-approximation lib
import libs.plotting_lib as pl                                                                                          # Plotting lib
import libs.eng_calcs_lib as ec                                                                                         # Engineering calcs lib

########
# VARS #
########

# Table files vars: (.csv) table files
csv_air_atmp_filepath = "../thermophys_vars/air_atmp_thermo_vars.csv"                                                   # Air at atm-pressure thermophysical vars data filepath (.csv table file)
csv_water_filepath = "../thermophys_vars/water_thermo_vars.csv"                                                         # Water thermophysical vars data filepath (.csv table file)
csv_aisi_316_filepath = "../thermophys_vars/aisi_316_thermo_vars.csv"                                                   # AISI-316 stainless-steel thermophysical vars data filepath (.csv table file)
# Table files structure vars
csv_sep_chr = ';'                                                                                                       # Sep chr in (.csv) table files
# Table files columns vars
temp_col = "T(degC)"                                                                                                    # Temp col in table files (Temperature [°C])
temp_col_k = "T(K)"                                                                                                     # Temp col in table files (Temperature [K])
rho_col = "Rho(kg/m3)"                                                                                                  # Rho col in table files (Density [kg/m^3])
cp_col = "Cp(kJ/kg*K)"                                                                                                  # Cp col in table files (Specific heat at constant pressure [kJ/(kg*K)])
lambda_col = "Lambda(W/m*K)"                                                                                            # Lambda col in table files (Thermal conductivity [W/(m*K)])
ni_col = "Ni(m2/s)"                                                                                                     # Ni col in table files (Kinematic viscosity [m^2/s])
beta_col = "Beta(1/K)"                                                                                                  # Beta col in table files (Thermal expansion coefficient [1/K])
pr_col = "Pr(-)"                                                                                                        # Pr col in table files (Prandtl number [adimensional])
# Air atmp-pressure thermophysical variables vs temp interpolation vars
air_atmp_intp_typ = "cubic"                                                                                             # Air thermophysical variables at atm pressure interpolation type
air_atmp_intp_fit_pts = 1000                                                                                            # Air thermophysical variables at atm pressure number of interpolation/fitting plotting points
# Water atmp-pressure thermophysical variables vs temp interpolation vars
water_intp_typ = "cubic"                                                                                                # Water thermophysical variables interpolation type
water_intp_fit_pts = 1000                                                                                               # Water thermophysical variables number of interpolation/fitting plotting points
# AISI-316 stainless-steel thermophysical variables vs temp interpolation vars
aisi_316_intp_typ = "cubic"                                                                                             # AISI-316 stainless-steel thermophysical variables interpolation type
aisi_316_intp_fit_pts = 1000                                                                                            # AISI-316 stainless-steel thermophysical variables number of interpolation/fitting plotting points

########
# DEFS #
########

# Thermophysical-variables polynomial approximation discarded/accepted result enum definition to plot graph text-box
class Poly_approx_res(en.Enum):                                                                                         # Poly approx result enum class
  disc = 0                                                                                                              # Poly approx discarded
  acc = 1                                                                                                               # Poly approx accepted
# Thermophysical-variables material type text plotting enum definition
class Mat_type(en.Enum):                                                                                                # Material type text enum class
  air_atmp = 0                                                                                                          # Air at atmospheric pressure
  water = 1                                                                                                             # Water
  aisi_316 = 2                                                                                                          # AISI-316-stainless-steel

##########
# FUNCTS #
##########

# Function definition load thermophysical vars data from (.csv) files and return DataFrames vars
def load_thermophys_vars_data():                                                                                        # load_thermophys_vars_data()
  air_atmp = pd.read_csv(csv_air_atmp_filepath, sep=csv_sep_chr, encoding="utf8")                                       # Import data from (.csv) table file and create a new panda DataFrame variable: 'air_atmp'
  water = pd.read_csv(csv_water_filepath, sep=csv_sep_chr, encoding="utf8")                                             # Import data from (.csv) table file and create a new panda DataFrame variable: 'water'
  aisi_316 = pd.read_csv(csv_aisi_316_filepath, sep=csv_sep_chr, encoding="utf8")                                       # Import data from (.csv) table file and create a new panda DataFrame variable: 'aisi_316'
  return air_atmp, water, aisi_316                                                                                      # Return 'air_atmp', 'water' and 'aisi_316' DataFrames vars containing the respective thermophysical variables

# Function definition to apply polynomial approximation and plot
# air's thermophysical variables vs temperature (at atmospheric pressure)
def poly_approx_plot_air_atmp_thermophys_vars(air_atmp, plt_flg):                                                       # poly_approx_plot_thermophys_vars_air_atmp(Air at atm-pressure DataFrame var, Plotting flag)
  air_atmp_temp = np.array(air_atmp[temp_col])                                                                          # Extract air interpolation/fitting temperatures array (at atm pressure) from DataFrame
  air_atmp_rho = np.array(air_atmp[rho_col])                                                                            # Extract air density array from DataFrame (thermophysic variable at atm pressure vs temp)
  air_atmp_cp = np.array(air_atmp[cp_col])                                                                              # Extract air specific heat at constant pressure array from DataFrame (thermophysic variable at atm pressure vs temp)
  air_atmp_lambda = np.array(air_atmp[lambda_col])                                                                      # Extract air thermal conductivity array from DataFrame (thermophysic variable at atm pressure vs temp)
  air_atmp_ni = np.array(air_atmp[ni_col])                                                                              # Extract air kinematic viscosity array from DataFrame (thermophysic variable at atm pressure vs temp)
  air_atmp_beta = np.array(air_atmp[beta_col])                                                                          # Extract air thermal expansion coefficient array from DataFrame (thermophysic variable at atm pressure vs temp)
  air_atmp_pr = np.array(air_atmp[pr_col])                                                                              # Extract air Prandtl number array from DataFrame(thermophysic variable at atm pressure vs temp)
  # Poly-fitting discarded, bad approximation
  f_air_atmp_rho, x, y = pa.poly_approx(air_atmp_temp, air_atmp_rho, pa.Poly_approx_md.fit, 4)                          # Air density fitting function (thermophysic variable at atm pressure vs temp)
  if (plt_flg):                                                                                                         # If plotting flag is ena
    pl.plot_tp_vars(x, y, f_air_atmp_rho, air_atmp_intp_fit_pts, Mat_type.air_atmp, pl.plt_temp_lbl,
                    pl.plt_rho_lbl, pl.plt_fit_rho_lbl, Poly_approx_res.disc, 0.6)                                      # Function call to plot air density and curve-fitting (thermophysic variable at atm pressure vs temp)
  # Poly-interpolation accepted, good approximation
  f_air_atmp_rho, x, y = pa.poly_approx(air_atmp_temp, air_atmp_rho, pa.Poly_approx_md.int, air_atmp_intp_typ)          # Air density interpolation function (thermophysic variable at atm pressure vs temp)
  if (plt_flg):                                                                                                         # If plotting flag is ena
    pl.plot_tp_vars(x, y, f_air_atmp_rho, air_atmp_intp_fit_pts, Mat_type.air_atmp, pl.plt_temp_lbl,
                    pl.plt_rho_lbl, pl.plt_intp_rho_lbl, Poly_approx_res.acc, 0.6)                                      # Function call to plot air density and interpolation (thermophysic variable at atm pressure vs temp)
  # Poly-interpolation discarded, bad approximation
  f_air_atmp_cp, x, y = pa.poly_approx(air_atmp_temp, air_atmp_cp, pa.Poly_approx_md.int, air_atmp_intp_typ)            # Air specific heat at constant pressure interpolation function (thermophysic variable at atm pressure vs temp)
  if (plt_flg):                                                                                                         # If plotting flag is ena
    pl.plot_tp_vars(x, y, f_air_atmp_cp, air_atmp_intp_fit_pts, Mat_type.air_atmp, pl.plt_temp_lbl,
                    pl.plt_cp_lbl, pl.plt_intp_cp_lbl, Poly_approx_res.disc, 1.12*1e-2)                                 # Function call to plot air specific heat at constant pressure and interpolation (thermophysic variable at atm pressure vs temp)
  # Poly-fitting accepted, good approximation
  f_air_atmp_cp, x, y = pa.poly_approx(air_atmp_temp, air_atmp_cp, pa.Poly_approx_md.fit, 4)                            # Air specific heat at constant pressure fitting function (thermophysic variable at atm pressure vs temp)
  if (plt_flg):                                                                                                         # If plotting flag is ena
    pl.plot_tp_vars(x, y, f_air_atmp_cp, air_atmp_intp_fit_pts, Mat_type.air_atmp, pl.plt_temp_lbl,
                    pl.plt_cp_lbl, pl.plt_fit_cp_lbl, Poly_approx_res.acc, 1.12*1e-2)                                   # Function call to plot air specific heat at constant pressure and curve-fitting (thermophysic variable at atm pressure vs temp)
  # Poly-interpolation discarded, bad approximation
  f_air_atmp_lambda, x, y = pa.poly_approx(air_atmp_temp, air_atmp_lambda, pa.Poly_approx_md.int, air_atmp_intp_typ)    # Air thermal conductivity interpolation function (thermophysic variable at atm pressure vs temp)
  if (plt_flg):                                                                                                         # If plotting flag is ena
    pl.plot_tp_vars(x, y, f_air_atmp_lambda, air_atmp_intp_fit_pts, Mat_type.air_atmp, pl.plt_temp_lbl,
                    pl.plt_lambda_lbl, pl.plt_intp_lambda_lbl, Poly_approx_res.disc, 1.17*1e-2)                         # Function call to plot air thermal conductivity and interpolation (thermophysic variable at atm pressure vs temp)
  # Poly-fitting accepted, good approximation
  f_air_atmp_lambda, x, y = pa.poly_approx(air_atmp_temp, air_atmp_lambda, pa.Poly_approx_md.fit, 4)                    # Air thermal conductivity fitting function (thermophysic variable at atm pressure vs temp)
  if (plt_flg):                                                                                                         # If plotting flag is ena
    pl.plot_tp_vars(x, y, f_air_atmp_lambda, air_atmp_intp_fit_pts, Mat_type.air_atmp, pl.plt_temp_lbl,
                    pl.plt_lambda_lbl, pl.plt_fit_lambda_lbl, Poly_approx_res.acc, 1.17*1e-2)                           # Function call to plot air thermal conductivity and curve-fitting (thermophysic variable at atm pressure vs temp)
  # Poly-interpolation discarded, bad approximation
  f_air_atmp_ni, x, y = pa.poly_approx(air_atmp_temp, air_atmp_ni, pa.Poly_approx_md.int, air_atmp_intp_typ)            # Air kinematic viscosity interpolation function (thermophysic variable at atm pressure vs temp)
  if (plt_flg):                                                                                                         # If plotting flag is ena
    pl.plot_tp_vars(x, y, f_air_atmp_ni, air_atmp_intp_fit_pts, Mat_type.air_atmp, pl.plt_temp_lbl,
                    pl.plt_ni_lbl, pl.plt_intp_ni_lbl, Poly_approx_res.disc, 1.38*1e-5)                                 # Function call to plot air kinematic viscosity and interpolation (thermophysic variable at atm pressure vs temp)
  # Poly-fitting accepted, good approximation
  f_air_atmp_ni, x, y = pa.poly_approx(air_atmp_temp, air_atmp_ni, pa.Poly_approx_md.fit, 4)                            # Air kinematic viscosity fitting function (thermophysic variable at atm pressure vs temp)
  if (plt_flg):                                                                                                         # If plotting flag is ena
    pl.plot_tp_vars(x, y, f_air_atmp_ni, air_atmp_intp_fit_pts, Mat_type.air_atmp, pl.plt_temp_lbl,
                    pl.plt_ni_lbl, pl.plt_fit_ni_lbl, Poly_approx_res.acc, 1.38*1e-5)                                   # Function call to plot air kinematic viscosity and curve-fitting (thermophysic variable at atm pressure vs temp)
  # Poly-fitting discarded, bad approximation
  f_air_atmp_beta, x, y = pa.poly_approx(air_atmp_temp, air_atmp_beta, pa.Poly_approx_md.fit, 4)                        # Air thermal expansion coefficient fitting function (thermophysic variable at atm pressure vs temp)
  if (plt_flg):                                                                                                         # If plotting flag is ena
    pl.plot_tp_vars(x, y, f_air_atmp_beta, air_atmp_intp_fit_pts, Mat_type.air_atmp, pl.plt_temp_lbl,
                    pl.plt_beta_lbl, pl.plt_fit_beta_lbl, Poly_approx_res.disc, 1.65*1e-3)                              # Function call to plot air thermal expansion coefficient and curve-fitting (thermophysic variable at atm pressure vs temp)
  # Poly-interpolation accepted, good approximation
  f_air_atmp_beta, x, y = pa.poly_approx(air_atmp_temp, air_atmp_beta, pa.Poly_approx_md.int, air_atmp_intp_typ)        # Air thermal expansion coefficient interpolation function (thermophysic variable at atm pressure vs temp)
  if (plt_flg):                                                                                                         # If plotting flag is ena
    pl.plot_tp_vars(x, y, f_air_atmp_beta, air_atmp_intp_fit_pts, Mat_type.air_atmp, pl.plt_temp_lbl,
                    pl.plt_beta_lbl, pl.plt_intp_beta_lbl, Poly_approx_res.acc, 1.65*1e-3)                              # Function call to plot air thermal expansion coefficient and interpolation (thermophysic variable at atm pressure vs temp)
  # Poly-interpolation discarded, bad approximation
  f_air_atmp_pr, x, y = pa.poly_approx(air_atmp_temp, air_atmp_pr, pa.Poly_approx_md.int, air_atmp_intp_typ)            # Air Prandtl number interpolation function (thermophysic variable at atm pressure vs temp)
  if (plt_flg):                                                                                                         # If plotting flag is ena
    pl.plot_tp_vars(x, y, f_air_atmp_pr, air_atmp_intp_fit_pts, Mat_type.air_atmp, pl.plt_temp_lbl,
                    pl.plt_pr_lbl, pl.plt_intp_pr_lbl, Poly_approx_res.disc, 1.6*1e-2)                                  # Function call to plot air Prandtl number and interpolation (thermophysic variable at atm pressure vs temp)
  # Poly-fitting accepted, good approximation
  f_air_atmp_pr, x, y = pa.poly_approx(air_atmp_temp, air_atmp_pr, pa.Poly_approx_md.fit, 7)                            # Air Prandtl number fitting function (thermophysic variable at atm pressure vs temp)
  if (plt_flg):                                                                                                         # If plotting flag is ena
    pl.plot_tp_vars(x, y, f_air_atmp_pr, air_atmp_intp_fit_pts, Mat_type.air_atmp, pl.plt_temp_lbl,
                    pl.plt_pr_lbl, pl.plt_fit_pr_lbl, Poly_approx_res.acc, 1.6*1e-2)                                    # Function call to plot air Prandtl number and curve-fitting (thermophysic variable at atm pressure vs temp)
  return f_air_atmp_rho, f_air_atmp_cp, f_air_atmp_lambda, f_air_atmp_ni, f_air_atmp_beta, f_air_atmp_pr                # Return air's thermophysical vars poly-approx functions (at atm-pressure)

# Function definition to apply polynomial approximation and plot water's thermophysical variables vs temperature
def poly_approx_plot_water_thermophys_vars(water, plt_flg):                                                             # poly_approx_plot_water_thermophys_vars(Water DataFrame var, Plotting flag)
  water_temp = np.array(water[temp_col])                                                                                # Extract water interpolation/fitting temperatures array from DataFrame
  water_rho = np.array(water[rho_col])                                                                                  # Extract water density array from DataFrame (thermophysic variable at atm pressure vs temp)
  water_cp = np.array(water[cp_col])                                                                                    # Extract water specific heat at constant pressure array from DataFrame (thermophysic variable vs temp)
  water_lambda = np.array(water[lambda_col])                                                                            # Extract water thermal conductivity array from DataFrame (thermophysic variable vs temp)
  water_ni = np.array(water[ni_col])                                                                                    # Extract water kinematic viscosity array from DataFrame (thermophysic variable vs temp)
  water_beta = np.array(water[beta_col])                                                                                # Extract water thermal expansion coefficient array from DataFrame (thermophysic variable vs temp)
  water_pr = np.array(water[pr_col])                                                                                    # Extract water Prandtl number array from DataFrame (thermophysic variable vs temp)
  # Poly-interpolation discarded, bad approximation
  f_water_rho, x, y = pa.poly_approx(water_temp, water_rho, pa.Poly_approx_md.int, water_intp_typ)                      # Water density interpolation function (thermophysic variable at atm pressure vs temp)
  if (plt_flg):                                                                                                         # If plotting flag is ena
    pl.plot_tp_vars(x, y, f_water_rho, water_intp_fit_pts, Mat_type.water, pl.plt_temp_lbl,
                    pl.plt_rho_lbl, pl.plt_intp_rho_lbl, Poly_approx_res.disc, -8.2)                                    # Function call to plot water density and interpolation (thermophysic variable at atm pressure vs temp)
  # Poly-fitting accepted, good approximation
  f_water_rho, x, y = pa.poly_approx(water_temp, water_rho, pa.Poly_approx_md.fit, 4)                                   # Water density fitting function (thermophysic variable at atm pressure vs temp)
  if (plt_flg):                                                                                                         # If plotting flag is ena
    pl.plot_tp_vars(x, y, f_water_rho, water_intp_fit_pts, Mat_type.water, pl.plt_temp_lbl,
                    pl.plt_rho_lbl, pl.plt_fit_rho_lbl, Poly_approx_res.acc, -8.2)                                      # Function call to plot water density and curve-fitting (thermophysic variable at atm pressure vs temp)
  # Poly-interpolation discarded, bad approximation
  f_water_cp, x, y = pa.poly_approx(water_temp, water_cp, pa.Poly_approx_md.int, water_intp_typ)                        # Water specific heat at constant pressure interpolation function (thermophysic variable at atm pressure vs temp)
  if (plt_flg):                                                                                                         # If plotting flag is ena
    pl.plot_tp_vars(x, y, f_water_cp, water_intp_fit_pts, Mat_type.water, pl.plt_temp_lbl,
                    pl.plt_cp_lbl, pl.plt_intp_cp_lbl, Poly_approx_res.disc, 2.45*1e-3)                                 # Function call to plot water specific heat at constant pressure and interpolation (thermophysic variable at atm pressure vs temp)
  # Poly-fitting accepted, good approximation
  f_water_cp, x, y = pa.poly_approx(water_temp, water_cp, pa.Poly_approx_md.fit, 4)                                     # Water specific heat at constant pressure fitting function (thermophysic variable at atm pressure vs temp)
  if (plt_flg):                                                                                                         # If plotting flag is ena
    pl.plot_tp_vars(x, y, f_water_cp, water_intp_fit_pts, Mat_type.water, pl.plt_temp_lbl,
                    pl.plt_cp_lbl, pl.plt_fit_cp_lbl, Poly_approx_res.acc, 2.45*1e-3)                                   # Function call to plot water specific heat at constant pressure and curve-fitting (thermophysic variable at atm pressure vs temp)
  # Poly-interpolation discarded, bad approximation
  f_water_lambda, x, y = pa.poly_approx(water_temp, water_lambda, pa.Poly_approx_md.int, water_intp_typ)                # Water thermal conductivity interpolation function (thermophysic variable at atm pressure vs temp)
  if (plt_flg):                                                                                                         # If plotting flag is ena
    pl.plot_tp_vars(x, y, f_water_lambda, water_intp_fit_pts, Mat_type.water, pl.plt_temp_lbl,
                    pl.plt_lambda_lbl, pl.plt_intp_lambda_lbl, Poly_approx_res.disc, -2.1*1e-2)                         # Function call to plot water thermal conductivity and interpolation (thermophysic variable at atm pressure vs temp)
  # Poly-fitting accepted, good approximation
  f_water_lambda, x, y = pa.poly_approx(water_temp, water_lambda, pa.Poly_approx_md.fit, 4)                             # Water thermal conductivity fitting function (thermophysic variable at atm pressure vs temp)
  if (plt_flg):                                                                                                         # If plotting flag is ena
    pl.plot_tp_vars(x, y, f_water_lambda, water_intp_fit_pts, Mat_type.water, pl.plt_temp_lbl,
                    pl.plt_lambda_lbl, pl.plt_fit_lambda_lbl, Poly_approx_res.acc, -2.1*1e-2)                           # Function call to plot water thermal conductivity and curve-fitting (thermophysic variable at atm pressure vs temp)
  # Poly-fitting discarded, bad approximation
  f_water_ni, x, y = pa.poly_approx(water_temp, water_ni, pa.Poly_approx_md.fit, 4)                                     # Water kinematic viscosity fitting function (thermophysic variable at atm pressure vs temp)
  if (plt_flg):                                                                                                         # If plotting flag is ena
    pl.plot_tp_vars(x, y, f_water_ni, water_intp_fit_pts, Mat_type.water, pl.plt_temp_lbl,
                    pl.plt_ni_lbl, pl.plt_fit_ni_lbl, Poly_approx_res.disc, 2.65*1e-7)                                  # Function call to plot water kinematic viscosity and curve-fitting (thermophysic variable at atm pressure vs temp)
  # Poly-interpolation accepted, good approximation
  f_water_ni, x, y = pa.poly_approx(water_temp, water_ni, pa.Poly_approx_md.int, water_intp_typ)                        # Water kinematic viscosity interpolation function (thermophysic variable at atm pressure vs temp)
  if (plt_flg):                                                                                                         # If plotting flag is ena
    pl.plot_tp_vars(x, y, f_water_ni, water_intp_fit_pts, Mat_type.water, pl.plt_temp_lbl,
                    pl.plt_ni_lbl, pl.plt_intp_ni_lbl, Poly_approx_res.acc, 2.65*1e-7)                                  # Function call to plot water kinematic viscosity and interpolation (thermophysic variable at atm pressure vs temp)
  # Poly-interpolation discarded, bad approximation
  f_water_beta, x, y = pa.poly_approx(water_temp, water_beta, pa.Poly_approx_md.int, water_intp_typ)                    # Water thermal expansion coefficient interpolation function (thermophysic variable at atm pressure vs temp)
  if (plt_flg):                                                                                                         # If plotting flag is ena
    pl.plot_tp_vars(x, y, f_water_beta, water_intp_fit_pts, Mat_type.water, pl.plt_temp_lbl,
                    pl.plt_beta_lbl, pl.plt_intp_beta_lbl, Poly_approx_res.disc, -1.55*1e-4)                            # Function call to plot water thermal expansion coefficient and interpolation (thermophysic variable at atm pressure vs temp)
  # Poly-fitting accepted, good approximation
  f_water_beta, x, y = pa.poly_approx(water_temp, water_beta, pa.Poly_approx_md.fit, 2)                                 # Water thermal expansion coefficient fitting function (thermophysic variable at atm pressure vs temp)
  if (plt_flg):                                                                                                         # If plotting flag is ena
    pl.plot_tp_vars(x, y, f_water_beta, water_intp_fit_pts, Mat_type.water, pl.plt_temp_lbl,
                    pl.plt_beta_lbl, pl.plt_fit_beta_lbl, Poly_approx_res.acc, -1.55*1e-4)                              # Function call to plot water thermal expansion coefficient and curve-fitting (thermophysic variable at atm pressure vs temp)
  # Poly-fitting discarded, bad approximation
  f_water_pr, x, y = pa.poly_approx(water_temp, water_pr, pa.Poly_approx_md.fit, 4)                                     # Water Prandtl number fitting function (thermophysic variable at atm pressure vs temp)
  if (plt_flg):                                                                                                         # If plotting flag is ena
    pl.plot_tp_vars(x, y, f_water_pr, water_intp_fit_pts, Mat_type.water, pl.plt_temp_lbl,
                    pl.plt_pr_lbl, pl.plt_fit_pr_lbl, Poly_approx_res.disc, 2.75)                                       # Function call to plot water Prandtl number and curve-fitting (thermophysic variable at atm pressure vs temp)
  # Poly-interpolation accepted, good approximation
  f_water_pr, x, y = pa.poly_approx(water_temp, water_pr, pa.Poly_approx_md.int, water_intp_typ)                        # Water Prandtl number interpolation function (thermophysic variable at atm pressure vs temp)
  if (plt_flg):                                                                                                         # If plotting flag is ena
    pl.plot_tp_vars(x, y, f_water_pr, water_intp_fit_pts, Mat_type.water, pl.plt_temp_lbl,
                    pl.plt_pr_lbl, pl.plt_intp_pr_lbl, Poly_approx_res.acc, 2.75)                                       # Function call to plot water Prandtl number and interpolation (thermophysic variable at atm pressure vs temp)    
  return f_water_rho, f_water_cp, f_water_lambda, f_water_ni, f_water_beta, f_water_pr                                  # Return water's thermophysical vars poly-approx functions

# Function definition to apply polynomial approximation and plot AISI-316-stainless-steel's
# thermophysical variables vs temperature
def poly_approx_plot_aisi_316_thermophys_vars(aisi_316, plt_flg):                                                       # poly_approx_plot_aisi_316_thermophys_vars(AISI-316 DataFrame var, Plotting flag)
  aisi_316_temp = ec.conv_temp_k_c(np.array(aisi_316[temp_col_k]))                                                      # Extract AISI-316 interpolation/fitting temperatures array from DataFrame and convert it from [K] into [°C]
  aisi_316_lambda = np.array(aisi_316[lambda_col])                                                                      # Extract AISI-316 thermal conductivity array from DataFrame (thermophysic variable vs temp)
  # Poly-interpolation discarded, bad approximation
  f_aisi_316_lambda, x, y = pa.poly_approx(aisi_316_temp, aisi_316_lambda, pa.Poly_approx_md.int, aisi_316_intp_typ)    # AISI-316 thermal conductivity interpolation function (thermophysic variable at atm pressure vs temp)
  if (plt_flg):                                                                                                         # If plotting flag is ena
    pl.plot_tp_vars(x, y, f_aisi_316_lambda, aisi_316_intp_fit_pts, Mat_type.aisi_316, pl.plt_temp_lbl,
                    pl.plt_lambda_lbl, pl.plt_intp_lambda_lbl, Poly_approx_res.disc, -2.45)                             # Function call to plot AISI-316 thermal conductivity and interpolation (thermophysic variable at atm pressure vs temp)
  # Poly-fitting accepted, good approximation
  f_aisi_316_lambda, x, y = pa.poly_approx(aisi_316_temp, aisi_316_lambda, pa.Poly_approx_md.fit, 4)                    # AISI-316 thermal conductivity fitting function (thermophysic variable at atm pressure vs temp)
  if (plt_flg):                                                                                                         # If plotting flag is ena
    pl.plot_tp_vars(x, y, f_aisi_316_lambda, aisi_316_intp_fit_pts, Mat_type.aisi_316, pl.plt_temp_lbl,
                    pl.plt_lambda_lbl, pl.plt_fit_lambda_lbl, Poly_approx_res.acc, -2.45)                               # Function call to plot AISI-316 thermal conductivity and curve-fitting (thermophysic variable at atm pressure vs temp)
  return f_aisi_316_lambda                                                                                              # Return AISI-316-stainless-steel's thermophysical vars poly-approx functions
