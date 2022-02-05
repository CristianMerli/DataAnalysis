########
# LIBS #
########

# Libraries import
import enum as en                                                                                                       # Enum lib
import numpy as np                                                                                                      # Data-analysis numpy lib
import math as mt                                                                                                       # Math lib
import numpy.polynomial.polynomial as poly                                                                              # Polynomials management lib
from scipy.interpolate import interp1d                                                                                  # 1D-interpolation lib

########
# DEFS #
########

# Polynomial approximation mode enum definition
class Poly_approx_md(en.Enum):                                                                                          # Poly approx mode enum class
  int = 1                                                                                                               # Poly approx with interpolation
  fit = 2                                                                                                               # Poly approx with curve-fitting

##########
# FUNCTS #
##########

# Function definition to clean data invalid values for polynomial approximation
def cln_invalid_data(x_arr, y_arr):                                                                                     # cln_invalid_data(X-array, Y-array)
  x_cln_arr = list()                                                                                                    # Define empy cleaned X-array list
  y_cln_arr = list()                                                                                                    # Define empy cleaned Y-array list
  for idx in range(len(x_arr)):                                                                                         # X and Y arrays elements scrollin' cycle
    if (not mt.isnan(y_arr[idx])):                                                                                      # Valid element detectin' cond
      x_cln_arr.append(x_arr[idx])                                                                                      # Add sel element from X-array into cleaned X-array list
      y_cln_arr.append(y_arr[idx])                                                                                      # Add sel element from Y-array into cleaned Y-array list
  return np.array(x_cln_arr), np.array(y_cln_arr)                                                                       # Return X and Y cleaned arrays lists converted into num-arrays

# Function definition to apply selected polynomial approximation (interpolation/cruve-fitting)
# and return approximating function
def poly_approx(x_arr, y_arr, approx_mode, typ_ord):                                                                    # poly_approx(X-array, Y-array, Approx-mode: poly-interpolation/poly-fitting, Poly-interpolation-type/Fitting-poly-order)
  x_cln_arr, y_cln_arr = cln_invalid_data(x_arr, y_arr)                                                                 # Function call to get X and Y cleand arrays deletin' invalid vals
  if (approx_mode == Poly_approx_md.int):                                                                               # In case of poly-interpolation approx mode sel
    return interp1d(x_cln_arr, y_cln_arr, kind=typ_ord), x_cln_arr, y_cln_arr                                           # Return polynomial interpolation function and X/Y cleand arrays
  else:                                                                                                                 # Else in case of poly-fitting approx mode sel
    return poly.Polynomial.fit(x_cln_arr, y_cln_arr, typ_ord), x_cln_arr, y_cln_arr                                     # Return polynomial fitting function and X/Y cleand arrays
