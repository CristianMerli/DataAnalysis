########
# LIBS #
########

# Libraries import
import enum as en                                                                                                       # Enum lib

########
# VARS #
########

# General output files vars
out_dir = "../../final_doc/code_exports/output/"                                                                        # Output files directory path def
out_file_ext = ".txt"                                                                                                   # Output files extension def
# Filenames vars
filenames = ["delim_idxs", "std_devs", "measures", "he", "measures_approx_calcs", "measures_calcs"]                     # Filenames strings array def

########
# DEFS #
########

# Output enum definition
class Output_typ(en.Enum):                                                                                              # Output type enum class
  delim_idxs = 0                                                                                                        # Delimiting indexes array output
  std_devs = 1                                                                                                          # Standard-deviations output
  measures = 2                                                                                                          # Measure-values output
  he = 3                                                                                                                # Heat-exchanger output
  measures_approx_calcs = 4                                                                                             # Measures calculated approximative values output
  measures_calcs = 5                                                                                                    # Measures calculated values output

##########
# FUNCTS #
##########

# Function definition to save text into an output file
def save_output(out_typ, txt):                                                                                          # save_output(Output save type, Text to save)
  with open(out_dir+filenames[out_typ.value]+out_file_ext, 'w') as out_file:                                            # Open output file in write mode
    out_file.write(txt)                                                                                                 # Write text into output file
  return                                                                                                                # Return nothing
