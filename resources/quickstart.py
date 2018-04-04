import UWG
import os

# Gets path of current directory
CURR_DIRECTORY = os.path.abspath(os.path.dirname(__file__))

# To run UWG provide the following inputs
epw_directory = os.path.join(CURR_DIRECTORY,"epw")  # EPW file directory
epw_filename = "SGP_Singapore.486980_IWEC.epw"      # EPW file name
uwg_param_directory = CURR_DIRECTORY                # .uwg file directory
uwg_param_filename = "initialize.uwg"               # .uwg file name

# Initialize the UWG object
uwg = UWG.UWG(epw_directory, epw_filename, uwg_param_directory, uwg_param_filename)

# Run the simulation
uwg.run()
