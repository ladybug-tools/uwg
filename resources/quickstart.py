from UWG import UWG

# Define the .epw, .uwg filenames to create an UWG object.
# UWG will look for the .epw file int the UWG/resources/epw folder,
# and the .uwg file in the UWG/resources/parameters folder.
epw_filename = "SGP_Singapore.486980_IWEC.epw"      # EPW file name
param_filename = "initialize_singapore.uwg"         # .uwg file name

# Initialize the UWG object and run the simulation
uwg = UWG(epw_filename, param_filename)
uwg.run()
