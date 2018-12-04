from uwg import uwg

# Define the .epw, .uwg filenames to create an UWG object.
# UWG will look for the .epw file in the UWG/resources/epw folder,
# and the .uwg file in the UWG/resources/parameters folder.
#epw_filename = "SGP_Singapore.486980_IWEC.epw"      # EPW file name
epw_filename = "ISR_Tel.Aviv-Bet.Dagan.401790_MSI.epw"
param_filename = "initialize_israel.uwg"         # .uwg file name

# Initialize the UWG object and run the simulation
uwg_ = uwg(epw_filename, param_filename)
uwg_.run()
