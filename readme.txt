###################### 1 #####################################
############# INPUT GENERATION/PREPARATION ###################

+-- monte_carlo
|   +-- inputs
|   +-- input_preparation.py
|   +-- dists.py

'input_preparation.py' generates a file (like 'MCsim_WTC.xls')
in the folder 'inputs' with data sampled form the inputs
distributions specified in the file 'dists.py'
##############################################################


#################### 2 #######################################
############### FLIGHTS SIMULATION (TRUE FLIGHTS) ############

+-- monte carlo
|   +-- descent
|   +-- flight
|   +-- fpls
|   +-- inputs
|   +-- mc_true_flights_sim.py
|   +-- simulate_true_descent_backwards.py
|   +-- simulate_true_flight.py
|   +-- res_log_WTC.txt

'mc_true_flights_sim.py' simulates the true flights for each WTC using
the inputs from the files like 'MCsim_WTC.xls' in the 'input'
folder. The generated flights (and descents) are stored
in the folders 'flight' (and 'descent').

Furthermore, the ids of the flights that were not simulated
successfully are recorded in the files like 'res_log_WTC.txt'.
These files are discarded from the further analysis.

###############################################################


################## 3 ##########################################
########### FILES SORTING AND CLEANING ########################

+-- cleaning files
|   +-- arrange_files_into_folders.py
|   +-- rename_failed_flights.py

'arrange_files_into_folders.py' moves the generated
true flights into the WTC folder according to the route number.

'rename_failed_flights.py' renames the failed flights
in the WTC/route folders based on the ids specified in the
'res_log_WTC.txt' files.

################################################################


################ 4 #############################################
########## PERFORMING PREDICTIONS FOR TRUE FLIGHTS #############

+-- prediction
|   +-- fpls
|   +-- errors_computation.py
|   +-- mc_prediction.py
|   +-- predict_descent_backwards.py
|   +-- predict_flight.py

'mc_prediction.py' performs predictions up to 20 minutes
look-ahead time for each (1) minute of a true flight progress time
and computes the ATE, CTE and AE for 5, 10, 15 and 20 minutes
look-ahead time and records the results into folders
like 'WTC/route/look-ahead time'.

################################################################


################ 5 #############################################
########## DATA ANALYSIS #######################################

+-- analysis
|   +-- lat_boxplot.py
|   +-- lat_distr.py
|   +-- lat_input_errors_corr_av_per_flight.py

'lat_boxplot.py' plots errors distribution as a boxplot for
look-ahead time of 5, 10, 15 and 20 minutes.

'lat_distr.py' plots errors distribution as a histogram for
look-ahead time of 5, 10, 15 and 20 minutes.

'lat_input_errors_corr_av_per_flight.py' - graphical or
tabular representation of the input-error correlation
coefficients.

################################################################