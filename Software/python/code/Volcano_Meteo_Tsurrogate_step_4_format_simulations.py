#!/usr/bin/python3

import numpy as np
import scipy.io as sio
import xarray as xr
from pathlib import Path as pth

print(' ')
print("-----------------------------")
print("Volcano-Meteo-Tsurrogate v1.0")
print("-----------------------------")
print("")
print("contact: Clea Denamiel")
print("email: clea.denamiel@live.fr")
print("")
print("Look at the User Manual for more information.")
print("")
print("----------------------------------------------------------")
print("STEP 4: format the user-provided deterministic simulations")
print("----------------------------------------------------------")	
print(' ')
print(' ')
print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
print('Before performing this step:')
print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
print('STEP 1: Edit and run Volcano_Meteo_Tsurrogate_step_1_user_input to generate the file: ../results/output_users.mat')
print('STEP 2: Run Volcano_Meteo_Tsurrogate_step_2_input_parameters to generate the file: ../results/output_param.mat')
print('STEP 3: Run the deterministic simulations outside Volcano-Meteo-Tsurrogate v1.0 and generate: ../data/input_simus.nc')
print('        input_simus.nc:  contains the maximum elevation, the maximum speed and the time of arrival of the meteotsunamis.')
print('        ---------------  ')
print('        netcdf input_simus {')
print('        dimensions:')
print('            volc = 11 ;  % number of volcanoes')
print('            city = 191 ;	% number of coastal locations')
print('            simu = 641 ;	% number of simulations') 
print('        variables:')
print('            double zeta_max(volc, city, simu) ;')
print('                zeta_max:units = "m" ;')
print('                zeta_max:long_name = "meteotsunami surge" ;')
print('                zeta_max:_FillValue = 9.96920996838687e+36 ;')
print('            double time_zeta_max(volc, city, simu) ;')
print('                time_zeta_max:units = "hours" ;')
print('                time_zeta_max:long_name = "time of arrival of the meteotsunami surge" ;')
print('                time_zeta_max:_FillValue = 9.96920996838687e+36 ;')
print('            double vel_max(volc, city, simu) ;')
print('                vel_max:units = "m/s" ;')
print('                vel_max:long_name = "meteotsunami maximum speed" ;')
print('                vel_max:_FillValue = 9.96920996838687e+36 ;')
print('        }')
print(' ')
print('        Provide the locations where to the surrogates are created and generate: ../data/surrogate_model_locations.mat')
print('        surrogate_model_locations.mat: contains loc_lon[nl x 1 double], loc_lat[nl x 1 double] and loc_name{nl x 1 cell},')
print('        ------------------------------  longitude, latitude and name of the locations and nl the number of locations')
print(' ')
print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
print('If the program and data are not edited this program formats the results from the 11 volcanoes of the test cases as') 
print('described in the User Manual.')
print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')

def Volcano_Meteo_Tsurrogate_step_4_format_simulations():	
				
	#-------------------------
	# Load general information
	#-------------------------
	# Read the user-defined input in output_users.mat
	users = pth('../results/output_users.mat')
	if users.is_file():
	   nvols = sio.loadmat(users)['nvols'].flatten()
	else:
	   print("The file: ",users," does not exist.")
	   print("Edit and run Volcano_Meteo_Tsurrogate_step_1_user_input.")
	   return
	# Read output_param.mat created by Volcano_Meteo_Tsurrogate_step_2_input_parameters
	my_file = pth('../results/output_param.mat')
	if my_file.is_file():
	   param_data = sio.loadmat(my_file)
	   param = np.squeeze(np.array(param_data['param']))
	   maxdeg=len(param.flatten())-1
	   nsim=len(param[maxdeg]['index'][0][0].flatten())
	else:
	   print("The file: ",my_file," does not exist.")
	   print("Run Volcano_Meteo_Tsurrogate_step_2_input_param.")
	   return	
	# Read the user provided file: surrogate_model_locations.mat, containing loc_lon, loc_lat and loc_name
	my_file = pth('../data/surrogate_model_locations.mat')
	if my_file.is_file():
           data = sio.loadmat(my_file)
           loc_lon = data['loc_lon'].flatten()
           loc_lat = data['loc_lat'].flatten()
           loc_name = data['loc_name'].flatten()
	else:
	   print("The file: ",my_file," does not exist.")
	   print("Produce this file as instructed above.")
	   return
	# Read the user provided file: input_simus.nc, containing the results of the deterministic simulations
	my_file = pth('../data/input_simus.nc')
	if my_file.is_file():
           ds = xr.open_dataset(my_file, decode_timedelta=False)
           zeta_max = ds['zeta_max'].values
           time_zeta_max = ds['time_zeta_max'].values
           velo_max = ds['vel_max'].values
	else:
	   print("The file: ",my_file," does not exist.")
	   print("Produce this file as instructed above.")
	   return	   

	#-------------------------
	# Reformat the simulations
	#-------------------------

	# Initialize the model list
	model = np.empty((nvols[0], nsim), dtype=object)
	# Loop through each volcano and simulation and extract data
	for j in range(nvols[0]):
	    for i in range(nsim):
	        #print(time_zeta_max[j, :, i])
	        model[j, i] = {
		            'Zeta_Z': zeta_max[j, :, i], 
			    'Velo_Z': velo_max[j, :, i],
			    'Time_Zeta_Z': time_zeta_max[j, :, i]
			    }
	# Save the model to a .mat file
	sio.savemat('../results/output_model.mat', {'model': model})
	
	print(" ")
	print("The results have been saved in: ../results/output_model.mat")
	print(" ")
	
if __name__ == '__main__':
    Volcano_Meteo_Tsurrogate_step_4_format_simulations()
