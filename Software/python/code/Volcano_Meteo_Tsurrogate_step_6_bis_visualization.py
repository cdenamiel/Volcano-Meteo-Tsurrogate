#!/usr/bin/python3
import scipy.io as sio
import numpy as np
from pathlib import Path as pth
import matplotlib.pyplot as plt

print(' ')
print("-----------------------------")
print("Volcano-Meteo-Tsurrogate v1.0")
print("-----------------------------")
print("")
print("contact: Clea Denamiel")
print("email: clea.denamiel@live.fr")
print("")
print("Look at the User Manual and the Article Draft for more information.")
print("")
print("-------------------------------------------")
print("*************  VISUALIZATION  *************")
print("-------------------------------------------")
print(' ')
print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
print('Before performing this step:')
print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
print('STEP 1: Edit and run Volcano_Meteo_Tsurrogate_step_1_user_input to generate the file: ../results/output_users.mat')
print('STEP 2: Run Volcano_Meteo_Tsurrogate_step_2_input_parameters to generate the file: ../results/output_param.mat')
print('STEP 3: Run the deterministic simulations outside Volcano-Meteo-Tsurrogate v1.0 and generate: ../data/input_simus.nc')
print('        input_simus.nc: contains the maximum elevation, the maximum speed and the time of arrival of the meteotsunamis.')
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
print('        surrogate_model_locations.mat: contains loc_lon[nl x 1 double], loc_lat[nl x 1 double] and loc_name{nl x 1 cell}, ')
print('        ------------------------------  longitude, latitude and name of the locations and nl the number of locations')
print('STEP 4: Run Volcano_Meteo_Tsurrogate_step_4_format_input to generate the file: ../results/output_model.mat')
print('STEP 5: Run Volcano_Meteo_Tsurrogate_step_5_coefficients to generate the file: ../results/output_coeff.mat')
print('STEP 6: Run Volcano_Meteo_Tsurrogate_step_6_evaluation to generate the files: ../results/output_evals.mat')
print('                                                                   and        ../results/output_sensi.mat')
print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
print(' ')

def Volcano_Meteo_Tsurrogate_step_6_bis_visualization():

        #--------------
	# Load the data
	#--------------
	
	users = pth('../results/output_users.mat')
	if users.is_file():
	   # Load the user-defined input
	   maxdeg = sio.loadmat(users)['maxdeg'].flatten()
	   nvols = sio.loadmat(users)['nvols'].flatten()
	else:
	   print("The file: ",users," does not exist.")
	   print("Edit and run Volcano_Meteo_Tsurrogate_user_input.")
	   return
	
	my_file = pth('../results/output_param.mat')
	if my_file.is_file():
	   param = sio.loadmat(my_file)['param']
	else:
	   print("The file: ",my_file," does not exist.")
	   print("Generate the file as instructed above.")
	   return

	my_file = pth('../results/output_evals.mat')
	if my_file.is_file():	 
	   evals = sio.loadmat(my_file)['evals']
	else:
	   print("The file: ",my_file," does not exist.")
	   print("Generate the file as instructed above.")
	   return
	
	my_file = pth('../results/output_sensi.mat')
	if my_file.is_file():	 
	   sensi = sio.loadmat(my_file)['sensi']
	else:
	   print("The file: ",my_file," does not exist.")
	   print("Generate the file as instructed above.")
	   return	

	user_input = input(f"Please enter a volcano number between 0 and {nvols[0]-1}: ")
	value = int(user_input)   
	
	print(" ")
	print("---------------------------------------------------------------------")
	print(f"Evaluating the gPCE surrogate models for the volcano number {value}.")
	print("---------------------------------------------------------------------")
	print(" ")

	#------------------------------------
	# Convergence of the surrogate models
	#------------------------------------
	
	nx = np.array(evals[value][0]['Zeta_Z'][0][0]).shape[0]
	maxdeg = maxdeg[0]
	norm_err_zeta = np.full((maxdeg+1, nx), np.nan)
	norm_err_velo = np.full((maxdeg+1, nx), np.nan)
	norm_err_time = np.full((maxdeg+1, nx), np.nan)
	den_time = np.full((nx), np.nan)
	
	den_zeta = np.sum((np.array(evals[value][0]['Zeta_Z'][0][0]) - np.array(evals[value][1]['Zeta_PCE'][0][0]))**2, axis=1)
	den_velo = np.sum((np.array(evals[value][0]['Velo_Z'][0][0]) - np.array(evals[value][1]['Velo_PCE'][0][0]))**2, axis=1)
	# time is treated differently as log transform was not applied and negative values must be eliminated
	for i in range(nx):
            temp_z   = np.array(evals[value][0]['Time_Z'][0][0])[i,:]
            temp_pce = np.array(evals[value][1]['Time_PCE'][0][0])[i,:]
            temp_z   = temp_z[np.isfinite(temp_pce)]
            temp_pce = temp_pce[np.isfinite(temp_pce)]
            den_time[i] = np.sum((temp_z-temp_pce)**2, axis=0)
	for d in range(1,maxdeg+2):
            num_zeta = np.sum((np.array(evals[value][0]['Zeta_Z'][0][0]) - np.array(evals[value][d]['Zeta_PCE'][0][0]))**2, axis=1)
            num_velo = np.sum((np.array(evals[value][0]['Velo_Z'][0][0]) - np.array(evals[value][d]['Velo_PCE'][0][0]))**2, axis=1)
            norm_err_zeta[d-1, :] = num_zeta / den_zeta
            norm_err_velo[d-1, :] = num_velo / den_velo
	    # time is treated differently as log transform was not applied and negative values must be eliminated
            num_time = np.full((nx), np.nan)
            for i in range(nx):
                temp_z   = np.array(evals[value][0]['Time_Z'][0][0])[i,:]
                temp_pce = np.array(evals[value][d]['Time_PCE'][0][0])[i,:]
                temp_z   = temp_z[np.isfinite(temp_pce)]
                temp_pce = temp_pce[np.isfinite(temp_pce)]
                num_time[i] = np.sum((temp_z-temp_pce)**2, axis=0)
            norm_err_time[d-1, :] = num_time / den_time
        
	plot_error(norm_err_zeta.T, 'Convergence - Maximum Meteo-Tsunami Elevation')
	plot_error(norm_err_velo.T, 'Convergence - Maximum Meteo-Tsunami Speed')
	plot_error(norm_err_time.T, 'Convergence - Meteo-Tsunami Time of Arrival')
	
	#------------------------------------
	# Accuracy of the surrogate models
	#------------------------------------	
	
	# Extract the two index arrays
	pm1=maxdeg-1
	index_pm1 = param[0][maxdeg-1]['index'][0][0].flatten()
	index_p = param[0][maxdeg]['index'][0][0].flatten()
	# Find those in index_p that are not in index_pm1
	ind_independent = ~np.isin(index_p, index_pm1)
	Z_true = np.array(evals[value][0]['Zeta_Z'][0][0])
	Z_pce  = np.array(evals[value][maxdeg-1]['Zeta_PCE'][0][0])
	V_true = np.array(evals[value][0]['Velo_Z'][0][0])
	V_pce  = np.array(evals[value][maxdeg-1]['Velo_PCE'][0][0])
	T_true  = np.array(evals[value][0]['Time_Z'][0][0])
	T_pce  = np.array(evals[value][maxdeg-1]['Time_PCE'][0][0])
	simu_test_z = Z_true[:, ind_independent].ravel()
	pce_test_z  = Z_pce[:, ind_independent].ravel()
	simu_test_v = V_true[:, ind_independent].ravel()
	pce_test_v  = V_pce[:, ind_independent].ravel()
	simu_test_t = T_true[:, ind_independent].ravel()
	pce_test_t  = T_pce[:, ind_independent].ravel()		
	
	plot_hex(simu_test_z, pce_test_z,lim=5,vmax=50,xlabel='Determinisric simulations',ylabel='Total Order '+str(pm1),title='Accuracy - Maximum Meteo-Tsunami Elevation')
	plot_hex(simu_test_v, pce_test_v,lim=3,vmax=50,xlabel='Deterministic simulations',ylabel='Total Order'+str(pm1),title='Accuracy - Maximum Meteo-Tsunami Speed')
	plot_hex(simu_test_t[np.isfinite(pce_test_t)], pce_test_t[np.isfinite(pce_test_t)],lim=40,vmax=5,xlabel='Deterministic simulations',ylabel='Total Order level '+str(pm1),title='Accuracy - Meteo-Tsunami Time of Arrival')
        
	#------------------------------------
	# Sensitivity of the surrogate models
	#------------------------------------

	sensi_dict = matlab_struct_to_dict(sensi[0, value])
	ST_Zeta = np.array(sensi_dict['ST_Zeta'])
	ST_Velo = np.array(sensi_dict['ST_Velo'])
	ST_Time = np.array(sensi_dict['ST_Time'])
	
	locations = np.arange(ST_Zeta.shape[0])
	labels = ['Pressure Amplitude', 'Wavelength', 'Speed', 'Time attenuation']
	
	plot_stacked(ST_Zeta, locations, labels, 'Sensitivity - Maximum Meteo-Tsunami Elevation')
	plot_stacked(ST_Velo, locations, labels, 'Sensitivity - Maximum Meteo-Tsunami Speed')
	plot_stacked(ST_Time, locations, labels, 'Sensitivity - Meteo-Tsunami Time of Arrival')


def plot_error(matrix, title):
    plt.figure()
    plt.pcolormesh(matrix, shading='gouraud', cmap='Reds', vmin=0, vmax=1)
    plt.xlabel('Total Order', fontsize=14)
    plt.ylabel('Locations', fontsize=14)
    plt.title(title, fontsize=16)
    plt.colorbar(label='Normalized Squared Error')
    plt.gca().tick_params(labelsize=14)
    plt.gcf().set_facecolor('w')
    plt.show()
    
def plot_hex(simu, pce, lim, vmax, cmap='Reds', xlabel='', ylabel='', title=''):
    plt.figure()
    plt.plot([0, lim], [0, lim], 'k--', linewidth=1)
    hb = plt.hexbin(simu, pce, gridsize=1000, xscale=[0, lim], yscale=[0, lim], cmap=cmap, vmin=0, vmax=vmax)
    plt.colorbar(hb, label='count')
    plt.axis('square')
    plt.xlim(0, lim)
    plt.ylim(0, lim)
    plt.xlabel(xlabel, fontsize=14)
    plt.ylabel(ylabel, fontsize=14)
    plt.title(title, fontsize=16)
    plt.gca().tick_params(labelsize=14)
    plt.tight_layout()
    plt.show()
    
def plot_stacked(ST, locations, labels, title):
    fig, ax = plt.subplots(figsize=(10, 6))
    bottom = np.zeros(ST.shape[0])    
    # For each variable, stack its bar segment
    for i in range(ST.shape[1]):
        ax.bar(locations, ST[:, i], bottom=bottom, label=labels[i])
        bottom += ST[:, i]    
    ax.set_facecolor('white')
    fig.patch.set_facecolor('white')
    ax.set_title(title, fontsize=20)
    ax.set_xlabel('Coastal Locations', fontsize=16)
    ax.set_ylabel('Total sensitivity', fontsize=16)
    ax.legend(fontsize=12)
    ax.tick_params(axis='both', which='major', labelsize=14)
    plt.tight_layout()
    plt.show()   

def matlab_struct_to_dict(mat_struct):
    """
    Recursively convert MATLAB structs loaded by scipy.io.loadmat to Python dicts.
    """
    if isinstance(mat_struct, np.ndarray):
        # Unwrap arrays of size 1
        if mat_struct.size == 1:
            return matlab_struct_to_dict(mat_struct[0])
        else:
            return [matlab_struct_to_dict(x) for x in mat_struct]
    elif isinstance(mat_struct, np.void):  # MATLAB struct
        return {name: matlab_struct_to_dict(mat_struct[name]) for name in mat_struct.dtype.names}
    else:
        return mat_struct     

if __name__ == '__main__':
    Volcano_Meteo_Tsurrogate_step_6_bis_visualization()
