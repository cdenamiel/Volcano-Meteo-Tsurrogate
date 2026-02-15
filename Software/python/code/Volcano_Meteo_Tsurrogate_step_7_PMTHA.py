#!/usr/bin/python3

import scipy.io as sio
import numpy as np
from math import comb, exp
from pathlib import Path as pth
import os
import time

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
print("--------------------------------------------------------------------")
print("******************************  PMTHA  *****************************")
print("--------------------------------------------------------------------")
print(' ')
print('***************  EDIT THIS SCRIPT BEFORE RUNNING IT  ***************')
print(' ')
print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
print('Before performing the PMTHA:')
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
print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
print(' ')

def Volcano_Meteo_Tsurrogate_step_7_PMTHA():

        #------------
	# User Inputs
	#------------

	#--------------------------------------------------------------------------------------------
	# In this section the user defines the limits of the distributions to generate the PMTHA 
	# The default values and setup are to reproduce the 2022 Hunga Tonga-Hunga Ha'apai eruption
	#--------------------------------------------------------------------------------------------
	        
	# chosen volcano 
	# for the test cases: 0=Askja, 1=Campi Flegrei, 2=Cotopaxi, 3=Hunga Tonga, 4=Katla, 5=Mt.Pinatubo, 
	#                     6=Popocatepetl, 7=Mt.St.Helens, 8=Sakurajima, 9=Vesuvius, 10=Yellowstone
	n = 3
	
	# Parameters of the Lamb wave
	a0=np.zeros(4)
	b0=np.zeros(4)
	# Amplitude of the pressure disturbance at the epicenter [Pa]
	a0[0] = 1600.0
	b0[0] = 1800.0
	# Wavelength of the pressure disturbance [m]
	a0[1] = 550000.0
	b0[1] = 650000.0
	# Speed of the pressure disturbance [m/s]
	a0[2] = 310.0
	b0[2] = 320.0
	# Attenuation in time [s]
	a0[3] = 4.0*3600.0
	b0[3] = 6.0*3600.0
	
	# chosen total order maxdeg0 such as 1 < maxdeg0 < 6
	maxdeg0 = 5
	
	# chosen number of samples for the PTHA
	nw0 = 1000

        #--------------
	# Load the data
	#--------------
	
	users = pth('../results/output_users.mat')
	if users.is_file():
	   # Load the user-defined input	   
	   nmodes = np.squeeze(np.array(sio.loadmat(users)['nmodes']))	 
	   if (len(a0) != nmodes) or (len(b0) != nmodes):
	      print("The defined input for PMTHA")
	      print("should be the same size than")
	      print("the number of user-defined")
	      print("stochastic variables.")
	      return
	   a = np.squeeze(np.array(sio.loadmat(users)['a']))
	   b = np.squeeze(np.array(sio.loadmat(users)['b']))
	else:
	   print("The file: ",users," does not exist.")
	   print("Edit and run Landslide_Tsurrogate_user_input.")
	   return
	my_file = pth('../results/output_coeff.mat')
	if my_file.is_file():
	   coeff = sio.loadmat(my_file)['coeff']
	else:
	   print("The file: ",my_file," does not exist.")
	   print("Generate the file as instructed above.")
	   return
	   
    #--------------------------
	# Generate the PMTHA
	#--------------------------
	
	start = time.time()	
	PMTHA_zeta, PMTHA_velo, PMTHA_time = surrogate_model_gauss_patterson_PMTHA(maxdeg0, a, b, coeff[n], nw0, a0, b0)
	end = time.time()
	length = end - start
	
	print(" ")
	print(" ---------------------------------------------------------------")
	print("PMTHA based on ",nw0," members produced in ",length, "seconds!")
	print(" ---------------------------------------------------------------")
	
	sio.savemat('../results/output_PMTHA.mat',{'PMTHA_zeta': PMTHA_zeta,'PMTHA_velo': PMTHA_velo,'PMTHA_time': PMTHA_time})
	
	print(" ")
	print("The results have been saved in: ../results/output_PMTHA.mat")
	print(" ")

def surrogate_model_gauss_patterson_PMTHA(maxdeg, a, b, coeff, nw, ar, br):
        
    nmodes = len(ar)
    nx = len(coeff[0]['zeta'][0][0].flatten())
    
    # Build Legendre polynomials up to maxdeg (Le[0] = constant 1, Le[1] = x, etc.)
    Le = [None] * (maxdeg + 1)
    Le[0] = np.array([1.0])  # L_0 = 1
    Le[1] = np.array([1.0, 0.0])  # L_1 = x
    for n in range(2, maxdeg + 1):
        Le[n] = ((2 * n - 1) / n) * np.concatenate((Le[n - 1], [0])) - ((n - 1) / n) * np.concatenate(([0, 0], Le[n - 2]))
        
    # generate random input  
    np.random.seed(0)
    seed_01 = np.random.rand(nmodes, nw)    
    # rescaling to fit the appropriate intervals
    Zw = np.tile(ar, (nw,1)).T + np.tile((br - ar), (nw,1)).T * seed_01
               
    # Preallocate true-model matrices
    zeta_temp = np.zeros((nx, nw))
    velo_temp = np.zeros((nx, nw))
    time_temp = np.zeros((nx, nw))
    
    for alpha_norm1 in range(max(0, maxdeg - nmodes + 1), maxdeg + 1):
        # Smolyak coefficient
        C_alpha = ((-1)**(maxdeg - alpha_norm1) * comb(nmodes - 1, maxdeg - alpha_norm1))            
        # Retrieve multi-indices and PCE coefficients for this alpha_norm1
        this = coeff[alpha_norm1]
        alpha    = np.array(this['alpha'][0][0]) 
        zeta_hat = np.array(this['zeta'][0][0])   
        velo_hat = np.array(this['velo'][0][0])
        time_hat = np.array(this['time'][0][0])
        nww = alpha.shape[0]            
        # For each multi-index row
        for l in range(nww):
            multH = np.ones(nw)
            for n in range(nmodes):
                x_scaled = (2 * Zw[n, :] - a[n] - b[n]) / (b[n] - a[n])
                multH *= np.polyval(Le[alpha[l,n]], x_scaled)
            for i in range(nx):
                zeta_temp[i, :] += C_alpha * zeta_hat[i, l] * multH
                velo_temp[i, :] += C_alpha * velo_hat[i, l] * multH
                time_temp[i, :] += C_alpha * time_hat[i, l] * multH

    # Exponentiate to undo the log transform and get the real values of the original quantities
    PMTHA_zeta = np.exp(zeta_temp)
    PMTHA_velo = np.exp(velo_temp)
    time_temp[time_temp < 0.0] =np.nan
    PMTHA_time = time_temp/3600.0
    
    return PMTHA_zeta, PMTHA_velo, PMTHA_time

if __name__ == '__main__':
    Volcano_Meteo_Tsurrogate_step_7_PMTHA()
