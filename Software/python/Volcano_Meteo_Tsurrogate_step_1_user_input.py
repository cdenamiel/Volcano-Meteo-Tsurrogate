#!/usr/bin/python3

import numpy as np
import scipy.io as sio
import os

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
print("-----------------------------")
print("STEP 1: user-defined criteria")
print("-----------------------------")
print(' ')
print(' ')
print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
print('Before performing this step:')
print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
print('Edit Volcano_Meteo_Tsurrogate_step_1_user_input to define:')
print('                                   - the stochastic vaiables: type, number and limits of the uniform distributions')
print('                                   - the volcanoes: number and coordinates')
print('                                   - the maximum total order')
print('                                   - the quadrature rule: Gauss-Patterson (GP) or Delayed Gauss-Patterson (DGP)')
print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
print('If the program is not edited: the user input below correspond to the 11 volcanoes described in the User Manual.')
print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')

def Volcano_Meteo_Tsurrogate_step_1_user_input():	
	
	"""
	USER SPECIFICATIONS
	"""
	#---------------------------------------------------------------------------------------------
	# In this section the user defines the number of stochastic variables and their distributions.
	# In version v1.0, all the distributions are assumed to be uniform. 
	#----------------------------------------------------------------------------------------------
		
	# Number of stochastic variables used to build the surrogate models
	nmodes = 4	

	# Limits [a,b] of the uniform distributions
	a=np.zeros(nmodes)
	b=np.zeros(nmodes)
	# uniforme distribution of the Lamb wave amplitude at the epicenter of the volcanic eruption [Pa] 
	a[0] = 1000.0
	b[0] = 20000.0
	# uniforme distribution of the Lamb wave wavelength [m]
	a[1] = 300000.0
	b[1] = 900000.0
	# uniforme distribution of the Lamb wave speed of propagation [m/s]
	a[2] = 280.0
	b[2] = 340.0
	# uniforme distribution of the Lamb wave time attenuation [s]
	a[3] = 3600.0
	b[3] = 43200.0
	
	# Number of volcanoes to simulate
	nvols = 11
	
	# Coordinates of the Volcano locations
	lon=np.zeros(nvols)
	lat=np.zeros(nvols)
	# Askja
	lon[0] = -16.7485
	lat[0] =  65.0111
	# Campi Flegrei
	lon[1] =  14.1390 
	lat[1] =  40.8270
	# Cotopaxi
	lon[2] = -78.4372
	lat[2] =  -0.6838 
	# Hunga Tonga-Hunga Ha'apai
	lon[3] =-175.3800
	lat[3] = -20.5700
	# Katla
	lon[4] = -19.1303
	lat[4] =  63.6467
	# Mt. Pinatubo
	lon[5] = 120.3496
	lat[5] =  15.1429
	# Popocatepetl
	lon[6] = -98.6279
	lat[6] =  19.0224 	
	# Mt. St. Helens
	lon[7] =-122.1956
	lat[7] =  46.1914
	# Mt. Sakurajima
	lon[8] = 130.6500
	lat[8] =  31.5833
	# Vesuvius
	lon[9] =  14.4289
	lat[9] =  40.8224 
	# Yellowstone
	lon[10] =-110.7232
	lat[10] =  44.4123	
	
	# Maximum total order of the Legendre polynomials that will be tested
	maxdeg = 6
	
	# type of nested grid option = 0 for Gauss-Patterson (GP) & option = 1 for Delayed Gauss-Patterson (DGP)
	option = 1
	
	"""
	END USER SPECIFICATIONS
	"""
	# save as mat file in order to be used later for the construction, convergence, evaluation, sensitivity of the surrogate models
	sio.savemat('../results/output_users.mat', {'nmodes': nmodes,'a': a,'b': b,'nvols': nvols,'lon': lon,'lat': lat,'maxdeg': maxdeg,'option': option})
	
	print(" ")
	print("The results have been saved in: ../results/output_users.mat")
	print(" ")

if __name__ == '__main__':
    Volcano_Meteo_Tsurrogate_step_1_user_input()
    
