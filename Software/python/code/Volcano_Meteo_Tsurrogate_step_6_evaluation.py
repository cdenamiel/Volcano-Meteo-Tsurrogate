#!/usr/bin/python3

import scipy.io as sio
import numpy as np
from pathlib import Path as pth
from math import comb, exp
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
print("----------------------------------------------------------------------------------")
print("STEP 6: evaluation of the surrogate model convergence, performance and sensitivity")
print("----------------------------------------------------------------------------------")
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
print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
print(' ')

def Volcano_Meteo_Tsurrogate_step_6_evaluation():

        #--------------
	# Load the data
	#--------------
	
	users = pth('../results/output_users.mat')
	if users.is_file():
	   # Load the user-defined input
	   maxdeg = sio.loadmat(users)['maxdeg'].flatten()
	   a = np.squeeze(np.array(sio.loadmat(users)['a']))
	   b = np.squeeze(np.array(sio.loadmat(users)['b']))
	   nvols  = sio.loadmat(users)['nvols'].flatten()
	else:
	   print("The file: ",users," does not exist.")
	   print("Edit and run Volcano_Meteo_Tsurrogate_step_1_user_input.")
	   return
	my_file = pth('../results/output_param.mat')
	if my_file.is_file():
	   param = sio.loadmat(my_file)['param']
	else:
	   print("The file: ",my_file," does not exist.")
	   print("Run Volcano_Meteo_Tsurrogate_step_2_input_parameters.")
	   return
	my_file = pth('../results/output_model.mat')
	if my_file.is_file():	 
	   model = sio.loadmat(my_file)['model']
	else:
	   print("The file: ",my_file," does not exist.")
	   print("Run Volcano_Meteo_Tsurrogate_step_4_format_input.")
	   return
	my_file = pth('../results/output_coeff.mat')
	if my_file.is_file():	 
	   coeff = sio.loadmat(my_file)['coeff']
	else:
	   print("The file: ",my_file," does not exist.")
	   print("Run Volcano_Meteo_Tsurrogate_step_5_coefficients.")
	   return
	   	   
        #---------------------------------------------
	# Generate data for evaluation and sensitivity
	#---------------------------------------------
	
	evals = np.empty((nvols[0], maxdeg[0]+2), dtype=object)
	for n in range(nvols[0]):
	    model_temp = model[n]
	    coeff_temp = coeff[n]
	    evals_temp = surrogate_model_gauss_patterson_evals(param,model_temp,coeff_temp,a,b,maxdeg[0])
	    evals[n] = np.array(evals_temp)
	sio.savemat('../results/output_evals.mat',{'evals': evals})
	
	print(" ")
	print("The results have been saved in: ../results/output_evals.mat")
	print(" ")
	
	sensi = np.empty((nvols[0]), dtype=object)
	for n in range(nvols[0]):
	    coeff_temp = coeff[n]
	    sensi_temp = surrogate_model_gauss_patterson_sensi(coeff_temp,maxdeg[0]-1)
	    sensi[n] = np.array(sensi_temp)
	sio.savemat('../results/output_sensi.mat',{'sensi': sensi})
	
	print(" ")
	print("The results have been saved in: ../results/output_sensi.mat")
	print(" ")

def surrogate_model_gauss_patterson_evals(param, model, coeff, a, b, maxdeg):
    """
    Evaluation of Pseudo Spectral Approximation with Gauss-Patterson Sparse grids.
    """
    
    nmodes = len(a)
    nx = len(model[0][0]['Zeta_Z'][0][0].flatten())
    
    # Build Legendre polynomials up to maxdeg (Le[0] = constant 1, Le[1] = x, etc.)
    Le = [None] * (maxdeg + 1)
    Le[0] = np.array([1.0])  # L_0 = 1
    Le[1] = np.array([1.0, 0.0])  # L_1 = x
    for n in range(2, maxdeg + 1):
        Le[n] = ((2 * n - 1) / n) * np.concatenate((Le[n - 1], [0])) - ((n - 1) / n) * np.concatenate(([0, 0], Le[n - 2]))
    # Prepare error structure as a list of dicts
    evals = []
    # Retrieve all simulations  
    Z = np.array(param[0][maxdeg]['Z'][0][0])
    index = param[0][maxdeg]['index'][0][0].flatten()
    ns = Z.shape[0]
    Z = Z.T
    # Preallocate true-model matrices
    Zeta_Z = np.empty((nx, ns))
    Velo_Z = np.empty((nx, ns))
    Time_Z = np.empty((nx, ns))
    for i in range(ns):
        entry = model[index[i]]
        Zeta_Z[:, i] = np.array(entry['Zeta_Z'][0][0])
        Velo_Z[:, i] = np.array(entry['Velo_Z'][0][0])
        Time_Z[:, i] = np.array(entry['Time_Zeta_Z'][0][0])
    # Store the true values as error[0]
    evals.append({
        'Zeta_Z': Zeta_Z.copy(),
        'Velo_Z': Velo_Z.copy(),
        'Time_Z': Time_Z.copy()
    })
    # Loop over increasing maximal gPCE degree
    for nmaxdeg in range(0, maxdeg+1):    
        # initialize PCE approximations
        zeta_PCE = np.zeros((nx, ns))
        velo_PCE = np.zeros((nx, ns))
        time_PCE = np.zeros((nx, ns))        
        # Smolyak summation over alpha norms
        for alpha_norm1 in range(max(0, nmaxdeg - nmodes + 1), nmaxdeg + 1):
            # Smolyak coefficient
            C_alpha = ((-1)**(nmaxdeg - alpha_norm1) 
                       * comb(nmodes-1, nmaxdeg - alpha_norm1))            
            # Retrieve multi-indices and PCE coefficients for this alpha_norm1
            this = coeff[alpha_norm1]
            alpha_mat = np.array(this['alpha'][0][0]) 
            zeta_hat = np.array(this['zeta'][0][0])   
            velo_hat = np.array(this['velo'][0][0])
            time_hat = np.array(this['time'][0][0])
            nw = alpha_mat.shape[0]            
            # For each multi-index row
            for l in range(nw):
                alpha_row = alpha_mat[l]           
                # Build the product of univariate polynomials
                # multH: shape (ns,)
                multH = np.ones(ns)
                for m in range(nmodes):
                    x_scaled = (2 * Z[m, :] - a[m] - b[m]) / (b[m] - a[m])
                    P = np.polyval(Le[alpha_row[m]], x_scaled)
                    multH *= P                
                # Accumulate PCE contributions
                # Note: zeta_hat[:, l] is shape (nx,), so we outer with multH
                zeta_PCE += C_alpha * np.outer(zeta_hat[:, l], multH)
                velo_PCE += C_alpha * np.outer(velo_hat[:, l], multH)
                time_PCE += C_alpha * np.outer(time_hat[:, l], multH) 
        # Make sure time is never negative
        time_PCE[time_PCE < 0.0] = np.nan
        print((np.isnan(time_PCE).sum()*100)/time_PCE.size,"% of negative time of arrival values were generated at total order p = ", alpha_norm1)	   
        # Exponentiate to undo the log transform and get the real values of the original quantities
        evals.append({
            'Zeta_PCE': np.exp(zeta_PCE),
            'Velo_PCE': np.exp(velo_PCE),
	    # no exponential and convert time of arrival in hours
            'Time_PCE': time_PCE/3600.0
        })
    
    return evals
    
def surrogate_model_gauss_patterson_sensi(coeff, maxdeg):
    import numpy as np
    from scipy.special import comb

    nmaxdeg = maxdeg - 1
    alpha0 = np.array(coeff[0]['alpha'][0][0])
    nx = np.array(coeff[0]['zeta'][0][0]).shape[0]
    nmodes = alpha0.shape[1]
    
    # Preallocate total sensitivity arrays
    sensi = []   
    ST_Zeta = np.full((nx, nmodes), np.nan)
    ST_Velo = np.full((nx, nmodes), np.nan)
    ST_Time = np.full((nx, nmodes), np.nan)

    # Gather all alpha and weighted coefficients
    alpha_list = []
    coeff_zeta_list = []
    coeff_velo_list = []
    coeff_time_list = []

    for alpha_norm1 in range(max(0, nmaxdeg - nmodes + 1), nmaxdeg + 1):
        C_alpha = (-1)**(nmaxdeg - alpha_norm1) * comb(nmodes - 1, nmaxdeg - alpha_norm1)
        entry = coeff[alpha_norm1]
        alpha_list.append(entry['alpha'][0][0])
        coeff_zeta_list.append(C_alpha * entry['zeta'][0][0])
        coeff_velo_list.append(C_alpha * entry['velo'][0][0])
        coeff_time_list.append(C_alpha * entry['time'][0][0])

    alpha = np.vstack(alpha_list)
    coeff_zeta = np.hstack(coeff_zeta_list)
    coeff_velo = np.hstack(coeff_velo_list)
    coeff_time = np.hstack(coeff_time_list)

    # Remove purely constant terms
    non_constant = ~(np.sum(alpha == 0, axis=1) == nmodes)
    alpha_nc = alpha[non_constant, :]
    coeff_zeta_nc = coeff_zeta[:, non_constant]
    coeff_velo_nc = coeff_velo[:, non_constant]
    coeff_time_nc = coeff_time[:, non_constant]

    # Total variances
    D_zeta = np.sum(coeff_zeta_nc**2, axis=1)
    D_velo = np.sum(coeff_velo_nc**2, axis=1)
    D_time = np.sum(coeff_time_nc**2, axis=1)

    # Compute total Sobol indices
    for i in range(nmodes):
        indT = np.where(alpha_nc[:, i] > 0)[0]
        ST_Zeta[:, i] = np.sum(coeff_zeta_nc[:, indT]**2, axis=1) / D_zeta
        ST_Velo[:, i] = np.sum(coeff_velo_nc[:, indT]**2, axis=1) / D_velo
        ST_Time[:, i] = np.sum(coeff_time_nc[:, indT]**2, axis=1) / D_time
    
    # Store final results
    sensi.append({'ST_Zeta': ST_Zeta.copy(),
                        'ST_Velo': ST_Velo.copy(),
                        'ST_Time': ST_Time.copy()
                       })
    return sensi

if __name__ == '__main__':
    Volcano_Meteo_Tsurrogate_step_6_evaluation()
    
