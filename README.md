# Volcano-Meteo-Tsurrogate-V1.0

**Contact**: Clea Denamiel <br>
**email (institutional)**: cdenami@irb.hr <br>
**email (permanent)**: clea.denamiel@live.fr <br>

The project contains three main folders, **Software**, **User_Manual** and **GUI**, that describe, document and share the codes of the Volcano-Meteo-Tsurrogate v1.0 model.  

## Software

Contains two subfolders with the Landslide-Tsurrogate v1.0 programs and routines written in **Python**.

The python folder contains three subfolders: **code**, **data** and **results**. 

<img src="https://github.com/user-attachments/assets/7070e830-97b9-43db-81c6-57c52b944673" width="300" height="300"/><br>

The code folder contains the main Volcano_Meteo_Tsurrogate programs. The workflow pipeline of the code follows the steps below:

**STEP 1**: Edit and run Volcano_Meteo_Tsurrogate_step_1_user_input to generate the file: ../results/output_users.mat <br>
**STEP 2**: Run Volcano_Meteo_Tsurrogate_step_2_input_parameters to generate the file: ../results/output_param.mat <br>
**STEP 3**: <br>
* Prepare the locations where to create the surrogates and generate: ../data/surrogate_model_locations.mat <br>
  surrogate_model_locations.mat: contains loc_lat[nl] double (latitude), loc_lon[nl] double (longitude) and loc_name[nl] cell (associated name) with nl the number of surrogate models to build <br>
* Run the deterministic simulations outside Volcano-Meteo-Tsurrogate v1.0 and generate: ../data/input_simus.nc <br>
  **netcdf input_simus {**
  
  dimensions:

        volc = 11 ;	 % number of volcanoes
        city = 191 ;	% number of coastal locations
        simu = 641 ;	% number of simulations 

  variables:

        double zeta_max(volc, city, simu) ;
                zeta_max:units = "m" ;
                zeta_max:long_name = "meteotsunami surge" ;
                zeta_max:_FillValue = 9.96920996838687e+36 ;
        double time_zeta_max(volc, city, simu) ;
                time_zeta_max:units = "hours" ;
                time_zeta_max:long_name = "time of arrival of the meteotsunami surge" ;
                time_zeta_max:_FillValue = 9.96920996838687e+36 ;
        double vel_max(volc, city, simu) ;
                vel_max:units = "m/s" ;
                vel_max:long_name = "meteotsunami maximum speed" ;
                vel_max:_FillValue = 9.96920996838687e+36 ;
        double time_vel_max(volc, city, simu) ;
                time_vel_max:units = "hours" ;
                time_vel_max:long_name = "time of arrival of the meteotsunami maximum speed" ;
                time_vel_max:_FillValue = 9.96920996838687e+36 ;
  **}** <br>

**STEP 4**: Run Volcano_Meteo_Tsurrogate_step_4_format_input to generate the file: ../results/output_model.mat <br>
**STEP 5**: Run Volcano_Meteo_Tsurrogate_step_5_coefficients to generate the file: ../results/output_coeff.mat <br>
**STEP 6**: Run Volcano_Meteo_Tsurrogate_step_6_evaluation to generate the files: ../results/output_evals.mat and ../results/output_sensi.mat  <br>
**STEP 7**: Run Volcano_Meteo_Tsurrogate_step_7_PTHA to generate the file: ../results/output_PTHA.mat <br>

All these steps can be run with the provided files for the test cases: ../data/surrogate_model_locations.mat is already provided and input_simus.nc should be downloaded separately (see Code/matlab/data/README.md or Code/python/data/README.md). 

The results are created in the ../results folder.

## User Manual

The  User Manual is a Jupyter Notebook that describes, and illustrate for the test cases, the different steps and functions needed to build 
surrogate models based on gPCE for volcano generated meteotsunamis.

The User Manual consists in:

- the jupyter notebook (JN): User_Manual.ipynb
- 3 different sub-folders: <br>
    * **data**: all the data used to build the surrogate models for the Mayotte test case <br>
    * **figures**: all the figures used in the JN <br>
    * **results**: an empty folder where the results from the JN are copied <br>

## GUI

### Python  

Jupyter Widget Notebook under GUI/python with possibility to deploy it as a web app through Voila for users not familiar with python.


