# -*- coding: utf-8 -*-
"""
Created on Tue Jun  9 17:23:33 2020

Script uses the fitting function in iFit to fit a single clear spectrum across
a large wavelength range, then stores the fitting parameters. These fitting 
parameters can then be altered to add light dilution and varying SO2 
quantities, creating a series of synthetic spectra which should replicate the 
real response spectra to light dilution. These spectra can then be reanalysed 
to see how the fitting function responds to the changes.

@author: Matthew Varnam - The University of Manchester
@email : matthew.varnam(-at-)manchester.ac.uk
@version: 2.0
"""

# =============================================================================
# Libraries
# =============================================================================

# Import numpy for numerical calculations
import numpy as np

# iFit reads and analyses UV spectra of SO2
from ifit.load_spectra import read_spectrum, average_spectra
from ifit.parameters import Parameters
from ifit.spectral_analysis import Analyser

# iFit_mod is alterations to iFit to allow analysis of light diluted spectra
from ifit_mod.synthetic_suite import Analyser_ld

import matplotlib.pyplot as plt

# =============================================================================
# User settings
# =============================================================================

mode = 0

# Define the directory containing all work files
prm_dir = ('F:/Black_Hole/Data/201801_Masaya/Ryu_spectra/'+
          '20180115_DOAS_for_camera_all/Ben_format_shifted/')

#prm_dir = ('F:/Black_Hole/Data/201801_Masaya/Ryu_spectra/'+
#           '20180114/Ben_format_shifted/')

#prm_dir = ('F:/Black_Hole/Data/2012_0318_EtnaLavaFountain/Model_Input/'+
#           'EMIL/Zenith_7.45/')

# Set spectrometer name and directory
spec_name = 'FLMS02101_2'
spec_dir  = prm_dir + ('spectra/spectrum_')
drk_dir   = prm_dir + ('dark/spectrum_')

# Create list of spectra numbers
spec_num = [511,521]
#spec_num = [210,220] 
drk_num  = [  0,  0]

# Set model parameters
#so2_grid_ppmm = np.arange(0,5010,20)
#ldf_grid = np.arange(0,1.0,0.002)

so2_grid_ppmm = np.arange(0,3020,100)
ldf_grid = np.arange(0,1.0,0.1)

# Set wavebands to view
pad  = 1
wav0 = [306,316]
wav1 = [312,322]

# =============================================================================
# Model Setup
# =============================================================================

# Create parameter dictionary
params = Parameters()

# Add the gases
params.add('SO2',  value=1.0e16, vary=True, xpath='Ref/SO2_295K.txt')
params.add('O3',   value=1.0e19, vary=True, xpath='Ref/O3_233K.txt')
params.add('Ring', value=0.1,    vary=True, xpath='Ref/Ring.txt')
    
# Add background polynomial parameters
params.add('bg_poly0', value=0.0, vary=True)
params.add('bg_poly1', value=0.0, vary=True)
params.add('bg_poly2', value=0.0, vary=True)
params.add('bg_poly3', value=1.0, vary=True)

# Add intensity offset parameters
#params.add('offset0', value=0.0, vary=True)

# Add wavelength shift parameters
params.add('shift0', value=0.002, vary=True)
params.add('shift1', value=-1.1, vary=True)

# Add ILS parameters
params.add('fwem', value = 0.68 , vary = False)
params.add('k',    value = 2.58 , vary = False)
params.add('a_w',  value = 0.064, vary = False)
params.add('a_k',  value = 0.71 , vary = False)

# Generate the first analyser for setup
analyser_prm = Analyser(params,
                        fit_window = [wav0[0]-pad, wav1[1]+pad],
                        frs_path   = 'Ref/sao2010.txt',
                        flat_flag  = True,
                        flat_path  = 'Spectrometer/flat_FLMS02101_2.txt',
                        stray_flag = True,
                        dark_flag  = True,
                        ils_type   = 'Manual',
                        ils_path   = 'Spectrometer/ils_params_FLMS02101.txt')

# =============================================================================
# Pre-processing
# =============================================================================

# List of names of all spectra files
spec_range = np.arange(spec_num[0],spec_num[1]+1)
spec_fnames = [spec_dir + "{:05d}".format(num) + '.txt'for num in spec_range]

# List of names of all dark files
drk_range  = np.arange(drk_num[0],drk_num[1]+1)
drk_fnames = [drk_dir + "{:05d}".format(num) + '.txt'for num in drk_range]

if len(spec_fnames) == 1:
    # Load single specrum
    x, y, spec_info, read_err = read_spectrum(spec_fnames[0])

else:
    # Load spectra and average them to produce a single spectrum
    x, y = average_spectra(spec_fnames)

if len(drk_fnames) == 1:
    # Load dark spectrum from file
    x, dark, spec_info, read_err = read_spectrum(drk_fnames[0])

else:
    # Load spectra and average them to produce a single spectrum
    x, dark = average_spectra(drk_fnames)

analyser_prm.dark_spec = dark

# =============================================================================
# Run broad wavelength analysis
# =============================================================================

fit_prm = analyser_prm.fit_spectrum([x,y], calc_od=['SO2', 'Ring', 'O3'])
print(fit_prm.params.pretty_print())

# =============================================================================
# Create suite of synthetic spectra
# =============================================================================

# Create second diluted ifit analyser for setup
analyser_ld = Analyser_ld(params,
                          fit_window = [wav0[0]-pad, wav1[1]+pad],
                          frs_path   = 'Ref/sao2010.txt',
                          flat_flag  = True,
                          flat_path  = 'Spectrometer/flat_FLMS02101_2.txt',
                          stray_flag = True,
                          dark_flag  = True,
                          ils_type   = 'Manual')

analyser_ld.params.update_values(fit_prm.params.popt_list())
analyser_ld.params.add('LDF', value = 0.0, vary = True)
analyser_ld.interp_method = 'cubic'

# Convert SO2 in ppmm to molecules/cm2
so2_grid = np.multiply(so2_grid_ppmm , 2.652e+15)

# Use shape of grid and so2 value to produce a single array
shape = (len(so2_grid),len(fit_prm.spec),len(ldf_grid))
spectra_suite = np.zeros(shape)

# Create synthetic spectra by updating parameters
for i, so2 in enumerate(so2_grid):
    for j, ldf in enumerate(ldf_grid):
        
        # Update parameters of synthetic spectra to generate        
        analyser_ld.params['SO2'].set(value = so2)
        analyser_ld.params['LDF'].set(value = ldf) 
        
        # Extract parameter list
        fit_params = analyser_ld.params.fittedvalueslist()
        
        # Create synthetic spectrum
        spectra_suite[i,...,j] = analyser_ld.fwd_model(fit_prm.grid,*fit_params)

# =============================================================================
# Analyse spectra in first waveband
# =============================================================================
        
print('Analyse synthetic spectra in waveband 0')

# Create first ifit analyser for results 
analyser0 = Analyser(params,
                     fit_window = wav0,
                     frs_path   = 'Ref/sao2010.txt',
                     flat_flag  = False,
                     stray_flag = False,
                     dark_flag  = False,
                     ils_type   = 'Manual')

# Create arrays to store answers
ifit_so2_0 = np.zeros((shape[0],shape[2]))
ifit_err_0 = np.zeros((shape[0],shape[2]))

# Loop through each synthetic spectrum
for i, so2 in enumerate(so2_grid):
    for j, ldf in enumerate(ldf_grid):

        #Extract syntheteic spectrum for suite of spectra
        spectrum = [fit_prm.grid,spectra_suite[i,...,j]]
        
        # Analyse spectrum
        fit = analyser0.fit_spectrum(spectrum, calc_od=['SO2', 'Ring', 'O3'])
        
        # Store SO2 fit parameters in array
        ifit_so2_0[i,j] = fit.params['SO2'].fit_val
        ifit_err_0[i,j] = fit.params['SO2'].fit_err
        
    print(i)

#Create new ifit_so2 with units in ppm.m
ifit_so2_ppmm0 = np.divide(ifit_so2_0,2.652e15)
ifit_err_ppmm0 = np.divide(ifit_err_0,2.652e15)

# =============================================================================
# Analyse spectra in second waveband
# =============================================================================

print('Analyse synthetic spectra in waveband 1')


# Create second ifit analyser for results 
analyser1 = Analyser(params,
                     fit_window = wav1,
                     frs_path   = 'Ref/sao2010.txt',
                     flat_flag  = False,
                     stray_flag = False,
                     dark_flag  = False,
                     ils_type   = 'Manual')

# Create arrays to store answers
ifit_so2_1 = np.zeros((shape[0],shape[2]))
ifit_err_1 = np.zeros((shape[0],shape[2]))

# Loop through each synthetic spectrum
for i, so2 in enumerate(so2_grid):
    for j, ldf in enumerate(ldf_grid):
        
        #Extract syntheteic spectrum for suite of spectra
        spectrum = [fit_prm.grid,spectra_suite[i,...,j]]
        
        fit = analyser1.fit_spectrum(spectrum, calc_od=['SO2', 'Ring', 'O3'])
        
        # Store SO2 fit parameters in array
        ifit_so2_1[i,j] = fit.params['SO2'].fit_val
        ifit_err_1[i,j] = fit.params['SO2'].fit_err
        
    print(i)
 
#Create new ifit_so2 with units in ppm.m
ifit_so2_ppmm1 = np.divide(ifit_so2_1,2.652e15)
ifit_err_ppmm1 = np.divide(ifit_err_1,2.652e15)

# =============================================================================
# Create plot to show light dilution curves
# =============================================================================
   
plt.figure()
#Create comparison curves for my data
for i in range(len(ldf_grid)):
   
    #Define the current light dilution being viewed
    ldf = ldf_grid[i]
    
    #Grab the corresponding rows of fitted so2 amounts
    waveband0 = ifit_so2_ppmm0[...,i]
    waveband1 = ifit_so2_ppmm1[...,i]
      
    plt.plot(waveband0,waveband1,alpha = 0.5)#,c='C0')
    
    plt.xlim(-50,1050)
    plt.ylim(-50,1050)
    plt.xlabel('Fitted SO$_2$ between 306-316nm (ppm.m)')
    plt.ylabel('Fitted SO$_2$ between 312-322nm (ppm.m)')
    
    plt.legend(title = 'Spectrometer')
