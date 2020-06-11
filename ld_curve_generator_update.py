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
from ifit.spectral_analysis import Analyser

# iFit_mod is alterations to iFit to allow analysis of light diluted spectra
from ifit_mod.settings import import_settings 
from ifit_mod.synthetic_suite import Analyser_ld


import matplotlib.pyplot as plt

# =============================================================================
# User settings
# =============================================================================

# Define the directory containing all work files
mst_dir = ('C:/Work/PhD/spectra/20180115-Stationary/Ben_format_shifted/')

# Set spectrometer name and directory
spec_name = 'FLMS02101_2'
spec_dir  = mst_dir + ('spectra/spectrum_')
drk_dir   = mst_dir + ('dark/spectrum_')

# Create list of spectra numbers
spec_num = [511,521]
drk_num  = [  0,  0]

# Set model parameters
so2_grid_ppmm = np.arange(0,3100,100)
ldf_grid = np.arange(0,1.0,0.1)

# Set wavebands to view
pad  = 1
wav0 = [306,316]
wav1 = [312,322]

# =============================================================================
# Pre-processing
# =============================================================================

# List of names of all spectra files
spec_range = np.arange(spec_num[0],spec_num[1]+1)
spec_fnames = [spec_dir + "{:05d}".format(num) + '.txt'for num in spec_range]

# List of names of all dark files
drk_range  = np.arange(drk_num[0],drk_num[1]+1)
drk_fnames = [drk_dir + "{:05d}".format(num) + '.txt'for num in drk_range]

# Create wavelength range of large fit containing two sub wavebands
wav_mst  = [wav0[0]-pad, wav1[1]+pad]

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

# Load model parameters
common, settings = import_settings(spec_name,wav_mst)

# Find the fit and stray windows
common['fit_idx']   = np.where(np.logical_and(x > settings['w_lo'],
                                              x < settings['w_hi']))
common['stray_idx'] = np.where(np.logical_and(x > settings['s_lo'],
                                              x < settings['s_hi']))

# Store wide measurement grid (trimmed spectrometer grid) for later use
meas_grid = x[common['fit_idx']]

# If no stray light pixels available, turn off the flag
if len(common['stray_idx'][0]) == 0:
    common['stray_flag'] = False     
else:
    common['stray_flag'] = True

#Copy dark spectrum to common
common['dark'] = dark

# Create ifit analyser
analyser_mst = Analyser(common)

# =============================================================================
# Run broad wavelength analysis
# =============================================================================

fit = analyser_mst.fit_spectrum([x,y], calc_od=['SO2', 'Ring', 'O3'])

print(fit.params.pretty_print())

# =============================================================================
# Create suite of synthetic spectra
# =============================================================================

# Create second ifit analyser 
analyser_ld = Analyser_ld(common)
analyser_ld.params.update_values(fit.params.popt_list())
analyser_ld.params.add('LDF', value = 0.0, vary = True)

# Convert SO2 in ppmm to molecules/cm2
so2_grid = np.multiply(so2_grid_ppmm , 2.463e+15)

# Use shape of grid and so2 value to produce a single array
shape = (len(so2_grid),len(common['fit_idx'][0]),len(ldf_grid))
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
        spectra_suite[i,...,j] = analyser_ld.fwd_model(meas_grid,*fit_params)

# =============================================================================
# Analyse spectra in first waveband
# =============================================================================
        
print('Analyse synthetic spectra in waveband 0')

# Load model parameters
common, settings = import_settings(spec_name,wav0)

# Remove dark and stray from process as not included in synthetic spectra        
common['dark_flag'] = False
common['stray_flag']= False
common['flat_flag'] = False

# Define maximum and minimum index within synthetic measurement grid
common['fit_idx']   = np.where(np.logical_and(meas_grid > settings['w_lo'],
                                              meas_grid < settings['w_hi']))

# Create new analyser for 
analyser0 = Analyser(common)

# Create arrays to store answers
ifit_so2_0 = np.zeros((shape[0],shape[2]))
ifit_err_0 = np.zeros((shape[0],shape[2]))

# Loop through each synthetic spectrum
for i, so2 in enumerate(so2_grid):
    for j, ldf in enumerate(ldf_grid):

        #Extract syntheteic spectrum for suite of spectra
        spectrum = [meas_grid,spectra_suite[i,...,j]]
        
        # Analyse spectrum
        fit = analyser0.fit_spectrum(spectrum, calc_od=['SO2', 'Ring', 'O3'])
        
        # Store SO2 fit parameters in array
        ifit_so2_0[i,j] = fit.params['SO2'].fit_val
        ifit_err_0[i,j] = fit.params['SO2'].fit_err
        
    print(i)

#Create new ifit_so2 with units in ppm.m
ifit_so2_ppmm0 = np.divide(ifit_so2_0,2.463e15)

# =============================================================================
# Analyse spectra in second waveband
# =============================================================================

print('Analyse synthetic spectra in waveband 1')

# Load model parameters
common, settings = import_settings(spec_name,wav1)

# Remove dark and stray from process as not included in synthetic spectra        
common['dark_flag'] = False
common['stray_flag']= False
common['flat_flag'] = False

# Define maximum and minimum index within synthetic measurement grid
common['fit_idx']   = np.where(np.logical_and(meas_grid > settings['w_lo'],
                                              meas_grid < settings['w_hi']))

# Create new analyser for 
analyser1 = Analyser(common)

# Create arrays to store answers
ifit_so2_1 = np.zeros((shape[0],shape[2]))
ifit_err_1 = np.zeros((shape[0],shape[2]))

# Loop through each synthetic spectrum
for i, so2 in enumerate(so2_grid):
    for j, ldf in enumerate(ldf_grid):
        
        #Extract syntheteic spectrum for suite of spectra
        spectrum = [meas_grid,spectra_suite[i,...,j]]
        
        fit = analyser1.fit_spectrum(spectrum, calc_od=['SO2', 'Ring', 'O3'])
        
        # Store SO2 fit parameters in array
        ifit_so2_1[i,j] = fit.params['SO2'].fit_val
        ifit_err_1[i,j] = fit.params['SO2'].fit_err
        
    print(i)
 
#Create new ifit_so2 with units in ppm.m
ifit_so2_ppmm1 = np.divide(ifit_so2_1,2.463e15)

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