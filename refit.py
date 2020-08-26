# -*- coding: utf-8 -*-
"""
Created on Thu Jul  9 11:23:51 2020

Script to allow refitting of spectra using best guess of ldf from lookup table

@author: Matthew Varnam - The University of Manchester
@email : matthew.varnam(-at-)manchester.ac.uk
"""

import string

# Import numpy for numerical calculations
import numpy as np

# matplotlib produces high quality graphs and figures
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

# iFit reads and analyses UV spectra of SO2
from ifit.load_spectra import read_spectrum, average_spectra, read_scan
from ifit.parameters import Parameters
from ifit.spectral_analysis import Analyser

# iFit_mod is alterations to iFit to allow analysis of light diluted spectra
from ifit_mod.synthetic_suite import Analyser_ld

# Glob allows listing of all files in a particular directory
from glob import glob

# Use pandas for the Dataframe class that makes exporting results table easy
import pandas as pd

# =============================================================================
# Define analysis files
# =============================================================================
   
# Define primary directory containing all files
prm_dir = ('F:/Black_Hole/Data/201801_Masaya/Ryu_spectra/' +
           '20180115_DOAS_for_camera_all/')

# Read in lookup table results
df = pd.read_csv(prm_dir + 'lookup_results.csv')

#Set flat and stray light toggles
flat_bool = True
stray_bool = True

# Set spectrometer name and directory
spec_name = 'FLMS02101_2'
spec_fpath = prm_dir + ('dark/spectra_')
dark_fpath = prm_dir + ('dark/spectra_')

stray_range = [280,290]

# =============================================================================
# Define some filepaths
# =============================================================================

# Load file of best fit ldf and list all spectra files
array = np.array(df['LDF'])
ldf_best = array
past_list = np.full((20),0.3)
index = 0

for i,ldf in enumerate(array):
    if not np.isnan(ldf) and ldf < 0.7:          
        past_list[index] = ldf
        if index<19:
            index += 1
        else:
            index = 0
    else:
        ldf_best[i] = np.mean(past_list)

# Create list of all spectra filenames
spec_fnames = glob(spec_fpath + '*.txt')

# =============================================================================
# Model Setup
# =============================================================================

flat_fpath = ('Spectrometer/flat_' + spec_name + '.txt')

# Create parameter dictionary
params = Parameters()

# Add the gases
params.add('SO2',  value=1.0e17, vary=True, xpath='Ref/SO2_295K.txt')
params.add('O3',   value=1.0e19, vary=True, xpath='Ref/O3_233K.txt')
params.add('Ring', value=0.2,    vary=True, xpath='Ref/Ring.txt')
    
# Add background polynomial parameters
params.add('bg_poly0', value=0.0, vary=True)
params.add('bg_poly1', value=0.0, vary=True)
params.add('bg_poly2', value=0.0, vary=True)
params.add('bg_poly3', value=1.0, vary=True)

# Add wavelength shift parameters
params.add('shift0', value=+0.004, vary=True)
params.add('shift1', value=-2.0 , vary=True)

# Add ILS parameters
#params.add('fwem', value = 0.65 , vary = False)
#params.add('k',    value = 2.32 , vary = False)
#params.add('a_w',  value = 0.054 , vary = False)
#params.add('a_k',  value = 0.72 , vary = False)

params.add('fwem', value = 0.60 , vary = False)
params.add('k',    value = 3.18 , vary = False)
params.add('a_w',  value = 0.126, vary = False)
params.add('a_k',  value = 1.62 , vary = False)

# Add LD parameter
params.add('LDF', value = 0.0, vary = False)

# Generate the analysers
analyser0 = Analyser_ld(params,
                        fit_window = [306.01, 316.01],
                        frs_path   = 'Ref/sao2010.txt',
                        flat_flag  = flat_bool,
                        flat_path  = flat_fpath,
                        stray_flag = stray_bool,
                        stray_window=stray_range,
                        dark_flag  = False,
                        ils_type   = 'Manual')

analyser1 = Analyser_ld(params,
                        fit_window = [312, 322],
                        frs_path   = 'Ref/sao2010.txt',
                        flat_flag  = flat_bool,
                        flat_path  = flat_fpath,
                        stray_flag = stray_bool,
                        stray_window=stray_range,
                        dark_flag  = False,
                        ils_type   = 'Manual')

print(params.pretty_print())

print('Done!\n')

# =============================================================================
# Run analysis
# =============================================================================

# Read in the dark spectrum and asign it in the analyser
#x, dark = np.loadtxt(dark_fname, unpack=True)
x, dark = average_spectra(dark_fnames,spec_type = 'Spectrasuite')
analyser0.dark_spec = dark
analyser1.dark_spec = dark

# Create results dataframe
col_names = ('num','time','LDF','SO2_0','SO2_0_err','SO2_1','SO2_1_err')
n_spec = np.arange(len(spec_fnames))
results_df = pd.DataFrame(index = n_spec, columns = col_names)

# Iterate through spectra 
for i,spec_fname in enumerate(spec_fnames):
    
    ldf = ldf_best[i]
    
    if np.isnan(ldf):
        # Change analyser values
        analyser0.params['LDF'].set(value = 0)
        analyser1.params['LDF'].set(value = 0)
    
    else:
        # Change analyser values
        analyser0.params['LDF'].set(value = ldf)
        analyser1.params['LDF'].set(value = ldf)
    
    print(analyser0.params['LDF'].value)
    
    # Read in the spectrum
    x, y, spec_info, read_err = read_spectrum(spec_fname, 'Spectrasuite')
    
    fit0 = analyser0.fit_spectrum([x,y], calc_od=['SO2', 'Ring', 'O3'],
                                  update_params=False)
    
    fit1 = analyser1.fit_spectrum([x,y], calc_od=['SO2', 'Ring', 'O3'],
                                  update_params=False)
    
    row = [i, spec_info['time'],ldf]
    row += [fit0.params['SO2'].fit_val,fit0.params['SO2'].fit_err]
    row += [fit1.params['SO2'].fit_val,fit1.params['SO2'].fit_err]
    results_df.loc[i] = row 
    
    print(i)

#results_df.to_csv('refit_results.csv')