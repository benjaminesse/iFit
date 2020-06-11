# -*- coding: utf-8 -*-
"""
Created on Wed Jun 10 11:36:27 2020

An adaptation of a fitting example written by Ben Esse to load in preset values
for fitting spectra normally

@author: Matthew Varnam - The University of Manchester
@email : matthew.varnam(-at-)manchester.ac.uk
"""
# =============================================================================
# Import libraries
# =============================================================================

import numpy as np
from ifit.model_setup import model_setup
from ifit.parameters import Parameters

# =============================================================================
# Model settings
# =============================================================================

def import_settings(spec_name,wavelength):
   
    # Set up the program settings
    settings = {'w_lo':          wavelength[0],
                'w_hi':          wavelength[1],
                's_lo':          280,
                's_hi':          290,
                'model_spacing': 0.01,
                'model_padding': 1.0,
                'flat_flag':     True,
                'stray_flag':    True,
                'dark_flag':     True,
                'ils_type':      'Fit',
                'frs_path':      'Ref/sao2010.txt',
                'flat_path':     f'Spectrometer/flat_{spec_name}.txt',
                'ils_path':      f'Spectrometer/ils_params_{spec_name}.txt',
                'generate_ils':  True
                }
    
    # Set the gas data
    gas_data = {}
    gas_data['SO2']  = ['Ref/SO2_295K.txt', 1.0e16, True]
    gas_data['O3']   = ['Ref/O3_223K.txt',  1.0e19, True]
    gas_data['Ring'] = ['Ref/Ring.txt',     0.1,    True]
    #gas_data['NO2']  = ['Ref/NO2_223K.txt',      1.0e14, True]
    
    # Add to the settings dictionary
    settings['gas_data'] = gas_data
    
    # Build the forward model
    common = model_setup(settings)
    
    # =========================================================================
    # Parameter Setup
    # =========================================================================
    
    # Create parameter dictionary
    params = Parameters()
        
    # Set other model parameters
    params.add('bg_poly0', value = -1.0e-13, vary = True)
    params.add('bg_poly1', value =  1.0e-10, vary = True)
    params.add('bg_poly2', value = -1.0e-8, vary = True)
    params.add('bg_poly3', value =  1.0e-6, vary = True)
    # params.add('bg_poly4', value = 1.0, vary = True)
    # params.add('offset0',  value = 0.0, vary = True)
    # params.add('offset1',  value = 0.0, vary = True)
    # params.add('offset2',  value = 0.0, vary = True)
    params.add('shift0',   value = 0.0, vary = True)
    params.add('shift1',   value = -0.10, vary = True)
    
    # Add ILS parameters
    params.add('fwem', value = 0.66, vary = False)
    params.add('k',    value = 2.2, vary = False)
    params.add('a_w',  value = 0.06, vary = False)
    params.add('a_k',  value = 0.50, vary = False)
    
    #params.add('fwem', value = 0.65, vary = True)
    #params.add('k',    value = 2.0, vary = True)
    #params.add('a_w',  value = 0.0, vary = True)
    #params.add('a_k',  value = 0.0, vary = True)
    
    # Add the gases
    for gas in gas_data:
        params.add(gas, value = gas_data[gas][1], vary = gas_data[gas][2])
    
    # Add to the common
    common['params'] = params
    
    return common,settings