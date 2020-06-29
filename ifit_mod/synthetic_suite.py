# -*- coding: utf-8 -*-
"""
Created on Wed Jun 10 14:43:07 2020

@author: Matthew Varnam - The University of Manchester
@email : matthew.varnam(-at-)manchester.ac.uk
"""

import numpy as np
from scipy.interpolate import griddata

from ifit.make_ils import make_ils
from ifit.spectral_analysis import Analyser

class Analyser_ld(Analyser):
    def __init__(self, params, fit_window, frs_path, model_padding=1.0, 
                 model_spacing=0.01, flat_flag=False, flat_path=None,
                 stray_flag=False, stray_window=[280, 290], dark_flag=False,
                 ils_type='Manual', ils_path=None):
        
        Analyser.__init__(self, params, fit_window, frs_path, model_padding=1.0, 
                          model_spacing=0.01, flat_flag=False, flat_path=None,
                          stray_flag=False, stray_window=[280, 290], dark_flag=False,
                          ils_type='Manual', ils_path=None)
        
    def fwd_model(self, x, *p0):

        '''
        iFit forward model to fit measured UV sky spectra:
        I(w) = ILS *conv* {I_off(w) + I*(w) x P(w) x exp( SUM[-xsec(w) . amt])}
        where w is the wavelength.
        Requires the following to be defined in the common dictionary:
            - params:       Parameters object holding the fit parameters
            - model_grid:   The wavelength grid on which the forward model is
                            built
            - frs:          The Fraunhofer reference spectrum interpolated onto
                            the model_grid
            - xsecs:        Dictionary of the absorber cross sections that have
                            been pre-interpolated onto the model grid.
                            Typically includes all gas spectra and the Ring
                            spectrum
            - generate_ils: Boolian flag telling the function whether to
                            build the ILS or not. If False then the ILS
                            must be predefined in the common
            - ils           The instrument line shape of the spectrometer. Only
                            used if generate ILS is False.
        Parameters
        ----------
        grid, array
            Measurement wavelength grid
        *x0, list
            Forward model state vector. Should consist of:
                - bg_polyx: Background polynomial coefficients
                - offsetx:  The intensity offset polynomial coefficients
                - shiftx:   The wavelength shift polynomial
                - gases:    Any variable with an associated cross section,
                            including absorbing gases and Ring. Each "gas" is
                            converted to transmittance through:
                                      gas_T = exp(-xsec . amt)
                For polynomial parameters x represents ascending intergers
                starting from 0 which correspond to the decreasing power of
                that coefficient
        Returns
        -------
        fit, array
            Fitted spectrum interpolated onto the spectrometer wavelength grid
        '''

        # Get dictionary of fitted parameters
        params = self.params
        p = params.valuesdict()

        # Update the fitted parameter values with those supplied to the forward
        # model
        i = 0
        for par in params.values():
            if par.vary:
                p[par.name] = p0[i]
                i+=1
            else:
                p[par.name] = par.value

        # Unpack polynomial parameters
        bg_poly_coefs = [p[n] for n in p if 'bg_poly' in n]
        offset_coefs  = [p[n] for n in p if 'offset'  in n]
        shift_coefs   = [p[n] for n in p if 'shift'   in n]
        ldf_coefs     = [p[n] for n in p if 'LDF'     in n]

        # Construct background polynomial
        bg_poly = np.polyval(bg_poly_coefs, self.model_grid)
        frs = np.multiply(self.frs, bg_poly)

        # Create empty array to hold optical depth spectra
        gas_T = np.zeros((len(self.xsecs),
                          len(self.model_grid)))

        # Calculate the gas optical depth spectra
        for n, gas in enumerate(self.xsecs):
            if gas == 'SO2':
                so2_T  = (np.multiply(self.xsecs[gas], p[gas]))
            else:
                gas_T[n] = (np.multiply(self.xsecs[gas], p[gas]))
                
        # Calculate gas optical depth spectra
        plm_T = np.vstack((so2_T,gas_T))
        
        # Sum the gas ODs
        sum_gas_T = np.sum(gas_T, axis=0)
        sum_plm_T = np.sum(plm_T, axis=0)

        # Build the exponent term
        exponent     = np.exp(-sum_gas_T)
        exponent_plm = np.exp(-sum_plm_T)

        # Build the complete model
        dil_F = np.multiply(frs, exponent)
        plm_F = np.multiply(frs, exponent_plm)
        
        #Add wavelength dependancy to light dilution factor
        ldf_const = - np.log(1-ldf_coefs[0])*(310**4)
        rayleigh  = self.model_grid**-4
        ldf = 1-np.exp(-ldf_const*rayleigh)
        
        # Calc light dilution effect
        dil_F = np.multiply(dil_F, ldf)

        # Combine with normal spectrum
        plm_F = np.multiply(plm_F, 1-ldf)
            
        # Build the baseline polynomial
        offset = np.polyval(offset_coefs, self.model_grid)
        
        # Add diluted and undilted light
        raw_F = np.add(dil_F, plm_F) + offset
        
        # Generate the ILS
        if self.generate_ils:

            # Unpack ILS params
            ils = make_ils(self.model_spacing,
                           p['fwem'],
                           p['k'],
                           p['a_w'],
                           p['a_k'])
        else:
            ils = self.ils

        # Apply the ILS convolution
        F_conv = np.convolve(raw_F, ils, 'same')

        # Apply shift and stretch to the model_grid
        wl_shift = np.polyval(shift_coefs, self.model_grid)
        shift_model_grid = np.add(self.model_grid, wl_shift)

        # Interpolate onto measurement wavelength grid
        fit = griddata(shift_model_grid, F_conv, x,
                       method=self.interp_method)

        return fit
