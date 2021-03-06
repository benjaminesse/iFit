Using iFit
###########

There are two main ways in which to use iFit. The first is to write your own program using ``Python`` by importing the relevant libraries. This is the most flexible method and allows for customisable analysis with minimal reliance on additional libraries.

The second is to use the Graphical User Interface, either by running ``iFitUI.py`` with Python or by using the executable from the latest release. More details on this can be found in the :ref:`gui` section.

Writing a Python Script
========================

iFit operates using two main python objects: :class:`~parameters.Parameters` and :class:`~spectral_analysis.Analyser`.

Parameters
-----------

A :class:`~parameters.Parameters` object contains all the information on the fit parameters used to analyse a spectrum. Each :class:`~parameters.Parameter` in the :class:`~parameters.Parameters` object has the following values:

* ``name`` identifies the :class:`~parameters.Parameter` (and so must be unique)

* ``value`` gives the initial guess for that :class:`~parameters.Parameter` in the fit

* ``vary`` controls whether that :class:`~parameters.Parameter` is allowed to be varied by the model

Additionally for gas cross-sections the ``xpath`` value is used to set the file path to the cross-section file.

So a :class:`~parameters.Parameters` object could be generated like this:

.. code-block:: python

  from ifit.parameters import Parameters

  # Create parameter dictionary
  params = Parameters()

  # Add the gases
  params.add('SO2',  value=1.0e16, vary=True, xpath='Ref/SO2.txt')
  params.add('O3',   value=1.0e19, vary=True, xpath='Ref/O3.txt')
  params.add('Ring', value=0.1,    vary=True, xpath='Ref/Ring.txt')

  # Add background polynomial parameters
  params.add('bg_poly0', value=0.0, vary=True)
  params.add('bg_poly1', value=0.0, vary=True)
  params.add('bg_poly2', value=0.0, vary=True)
  params.add('bg_poly3', value=1.0, vary=True)

  # Add intensity offset parameters
  params.add('offset0', value=0.0, vary=True)

  # Add wavelength shift parameters
  params.add('shift0', value=0.0, vary=True)
  params.add('shift1', value=0.1, vary=True)


This defines three gas :class:`~parameters.Parameter` objects for SO2, O3 and Ring, as well as the polynomial coefficients for the background polynomial, a constant intensity offset and a wavelength shift and squeeze. Not that the naming convention for the polynomial parameters (``bg_poly{n}``, ``offset{n}`` and ``shift{n}``) is fixed. Once the :class:`~parameters.Parameter` is defined the :class:`~spectral_analysis.Analyser` can be generated.

Analyser
---------

The :class:`~spectral_analysis.Analyser` handles the actual analysis of the spectra. It must be generated first, defining certain settings for the analysis, as well as the :class:`~parameters.Parameters` already defined:

.. code-block:: python

  from ifit.spectral_analysis import Analyser

  # Generate the analyser
  analyser = Analyser(params=params,
                      fit_window=[310, 320],
                      frs_path='Ref/sao2010.txt',
                      stray_flag=True,
                      stray_window=[280, 290])


This will generate an analyser that will fit the measured spectra between 310 - 320 nm, performing a stray light correction using the measured intensities between 280 - 290 nm.

Measured spectra can then be analysed by using :meth:`~spectral_analysis.Analyser.fit_spectrum`:

.. code-block:: python

  fit = analyser.fit_spectrum([x, y])

In this case x and y are the measured spectrum wavelengths and intensities respectively. This returns a :class:`~spectral_analysis.FitResult` object which holds the fit data and useful information.

The :class:`~spectral_analysis.FitResult` object contains a copy of the :class:`~parameters.Parameters` object that was passed to the :class:`~spectral_analysis.Analyser` with the ``fit_val`` and ``fit_err`` values populated with the optimised value and associated error for each :class:`~parameters.Parameter`. It also contains:

* ``grid`` the wavelength grid of the fit window

* ``spec`` the measured spectrum (after pre-processing) in the fit window

* ``fit`` the optimised model spectrum

* ``resid`` the residual between the measurement and the model

An example script is given in ``iFit.py``.
