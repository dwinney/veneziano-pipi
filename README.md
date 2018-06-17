# Pi-Pi Scattering

Toolkit for plotting and fitting pion-pion scattering amplitudes.

To install and use features in the core program simply clone and/or unpack and in main directory use ``make``.

Requires a working distribution of ROOT to both compile and use. Tested only with ROOT 6.12.06.

Written by Daniel Winney [contact: dwinney@iu.edu] for the JPAC Collaboration.

http://cgl.soic.indiana.edu/jpac/index.php

#### PLOTTING
Use ``./pipi plot (MODEL) (AMPLITUDE) (OPTIONS) `` to generate data sets of specified amplitudes as a function of energy and create a PDF plot.

Two MODELS are currently implemented:

 * ``GKPY`` - Based on the parameterization of [[1]](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.83.074004) based on a dispersive Roy-like description.
 * ``VENEZ`` - Based on a Veneziano amplitude [[2]](https://www.sciencedirect.com/science/article/pii/S0370269314006376). Extension of the Lovelace-Shapiro model, with linear (complex) Regge Trajectories.

For ``GKPY`` available amplitudes include :
 * ``isospin (0, 1, 2)`` - Produces the scattering amplitude of definite isospin.
 * ``total`` - Produces the total scattering amplitude for pion scattering.
 * ``partial (WAVE)`` - Produces desired partial-wave with definite angular momentum and isospin. The ``(WAVE)`` option should be input using the usual spectroscopic notation (S0, S2, P1, D0, D2, or F1).
 * ``phaseshift (WAVE)`` - Special to the GKPY parameterization, produces the phase shifts and inelasticities for a specified partial wave.

For ``VENEZ`` all amplitudes require a text file containing the coupling constants and linear Regge trajectory parameters denoted by ``(INPUT)`` which should be given as a relative path from the main directory.
 * ``isospin (0, 1 ,2) (INPUT)``- same as above.
 * ``total (INPUT)`` - self explanatory.
 * ``partial (WAVE) (INPUT)`` - note here as well ``(WAVE)`` should be a char string with spectroscopic notation.

Additional MODEL option, ``FILE``, to plot an already existing text file containing the real and imaginary parts of an amplitude as a function of energy.
To use: ``./pipi FILE (PATH TO FILE) (OUTPUT FILE NAME)``.

###### EXAMPLE
To produce the P1-wave Veneziano partial wave from couplings in ./example_coup.dat:
```
./pipi plot VENEZ partial P1 ./example_cout.dat
```

#### FITTING
UNDER CONSTRUCTION
