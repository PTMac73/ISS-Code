# ISS-Code
Code used to analyse data from the [ISOLDE Solenoidal Spectrometer](https://isolde-solenoidal-spectrometer.web.cern.ch/). A number of different codes used here:


##analysis-codes
#### AlphaCalibration
Performs an alpha calibration on the data that you feed in. Used a quadruple alpha source, so selects for those peaks (pre-defined positions as working with a small data set), and then fits a series of Gaussian peaks for extracting centroids and doing an energy calibration.

#### Doublet
Numerically tries to fit a doublet with two distributions defined by two scaling factors.

#### ELUMSolidAngle
Simulation to calculate the ELUM solid angle by modelling trajectories for a range of CM angles.

#### ELUMYield
Extracts yields from the ELUM data by fitting a Gaussian peak on a quadratic background.

#### Plotter
\[OLD\] Plots quantities using custom formulae applied to a TTree. Superceded by analyse-tree.

#### analyse-tree
Runs a TSelector code through the TTree that encodes the ISS data, and selects and plots the elements that I want.

#### array-geometry
Simulates ejectile and residual nucleus trajectories within ISS.

#### cm-angle-calculator
Calculates the centre-of-mass angles within ISS for a given reaction.

#### energy-calibration
Takes a number of data points and fits a straight line through them for the energy calibration.

#### extract-yields
Fits the excitation spectrum with a series of fitted Gaussians.

#### fit-background-yields
\[OLD\] Fits the background using a TTree.

#### super-unbound-fit
\[OLD\] Specifically fits 3 very unbound states. Superceded by extract-yields.

#### unbound-doublet-fit
\[OLD\] Specifically fits an unbound doublet. Superceded by extract-yields.
