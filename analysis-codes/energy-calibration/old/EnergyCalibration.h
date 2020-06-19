// EnergyCalibration.h
// Header for EnergyCalibration.C
// ============================================================================================= //
// Patrick MacGregor
// Nuclear Physics Research Group
// Department of Physics and Astronomy
// The University of Manchester
// ============================================================================================= //
#ifndef ENERGY_CALIBRATION_H_
#define ENERGY_CALIBRATION_H_

const Int_t NUM_POINTS = 5;
Int_t C_WIDTH = 1200;
Int_t C_HEIGHT = 900;

Double_t monitors_data[NUM_POINTS][2] = {
	{ 1.0702, 0.0034 },
	{ 1.4139, 0.0023 },
	{ 2.2325, 0.0102 },
	{ 2.4665, 0.0053 },
	{ 4.3224, 0.0044 }
};

// -1 if not used in the fitting
Double_t ensdf_data[NUM_POINTS][2] = {
	{ 1.0946, 0.0002 },
	{ 1.4307, 0.0004 },
	{ 2.2659, 0.0007 },
	{ 2.4999, 0.0010 },
	{ 4.2800, 0.0400 }
};






#endif
