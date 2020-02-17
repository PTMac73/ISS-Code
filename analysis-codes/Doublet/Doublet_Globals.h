// Doublet_Globals.h
// Contains all of the global variables
// ============================================================================================= //
// Patrick MacGregor
// Nuclear Physics Research Group
// School of Physics and Astronomy
// The University of Manchester
// ============================================================================================= //
#ifndef DOUBLET_GLOBAL_H
#define DOUBLET_GLOBAL_H

// GLOBAL VARIABLES
TString input_file_dir = "/home/ptmac/Documents/07-CERN-ISS-Mg/analysis/analysis-codes/Doublet/doublet_data.dat";
const int NUM_DATA_POINTS = 12;
const int NUM_L = 4;

// SWITCHES
const Bool_t SWITCH_DRAW_FITS = 0;
const Bool_t SWITCH_BRUTE_FORCE = 1;
const Bool_t SWITCH_VERBOSE = 0;

// Define storage containers
Double_t X[NUM_DATA_POINTS];
Double_t Y[NUM_DATA_POINTS];
Double_t E[NUM_DATA_POINTS];
Double_t PT[NUM_L][NUM_DATA_POINTS];

// Define spacing and upper and lower bounds (max and min - they cannot get higher or lower than these)
const Double_t STEP_SIZE_I = 1e-3;
const Double_t LB_0 = 0.0;
const Double_t UB_0 = 1.0;
const Double_t LB_1 = 0.0;
const Double_t UB_1 = 1.0;

#endif
