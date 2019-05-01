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
TString input_file_dir = "/home/ptmac/Documents/07-CERN-ISS-Mg/analysis/analysis_codes/Doublet/doublet_data.dat";
const int NUM_DATA_POINTS = 12;
const int NUM_L = 4;

// SWITCHES
Bool_t SWITCH_DRAW_FITS = 0;
Bool_t SWITCH_BRUTE_FORCE = 1;

// Define storage containers
Double_t X[NUM_DATA_POINTS];
Double_t Y[NUM_DATA_POINTS];
Double_t E[NUM_DATA_POINTS];
Double_t PT[NUM_L][NUM_DATA_POINTS];

#endif
