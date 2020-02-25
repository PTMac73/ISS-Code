// AT_Settings.h
// Contains settings and global variables for the AnalyseTree.C script
// ============================================================================================= //
// Patrick MacGregor
// Nuclear Physics Research Group
// Department of Physics and Astronomy
// The University of Manchester
// ============================================================================================= //
#ifndef AT_SETTINGS_H_
#define AT_SETTINGS_H_

#include <TString.h>
#include <iostream>

// CONSTANTS
// Canvas
const Int_t C_WIDTH = 1200;
const Int_t C_HEIGHT = 900;

// Directories
TString print_dir = "/home/ptmac/Documents/07-CERN-ISS-Mg/Mg-Analysis/SPE-Files";

// Select row number to look at (-1 means do them all)
const Int_t ROW_NUMBER = 3;

// Select angle cuts
const Double_t THETA_MIN = 11.0;
const Double_t THETA_LB = 17.0;
const Double_t THETA_UB = 30.0;

// PRINT OPTIONS
const Bool_t PRINT_PDF = 1;
const Bool_t PRINT_PNG = 0;
const Bool_t PRINT_ROOT = 1;
TFile* out_root_file;


// SWITCHES
/* Form is:
	(0) Make the histograms for this
	(1) Print the histograms (pdf, root etc)
	(2) Write SPE files
*/
const Bool_t  SW_EX_COMPARE[3] = { 1, 0, 0 };
const Bool_t    SW_RDT_CUTS[3] = { 0, 0, 0 };
const Bool_t SW_EVZ_COMPARE[3] = { 1, 0, 0 };
const Bool_t         SW_EVZ[3] = { 0, 0, 0 };

// GLOBAL VARIABLES
TObjArray* cut_list;

// Cut booleans
Bool_t is_in_used_det;
Bool_t is_in_rdt;
Bool_t is_in_td;
Bool_t is_in_xcal;
Bool_t is_in_theta_min;
Bool_t is_in_theta_range;
Bool_t is_in_row;




#endif
