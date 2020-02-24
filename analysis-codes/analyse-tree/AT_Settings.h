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

#include <iostream>

// CONSTANTS
// Canvas
const Int_t C_WIDTH = 1200;
const Int_t C_HEIGHT = 900;

// Select row number to look at (-1 means do them all)
const Int_t ROW_NUMBER = 3;

// Select angle cuts
const Double_t THETA_MIN = 11.0;
const Double_t THETA_LB = 12.5;
const Double_t THETA_UB = 20.0;


// SWITCHES
/* Form is:
	(0) Make the histograms for this
	(1) Print the histograms
	(2) Write SPE files
*/
const Bool_t  SW_EX_COMPARE[3] = { 1, 0, 1 };
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
