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

/* TODO LIST
	* Change naming convention for files so it sorts well alphabetically:
		[HIST_NAME]-[POS]-[ROW]-[ADDITIONAL INFO].extension
	* Add more plots that will be useful with all the correct formatting - look at Plotter
	* Make it create all the spectra, superseding CreateMgSpectra
	* Make it do the theoretical lines on the E v.s. z plot.
	* Make a common axis labelling and histogram boundaries function for similar plots.
	* Add full excitation plots in last slot i.e. make the RBR things have TH1F* h[7], with h[6] = full spec.

*/

// CONSTANTS
// Canvas
const Int_t C_WIDTH = 1200;
const Int_t C_HEIGHT = 900;

// Directories
TString print_dir = "/home/ptmac/Documents/07-CERN-ISS-Mg/Mg-Analysis/SPE-Files";
TString cut_dir = "/home/ptmac/Documents/07-CERN-ISS-Mg/analysis/working/ALL-MgCuts3.root";
TString cut_dir_si = "/home/ptmac/Documents/07-CERN-ISS-Mg/analysis/analysis-codes/analyse-tree/cuttlefish.root";

// Select row number to look at (-1 means do them all)
const Int_t ROW_NUMBER = -1;

// Select angle cuts
const Double_t THETA_MIN = 11.0;
const Double_t THETA_LB = 12.0;
const Double_t THETA_UB = 19.0;

// PRINT OPTIONS
const Bool_t PRINT_PDF = 1;
const Bool_t PRINT_PNG = 0;
const Bool_t PRINT_ROOT = 0;
const Bool_t CANVAS_COMBINE = 0;

// Cut creator
const Bool_t DRAW_NEW_CUTS = 0;
TFile* out_root_file;

// SWITCHES
/* Form is:
	(0) Make the histograms for this
	(1) Print the histograms (pdf, root etc)
	(2) Write SPE files
*/
const Bool_t  SW_EX_COMPARE[3] = { 0, 0, 0 };
const Bool_t    SW_RDT_CUTS[3] = { 0, 0, 0 };
const Bool_t      SW_EVZ_SI[3] = { 1, 1, 0 };
const Bool_t       SW_EX_SI[3] = { 1, 1, 0 };
const Bool_t SW_EVZ_COMPARE[3] = { 0, 0, 0 };
const Bool_t         SW_EVZ[3] = { 0, 0, 0 };

// GLOBAL VARIABLES
TObjArray* cut_list;
TObjArray* cut_list_si;		// Si cuts to read in
TObjArray* cuttlefish;		// Cuts to be read out

// Cut booleans
Bool_t is_in_used_det = 0;
Bool_t is_in_rdt = 0;
Bool_t is_in_rdt_si = 0;
Bool_t is_in_td = 0;
Bool_t is_in_xcal = 0;
Bool_t is_in_theta_min = 0;
Bool_t is_in_theta_range = 0;
Bool_t is_in_row = 0;

// Other booleans
Bool_t found_si_cuts = 0;




#endif
