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

// Decide how to plot excitation spectra
const Bool_t ALL_ROWS = 0;
const Bool_t ROW_BY_ROW = 1;
const Bool_t DET_BY_DET = 0;

// Select row/det number to look at (-1 means do them all)
const Int_t DET_NUMBER = -1;
const Int_t ROW_NUMBER = -1;


// Select angle cuts
const Double_t THETA_MIN = 11.0; //16.6278;
const Double_t THETA_LB = 11.0;
const Double_t THETA_UB = 14.5;

// PRINT OPTIONS
const Bool_t DISPLAY_CANVAS = 0;
const Bool_t PRINT_PDF = 1;
const Bool_t PRINT_PNG = 0;
const Bool_t PRINT_ROOT = 1;
const Bool_t CANVAS_COMBINE = 0;

// Cut creator
const Bool_t DRAW_NEW_CUTS = 0;
TFile* out_root_file;

// Other constant booleans


// SWITCHES
/* Form is:
	(0) Make the histograms for this
	(1) Print the histograms (pdf, root etc)
	(2) Write SPE files
*/
const Bool_t  SW_EX_COMPARE[3] = { 0, 1, 1 };
const Bool_t    SW_RDT_CUTS[3] = { 0, 0, 0 };
const Bool_t      SW_EVZ_SI[3] = { 0, 0, 0 };
const Bool_t          SW_EX[3] = { 0, 0, 0 };
const Bool_t       SW_EX_SI[3] = { 0, 0, 0 };
const Bool_t SW_EVZ_COMPARE[3] = { 0, 1, 0 };
const Bool_t         SW_EVZ[3] = { 0, 0, 0 };
const Bool_t        SW_XNXF[3] = { 1, 1, 0 };
const Bool_t        SW_XCAL[3] = { 0, 1, 0 };

// GLOBAL VARIABLES
TObjArray* cut_list;
TObjArray* cut_list_si;		// Si cuts to read in
TObjArray* cut_list_xnxf;
TObjArray* cuttlefish;		// Cuts to be read out

// Cut booleans
Bool_t is_in_used_det = 0;
Bool_t is_in_rdt = 0;
Bool_t is_in_rdt_si = 0;
Bool_t is_in_td = 0;
Bool_t is_in_xcal = 0;
Bool_t is_in_xcal_mid = 0;
Bool_t is_in_theta_min = 0;
Bool_t is_in_theta_range = 0;
Bool_t is_in_row = 0;
Bool_t is_in_xnxf_cut = 0;

// Other booleans
Bool_t found_si_cuts = 0;

// Detector array
/*Bool_t det_array[6][4] = {
	{ 1, 1, 1, 1 },
	{ 1, 0, 0, 1 },
	{ 0, 1, 1, 1 },
	{ 1, 1, 0, 1 },
	{ 1, 1, 1, 1 },
	{ 1, 1, 0, 1 }
};*/

Bool_t det_array[6][4] = {
	{ 0, 1, 0, 1 },
	{ 1, 1, 1, 1 },
	{ 1, 1, 1, 1 },
	{ 1, 1, 1, 1 },
	{ 1, 1, 1, 1 },
	{ 1, 0, 1, 1 }
};
/*
// Form is xlow, xup, ylow, yup ----OLD ONES BASED ON fin1.root
Double_t xnxf_lims[24][4] = {
	{ 0200.0, 0700.0, 0700.0, 1350.0 }, // 00
	{ 0650.0, 1200.0, 0060.0, 0910.0 }, // 01
	{ 0400.0, 1000.0, 0300.0, 1000.0 }, // 02
	{ 0500.0, 1000.0, 0200.0, 0800.0 }, // 03
	{ 0400.0, 1200.0, 0200.0, 1100.0 }, // 04
	{ 0300.0, 1100.0, 0200.0, 1100.0 }, // 05
	{ 0250.0, 0800.0, 0500.0, 1100.0 }, // 06
	{ 0250.0, 0750.0, 0500.0, 1100.0 }, // 07
	{ 0300.0, 0850.0, 0400.0, 0800.0 }, // 08
	{ 0850.0, 1250.0, 0250.0, 0650.0 }, // 09
	{ 0200.0, 0600.0, 0700.0, 1300.0 }, // 10
	{ 0000.0, 0000.0, 0000.0, 0000.0 }, // 11 - NOT USED
	{ 0000.0, 0000.0, 0000.0, 0000.0 }, // 12 - NOT USED
	{ 0250.0, 1200.0, 0200.0, 1200.0 }, // 13
	{ 0250.0, 1100.0, 0200.0, 1050.0 }, // 14
	{ 0250.0, 1150.0, 0200.0, 1150.0 }, // 15
	{ 0250.0, 1200.0, 0200.0, 1100.0 }, // 16
	{ 0200.0, 1150.0, 0200.0, 1200.0 }, // 17
	{ 0600.0, 1150.0, 0200.0, 0800.0 }, // 18
	{ 0250.0, 1200.0, 0200.0, 1200.0 }, // 19
	{ 0250.0, 1100.0, 0200.0, 1100.0 }, // 20
	{ 0250.0, 1250.0, 0300.0, 1200.0 }, // 21
	{ 0250.0, 1150.0, 0200.0, 1000.0 }, // 22
	{ 0200.0, 1100.0, 0300.0, 1000.0 }  // 23
};
*/
// Form is xlow, xup, ylow, yup
Double_t xnxf_lims[24][4] = {
	{ 0000.0, 2000.0, 0000.0, 2000.0 }, // 00 - NOT USED
	{ 0650.0, 1200.0, 0100.0, 0850.0 }, // 01
	{ 0300.0, 1100.0, 0150.0, 1200.0 }, // 02
	{ 0450.0, 1100.0, 0150.0, 1000.0 }, // 03
	{ 0400.0, 1200.0, 0150.0, 1100.0 }, // 04
	{ 0350.0, 0900.0, 0350.0, 1100.0 }, // 05
	{ 0200.0, 1100.0, 0200.0, 1200.0 }, // 06
	{ 0250.0, 1100.0, 0200.0, 1200.0 }, // 07
	{ 0200.0, 0850.0, 0200.0, 1150.0 }, // 08
	{ 0250.0, 1150.0, 0150.0, 1250.0 }, // 09
	{ 0150.0, 1150.0, 0200.0, 1200.0 }, // 10
	{ 0000.0, 0000.0, 0000.0, 0000.0 }, // 11 - NOT USED
	{ 0000.0, 0000.0, 0000.0, 0000.0 }, // 12 - NOT USED
	{ 0250.0, 1200.0, 0150.0, 1300.0 }, // 13
	{ 0200.0, 1150.0, 0200.0, 1300.0 }, // 14
	{ 0200.0, 1100.0, 0150.0, 1300.0 }, // 15
	{ 0200.0, 1150.0, 0200.0, 1300.0 }, // 16
	{ 0200.0, 1150.0, 0200.0, 1300.0 }, // 17
	{ 0450.0, 1200.0, 0150.0, 1100.0 }, // 18
	{ 0300.0, 1150.0, 0150.0, 1200.0 }, // 19
	{ 0200.0, 1100.0, 0150.0, 1300.0 }, // 20
	{ 0250.0, 1200.0, 0150.0, 1300.0 }, // 21
	{ 0200.0, 1150.0, 0200.0, 1300.0 }, // 22
	{ 0300.0, 0550.0, 0200.0, 0600.0 }  // 23
};

Double_t xnCorr[24] = {
	0.000000,	// 00
	0.997927,	// 01
	0.989599,	// 02
	0.926229,	// 03
	1.047585,	// 04
	0.950979,	// 05
	0.929986,	// 06
	0.967946,	// 07
	1.024587,	// 08
	0.962987,	// 09
	1.007415,	// 10
	0.000000,	// 11
	0.000000,	// 12
	1.001459,	// 13
	1.039200,	// 14
	0.971375,	// 15
	1.031104,	// 16
	0.927420,	// 17
	0.951143,	// 18
	1.003948,	// 19
	0.939985,	// 20
	0.991204,	// 21
	1.056439,	// 22
	1.099749	// 23
};

Double_t xfxneCorr[24][2] = {
	{   0.000000, 0.000000 },	// 00
	{  10.205986, 0.916083 },	// 01
	{   0.145663, 1.039639 },	// 02
	{   2.554639, 1.064789 },	// 03
	{   9.310831, 0.880631 },	// 04
	{  -3.812114, 0.957250 },	// 05
	{   1.975420, 1.011143 },	// 06
	{   3.681734, 0.961259 },	// 07
	{   1.522669, 0.982226 },	// 08
	{  -4.465927, 0.978012 },	// 09
	{  48.029516, 0.810101 },	// 10
	{   0.000000, 0.000000 },	// 11
	{   0.000000, 0.000000 },	// 12
	{  -0.498766, 0.919494 },	// 13
	{   1.744128, 0.926475 },	// 14
	{  -0.984655, 0.985022 },	// 15
	{  -2.544846, 0.990746 },	// 16
	{  11.467507, 0.956446 },	// 17
	{ -34.989082, 0.954787 },	// 18
	{   2.315599, 0.950286 },	// 19
	{  -0.466768, 1.032458 },	// 20
	{  17.498454, 0.866784 },	// 21
	{  21.096895, 1.044052 },	// 22
	{ -18.912714, 0.908218 }	// 23
};




Double_t xnxfE_lims[24][4] = {
	{ 0200.0, 1900.0, 0200.0, 1900.0 }, // 00
	{ 0200.0, 1900.0, 0200.0, 1900.0 }, // 01
	{ 0200.0, 1900.0, 0200.0, 1900.0 }, // 02
	{ 0200.0, 1900.0, 0200.0, 1900.0 }, // 03
	{ 0200.0, 1900.0, 0200.0, 1900.0 }, // 04
	{ 0200.0, 1900.0, 0200.0, 1900.0 }, // 05
	{ 0200.0, 1900.0, 0200.0, 1900.0 }, // 06
	{ 0200.0, 1900.0, 0200.0, 1900.0 }, // 07
	{ 0200.0, 1900.0, 0200.0, 1900.0 }, // 08
	{ 0200.0, 1900.0, 0200.0, 1900.0 }, // 09
	{ 0200.0, 1900.0, 0200.0, 1900.0 }, // 10
	{ 0200.0, 1900.0, 0200.0, 1900.0 }, // 11
	{ 0200.0, 1900.0, 0200.0, 1900.0 }, // 12
	{ 0200.0, 1900.0, 0200.0, 1900.0 }, // 13
	{ 0200.0, 1900.0, 0200.0, 1900.0 }, // 14
	{ 0200.0, 1900.0, 0200.0, 1900.0 }, // 15
	{ 0200.0, 1900.0, 0200.0, 1900.0 }, // 16
	{ 0200.0, 1900.0, 0200.0, 1900.0 }, // 17
	{ 0200.0, 1900.0, 0200.0, 1900.0 }, // 18
	{ 0200.0, 1900.0, 0200.0, 1900.0 }, // 19
	{ 0200.0, 1900.0, 0200.0, 1900.0 }, // 20
	{ 0200.0, 1900.0, 0200.0, 1900.0 }, // 21
	{ 0200.0, 1900.0, 0200.0, 1900.0 }, // 22
	{ 0200.0, 1900.0, 0200.0, 1900.0 }  // 23
};



#endif
