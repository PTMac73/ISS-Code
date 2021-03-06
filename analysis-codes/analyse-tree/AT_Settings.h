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
//TString print_dir = "/home/ptmac/Documents/07-CERN-ISS-Mg/Mg-Analysis/SPE-Files";
TString print_dir = "/home/ptmac/Documents/07-CERN-ISS-Mg/analysis/analysis-codes/extract-yields/print";
//TString cut_dir = "/home/ptmac/Documents/07-CERN-ISS-Mg/analysis/working/ALL-MgCuts3.root";
TString cut_dir = "/home/ptmac/Documents/07-CERN-ISS-Mg/analysis/analysis-codes/analyse-tree/mg_cuts.root";
TString cut_dir_new = "/home/ptmac/Documents/07-CERN-ISS-Mg/analysis/analysis-codes/analyse-tree/cuttlefish.root";
TString cut_dir_si = "/home/ptmac/Documents/07-CERN-ISS-Mg/analysis/analysis-codes/analyse-tree/si_cuts.root";

// Decide how to plot excitation spectra
const Bool_t ALL_ROWS = 1;
const Bool_t ROW_BY_ROW = 1;
const Bool_t DET_BY_DET = 1;

// Select row/det number to look at (-1 means do them all)
const Int_t DET_NUMBER = -1;
const Int_t ROW_NUMBER = -1;

// Array position - can only be one or two
const Int_t ARR_POSITION = 1;

// Select angle cuts
const Double_t THETA_MIN = 11.0;//19.21; //16.6278;
const Double_t THETA_LB = 12.5;
const Double_t THETA_UB = 20.0;

// PRINT OPTIONS
const Bool_t DISPLAY_CANVAS = 0;
const Bool_t PRINT_PDF = 0;
const Bool_t PRINT_PNG = 0;
const Bool_t PRINT_TEX = 0;
const Bool_t PRINT_ROOT = 0;
const Bool_t CANVAS_COMBINE = 0;

// Cut creator
const Bool_t DRAW_NEW_CUTS = 0;

// Other constant booleans


// SWITCHES
/* Form is:
	(0) Make the histograms for this
	(1) Print the histograms (pdf, root etc)
	(2) Write SPE files
*/
const Bool_t  SW_EX_COMPARE[3] = { 0, 1, 0 };
const Bool_t    SW_RDT_CUTS[3] = { 0, 1, 0 };
const Bool_t      SW_EVZ_SI[3] = { 0, 0, 0 };
const Bool_t          SW_EX[3] = { 1, 1, 1 };
const Bool_t       SW_EX_SI[3] = { 0, 1, 0 };
const Bool_t SW_EVZ_COMPARE[3] = { 0, 1, 0 };
const Bool_t         SW_EVZ[3] = { 0, 1, 0 };
const Bool_t        SW_XNXF[3] = { 0, 1, 0 };
const Bool_t        SW_XCAL[3] = { 0, 1, 0 };
const Bool_t          SW_TD[3] = { 0, 1, 0 };
const Bool_t     SW_SIGTIME[3] = { 0, 1, 0 };

// GLOBAL VARIABLES
TObjArray* cut_list;
TObjArray* cut_list_si;		// Si cuts to read in
TObjArray* cut_list_xnxf;
TObjArray* cuttlefish;		// Cuts to be read out

// Cut booleans
Bool_t is_in_used_det = 0;
Bool_t is_in_best_det = 0;

Bool_t is_in_rdt[4] = { 0, 0, 0, 0 };
Bool_t is_in_rdt_total = 0;

Bool_t is_in_rdt_si[4] = { 0, 0, 0, 0 };
Bool_t is_in_rdt_si_total = 0;

Bool_t is_in_td[4] = { 0, 0, 0, 0 };
Bool_t is_in_td_total = 0;

Bool_t is_in_TD[4] = { 0, 0, 0, 0 };			// Local version of time difference cuts
Bool_t is_in_TD_total = 0;

Bool_t is_in_rdt_and_td[4] = { 0, 0, 0, 0 };
Bool_t is_in_rdt_and_td_total = 0;

Bool_t is_in_rdt_si_and_td[4] = { 0, 0, 0, 0 };
Bool_t is_in_rdt_si_and_td_total = 0;

Bool_t is_in_rdt_and_TD[4] = { 0, 0, 0, 0 };
Bool_t is_in_rdt_and_TD_total = 0;

Bool_t is_in_xcal = 0;
//Bool_t is_in_xcal_mid = 0;
//Bool_t is_in_XCAL = 0;			// Local version of xcal cuts
Bool_t is_in_theta_min = 0;
Bool_t is_in_theta_range = 0;
Bool_t is_in_theta_custom = 0;
Bool_t is_in_theta_singles = 0;
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
	{ 1, 1, 1, 1 },
	{ 1, 1, 1, 1 },
	{ 1, 1, 1, 1 },
	{ 1, 1, 1, 1 },
	{ 1, 1, 1, 1 },
	{ 1, 0, 1, 1 }
};

Bool_t best_det_array[6][4] = {
	{ 0, 1, 0, 0 },
	{ 1, 0, 0, 0 },
	{ 1, 0, 0, 0 },
	{ 1, 0, 0, 0 },
	{ 0, 1, 0, 0 },
	{ 0, 0, 1, 0 }
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
	{ 0300.0, 0750.0, 0700.0, 1450.0 }, // 00
	{ 0650.0, 1200.0, 0100.0, 0850.0 }, // 01
	{ 0300.0, 1100.0, 0150.0, 1200.0 }, // 02
	{ 0450.0, 1100.0, 0150.0, 1000.0 }, // 03
	{ 0400.0, 1200.0, 0150.0, 1100.0 }, // 04
	{ 0350.0, 0800.0, 0350.0, 1100.0 }, // 05
	{ 0200.0, 1100.0, 0200.0, 1200.0 }, // 06
	{ 0250.0, 1100.0, 0200.0, 1200.0 }, // 07
	{ 0200.0, 0850.0, 0200.0, 1150.0 }, // 08
	{ 0250.0, 1150.0, 0150.0, 1250.0 }, // 09
	{ 0150.0, 1100.0, 0150.0, 1300.0 }, // 10
	{ 0200.0, 1800.0, 0200.0, 1800.0 }, // 11 - NOT USED
	{ 0150.0, 1150.0, 0150.0, 1400.0 }, // 12
	{ 0250.0, 1200.0, 0150.0, 1300.0 }, // 13
	{ 0200.0, 1150.0, 0100.0, 1200.0 }, // 14
	{ 0200.0, 1100.0, 0100.0, 1300.0 }, // 15
	{ 0300.0, 1150.0, 0150.0, 1200.0 }, // 16
	{ 0200.0, 1150.0, 0100.0, 1300.0 }, // 17
	{ 0500.0, 1150.0, 0150.0, 1100.0 }, // 18
	{ 0300.0, 1150.0, 0100.0, 1200.0 }, // 19
	{ 0300.0, 1050.0, 0150.0, 1200.0 }, // 20
	{ 0250.0, 1150.0, 0150.0, 1300.0 }, // 21
	{ 0250.0, 1050.0, 0150.0, 1200.0 }, // 22
	{ 0300.0, 0550.0, 0150.0, 0600.0 }  // 23
};

Double_t xnCorr[24] = {
	0.966727,	// 00
	0.997927,	// 01
	0.989599,	// 02
	0.926229,	// 03
	1.047585,	// 04
	0.909146,	// 05
	0.929986,	// 06
	0.967946,	// 07
	1.024587,	// 08
	0.962987,	// 09
	1.004322,	// 10
	0.000000,	// 11
	0.990265,	// 12
	1.001459,	// 13
	1.030865,	// 14
	0.971375,	// 15
	1.030609,	// 16
	0.925257,	// 17
	0.947191,	// 18
	1.003948,	// 19
	0.922519,	// 20
	0.992650,	// 21
	1.048827,	// 22
	1.099749,	// 23
};

Double_t thetaCM_cuts[6][2] = {
	{ 17.0, 17.0 },	// ROW 0
	{ 17.7, 17.7 },	// ROW 1
	{ 18.3, 18.8 },	// ROW 2
	{ 18.8, 19.4 },	// ROW 3
	//{ 19.4, 19.4 },	// ROW 4
	//{ 19.4, 19.4 },	// ROW 5
	{ 20.63, 20.63 }, 	// ROW 4
	{ 20.63, 20.63 } 	// ROW 5
};

// 0 is intercept, 1 is gradient
Double_t xfxneCorr[24][2] = {
	{   3.375253, 0.887861 },	// 00
	{   9.919062, 0.918479 },	// 01
	{   0.271875, 1.039087 },	// 02
	{   3.018420, 1.063643 },	// 03
	{   0.694993, 0.882137 },	// 04
	{  12.472202, 0.972475 },	// 05
	{   4.401876, 1.008029 },	// 06
	{   2.290693, 0.964587 },	// 07
	{   0.896536, 0.986746 },	// 08
	{   0.066819, 0.970263 },	// 09
	{  -2.399563, 0.909507 },	// 10
	{   0.000000, 0.000000 },	// 11
	{   0.699133, 0.783315 },	// 12
	{  -0.418725, 0.920520 },	// 13
	{  -0.859874, 0.933119 },	// 14
	{   0.475548, 0.984137 },	// 15
	{  -2.388185, 0.991268 },	// 16
	{   8.690010, 0.966768 },	// 17
	{   0.462283, 0.910253 },	// 18
	{  20.404622, 0.924948 },	// 19
	{   0.399083, 1.036339 },	// 20
	{   1.132734, 0.888264 },	// 21
	{  27.645136, 1.051020 },	// 22
	{  -1.097706, 0.874151 },	// 23
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

Double_t e_calibration_range = 32;

Double_t rawE_pos[24][4] = {
	{ 0815.0,  1325.0, 1400.0, 1480.0 }, // 00
	{ 0800.0,  1310.0, 1390.0, 1475.0 }, // 01
	{ 0875.0,  1425.0, 1525.0, 1610.0 }, // 02
	{ 0850.0,  1375.0, 1475.0, 1560.0 }, // 03
	{ 0800.0,  1290.0, 1375.0, 1450.0 }, // 04
	{ 0775.0,  1275.0, 1355.0, 1430.0 }, // 05
	{ 0835.0,  1350.0, 1440.0, 1525.0 }, // 06
	{ 0800.0,  1300.0, 1390.0, 1465.0 }, // 07
	{ 0830.0,  1355.0, 1420.0, 1520.0 }, // 08
	{ 0810.0,  1315.0, 1400.0, 1480.0 }, // 09
	{ 0760.0,  1240.0, 1310.0, 1390.0 }, // 10
	{ 0800.0,  1275.0, 1350.0, 1450.0 }, // 11
	{ 0675.0,  1100.0, 1175.0, 1250.0 }, // 12
	{ 0800.0,  1285.0, 1375.0, 1450.0 }, // 13
	{ 0775.0,  1260.0, 1340.0, 1415.0 }, // 14
	{ 0800.0,  1300.0, 1390.0, 1465.0 }, // 15
	{ 0840.0,  1355.0, 1445.0, 1525.0 }, // 16
	{ 0815.0,  1325.0, 1410.0, 1490.0 }, // 17
	{ 0780.0,  1280.0, 1365.0, 1440.0 }, // 18
	{ 0815.0,  1330.0, 1420.0, 1505.0 }, // 19
	{ 0840.0,  1375.0, 1465.0, 1540.0 }, // 20
	{ 0775.0,  1275.0, 1360.0, 1440.0 }, // 21
	{ 0925.0,  1500.0, 1600.0, 1700.0 }, // 22
	{ 0760.0,  1240.0, 1310.0, 1390.0 }  // 23

};
Float_t XCAL_cuts[24][2] = {
	{ 0.01, 0.96 }, // 00
	{ 0.00, 0.96 }, // 01
	{ 0.04, 1.00 }, // 02
	{ 0.03, 0.97 }, // 03
	{ 0.03, 1.00 }, // 04
	{ 0.05, 1.00 }, // 05
	{ 0.13, 0.87 }, // 06
	{ 0.06, 0.89 }, // 07
	{ 0.00, 0.87 }, // 08
	{ 0.01, 0.99 }, // 09
	{ 0.05, 0.94 }, // 10
	{ 0.00, 1.00 }, // 11
	{ 0.00, 1.00 }, // 12
	{ 0.00, 0.85 }, // 13
	{ 0.08, 0.88 }, // 14
	{ 0.04, 0.89 }, // 15
	{ 0.00, 0.96 }, // 16
	{ 0.02, 1.00 }, // 17
	{ 0.05, 1.00 }, // 18
	{ 0.00, 0.96 }, // 19
	{ 0.04, 1.00 }, // 20
	{ 0.06, 1.00 }, // 21
	{ 0.06, 0.99 }, // 22
	{ 0.00, 0.95 }  // 23
};
/*
Float_t XCAL_cuts[24][2] = {
	{ 0.00, 1.00 }, // 00
	{ 0.00, 1.00 }, // 01
	{ 0.00, 1.00 }, // 02
	{ 0.00, 1.00 }, // 03
	{ 0.00, 1.00 }, // 04
	{ 0.00, 1.00 }, // 05
	{ 0.00, 1.00 }, // 06
	{ 0.00, 1.00 }, // 07
	{ 0.00, 1.00 }, // 08
	{ 0.00, 1.00 }, // 09
	{ 0.00, 1.00 }, // 10
	{ 0.00, 1.00 }, // 11
	{ 0.00, 1.00 }, // 12
	{ 0.00, 1.00 }, // 13
	{ 0.00, 1.00 }, // 14
	{ 0.00, 1.00 }, // 15
	{ 0.00, 1.00 }, // 16
	{ 0.00, 1.00 }, // 17
	{ 0.00, 1.00 }, // 18
	{ 0.00, 1.00 }, // 19
	{ 0.00, 1.00 }, // 20
	{ 0.00, 1.00 }, // 21
	{ 0.00, 1.00 }, // 22
	{ 0.00, 1.00 }  // 23
};*/

Int_t TD_rdt_e_cuts[24][2] = {
	{  -7,  4 },	// 00
	{ -11,  4 },	// 01
	{  -9,  4 },	// 02
	{ -14,  3 },	// 03
	{ -16,  4 },	// 04
	{ -20,  3 },	// 05
	{  -5,  6 },	// 06
	{ -10,  3 },	// 07
	{  -10,  5 },	// 08
	{ -13,  4 },	// 09
	{ -19,  4 },	// 10
	{-200,200 },	// 11
	{ -12,  5 },	// 12
	{  -7,  7 },	// 13
	{  -9,  7 },	// 14
	{  -8,  5 },	// 15
	{ -20,  5 },	// 16
	{ -20,  6 },	// 17
	{ -11,  6 },	// 18
	{ -15,  6 },	// 19
	{ -11,  8 },	// 20
	{ -18,  5 },	// 21
	{ -22,  6 },	// 22
	{ -16,  6 } 	// 23
};


Float_t thetaCM_singles_cuts[6][2][2] = {
	{ { 11.0, 14.5 }, { 11.0, 90.0 } },
	{ { 11.0, 17.0 }, { 11.0, 90.0 } },
	{ { 12.0, 19.0 }, { 14.0, 90.0 } },
	{ { 12.5, 20.0 }, { 16.5, 24.0 } },
	{ { 18.5, 23.0 }, { 19.0, 25.0 } },
	{ { 20.5, 24.5 }, { 21.0, 26.0 } }
};


Double_t gain_match_pars[2][2] = {
	{ 1.00756, 0.0137467 },
	{ 1.00611, 0.0150445 }
};






#endif
