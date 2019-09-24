// PTPlotterINIT.h
// Initialises all the variables for the PTPlotter.C script
// ============================================================================================= //
// Patrick MacGregor
// Nuclear Physics Research Group
// School of Physics and Astronomy
// The University of Manchester
// ============================================================================================= //
#ifndef PTPLOTTER_INIT_H_
#define PTPLOTTER_INIT_H_

#include <TFile.h>
#include <TTree.h>
#include <TCutG.h>
#include <TString.h>
#include <TROOT.h>
#include <TSystem.h>
#include <THStack.h>
#include <TStyle.h>
#include <iostream>

// DEFINE GLOBAL VARIABLES
// Switches for functions
const Bool_t SWITCH_EVZ = 1;
const Bool_t SWITCH_GAMMA_BRANCH = 0;
const Bool_t SWITCH_XCAL_24 = 0;
const Bool_t SWITCH_TCUTS_24 = 0;
const Bool_t SWITCH_EX_6 = 0;
const Bool_t SWITCH_EX_FULL = 0;
const Bool_t SWITCH_CUT_EXAMPLES = 0;
const Bool_t SWITCH_DRAW_BEST_RESOLUTION = 0;

// Global switches
const Bool_t SWITCH_ITERATE_PLOTS = 0;			// Iterates plots over cuts
const Bool_t SWITCH_COMPARE_PLOTS = 0;			// Compares plots between cuts
const Bool_t SWITCH_DISPLAY_CANVAS = 1;			// Displays the canvas before the program ends
const Bool_t SWITCH_PRINT_CANVAS = 0;			// Prints to file
const Bool_t SWITCH_PLOT_DETAILS = 0;			// Adds OptStat box

// EVZ plotting options
const Bool_t DRAW_SI_STRIP_DIVIDERS = 1;
const Bool_t DRAW_THEORETICAL_LINES_NR = 1;

// Number of cut iterations (4 cuts = 5 different states)
const Int_t NUM_CUT_ITER = 5;

// Choose which excitation energy to use (0 = Ex, 1 = Ex_corrected)
const Int_t WHICH_EXCITATION = 1;
const TString excitation_mode[2] = { "Ex", "Ex_corrected" };

// Canvas size (for 4:3 screen)
const Int_t C_WIDTH = 1200;
const Int_t C_HEIGHT = 900;

// Z OFFSET POSITIONS
Double_t Z_OFF[2]= {9.498,6.5};
Double_t Z_ARRAY_POS[6] = {35.868, 29.987, 24.111, 18.248, 12.412, 6.676};

// Print format
const TString PRINT_FORMAT = ".pdf";

// Declare options struct
struct plotterOptions{
	TString printCanvasOption;
	TString cutList[NUM_CUT_ITER];
	TString cutListDesc[NUM_CUT_ITER];
	Int_t numIter;
	TString printDir;
};

// Declare global TStyle pointer
TStyle *ptm_style;

void SetPadMargins( TStyle *st, Int_t opt = 1, Int_t ratio = 43 ){
	if ( ratio == 11 ){
		if ( opt == 1 ){
			// 1D HIST
			st->SetPadLeftMargin(0.12);
			st->SetPadRightMargin(0.04);
		}
		else if ( opt == 2 ){
			// 2D HIST
			st->SetPadLeftMargin(0.12);
			st->SetPadRightMargin(0.12);
		}
	}
	else if ( ratio == 43 ){
		if ( opt == 1 ){
			// 1D HIST
			st->SetPadLeftMargin(0.08);
			st->SetPadRightMargin(0.02);
		}
		else if ( opt == 2 ){
			// 2D HIST
			//st->SetPadLeftMargin(0.12);
			//st->SetPadRightMargin(0.12);
		}
	}
	return;
}

Int_t GetRatio( Int_t width, Int_t height ){
	if ( width == height ){ return 11; }
	else if ( (Double_t)height/(Double_t)width == 0.75 ){ return 43; }
	else{ return 0; }
}


/* ********************************************************************************************* //
 * NOTES ON CUT ORDER
 * The order of the cuts is as follows:
 *    (i) Recoil cuts, to ensure you get Mg
 *   (ii) Theta cuts, to eliminate knees
 *  (iii) Timing cuts, to ensure a good timing resolution
 *   (iv) xcal cuts, to ensure good position resolution
 */


// FUNCTIONS =================================================================================== //
// Inititalise options
void initialiseOptions( plotterOptions &opt_s ){
	// CREATE OPTIONS FROM SWITCHES
	// Printing to a canvas
	if ( SWITCH_DISPLAY_CANVAS == 1 ){
		gROOT->SetBatch(kFALSE);
	}
	else{
		gROOT->SetBatch(kTRUE);
	}
	
	// Iterating cuts
	opt_s.cutList[4] = "";
	opt_s.cutList[3] = "( cut0 || cut1 || cut2 || cut3 )";
	opt_s.cutList[2] = opt_s.cutList[3] + " && thetaCM > 11";
	opt_s.cutList[1] = opt_s.cutList[2] + " && td_rdt_e[] > td_rdt_e_cuts[][0] && td_rdt_e[] < td_rdt_e_cuts[][1]";
	opt_s.cutList[0] = opt_s.cutList[1] + " && xcal[] > xcal_cuts[][0] && xcal[] < xcal_cuts[][1]";
	
	// Iterating cut descriptions
	opt_s.cutListDesc[4] = "no";
	opt_s.cutListDesc[3] = "Recoil";
	opt_s.cutListDesc[2] = "Angle";
	opt_s.cutListDesc[1] = "Timing";
	opt_s.cutListDesc[0] = "Position";

	// Iterate plots
	if ( SWITCH_ITERATE_PLOTS == 1 ){
		opt_s.numIter = NUM_CUT_ITER;
	}
	else{
		opt_s.numIter = 1;
	}
	
	// Directory for saving files
	opt_s.printDir = "/home/ptmac/Documents/07-CERN-ISS-Mg/analysis/PLOTS";
	
	ptm_style = (TStyle*)gStyle->Clone();
	ptm_style->SetName("ptm_style");
	
	// Style Options
	if ( SWITCH_PLOT_DETAILS == 1 ){
		ptm_style->SetOptStat("meni");
		ptm_style->SetPadTopMargin(0.08);
	}
	else{
		ptm_style->SetOptStat(0);
		ptm_style->SetPadTopMargin(0.02);
	}
	ptm_style->SetPadBottomMargin(0.08);
	
	// Set Pad margins for left and right
	SetPadMargins( ptm_style, 1, GetRatio( C_WIDTH, C_HEIGHT ) );
	
	ptm_style->SetTitleOffset(1.1,"y");
	ptm_style->SetTitleOffset(1.1,"xz");


	//ptm_style->SetPalette(kRainBow);
	ptm_style->SetLabelFont(62, "xyz");
	ptm_style->SetTitleFont(62, "xyz");
	ptm_style->SetTitleFont(62, "w");
	ptm_style->SetLegendFont(62);
	ptm_style->SetStatFont(62);
	
	gROOT->SetStyle("ptm_style");
	ptm_style->cd();

	// Clear the terminal
	gSystem->Exec("clear");
	
}


// --------------------------------------------------------------------------------------------- //
// Draw a dividing line
void printDiv(){
	printf("================================================================================\n");
}



#endif  // PTPLOTTER_H_

