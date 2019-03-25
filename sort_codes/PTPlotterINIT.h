// PTPlotterINIT.h
// Initialises all the variables for the PTPlotter.C script
// ============================================================================================= //
// Patrick MacGregor
// Nuclear Physics Research Group
// School of Physics and Astronomy
// The University of Manchester
// LAST EDITED: 21/03/19
// ============================================================================================= //
#ifndef PTPLOTTER_INIT_H_
#define PTPLOTTER_INIT_H_

#include "PTMonitors.h"
#include <TFile.h>
#include <TTree.h>
#include <TCutG.h>
#include <TString.h>
#include <TROOT.h>
#include <TSystem.h>
#include <THStack.h>
#include <iostream>

// DEFINE GLOBAL VARIABLES
// Switches for functions
const Bool_t SWITCH_EVZ = 0;
const Bool_t SWITCH_XCAL_24 = 0;
const Bool_t SWITCH_TCUTS_24 = 0;
const Bool_t SWITCH_EX_6 = 1;
const Bool_t SWITCH_EX_FULL = 0;
const Bool_t SWITCH_CUT_EXAMPLES = 0;

// Global switches
const Bool_t SWITCH_ITERATE_PLOTS = 1;
const Bool_t SWITCH_DISPLAY_CANVAS = 0;
const Bool_t SWITCH_PRINT_CANVAS = 1;
const Int_t NUM_CUT_ITER = 5;

// Canvas size (for 4:3 screen)
const Int_t C_WIDTH = 1200;
const Int_t C_HEIGHT = 900;

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
		opt_s.numIter = 0;
	}
	
	// Directory for saving files
	opt_s.printDir = "../PLOTS/";
	
	// Style Options
	gStyle->SetOptStat("meni");
	gStyle->SetPalette(kRainBow);

	// Clear the terminal
	gSystem->Exec("clear");
	
}


// --------------------------------------------------------------------------------------------- //
// Draw a dividing line
void printDiv(){
	printf("================================================================================\n");
}
// --------------------------------------------------------------------------------------------- //
// // Gets the position number from a file name (character before the . in "XYZ-2.root")
Int_t getPosNumber( TString s ){
	TString file_ending = ".root";
	Int_t pos_number;

	// Remove the file suffix
	s.Remove( s.Length() - file_ending.Length(), file_ending.Length() );
	
	// Remove all but the last character
	s.Remove( 0, s.Length() - 1 );

	// Calculate the number
	pos_number = s.Atoi();

	// Return it
	return pos_number;
}


#endif  // PTPLOTTER_H_

