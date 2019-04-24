// PTPlotter.C
// Creates the 6 spectra for the array, and also plots the iterative change in the plots for
// plots where appropriate
// ============================================================================================= //
// Patrick MacGregor
// Nuclear Physics Research Group
// School of Physics and Astronomy
// The University of Manchester
// ============================================================================================= //
#include "PTPlotter.h"

#include <TFile.h>
#include <TTree.h>
#include <TCutG.h>
#include <TString.h>
#include <TROOT.h>

// MAIN FUNCTION =============================================================================== //
void PTPlotter( TFile *f ){
	// Declare options structure
	plotterOptions opt_s;
	
	// Initialise options
	initialiseOptions( opt_s );
	Int_t pos_number = GetPosNumber( (TString)f->GetName() );


	// Get the TTree and the TCutG's
	TTree *t = (TTree*)f->Get("fin_tree");
	TCutG *cut0 = (TCutG*)f->Get("cut0");
	TCutG *cut1 = (TCutG*)f->Get("cut1");
	TCutG *cut2 = (TCutG*)f->Get("cut2");
	TCutG *cut3 = (TCutG*)f->Get("cut3");

	// PLOT THE THINGS
	if ( SWITCH_EVZ == true )DrawEVZ( t, opt_s, pos_number );
	if ( SWITCH_XCAL_24 == true )drawXcal24( t, opt_s );
	if ( SWITCH_TCUTS_24 == true )drawTCuts24( t, opt_s );
	if ( SWITCH_EX_6 == true )drawEx6( t, opt_s, pos_number );
	if ( SWITCH_EX_FULL == true )drawExFull( t, opt_s );
	if ( SWITCH_CUT_EXAMPLES == true)drawCutExamples( t, opt_s, cut0 );
	if ( SWITCH_DRAW_BEST_RESOLUTION == true)DrawBestResolution( t, opt_s );
}

