// PTP_detector-by-detector.h
// Draws different things detector by detector
// ============================================================================================= //
// Patrick MacGregor
// Nuclear Physics Research Group
// School of Physics and Astronomy
// The University of Manchester
// ============================================================================================= //
#ifndef PTP_DRAW_DET_BY_DET_
#define PTP_DRAW_DET_BY_DET_

#include "PTPlotterINIT.h"
#include "PTP_DBD_classH24.h"
#include "PTP_DBD_SetType.h"
#include <TFile.h>
#include <TTree.h>
#include <TCutG.h>
#include <TString.h>
#include <TROOT.h>
#include <TSystem.h>
#include <THStack.h>
#include <iostream>

// GLOBALS


// FUNCTION - Draw all of the plots
void Draw24( TTree *t, plotterOptions &opt_s, Int_t type = 0 ){
	// Print welcome message
	printDiv(); printf("PLOTTING 24 PLOTS\n"); printDiv();

	// Define a global hist to be used
	H24 glob_h24 = SetType( type );

	// Create canvases
	TCanvas *c1;
	TCanvas *c_ind[24];

	if ( SWITCH_INDIVIDUAL_CANVASES == 0 ){
		// Create a TCanvas with 24 slots to hold the xcal cut plots
		c1 = new TCanvas( "canvas", "xcal plots", 1350, 900);
		c1->Divide(6,4);
	}

	// Start populating the TCanvases
	for ( Int_t i = 0; i < 6; i++ ){
		for ( Int_t j = 0; j < 4; j++ ){
			
			if ( SWITCH_INDIVIDUAL_CANVASES == 0 ){
				// Change to the right pad
				c1->cd(6*j + i + 1);
			}
			else{
				c_ind[6*j + i] = new TCanvas( Form( "c_det%i", 6*j + i ), Form( "XCAL PLOT DET #%i", 6*j + i ), C_WIDTH, C_HEIGHT );
				c1 = c_ind[6*j + i];
				c1->cd();
			}

			// Draw using the TTree, but don't draw graphically
			t->Draw( Form( glob_h24.MakeHistStr().Data(), 6*j+i ), Form( glob_h24.GetCutString(0), 6*j+i ), "goff" );

			// Get the histogram
			glob_h24.SetHFG( (TH1*)gDirectory->Get( Form( glob_h24.GetHistName().Data(), 6*j+i ) ), 6*j + i );


			// Do the same again if there is a background spectrum
			if ( glob_h24.GetBG() == true ){
				t->Draw( Form( glob_h24.MakeHistStr(1).Data(), 6*j+i ), Form( glob_h24.GetCutString(1), 6*j+i ), "goff" );
				glob_h24.SetHBG( (TH1*)gDirectory->Get( Form( glob_h24.GetHistName(1).Data(), 6*j+i ) ), 6*j + i );
			}

			// Format the histograms
			FormatHistograms( type, glob_h24.GetHFG( 6*j+i ), glob_h24.GetHBG( 6*j+i ) );

			// Draw the histograms
			if ( glob_h24.GetBG() == true && glob_h24.GetHBG( 6*j+i ) != NULL ){ glob_h24.GetHBG( 6*j+i )->Draw(""); }
			else{ std::cout << "BG is NULL\n"; }
			glob_h24.GetHFG( 6*j+i )->Draw("SAME");

			// Print message
			printf( "Drawn plot %02i/24\n", 6*j + i + 1 );

			if ( SWITCH_PRINT_CANVAS == 1 && SWITCH_INDIVIDUAL_CANVASES == 1 ){
				c1->Print( Form( "%s/%s/%s%s", opt_s.printDir.Data(), glob_h24.GetPrintFolder().Data(), Form( glob_h24.GetHistName().Data(), 6*j + i ), PRINT_FORMAT.Data() ) );
			}
			
		}
	}

	// Print the canvas
	if ( SWITCH_PRINT_CANVAS == 1 && SWITCH_INDIVIDUAL_CANVASES == 0){
		c1->Print( Form( "%s/xcal24%s", opt_s.printDir.Data(), PRINT_FORMAT.Data() ) );
	}
}


#endif  // PTP_DRAW_XCAL24_H_
