// AlphaDetByDet
// Plots the alpha calibration data detector by detector with the raw energy. Then writes to SPE
// file for analysis in gf3
// ============================================================================================= //
// Patrick MacGregor
// Nuclear Physics Research Group
// School of Physics and Astronomy
// The University of Manchester
// ============================================================================================= //
#ifndef ALPHA_DET_BY_DET_H_
#define ALPHA_DET_BY_DET_H_

#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TStyle.h>

#include "AlphaGlobals.h"
#include "WriteSPE.h"

void AlphaDetByDet( TTree *t ){
	
	// Declare variables
	const Int_t NUM_DETECTORS = 24;
	TH1F *h_det[NUM_DETECTORS];
	TCanvas *c_det[NUM_DETECTORS];

	// Loop over the detectors
	for ( Int_t i = 0; i < NUM_DETECTORS; i++ ){
		// Plot the detectors
		c_det[i] = new TCanvas( Form( "c_det%i", i ), Form( "Detector Canvas (%i/%i)", i + 1, NUM_DETECTORS ), 1200, 900 );
		t->Draw( Form( "e>>h_det%i(1000, 0, 2000)", i ), Form( "detID == %i && xcal[%i] > xcal_cuts[%i][0] && xcal[%i] < xcal_cuts[%i][1]", i, i, i, i, i ), "goff" );
		
		// Grab the histograms
		h_det[i] = (TH1F*)gDirectory->Get( Form ( "h_det%i", i ) );

		// Format and display the canvas (if option turned on)
		h_det[i]->SetTitle( Form( "Det #%i raw energy", i ) );
		h_det[i]->GetXaxis()->SetTitle("Raw Energy / arb.");
		h_det[i]->GetYaxis()->SetTitle("#");
		gStyle->SetOptStat("meni");
		h_det[i]->Draw();

		// Print the canvas if desired
		if ( SWITCH_PRINT_CANVAS == 1 ){
			c_det[i]->Print( Form( "%sAlpha_Raw_Energy_%i.pdf", PRINT_DIRECTORY.Data(), i ) );
		}

		// Write spectrum file if desired
		if ( SWITCH_WRITE_SPE == 1 ){
			WriteSPE( h_det[i]->GetName(), Form( "%sAlpha_Raw_Energy_%i", PRINT_DIRECTORY.Data(), i ) );
		}
		printf(" >>> Detector %02d complete\n\n", i );
	}
	
}

#endif
