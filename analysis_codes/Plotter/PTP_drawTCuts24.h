// PTPlotter.h
// Draws the 24 timing cuts on one canvas
// ============================================================================================= //
// Patrick MacGregor
// Nuclear Physics Research Group
// School of Physics and Astronomy
// The University of Manchester
// LAST EDITED: 21/03/19
// ============================================================================================= //
#ifndef PTP_DRAW_TCUTS24_H_
#define PTP_DRAW_TCUTS24_H_

#include "PTPlotterINIT.h"
#include <TFile.h>
#include <TTree.h>
#include <TCutG.h>
#include <TString.h>
#include <TROOT.h>
#include <TSystem.h>
#include <THStack.h>
#include <iostream>

struct plotterOptions;

// Draw the 24 timing cuts between array and recoils
void drawTCuts24( TTree *t, plotterOptions &opt_s ){
	// Print welcome message
	printDiv(); printf("PLOTTING THE 24 TIMING CUT PLOTS\n"); printDiv();

	// Create a TCanvas with 24 slots to hold the xcal cut plots
	TCanvas *c1;
	TH1F *a[24], *b[24];
	c1 = new TCanvas( "canvas", "Timing cut plots", 1350, 900);
	c1->Divide(6,4);

	// Start populating the TCanvases
	for ( Int_t i = 0; i < 6; i++ ){
		for ( Int_t j = 0; j < 4; j++ ){
			// Change to the right pad
			c1->cd(6*j + i + 1);

			// Draw using the TTree, but don't draw graphically
			t->Draw( Form( "td_rdt_e[%i]>>a%i(60,-30,30)", 6*j + i, 6*j + i ), Form( "(cut0 || cut1 || cut2 || cut3) && td_rdt_e[%i] > td_rdt_e_cuts[%i][0] && td_rdt_e[%i] < td_rdt_e_cuts[%i][1]", 6*j + i, 6*j + i, 6*j + i, 6*j + i ), "goff" );		
			t->Draw( Form( "td_rdt_e[%i]>>b%i(60,-30,30)", 6*j + i, 6*j + i ), "( cut0 || cut1 || cut2 || cut3 )", "goff" );

			// Store the histogram
			a[6*j + i] = (TH1F*)gDirectory->Get( Form("a%i", 6*j + i ) );
			b[6*j + i] = (TH1F*)gDirectory->Get( Form("b%i", 6*j + i ) );

			// Draw the histogram
			b[6*j + i]->SetFillColor(2);
			a[6*j + i]->SetFillColor(5);
			a[6*j + i]->SetTitle( Form( "TCut[%i]", 6*j + i ) );
			b[6*j + i]->SetTitle( Form( "TCut[%i]", 6*j + i ) );

			b[6*j + i]->Draw();
			a[6*j + i]->Draw("SAME");

			// Print message
			printf( "Drawn plot %02i/24\n", 6*j + i );
		}
	}

	// Print the canvas
	if ( SWITCH_PRINT_CANVAS == 1 ){
		c1->Print( Form( "%sTD_rdt-e_24%s", opt_s.printDir.Data(), PRINT_FORMAT.Data() ) );
	}
}

#endif  // PTP_DRAW_TCUTS24_H_
