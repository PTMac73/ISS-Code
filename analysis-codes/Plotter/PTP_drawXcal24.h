// PTP_drawXcal24.h
// Draws the 24 xcal plots on one canvas
// ============================================================================================= //
// Patrick MacGregor
// Nuclear Physics Research Group
// School of Physics and Astronomy
// The University of Manchester
// LAST EDITED: 21/03/19
// ============================================================================================= //
#ifndef PTP_DRAW_XCAL24_H_
#define PTP_DRAW_XCAL24_H_

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

// Draw the 24 xcal plots
void drawXcal24( TTree *t, plotterOptions &opt_s ){
	// Print welcome message
	printDiv(); printf("PLOTTING THE 24 XCAL PLOTS\n"); printDiv();

	// Create variables
	TH1F *a[24], *b[24];

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
			t->Draw( Form( "xcal[%i]>>a%i(201,-0.5,1.5)", 6*j + i, 6*j + i ), "(cut0 || cut1 || cut2 || cut3) && thetaCM > 11 && td_rdt_e[] > td_rdt_e_cuts[][0] && td_rdt_e[] < td_rdt_e_cuts[][1] && xcal[] > xcal_cuts[][0] && xcal[] < xcal_cuts[][1]", "goff" );		
			t->Draw( Form( "xcal[%i]>>b%i(201,-0.5,1.5)", 6*j + i, 6*j + i ), "(cut0 || cut1 || cut2 || cut3) && thetaCM > 11 && td_rdt_e[] > td_rdt_e_cuts[][0] && td_rdt_e[] < td_rdt_e_cuts[][1]", "goff" );
			//t->Draw( Form( "xcal[%i]>>a%i(201,-0.5,1.5)", 6*j + i, 6*j + i ), "xcal[] > xcal_cuts[][0] && xcal[] < xcal_cuts[][1]", "goff" );		
			//t->Draw( Form( "xcal[%i]>>b%i(201,-0.5,1.5)", 6*j + i, 6*j + i ), "", "goff" );


			// Store the histogram
			a[6*j + i] = (TH1F*)gDirectory->Get( Form("a%i", 6*j + i ) );
			b[6*j + i] = (TH1F*)gDirectory->Get( Form("b%i", 6*j + i ) );

			// Draw the histogram
			b[6*j + i]->SetFillColor(2);
			a[6*j + i]->SetFillColor(5);
			//a[6*j + i]->SetTitle( Form( "xcal[%i]", 6*j + i ) );
			a[6*j + i]->SetTitle("");
			b[6*j + i]->SetTitle("");
	
			a[6*j + i]->GetXaxis()->SetTitle("xcal");
			a[6*j + i]->GetYaxis()->SetTitle("Counts");

			b[6*j + i]->Draw();
			a[6*j + i]->Draw("SAME");

			// Print message
			printf( "Drawn plot %02i/24\n", 6*j + i );

			if ( SWITCH_PRINT_CANVAS == 1 && SWITCH_INDIVIDUAL_CANVASES == 1 ){
				c1->Print( Form( "%s/xcal_det%i%s", opt_s.printDir.Data(), 6*j + i, PRINT_FORMAT.Data() ) );
			}
		}
	}

	// Print the canvas
	if ( SWITCH_PRINT_CANVAS == 1 && SWITCH_INDIVIDUAL_CANVASES == 0){
		c1->Print( Form( "%s/xcal24%s", opt_s.printDir.Data(), PRINT_FORMAT.Data() ) );
	}
}


#endif  // PTP_DRAW_XCAL24_H_
