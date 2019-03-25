// PTP_drawEVZ.h
// Draw the E v.s. z plot
// ============================================================================================= //
// Patrick MacGregor
// Nuclear Physics Research Group
// School of Physics and Astronomy
// The University of Manchester
// LAST EDITED: 21/03/19
// ============================================================================================= //
#ifndef PTP_DRAW_EVZ_H_
#define PTP_DRAW_EVZ_H_

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

// --------------------------------------------------------------------------------------------- //
// Draw EVZ plot
void drawEVZ( TTree *t, plotterOptions &opt_s ){
	gStyle->SetOptStat(kFALSE);
	
	// Print welcome message
	printDiv(); printf("PLOTTING THE E v.s. z PLOTS\n"); printDiv();

	// Work out whether to iterate plots
	Int_t first_index;
	if ( SWITCH_ITERATE_PLOTS == 1 ){
		first_index = 0;
	}
	else{
		first_index = opt_s.numIter - 1;
	}

	// Define variables
	TCanvas *c_evz[opt_s.numIter];
	TH2F *h_evz[opt_s.numIter];

	// Loop to create graphs
	TString plotTitle;
	for ( Int_t i = first_index; i < opt_s.numIter; i++ ){
		// Call the canvas
		c_evz[i] = new TCanvas( Form( "c_evz%i", i ), Form( "E v.s. z (%i/%i)", i+1, opt_s.numIter - 1 ), C_WIDTH, C_HEIGHT );

		// Define the histogram
		t->Draw( Form( "ecrr:z>>evz_%i(400, -50, -10, 900, 0, 9)", i), opt_s.cutList[i], "colz" );

		// Store the histogram
		h_evz[i] = (TH2F*)gDirectory->Get( Form( "evz_%i", i ) );

		// Format the histogram
		h_evz[i]->GetXaxis()->SetTitle("z / cm");
		h_evz[i]->GetYaxis()->SetTitle("E / MeV");
		
		// Produce title of histogram
		plotTitle = "E v.s. z with ";

		// Create the cut names
		if ( i == opt_s.numIter - 1 ){ plotTitle += opt_s.cutListDesc[i].Data(); }
		else{
			for ( Int_t j = i; j <= opt_s.numIter - 2; j++ ){
				plotTitle += opt_s.cutListDesc[j].Data();
				if ( j < opt_s.numIter - 2 ){ plotTitle += ", "; }
				if ( j == opt_s.numIter - 3 ){ plotTitle += "and "; }
			}
		}
		plotTitle +=" cut";
		if (i != opt_s.numIter - 2){ plotTitle +="s"; }
		printf("%i: %s\n", i, plotTitle.Data() );
		
		h_evz[i]->SetTitle( plotTitle.Data() );

		// Update the canvas
		c_evz[i]->Modified(); c_evz[i]->Update();

		// Print if desired
		// Print the canvas
		if ( SWITCH_PRINT_CANVAS == 1 ){
			c_evz[i]->Print( Form( "%sevz_%i%s", opt_s.printDir.Data(), i, PRINT_FORMAT.Data() ) );
		}
	}
	gStyle->SetOptStat(kTRUE);
	
}

#endif  // PTP_DRAW_EVZ_H_
