// PTP_drawEVZ.h
// Draw the E v.s. z plot
// ============================================================================================= //
// Patrick MacGregor
// Nuclear Physics Research Group
// School of Physics and Astronomy
// The University of Manchester
// ============================================================================================= //
#ifndef PTP_DRAW_EVZ_H_
#define PTP_DRAW_EVZ_H_

#include "PTPlotterINIT.h"
#include "PTF_GetPosNumber.h"
#include <TFile.h>
#include <TTree.h>
#include <TCutG.h>
#include <TString.h>
#include <TROOT.h>
#include <TSystem.h>
#include <THStack.h>
#include <TLine.h>
#include <TMath.h>
#include <iostream>

struct plotterOptions;

// --------------------------------------------------------------------------------------------- //
// Draw EVZ plot
void drawEVZ( TTree *t, plotterOptions &opt_s, Int_t pos_number ){
	gStyle->SetOptStat(kFALSE);
	
	// Print welcome message
	printDiv(); printf("PLOTTING THE E v.s. z PLOTS IN POSITION %i\n", pos_number ); printDiv();

	// Work out whether to iterate plots
	Int_t last_index;
	if ( SWITCH_ITERATE_PLOTS == 1 ){
		last_index = opt_s.numIter;
	}
	else{
		last_index = 1;
	}

	// Define variables
	TCanvas *c_evz[opt_s.numIter];
	TH2F *h_evz[opt_s.numIter];

	// Define axis limits
	Double_t x[4] = {-50,0,-10,9};


	// Loop to create graphs
	TString plotTitle;
	for ( Int_t i = 0; i < last_index; i++ ){
		// Call the canvas
		c_evz[i] = new TCanvas( Form( "c_evz%i", i ), Form( "E v.s. z (%i/%i)", i+1, opt_s.numIter - 1 ), C_WIDTH, C_HEIGHT );

		// Define the histogram
		t->Draw( Form( "ecrr:z>>evz_%i(400, %f, %f, 900, %f, %f)", i, x[0], x[2], x[1], x[3] ), opt_s.cutList[i], "colz" );

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

		// DRAW SI STRIP DIVIDING LINES
		if ( DRAW_SI_STRIP_DIVIDERS == 1 ){
			// Define dividers
			TLine *l_divider[6][2];
			for ( Int_t j = 0; j < 6; j++ ){
				for ( Int_t k = 0; k < 2; k++ ){
					// Define positions for each line
					l_divider[j][k] = new TLine();
					l_divider[j][k]->SetX1( TMath::Power(-1.0,k)*2.5 - Z_OFF[pos_number-1] - Z_ARRAY_POS[j] );
					l_divider[j][k]->SetX2( l_divider[j][k]->GetX1() );
					l_divider[j][k]->SetY1( x[1] );
					l_divider[j][k]->SetY2( x[3] );

					// Define styles
					l_divider[j][k]->SetLineWidth(2);
					if ( j % 2 == 0 ){
						l_divider[j][k]->SetLineColor(kBlack);
					}
					else{
						l_divider[j][k]->SetLineColor(kBlue);
					}
					// Draw the line
					l_divider[j][k]->Draw("same");
				}
			}		
		}

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
