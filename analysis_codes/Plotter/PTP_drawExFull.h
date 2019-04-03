// PTPlotter.h
// Draws the full excitation spectrum for all of the rows on the array
// ============================================================================================= //
// Patrick MacGregor
// Nuclear Physics Research Group
// School of Physics and Astronomy
// The University of Manchester
// LAST EDITED: 21/03/19
// ============================================================================================= //
#ifndef PTP_DRAW_EXFULL_H_
#define PTP_DRAW_EXFULL_H_

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

// Draw the full excitation spectrum
void drawExFull( TTree *t, plotterOptions &opt_s ){
	// Print welcome message
	printDiv(); printf("PLOTTING THE FULL EXCITATION SPECTRUM\n"); printDiv();
	
	// Define variables
	TCanvas *cExFull[opt_s.numIter];
	TH1F * hExFull[opt_s.numIter];

	// Work out whether to iterate plots
	Int_t first_index;
	if ( SWITCH_ITERATE_PLOTS == 1 ){
		first_index = 0;
	}
	else{
		first_index = opt_s.numIter;
	}
	
	// Declare a plot title for the histograms
	TString plotTitle;

	// Iterate over the cuts applied (cumulatively)
	for ( Int_t i = first_index; i < opt_s.numIter; i++ ){
		// Draw the plot
		cExFull[i] = new TCanvas( Form( "Ex-Canvas%i", i ), Form( "Full Excitation Spectrum (%i/%i)", i+1, opt_s.numIter ), C_WIDTH, C_HEIGHT );
		t->Draw( Form( "Ex>>ex_%i(450, -1, 8)", i ), opt_s.cutList[opt_s.numIter - i - 1].Data() );
		
		// Store the plot
		hExFull[i] = (TH1F*)gDirectory->Get( Form( "ex_%i", i ) );

		// Format the plot
		plotTitle = "Ex with ";
		if (i == 0){plotTitle += opt_s.cutListDesc[opt_s.numIter - 1];}
		else{
			for ( Int_t j = 1; j < i + 1; j++ ){
				plotTitle += opt_s.cutListDesc[opt_s.numIter - j - 1].Data();
				if ( j != i ){plotTitle += ", ";}
				if ( j == i - 1 ){plotTitle += "and ";}
			}
		}
		plotTitle +=" cut";
		if (i != 1){ plotTitle +="s"; }

		// Set the plot title and axes names
		hExFull[i]->SetTitle( plotTitle.Data() );
		hExFull[i]->GetXaxis()->SetTitle("Ex / MeV");
		hExFull[i]->GetYaxis()->SetTitle("#");

		// Update the canvas
		cExFull[i]->Modified(); cExFull[i]->Update();
		
		// Print the plot if desired
		if ( SWITCH_PRINT_CANVAS == 1 ){
			cExFull[i]->Print( Form( "%sExFull_%i%s", opt_s.printDir.Data(), i, PRINT_FORMAT.Data() ) );
		}
		printf("===> PLOTTED SPECTRUM %i/%i\n\n", i + 1, opt_s.numIter);
	}


	// NOW PRINT COMPARITIVE PLOTS
	if (SWITCH_ITERATE_PLOTS == 1){
		// Define some new canvases
		TCanvas *cComp[opt_s.numIter-1];
		
		// Clone the histograms so that they can be changed
		TH1F* h_Comp[opt_s.numIter];
		for ( Int_t l = 0; l < opt_s.numIter; l++ ){
			h_Comp[l] = (TH1F*)hExFull[l]->Clone( Form( "h_Comp%i", l ) );
		}

		// Define a histogram stack
		THStack *hs[opt_s.numIter - 1];
		gDirectory->ls();

		TString cutString = "";
		for (Int_t k = 0; k < opt_s.numIter - 1; k++ ){
			// Define the canvas
			cComp[k] = new TCanvas( Form( "Ex-Comp-Canvas%i", k ), Form( "Comparitive Excitation Spectrum (%i/%i)", k + 1, opt_s.numIter - 1 ), C_WIDTH, C_HEIGHT );

		
			// Define the histogram stack (different based on title)
			if ( k == 0 ){
				hs[k] = new THStack( Form( "hs-%i", k ), Form( "Difference when adding %s cut", opt_s.cutListDesc[opt_s.numIter - k - 2].Data() ) );
			}
			else{
				// Append to the title string
				cutString += opt_s.cutListDesc[opt_s.numIter - k - 1].Data();

				// Define the histogram stack
				hs[k] = new THStack( Form( "hs-%i", k ), Form( "Difference between {%s} when adding %s cuts", cutString.Data(), opt_s.cutListDesc[opt_s.numIter - k - 2].Data() ) );
			}
			
			h_Comp[k]->SetFillColor(kAzure);
			h_Comp[k+1]->SetFillColor(kOrange);
			
			// Draw the two histograms
			hs[k]->Add(h_Comp[k]);		
			hs[k]->Add(h_Comp[k+1]);
			hs[k]->Draw("nostack");
			hs[k]->GetXaxis()->SetTitle("Ex / MeV");
			hs[k]->GetYaxis()->SetTitle("#");
			cComp[k]->Modified(); cComp[k]->Update();	
			
			// Print the plot if desired
			if ( SWITCH_PRINT_CANVAS == 1 ){
				cComp[k]->Print( Form( "%sExCompCuts_%i%s", opt_s.printDir.Data(), k, PRINT_FORMAT.Data() ) );
			}
			
			// Append a comma and space to cutString
			if( k > 0)cutString += ", ";

			// Print friendly message
			printf("===> PLOTTED COMPARITIVE SPECTRUM %i/%i\n\n", k + 1, opt_s.numIter - 1);
		}
	}
}

#endif  // PTP_DRAW_EXFULL_H_
