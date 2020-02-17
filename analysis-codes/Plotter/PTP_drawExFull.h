// PTPlotter.h
// Draws the full excitation spectrum for all of the rows on the array
// ============================================================================================= //
// Patrick MacGregor
// Nuclear Physics Research Group
// Department of Physics and Astronomy
// The University of Manchester
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
#include <TCanvas.h>
#include <TPad.h>
#include <iostream>

struct plotterOptions;

// Draw the full excitation spectrum
void DrawExFull( TTree *t, plotterOptions &opt_s ){
	// Print welcome message
	printDiv(); printf("PLOTTING THE FULL EXCITATION SPECTRUM\n"); printDiv();

	// Define variables
	TCanvas *c_ex_full[opt_s.numIter];
	TH1F* h_ex_full[opt_s.numIter];
	TPad* pad;
	
	// Declare a plot title for the histograms
	TString plot_title = "";
	
	// Iterate over the cuts applied (cumulatively)
	for ( Int_t i = 0; i < opt_s.numIter; i++ ){

		// Draw the plot
		c_ex_full[i] = new TCanvas( Form( "c_ex_full_%i", i ), Form( "FULL-EX-%i/%i", i+1, opt_s.numIter ), C_WIDTH, C_HEIGHT );

		pad = (TPad*)c_ex_full[i]->GetPad(0);
		pad->SetLeftMargin(0.08);
		pad->SetRightMargin(0.01);
		pad->SetTopMargin(0.07);
		pad->SetBottomMargin(0.09);

		t->Draw( Form( "%s>>h_ex_full_%i(450, -1, 8)", excitation_mode[WHICH_EXCITATION].Data(), i ), opt_s.cutList[i].Data() );

		// Store the plot
		h_ex_full[i] = (TH1F*)gDirectory->Get( Form( "h_ex_full_%i", i ) );

		// Format the plot
		if ( SWITCH_PLOT_DETAILS == 1 ){
			plot_title = opt_s.cutListDesc[i].Data();
		}

		// Set the plot title and axes names
		h_ex_full[i]->SetTitle( plot_title.Data() );
		h_ex_full[i]->GetXaxis()->SetTitle("Excitation Energy (MeV)");
		h_ex_full[i]->GetYaxis()->SetTitle("Counts per 20 keV");

		h_ex_full[i]->SetLineColor(kRed);

		// Update the canvas
		c_ex_full[i]->Modified(); c_ex_full[i]->Update();

		// Print the plot if desired
		if ( SWITCH_PRINT_CANVAS == 1 ){
			c_ex_full[i]->Print( Form( "%s/EX/%s_full_%i%s", opt_s.printDir.Data(), excitation_mode[WHICH_EXCITATION].Data(), i, PRINT_FORMAT.Data() ) );
		}
		printf("===> PLOTTED SPECTRUM %i/%i\n\n", i + 1, opt_s.numIter);
	}



	// NOW PRINT COMPARITIVE PLOTS
	if (SWITCH_COMPARE_PLOTS == 1){
		// Define some new canvases
		TCanvas *cComp[opt_s.numIter-1];
		
		// Clone the histograms so that they can be changed
		TH1F* h_Comp[opt_s.numIter];
		for ( Int_t l = 0; l < opt_s.numIter; l++ ){
			h_Comp[l] = (TH1F*)h_ex_full[l]->Clone( Form( "h_Comp%i", l ) );
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
				cComp[k]->Print( Form( "%s/EX/%sCompCuts_%i%s", opt_s.printDir.Data(), excitation_mode[WHICH_EXCITATION].Data(), k, PRINT_FORMAT.Data() ) );
			}
			
			// Append a comma and space to cutString
			if( k > 0)cutString += ", ";

			// Print friendly message
			printf("===> PLOTTED COMPARITIVE SPECTRUM %i/%i\n\n", k + 1, opt_s.numIter - 1);
		}
	}
}

#endif  // PTP_DRAW_EXFULL_H_
