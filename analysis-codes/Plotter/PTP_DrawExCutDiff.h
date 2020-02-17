// PTP_DrawExCutDiff.h
// Draws the 6 excitation spectra with the cut difference between singles and coincidences
// ============================================================================================= //
// Patrick MacGregor
// Nuclear Physics Research Group
// Department of Physics and Astronomy
// The University of Manchester
// ============================================================================================= //
#ifndef PTP_DRAW_EX_CUT_DIFF_H_
#define PTP_DRAW_EX_CUT_DIFF_H_

#include "PTPlotterINIT.h"
#include "WriteSPE.h"

#include <TFile.h>
#include <TTree.h>
#include <THStack.h>
#include <TH1.h>

struct plotterOptions;

// This determines the detectors to use when making the excitation spectra.
Bool_t USE_DETECTOR[24] = {
	1,  // DET #0
	1,  // DET #1
	0,  // DET #2
	1,  // DET #3
	1,  // DET #4
	1,  // DET #5
	1,  // DET #6
	0,  // DET #7
	0,  // DET #8
	1,  // DET #9
	1,  // DET #10
	1,  // DET #11
	1,  // DET #12
	0,  // DET #13
	1,  // DET #14
	0,  // DET #15
	1,  // DET #16
	0,  // DET #17
	1,  // DET #18
	1,  // DET #19
	1,  // DET #20
	1,  // DET #21
	1,  // DET #22
	1  // DET #23
};

// Function to get the Det ID cut string
TString GetDetIDCutString( Int_t row_number ){
	TString cut_string = "";

	// Get rid of nonsense
	if ( row_number < 0 || row_number > 6 ){
		return cut_string;
	}

	// Add first bracket
	cut_string += "( ";

	// Do row-by-row
	if ( row_number < 6 ){
		// Add detector if it's needed
		for ( Int_t i = 0; i < 4; i++ ){
			if ( USE_DETECTOR[6*i + row_number] == 1 ){
				if ( cut_string != "( " ){ cut_string += " || "; }
				cut_string += Form( "detID == %i", 6*i + row_number );
			}
		}
	}

	// Do full excitation spectrum
	else{
		for ( Int_t i = 0; i < 24; i++ ){
			if( USE_DETECTOR[i] == 1 ){
				if ( cut_string != "( " ){ cut_string += " || "; }
				cut_string += Form( "detID == %i", i );
			}
		}
	}
	
	// Add final bracket
	cut_string += " )";

	// Return final cut string
	return cut_string;
}

void DrawExCutDiff( TTree *t, plotterOptions &opt_s ){
	// Print welcome message
	printDiv(); printf("PLOTTING THE DIFFERENCE BETWEEN EXCITATION SPECTRA\n"); printDiv();

	// Define pointer arrays to everything
	TH1D* h_all_cuts[7];
	TH1D* h_red_cuts[7];
	THStack* h_stack[7];
	TCanvas* c_row[7];

	// Loop over each of the detector rows
	for ( Int_t i = 0; i < 6; i++ ){

		// Make canvas
		c_row[i] = new TCanvas( Form( "c_row%i", i ), Form( "CANVAS ROW %i", i ), C_WIDTH, C_HEIGHT );
		std::cout << "About to draw hists\n";

		// Get pointers to histograms
		t->Draw( Form( "%s>>ex_all%i(450,-1,8)", excitation_mode[WHICH_EXCITATION].Data(), i ), Form( "%s", opt_s.cutList[0].Data() ), "goff" );
		h_all_cuts[i] = (TH1D*)gDirectory->Get( Form( "ex_all%i", i ) );
		std::cout << "Drawn 1 hist\n";
		t->Draw( Form( "%s>>ex_red%i(450,-1,8)", excitation_mode[WHICH_EXCITATION].Data(), i ), Form( "%s && %s", CUTS[1].Data(), CUTS[3].Data() ), "goff" );
		h_red_cuts[i] = (TH1D*)gDirectory->Get( Form( "ex_red%i", i ) );
		std::cout << "Drawn all hists\n";

		// Format the histograms
		

		// Now stack the histograms
		h_stack[i]->Add( h_all_cuts[i] );
		h_stack[i]->Add( h_red_cuts[i] );

		// Draw the histograms
		h_stack[i]->Draw();

		// Print the canvases
		if ( SWITCH_PRINT_CANVAS == 1 ){
			c_row[i]->Print( Form( "/home/ptmac/Desktop/Test/Plots/canvas%i.pdf", i ) );
		}
	
		std::cout << "Drawn stack " << i << "\n";
	}








}
#endif
