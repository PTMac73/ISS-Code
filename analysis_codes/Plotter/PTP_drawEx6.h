// PTP_drawEx6.h
// Draws the 6 excitation spectra
// ============================================================================================= //
// Patrick MacGregor
// Nuclear Physics Research Group
// School of Physics and Astronomy
// The University of Manchester
// LAST EDITED: 21/03/19
// ============================================================================================= //
#ifndef PTP_DRAW_EX6_H_
#define PTP_DRAW_EX6_H_

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

// Draw (and save as SPE files) the 6 excitation spectra
void drawEx6( TTree *t, plotterOptions &opt_s, Int_t pos_number ){
	// Print welcome message
	printDiv(); printf("PLOTTING AND SAVING THE 6 EXCITATION SPECTRA\n"); printDiv();

	// Define some variables
	TH1F *ex6[6];
	TCanvas *cEx6[6];

	// Loop over each of the detector rows
	for ( Int_t i = 0; i < 6; i++ ){
		// Draw the excitation spectra
		cEx6[i] = new TCanvas( Form( "cEx6_%i", i ), Form( "Detector Row %i", i ), C_WIDTH, C_HEIGHT );
		t->Draw( Form( "Ex>>ex6_%i(450, -1, 8)", i ), Form( "%s && detID %% 6 == %i", opt_s.cutList[0].Data(), i ) );
		
		// Get the excitation spectra
		ex6[i] = (TH1F*)gDirectory->Get( Form( "ex6_%i", i ) );
		
		// Save the excitation spectra as pdf and spe file
		if ( SWITCH_PRINT_CANVAS == 1 ){
			TString fileName = Form( "%sEX%i_%i|6", opt_s.printDir.Data(), pos_number, i );
			writespe( Form("ex6_%i", i ), Form( "%s", fileName.Data() ) );
			cEx6[i]->Print( Form( "%s%s", fileName.Data(), PRINT_FORMAT.Data() ) );

			// Move the spe file to the PLOTS folder
			gSystem->Exec( Form( "mv %s.spe %s%s.spe", fileName.Data(), opt_s.printDir.Data(), fileName.Data() ) );
		}
		
		// Write a message
		printf("===> PLOTTED SPECTRUM %i\n\n", i + 1);
	}
	
}

#endif  // PTP_DRAW_EX6_H_
