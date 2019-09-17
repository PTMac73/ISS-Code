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
#include "WriteSPE.h"
#include <iostream>

struct plotterOptions;

TString ID_str[6] = {
	"( detID == 6 || detID == 12 || detID == 18 )",
	"( detID == 7 || detID == 13 || detID == 19 )",
	"( detID == 14 || detID == 20 )",
	"( detID == 9 || detID == 15 )",
	"( detID == 10 )",
	"( detID == 5 || detID == 11 )"
};

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
		//t->Draw( Form( "Ex>>ex6_%i(450, -1, 8)", i ), Form( "%s && detID %% 6 == %i", opt_s.cutList[0].Data(), i ) );
		std::cout << "Canvas " << i+1 << " created" << std::endl;
		t->Draw( Form( "%s>>ex6_%i(450, -1, 8)", excitation_mode[WHICH_EXCITATION].Data(), i ), Form( "%s && %s", opt_s.cutList[0].Data(), ID_str[i].Data() ) );
		std::cout << "Histogram " << i+1 << " drawn" << std::endl;
		
		// Get the excitation spectra
		ex6[i] = (TH1F*)gDirectory->Get( Form( "ex6_%i", i ) );
		
		// Save the excitation spectra as pdf and spe file
		if ( SWITCH_PRINT_CANVAS == 1 ){
			TString file_name = Form( "%s%i_%i-6", excitation_mode[WHICH_EXCITATION].Data(), pos_number, i );
			WriteSPE( Form("ex6_%i", i ), Form( "%s", file_name.Data() ) );
			cEx6[i]->Print( Form( "%s%s", file_name.Data(), PRINT_FORMAT.Data() ) );

			// Move the spe file to the PLOTS folder
			gSystem->Exec( Form( "mv %s.spe %s/28MgDP_ReducedDetectors/%s.spe", file_name.Data(), opt_s.printDir.Data(), file_name.Data() ) );
			gSystem->Exec( Form( "mv %s%s %s/28MgDP_ReducedDetectors/%s%s", file_name.Data(), PRINT_FORMAT.Data(), opt_s.printDir.Data(), file_name.Data(), PRINT_FORMAT.Data() ) );
		}
		
		// Write a message
		printf("===> PLOTTED SPECTRUM %i\n\n", i + 1);
	}
	
}

#endif  // PTP_DRAW_EX6_H_
