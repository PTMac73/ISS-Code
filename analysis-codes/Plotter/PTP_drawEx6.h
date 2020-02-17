// PTP_drawh_ex.h
// Draws the 6 excitation spectra
// ============================================================================================= //
// Patrick MacGregor
// Nuclear Physics Research Group
// Department of Physics and Astronomy
// The University of Manchester
// ============================================================================================= //
#ifndef PTP_DRAW_h_ex_H_
#define PTP_DRAW_h_ex_H_

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
	TH1F* h_ex[6];
	TCanvas* c_ex[6];

	// Loop over each of the detector rows
	for ( Int_t i = 0; i < 6; i++ ){
		// Draw the excitation spectra
		c_ex[i] = new TCanvas( Form( "c_ex_%i", i ), Form( "Detector Row %i", i ), C_WIDTH, C_HEIGHT );
		std::cout << "Canvas " << i+1 << " created" << std::endl;
		t->Draw( Form( "%s>>h_ex_%i(450, -1, 8)", excitation_mode[WHICH_EXCITATION].Data(), i ), Form( "%s && %s", opt_s.cutList[0].Data(), ID_str[i].Data() ) );
		std::cout << "Histogram " << i+1 << " drawn" << std::endl;
		
		// Get the excitation spectra
		h_ex[i] = (TH1F*)gDirectory->Get( Form( "h_ex_%i", i ) );
		
		// Save the excitation spectra as pdf and spe file
		if ( SWITCH_PRINT_CANVAS == 1 ){
			TString file_name = Form( "%s%i_%i-6", excitation_mode[WHICH_EXCITATION].Data(), pos_number, i );
			c_ex[i]->Print( Form( "%s%s", file_name.Data(), PRINT_FORMAT.Data() ) );

			// Move the spe file to the PLOTS folder
			gSystem->Exec( Form( "mv %s%s %s/28MgDP_ReducedDetectors/%s%s", file_name.Data(), PRINT_FORMAT.Data(), opt_s.printDir.Data(), file_name.Data(), PRINT_FORMAT.Data() ) );
		}
		
		// Write a message
		printf("===> PLOTTED SPECTRUM %i\n\n", i + 1);
	}
	
}

#endif  // PTP_DRAW_h_ex_H_
