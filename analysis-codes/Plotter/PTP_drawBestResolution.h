// PTP_drawBestResolution.h
// Header file containing the best resolution detectors
// ============================================================================================= //
// Patrick MacGregor
// Nuclear Physics Research Group
// School of Physics and Astronomy
// The University of Manchester
// ============================================================================================= //
#ifndef PTP_DRAW_BEST_RESOLUTION_H_
#define PTP_DRAW_BEST_RESOLUTION_H_

#include "PTPlotterINIT.h"
#include <TFile.h>
#include <TTree.h>
#include <TCutG.h>
#include <TString.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TH1.h>
#include <TSystem.h>
#include <iostream>
#include "WriteSPE.h"

void DrawBestResolution( TTree *t, plotterOptions &opt_s ){
	// Define the best detector numbers (based on width of peaks)
	const int NUM_DETS = 6;
	int best_detectors[NUM_DETS] = {0,	13,	8,	3,	10,	23};
	
	// Print welcome message
	printDiv(); printf("PLOTTING THE BEST RESOLUTION EXCITATION SPECTRUM\n" ); printDiv();
	
	// Draw the plot
	TString cut_str = opt_s.cutList[0];
	cut_str += " && ( ";
	for ( Int_t i = 0; i < NUM_DETS; i++ ){
		cut_str += Form( "detID == %d", best_detectors[i] );
		if ( i < NUM_DETS - 1 ){ cut_str += " || "; }
	}
	cut_str += " )";
	
	t->Draw("Ex>>h_BR(450,-1,8)", cut_str.Data(), "goff" );
	
	// Get the plot
	TH1F *h_BR = (TH1F*)gDirectory->Get("h_BR");
	
	// Format the histogram
	h_BR->SetTitle("Excitation Spectrum");
	h_BR->GetXaxis()->SetTitle("Excitation energy (MeV)");
	h_BR->GetYaxis()->SetTitle("Counts / 20 keV");
	gStyle->SetOptStat(0);
	
	// Draw the canvas
	TCanvas *c_BR = new TCanvas( "c_BR", "Best resolution excitation spectrum", C_WIDTH, C_HEIGHT );
	h_BR->Draw();
	
	if ( SWITCH_PRINT_CANVAS == 1 ){
		WriteSPE( h_BR->GetName(), Form( "BestResolution" ) );
		c_BR->Print( Form( "%sBestResolution%s", opt_s.printDir.Data(), PRINT_FORMAT.Data() ) );
		
		// Move the spe file to the PLOTS folder
		gSystem->Exec( Form( "mv BestResolution.spe %sBestResolution.spe", opt_s.printDir.Data() ) );
		
	}
}



#endif
