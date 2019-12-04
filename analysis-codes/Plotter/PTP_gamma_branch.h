// PTP_gamma_branch.h
// Draw the gamma branch plots and saves them
// ============================================================================================= //
// Patrick MacGregor
// Nuclear Physics Research Group
// School of Physics and Astronomy
// The University of Manchester
// ============================================================================================= //
#ifndef PTP_GAMMA_BRANCH_H
#define PTP_GAMMA_BRANCH_H

#include "PTPlotterINIT.h"
#include "PTF_GetPosNumber.h"
#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>
#include <TStyle.h>
#include <TCutG.h>
#include <iostream>

void DrawGammaBranch( TTree *t, plotterOptions &opt_s, Int_t pos_number ){
	gStyle->SetOptStat(kFALSE);

	// Print welcome message
	printDiv(); printf("PLOTTING THE GAMMA BRANCH GRAPHS FOR POSITION %i\n", pos_number ); printDiv();
	
	// Draw the four graphs
	TCanvas *c_rdt = new TCanvas( "TITLE", "TITLE", C_WIDTH, C_HEIGHT );
	t->Draw("rdt[0]:rdt[4]>>h_rdt(300,0,9000,167,0,5000)", "");
	TH2F *h_rdt = (TH2F*)gDirectory->Get("h_rdt");
	gStyle->SetPalette(kRainBow);
	h_rdt->GetXaxis()->SetTitle("E");
	h_rdt->GetYaxis()->SetTitle("#Delta E");
	c_rdt->Modified(); c_rdt->Update();
	
}



#endif
