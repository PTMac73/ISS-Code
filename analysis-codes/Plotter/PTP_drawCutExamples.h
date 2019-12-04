// PTP_drawCutExamples.h
// Draws examples of the cuts
// ============================================================================================= //
// Patrick MacGregor
// Nuclear Physics Research Group
// School of Physics and Astronomy
// The University of Manchester
// LAST EDITED: 21/03/19
// ============================================================================================= //

#ifndef PTP_DRAW_CUT_EXAMPLES_H_
#define PTP_DRAW_CUT_EXAMPLES_H_

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
// Draw the recoil detector cuts
void drawCutExamples( TTree *t, plotterOptions &opt_s, TCutG* cut ){
	// Print welcome message
	printDiv(); printf("PLOTTING EXAMPLES OF CUTS\n"); printDiv();


	// Define stuff
	TCanvas *c_recoil; TH2F *h_recoil;
	TCanvas *c_timing; TH1F *h_timing0, *h_timing1;
	TCanvas *c_pos; TH1F *h_pos0, *h_pos1;
	THStack *hs_timing, *hs_pos;

	// RECOIL CUT
	gStyle->SetOptStat(kFALSE);
	SetPadMargins( ptm_style, 2 );
	c_recoil = new TCanvas( "c_recoil", "Recoil Cuts", C_WIDTH, C_HEIGHT );
	t->Draw("rdt[0]:rdt[4]>>h_recoil(500, 2000, 7000, 350, 1000, 3500)", "", "colz");
	h_recoil = (TH2F*)gDirectory->Get("h_recoil");
	if ( SWITCH_PLOT_DETAILS == 1 ){
		h_recoil->SetTitle("Recoil Detector Cuts");
	}
	else{
		h_recoil->SetTitle("");
	}
	h_recoil->GetXaxis()->SetTitle("dE");
	h_recoil->GetXaxis()->SetNdivisions(505);
	h_recoil->GetXaxis()->SetLabelOffset(0.015);
	h_recoil->GetYaxis()->SetTitle("E");

	// Print the plot if desired
	if ( SWITCH_PRINT_CANVAS == 1 ){
		c_recoil->Print( Form( "%s/cutRecoils0%s", opt_s.printDir.Data(), PRINT_FORMAT.Data() ), "EmbedFonts" );
	}

	// Now draw the cut
	cut->SetLineColor(1);
	cut->SetLineWidth(2);
	cut->Draw("same");

	// Print the plot if desired
	if ( SWITCH_PRINT_CANVAS == 1 ){
		c_recoil->Print( Form( "%s/cutRecoils1%s", opt_s.printDir.Data(), PRINT_FORMAT.Data() ), "EmbedFonts" );
	}

	// TIMING CUT
	SetPadMargins( ptm_style );
	c_timing = new TCanvas( "c_timing", "Timing Cuts", C_WIDTH, C_HEIGHT );
	Int_t cutNumber = 17;
	t->Draw( Form( "td_rdt_e[%i]>>h_timing0(100, -50, 50)", cutNumber ) );
	t->Draw( Form( "td_rdt_e[%i]>>h_timing1(100, -50, 50)", cutNumber ), Form( "td_rdt_e[%i] > td_rdt_e_cuts[%i][0] && td_rdt_e[%i] < td_rdt_e_cuts[%i][1]", cutNumber, cutNumber, cutNumber, cutNumber ) );
	h_timing0 = (TH1F*)gDirectory->Get("h_timing0");
	h_timing1 = (TH1F*)gDirectory->Get("h_timing1");
	h_timing0->SetFillColor(kRed);
	h_timing1->SetFillColor(kAzure);
	if ( SWITCH_PLOT_DETAILS == 1 ){
		hs_timing = new THStack("hs_timing","Timing Cut");
		gStyle->SetOptStat(kTRUE);
	}
	else{
		hs_timing = new THStack("hs_timing","");
	}
	hs_timing->Add(h_timing0);
	hs_timing->Add(h_timing1);
	hs_timing->Draw("nostack");
	hs_timing->GetYaxis()->SetTitle("#");
	hs_timing->GetXaxis()->SetTitle("t / 10^{-8} s");
	c_timing->Modified(); c_timing->Update();
	
	if ( SWITCH_PRINT_CANVAS == 1 ){
		c_timing->Print( Form( "%s/cutTiming%s", opt_s.printDir.Data(), PRINT_FORMAT.Data() ), "EmbedFonts" );
	}


	// POSITION CUT
	c_pos = new TCanvas( "c_pos", "Position Cuts", C_WIDTH, C_HEIGHT );
	SetPadMargins( ptm_style );
	Int_t cutNumber2 = 16;
	
	t->Draw( Form( "xcal[%i]>>h_pos0(200, -0.5, 1.5)", cutNumber2 ), "(cut0 || cut1 || cut2 || cut3)" );
	t->Draw( Form( "xcal[%i]>>h_pos1(200, -0.5, 1.5)", cutNumber2 ), Form( "(cut0 || cut1 || cut2 || cut3) && xcal[%i] > xcal_cuts[%i][0] && xcal[%i] < xcal_cuts[%i][1]", cutNumber2, cutNumber2, cutNumber2, cutNumber2 ) );
	h_pos0 = (TH1F*)gDirectory->Get("h_pos0");
	h_pos1 = (TH1F*)gDirectory->Get("h_pos1");
	h_pos0->SetFillColor(kRed);
	h_pos1->SetFillColor(kAzure);
	if ( SWITCH_PLOT_DETAILS == 1 ){
		hs_pos = new THStack("hs_pos","Position Cut");
	}
	else{
		hs_pos = new THStack("hs_pos","");
	}
	hs_pos->Add(h_pos0);
	hs_pos->Add(h_pos1);
	hs_pos->Draw("nostack");
	hs_pos->GetYaxis()->SetTitle("#");
	hs_pos->GetXaxis()->SetTitle("x");
	c_pos->Modified(); c_pos->Update();

	if ( SWITCH_PRINT_CANVAS == 1 ){
		c_pos->Print( Form( "%s/cutPosition%s", opt_s.printDir.Data(), PRINT_FORMAT.Data() ), "EmbedFonts" );
	}

}


#endif
