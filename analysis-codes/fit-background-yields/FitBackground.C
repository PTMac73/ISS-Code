// FitBackground.C
// Fits the background of the singles spectrum
// ============================================================================================= //
// Patrick MacGregor
// Nuclear Physics Research Group
// Department of Physics and Astronomy
// The University of Manchester
// ============================================================================================= //

#include <TFile.h>
#include <TTree.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TSpectrum.h>
#include <TCutG.h>
#include <iostream>


void FitBackground( TFile* f ){

	// Formatting options
	gStyle->SetOptStat(0);
	gStyle->SetPadLeftMargin( 0.1 );
	gStyle->SetPadRightMargin( 0.02 );
	gStyle->SetPadTopMargin( 0.02 );
	gStyle->SetPadBottomMargin( 0.1 );

	// Constants
	const Int_t ROW_NUM = 2;
	TString cut_start = "( cut0 || cut1 || cut2 || cut3 ) && td_rdt_e[] > td_rdt_e_cuts[][0] && td_rdt_e[] < td_rdt_e_cuts[][1]";
	TString cut_end = "xcal[] >= xcal_cuts[][0] && xcal[] <= xcal_cuts[][1] && ( xcal[] <= xcal_cuts[][2] || xcal[] >= xcal_cuts[][3]) && thetaCM >= 16.6278 && %s";
	TString cut_tot = cut_start + " && " + cut_end;

	TString det_id_str[6] = {
		"( detID == 00 || detID == 06 || detID == 12 || detID == 18 )",
		"( detID == 01 || detID == 19 )",
		"( detID == 14 || detID == 20 )",
		"( detID == 03 || detID == 09 || detID == 21 )",
		"( detID == 04 || detID == 10 || detID == 16 || detID == 22 )",
		"( detID == 05 || detID == 11 || detID == 23 )"
	};

	// Get the tree
	TTree* t = (TTree*)f->Get("fin_tree");

	// Get the singles spectrum and draw it on TCanvas
	const Int_t nbins = 450;
	const Double_t lb = -1.0;
	const Double_t ub = 8.0;
	TCanvas* c0 = new TCanvas( "c0", "CANVAS-00", 1200, 900 );
	t->Draw( Form( "Ex>>h_raw(%i,%f,%f)", nbins, lb, ub ), Form( cut_end, det_id_str[ROW_NUM].Data() ) );

	// Get the raw spectrum and format it
	TH1D* h_raw = (TH1D*)gDirectory->Get("h_raw");
	h_raw->SetTitle(";Ex (MeV);Counts");
	
	// Define a TSpectrum
	TSpectrum* s = new TSpectrum();

	// Define an array to hold the raw spectrum
	Double_t source[nbins];
	for ( Int_t i = 0; i < nbins; i++ ){
		source[i] = h_raw->GetBinContent(i+1);
	}
	
	// Get the background spectrum
	s->Background( source, nbins, 10, TSpectrum::kBackIncreasingWindow, TSpectrum::kBackOrder2, kTRUE,  TSpectrum::kBackSmoothing3, kFALSE );

	// Define a background histogram and set the values
	TH1D* h_bg = new TH1D("h_bg", "", nbins, lb, ub);
	for ( Int_t i = 0; i < nbins; i++ ){
		h_bg->SetBinContent(i+1, source[i]);
	}
	h_bg->SetLineColor(kRed);
	h_bg->Draw("SAME L");

	TH1D* h_spec = (TH1D*)h_raw->Clone();
	h_spec->Add(h_bg, -1 );
	h_spec->SetLineWidth(2);
	h_spec->SetLineColor( kBlue );
	h_spec->Draw("SAME");

	c0->Modified(); c0->Update();

	// ----------------------------------------------------------------------------------------- //
	TCutG *cut = (TCutG*)f->Get("cut0");
	cut = (TCutG*)f->Get("cut1");
	cut = (TCutG*)f->Get("cut2");
	cut = (TCutG*)f->Get("cut3");
	
	// Draw the new spectrum on top of the old one
	TCanvas* c1 = new TCanvas( "c1", "CANVAS-01", 1200, 900 );
	c1->cd();	
	
	// Get the total (proper) spectrum
	t->Draw( Form( "Ex>>h_prop(%i,%f,%f)", nbins, lb, ub ), Form( cut_tot, det_id_str[ROW_NUM].Data() ) );
	TH1D* h_prop = (TH1D*)gDirectory->Get("h_prop");
	h_prop->SetLineColor(kRed);
	h_prop->SetLineWidth(2);
	h_prop->Draw();

	// Add the singles spectrum
	h_spec->Draw("SAME");
	
	// Set the maximum
	h_prop->GetYaxis()->SetRangeUser( 0, TMath::Max( h_spec->GetBinContent( h_spec->GetMaximumBin() ), h_prop->GetBinContent( h_prop->GetMaximumBin() ) ) );


	c1->Modified(); c1->Update();



}
