// MassVariations.h
// Functions for the mass variation stuff
// ============================================================================================= //
// Patrick MacGregor
// Nuclear Physics Research Group
// School of Physics and Astronomy
// The University of Manchester
// ============================================================================================= //
#ifndef MASS_VARIATIONS_H_
#define MASS_VARIATIONS_H_

#include <TAxis.h>
#include <TCanvas.h>
#include <TGraph.h>

// CONSTANTS
const Int_t C_WIDTH = 1200;
const Int_t C_HEIGHT = 900;

// Mass bases
const Int_t NUM_MASS = 3;
Double_t mass_base[NUM_MASS] = { 26066.8, 1876.12, 27002.7 };	// Define the base masses from which to start

// GLOBAL FUNCTIONS
// Set canvas margins
void GlobSetCanvasMargins( TCanvas *c, Double_t l = 0.1, Double_t r = 0.02, Double_t t = 0.02, Double_t b = 0.1 ){
	TPad* pad = (TPad*)c;
	pad->SetLeftMargin( l );
	pad->SetRightMargin( r );
	pad->SetTopMargin( t );
	pad->SetBottomMargin( b );
	return;
}

void GlobSetHistFonts( TGraph* g ){
	g->GetXaxis()->SetTitleFont(62);
	g->GetXaxis()->SetLabelFont(62);
	g->GetXaxis()->CenterTitle();
	g->GetYaxis()->SetTitleFont(62);
	g->GetYaxis()->SetLabelFont(62);
	g->GetYaxis()->CenterTitle();
	return;
}

void CreateMassGraph( TGraph* g, Int_t i ){
	GlobSetHistFonts(g);
	g->SetTitle("");
	g->GetXaxis()->SetTitle( Form( "m_{%i} (MeV)", i + 1 + ( i == 2 ? 1 : 0 ) ) );
	g->GetXaxis()->SetRangeUser( mass_base[i] -13*0.511, mass_base[i] + 0.511 );
	g->GetYaxis()->SetTitle( "#theta_{cm} (#circ)" );
	g->GetYaxis()->SetRangeUser(0,60);
	g->SetMarkerStyle(20);
	g->SetMarkerSize(0.1);
	g->SetMarkerColor(kRed);
	g->SetLineColor(kOrange);
	return;
}

void CreateMassDerivGraph( TGraph* g, Int_t i ){
	GlobSetHistFonts(g);
	g->SetTitle("");
	g->GetXaxis()->SetTitle( Form( "m_{%i} (MeV)", i + 1 + ( i == 2 ? 1 : 0 ) ) );
	g->GetXaxis()->SetRangeUser( mass_base[i] -13*0.511, mass_base[i] + 0.511 );
	g->GetYaxis()->SetTitle( Form( "#frac{d#theta_{cm}}{dm_{%i}} (#circ/MeV)", i + 1 + ( i == 2 ? 1 : 0 ) ) );
	g->GetYaxis()->SetRangeUser(-8,6);
	g->SetMarkerStyle(20);
	g->SetMarkerSize(0.1);
	g->SetMarkerColor(kRed);
	g->SetLineColor(kOrange);
	return;
}

#endif
