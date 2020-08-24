// AG_style.h
// TStyle definitions for the file ArrayGeometry.C
// ============================================================================================= //
// Patrick MacGregor
// Nuclear Physics Research Group
// Department of Physics and Astronomy
// The University of Manchester
// ============================================================================================= //
#ifndef AG_STYLE_H_
#define AG_STYLE_H_

#include <TCanvas.h>
#include <TColor.h>
#include <TH1.h>
#include <TH2.h>
#include <TStyle.h>
#include <TROOT.h>

// PRINTING OPTIONS
Bool_t PRINT_SVG = 1;
Bool_t PRINT_PDF = 1;
Bool_t PRINT_TEX = 1;
Bool_t PRINT_PNG = 1;

TStyle* ptm_style = (TStyle*)gStyle->Clone();

// Define margins and axis positions
Double_t marg_l = 0.100;
Double_t marg_r = 0.020;
Double_t marg_t = 0.020;
Double_t marg_b = 0.087;

Double_t h = 4.0;
Double_t w = 3.0;


void CreateStyle( TStyle* st ){
	st->SetName("ptm_style");

	// Set margins
	st->SetPadLeftMargin(marg_l);
	st->SetPadRightMargin(marg_r);
	st->SetPadTopMargin(marg_t);
	st->SetPadBottomMargin(marg_b);

	// Fonts
	st->SetTitleOffset( 1.0, "xyz" );
	st->SetOptStat(0);
}

void FormatFrame( TH1* h, TString xlabel, TString ylabel ){
	h->SetTitle( "" );
	h->GetXaxis()->SetTitle( xlabel );
	h->GetYaxis()->SetTitle( ylabel );
	
	h->GetXaxis()->CenterTitle();
	h->GetYaxis()->CenterTitle();
	
	h->GetXaxis()->SetTitleFont(62);
	h->GetYaxis()->SetTitleFont(62);
	h->GetZaxis()->SetTitleFont(62);
	
	h->GetXaxis()->SetLabelFont(62);
	h->GetYaxis()->SetLabelFont(62);
	h->GetZaxis()->SetLabelFont(62);
	return;	
}

void FormatFrame( TH2* h, TString xlabel, TString ylabel ){
	h->SetTitle( "" );
	h->GetXaxis()->SetTitle( xlabel );
	h->GetYaxis()->SetTitle( ylabel );
	
	h->GetXaxis()->CenterTitle();
	h->GetYaxis()->CenterTitle();
	
	h->GetXaxis()->SetTitleFont(62);
	h->GetYaxis()->SetTitleFont(62);
	h->GetZaxis()->SetTitleFont(62);
	
	h->GetXaxis()->SetLabelFont(62);
	h->GetYaxis()->SetLabelFont(62);
	h->GetZaxis()->SetLabelFont(62);
	return;	
}


void PrintCanvas( TCanvas* c, TString base_name ){
	if ( PRINT_SVG ){ c->Print( Form( "%s.svg", base_name.Data() ) ); }
	if ( PRINT_PDF ){ c->Print( Form( "%s.pdf", base_name.Data() ) ); }
	if ( PRINT_TEX ){ c->Print( Form( "%s.tex", base_name.Data() ) ); }
	if ( PRINT_PNG ){ c->Print( Form( "%s.png", base_name.Data() ) ); }
	return;
}






#endif
