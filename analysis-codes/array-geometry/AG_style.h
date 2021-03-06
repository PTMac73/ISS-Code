// AG_style.h
// TStyle definitions for the file ArrayGeometry.C
// ============================================================================================= //
// Patrick MacGregor
// Nuclear Physics Research Group
// School of Physics and Astronomy
// The University of Manchester
// ============================================================================================= //
#ifndef AG_STYLE_H_
#define AG_STYLE_H_


#include <TColor.h>
#include <TH1F.h>
#include <TROOT.h>
#include <TString.h>
#include <TStyle.h>

TStyle* ptm_style = (TStyle*)gStyle->Clone();

// Define margins and axis positions
Double_t marg_l = 0.100;
Double_t marg_r = 0.020;
Double_t marg_t = 0.020;
Double_t marg_b = 0.087;

Double_t h = 4.0;
Double_t w = 3.0;

Double_t X1 = -50.0;
Double_t X2 = 5.0;
Double_t DY = (X2-X1)*(1 - marg_t - marg_b)*h/( (1 - marg_r - marg_l)*w );
Double_t Y1 = -0.5*DY;
Double_t Y2 = 0.5*DY;


void CreateStyle( TStyle* st ){
	st->SetName("ptm_style");

	// Set margins
	st->SetPadLeftMargin(marg_l);
	st->SetPadRightMargin(marg_r);
	st->SetPadTopMargin(marg_t);
	st->SetPadBottomMargin(marg_b);
}



void FormatFrame( TH1F* h, TString xlabel, TString ylabel ){
	h->SetTitle( "" );
	h->GetXaxis()->SetTitle( xlabel );
	h->GetYaxis()->SetTitle( ylabel );
	
	h->GetXaxis()->CenterTitle();
	h->GetYaxis()->CenterTitle();
	
	h->GetXaxis()->SetTitleFont(62);
	h->GetYaxis()->SetTitleFont(62);
	
	h->GetXaxis()->SetLabelFont(62);
	h->GetYaxis()->SetLabelFont(62);
	return;	
}



#endif
