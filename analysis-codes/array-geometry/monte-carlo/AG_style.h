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

#include <TStyle.h>
#include <TColor.h>
#include <TROOT.h>

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

	// Fonts
	st->SetTitleFont( 62, "xyz" );
	st->SetTitleOffset( 1.0, "xyz" );
	st->SetLabelFont( 62, "xyz" );



	st->SetOptStat(0);
}


#endif
