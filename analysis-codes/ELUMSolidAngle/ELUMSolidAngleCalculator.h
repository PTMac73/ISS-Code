// ELUMSolidAngleCalculator.h
// Header for ELUMSolidAngleCalculator.h
// ============================================================================================= //
// Patrick MacGregor
// Nuclear Physics Research Group
// School of Physics and Astronomy
// The University of Manchester
// ============================================================================================= //
#ifndef ELUM_SOLID_ANGLE_CALCULATOR_H_
#define ELUM_SOLID_ANGLE_CALCULATOR_H_

#include <TCanvas.h>
#include <TColor.h>
#include <TH1.h>
#include <THStack.h>
#include <TMath.h>

#include "conversion.h"

// Set colours
Int_t pcb_green_i = TColor::GetFreeColorIndex();
TColor *pcb_green = new TColor( pcb_green_i, 0.0, 111.0/255.0, 75.0/255.0 );

// Define margins and axis positions
Double_t marg_l = 0.080;
Double_t marg_r = 0.010;
Double_t marg_t = 0.010;
Double_t marg_b = 0.080;

// DEFINE GLOBAL SWITCHES
Bool_t SWITCH_DRAW_CANVAS = 0;
Bool_t SWITCH_PRINT_CANVAS = 0;
Bool_t SWITCH_VERBOSE = 0;
Int_t C_HEIGHT = 900;

// Define global storage for angle things
// [REV NUMBER][TYPE][MIN/MAX]
const Int_t CHECK_NUM_REV = 2;
const Int_t CHECK_TYPES = 3;
const Int_t CHECK_MINMAX = 2;



// FUNCTIONS =================================================================================== //
// Define radius function [ mm ] - > units are [ [s^-1], [mm], [1/c], [1/c] ]
// For maximum vertical radius in y direction, y = r(1-cos(wt)), x = rsin(wt)
Double_t CalculateRadius( Double_t cyclotron_freq, Double_t z, Double_t v_para, Double_t v_perp ){
	if ( v_para == 0 ){
		return TMath::QuietNaN();
	}
	else{
		return ( v_perp*TMath::C()*MToMM()/cyclotron_freq )*( 1.0 - TMath::Cos( cyclotron_freq*z/( v_para*TMath::C()*MToMM() ) ) );
	}
}

Double_t CalculateSolidAngle( Double_t theta1, Double_t theta2 ){
	return 2*TMath::Pi()*( TMath::Cos( theta1 ) - TMath::Cos( theta2 ) );
}

struct Traj_t{
	Int_t index;	// theta_cm index
	Int_t status;	// status of trajectory - 0 = green (hits elum, miss shield), 1 = orange (hits elum, hits shield), 2 = red (misses elum, hits shield), 3 = blue (misses elum, misses shield)
	Int_t rev;		// number of turns
	Bool_t above_shield;	// Bool that measures whether the trajectory passed above or below the shield.
};


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

void FormatFrame( THStack* hs, TString xlabel, TString ylabel ){
	hs->GetXaxis()->SetTitle( xlabel );
	hs->GetYaxis()->SetTitle( ylabel );

	hs->GetXaxis()->CenterTitle();
	hs->GetYaxis()->CenterTitle();

	hs->GetXaxis()->SetTitleFont(62);
	hs->GetYaxis()->SetTitleFont(62);

	hs->GetXaxis()->SetLabelFont(62);
	hs->GetYaxis()->SetLabelFont(62);

	return;	
}


void GlobSetCanvasMargins( TCanvas *c ){
	TPad* pad = (TPad*)c;
	pad->SetLeftMargin( marg_l );
	pad->SetRightMargin( marg_r );
	pad->SetTopMargin( marg_t );
	pad->SetBottomMargin( marg_b );
	return;
}

Int_t GetCanvasWidth( Int_t h, Double_t x1, Double_t y1, Double_t x2, Double_t y2 ){
	return (Int_t)( 4*h*( ( marg_l + marg_r + x2 - x1 )/( marg_b + marg_t + y2 - y1 ) ) );
}

#endif
