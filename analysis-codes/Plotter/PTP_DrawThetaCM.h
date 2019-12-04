// PTP_DrawThetaCM.h
// Draws the theta_cm histogram for every row/detector-by-detector
// ============================================================================================= //
// Patrick MacGregor
// Nuclear Physics Research Group
// Department of Physics and Astronomy
// The University of Manchester
// ============================================================================================= //
#ifndef PTP_DRAW_THETA_CM_H_
#define PTP_DRAW_THETA_CM_H_

#include "PTPlotterINIT.h"
#include "PTF_GetPosNumber.h"

#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TPad.h>

#include <iostream>

struct plotterOptions;

const Int_t NUM_STATES = 9;

Double_t ex_lims[NUM_STATES][2] = {
	{ -0.3, 0.30 },
	{ 0.80, 1.25 },
	{ 1.26, 1.65 },
	{ 2.00, 2.37 },
	{ 2.38, 2.70 },
	{ 2.77, 3.02 },
	{ 3.03, 3.40 },
	{ 3.60, 4.13 },
	{ 4.14, 4.66 }
};

// ThetaCM limits (based on eyeballing, with array size 2 for position)
Double_t thetaCM_lims[NUM_STATES][2] = {
	{ 15.0, 19.0 },
	{ 12.5, 12.5 },
	{ 15.5, 16.0 },
	{ 15.0, 15.0 },
	{ 13.5, 13.0 },
	{ 0.00, 0.00 },
	{ 13.0, 13.5 },
	{ 13.0, 13.0 },
	{ 13.0, 10.5 }
};



// FUNCTIONS ----------------------------------------------------------------------------------- //
TString GenerateXCALCutString( Double_t EX_LB = -0.5, Double_t EX_UB = 0.5 ){
	return Form( "( cut0 || cut1 || cut2 || cut3 ) && td_rdt_e[] >= td_rdt_e_cuts[][0] && td_rdt_e[] <= td_rdt_e_cuts[][1] && xcal[] >= xcal_cuts[][0] && xcal[] <= xcal_cuts[][1] && ( xcal[] <= xcal_cuts[][2] || xcal[] >= xcal_cuts[][3] ) && Ex >= %f && Ex <= %f", EX_LB, EX_UB );
}

TString GenerateXCALCutStringExclusion( Double_t LB, Double_t UB, Double_t EX_LB = -0.5, Double_t EX_UB = 0.5 ){
	return GenerateXCALCutString( EX_LB, EX_UB ) + (TString)Form( "&& thetaCM >= %f && thetaCM <= %f", LB, UB );
} 


// MAIN ======================================================================================== //
void DrawThetaCM( TTree *t, plotterOptions &opt_s, Int_t pos_number ){
	gStyle->SetOptStat(kFALSE);
	printDiv(); printf( "PLOTTING THE THETA_CM GRAPHS WITH ANGLE\n" ); printDiv();

	TCanvas* c_theta[NUM_STATES];
	TH1F* h_theta[NUM_STATES][2];
	TPad* pad;
	TFile* f = new TFile( Form( "%s/thetaCM-plots/thetaCMex_pos%i_states.root", opt_s.printDir.Data(), pos_number ), "RECREATE" );

	// Draw the graphs
	for ( Int_t i = 0; i < NUM_STATES; i++ ){
		c_theta[i] = new TCanvas( Form( "c_theta%i", i ), Form( "CANVAS %i", i ), C_WIDTH, C_HEIGHT );		
		t->Draw( Form( "thetaCM>>h_theta%i( 100, 0, 50 )", i ), GenerateXCALCutStringExclusion( thetaCM_lims[i][pos_number-1], 50.0, ex_lims[i][0], ex_lims[i][1] ), "goff" );
		h_theta[i][0] = (TH1F*)gDirectory->Get( Form( "h_theta%i", i ) );
		h_theta[i][0]->SetTitle("");
		h_theta[i][0]->GetXaxis()->SetTitle("#theta_{cm}");
		h_theta[i][0]->GetYaxis()->SetTitle("Counts");
		h_theta[i][0]->GetYaxis()->SetTitleOffset(1.4);
		h_theta[i][0]->GetXaxis()->SetTitleOffset(1.0);
		h_theta[i][0]->SetFillColor(kGreen);
		h_theta[i][0]->Draw();

		t->Draw( Form( "thetaCM>>h_theta_ex%i( 100, 0, 50 )", i ), GenerateXCALCutStringExclusion( 0.0, thetaCM_lims[i][pos_number-1], ex_lims[i][0], ex_lims[i][1] ), "goff" );
		h_theta[i][1] = (TH1F*)gDirectory->Get( Form( "h_theta_ex%i", i ) );
		h_theta[i][1]->SetFillColor(kRed);
		h_theta[i][1]->Draw("SAME");
		

		// Print the graphs
		if ( SWITCH_PRINT_CANVAS == 1 ){
			c_theta[i]->Print( Form( "%s/thetaCM-plots/thetaCMex_pos%i_state%i%s", opt_s.printDir.Data(), pos_number, i+1, PRINT_FORMAT.Data() ) );
		}

		// Write to ROOT file
		h_theta[i][0]->Write();
		h_theta[i][1]->Write();

		std::cout << "State " << i << " plotted\n";

	}


	// Close ROOT file
	f->Close();



}















#endif
