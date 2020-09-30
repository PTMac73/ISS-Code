// Defines all the histogram formatting and functions etc for the AnalyseTree.C script
// ============================================================================================= //
// Patrick MacGregor
// Nuclear Physics Research Group
// Department of Physics and Astronomy
// The University of Manchester
// ============================================================================================= //

#include <TCanvas.h>
#include <TCutG.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TLine.h>
#include <TObjArray.h>
#include <TString.h>

#include <iostream>

#include "AT_HistogramGlobals.h"
#include "../AT_Settings.h"
#include "WriteSPE.h"


#ifndef ATH_EX_SI_H_
#define ATH_EX_SI_H_

// Switch for this is SW_EX_SI
// 0 is the Mg spectrum
// 1 is the Si spectrum
// 2 is the Mg spectrum in terms of Si excitation
// 3 is the Si spectrum in terms of Si excitation
TH1F* h_ex_si[4];		// In terms of Mg excitation
TFile *f;
const Int_t NUM_SI_STATES= 10;
TLine *si_state_line[NUM_SI_STATES];

// Define Si states and colours
Double_t si_state_ensdf[NUM_SI_STATES] = {
	3.06713,
	3.62349,
	6.69593,
	7.9968,
	9.252,
	9.298,
	9.326,
	9.392,
	9.413,
	9.518
};

Int_t si_state_ensdf_colour[NUM_SI_STATES] = {
	kRed,
	kOrange - 9,
	kBlack,
	kRed,
	kMagenta,
	kMagenta,
	kMagenta,
	kMagenta,
	kMagenta,
	kMagenta
};




void HCreateExSi(){
	for ( Int_t i = 0; i < 2; i++ ){
		CreateExSpectrum( h_ex_si[i], Form( "h_ex_si_%i", i ) );
		CreateExSpectrum( h_ex_si[i+2], Form( "h_ex_si_%i", i+2 ), 5, 15 );
	}
	for ( Int_t i = 0; i < 2; i++ ){
		h_ex_si[2*i]->SetFillColor( kGray );
		h_ex_si[2*i]->SetLineColor( kGray );
		h_ex_si[2*i+1]->SetFillColor( kRed );
		h_ex_si[2*i+1]->SetLineColor( kRed );
	}
	
	return;
}

void HDrawExSi(){
	TCanvas *c_ex_si[2];
	
	if ( SW_EX_SI[1] == 1 && PRINT_ROOT == 1 ){
		f = new TFile( Form( "%s/pos%i_exsi.root", print_dir.Data(), ARR_POSITION ), "RECREATE" );
	}
	
	// Create Si state lines
	for ( Int_t i = 0; i < NUM_SI_STATES; i++ ){
		si_state_line[i] = new TLine( si_state_ensdf[i], 0, si_state_ensdf[i], 1.05*h_ex_si[2]->GetMaximum() );
		si_state_line[i]->SetLineColor( si_state_ensdf_colour[i] );
	}
	
	
	
	
	// In terms of Mg
	c_ex_si[0] = new TCanvas( "c_ex_si_0", "EX-Si-Mg", C_WIDTH, C_HEIGHT );
	GlobSetCanvasMargins( c_ex_si[0] );
	for ( Int_t i = 0; i < 2; i++ ){
		h_ex_si[i]->Draw("SAME");
	}
	
	// In terms of Si
	c_ex_si[1] = new TCanvas( "c_ex_si_1", "EX-Si", C_WIDTH, C_HEIGHT );
	GlobSetCanvasMargins( c_ex_si[1] );
	SetPadTicks( c_ex_si[1] );
	h_ex_si[2]->Draw();
	h_ex_si[3]->Draw("SAME");
	
	
	/*TObjArray* si_lines = new TObjArray();
	for ( Int_t i = 0; i < NUM_SI_STATES; i++ ){
		si_state_line[i]->Draw("SAME");
		si_lines->Add( si_state_line[i] );
	}*/
	
	TString spec_name = Form( "%s/pos%i_exsi", print_dir.Data(), ARR_POSITION );

	// Print spectrum if desired
	if ( SW_EX_SI[1] == 1 ){
		PrintAll( c_ex_si[0], spec_name + "_mg" );
		PrintAll( c_ex_si[1], spec_name );
		
		if ( PRINT_ROOT == 1 ){
			f->cd();
			for ( Int_t i = 0; i < 4; i++ ){
				h_ex_si[i]->Write();
			}
			//si_lines->Write("si_lines", TObject::kSingleKey );
		}
		
	}
	
	// Write SPE file if desired
	if ( SW_EX_SI[2] == 1 ){
		WriteSPE( h_ex_si[2]->GetName(), Form( "%s", ( spec_name + "_mg" ).Data() ) );	// Mg in terms of Si
		WriteSPE( h_ex_si[3]->GetName(), Form( "%s", ( spec_name + "_si" ).Data() ) );	// Si in terms of Si
	}

	// Close ROOT file
	if ( SW_EX_SI[1] == 1 && PRINT_ROOT == 1 && f != NULL ){ if ( f->IsOpen() ){ f->Close(); } }

	return;
}


#endif



























