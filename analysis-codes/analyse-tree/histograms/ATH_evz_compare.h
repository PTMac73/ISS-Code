// ATH_evz_compare.h
// Compare E v.s. z plots for different angle cuts
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
#include <TString.h>

#include <iostream>

#include "AT_HistogramGlobals.h"
#include "../AT_Settings.h"
#include "WriteSPE.h"


#ifndef ATH_EVZ_COMPARE_H_
#define ATH_EVZ_COMPARE_H_

// Switch for this is SW_EVZ_COMPARE
TH2F* h_evz_compare[4];


void HCreateEVZCompare(){
	for ( Int_t i = 0; i < 4; i++ ){
		h_evz_compare[i] = new TH2F( Form( "h_evz_compare_%i", i ), "", 450, -50, -10, 900, 0, 9 );
		h_evz_compare[i]->SetMarkerStyle(20);
		h_evz_compare[i]->SetMarkerSize(0.5);
		h_evz_compare[i]->SetMarkerColor( i % 2 == 0 ? kRed : kBlack );
		h_evz_compare[i]->GetYaxis()->SetTitle( "Energy (MeV)" );
		h_evz_compare[i]->GetXaxis()->SetTitle( "z (cm)" );
		GlobSetHistFonts( h_evz_compare[i] );
	}
	
	return;
}


void HDrawEVZCompare(){
	TCanvas *c_evz_compare_sing = new TCanvas( "c_evz_compare_sing", "COMPARING EVZ AT DIFFERENT ANGLES -- SINGLES", C_WIDTH, C_HEIGHT );
	GlobSetCanvasMargins( c_evz_compare_sing );
	h_evz_compare[0]->Draw();
	h_evz_compare[1]->Draw("SAME");
	
	TCanvas *c_evz_compare_clean = new TCanvas( "c_evz_compare_clean", "COMPARING EVZ AT DIFFERENT ANGLES -- CLEAN", C_WIDTH, C_HEIGHT );
	GlobSetCanvasMargins( c_evz_compare_clean );
	h_evz_compare[2]->Draw();
	h_evz_compare[3]->Draw("SAME");
	
	
	
	// Print spectra if desired
	if ( SW_EVZ_COMPARE[1] == 1 ){
		TString base = "%s/evz_comp%s_th_%3.2f-%3.2f_posXXX";
		TString root_name = Form( base.Data(), print_dir.Data(), "", THETA_LB, THETA_UB );
		TString spec_name_sing = Form( base.Data(), print_dir.Data(), "_sing", THETA_LB, THETA_UB );
		PrintAll( c_evz_compare_sing, spec_name_sing );
		
		TString spec_name_clean = Form( base.Data(), print_dir.Data(), "_clean", THETA_LB, THETA_UB );
		PrintAll( c_evz_compare_clean, spec_name_clean );
		
		if ( PRINT_ROOT == 1 ){
			TFile* f = new TFile( ( root_name + ".root" ).Data(), "RECREATE" );
			for ( Int_t i = 0; i < 4; i++ ){
				f->cd();
				h_evz_compare[i]->Write();
			}
			f->Close();
		}
	}
	
	// Write SPE file if desired
	if ( SW_EVZ_COMPARE[2] == 1 ){
		ErrorSPE( "c_evz_compare_sing and c_evz_compare_clean" );
	}

	return;
}


#endif
