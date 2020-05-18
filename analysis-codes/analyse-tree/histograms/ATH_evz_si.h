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
#include <TString.h>

#include <iostream>

#include "AT_HistogramGlobals.h"
#include "../AT_Settings.h"
#include "WriteSPE.h"


#ifndef ATH_EVZ_SI_H_
#define ATH_EVZ_SI_H_

// Switch for this is SW_EVZ_SI
// 0 is the singles spectrum
// 1 is the Mg spectrum
// 2 is the Si spectrum
TH2F* h_evz_si[3];


void HCreateEVZSi(){
	for ( Int_t i = 0; i < 3; i++ ){
		h_evz_si[i] = new TH2F( Form( "h_evz_si_%i", i ), "", 400, -50, -10, 900, 0, 9 );
		h_evz_si[i]->SetMarkerStyle(20);
		h_evz_si[i]->SetMarkerSize(0.5);
		GlobSetHistFonts( h_evz_si[i] );
		
		if ( i == 0 ){
			h_evz_si[i]->SetMarkerColorAlpha( kBlack, 0.1 );
			h_evz_si[i]->GetYaxis()->SetTitle( "Energy (MeV)" );
			h_evz_si[i]->GetXaxis()->SetTitle( "z (cm)" );
		}
		else if ( i == 1 ){
			h_evz_si[i]->SetMarkerColorAlpha( kRed, 0.5 );
		}
		else if ( i == 2 ) {
			h_evz_si[i]->SetMarkerColor( kBlue );
		}
		
	}
	
}


void HDrawEVZSi(){
	TCanvas *c_evz_si = new TCanvas( "c_evz_si", "EVZ-Si", C_WIDTH, C_HEIGHT );
	GlobSetCanvasMargins( c_evz_si );
	for ( Int_t i = 0; i < 3; i++ ){
		h_evz_si[i]->Draw("SAME");
	}
	
	// Print spectrum if desired
	if ( SW_EVZ_SI[1] == 1 ){
		TString spec_name = Form( "%s/posXXX_evzsi", print_dir.Data() );
		PrintAll( c_evz_si, spec_name );
		
		/*if ( PRINT_ROOT == 1 ){
			out_root_file->cd();
			for ( Int_t i = 0; i < 3; i++ ){
				h_evz_si[i]->Write();
			}
		}*/
		
	}
	
	// Write SPE file if desired
	if ( SW_EVZ_SI[2] == 1 ){
		ErrorSPE( "c_evz_si" );
	}
	
	return;
}


#endif
