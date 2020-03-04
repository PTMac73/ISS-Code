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


#ifndef ATH_EX_COMPARE_H_
#define ATH_EX_COMPARE_H_

// Switch for this is SW_EVZ
TH2F* h_evz_si;


void HCreateEVZSi(){
	h_evz_si = new TH2F( "h_evz_si", "", 450, -50, -10, 900, 0, 9 );
	h_evz_si->SetMarkerStyle(20);
	h_evz_si->SetMarkerSize(0.5);
	h_evz_si->SetMarkerColor( kRed );
	h_evz_si->GetYaxis()->SetTitle( "Energy (MeV)" );
	h_evz_si->GetXaxis()->SetTitle( "z (cm)" );
	GlobSetHistFonts( h_evz_si );
}


void HDrawEVZSi(){
	TCanvas *c_evz_si = new TCanvas( "c_evz_si", "EVZ-Si", C_WIDTH, C_HEIGHT );
	GlobSetCanvasMargins( c_evz_si );
	h_evz_si->Draw();
	
	// Print spectrum if desired
	if ( SW_EVZ[1] == 1 ){
		TString spec_name = Form( "%s/posXXX_evzsi", print_dir.Data() );
		PrintAll( c_evz_si, spec_name );
		
		if ( PRINT_ROOT == 1 ){
			out_root_file->cd();
			h_evz_si->Write();
		}
		
	}
	
	// Write SPE file if desired
	if ( SW_EVZ[2] == 1 ){
		ErrorSPE( "c_evz_si" );
	}
	
	return;
}


#endif
