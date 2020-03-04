// ATH_evz.h
// Full E v.s. z plot
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


#ifndef ATH_EVZ_H_
#define ATH_EVZ_H_

// Switch for this is SW_EVZ
TH2F* h_evz;


void HCreateEVZ(){
	h_evz = new TH2F( "h_evz", "", 450, -50, -10, 900, 0, 9 );
	h_evz->SetMarkerStyle(20);
	h_evz->SetMarkerSize(0.5);
	h_evz->SetMarkerColor( kRed );
	h_evz->GetYaxis()->SetTitle( "Energy (MeV)" );
	h_evz->GetXaxis()->SetTitle( "z (cm)" );
	GlobSetHistFonts( h_evz );
}


void HDrawEVZ(){
	TCanvas *c_evz = new TCanvas( "c_evz", "EVZ", C_WIDTH, C_HEIGHT );
	GlobSetCanvasMargins( c_evz );
	h_evz->Draw();
	
	// Print spectrum if desired
	if ( SW_EVZ[1] == 1 ){
		TString spec_name = Form( "%s/evz_posXXX", print_dir.Data() );
		PrintAll( c_evz, spec_name );
		
		if ( PRINT_ROOT == 1 ){
			out_root_file->cd();
			h_evz->Write();
		}
		
	}
	
	// Write SPE file if desired
	if ( SW_EVZ[2] == 1 ){
		ErrorSPE( "c_evz" );
	}
	
	return;
}



#endif
