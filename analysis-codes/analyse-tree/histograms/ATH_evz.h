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
TH2F* h_evz_evolution[5];


void HCreateEVZ(){
	CreateEVZSpectrum( h_evz, "h_evz" );
	for ( Int_t i = 0; i < 5; i++ ){
		CreateEVZSpectrum( h_evz_evolution[i], Form( "h_evz_evolution_%i", i ) );
	}
	
}

void HDrawEVZ(){
	TCanvas* c_evz;
	TCanvas* c_evz_evolution[5];
	
	c_evz = new TCanvas( "c_evz", "EVZ", C_WIDTH, C_HEIGHT );
	GlobSetCanvasMargins( c_evz );
	h_evz->Draw();
	
	TString root_name = Form( "%s/posXXX_evz", print_dir.Data() );
	TFile* f;
	
	// Open root file if desired
	if ( SW_EVZ[1] == 1 && PRINT_ROOT == 1 ){
		f = new TFile( ( root_name + ".root" ).Data(), "RECREATE" );
	}
	
	for ( Int_t i = 0; i < 5; i++ ){
		c_evz_evolution[i] = new TCanvas( Form( "c_evz_evolution_%i", i ), Form( "EVOLUTION OF EVZ | Case %i", i ), C_WIDTH, C_HEIGHT );
		GlobSetCanvasMargins( c_evz_evolution[i] );
		h_evz_evolution[i]->Draw();
	}
	
	// Print spectrum if desired
	if ( SW_EVZ[1] == 1 ){
		PrintAll( c_evz, Form( "%s/evz_posXXX", print_dir.Data() ) );
		for ( Int_t i = 0; i < 5; i++ ){
			PrintAll( c_evz_evolution[i], Form( "%s/evz_posXXX_evolution_%i", print_dir.Data(), i ) );
		}
	}
	
	// Write ROOT file if desired
	if ( PRINT_ROOT == 1 && SW_EVZ[1] == 1 ){ 
		f->cd(); h_evz->Write();
		for ( Int_t i = 0; i < 5; i++ ){ h_evz_evolution[i]->Write(); } 
	}
	
	// Write SPE file if desired
	if ( SW_EVZ[2] == 1 ){
		ErrorSPE( "EVZ Spectra" );
	}
	
	return;
}



#endif
