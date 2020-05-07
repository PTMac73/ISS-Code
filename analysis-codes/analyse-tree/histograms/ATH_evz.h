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
TH2F* h_evz_bands[5];


void HCreateEVZ(){
	CreateEVZSpectrum( h_evz, "h_evz" );
	for ( Int_t i = 0; i < 5; i++ ){
		CreateEVZSpectrum( h_evz_evolution[i], Form( "h_evz_evolution_%i", i ) );
		CreateEVZSpectrum( h_evz_bands[i], Form( "h_evz_bands_%i", i ) );
	}
	h_evz_bands[0]->SetMarkerColor( kBlack );
	h_evz_bands[1]->SetMarkerColor( kRed );
	h_evz_bands[2]->SetMarkerColor( kGreen );
	h_evz_bands[3]->SetMarkerColor( kBlue );
	h_evz_bands[4]->SetMarkerColor( kCyan );
	
}

void HDrawEVZ(){
	TCanvas* c_evz;
	TCanvas* c_evz_evolution[5];
	TCanvas* c_evz_bands;
	
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
	
	c_evz_bands = new TCanvas( "c_evz_bands", "EVZ 2\u00b0 BANDS", C_WIDTH, C_HEIGHT );
	GlobSetCanvasMargins( c_evz_bands );
	h_evz_bands[0]->Draw();
	h_evz_bands[1]->Draw("SAME");
	h_evz_bands[2]->Draw("SAME");
	h_evz_bands[3]->Draw("SAME");
	h_evz_bands[4]->Draw("SAME");
	
	// Print spectrum if desired
	if ( SW_EVZ[1] == 1 ){
		PrintAll( c_evz, Form( "%s/posXXX_evz", print_dir.Data() ) );
		for ( Int_t i = 0; i < 5; i++ ){
			PrintAll( c_evz_evolution[i], Form( "%s/posXXX_evz_evolution_%i", print_dir.Data(), i ) );
		}
		PrintAll( c_evz_bands, Form( "%s/posXXX_evz_bands", print_dir.Data() ) );
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
