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
TH2F* h_evz_custom;
TH2F* h_evz_evolution[5];
TH2F* h_evz_bands[5];
TH2F* h_evz_sides[4];

Bool_t evz_print_opt[5] = {
	0, // (0) EVZ Spectrum  --> Standard spectrum
	0, // (1) EVZ Custom    --> Custom thetaCM cuts to avoid the lack of recoil-coincidences
	1, // (2) EVZ Evolution --> Evolves through the cuts
	0, // (3) EVZ Bands     --> Splits the spectrum up into 2 DEG bands
	1  // (4) EVZ Sides     --> Does each side of the array individually
};


void HCreateEVZ(){
	// (0) EVZ Spectrum
	CreateEVZSpectrum( h_evz, "h_evz" );
	
	// (1) EVZ Custom
	CreateEVZSpectrum( h_evz_custom, "h_evz_custom" );
	
	// (2) EVZ Evolution
	for ( Int_t i = 0; i < 5; i++ ){
		CreateEVZSpectrum( h_evz_evolution[i], Form( "h_evz_evolution_%i", i ) );
	}
	
	// (3) EVZ Bands
	for ( Int_t i = 0; i < 5; i++ ){
		CreateEVZSpectrum( h_evz_bands[i], Form( "h_evz_bands_%i", i ) ); 
	}
	h_evz_bands[0]->SetMarkerColor( kBlack );
	h_evz_bands[1]->SetMarkerColor( kRed );
	h_evz_bands[2]->SetMarkerColor( kGreen );
	h_evz_bands[3]->SetMarkerColor( kBlue );
	h_evz_bands[4]->SetMarkerColor( kCyan );
	
	// (4) EVZ Sides
	for ( Int_t i = 0; i < 4; i++ ){
		CreateEVZSpectrum( h_evz_sides[i], Form( "h_evz_%s", SideString( 6*i ).Data() ) );
	}
}

void HDrawEVZ(){
	TCanvas* c_evz;
	TCanvas* c_evz_evolution[5];
	TCanvas* c_evz_bands;
	TCanvas* c_evz_custom;
	TCanvas* c_evz_sides[4];
	
	TString root_name = Form( "%s/pos%i_evz", print_dir.Data(), ARR_POSITION );
	TFile* f;
	
	// Open root file if desired
	if ( SW_EVZ[1] == 1 && PRINT_ROOT == 1 ){
		f = new TFile( ( root_name + ".root" ).Data(), "RECREATE" );
	}
	
	
	// (0) EVZ Spectrum
	c_evz = new TCanvas( "c_evz", "EVZ", C_WIDTH, C_HEIGHT );
	GlobSetCanvasMargins( c_evz );
	h_evz->Draw();
	
	// (1) EVZ Custom
	c_evz_custom = new TCanvas( "c_evz_custom", "EVZ custom", C_WIDTH, C_HEIGHT );
	GlobSetCanvasMargins( c_evz_custom );
	h_evz_custom->Draw();
	
	// (2) EVZ Evolution
	for ( Int_t i = 0; i < 5; i++ ){
		c_evz_evolution[i] = new TCanvas( Form( "c_evz_evolution_%i", i ), Form( "EVOLUTION OF EVZ | Case %i", i ), C_WIDTH, C_HEIGHT );
		GlobSetCanvasMargins( c_evz_evolution[i] );
		h_evz_evolution[i]->Draw();
	}
	
	// (3) EVZ Bands
	c_evz_bands = new TCanvas( "c_evz_bands", "EVZ 2\u00b0 BANDS", C_WIDTH, C_HEIGHT );
	GlobSetCanvasMargins( c_evz_bands );
	h_evz_bands[0]->Draw();
	h_evz_bands[1]->Draw("SAME");
	h_evz_bands[2]->Draw("SAME");
	h_evz_bands[3]->Draw("SAME");
	h_evz_bands[4]->Draw("SAME");

	
	// (4) EVZ Sides
	for ( Int_t i = 0; i < 4; i++ ){
		c_evz_sides[i] = new TCanvas( Form( "c_evz_%s", SideString(6*i).Data() ), Form( "EVZ FOR %s SIDE OF ARRAY", SideString(6*i).Data() ) );
		GlobSetCanvasMargins( c_evz_sides[i] );
		h_evz_sides[i]->Draw();
	}
	
	
	// Print spectrum if desired
	if ( SW_EVZ[1] == 1 ){
		if ( evz_print_opt[0] == 1 ){ PrintAll( c_evz, Form( "%s/pos%i_evz", print_dir.Data(), ARR_POSITION ) ); }		// (0) EVZ Spectrum
		if ( evz_print_opt[1] == 1 ){ PrintAll( c_evz_custom, Form( "%s/pos%i_evz_custom", print_dir.Data(), ARR_POSITION ) ); }							// (1) EVZ Custom
		if ( evz_print_opt[2] == 1 ){ 
			for ( Int_t i = 0; i < 5; i++ ){ PrintAll( c_evz_evolution[i], Form( "%s/pos%i_evz_evolution_%i", print_dir.Data(), ARR_POSITION, i ) ); }	// (2) EVZ Evolution
		}
		if ( evz_print_opt[3] == 1 ){ PrintAll( c_evz_bands, Form( "%s/pos%i_evz_bands", print_dir.Data(), ARR_POSITION ) ); }	// (3) EVZ Bands
		if ( evz_print_opt[4] == 1 ){
			for ( Int_t i = 0; i < 4; i++ ){ PrintAll( c_evz_sides[i], Form( "%s/pos%i_evz_sides_%s", print_dir.Data(), ARR_POSITION, SideString(6*i).Data() ) ); }	// (4) EVZ Sides
		}
	}
	
	// Write ROOT file if desired
	if ( PRINT_ROOT == 1 && SW_EVZ[1] == 1 ){ 
		f->cd(); 
		if ( evz_print_opt[0] == 1 ){ h_evz->Write(); }						// (0) EVZ Spectrum
		if ( evz_print_opt[1] == 1 ){ h_evz_custom->Write(); } 				// (1) EVZ Custom
		if ( evz_print_opt[2] == 1 ){ 
			for ( Int_t i = 0; i < 5; i++ ){ h_evz_evolution[i]->Write(); }	// (2) EVZ Evolution
		}
		if ( evz_print_opt[3] == 1 ){ 
			for ( Int_t i = 0; i < 5; i++ ){ h_evz_bands[i]->Write(); }		// (3) EVZ Bands
		}
		if ( evz_print_opt[4] == 1 ){ 
			for ( Int_t i = 0; i < 4; i++ ){ h_evz_sides[i]->Write(); }		// (4) EVZ Sides
		}
	}
	
	// Write SPE file if desired
	if ( SW_EVZ[2] == 1 ){
		ErrorSPE( "EVZ Spectra" );
	}
	
	return;
}



#endif
