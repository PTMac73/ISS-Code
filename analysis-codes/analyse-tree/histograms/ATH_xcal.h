// ATH_xcal.h
// Xcal spectrum on a detector-by-detector basis
// ============================================================================================= //
// Patrick MacGregor
// Nuclear Physics Research Group
// Department of Physics and Astronomy
// The University of Manchester
// ============================================================================================= //
#ifndef ATH_XCAL_H_
#define ATH_XCAL_H_

#include <TCanvas.h>
#include <TCutG.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <THStack.h>
#include <TLine.h>
#include <TMath.h>
#include <TString.h>
#include <TStyle.h>


#include <iostream>

#include "AT_HistogramGlobals.h"
#include "../AT_Settings.h"
#include "WriteSPE.h"

// Switch for this is SW_XCAL
TH1F* h_xcal[24][2];		// Full spectrum
TH2F* h_xcal_e[24][2];		// E v.s. xcal plot
TH1F* h_xcal_full_comp[3];	// Divide xcal up by colour
THStack* hs_xcal_full_comp;

Bool_t xcal_print_opt[3] = {
	1,	// (0) xcal
	1,	// (1) xcal-energy
	0,	// (2) full xcal spectrum (colour-coded)
};


void HCreateXCAL(){
	// *LOOP* over detectors
	for ( Int_t i = 0; i < 24; i++ ){
		// Detector by detector
		if ( ( i == DET_NUMBER || ( DET_NUMBER == -1 ) ) && det_array[i % 6][(Int_t)TMath::Floor(i/6)] != 0 ){
			for ( Int_t j = 0; j < 2; j++ ){
				h_xcal[i][j] = new TH1F( Form( "h_xcal_%s%i", ( j == 0 ? "" : "cut_" ), i ), "", 200, -0.5, 1.5 );
				h_xcal[i][j]->SetTitle("");
				h_xcal[i][j]->GetXaxis()->SetTitle( ( j == 0 ? "X_{psd}" : "" ) );
				h_xcal[i][j]->GetYaxis()->SetTitle( ( j == 0 ? "Counts" : "" ) );	
				h_xcal[i][j]->SetFillColor( ( j == 0 ? kRed : kYellow ) );	
				h_xcal[i][j]->SetLineWidth(1);
				h_xcal[i][j]->SetLineColor(kBlack);
				GlobSetHistFonts( h_xcal[i][j] );
				
				h_xcal_e[i][j] = new TH2F( Form( "h_xcal_e_%s%i", ( j == 0 ? "" : "cut_" ), i ), "", 200, -0.5, 1.5, 900, 0, 9 );
				GlobCreate2DHists( h_xcal_e[i][j], "xcal", "Energy (MeV)" );
			}
			h_xcal_e[i][1]->SetMarkerSize(0.3);
			h_xcal_e[i][1]->SetMarkerColor(kYellow);
		
		}
	}
	
	for ( Int_t i = 0; i < 3; i++ ){
		h_xcal_full_comp[i] = new TH1F( Form( "h_xcal_full_comp_%i", i ), "", 200, -0.5, 1.5 );
	}
	h_xcal_full_comp[0]->SetFillColor(kGreen-7);
	h_xcal_full_comp[1]->SetFillColor(kRed-7);
	h_xcal_full_comp[2]->SetFillColor(kBlue-7);
	
	h_xcal_full_comp[0]->SetLineColor(kGreen-7);
	h_xcal_full_comp[1]->SetLineColor(kRed-7);
	h_xcal_full_comp[2]->SetLineColor(kBlue-7);
	
	return;
}


void HDrawXCAL(){
	// Define some local variables
	TCanvas* c_xcal_comb;
	TCanvas* c_xcal[24];
	TCanvas* c_xcal_e[24];
	TCanvas* c_xcal_full_comp;
	
	if ( CANVAS_COMBINE == 1 ){
		c_xcal_comb = new TCanvas( "c_xcal_comb", "XCAL Spectrum Combined", 2*C_WIDTH, 2*C_HEIGHT );
		c_xcal_comb->Divide(6,4);
	}
	
	TString root_name = Form( "%s/posXXX_xcal", print_dir.Data() );
	TFile* f;
	TString spec_name;
	
	// Open root file if desired
	if ( SW_XCAL[1] == 1 && PRINT_ROOT == 1 ){
		f = new TFile( ( root_name + ".root" ).Data(), "RECREATE" );
	}
	
	
	// *LOOP* over detectors (i)
	for ( Int_t i = 0; i < 24; i++ ){

		// Detector by detector
		if ( ( i == DET_NUMBER || DET_NUMBER == -1 ) && det_array[i % 6][(Int_t)TMath::Floor(i/6)] != 0 ){
			spec_name = Form( "%s/posXXX_xcal_%i", print_dir.Data(), i );
			
			TLine* cut_line[2];
			for ( Int_t j = 0; j < 2; j++ ){
				cut_line[j] = new TLine(
					XCAL_cuts[i][j],				// X1
					0.0,							// Y1
					XCAL_cuts[i][j],				// X2
					1.05*h_xcal[i][0]->GetMaximum()	// Y2
				);
				cut_line[j]->SetLineWidth(1);
				cut_line[j]->SetLineStyle(2);
				cut_line[j]->SetLineColor(kBlack);
			}
			
			
			// Plot spectrum
			c_xcal[i] = new TCanvas( Form( "c_xcal_%i", i ), Form( "xcal | Det %i", i ), C_WIDTH, C_HEIGHT );
			GlobSetCanvasMargins( c_xcal[i] );
			TPad* pad = (TPad*)c_xcal[i]->GetPad(0);
			pad->SetTicks(1,1);
			h_xcal[i][0]->Draw();
			h_xcal[i][1]->Draw("SAME");
			cut_line[0]->Draw("SAME");
			cut_line[1]->Draw("SAME");
			gStyle->SetTitleFont(62);
			
			c_xcal_e[i] = new TCanvas( Form( "c_xcal_e_%i", i ), Form( "xcal-E | Det %i", i ), C_WIDTH, C_HEIGHT );
			GlobSetCanvasMargins( c_xcal_e[i] );
			h_xcal_e[i][0]->Draw();
			h_xcal_e[i][1]->Draw("SAME");
			
			// Print spectrum if desired
			if ( SW_XCAL[1] == 1 ){
				if ( xcal_print_opt[0] == 1 ){ PrintAll( c_xcal[i], spec_name ); }
				if ( xcal_print_opt[1] == 1 ){ PrintAll( c_xcal_e[i], Form( "%s/posXXX_xcal_e_%i", print_dir.Data(), i ) ); }
			}
				
			// Write ROOT file if desired
			if ( PRINT_ROOT == 1 && SW_XCAL[1] == 1 ){
				f->cd();
				if ( xcal_print_opt[0] == 1 ){ h_xcal[i][0]->Write();h_xcal[i][1]->Write(); }
				if ( xcal_print_opt[1] == 1 ){ h_xcal_e[i][0]->Write();h_xcal_e[i][1]->Write(); }
			}
			
			// Write SPE file if desired
			if ( SW_XCAL[2] == 1 ){ ErrorSPE( "XCAL spectra" ); }
		}
		
	} // *LOOP* (i)
		
	// Draw the xcal stack
	std::cout << h_xcal_full_comp[0]->GetSize() << "\n";
	hs_xcal_full_comp = new THStack( "hs", "" );
	for ( Int_t i = 0; i < 3; i++ ){
		hs_xcal_full_comp->Add( h_xcal_full_comp[i] );
	}
	
	c_xcal_full_comp = new TCanvas( "c_xcal_full_comp", "", C_WIDTH, C_HEIGHT );
	GlobSetCanvasMargins( c_xcal_full_comp );
	hs_xcal_full_comp->Draw();
	
	TPad* pad = (TPad*)c_xcal_full_comp->GetPad(0);
	pad->SetTicks(1,1);
	
	GlobSetHistFonts(hs_xcal_full_comp);
	hs_xcal_full_comp->GetXaxis()->SetTitle("X_{cal}");
	hs_xcal_full_comp->GetYaxis()->SetTitle("Counts");
	
	c_xcal_full_comp->Modified(); c_xcal_full_comp->Update();
	
	
	if ( SW_XCAL[1] == 1 && xcal_print_opt[2] == 1  ){
		PrintAll( c_xcal_full_comp, Form( "%s/posXXX_xcal_full_comp", print_dir.Data() ) );
		if ( PRINT_ROOT == 1 ){
			f->cd(); hs_xcal_full_comp->Write();
		}
	}
	
	
	
	// Close the root file
	if ( SW_XCAL[1] == 1 && PRINT_ROOT == 1 && f != NULL ){ if ( f->IsOpen() ){ f->Close(); } }
	return;
}

#endif
