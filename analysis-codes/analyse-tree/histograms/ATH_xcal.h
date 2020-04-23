// ATH_xcal.h
// Xcal spectrum on a detector-by-detector basis
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
#include <TMath.h>
#include <TString.h>
#include <TStyle.h>


#include <iostream>

#include "AT_HistogramGlobals.h"
#include "../AT_Settings.h"
#include "WriteSPE.h"


#ifndef ATH_XCAL_H_
#define ATH_XCAL_H_

// Switch for this is SW_XCAL
TH1F* h_xcal[24][2];		// Full spectrum


void HCreateXCAL(){
	// *LOOP* over detectors
	for ( Int_t i = 0; i < 24; i++ ){
		// Detector by detector
		if ( ( i == DET_NUMBER || ( DET_NUMBER == -1 ) ) && det_array[i % 6][(Int_t)TMath::Floor(i/6)] != 0 ){
			h_xcal[i][0] = new TH1F( Form( "h_xcal_%i", i ), "", 200, -0.5, 1.5 );
			h_xcal[i][0]->SetTitle("");
			h_xcal[i][0]->GetXaxis()->SetTitle("xcal");
			h_xcal[i][0]->GetYaxis()->SetTitle("Counts");	
			h_xcal[i][0]->SetFillColor(kRed);	
			h_xcal[i][0]->SetLineWidth(1);
			h_xcal[i][0]->SetLineColor(kBlack);
			GlobSetHistFonts( h_xcal[i][0] );
			
			h_xcal[i][1] = new TH1F( Form( "h_xcal_cut_%i", i ), "", 200, -0.5, 1.5 );
			h_xcal[i][1]->SetTitle("");
			h_xcal[i][1]->GetXaxis()->SetTitle("");
			h_xcal[i][1]->GetYaxis()->SetTitle("");	
			h_xcal[i][1]->SetFillColor(kYellow);
			h_xcal[i][1]->SetLineWidth(1);
			h_xcal[i][1]->SetLineColor(kBlack);
			GlobSetHistFonts( h_xcal[i][1] );
		}
	}
	
	return;
}


void HDrawXCAL(){
	// Define some local variables
	TCanvas* c_xcal_comb;
	TCanvas* c_xcal[24];
	
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
			
			if ( CANVAS_COMBINE == 0 ){
				// Plot spectrum
				c_xcal[i] = new TCanvas( Form( "c_xcal_%i", i ), Form( "xcal | Det %i", i ), C_WIDTH, C_HEIGHT );
				GlobSetCanvasMargins( c_xcal[i] );
				h_xcal[i][0]->Draw();
				h_xcal[i][1]->Draw("SAME");
				gStyle->SetTitleFont(62);
				
				// Print spectrum if desired
				if ( SW_XCAL[1] == 1 ){
					PrintAll( c_xcal[i], spec_name );
				}
				
			}
			else{
				c_xcal_comb->cd(i+1);
				h_xcal[i][0]->SetTitle( Form( "Det #%i", i ) );
				h_xcal[i][0]->Draw();
				h_xcal[i][1]->Draw("SAME");
				//TPad* pad = (TPad*)c_xcal_comb->GetPad(i+1);
				//SetCanvasTitleFont( pad );
				
				// Print spectrum if desired
				if ( SW_XCAL[1] == 1 && i == 23 ){
					PrintAll( c_xcal_comb, Form( "%s/posXXX_xcal_comb", print_dir.Data() ) );
				}
				
				
			}
			
			// Write ROOT file if desired
			if ( PRINT_ROOT == 1 && SW_XCAL[1] == 1 ){ f->cd(); h_xcal[i][0]->Write(); h_xcal[i][1]->Write(); }
			
			// Write SPE file if desired
			if ( SW_XCAL[2] == 1 ){ ErrorSPE( "XCAL spectra" ); }
		} 
	
	
	}
	// Close the root file
	if ( SW_XCAL[1] == 1 && PRINT_ROOT == 1 && f != NULL ){ if ( f->IsOpen() ){ f->Close(); } }
	return;
}

#endif
