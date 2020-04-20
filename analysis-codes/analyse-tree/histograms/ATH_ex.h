// ATH_ex.h
// Excitation spectrum on a row-by-row/full basis
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
#include <TStyle.h>


#include <iostream>

#include "AT_HistogramGlobals.h"
#include "../AT_Settings.h"
#include "WriteSPE.h"


#ifndef ATH_EX_H_
#define ATH_EX_H_

// Switch for this is SW_EX
// Number of histograms = 24 (DBD) + 6 (RBR) + 1 (Full)
TH1F* h_ex_full = NULL;		// Full spectrum
TH1F* h_ex_rbr[6];		// Full spectrum
TH1F* h_ex_dbd[24];		// Full spectrum


void HCreateEx(){
	// *LOOP* over detectors
	for ( Int_t i = 0; i < 24; i++ ){
	
		// Full spectrum
		if ( i == 0 && ALL_ROWS == 1 ){
			CreateExSpectrum( h_ex_full, "h_ex_full" );
		}
		
		// Row by row
		if ( ( i == ROW_NUMBER || ( ROW_NUMBER == -1 && i < 6 ) ) && ROW_BY_ROW == 1 ){
			CreateExSpectrum( h_ex_rbr[i], Form( "h_ex_rbr_%i", i ) );
		}
		
		// Detector by detector
		if ( ( i == DET_NUMBER || ( DET_NUMBER == -1 ) ) && DET_BY_DET == 1 ){
			CreateExSpectrum( h_ex_dbd[i], Form( "h_ex_dbd_%i", i ) );
		}
	}
	
	return;
}


void HDrawEx(){
	// Define some local variables
	TCanvas* c_ex_full;
	TCanvas* c_ex_rbr[6];
	TCanvas* c_ex_dbd[24];
	TCanvas* c_ex_dbd_comb;
	if ( CANVAS_COMBINE == 1 ){
		c_ex_dbd_comb = new TCanvas( "c_ex_dbd_comb", "DBD Ex Spectrum", C_WIDTH, C_HEIGHT );
		c_ex_dbd_comb->Divide(6,4);
	}
	
	TString root_name = Form( "%s/posXXX_ex", print_dir.Data() );
	TFile* f;
	TString spec_name;
	
	// Open root file if desired
	if ( SW_EX[1] == 1 && PRINT_ROOT == 1 ){
		f = new TFile( ( root_name + ".root" ).Data(), "RECREATE" );
	}
	
	
	// *LOOP* over detectors (i)
	for ( Int_t i = 0; i < 24; i++ ){
	
	
		// Full spectrum
		if ( i == 0 && ALL_ROWS == 1 ){
			// Plot spectrum
			c_ex_full = new TCanvas( "c_ex_full", "Full excitation spectrum",  C_WIDTH, C_HEIGHT );
			GlobSetCanvasMargins( c_ex_full );
			h_ex_full->Draw();
			spec_name = Form( "%s/posXXX_ex_full", print_dir.Data() );
			
			// Print spectrum if desired
			if ( SW_EX[1] == 1 ){
				PrintAll( c_ex_full, spec_name );
				if ( PRINT_ROOT == 1 ){ f->cd(); h_ex_full->Write(); }
			}
			
			// Write SPE file if desired
			if ( SW_EX[2] == 1 ){ WriteSPE( h_ex_full->GetName(), Form( "%s", spec_name.Data() ) ); }
			
		}
		
		
		// Row by row
		if ( ( i == ROW_NUMBER || ( ROW_NUMBER == -1 && i < 6 ) ) && ROW_BY_ROW == 1 ){
			// Plot spectrum
			c_ex_rbr[i] = new TCanvas( Form( "c_ex_rbr_%i", i ), Form( "Ex | Row %i", i ), C_WIDTH, C_HEIGHT );
			GlobSetCanvasMargins( c_ex_rbr[i] );
			h_ex_rbr[i]->Draw();
			spec_name = Form( "%s/posXXX_ex_rbr_%i", print_dir.Data(), i );
			
			// Print spectrum if desired
			if ( SW_EX[1] == 1 ){
				PrintAll( c_ex_rbr[i], spec_name );
				if ( PRINT_ROOT == 1 ){ f->cd(); h_ex_rbr[i]->Write(); }
			}
			
			// Write SPE file if desired
			if ( SW_EX[2] == 1 ){ WriteSPE( h_ex_rbr[i]->GetName(), Form( "%s", spec_name.Data() ) ); }
			
		}
		
		
		// Detector by detector
		if ( ( i == DET_NUMBER || ( DET_NUMBER == -1 ) ) && DET_BY_DET == 1 ){
			spec_name = Form( "%s/posXXX_ex_dbd_%i", print_dir.Data(), i );
			if ( CANVAS_COMBINE == 0 ){
				// Plot spectrum
				c_ex_dbd[i] = new TCanvas( Form( "c_ex_dbd_%i", i ), Form( "Ex | Det %i", i ), C_WIDTH, C_HEIGHT );
				GlobSetCanvasMargins( c_ex_dbd[i] );
				h_ex_dbd[i]->Draw();
				gStyle->SetTitleFont(62);
				
				// Print spectrum if desired
				if ( SW_EX[1] == 1 ){
					PrintAll( c_ex_dbd[i], spec_name );
				}
				
			}
			else{
				c_ex_dbd_comb->cd(i+1);
				h_ex_dbd[i]->SetTitle( Form( "Det #%i", i ) );
				h_ex_dbd[i]->Draw();
				TPad* pad = (TPad*)c_ex_dbd_comb->GetPad(i+1);
				SetCanvasTitleFont( pad );
				
				// Print spectrum if desired
				if ( SW_EX[1] == 1 && i == 23 ){
					PrintAll( c_ex_dbd_comb, Form( "%s/posXXX_ex_dbd_comb", print_dir.Data() ) );
				}
				
				
			}
			
			// Write ROOT file if desired
			if ( PRINT_ROOT == 1 && SW_EX[1] == 1 ){ f->cd(); h_ex_dbd[i]->Write(); }
			
			// Write SPE file if desired
			if ( SW_EX[2] == 1 ){ WriteSPE( h_ex_dbd[i]->GetName(), Form( "%s", spec_name.Data() ) ); }
		} 
	
	
	}
	// Close the root file
	if ( SW_EX[1] == 1 && PRINT_ROOT == 1 && f != NULL ){ if ( f->IsOpen() ){ f->Close(); } }
	return;
}

#endif
