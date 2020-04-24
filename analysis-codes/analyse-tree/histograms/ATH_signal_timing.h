// ATH_signal_timing.h
// Signals plotted against time
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


#ifndef ATH_SIGNAL_TIMING_H_
#define ATH_SIGNAL_TIMING_H_

// Switch for this is SW_SIGTIME
TH2F* h_sigtime_e[24];		// Full spectrum
// *TODO* the rest of the spectra

// Define which ones to print
const Int_t NUM_HISTS_SIGTIME = 7;
Bool_t sigtime_print_opt[NUM_HISTS_SIGTIME] = {
	1,	// (0) e
	0,	// (1) xn
	0,	// (2) xf
	0,	// (3) rdt
	0,	// (4) tac
	0,	// (5) elum
	0	// (6) ezero
};

// ------------------------------------{ 00, 01, 02, 03, 04, 05, 06 };
TString sigtime_spec_name[NUM_HISTS_SIGTIME] = { "", "", "", "", "", "", "" };

// Populate TString array of names
void MakeSpecNameSIGTIME( Int_t i ){
	sigtime_spec_name[0]  = Form( "%s/posXXX_sigtime_e_%i", print_dir.Data(), i );		// (0) e
	//sigtime_spec_name[1]  = Form( "%s/posXXX_xnxf_p_%i", print_dir.Data(), i );		// (1) xn
	//sigtime_spec_name[2]  = Form( "%s/posXXX_xnxfE_%i", print_dir.Data(), i );		// (2) xf
	//sigtime_spec_name[3]  = Form( "%s/posXXX_xnxfE_p_%i", print_dir.Data(), i );	// (3) rdt
	//sigtime_spec_name[4]  = Form( "%s/posXXX_xnE_%i", print_dir.Data(), i );		// (4) tac
	//sigtime_spec_name[5]  = Form( "%s/posXXX_xfE_%i", print_dir.Data(), i );		// (5) elum
	//sigtime_spec_name[6]  = Form( "%s/posXXX_xnxf_col_%i", print_dir.Data(), i );	// (6) ezero
}




void HCreateSIGTIME(){
	// *LOOP* over detectors
	for ( Int_t i = 0; i < 24; i++ ){
		
		// Initially declare as NULL
		h_sigtime_e[i] = NULL;	// (0) e
	
	
		// Detector by detector
		if ( ( i == DET_NUMBER || ( DET_NUMBER == -1 ) ) && det_array[i % 6][(Int_t)TMath::Floor(i/6)] != 0 ){
			// (0) e
			h_sigtime_e[i] = new TH2F( Form( "h_sigtime_e_%i", i ), Form( "E SIGTIME | Det %i", i ), 500, 3E11, 8E11, 2000, 0, 2000 );
			GlobCreate2DHists( h_sigtime_e[i], "Time", "Raw Energy Signal" );
			h_sigtime_e[i]->SetMarkerSize(0.1);
		}
	}
	
	return;
}


void HDrawSIGTIME(){
	// Define some local variables
	TCanvas* c_sigtime_e_comb;
	TCanvas* c_sigtime_e[24];
	
	if ( CANVAS_COMBINE == 1 ){
		c_sigtime_e_comb = new TCanvas( "c_sigtime_e_comb", "SIGTIME E Spectrum Combined", 2*C_WIDTH, 2*C_HEIGHT );
		c_sigtime_e_comb->Divide(6,4);
	}
	
	TString root_name = Form( "%s/posXXX_sigtime", print_dir.Data() );
	TFile* f;

	
	// Open root file if desired
	if ( SW_SIGTIME[1] == 1 && PRINT_ROOT == 1 ){
		f = new TFile( ( root_name + ".root" ).Data(), "RECREATE" );
	}
	
	
	// *LOOP* over detectors (i)
	for ( Int_t i = 0; i < 24; i++ ){

		// Detector by detector
		if ( ( i == DET_NUMBER || DET_NUMBER == -1 ) && det_array[i % 6][(Int_t)TMath::Floor(i/6)] != 0 ){
			
			MakeSpecNameSIGTIME( i );
			
			if ( CANVAS_COMBINE == 0 ){
				// (0) e
				c_sigtime_e[i] = new TCanvas( Form( "c_sigtime_%i", i ), Form( "sigtime | Det %i", i ), C_WIDTH, C_HEIGHT );
				GlobSetCanvasMargins( c_sigtime_e[i] );
				h_sigtime_e[i]->Draw();
				
				// Print spectrum if desired
				if ( SW_SIGTIME[1] == 1 ){
					if ( sigtime_print_opt[0] == 1 ){ PrintAll( c_sigtime_e[i], sigtime_spec_name[0] ); }
				}
				
			}
			else{
				c_sigtime_e_comb->cd(i+1);
				h_sigtime_e[i]->SetTitle( Form( "Det #%i", i ) );
				h_sigtime_e[i]->Draw();
				//TPad* pad = (TPad*)c_sigtime_e_comb->GetPad(i+1);
				//SetCanvasTitleFont( pad );
				
				// Print spectrum if desired
				if ( SW_SIGTIME[1] == 1 && i == 23 ){
					if ( sigtime_print_opt[0] == 1 ){ PrintAll( c_sigtime_e_comb, Form( "%s/posXXX_sigtime_comb", print_dir.Data() ) ); }
				}
				
				
			}
			
			// Write ROOT file if desired
			if ( PRINT_ROOT == 1 && SW_SIGTIME[1] == 1 ){ f->cd(); h_sigtime_e[i]->Write(); }
			
			// Write SPE file if desired
			if ( SW_SIGTIME[2] == 1 ){ ErrorSPE( "SIGTIME spectra" ); }
		} 
	
	
	}
	// Close the root file
	if ( SW_SIGTIME[1] == 1 && PRINT_ROOT == 1 && f != NULL ){ if ( f->IsOpen() ){ f->Close(); } }
	return;
}

#endif
