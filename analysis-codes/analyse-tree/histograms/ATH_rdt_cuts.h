// ATH_rdt_cuts.h
// Recoil detector cuts
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

#include "AT_CutCreator.cxx"
#include "AT_HistogramGlobals.h"
#include "../AT_Settings.h"
#include "WriteSPE.h"


#ifndef ATH_RDT_CUTS_H_
#define ATH_RDT_CUTS_H_

// Switch for this is SW_RDT_CUTS
TH2F* h_rdt_cuts[4];


void HCreateRDTCuts(){
	for ( Int_t i = 0; i < 4; i++ ){
		h_rdt_cuts[i] = new TH2F( Form( "h_rdt_cuts_%i", i ), "", 900, 0, 9000, 1000, 0, 8500 );
		h_rdt_cuts[i]->SetMarkerStyle(20);
		h_rdt_cuts[i]->SetMarkerSize(0.5);
		h_rdt_cuts[i]->SetMarkerColor(kRed);
		h_rdt_cuts[i]->GetYaxis()->SetTitle( Form( "rdt[%i]", i ) );
		h_rdt_cuts[i]->GetXaxis()->SetTitle( Form( "rdt[%i]", i + 4 ) );
		
		GlobSetHistFonts( h_rdt_cuts[i] );
	}
	
	return;
}


void HDrawRDTCuts(){
	TCanvas* c_rdt_cuts_comb;
	TCanvas* c_rdt_cuts[4];
	TObjArray* cuttlefish = new TObjArray();
	TCutG* cuttle;

	// If combining canvases, then divide up
	if ( CANVAS_COMBINE == 1 ){
		TCanvas* c_rdt_cuts_comb = new TCanvas( "c_rdt_cuts_comb", "c_rdt_cuts_comb", C_WIDTH, C_HEIGHT );
		c_rdt_cuts_comb->Divide(2,2);
	}
	
	// *LOOP* over number of spectra
	for ( Int_t i = 0; i < 4; i++ ){
		// Change to the relevant canvas
		if ( CANVAS_COMBINE == 1 ){ c_rdt_cuts_comb->cd(i+1); }
		else{
			c_rdt_cuts[i] = new TCanvas( Form( "c_rdt_cuts_%i", i ), Form( "c_rdt_cuts_%i", i ), C_WIDTH, C_HEIGHT );
			GlobSetCanvasMargins( c_rdt_cuts[i] );	
		}
		
		// Draw the spectrum and the cuts
		if ( DRAW_NEW_CUTS == 1 && CANVAS_COMBINE == 0 ){
			cuttle = CreateCut( h_rdt_cuts[i], c_rdt_cuts[i], Form("cuttle%i", i ) );
			cuttlefish->Add( cuttle );
		}
		else{
			h_rdt_cuts[i]->Draw();
		}
		
		if ( cut_list != NULL ){
			TCutG* cut = (TCutG*)cut_list->At(i);
			cut->Draw("SAME");
		}
	}
	
	// Write to file if new cuts defined
	if ( DRAW_NEW_CUTS == 1 ){
		cuttlefish->SetName("cuttlefish");
		WriteCutFile( cuttlefish );
	}
	
	// Print spectrum if desired
	if ( SW_RDT_CUTS[1] == 1 ){
		TString spec_name = "";
		
		// Combine canvas print settings
		if ( CANVAS_COMBINE == 1 ){
			spec_name = Form( "%s/posXXX_rdt_cuts_all", print_dir.Data() );
			PrintAll( c_rdt_cuts_comb, spec_name );
		}
		
		// Separate canvas print settings
		else{
			for ( Int_t i = 0; i < 4; i++ ){
				spec_name = Form( "%s/posXXX_rdt_cuts_%i", print_dir.Data(), i );
				PrintAll( c_rdt_cuts[i], spec_name );;
			}
		}
		
		// Root file options
		if ( PRINT_ROOT == 1 ){
			for ( Int_t i = 0; i < 4; i++ ){
				out_root_file->cd();
				h_rdt_cuts[i]->Write();
			}
		}
	}
	
	// Write SPE file if desired
	if ( SW_RDT_CUTS[2] == 1 ){
		ErrorSPE( "c_rdt_cuts_comb" );
	}
	
	return;
}



#endif
