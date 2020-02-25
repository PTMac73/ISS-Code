// AT_Histograms.h
// Defines all the histogram formatting and functions etc for the AnalyseTree.C script
// ============================================================================================= //
// Patrick MacGregor
// Nuclear Physics Research Group
// Department of Physics and Astronomy
// The University of Manchester
// ============================================================================================= //
#ifndef AT_HISTOGRAMS_H_
#define AT_HISTOGRAMS_H_

#include <TCanvas.h>
#include <TCutG.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TString.h>

#include <iostream>

#include "AT_Settings.h"
#include "AT_HistogramGlobals.h"
#include "WriteSPE.h"

// --------------------------------------------------------------------------------------------- //
// Compare excitation energy spectra for given cuts on a row-by-row basis
// Switch for this is SW_EX_COMPARE
TH1F* h_ex_compare1[6];		// Singles spectrum
TH1F* h_ex_compare2[6];		// Clean spectrum


void HCreateExCompare(){
	// *LOOP* over rows
	for ( Int_t i = 0; i < 6; i++ ){
	
		// Check row number
		if ( i == ROW_NUMBER || ROW_NUMBER == -1 ){
		
			// Create histograms
			h_ex_compare1[i] = new TH1F( Form( "h_ex_compare1_%i", i ), "", 450, -1, 8 );
			h_ex_compare2[i] = new TH1F( Form( "h_ex_compare2_%i", i ), "", 450, -1, 8 );
			
			// Format histograms
			h_ex_compare1[i]->SetTitle("");
			h_ex_compare1[i]->GetXaxis()->SetTitle("Excitation Energy (MeV)");
			h_ex_compare1[i]->GetYaxis()->SetTitle("Counts per 20 keV");
			
			h_ex_compare1[i]->SetLineColor(kBlack);
			h_ex_compare2[i]->SetLineColor(kRed);
			
			GlobSetHistFonts( h_ex_compare1[i] );
		}
	}
	
	return;
}


void HDrawExCompare(){
	// Define some local variables
	TCanvas* c_ex_compare[6];
	
	TString root_name = Form( "%s/ex_comp_th%s-%s_posXXX%s", print_dir.Data(), DoubleToString( THETA_LB ).Data(), DoubleToString( THETA_UB ).Data(), ( ROW_NUMBER == -1 ? "" : Form( "_row%i", ROW_NUMBER ) ) );
	TFile* f;
	
	// Open root file if desired
	if ( SW_EX_COMPARE[1] == 1 && PRINT_ROOT == 1 ){
		f = new TFile( ( root_name + ".root" ).Data(), "RECREATE" );
	}
	
	
	// *LOOP* over rows (i)
	for ( Int_t i = 0; i < 6; i++ ){
	
		// Check row number
		if ( i == ROW_NUMBER || ROW_NUMBER == -1 ){
		
			// Define canvas and draw
			c_ex_compare[i] = new TCanvas( Form( "c_ex_compare_%i", i ), Form( "COMPARE HISTOGRAMS | ROW %i", i ),  C_WIDTH, C_HEIGHT );
			GlobSetCanvasMargins( c_ex_compare[i] );
			
			h_ex_compare1[i]->Draw();
			h_ex_compare2[i]->Draw("SAME");
			
			// Define a file name for printing
			TString spec_name = ( ROW_NUMBER == -1 ? root_name + Form( "_row%i", i ) : root_name );
			
			
			// Print spectrum if desired
			if ( SW_EX_COMPARE[1] == 1 ){
				PrintAll( c_ex_compare[i], spec_name );
				if ( PRINT_ROOT == 1 ){
					f->cd();
					h_ex_compare1[i]->Write();
					h_ex_compare2[i]->Write();
				}
			}
			
			// Write SPE file if desired
			if ( SW_EX_COMPARE[2] == 1 ){
				WriteSPE( h_ex_compare1[i]->GetName(), Form( "%s", spec_name.Data() ) );
			}
			
		}	// *LOOP* over rows (i)
	
	}
	// Close the root file
	if ( SW_EX_COMPARE[1] == 1 && PRINT_ROOT == 1 && f != NULL ){ if ( f->IsOpen() ){ f->Close(); } }
	return;
}



// --------------------------------------------------------------------------------------------- //
// Recoil detector cuts
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
	TCanvas* c_rdt_cuts = new TCanvas( "c_rdt_cuts", "c_rdt_cuts", C_WIDTH, C_HEIGHT );
	c_rdt_cuts->Divide(2,2);
	for ( Int_t i = 0; i < 4; i++ ){
		c_rdt_cuts->cd(i+1);
		h_rdt_cuts[i]->Draw();
		if ( cut_list != NULL ){
			TCutG* cut = (TCutG*)cut_list->At(i);
			cut->Draw("SAME");
		}
	}
	
	// Print spectrum if desired
	if ( SW_RDT_CUTS[1] == 1 ){
		TString spec_name = Form( "%s/rdt_cuts_posXXX", print_dir.Data() );
		PrintAll( c_rdt_cuts, spec_name );
		if ( PRINT_ROOT == 1 ){
			for ( Int_t i = 0; i < 4; i++ ){
				out_root_file->cd();
				h_rdt_cuts[i]->Write();
			}
		}
	}
	
	// Write SPE file if desired
	if ( SW_RDT_CUTS[2] == 1 ){
		ErrorSPE( "c_rdt_cuts" );
	}
	
	return;
}

// --------------------------------------------------------------------------------------------- //
// Compare E v.s. z plots for different angle cuts
// Switch for this is SW_EVZ_COMPARE
TH2F* h_evz_compare[4];


void HCreateEVZCompare(){
	for ( Int_t i = 0; i < 4; i++ ){
		h_evz_compare[i] = new TH2F( Form( "h_evz_compare_%i", i ), "", 450, -50, -10, 900, 0, 9 );
		h_evz_compare[i]->SetMarkerStyle(20);
		h_evz_compare[i]->SetMarkerSize(0.5);
		h_evz_compare[i]->SetMarkerColor( i % 2 == 0 ? kRed : kBlack );
		h_evz_compare[i]->GetYaxis()->SetTitle( "Energy (MeV)" );
		h_evz_compare[i]->GetXaxis()->SetTitle( "z (cm)" );
		GlobSetHistFonts( h_evz_compare[i] );
	}
	
	return;
}


void HDrawEVZCompare(){
	TCanvas *c_evz_compare_sing = new TCanvas( "c_evz_compare_sing", "COMPARING EVZ AT DIFFERENT ANGLES -- SINGLES", C_WIDTH, C_HEIGHT );
	GlobSetCanvasMargins( c_evz_compare_sing );
	h_evz_compare[0]->Draw();
	h_evz_compare[1]->Draw("SAME");
	
	TCanvas *c_evz_compare_clean = new TCanvas( "c_evz_compare_clean", "COMPARING EVZ AT DIFFERENT ANGLES -- CLEAN", C_WIDTH, C_HEIGHT );
	GlobSetCanvasMargins( c_evz_compare_clean );
	h_evz_compare[2]->Draw();
	h_evz_compare[3]->Draw("SAME");
	
	
	
	// Print spectra if desired
	if ( SW_EVZ_COMPARE[1] == 1 ){
		TString base = "%s/evz_comp%s_th_%3.2f-%3.2f_posXXX";
		TString root_name = Form( base.Data(), print_dir.Data(), "", THETA_LB, THETA_UB );
		TString spec_name_sing = Form( base.Data(), print_dir.Data(), "_sing", THETA_LB, THETA_UB );
		PrintAll( c_evz_compare_sing, spec_name_sing );
		
		TString spec_name_clean = Form( base.Data(), print_dir.Data(), "_clean", THETA_LB, THETA_UB );
		PrintAll( c_evz_compare_clean, spec_name_clean );
		
		if ( PRINT_ROOT == 1 ){
			TFile* f = new TFile( ( root_name + ".root" ).Data(), "RECREATE" );
			for ( Int_t i = 0; i < 4; i++ ){
				f->cd();
				h_evz_compare[i]->Write();
			}
			f->Close();
		}
	}
	
	// Write SPE file if desired
	if ( SW_EVZ_COMPARE[2] == 1 ){
		ErrorSPE( "c_evz_compare_sing and c_evz_compare_clean" );
	}

	return;
}

// --------------------------------------------------------------------------------------------- //
// Full E v.s. z plot
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
