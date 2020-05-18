// ATH_ex_compare.h
// Compare excitation energy spectra for given cuts on a row-by-row basis
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


#ifndef ATH_EX_COMPARE_H_
#define ATH_EX_COMPARE_H_

// Switch for this is SW_EX_COMPARE
TH1F* h_ex_compare1[6];		// Singles spectrum
TH1F* h_ex_compare2[6];		// Clean spectrum

TH1F* h_ex_dbd_singl[24];		// Singles spectrum
TH1F* h_ex_dbd_clean[24];		// Clean spectrum


void HCreateExCompare(){
	// *LOOP* over rows
	for ( Int_t i = 0; i < 6; i++ ){
	
		// Check row number
		if ( ROW_BY_ROW == 1 && ( i == ROW_NUMBER || ROW_NUMBER == -1 ) ){
		
			// Create histograms
			CreateExSpectrum( h_ex_compare1[i], Form( "h_ex_sing1_%i", i ) );
			CreateExSpectrum( h_ex_compare2[i], Form( "h_ex_sing2_%i", i ) );
			h_ex_compare1[i]->SetLineColor(kBlack);
			h_ex_compare2[i]->SetLineColor(kRed);
			
			GlobSetHistFonts( h_ex_compare1[i] );
			GlobSetHistFonts( h_ex_compare2[i] );
		}
	}
	
	for ( Int_t i = 0; i < 24; i++ ){
	
		// Check detector number
		if ( DET_BY_DET == 1 && ( i == DET_NUMBER || DET_NUMBER == -1 ) ){
			CreateExSpectrum( h_ex_dbd_singl[i], Form( "h_ex_dbd_singl_%i", i ) );
			CreateExSpectrum( h_ex_dbd_clean[i], Form( "h_ex_dbd_clean_%i", i ) );
			h_ex_dbd_singl[i]->SetLineColor(kBlack);
			h_ex_dbd_clean[i]->SetLineColor(kRed);
			GlobSetHistFonts( h_ex_dbd_singl[i] );
			GlobSetHistFonts( h_ex_dbd_clean[i] );
		}
	}
	
	
	return;
}


void HDrawExCompare(){
	// Define some local variables
	TCanvas* c_ex_compare[6];
	TCanvas* c_ex_dbd_compare[24];
	
	TString root_name = Form( "%s/pos%i%s_ex_sing_th%s-%s", print_dir.Data(), ARR_POSITION, ( ROW_NUMBER == -1 ? "" : Form( "_row%i", ROW_NUMBER ) ),  DoubleToString( THETA_LB ).Data(), DoubleToString( THETA_UB ).Data() );
	TFile* f;
	
	// Open root file if desired
	if ( SW_EX_COMPARE[1] == 1 && PRINT_ROOT == 1 ){
		f = new TFile( ( root_name + ".root" ).Data(), "RECREATE" );
	}
	
	
	// *LOOP* over rows (i)
	for ( Int_t i = 0; i < 6; i++ ){
	
		// Check row number
		if ( ROW_BY_ROW == 1 && ( i == ROW_NUMBER || ROW_NUMBER == -1 ) ){
		
			// Define canvas and draw
			c_ex_compare[i] = new TCanvas( Form( "c_ex_sing_%i", i ), Form( "COMPARE HISTOGRAMS | ROW %i", i ),  C_WIDTH, C_HEIGHT );
			GlobSetCanvasMargins( c_ex_compare[i] );
			
			h_ex_compare1[i]->Draw();
			h_ex_compare2[i]->Draw("SAME");
			
			// Define a file name for printing
			TString spec_name = Form( "%s/pos%i_ex_compare_row%02i_th%s-%s", print_dir.Data(), ARR_POSITION, i, DoubleToString( thetaCM_singles_cuts[i][ARR_POSITION-1][0] ).Data(), DoubleToString( thetaCM_singles_cuts[i % 6][ARR_POSITION-1][1] ).Data() );
			
			
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
			
		}	
	
	}	// *LOOP* over rows (i)
	
	
	// *LOOP* over detectors (i)
	for ( Int_t i = 0; i < 24; i++ ){
		
		if ( DET_BY_DET == 1 && ( i == DET_NUMBER || DET_NUMBER == -1 ) ){
			c_ex_dbd_compare[i] = new TCanvas( Form( "c_ex_dbd_compare_%i", i ), Form( "COMPARE EXCITATION SPECTRA | DET %i", i ), C_WIDTH, C_HEIGHT );
			GlobSetCanvasMargins( c_ex_dbd_compare[i] );
			h_ex_dbd_singl[i]->Draw();
			h_ex_dbd_clean[i]->Draw("SAME");
			
			TString spec_name = Form( "%s/pos%i_ex_compare_det%02i_th%s-%s", print_dir.Data(), ARR_POSITION, i, DoubleToString( thetaCM_singles_cuts[i % 6][ARR_POSITION-1][0] ).Data(), DoubleToString( thetaCM_singles_cuts[i % 6][ARR_POSITION-1][1] ).Data() );
			
			// Print spectrum if desired
			if ( SW_EX_COMPARE[1] == 1 ){
				PrintAll( c_ex_dbd_compare[i], spec_name );
				if ( PRINT_ROOT == 1 ){
					f->cd();
					h_ex_dbd_singl[i]->Write();
					h_ex_dbd_clean[i]->Write();
				}
			}
			
			if ( SW_EX_COMPARE[2] == 1 ){
				WriteSPE( h_ex_dbd_singl[i]->GetName(), Form( "%s", spec_name.Data() ) );
			}
		}
	}
	
	
	// Close the root file
	if ( SW_EX_COMPARE[1] == 1 && PRINT_ROOT == 1 && f != NULL ){ if ( f->IsOpen() ){ f->Close(); } }
	return;
}

#endif
