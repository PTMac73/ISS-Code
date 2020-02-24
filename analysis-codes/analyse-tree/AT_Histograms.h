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
#include <TH1.h>
#include <TH2.h>
#include <TString.h>

#include <iostream>

#include "AT_Settings.h"
#include "WriteSPE.h"


// --------------------------------------------------------------------------------------------- //
// GLOBAL FUNCTIONS
// Set canvas margins
void GlobSetCanvasMargins( TCanvas *c, Double_t l = 0.1, Double_t r = 0.02, Double_t t = 0.02, Double_t b = 0.1 ){
	TPad* pad = (TPad*)c;
	pad->SetLeftMargin( l );
	pad->SetRightMargin( r );
	pad->SetTopMargin( t );
	pad->SetBottomMargin( b );
	return;
}

// Set histogram fonts
void GlobSetHistFonts( TH1* h ){
	h->GetXaxis()->SetTitleFont(62);
	h->GetXaxis()->SetLabelFont(62);
	h->GetYaxis()->SetTitleFont(62);
	h->GetYaxis()->SetLabelFont(62);
	return;
}

void GlobSetHistFonts( TH2* h ){
	h->GetXaxis()->SetTitleFont(62);
	h->GetXaxis()->SetLabelFont(62);
	h->GetYaxis()->SetTitleFont(62);
	h->GetYaxis()->SetLabelFont(62);
	return;
}

// Write error message about spe files
void ErrorSPE( TString str = "This mode" ){
	std::cout << str << " cannot have a valid .spe file" << "\n";
}


// --------------------------------------------------------------------------------------------- //
// Compare excitation energy spectra for given cuts on a row-by-row basis
// Switch for this is SW_EX_COMPARE
TH1F* h_ex_compare1[6];
TH1F* h_ex_compare2[6];


void HCreateExCompare(){
	// *LOOP* over rows
	for ( Int_t i = 0; i < 6; i++ ){
		if ( i == ROW_NUMBER || ROW_NUMBER == -1 ){
			h_ex_compare1[i] = new TH1F( Form( "h_ex_compare1_%i", i ), "", 450, -1, 8 );
			h_ex_compare2[i] = new TH1F( Form( "h_ex_compare2_%i", i ), "", 450, -1, 8 );
			
			h_ex_compare1[i]->SetTitle("");
			h_ex_compare1[i]->GetXaxis()->SetTitle("Excitation Energy (MeV)");
			h_ex_compare1[i]->GetYaxis()->SetTitle("Counts per 20 keV");
			h_ex_compare1[i]->SetLineColor(kBlack);
			
			h_ex_compare2[i]->SetLineColor(kRed);
			
			GlobSetHistFonts( h_ex_compare1[i] );
		}
	}
}


void HDrawExCompare(){
	TCanvas *c_ex_compare[6];
	for ( Int_t i = 0; i < 6; i++ ){
		if ( i == ROW_NUMBER || ROW_NUMBER == -1 ){
			c_ex_compare[i] = new TCanvas( Form( "c_ex_compare_%i", i ), Form( "COMPARE HISTOGRAMS | ROW %i", i ),  C_WIDTH, C_HEIGHT );
			GlobSetCanvasMargins( c_ex_compare[i] );
			h_ex_compare1[i]->Draw();
			h_ex_compare2[i]->Draw("SAME");
			
			if ( SW_EX_COMPARE[1] == 1 ){
				// PRINT ME
				std::cout << "DEVELOP PRINTING" << "\n"; /*TODO*/
			}
			if ( SW_EX_COMPARE[2] == 1 ){
				// SPE ME
				//WriteSPE( h_ex_compare1[i], "" ); /*TODO*/
			}
		}
	}
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
	
	if ( SW_RDT_CUTS[1] == 1 ){
		// PRINT ME
		std::cout << "DEVELOP PRINTING" << "\n"; /*TODO*/
	}
	if ( SW_RDT_CUTS[2] == 1 ){
		ErrorSPE( "c_rdt_cuts" );
	}
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
	
	if ( SW_EVZ_COMPARE[1] == 1 ){
		// PRINT ME
		std::cout << "DEVELOP PRINTING" << "\n"; /*TODO*/
	}
	if ( SW_EVZ_COMPARE[2] == 1 ){
		ErrorSPE( "c_evz_compare_sing and c_evz_compare_clean" );
	}

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
	
	if ( SW_EVZ[1] == 1 ){
		// PRINT ME
		std::cout << "DEVELOP PRINTING" << "\n"; /*TODO*/
	}
	if ( SW_EVZ[2] == 1 ){
		ErrorSPE( "c_evz" );
	}
}



































#endif
