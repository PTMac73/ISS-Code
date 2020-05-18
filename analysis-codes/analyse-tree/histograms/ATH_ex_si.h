// Defines all the histogram formatting and functions etc for the AnalyseTree.C script
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


#ifndef ATH_EX_SI_H_
#define ATH_EX_SI_H_

// Switch for this is SW_EX_SI
// 0 is the Mg spectrum
// 1 is the Si spectrum
TH1F* h_ex_si[2];

void HCreateExSi(){
	for ( Int_t i = 0; i < 2; i++ ){
		h_ex_si[i] = new TH1F( Form( "h_ex_si_%i", i ), "", 450, -1, 8 );
		GlobSetHistFonts( h_ex_si[i] );
		
		if ( i == 0 ){
			h_ex_si[i]->SetFillColor( kGray );
			h_ex_si[i]->SetLineColor( kGray );
			h_ex_si[i]->GetXaxis()->SetTitle( "Excitation Energy (MeV)" );
			h_ex_si[i]->GetYaxis()->SetTitle( "Counts per 20 keV" );
		}
		else if ( i == 1 ){
			h_ex_si[i]->SetFillColor( kRed );
			h_ex_si[i]->SetLineColor( kBlack );
		}
		
	}
	
	return;
}

void HDrawExSi(){
	TCanvas *c_ex_si = new TCanvas( "c_ex_si", "EX-Si", C_WIDTH, C_HEIGHT );
	GlobSetCanvasMargins( c_ex_si );
	std::cout << "Banana0" << "\n";
	for ( Int_t i = 0; i < 2; i++ ){
		std::cout << "Banana1." << i << "\n";
		h_ex_si[i]->Draw("SAME");
	}
	TString spec_name = Form( "%s/posXXX_exsi", print_dir.Data() );

	// Print spectrum if desired
	if ( SW_EX_SI[1] == 1 ){
		PrintAll( c_ex_si, spec_name );
		
		/*if ( PRINT_ROOT == 1 ){
			out_root_file->cd();
			for ( Int_t i = 0; i < 2; i++ ){
				h_ex_si[i]->Write();
			}
		}*/
		
	}
	
	// Write SPE file if desired
	if ( SW_EX_SI[2] == 1 ){
		WriteSPE( h_ex_si[2]->GetName(), Form( "%s", spec_name.Data() ) );
	}

	return;
}


#endif



























