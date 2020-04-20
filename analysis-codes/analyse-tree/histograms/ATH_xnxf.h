// ATH_xnxf.h
// XNXF plots
// ============================================================================================= //
// Patrick MacGregor
// Nuclear Physics Research Group
// Department of Physics and Astronomy
// The University of Manchester
// ============================================================================================= //

#include <TCanvas.h>
#include <TCutG.h>
#include <TF1.h>
#include <TFile.h>
#include <TFitResult.h>
#include <TH1.h>
#include <TH2.h>
#include <TLine.h>
#include <TProfile.h>
#include <TString.h>
#include <TStyle.h>

#include <iostream>

#include "AT_HistogramGlobals.h"
#include "../AT_Settings.h"
#include "WriteSPE.h"


#ifndef ATH_XNXF_H_
#define ATH_XNXF_H_

// Switch for this is SW_XNXF
// DEFINE HISTOGRAMS/PROFILES
const Int_t NUM_DIFF_HISTS = 10;
const Int_t NUM_DETS = 24;
const Double_t XNXF_FRAC = 0.1;
TH2F* h_xnxf[NUM_DETS];					// (0) XN-XF hist monochrome
TProfile* p_xnxf[NUM_DETS];				// (1) XN-XF profile
TH2F* h_xnxfE[NUM_DETS];				// (2) XNXF-E hist monochrome
TProfile* p_xnxfE[NUM_DETS];			// (3) XNXF-E profile
TH2F* h_xnE[NUM_DETS];					// (4) XN-E hist monochrome
TH2F* h_xfE[NUM_DETS];					// (5) XF-E hist monochrome
TH2F* h_xnxf_colour[NUM_DETS][4];		// (6) XN-XF hist coloured
TH2F* h_xnE_colour[NUM_DETS][2];		// (7) XN-E hist coloured
TH2F* h_xfE_colour[NUM_DETS][2];		// (8) XF-E hist coloured
TH2F* h_xnxfE_colour[NUM_DETS][3];		// (9) XNXF-E hist coloured

// Define which ones to print
Bool_t xnxf_print_opt[NUM_DIFF_HISTS] = {
	0,	// (0) XN-XF hist monochrome
	1,	// (1) XN-XF profile
	1,	// (2) XNXF-E hist monochrome
	1,	// (3) XNXF-E profile
	0,	// (4) XN-E hist monochrome
	0,	// (5) XF-E hist monochrome
	1,	// (6) XN-XF hist coloured
	1,	// (7) XN-E hist coloured
	1,	// (8) XF-E hist coloured
	1	// (9) XNXF-E hist coloured
};

// -------------------------------- { 00, 01, 02, 03, 04, 05, 06, 07, 08, 09 };
TString spec_name[NUM_DIFF_HISTS] = { "", "", "", "", "", "", "", "", "", "" };

// Additional variables
TString xnxfcut_dir = "/home/ptmac/Documents/07-CERN-ISS-Mg/Mg-Analysis/SPE-Files/xnxf-calibration/pos1_xnxf_cuts.root";
TCutG* xnxf_cut;

// Store the correction parameters
Double_t xn_corr[NUM_DETS];
Double_t xfxne_corr[NUM_DETS][2];



// FUNCTIONS ----------------------------------------------------------------------------------- //
// Generate XNXF CUT NAME
TString XNXFCutName( Int_t i ){
	return Form( "xnxf_cut_%i", i );
}

// Generate cut array from a file
TObjArray* GetXNXFCutArray( TFile* f ){
	TCutG* cut;
	TObjArray* arr = new TObjArray();
	for ( Int_t i = 0; i < NUM_DETS; i++ ){
		cut = (TCutG*)f->Get( XNXFCutName(i) );
		arr->Add( cut );
	}
	return arr;
}

// Print parameters XN-XF
void PrintXNXFFitPars( TFitResultPtr r ){
	std::cout << "--------------------------------------------------------------------------------\n";
	std::cout << "Fit: XN = A * XF + B" << "\n";
	std::cout << Form( "A = %8.6f +/- %8.6f", r->Parameter(1), r->ParError(1) ) << "\n";
	std::cout << Form( "B = %8.6f +/- %8.6f", r->Parameter(0), r->ParError(0) ) << "\n";
	std::cout << Form( "Chi^2 = %8.6f", r->Chi2() ) << "\n";
	std::cout << Form( "\u2192 %8.6f * XN = - XF + %8.6f", -1.0/r->Parameter(1), -1.0*r->Parameter(0)/r->Parameter(1) ) << "\n";
}

// Print parameters XNXF-E
void PrintXNXFEFitPars( TFitResultPtr r ){
	std::cout << "--------------------------------------------------------------------------------\n";
	std::cout << "Fit: E = A * X + B, where X is calibrated sum of XN and XF" << "\n";
	std::cout << Form( "A = %8.6f +/- %8.6f", r->Parameter(1), r->ParError(1) ) << "\n";
	std::cout << Form( "B = %8.6f +/- %8.6f", r->Parameter(0), r->ParError(0) ) << "\n";
	std::cout << Form( "Chi^2 = %8.6f", r->Chi2() ) << "\n";
}

// Print the correction parameters XN-XF
void PrintXNCorr(){
	std::cout << "Double_t xnCorr[24] = {" << "\n";
	for ( Int_t i = 0; i < NUM_DETS; i++ ){
		std::cout << "\t" << std::right << std::fixed << std::setw(8) << std::setprecision(6) << xn_corr[i] << ",\t// " << Form("%02d", i ) << "\n";
	}	
	std::cout << "};" << "\n";
	return;
}

// Print the correction parameters XNXF-E
void PrintXNXFECorr(){
	std::cout << "Double_t xfxneCorr[24][2] = {" << "\n";
	for ( Int_t i = 0; i < NUM_DETS; i++ ){
		std::cout << "\t{ " << std::right << std::fixed << std::setw(10) << std::setprecision(6) << xfxne_corr[i][0] << ", " << std::fixed << std::setw(8) << std::setprecision(6) << xfxne_corr[i][1] << " },\t// " << Form("%02d", i ) << "\n";
	}	
	std::cout << "};" << "\n";
	return;
}

// Format 2D histograms
void Create2DHists( TH2F* h, TString x_label, TString y_label ){
	h->SetTitle("");
	h->GetXaxis()->SetTitle( x_label.Data() );
	h->GetYaxis()->SetTitle( y_label.Data() );
	h->SetMarkerStyle(20);
	h->SetMarkerSize(0.5);
	h->SetMarkerColor(kRed);	
	GlobSetHistFonts( h );
	return;
}

// Format profiles
void CreateProfile( TProfile* p, TString x_label, TString y_label ){
	p->SetTitle("");
	p->GetXaxis()->SetTitle( x_label.Data() );
	p->GetYaxis()->SetTitle( y_label.Data() );
	GlobSetHistFonts( p );
	return;
}


// Populate TString array of names
void MakeSpecName( Int_t i ){
	spec_name[0] = Form( "%s/posXXX_xnxf_%i", print_dir.Data(), i );		// (0) XN-XF hist monochrome
	spec_name[1] = Form( "%s/posXXX_xnxf_p_%i", print_dir.Data(), i );		// (1) XN-XF profile
	spec_name[2] = Form( "%s/posXXX_xnxfE_%i", print_dir.Data(), i );		// (2) XNXF-E hist monochrome
	spec_name[3] = Form( "%s/posXXX_xnxfE_p_%i", print_dir.Data(), i );		// (3) XNXF-E profile
	spec_name[4] = Form( "%s/posXXX_xnE_%i", print_dir.Data(), i );			// (4) XN-E hist monochrome
	spec_name[5] = Form( "%s/posXXX_xfE_%i", print_dir.Data(), i );			// (5) XF-E hist monochrome
	spec_name[6] = Form( "%s/posXXX_xnxf_col_%i", print_dir.Data(), i );	// (6) XN-XF hist coloured
	spec_name[7] = Form( "%s/posXXX_xnE_col_%i", print_dir.Data(), i );		// (7) XN-E hist coloured
	spec_name[8] = Form( "%s/posXXX_xfE_col_%i", print_dir.Data(), i );		// (8) XF-E hist coloured
	spec_name[9] = Form( "%s/posXXX_xnxfE_col_%i", print_dir.Data(), i );	// (9) XNXF-E hist coloured
}

// MAIN FUNCTIONS ------------------------------------------------------------------------------ //
void HCreateXNXF(){

	// *LOOP* over detectors
	for ( Int_t i = 0; i < 24; i++ ){
	
		// Initially declare them as NULL to help with error detection
		h_xnxf[i] = NULL;		// (0) XN-XF hist monochrome
		p_xnxf[i] = NULL;		// (1) XN-XF profile
		h_xnxfE[i] = NULL;		// (2) XNXF-E hist monochrome
		p_xnxfE[i] = NULL;		// (3) XNXF-E profile
		h_xnE[i] = NULL;		// (4) XN-E hist monochrome
		h_xfE[i] = NULL;		// (5) XF-E hist monochrome
		for ( Int_t j = 0; j < 4; j++ ){ h_xnxf_colour[i][j] = NULL; }	// (6) XN-XF hist coloured
		for ( Int_t j = 0; j < 2; j++ ){ h_xnE_colour[i][j] = NULL; } 	// (7) XN-E hist coloured
		for ( Int_t j = 0; j < 2; j++ ){ h_xfE_colour[i][j] = NULL; }	// (8) XF-E hist coloured
		for ( Int_t j = 0; j < 3; j++ ){ h_xnxfE_colour[i][j] = NULL; }	// (9) XNXF-E hist coloured
		
		// Detector by detector
		if ( i == DET_NUMBER || ( DET_NUMBER == -1 ) ){
		
			// (0) *HIST* XN-XF hist monochrome
			h_xnxf[i] = new TH2F( Form( "h_xnxf_%i", i ), Form( "XN-XF SPECTRUM %i", i ), 200, 0, 2000, 200, 0, 2000 );
			Create2DHists( h_xnxf[i], "XF", "XN" );
			
			// (1) *HIST* XN-XF profile
			p_xnxf[i] = new TProfile( Form( "p_xnxf_%i", i ), Form( "XN-XF SPECTRUM %i", i ), 200, xnxf_lims[i][0], xnxf_lims[i][1], xnxf_lims[i][2], xnxf_lims[i][3] );
			CreateProfile( p_xnxf[i], "XF", "XN" );
			
			// (2) *HIST* XNXF-E hist monochrome
			h_xnxfE[i] = new TH2F( Form( "h_xnxfE_%i", i ), Form( "XNXF-E SPECTRUM %i", i ), 200, 0, 2000, 200, 0, 2000 );
			Create2DHists( h_xnxfE[i], "XNXF", "E" );
			
			// (3) *HIST* XNXF-E profile
			p_xnxfE[i] = new TProfile( Form( "p_xnxfE_%i", i ), Form( "XNXF-E SPECTRUM %i", i ), 200, xnxfE_lims[i][0], xnxfE_lims[i][1], xnxfE_lims[i][2], xnxfE_lims[i][3] );
			CreateProfile( p_xnxfE[i], "XNXF", "E" );
			
			// (4) *HIST* XN-E hist monochrome
			h_xnE[i] = new TH2F( Form( "h_xnE_%i", i ), Form( "XN-E SPECTRUM %i", i ), 200, 0, 2000, 200, 0, 2000 );
			Create2DHists( h_xnE[i], "XN", "E" );
			
			// (5) *HIST* XF-E hist monochrome
			h_xfE[i] = new TH2F( Form( "h_xfE_%i", i ), Form( "XF-E SPECTRUM %i", i ), 200, 0, 2000, 200, 0, 2000 );
			Create2DHists( h_xfE[i], "XF", "E" );
			
			// (6) *HIST* XN-XF hist coloured
			for ( Int_t j = 0; j < 4; j++ ){
				h_xnxf_colour[i][j] = new TH2F( Form( "h_xnxf_%i_%i", i, j ), Form( "XN-XF COLOUR SPECTRUM | Det %i | Case %i", i, j ), 200, 0, 2000, 200, 0, 2000 );
				Create2DHists( h_xnxf_colour[i][j], "XF", "XN" );
			}
			h_xnxf_colour[i][0]->SetMarkerColor( kBlue );
			h_xnxf_colour[i][1]->SetMarkerColor( kBlack );
			h_xnxf_colour[i][2]->SetMarkerColor( kRed );
			h_xnxf_colour[i][3]->SetMarkerColor( kGreen );
			
			// (7) *HIST*  XN-E hist coloured
			for ( Int_t j = 0; j < 2; j++ ){
				h_xnE_colour[i][j] = new TH2F( Form( "h_xnE_%i_%i", i, j ), Form( "XN-E COLOUR SPECTRUM | Det %i | XN %s E", i, ( j == 0 ? "<" : ">" ) ), 200, 0, 2000, 200, 0, 2000 );
				Create2DHists( h_xnE_colour[i][j], "XN", "E" );
			}
			h_xnE_colour[i][0]->SetMarkerColor( kBlue );
			h_xnE_colour[i][1]->SetMarkerColor( kGreen );
			
			// (8) *HIST* XF-E hist coloured
			for ( Int_t j = 0; j < 2; j++ ){
				h_xfE_colour[i][j] = new TH2F( Form( "h_xfE_%i_%i", i, j ), Form( "XF-E COLOUR SPECTRUM | Det %i | XF %s E", i, ( j == 0 ? "<" : ">" ) ), 200, 0, 2000, 200, 0, 2000 );
				Create2DHists( h_xfE_colour[i][j], "XF", "E" );
			}
			h_xfE_colour[i][0]->SetMarkerColor( kGreen );
			h_xfE_colour[i][1]->SetMarkerColor( kBlue );
			
			// (9) XNXF-E hist coloured
			for ( Int_t j = 0; j < 3; j++ ){
				h_xnxfE_colour[i][j] = new TH2F( Form( "h_xnxfE_%i_%i", i, j ), Form( "XNXF-E COLOUR SPECTRUM | Det %i | CASE %i", i, j ), 200, 0, 2000, 200, 0, 2000 );
				Create2DHists( h_xnxfE_colour[i][j], "XNXF", "E" );
			}
			h_xnxfE_colour[i][0]->SetMarkerColor( kRed );
			h_xnxfE_colour[i][1]->SetMarkerColor( kGreen );
			h_xnxfE_colour[i][2]->SetMarkerColor( kBlue );

			// Correction parameters		
			xn_corr[i] = 0.0;
			xfxne_corr[i][0] = 0.0;
			xfxne_corr[i][1] = 0.0;
		}
	}
	
	return;
}


void HDrawXNXF(){
	// DEFINE CANVASES
	TCanvas* c_xnxf[24];			// (0) XN-XF hist monochrome
	TCanvas* c_xnxf_comb;
	
	TCanvas* c_xnxf_p[24];			// (1) XN-XF profile
	TCanvas* c_xnxf_p_comb;
	
	TCanvas* c_xnxfE[24];			// (2) XNXF-E hist monochrome
	TCanvas* c_xnxfE_comb;
	
	TCanvas* c_xnxfE_p[24];			// (3) XNXF-E profile
	TCanvas* c_xnxfE_p_comb;
	
	TCanvas* c_xnE[24];				// (4) XN-E hist monochrome
	TCanvas* c_xnE_comb;
	
	TCanvas* c_xfE[24];				// (5) XF-E hist monochrome
	TCanvas* c_xfE_comb;
	
	TCanvas* c_xnxf_colour[24];		// (6) XN-XF hist coloured
	TCanvas* c_xnxf_colour_comb;

	TCanvas* c_xnE_colour[24];		// (7) XN-E hist coloured
	TCanvas* c_xnE_colour_comb;
	
	TCanvas* c_xfE_colour[24];		// (8) XF-E hist coloured
	TCanvas* c_xfE_colour_comb;
	
	TCanvas* c_xnxfE_colour[24];	// (9) XNXF-E hist coloured
	TCanvas* c_xnxfE_colour_comb;
	
	// DEFINE COMBINED CANVASES
	if ( CANVAS_COMBINE == 1 ){
		c_xnxf_comb = new TCanvas( "c_xnxf_comb", "XNXF PLOTS", C_WIDTH, C_HEIGHT );
		c_xnxf_comb->Divide(6,4);
		
		// *TODO* Finish this
	}
	
	
	// Define boundary lines
	TLine* bline_xnxf[24][4];
	TLine* bline_xnxfE[24][4];
	
	// Define fit lines for the different plots
	TF1* fitline_xnxf = new TF1( "fitline_xnxf", "[0] + [1]*x", 0, 100 );
	fitline_xnxf->SetLineWidth(2);
	fitline_xnxf->SetLineColor(kRed);
	
	TF1* fitline_xnxfE = new TF1( "fitline_xnxfE", "[0] + [1]*x", 0, 100 );
	fitline_xnxfE->SetLineWidth(2);
	fitline_xnxfE->SetLineColor(kRed);
	
	
	// Create names and files
	TString root_name = Form( "%s/posXXX_xnxf", print_dir.Data() );
	TFile* f;
	
	// Open root file if desired
	if ( SW_XNXF[1] == 1 && PRINT_ROOT == 1 ){
		f = new TFile( ( root_name + ".root" ).Data(), "RECREATE" );
	}
	
	
	// *LOOP* over detectors (i)
	for ( Int_t i = 0; i < 24; i++ ){

		// Detector by detector
		if ( ( i == DET_NUMBER || ( DET_NUMBER == -1 ) ) && det_array[i % 6][(Int_t)TMath::Floor(i/6)] != 0 ){
			
			// Define names for the spectra
			MakeSpecName(i);
			
			// Prepare boundary lines
			for ( Int_t j = 0; j < 4; j++ ){
				// For XN-XF plots
				bline_xnxf[i][j] = new TLine( 
					( j == 3 ? xnxf_lims[i][1] : xnxf_lims[i][0] ),
					( j == 2 ? xnxf_lims[i][3] : xnxf_lims[i][2] ),
					( j == 1 ? xnxf_lims[i][0] : xnxf_lims[i][1] ),
					( j == 0 ? xnxf_lims[i][2] : xnxf_lims[i][3] )
				);
				bline_xnxfE[i][j] = new TLine( 
					( j == 3 ? xnxfE_lims[i][1] : xnxfE_lims[i][0] ),
					( j == 2 ? xnxfE_lims[i][3] : xnxfE_lims[i][2] ),
					( j == 1 ? xnxfE_lims[i][0] : xnxfE_lims[i][1] ),
					( j == 0 ? xnxfE_lims[i][2] : xnxfE_lims[i][3] )
				);
			}
			
			// TPROFILE FIT LINES
			// Fit some lines to the XN:XF TProfile plots
			TFitResultPtr fit_ptr = p_xnxf[i]->Fit("pol1", "SQ0", "", xnxf_lims[i][0], xnxf_lims[i][1] );
			fitline_xnxf->SetRange( xnxf_lims[i][0], xnxf_lims[i][1] );
			fitline_xnxf->SetParameters( fit_ptr->Parameter(0), fit_ptr->Parameter(1) );
			xn_corr[i] = -1.0/fit_ptr->Parameter(1);
			//PrintXNXFFitPars( fit_ptr );
			
			// Fit some lines to the XNXF:E TProfile plots
			TFitResultPtr fit_ptrE = p_xnxfE[i]->Fit("pol1", "SQ0", "", xnxfE_lims[i][0], xnxfE_lims[i][1] );
			fitline_xnxfE->SetRange( xnxfE_lims[i][0], xnxfE_lims[i][1] );
			fitline_xnxfE->SetParameters( fit_ptrE->Parameter(0), fit_ptrE->Parameter(1) );
			xfxne_corr[i][0] = 0.5*fit_ptrE->Parameter(0); 	// Intercept
			xfxne_corr[i][1] = fit_ptrE->Parameter(1); 		// Gradient
			//PrintXNXFEFitPars( fit_ptrE );
			
			// Get the XNXF cut
			if ( cut_list_xnxf != NULL ){ xnxf_cut = (TCutG*)cut_list_xnxf->At(i); }
			
			// DRAW HISTOGRAMS ON CANVASES
			if ( CANVAS_COMBINE == 0 ){
			
				// (0) *CANVAS* XN-XF hist monochrome
				c_xnxf[i] = new TCanvas( Form( "c_xnxf_%i", i ), Form( "XN-XF | Det %i", i ), C_WIDTH, C_HEIGHT );
				GlobSetCanvasMargins( c_xnxf[i] );
				h_xnxf[i]->Draw();
				for ( Int_t j = 0; j < 4; j++ ){
					bline_xnxf[i][j]->Draw();
				}
				if ( xnxf_cut != NULL ){ xnxf_cut->Draw("SAME"); }
				
				
				// (1) *CANVAS* XN-XF profile
				c_xnxf_p[i] = new TCanvas( Form( "c_xnxf_p_%i", i ), Form( "XN-XF Profile | Det %i", i ), C_WIDTH, C_HEIGHT );
				GlobSetCanvasMargins( c_xnxf_p[i] );
				p_xnxf[i]->Draw();
				fitline_xnxf->Draw("SAME");
				
				
				// (2) *CANVAS* XNXF-E hist monochrome
				c_xnxfE[i] = new TCanvas( Form( "c_xnxfE_%i", i ), Form( "XNXF-E | Det %i", i ), C_WIDTH, C_HEIGHT );
				GlobSetCanvasMargins( c_xnxfE[i] );
				h_xnxfE[i]->Draw();
				for ( Int_t j = 0; j < 4; j++ ){
					bline_xnxfE[i][j]->Draw();
				}
				
				
				// (3) *CANVAS* XNXF-E profile
				c_xnxfE_p[i] = new TCanvas( Form( "c_xnxf_pE_%i", i ), Form( "XNXF-E Profile | Det %i", i ), C_WIDTH, C_HEIGHT );
				GlobSetCanvasMargins( c_xnxfE_p[i] );
				p_xnxfE[i]->Draw();
				fitline_xnxfE->Draw("SAME");
				
				
				// (4) *CANVAS* XN-E hist monochrome
				c_xnE[i] = new TCanvas( Form( "c_xnE_%i", i ), Form( "XN-E | Det %i", i ), C_WIDTH, C_HEIGHT );
				GlobSetCanvasMargins( c_xnE[i] );
				h_xnE[i]->Draw();
				
				
				// (5) *CANVAS* XF-E hist monochrome
				c_xfE[i] = new TCanvas( Form( "c_xfE_%i", i ), Form( "XF-E | Det %i", i ), C_WIDTH, C_HEIGHT );
				GlobSetCanvasMargins( c_xfE[i] );
				h_xfE[i]->Draw();
				
				
				// (6) *CANVAS* XN-XF hist coloured
				c_xnxf_colour[i] = new TCanvas( Form( "c_xnxf_col_%i", i ), Form( "XN-XF Colour Plot | Det %i", i ), C_WIDTH, C_HEIGHT );
				GlobSetCanvasMargins( c_xnxf_colour[i] );
				h_xnxf_colour[i][0]->Draw();
				h_xnxf_colour[i][1]->Draw("SAME");
				h_xnxf_colour[i][2]->Draw("SAME");
				h_xnxf_colour[i][3]->Draw("SAME");
				
				for ( Int_t j = 0; j < 4; j++ ){
					bline_xnxf[i][j]->Draw();
				}
				if ( xnxf_cut != NULL ){ xnxf_cut->Draw("SAME"); }
				
				// (7) *CANVAS* XN-E hist coloured
				c_xnE_colour[i] = new TCanvas( Form( "c_xnE_col_%i", i ), Form( "XN-E Colour Plot | Det %i", i ), C_WIDTH, C_HEIGHT );
				GlobSetCanvasMargins( c_xnE_colour[i] );
				h_xnE_colour[i][0]->Draw();
				h_xnE_colour[i][1]->Draw("SAME");
				
				// (8) *CANVAS* XF-E hist coloured
				c_xfE_colour[i] = new TCanvas( Form( "c_xfE_col_%i", i ), Form( "XF-E Colour Plot | Det %i", i ), C_WIDTH, C_HEIGHT );
				GlobSetCanvasMargins( c_xfE_colour[i] );
				h_xfE_colour[i][0]->Draw();
				h_xfE_colour[i][1]->Draw("SAME");
				
				// (9) XNXF-E hist coloured
				c_xnxfE_colour[i] = new TCanvas( Form( "c_xnxfE_col_%i", i ), Form( "XNXF-E Colour Plot | Det %i", i ), C_WIDTH, C_HEIGHT );
				GlobSetCanvasMargins( c_xnxfE_colour[i] );
				h_xnxfE_colour[i][0]->Draw();
				h_xnxfE_colour[i][1]->Draw("SAME");
				h_xnxfE_colour[i][2]->Draw("SAME");
				
				// PRINT SPECTRA
				if ( SW_XNXF[1] == 1 ){
					if( xnxf_print_opt[0] == 1 )PrintAll( c_xnxf[i],         spec_name[0] );		// (0) XN-XF hist monochrome
					if( xnxf_print_opt[1] == 1 )PrintAll( c_xnxf_p[i],       spec_name[1] );		// (1) XN-XF profile
					if( xnxf_print_opt[2] == 1 )PrintAll( c_xnxfE[i],        spec_name[2] );		// (2) XNXF-E hist monochrome
					if( xnxf_print_opt[3] == 1 )PrintAll( c_xnxfE_p[i],      spec_name[3] );		// (3) XNXF-E profile
					if( xnxf_print_opt[4] == 1 )PrintAll( c_xnE[i],          spec_name[4] );		// (4) XN-E hist monochrome
					if( xnxf_print_opt[5] == 1 )PrintAll( c_xfE[i],          spec_name[5] );		// (5) XF-E hist monochrome
					if( xnxf_print_opt[6] == 1 )PrintAll( c_xnxf_colour[i],  spec_name[6] );		// (6) XN-XF hist coloured
					if( xnxf_print_opt[7] == 1 )PrintAll( c_xnE_colour[i],   spec_name[7] );		// (7) XN-E hist coloured
					if( xnxf_print_opt[8] == 1 )PrintAll( c_xfE_colour[i],   spec_name[8] );		// (8) XF-E hist coloured
					if( xnxf_print_opt[9] == 1 )PrintAll( c_xnxfE_colour[i], spec_name[9] );		// (9) XNXF-E hist coloured
				}
				
			}
			else{
				// ***************************************************************************** //
				// TODO - Profile combine canvas
				/*
				c_xnxf_comb->cd(i+1);
				c_xnxf_comb->SetTitle( Form( "XN:XF #%i", i ) );
				h_xnxf[i]->Draw();
				TPad* pad = (TPad*)c_xnxf_comb->GetPad(i+1);
				//SetCanvasTitleFont( pad );
				
				for ( Int_t j = 0; j < 4; j++ ){
					bline_xnxf[i][j]->Draw();
				}
				
				
				
				// Print spectrum if desired
				if ( SW_EX[1] == 1 && i == 23 ){
					PrintAll( c_xnxf_comb, Form( "%s/posXXX_xnxf_comb", print_dir.Data() ) );
				}
				*/
				// ***************************************************************************** //
			}
			
			// Write ROOT file if desired
			if ( PRINT_ROOT == 1 && SW_XNXF[1] == 1 ){ f->cd(); h_xnxf[i]->Write(); }
			
			// Write SPE file if desired
			if ( SW_XNXF[2] == 1 ){ ErrorSPE("XN:XF"); }
		} 
	
	
	}
	PrintXNCorr();
	PrintXNXFECorr();
	// Close the root file
	if ( SW_XNXF[1] == 1 && PRINT_ROOT == 1 && f != NULL ){ if ( f->IsOpen() ){ f->Close(); } }
	return;
}

#endif