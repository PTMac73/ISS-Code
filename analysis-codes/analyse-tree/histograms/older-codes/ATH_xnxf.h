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
const Int_t NUM_DIFF_HISTS = 11;
const Int_t NUM_DETS = 24;
const Double_t XNXF_FRAC = 0.1;
TH2F* h_xnxf[NUM_DETS];					// (0) XN-XF hist monochrome
TProfile* p_xnxf[NUM_DETS];				// (1) XN-XF profile
TH2F* h_xnxfE[NUM_DETS];				// (2) XNXF-E hist monochrome
TProfile* p_xnxfE[NUM_DETS];			// (3) XNXF-E profile
TH2F* h_xnE[NUM_DETS];					// (4) XN-E hist monochrome
TH2F* h_xfE[NUM_DETS];					// (5) XF-E hist monochrome
TH2F* h_xnxf_colour[NUM_DETS][5];		// (6) XN-XF hist coloured
TH2F* h_xnE_colour[NUM_DETS][5];		// (7) XN-E hist coloured
TH2F* h_xfE_colour[NUM_DETS][5];		// (8) XF-E hist coloured
TH2F* h_xnxfE_colour[NUM_DETS][4];		// (9) XNXF-E hist coloured
TH1F* h_ecalibration[NUM_DETS][2];		//(10) E Calibration hist

// Define which ones to print
Bool_t xnxf_print_opt[NUM_DIFF_HISTS] = {
	0,	// (0) XN-XF hist monochrome (OBSOLETE)
	0,	// (1) XN-XF profile
	0,	// (2) XNXF-E hist monochrome
	0,	// (3) XNXF-E profile
	0,	// (4) XN-E hist monochrome (OBSOLETE)
	0,	// (5) XF-E hist monochrome (OBSOLETE)
	0,	// (6) XN-XF hist coloured
	0,	// (7) XN-E hist coloured
	0,	// (8) XF-E hist coloured
	0,	// (9) XNXF-E hist coloured
	1	//(10) E Calibration hist
};

// -------------------------------- { 00, 01, 02, 03, 04, 05, 06, 07, 08, 09, 10 };
TString spec_name[NUM_DIFF_HISTS] = { "", "", "", "", "", "", "", "", "", "", "" };

// Additional variables
TString xnxfcut_dir = "/home/ptmac/Documents/07-CERN-ISS-Mg/Mg-Analysis/SPE-Files/xnxf-calibration/pos1_xnxf_cuts.root";
TCutG* xnxf_cut;

// Store the correction parameters
Double_t xn_corr[NUM_DETS];
Double_t xfxne_corr[NUM_DETS][2];
Double_t e_corr[NUM_DETS][2];


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


// Populate TString array of names
void MakeSpecName( Int_t i ){
	spec_name[0]  = Form( "%s/posXXX_xnxf_%i", print_dir.Data(), i );		// (0) XN-XF hist monochrome
	spec_name[1]  = Form( "%s/posXXX_xnxf_p_%i", print_dir.Data(), i );		// (1) XN-XF profile
	spec_name[2]  = Form( "%s/posXXX_xnxfE_%i", print_dir.Data(), i );		// (2) XNXF-E hist monochrome
	spec_name[3]  = Form( "%s/posXXX_xnxfE_p_%i", print_dir.Data(), i );	// (3) XNXF-E profile
	spec_name[4]  = Form( "%s/posXXX_xnE_%i", print_dir.Data(), i );		// (4) XN-E hist monochrome
	spec_name[5]  = Form( "%s/posXXX_xfE_%i", print_dir.Data(), i );		// (5) XF-E hist monochrome
	spec_name[6]  = Form( "%s/posXXX_xnxf_col_%i", print_dir.Data(), i );	// (6) XN-XF hist coloured
	spec_name[7]  = Form( "%s/posXXX_xnE_col_%i", print_dir.Data(), i );	// (7) XN-E hist coloured
	spec_name[8]  = Form( "%s/posXXX_xfE_col_%i", print_dir.Data(), i );	// (8) XF-E hist coloured
	spec_name[9]  = Form( "%s/posXXX_xnxfE_col_%i", print_dir.Data(), i );	// (9) XNXF-E hist coloured
	spec_name[10] = Form( "%s/posXXX_ecalibration_%i", print_dir.Data(), i );	//(10) E Calibration hist
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
		for ( Int_t j = 0; j < 5; j++ ){ h_xnxf_colour[i][j] = NULL; }	// (6) XN-XF hist coloured
		for ( Int_t j = 0; j < 5; j++ ){ h_xnE_colour[i][j] = NULL; } 	// (7) XN-E hist coloured
		for ( Int_t j = 0; j < 5; j++ ){ h_xfE_colour[i][j] = NULL; }	// (8) XF-E hist coloured
		for ( Int_t j = 0; j < 4; j++ ){ h_xnxfE_colour[i][j] = NULL; }	// (9) XNXF-E hist coloured
		for ( Int_t j = 0; j < 2; j++ ){ h_ecalibration[i][j] = NULL; }	//(10) E Calibration hist
		
		// Detector by detector
		if ( i == DET_NUMBER || ( DET_NUMBER == -1 ) ){
		
			// (0) *HIST* XN-XF hist monochrome
			h_xnxf[i] = new TH2F( Form( "h_xnxf_%i", i ), Form( "XN-XF SPECTRUM %i", i ), 200, 0, 2000, 200, 0, 2000 );
			GlobCreate2DHists( h_xnxf[i], "XF", "XN" );
			
			// (1) *HIST* XN-XF profile
			p_xnxf[i] = new TProfile( Form( "p_xnxf_%i", i ), Form( "XN-XF SPECTRUM %i", i ), 200, xnxf_lims[i][0], xnxf_lims[i][1], xnxf_lims[i][2], xnxf_lims[i][3] );
			GlobCreateProfile( p_xnxf[i], "XF", "XN" );
			
			// (2) *HIST* XNXF-E hist monochrome
			h_xnxfE[i] = new TH2F( Form( "h_xnxfE_%i", i ), Form( "XNXF-E SPECTRUM %i", i ), 200, 0, 2000, 200, 0, 2000 );
			GlobCreate2DHists( h_xnxfE[i], "XNXF", "E" );
			
			// (3) *HIST* XNXF-E profile
			p_xnxfE[i] = new TProfile( Form( "p_xnxfE_%i", i ), Form( "XNXF-E SPECTRUM %i", i ), 200, xnxfE_lims[i][0], xnxfE_lims[i][1], xnxfE_lims[i][2], xnxfE_lims[i][3] );
			GlobCreateProfile( p_xnxfE[i], "XNXF", "E" );
			
			// (4) *HIST* XN-E hist monochrome
			h_xnE[i] = new TH2F( Form( "h_xnE_%i", i ), Form( "XN-E SPECTRUM %i", i ), 200, 0, 2000, 200, 0, 2000 );
			GlobCreate2DHists( h_xnE[i], "XN", "E" );
			
			// (5) *HIST* XF-E hist monochrome
			h_xfE[i] = new TH2F( Form( "h_xfE_%i", i ), Form( "XF-E SPECTRUM %i", i ), 200, 0, 2000, 200, 0, 2000 );
			GlobCreate2DHists( h_xfE[i], "XF", "E" );
			
			// (6) *HIST* XN-XF hist coloured
			for ( Int_t j = 0; j < 5; j++ ){
				h_xnxf_colour[i][j] = new TH2F( Form( "h_xnxf_%i_%i", i, j ), Form( "XN-XF COLOUR SPECTRUM | Det %i | Case %i", i, j ), 200, 0, 2000, 200, 0, 2000 );
				GlobCreate2DHists( h_xnxf_colour[i][j], "XF", "XN" );
			}
			h_xnxf_colour[i][0]->SetMarkerColor( kBlue );
			h_xnxf_colour[i][1]->SetMarkerColor( kBlack );
			h_xnxf_colour[i][2]->SetMarkerColor( kRed );
			h_xnxf_colour[i][3]->SetMarkerColor( kGreen );
			h_xnxf_colour[i][4]->SetMarkerColor( kMagenta+2 );
			
			// (7) *HIST*  XN-E hist coloured
			for ( Int_t j = 0; j < 5; j++ ){
				h_xnE_colour[i][j] = new TH2F( Form( "h_xnE_%i_%i", i, j ), Form( "XN-E COLOUR SPECTRUM | Det %i | CASE %i", i, j ), 200, 0, 2000, 200, 0, 2000 );
				GlobCreate2DHists( h_xnE_colour[i][j], "XN", "E" );
			}
			h_xnE_colour[i][0]->SetMarkerColor( kBlue );
			h_xnE_colour[i][1]->SetMarkerColor( kBlack );
			h_xnE_colour[i][2]->SetMarkerColor( kRed );
			h_xnE_colour[i][3]->SetMarkerColor( kGreen );
			h_xnE_colour[i][4]->SetMarkerColor( kMagenta+2 );
			
			// (8) *HIST* XF-E hist coloured
			for ( Int_t j = 0; j < 5; j++ ){
				h_xfE_colour[i][j] = new TH2F( Form( "h_xfE_%i_%i", i, j ), Form( "XF-E COLOUR SPECTRUM | Det %i | CASE %i", i, j ), 200, 0, 2000, 200, 0, 2000 );
				GlobCreate2DHists( h_xfE_colour[i][j], "XF", "E" );
			}
			h_xfE_colour[i][0]->SetMarkerColor( kBlue );
			h_xfE_colour[i][1]->SetMarkerColor( kBlack );
			h_xfE_colour[i][2]->SetMarkerColor( kRed );
			h_xfE_colour[i][3]->SetMarkerColor( kGreen );
			h_xfE_colour[i][4]->SetMarkerColor( kMagenta+2 );
			
			// (9) XNXF-E hist coloured
			for ( Int_t j = 0; j < 4; j++ ){
				h_xnxfE_colour[i][j] = new TH2F( Form( "h_xnxfE_%i_%i", i, j ), Form( "XNXF-E COLOUR SPECTRUM | Det %i | CASE %i", i, j ), 200, 0, 2000, 200, 0, 2000 );
				GlobCreate2DHists( h_xnxfE_colour[i][j], "XNXF", "E" );
			}
			h_xnxfE_colour[i][0]->SetMarkerColor( kBlue );
			h_xnxfE_colour[i][1]->SetMarkerColor( kGreen );
			h_xnxfE_colour[i][2]->SetMarkerColor( kRed );
			h_xnxfE_colour[i][3]->SetMarkerColor( kMagenta+2 );
			
			//(10) *HIST* E Calibration
			h_ecalibration[i][0] = new TH1F( Form( "h_ecalibration_%i", i ), Form( "E CALIBRATION SPECTRUM | DET %i", i ), 2000, 0, 2000 );
			h_ecalibration[i][1] = new TH1F( Form( "h_ecalibration_cuts%i", i ), Form( "E CALIBRATION SPECTRUM CUTS | DET %i", i ), 2000, 0, 2000 );
			
			h_ecalibration[i][0]->SetLineColorAlpha( kBlue, 0.5 );
			h_ecalibration[i][0]->SetTitle("");
			h_ecalibration[i][0]->GetXaxis()->SetTitle( "Raw E" );
			h_ecalibration[i][0]->GetYaxis()->SetTitle( "Counts" );
			GlobSetHistFonts( h_ecalibration[i][0] );
			
			h_ecalibration[i][1]->SetLineColor( kRed );
			h_ecalibration[i][1]->SetTitle("");
			h_ecalibration[i][1]->GetXaxis()->SetTitle( "" );
			h_ecalibration[i][1]->GetYaxis()->SetTitle( "" );
			GlobSetHistFonts( h_ecalibration[i][1] );
			

			// Correction parameters		
			xn_corr[i] = 0.0;
			xfxne_corr[i][0] = 0.0;
			xfxne_corr[i][1] = 0.0;
			e_corr[i][0] = 0.0;
			e_corr[i][1] = 0.0;
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
	
	TCanvas* c_ecalibration[24];	// (10) E calibration
	TCanvas* c_ecalibration_comb;
	
	
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
	fitline_xnxf->SetLineColorAlpha(kBlack,0.5);
	fitline_xnxf->SetParameters( 0, 0 );
	fitline_xnxf->SetRange( 200, 1800 );
	
	TF1* fitline_xnxfE = new TF1( "fitline_xnxfE", "[0] + [1]*x", 0, 100 );
	fitline_xnxfE->SetLineWidth(2);
	fitline_xnxfE->SetLineColorAlpha(kBlack,0.5);
	fitline_xnxfE->SetParameters( 0, 0 );
	fitline_xnxfE->SetRange( 200, 1800 );
	
	// Define peak markers and fitting functions
	TLine* peak_mark[4];
	TF1* fitline_ecalibration_gaus_fit[5];
	
	
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
			for ( Int_t j = 0; j < 4; j++ ){
				peak_mark[j] = new TLine( rawE_pos[i][j], 0.0, rawE_pos[i][j],1000.0 );
				peak_mark[j]->SetLineColorAlpha( kBlue, 0.5 );
			}
			
			// TPROFILE FIT LINES
			// Fit some lines to the XN:XF TProfile plots
			TFitResultPtr fit_ptr = p_xnxf[i]->Fit("pol1", "SQ0", "", xnxf_lims[i][0], xnxf_lims[i][1] );
			if ( fit_ptr.Get() != NULL ){
				fitline_xnxf->SetRange( xnxf_lims[i][0], xnxf_lims[i][1] );
				fitline_xnxf->SetParameters( fit_ptr->Parameter(0), fit_ptr->Parameter(1) );
				xn_corr[i] = -1.0/fit_ptr->Parameter(1);
				//PrintXNXFFitPars( fit_ptr );
			}
			
			// Fit some lines to the XNXF:E TProfile plots
			TFitResultPtr fit_ptrE = p_xnxfE[i]->Fit("pol1", "SQ0", "", xnxfE_lims[i][0], xnxfE_lims[i][1] );
			if ( fit_ptrE.Get() != NULL ){
				fitline_xnxfE->SetRange( xnxfE_lims[i][0], xnxfE_lims[i][1] );
				fitline_xnxfE->SetParameters( fit_ptrE->Parameter(0), fit_ptrE->Parameter(1) );
				xfxne_corr[i][0] = 0.5*fit_ptrE->Parameter(0); 	// Intercept
				xfxne_corr[i][1] = fit_ptrE->Parameter(1); 		// Gradient
				//PrintXNXFEFitPars( fit_ptrE );
			}
			// PEAK FITTING FOR ENERGY CALIBRATION (0-3 holds peak info, 4 is the combined triplet)
			for ( Int_t j = 0; j < 4; j++ ){
				fitline_ecalibration_gaus_fit[j] = new TF1( Form( "fitline_ecalibration_%i", j ), "[0] + [1]*x + [2]*exp(-0.5*((x-[3])/[4])^2)", rawE_pos[i][j] - e_calibration_range, rawE_pos[i][j] + e_calibration_range );
			}
			fitline_ecalibration_gaus_fit[4] = new TF1( "fitline_ecalibration_4", "[0] + [1]*x + [2]*exp(-0.5*((x-[3])/[8])^2) + [4]*exp(-0.5*((x-[5])/[8])^2) + [6]*exp(-0.5*((x-[7])/[8])^2)", rawE_pos[i][1] - e_calibration_range, rawE_pos[i][3] + e_calibration_range );
		
			
			// Set parameter limits
			fitline_ecalibration_gaus_fit[0]->SetParLimits( 0, 0, 100 );		// Constant BG
			fitline_ecalibration_gaus_fit[0]->SetParLimits( 1, -100, 100 );		// Gradient BG
			fitline_ecalibration_gaus_fit[0]->SetParLimits( 2, 0, 2000 );		// Amplitude
			fitline_ecalibration_gaus_fit[0]->SetParLimits( 3, rawE_pos[i][0] - e_calibration_range, rawE_pos[i][0] + e_calibration_range );	// Mu
			fitline_ecalibration_gaus_fit[0]->SetParLimits( 4, 0, 100 );		// Sigma
			
			fitline_ecalibration_gaus_fit[4]->SetParLimits( 0, 0, 100 );		// Constant BG
			fitline_ecalibration_gaus_fit[4]->SetParLimits( 1, -2, 2 );			// Gradient BG
			fitline_ecalibration_gaus_fit[4]->SetParLimits( 2, 0, 2000 );		// Amplitude
			fitline_ecalibration_gaus_fit[4]->SetParLimits( 3, rawE_pos[i][1] - e_calibration_range, rawE_pos[i][1] + e_calibration_range );	// Mu
			fitline_ecalibration_gaus_fit[4]->SetParLimits( 4, 0, 2000 );		// Amplitude
			fitline_ecalibration_gaus_fit[4]->SetParLimits( 5, rawE_pos[i][2] - e_calibration_range, rawE_pos[i][2] + e_calibration_range );	// Mu
			fitline_ecalibration_gaus_fit[4]->SetParLimits( 6, 0, 2000 );		// Amplitude
			fitline_ecalibration_gaus_fit[4]->SetParLimits( 7, rawE_pos[i][3] - e_calibration_range, rawE_pos[i][3] + e_calibration_range );	// Mu
			fitline_ecalibration_gaus_fit[4]->SetParLimits( 8, 0, 50 );			// Sigma
			
			
			// Set formatting
			for ( Int_t j = 0; j < 5; j++ ){
				fitline_ecalibration_gaus_fit[j]->SetLineWidth(1);
				fitline_ecalibration_gaus_fit[j]->SetLineColor( ( j == 0 || j == 4 ? kBlack : kBlue ) );
			}
			
			// Set initial values
			fitline_ecalibration_gaus_fit[0]->SetParameters( 0.0, 0.0, 500, rawE_pos[i][0], 20 );
			fitline_ecalibration_gaus_fit[4]->SetParameter( 0, 0.0 );					// Const. BG
			fitline_ecalibration_gaus_fit[4]->SetParameter( 1, 0.0 );					// Gradient BG
			fitline_ecalibration_gaus_fit[4]->SetParameter( 2, 300 );					// Amp 1
			fitline_ecalibration_gaus_fit[4]->SetParameter( 3, rawE_pos[i][1] );		// Mu 1
			fitline_ecalibration_gaus_fit[4]->SetParameter( 4, 300 );					// Amp 2
			fitline_ecalibration_gaus_fit[4]->SetParameter( 5, rawE_pos[i][2] );		// Mu 2
			fitline_ecalibration_gaus_fit[4]->SetParameter( 6, 300 );					// Amp 3
			fitline_ecalibration_gaus_fit[4]->SetParameter( 7, rawE_pos[i][3] );		// Mu 3
			fitline_ecalibration_gaus_fit[4]->SetParameter( 8, 40 );					// Sigma 1
			
			// Set ranges
			/*
			for ( Int_t j = 0; j < 4; j++ ){
				fitline_ecalibration_gaus_fit[j]->SetRange( rawE_pos[i][j] - 3*e_calibration_range, rawE_pos[i][j] + 3*e_calibration_range );
			}
			fitline_ecalibration_gaus_fit[4]->SetRange( rawE_pos[i][1] - 3*e_calibration_range, rawE_pos[i][3] + 3*e_calibration_range );
			*/
			
			h_ecalibration[i][1]->Fit( "fitline_ecalibration_0", "0" );
			h_ecalibration[i][1]->Fit( "fitline_ecalibration_4", "0" );
			
			// Set the fit parameters for peaks 1--3 here
			for ( Int_t j = 1; j < 4; j++ ){
				fitline_ecalibration_gaus_fit[j]->SetParameter( 0, fitline_ecalibration_gaus_fit[4]->GetParameter(0) );
				fitline_ecalibration_gaus_fit[j]->SetParameter( 1, fitline_ecalibration_gaus_fit[4]->GetParameter(1) );
				fitline_ecalibration_gaus_fit[j]->SetParameter( 2, fitline_ecalibration_gaus_fit[4]->GetParameter(2*j) );
				fitline_ecalibration_gaus_fit[j]->SetParameter( 3, fitline_ecalibration_gaus_fit[4]->GetParameter(2*j+1) );
				fitline_ecalibration_gaus_fit[j]->SetParameter( 4, fitline_ecalibration_gaus_fit[4]->GetParameter(8) );
			}
			
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
				p_xnxf[i]->Draw("E1");
				fitline_xnxf->Draw("SAME");
				
				
				// (2) *CANVAS* XNXF-E hist monochrome
				c_xnxfE[i] = new TCanvas( Form( "c_xnxfE_%i", i ), Form( "XNXF-E | Det %i", i ), C_WIDTH, C_HEIGHT );
				GlobSetCanvasMargins( c_xnxfE[i] );
				h_xnxfE[i]->Draw();
				/*for ( Int_t j = 0; j < 4; j++ ){
					bline_xnxfE[i][j]->Draw();
				}*/
				fitline_xnxfE->Draw("SAME");
				
				// (3) *CANVAS* XNXF-E profile
				c_xnxfE_p[i] = new TCanvas( Form( "c_xnxf_pE_%i", i ), Form( "XNXF-E Profile | Det %i", i ), C_WIDTH, C_HEIGHT );
				GlobSetCanvasMargins( c_xnxfE_p[i] );
				p_xnxfE[i]->Draw("E1");
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
				h_xnxf_colour[i][4]->Draw("SAME");
				
				for ( Int_t j = 0; j < 4; j++ ){
					bline_xnxf[i][j]->Draw();
				}
				if ( xnxf_cut != NULL ){ xnxf_cut->Draw("SAME"); }
				
				// (7) *CANVAS* XN-E hist coloured
				c_xnE_colour[i] = new TCanvas( Form( "c_xnE_col_%i", i ), Form( "XN-E Colour Plot | Det %i", i ), C_WIDTH, C_HEIGHT );
				GlobSetCanvasMargins( c_xnE_colour[i] );
				h_xnE_colour[i][0]->Draw();
				h_xnE_colour[i][1]->Draw("SAME");
				h_xnE_colour[i][2]->Draw("SAME");
				h_xnE_colour[i][3]->Draw("SAME");
				h_xnE_colour[i][4]->Draw("SAME");
				
				// (8) *CANVAS* XF-E hist coloured
				c_xfE_colour[i] = new TCanvas( Form( "c_xfE_col_%i", i ), Form( "XF-E Colour Plot | Det %i", i ), C_WIDTH, C_HEIGHT );
				GlobSetCanvasMargins( c_xfE_colour[i] );
				h_xfE_colour[i][0]->Draw();
				h_xfE_colour[i][1]->Draw("SAME");
				h_xfE_colour[i][2]->Draw("SAME");
				h_xfE_colour[i][3]->Draw("SAME");
				h_xfE_colour[i][4]->Draw("SAME");
				
				// (9) XNXF-E hist coloured
				c_xnxfE_colour[i] = new TCanvas( Form( "c_xnxfE_col_%i", i ), Form( "XNXF-E Colour Plot | Det %i", i ), C_WIDTH, C_HEIGHT );
				GlobSetCanvasMargins( c_xnxfE_colour[i] );
				h_xnxfE_colour[i][3]->Draw();
				h_xnxfE_colour[i][0]->Draw("SAME");
				h_xnxfE_colour[i][1]->Draw("SAME");
				h_xnxfE_colour[i][2]->Draw("SAME");
				
				fitline_xnxfE->Draw("SAME");
				
				// (10) E Calibration
				c_ecalibration[i] = new TCanvas( Form( "c_ecalibration_%i", i ), Form( "Energy Calibration | Det %i", i ), C_WIDTH, C_HEIGHT );
				h_ecalibration[i][0]->Draw();
				h_ecalibration[i][1]->Draw("SAME");
				/*fitline_ecalibration_gaus_fit[0]->Draw("SAME");
				fitline_ecalibration_gaus_fit[4]->Draw("SAME");
				fitline_ecalibration_gaus_fit[1]->Draw("SAME");
				fitline_ecalibration_gaus_fit[2]->Draw("SAME");
				fitline_ecalibration_gaus_fit[3]->Draw("SAME");*/

				for ( Int_t j = 0; j < 4; j++ ){
					peak_mark[j]->Draw("SAME");
				}
				
				// PRINT SPECTRA
				if ( SW_XNXF[1] == 1 ){
					if( xnxf_print_opt[0] ==  1 )PrintAll( c_xnxf[i],         spec_name[0]  );		// (0) XN-XF hist monochrome
					if( xnxf_print_opt[1] ==  1 )PrintAll( c_xnxf_p[i],       spec_name[1]  );		// (1) XN-XF profile
					if( xnxf_print_opt[2] ==  1 )PrintAll( c_xnxfE[i],        spec_name[2]  );		// (2) XNXF-E hist monochrome
					if( xnxf_print_opt[3] ==  1 )PrintAll( c_xnxfE_p[i],      spec_name[3]  );		// (3) XNXF-E profile
					if( xnxf_print_opt[4] ==  1 )PrintAll( c_xnE[i],          spec_name[4]  );		// (4) XN-E hist monochrome
					if( xnxf_print_opt[5] ==  1 )PrintAll( c_xfE[i],          spec_name[5]  );		// (5) XF-E hist monochrome
					if( xnxf_print_opt[6] ==  1 )PrintAll( c_xnxf_colour[i],  spec_name[6]  );		// (6) XN-XF hist coloured
					if( xnxf_print_opt[7] ==  1 )PrintAll( c_xnE_colour[i],   spec_name[7]  );		// (7) XN-E hist coloured
					if( xnxf_print_opt[8] ==  1 )PrintAll( c_xfE_colour[i],   spec_name[8]  );		// (8) XF-E hist coloured
					if( xnxf_print_opt[9] ==  1 )PrintAll( c_xnxfE_colour[i], spec_name[9]  );		// (9) XNXF-E hist coloured
					if( xnxf_print_opt[10] == 1 )PrintAll( c_ecalibration[i], spec_name[10] );		//(10) E Calibration
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
			if ( PRINT_ROOT == 1 && SW_XNXF[1] == 1 ){
				f->cd();
				if( xnxf_print_opt[0] ==  1 ){ h_xnxf[i]->Write(); }			// (0) XN-XF hist monochrome
				if( xnxf_print_opt[1] ==  1 ){ p_xnxf[i]->Write(); }			// (1) XN-XF profil
				if( xnxf_print_opt[2] ==  1 ){ h_xnxfE[i]->Write(); }			// (2) XNXF-E hist monochrome
				if( xnxf_print_opt[3] ==  1 ){ p_xnxfE[i]->Write(); }			// (3) XNXF-E profile
				if( xnxf_print_opt[4] ==  1 ){ h_xnE[i]->Write(); }				// (4) XN-E hist monochrome
				if( xnxf_print_opt[5] ==  1 ){ h_xfE[i]->Write(); }				// (5) XF-E hist monochrome
				if( xnxf_print_opt[6] ==  1 ){ for( Int_t j = 0; j < 5; j++ ){ h_xnxf_colour[i][j]->Write(); } }		// (6) XN-XF hist coloured
				if( xnxf_print_opt[7] ==  1 ){ for( Int_t j = 0; j < 5; j++ ){ h_xnE_colour[i][j]->Write(); } }			// (7) XN-E hist coloured
				if( xnxf_print_opt[8] ==  1 ){ for( Int_t j = 0; j < 5; j++ ){ h_xfE_colour[i][j]->Write(); } }			// (8) XF-E hist coloured
				if( xnxf_print_opt[9] ==  1 ){ for( Int_t j = 0; j < 4; j++ ){ h_xnxfE_colour[i][j]->Write(); } }		// (9) XNXF-E hist coloured
				if( xnxf_print_opt[10] == 1 ){ for( Int_t j = 0; j < 2; j++ ){ h_ecalibration[i][j]->Write(); } }		//(10) E Calibration
			}
			
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
