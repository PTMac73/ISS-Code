// ExtractYields.C
// Extracts the yields for a given histogram file provided when creating the Mg spectra
// ============================================================================================= //
// Patrick MacGregor
// Nuclear Physics Research Group
// Department of Physics and Astronomy
// The University of Manchester
// ============================================================================================= //

#include "ExtractYields.h"
#include "EY_GenerateFitName.h"
#include "EY_HistogramFunctions.h"
#include "EY_PeakFunctions.h"
#include "EY_FitPeaks.h"

#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TH1.h>
#include <TLatex.h>
#include <TLine.h>
#include <TPad.h>
#include <TROOT.h>
#include <TString.h>
#include <TStyle.h>
#include <TText.h>

#include <iostream>
#include <fstream>

struct FitPeakOptions_t;

void ExtractYieldsHist( TH1D* h, Int_t row_num, Int_t pos, Int_t spectrum_type ){
	// DECLARE SOME VARIABLES ------------------------------------------------------------------ //
	TCanvas* c_spec;			// Canvas pointers
	TCanvas* c_pre_fit;
	Double_t sig_est = 0;
	Double_t sig_range = 0;
	TFitResultPtr fit_result_ptr;

	// Format the histogram
	if ( SW_HIST_TITLE == 1 ){ h->SetTitle( GenerateHistTitle( pos, row_num, spectrum_type ) ); }
	else{ h->SetTitle(""); }
	FormatHistogram( h );

	
	// FIT THE TWO STRONG PEAKS IF POSSIBLE ------------------------------------------------ //
	if ( num_peaks_in_spectrum[row_num][pos-1] > 1 ){
	
		// Estimate the peak parameters
		FitPeakOptions_t pre_fit_opt;
		SetFitOptions( pre_fit_opt, 2, 2 );
		Int_t** pre_var_type_arr;
		TF1* pre_fit_func = EstimatePeakParameters( h, pre_fit_opt, pre_var_type_arr );
		
		// Fit the spectrum
		h->Fit( pre_fit_func, "RB" );
		
		// Draw if desired
		if ( SW_DRAW_PRE_FIT == 1 ){
			c_pre_fit = new TCanvas( "c_pre_fit", "CANVAS", CANVAS_WIDTH, CANVAS_HEIGHT );
			SetPadMargins( (TPad*)c_pre_fit->GetPad(0) );
			h->Draw();
			c_pre_fit->Modified(); c_pre_fit->Update();
			c_pre_fit->Print( GenerateCanvasPrintName( pos, row_num, spectrum_type, 1 ) );
		}
		
		// Calculate the width
		sig_est = pre_fit_func->GetParameter( pre_var_type_arr[1][2] );
		sig_range = pre_fit_func->GetParError( pre_var_type_arr[1][2] );
	}
	else{
		// Estimate the widths from this
		sig_est = 0.25*PEAK_WIDTH_ESTIMATE;
		sig_range = 0.2*PEAK_WIDTH_ESTIMATE;
	}
	std::cout << "SIG-EST  : " << std::setprecision(8) << sig_est << "\n";
	std::cout << "SIG-RANGE: " << std::setprecision(8) << sig_range << "\n";

	// DRAW A CANVAS AND DRAW THE INITIAL SPECTRUM ----------------------------------------- //
	// Define a canvas and set margins

	c_spec = new TCanvas( GenerateCanvasName( row_num ), "CANVAS", CANVAS_WIDTH, CANVAS_HEIGHT );
	SetPadMargins( (TPad*)c_spec->GetPad(0) );

	// Draw the canvas
	h->Draw();

	// Remove the stats box
	gStyle->SetOptStat(0);
	c_spec->Modified(); c_spec->Update();
	

	// START FULL FITTING PROCESS HERE ----------------------------------------------------- //
	// Estimate the peak parameters
	Int_t num_full_peaks = num_peaks_in_spectrum[row_num][pos-1];
	FitPeakOptions_t fit_opt;
	SetFitOptions( fit_opt, num_full_peaks, 0, sig_est, sig_range );
	Int_t** var_type_arr;
	TF1* fit_func = EstimatePeakParameters( h, fit_opt, var_type_arr );
	
	// Format the fit function
	fit_func->SetLineColor( kBlack );
	fit_func->SetLineWidth(3);
	
	// Fit the histogram
	fit_result_ptr = h->Fit( fit_func, "S" );
	PrintVTA( var_type_arr, num_full_peaks );

	// LOOP OVER THE PEAKS TO DO USEFUL PEAK-SPECIFIC THINGS ------------------------------- //
	// Variables
	Peak_t peaks[num_full_peaks];
	RedPeak_t gsdoub;
	TF1* f_peaks[num_full_peaks];
		
	// Delete the old data file if it exists
	if ( std::remove( GetDataFileName( pos, spectrum_type, row_num ) ) != 0 ){
		std::cout << "Error deleting data file.\n";
	}

	// Open the file
	std::ofstream out_file;
	out_file.open( GetDataFileName( pos, spectrum_type, row_num ), std::ofstream::app );
	
	
	// Start the loop
	for ( Int_t i = 0; i < num_full_peaks; i++ ){

		// Calculate just the yield of the ground state doublet
		/*if ( i == 0 ){
			CalculateRedPeakQuantities( fit_result_ptr, BG_DIM, GS_DOUB_LB, GS_DOUB_UB, h, gsdoub );
			PrintRedPeak( gsdoub, out_file );
		}*/
		
		// Calculate peak quantities and print the peaks to file ( including the GS doublet for comparison)
		CalculatePeakQuantities( fit_result_ptr, BG_DIM, var_type_arr[i+1][0], var_type_arr[i+1][2], peaks[i] );
		PrintPeak( peaks[i], out_file );
			
		
		// Populate TF1's to hold individual peak info
		/*if ( i == 0 ){
			f_peaks[i] = new TF1( Form( "fit_func%i", i ), GetPolynomialString( 0, BG_DIM ), gsdoub.lb, gsdoub.ub );
		}
		else{*/
			f_peaks[i] = new TF1( Form( "fit_func%i", i ), GetFitStringSimple( 1, BG_DIM ), peak_energies[i] - 1.5*PEAK_WIDTH_ESTIMATE, peak_energies[i] + 1.5*PEAK_WIDTH_ESTIMATE );
		//}
		f_peaks[i]->SetNpx(2000);
		
		// Set individual backgrounds
		for ( Int_t j = 0; j < BG_DIM + 1; j++ ){
			f_peaks[i]->FixParameter( j, fit_func->GetParameter( var_type_arr[0][j] ) );
		}
		
		// Set individual peak shapes
		//if ( i > 0 ){
			for ( Int_t j = 0; j < 3; j++ ){
				f_peaks[i]->FixParameter( j + BG_DIM + 1, fit_func->GetParameter(var_type_arr[i+1][j]) );
			}
		//}
		
		// Format individual peaks (based on whether width is fixed or not)
		if ( peak_fix_widths[i] == 1 || i == 0 ){ f_peaks[i]->SetLineColor(kRed); }		// Fixed width peaks
		else{ f_peaks[i]->SetLineColor(kBlue); }
												// Unbound states
		f_peaks[i]->SetLineWidth(2);
		f_peaks[i]->Draw("SAME");
		
		
	} // LOOP OVER PEAKS (i)
	
	// Add a line for the neutron separation energy
	std::cout << "MAX: " << h->GetYaxis()->GetXmax() << "\n";
	TLine* lNSE = new TLine( NEUTRON_SEP_ENERGY, 0, NEUTRON_SEP_ENERGY, c_spec->GetUymax() );
	lNSE->SetLineColor(kBlack);
	lNSE->SetLineStyle(2);
	lNSE->Draw("SAME");
	
	TLatex* tNSE = new TLatex( 0.95*NEUTRON_SEP_ENERGY, 0.9*h->GetMaximum(), Form( "S_{n} = %4.3f MeV", NEUTRON_SEP_ENERGY ) );
	tNSE->SetTextAlign(32);
	tNSE->Draw("SAME");
	
	// Print the last line
	PrintFitPars( fit_result_ptr, BG_DIM, out_file );
		
	// Append new line and close the file
	out_file.close();
	
	// Print the spectrum
	c_spec->Print( GenerateCanvasPrintName( pos, row_num, spectrum_type ) );
	
	
	
	
	
	
	/*
	// Delete arrays that are no longer needed and free them up for the next loop around
	for ( Int_t j = 0; j < 3; j++ ){
		delete[] pre_var_type_arr[j];
	}
	delete[] pre_var_type_arr;
	
	
	for ( Int_t j = 0; j < num_full_peaks+1; j++ ){
		delete[] var_type_arr[j];
		if ( j < num_full_peaks ){
			if ( f_peaks[j]->IsOnHeap() ){ f_peaks[j]->Delete(); }
		}
	}
	delete[] var_type_arr;
		
	
	// Delete other pointers
	if ( SW_BATCH_MODE == 1 ){
		if ( pre_fit_func->IsOnHeap() ){ pre_fit_func->Delete(); }
		if ( fit_func->IsOnHeap() ){ fit_func->Delete(); }
		if ( c_pre_fit->IsOnHeap() ){ delete c_pre_fit; }
		if ( c_spec->IsOnHeap() ){ delete c_spec; }
	}
*/
	return;
}
	
	
void ExtractYields( Int_t pos = 1, Int_t spectrum_type = 0 ){
	// Declare variables
	TFile* f;						// File containing the histograms
	
	// Batch mode
	if ( SW_BATCH_MODE == 1 ){ gROOT->SetBatch(kTRUE); }
	
	// Open the relevant TFile
	f = new TFile( GetROOTFileName( pos, spectrum_type ), "READ" );

	// Check that it opened
	if ( f->IsOpen() ){ std::cout << "File opened successfully!\n"; }
	else{
		std::cout << "File " << GetROOTFileName( pos, spectrum_type ) << " failed to open. Exiting process...\n";
		exit(1);
	}
	
	// Get a list of histograms in that TFile
	Int_t num_hists = GetNumHists(f);
	
	// Get the names in a TString array
	TString* hist_name_arr = new TString[num_hists];
	GetHistNames( (TList*)f->GetListOfKeys(), hist_name_arr );
	
	// Define a histogram pointer
	TH1D* h_spec;
	
	
	
	// EXTRACT YIELDS -------------------------------------------------------------------------- //
	if ( spectrum_type < 2 ){
		for ( Int_t i = 0; i < num_hists; i++ ){
			h_spec = (TH1D*)f->Get( hist_name_arr[i] );
			ExtractYieldsHist( h_spec, i, pos, spectrum_type );
		}
	}
	else if ( spectrum_type == 2 ){
		h_spec = (TH1D*)f->Get( hist_name_arr[0] );
		ExtractYieldsHist( h_spec, 5, pos, spectrum_type );
	}
	else{
		std::cout << "Invalid spectrum type" << "\n";
	}
	
	
	
	// Close the TFile and delete all from memory
	if ( SW_BATCH_MODE == 1 ){
	
		f->Close();
		if ( f->IsOnHeap() ){ f->Delete(); }
	}
	
	
	return;
}

		

/* TODO
 * Make axis label font bigger
 * Start from -0.5
*/



















































