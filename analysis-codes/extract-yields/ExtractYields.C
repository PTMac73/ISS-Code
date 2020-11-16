// ExtractYields.C
// Extracts the yields for a given histogram file provided when creating the Mg spectra
// 3 regions done here:
//	> Bound states
//	> Unbound doublet + 1
//	> Super unbound states
//
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

void ExtractYieldsHist( TH1D* h, Int_t pos ){
	// DECLARE SOME VARIABLES ------------------------------------------------------------------ //
	TCanvas* c_spec;			// Canvas pointers
	TCanvas* c_pre_fit;
	Double_t sig_est = 0;
	Double_t sig_range = 0;
	TFitResultPtr fit_result_ptr;
	
	// Remove the stats box
	gStyle->SetOptStat(0);
	
	// Define canvas
	c_pre_fit = new TCanvas( "c_pre_fit", "PRE-FIT PEAKS", CANVAS_WIDTH, CANVAS_HEIGHT );
	c_pre_fit->cd();
	SetPadMargins( (TPad*)c_pre_fit->GetPad(0) );

	// Format the histogram
	if ( SW_HIST_TITLE == 1 ){ h->SetTitle( GenerateHistTitle( pos ) ); }
	else{ h->SetTitle(""); }
	FormatHistogram( h );

	
	// FIT THE TWO STRONG PEAKS IF POSSIBLE ------------------------------------------------ //
	// Estimate the peak parameters
	FitPeakOptions_t pre_fit_opt;
	SetFitOptions( pre_fit_opt, 2, PRE_FIX_WID, PRE_FIX_POS, pre_peak_energies );
	Int_t** pre_var_type_arr;
	TF1* pre_fit_func = EstimatePeakParameters( h, pre_fit_opt, pre_var_type_arr, "pre_fit_func", PRE_MIN, PRE_MAX );

	pre_fit_func->SetLineWidth(2);
	pre_fit_func->SetLineColor(kBlack);
	
	// Fit the spectrum
	h->Fit( pre_fit_func, "RBQ" );
	
	// Draw if desired
	h->Draw();
	c_pre_fit->Modified(); c_pre_fit->Update();
	c_pre_fit->Print( GenerateFileName( pdf, pos, 1 ) );
	
	// Calculate the width
	sig_est = pre_fit_func->GetParameter( pre_var_type_arr[1][2] );
	sig_range = pre_fit_func->GetParError( pre_var_type_arr[1][2] );
	
	std::cout << "SIG-EST  : " << std::setprecision(8) << sig_est << "\n";
	std::cout << "SIG-RANGE: " << std::setprecision(8) << sig_range << "\n";


	// DEFINE CONSTANTS ------------------------------------------------------------------------ //
	FitPeakOptions_t fit_opt;	// Holds the fitting options for the peaks
	Int_t** var_type_arr;		// Holds the type for each variable number (e.g. var 2 is a amplitude etc.)
	TF1* fit_func;				// Fitting function for the whole region
	Peak_t peaks[NUM_PEAKS];	// Holds info on the peaks
	TF1* f_peaks[NUM_PEAKS];	// Fitting function for each peak
	
	
	
	// Open the output file
	std::ofstream out_file;
	
	// Delete the old data file if it exists
	if ( std::remove( GenerateFileName( dat, pos ) ) != 0 ){
		std::cout << "Error deleting data file.\n";
	}
	out_file.open( GenerateFileName( dat, pos ), std::ofstream::app );
	
	// DRAW A CANVAS AND DRAW THE INITIAL SPECTRUM --------------------------------------------- //
	// Define a canvas and set margins
	c_spec = new TCanvas( GenerateCanvasName( pos ), "FITTED SPECTRUM", CANVAS_WIDTH, CANVAS_HEIGHT );
	SetPadMargins( (TPad*)c_spec->GetPad(0) );
	c_spec->cd();
	
	// Draw the histogram
	h->Draw();
	

	// START FULL FITTING PROCESS HERE ----------------------------------------------------- //
	// Estimate the peak parameters
//	SetFitOptions( FitPeakOptions_t &opt, Int_t num_peaks, Bool_t* fix_widths, Bool_t* fix_positions, const Double_t* pe, Double_t sig_est, Double_t sig_range )
	SetFitOptions( fit_opt, NUM_PEAKS, fix_widths, fix_positions, peak_energies, sig_est, sig_range );
	
	// Define fitting function
	fit_func = EstimatePeakParameters( h, fit_opt, var_type_arr, "fit_func", E_LIMITS[0], E_LIMITS[1] );

	// Format the fit function
	fit_func->SetLineColor( kBlack );
	fit_func->SetLineWidth(2);
	
	// Fit the histogram
	fit_result_ptr = h->Fit( fit_func, "S" );
	PrintFF( fit_func, var_type_arr, NUM_PEAKS );
	
	// Print the peak information
	PrintFF( fit_func, var_type_arr, NUM_PEAKS, out_file );
	PrintPeakHeader( out_file );

	// LOOP OVER THE PEAKS TO DO USEFUL PEAK-SPECIFIC THINGS ------------------------------- //
	for ( Int_t j = 0; j < NUM_PEAKS; j++ ){
	
		// Calculate peak quantities and print the peaks to file ( including the GS doublet for comparison)
		CalculatePeakQuantities( fit_result_ptr, BG_DIM, var_type_arr[j+1][0], var_type_arr[j+1][2], peaks[j] ); 
		PrintPeak( peaks[j], out_file );
			
		f_peaks[j] = new TF1( Form( "fit_func_%i", j ), GetFitStringSimple( 1, BG_DIM ), peak_energies[j] - 4*PEAK_WIDTH_ESTIMATE, peak_energies[j] + 4*PEAK_WIDTH_ESTIMATE );
		f_peaks[j]->SetNpx(2000);
		
		// Set individual backgrounds
		for ( Int_t k = 0; k < BG_DIM + 1; k++ ){
			f_peaks[j]->FixParameter( k, fit_func->GetParameter( var_type_arr[0][k] ) );
		}
		
		// Set individual peak shapes
		for ( Int_t k = 0; k < 3; k++ ){
			f_peaks[j]->FixParameter( k + BG_DIM + 1, fit_func->GetParameter(var_type_arr[j+1][k]) );
		}
		
		// Format individual peaks
		f_peaks[j]->SetLineColor(peak_colours[j]);
		f_peaks[j]->SetLineWidth(1);
		f_peaks[j]->Draw("SAME");

		// Update the canvas
		c_spec->Modified(); c_spec->Update();
		
		
	} // LOOP OVER PEAKS (j)

	// Print the last line
	PrintFitPars( fit_result_ptr, BG_DIM, out_file );

	
	// Add a line for the neutron separation energy
	TLine* lNSE = new TLine( NEUTRON_SEP_ENERGY, 0, NEUTRON_SEP_ENERGY, c_spec->GetUymax() );
	lNSE->SetLineColor(kBlack);
	lNSE->SetLineStyle(2);
	lNSE->Draw("SAME");
	
	TLatex* tNSE = new TLatex( 1.05*NEUTRON_SEP_ENERGY, 0.9*h->GetMaximum(), Form( "S_{n} = %4.3f MeV", NEUTRON_SEP_ENERGY ) );
	tNSE->SetTextAlign(12);
	tNSE->SetTextFont(132);
	tNSE->Draw("SAME");
	
	TText* peak_label_gs[5];
	TText* peak_label_ubd[4];
	
	// Number the peaks
	TText* peak_labels[NUM_PEAKS][2];
	for ( Int_t i = 0; i < NUM_PEAKS; i++ ){
		if ( i == 0 ){
			for ( Int_t j = 0; j < 5; j++ ){
				peak_label_gs[j] = new TText(peaks[i].mu + peak_label_pos_offset[i][0], peaks[i].amp + peak_label_pos_offset[i][1] + 10*(4-j), peak_label_gs_str[j] );
				peak_label_gs[j]->SetTextAlign(22);
				peak_label_gs[j]->SetTextFont(132);
				peak_label_gs[j]->SetTextSize(0.03);
				peak_label_gs[j]->Draw("SAME");
			}
		}
		else if ( i == 7 ){
			for ( Int_t j = 0; j < 4; j++ ){
				peak_label_ubd[j] = new TText(peaks[i].mu + peak_label_pos_offset[i][0], peaks[i].amp + peak_label_pos_offset[i][1] + 10*(3-j), peak_label_ubd_str[j] );
				peak_label_ubd[j]->SetTextAlign(22);
				peak_label_ubd[j]->SetTextFont(132);
				peak_label_ubd[j]->SetTextSize(0.03);
				peak_label_ubd[j]->Draw("SAME");
			}
		}
		else if ( i == 8 ){
			// Do nothing
		}
		else{
			for ( Int_t j = 0; j < 2; j++ ){
				peak_labels[i][j] = new TText(peaks[i].mu + peak_label_pos_offset[i][0], peaks[i].amp + peak_label_pos_offset[i][1] + 10*(1-j), peak_label[i][j] );
				peak_labels[i][j]->SetTextAlign(22);
				peak_labels[i][j]->SetTextFont(132);
				peak_labels[i][j]->SetTextSize(0.03);
				peak_labels[i][j]->Draw("SAME");
			}
		}
		
		
	}
	
	
	// Update the canvas
	TPad* pad = (TPad*)c_spec->GetPad(0);
	pad->SetTicks(1,1);
	c_spec->Modified(); c_spec->Update();

	// Print the spectrum
	c_spec->Print( GenerateFileName( pdf, pos ) );
	c_spec->Print( GenerateFileName( tex, pos ) );
	
	// Append new line and close the file
	out_file.close();
	
	
	// Delete arrays that are no longer needed and free them up for the next loop around
	for ( Int_t j = 0; j < 3; j++ ){
		delete[] pre_var_type_arr[j];
	}
	delete[] pre_var_type_arr;
	
	
	for ( Int_t j = 0; j < NUM_PEAKS+1; j++ ){
		delete[] var_type_arr[j];
		if ( j < NUM_PEAKS ){
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
	
	return;
}
	
// ---===---===---===---===---===---===---===---===---===---===---===---===---===---===---===--- //
void ExtractYields( Int_t pos = 1 ){
	// Declare variables
	TFile* f;						// File containing the histograms
	
	// Batch mode
	if ( SW_BATCH_MODE == 1 ){ gROOT->SetBatch(kTRUE); }
	
	// Open the relevant TFile
	f = new TFile( GetROOTFileName( pos ), "READ" );

	// Check that it opened
	if ( f->IsOpen() ){ std::cout << "File opened successfully!\n"; }
	else{
		std::cout << "File " << GetROOTFileName( pos ) << " failed to open. Exiting process...\n";
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
	h_spec = (TH1D*)f->Get( hist_name_arr[0] );
	ExtractYieldsHist( h_spec, pos );
	
	// Close the TFile and delete all from memory
	if ( SW_BATCH_MODE == 1 ){
	
		f->Close();
		if ( f->IsOnHeap() ){ f->Delete(); }
	}
	
	
	return;
}

