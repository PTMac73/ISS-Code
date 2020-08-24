// EY_FitPeaks.h
// Fits peaks for the ExtractYields.C script
// ============================================================================================= //
// Patrick MacGregor
// Nuclear Physics Research Group
// Department of Physics and Astronomy
// The University of Manchester
// ============================================================================================= //

#ifndef EY_FIT_PEAKS_H_
#define EY_FIT_PEAKS_H_

#include "EY_HistogramFunctions.h"

#include <TF1.h>
#include <TH1.h>

struct FitPeakOptions_t{
	Int_t num_peaks;
	Int_t peak_num;
	Double_t sig_est;
	Double_t sig_range;
	Bool_t fix_widths;
	Bool_t fix_positions;
	const Double_t* pe;
};

void SetFitOptions( FitPeakOptions_t &opt, Int_t num_peaks, Int_t peak_num, Bool_t fix_widths, Bool_t fix_positions = 0, const Double_t* pe = peak_energies, Double_t sig_est = 0.25*PEAK_WIDTH_ESTIMATE, Double_t sig_range = 0.2*PEAK_WIDTH_ESTIMATE ){
	opt.num_peaks = num_peaks;
	opt.peak_num = peak_num;
	opt.sig_est = sig_est;
	opt.sig_range = sig_range;
	opt.fix_widths = fix_widths;
	opt.fix_positions = fix_positions;
	opt.pe = pe;
	return;
}

TF1* EstimatePeakParameters( TH1D* h, FitPeakOptions_t &opt, Int_t**& var_type_arr, TString fit_name = "fit_func", Double_t lb = 0, Double_t ub = 10 ){
	// Declare variables
	Double_t amp_est;

	// Populate variable type array with -1
	var_type_arr = new Int_t*[opt.num_peaks + 1];
	for ( Int_t i = 0; i < opt.num_peaks + 1; i++ ){
		var_type_arr[i] = new Int_t[3];
		for ( Int_t k = 0; k < 3; k++ ){
			var_type_arr[i][k] = -1;
		}
	}
	
	// Generate the fit
	TString fit_str = GetFitString( opt.num_peaks, BG_DIM, var_type_arr, opt.peak_num, opt.fix_widths );
	TF1* fit_func = new TF1( fit_name.Data(), fit_str, lb, ub );
	fit_func->SetLineColor(kBlack);
	fit_func->SetNpx(2000);

	// Loop over number of peaks + 1 (i)
	for ( Int_t i = 0; i < opt.num_peaks + 1; i++ ){
	
		// Loop over number of parameters (j)
		for ( Int_t j = 0; j < 3; j++ ){
		
			// Backgrounds
			if ( i == 0 && var_type_arr[i][j] >= 0 ){
				fit_func->SetParameter( var_type_arr[i][j], 0.01 );
				fit_func->SetParLimits( var_type_arr[i][j], 0, TMath::Power( 0.1*h->GetMaximum(), j+1) );
			}
			
			// Amplitudes
			else if ( i > 0 && j == 0 ){
				amp_est = GetAmpEstimate( h, opt.pe[i + opt.peak_num - 1] - 0.5*PEAK_WIDTH_ESTIMATE, opt.pe[i + opt.peak_num - 1] + 0.5*PEAK_WIDTH_ESTIMATE );
				fit_func->SetParameter( var_type_arr[i][j], amp_est );
				fit_func->SetParLimits( var_type_arr[i][j], 0.0*amp_est, 1.2*amp_est );
			}
			
			// Mus
			else if ( i > 0 && j == 1 ){
				if ( opt.fix_positions == 1 ){
					fit_func->FixParameter( var_type_arr[i][j], opt.pe[i + opt.peak_num - 1]);
				}
				else{
					fit_func->SetParameter( var_type_arr[i][j], opt.pe[i + opt.peak_num - 1]);
					fit_func->SetParLimits( var_type_arr[i][j], opt.pe[i + opt.peak_num - 1] - PEAK_WIDTH_ESTIMATE, opt.pe[i + opt.peak_num - 1] + PEAK_WIDTH_ESTIMATE );
				}
			}
			
			// Sigmas
			else if ( i > 0 && j == 2 ){
			
				// Ensure the fixed width has not already been set
				// Fixed-width peaks:
				if ( var_type_arr[i][j] > var_type_arr[i][j-1] && opt.fix_widths == 1 ){
					fit_func->SetParameter( var_type_arr[i][j], opt.sig_est );
					fit_func->SetParLimits( var_type_arr[i][j], opt.sig_est, opt.sig_est + 5*opt.sig_range );
				}
				
				// Variable-width peaks
				else if ( var_type_arr[i][j] > var_type_arr[i][j-1] && opt.fix_widths == 0 ){
					fit_func->SetParameter( var_type_arr[i][j], opt.sig_est );
					fit_func->SetParLimits( var_type_arr[i][j], opt.sig_est, TMath::Sqrt( TMath::Power( opt.sig_est, 2 ) + TMath::Power( MAX_SF_ESTIMATE*1.21, 2 ) ) );
				}
					
			}
			
		}
		
	}
	return fit_func;
}

void PrintFF( TF1* fit_func, Int_t** vta, Int_t num_peaks, std::ostream& f = std::cout ){
	// Define first-column width
	Int_t w = 10;
	Int_t prec = 6;
	TString sp = "    ";
	TString bl = "-";
	Double_t lb = 0;
	Double_t ub = 0;
	
	// Print header
	f << "\n>>> Fitting function parameters <<<\n";
	std::cout << "FORMULA: " << fit_func->GetExpFormula() << "\n";
	f << std::left << std::setw(w) << std::right << "KEY :" << std::left << sp
		<< std::setw(w) << "AMP" << sp << std::setw(w) << "[lb]" << sp << std::setw(w) << "[ub]" << sp
		<< std::setw(w) << "MU"  << sp << std::setw(w) << "[lb]" << sp << std::setw(w) << "[ub]" << sp
		<< std::setw(w) << "SIG"  << sp << std::setw(w) << "[lb]" << sp << std::setw(w) << "[ub]" << sp << "\n";
	
	for ( Int_t i = 0; i < num_peaks + 1; i++ ){
		// Print first column
		f << std::setw(w) << std::right << ( i == 0 ? "BG :" : Form( "P%i :", i-1 ) ) << std::left << sp;
	
		// Print other columns
		for ( Int_t j = 0; j < 3; j++ ){
			// Get parameter limits
			fit_func->GetParLimits( vta[i][j], lb, ub );
			
			// Print the section
			if ( vta[i][j] != -1 ){
				f << std::setprecision(prec) << std::setw(w) << fit_func->GetParameter( vta[i][j] ) << sp
					<< std::setprecision(prec) << std::setw(w) << lb << sp
					<< std::setprecision(prec) << std::setw(w) << ub << sp;
			}
			else{ 
				f << std::setprecision(prec) << std::setw(w) << bl << sp
					<< std::setprecision(prec) << std::setw(w) << bl << sp
					<< std::setprecision(prec) << std::setw(w) << bl << sp;
			}
		}
		
		// Finish the line
		f << "\n";
		
	}
	f << "\n";
	
	return;
}















#endif
