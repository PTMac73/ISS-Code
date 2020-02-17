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
};

void SetFitOptions( FitPeakOptions_t &opt, Int_t num_peaks, Int_t peak_num, Double_t sig_est = 0.25*PEAK_WIDTH_ESTIMATE, Double_t sig_range = 0.2*PEAK_WIDTH_ESTIMATE ){
	opt.num_peaks = num_peaks;
	opt.peak_num = peak_num;
	opt.sig_est = sig_est;
	opt.sig_range = sig_range;
	return;
}

TF1* EstimatePeakParameters( TH1D* h, FitPeakOptions_t &opt, Int_t**& var_type_arr ){
	// Declare variables
	Double_t amp_est;

	// Populate variable type array
	var_type_arr = new Int_t*[opt.num_peaks + 1];
	for ( Int_t i = 0; i < opt.num_peaks + 1; i++ ){
		var_type_arr[i] = new Int_t[3];
		for ( Int_t k = 0; k < 3; k++ ){
			var_type_arr[i][k] = -1;
		}
	}
	
	// Generate the fit
	TString fit_str = GetFitString( opt.num_peaks, BG_DIM, var_type_arr, opt.peak_num );
	TF1* fit_func = new TF1( "fit_func", fit_str, PRE_MIN, PRE_MAX );
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
				amp_est = GetAmpEstimate( h, peak_energies[i + opt.peak_num - 1] - 0.5*PEAK_WIDTH_ESTIMATE, peak_energies[i + opt.peak_num - 1] + 0.5*PEAK_WIDTH_ESTIMATE );
				std::cout << "AMP-EST (" << i-1 << "):" << amp_est << "\n";
				fit_func->SetParameter( var_type_arr[i][j], amp_est );
				fit_func->SetParLimits( var_type_arr[i][j], 0.0*amp_est, 1.2*amp_est );
			}
			
			// Mus
			else if ( i > 0 && j == 1 ){
				if ( peak_fix_positions[i + opt.peak_num - 1] == 1 ){
					fit_func->FixParameter( var_type_arr[i][j], peak_energies[i + opt.peak_num - 1]);
				}
				else{
					fit_func->SetParameter( var_type_arr[i][j], peak_energies[i + opt.peak_num - 1]);
					fit_func->SetParLimits( var_type_arr[i][j], peak_energies[i + opt.peak_num - 1] - PEAK_WIDTH_ESTIMATE, peak_energies[i + opt.peak_num - 1] + PEAK_WIDTH_ESTIMATE );
				}
			}
			
			// Sigmas
			else if ( i > 0 && j == 2 ){
			
				// Ensure the fixed width has not already been set
				// Fixed-width peaks:
				if ( var_type_arr[i][j] > var_type_arr[i][j-1] && peak_fix_widths[i + opt.peak_num - 1] == 1 ){
					fit_func->SetParameter( var_type_arr[i][j], opt.sig_est );
					fit_func->SetParLimits( var_type_arr[i][j], opt.sig_est - opt.sig_range, opt.sig_est + 5*opt.sig_range );
				}
				
				// Variable-width peaks
				else if ( var_type_arr[i][j] > var_type_arr[i][j-1] && peak_fix_widths[i + opt.peak_num - 1] == 0 ){
					fit_func->SetParameter( var_type_arr[i][j], opt.sig_est );
					fit_func->SetParLimits( var_type_arr[i][j], opt.sig_est, opt.sig_est + PEAK_WIDTH_ESTIMATE );
				}
					
			}
			
		}
		
	}
	return fit_func;
}
#endif
