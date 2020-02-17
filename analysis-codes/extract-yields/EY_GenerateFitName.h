// EY_GenerateFitName.h
// Generates the names of fits
// ============================================================================================= //
// Patrick MacGregor
// Nuclear Physics Research Group
// Department of Physics and Astronomy
// The University of Manchester
// ============================================================================================= //

#ifndef EY_GENERATE_FIT_NAME_H_
#define EY_GENERATE_FIT_NAME_H_

#include "ExtractYields.h"
#include <iostream>
#include <TString.h>


// GENERATE GAUSSIAN STRINGS ------------------------------------------------------------------- //
TString GetGaussianString( Int_t par_num_begin = 0 ){
	return Form("[%d]*exp(-0.5*((x-[%d])/[%d])**2)", par_num_begin, par_num_begin + 1, par_num_begin + 2 );
}

TString GetGaussianStringFW( Int_t par_num_begin = 0, Int_t width_var_num = 2 ){
	return Form("[%d]*exp(-0.5*((x-[%d])/[%d])**2)", par_num_begin, par_num_begin + 1, width_var_num );
}



// GENERATE POLYNOMIAL STRINGS ----------------------------------------------------------------- //
TString GetPolynomialString( Int_t par_num_begin = 0, Int_t dim = 0 ){
	TString poly_str = "";

	for ( Int_t i = 0; i < dim + 1; i++ ){
		if ( i == 0 ){ poly_str += Form( "[%d]", par_num_begin + i ); }
		if ( i == 1 ){ poly_str += Form( "[%d]*x", par_num_begin + i ); }
		if ( i >= 2 ){ poly_str += Form( "[%d]*x**%d", par_num_begin + i, i ); }
		if ( i != dim ){ poly_str += " + "; }
	}
	return poly_str;
}



// GET TOTAL NUMBER OF VARIABLES TO USE IN FIT WITH RELATIVE WIDTHS ---------------------------- //
Int_t GetNumFitVars( Int_t num_peaks, Int_t bg_dim = 0 ){
	// Calculate the number of variables
	Int_t num_pars = bg_dim + 1;
	for ( Int_t i = 0; i < num_peaks; i++ ){
		if ( peak_fix_widths[i] == 0 ){ num_pars += 3; }
		else{ num_pars += 2; }
	}
	return num_pars;
}



// GET INDICES OF DIFFERENT FEATURES IN ARRAY -------------------------------------------------- //
Int_t CountFeatures( Int_t feature_num, Int_t num_vars, Int_t* in_arr ){
	// Count the number that appear
	Int_t num_features = 0;
	for ( Int_t i = 0; i < num_vars; i++ ){
		if ( in_arr[i] == feature_num ){
			num_features++;
		}
	}
	return num_features;
}

void GetFeatureIndices( Int_t feature_num, Int_t num_vars, Int_t* in_arr, Int_t* out_arr ){
	Int_t out_ctr = 0;
	for ( Int_t i = 0; i < num_vars; i++ ){
		if ( in_arr[i] == feature_num ){
			out_arr[out_ctr] = i;
			out_ctr++;
		}
	}

}


// GENERATE TOTAL FIT STRINGS (SIMPLE) --------------------------------------------------------- //
TString GetFitStringSimple( Int_t num_peaks, Int_t bg_dim = 0 ){
	TString fit_str = "";

	// Add background
	fit_str += GetPolynomialString( 0, bg_dim );
	fit_str += " + ";

	// Add number of Gaussian peaks
	for ( Int_t i = 0; i < num_peaks; i++ ){
		fit_str += GetGaussianString(3*i + bg_dim + 1);
		if ( i != num_peaks - 1 ){ fit_str += " + "; }
	}

	return fit_str;
}



/* VARIABLE TYPE ARRAY KEY:
	- 0 --> Background
	- 1 --> Amplitude
	- 2 --> Mean
	- 3 --> Sigma
	- 4 --> Fixed sigma
	
	SIZE = ( n + 1 ) x 3, where n is number of peaks:
	- BG stored in slot 0 (quadratic at most)
	- P1 stored in slot 1
	- P2 stored in slot 2
	  ...
	- Pn stored in slot n
	
	STORES: indices of the variables, not the variables themselves!
*/

// PRINT VARIABLE TYPE ARRAY ------------------------------------------------------------------- //
void PrintVTA( Int_t** arr, Int_t num_peaks, std::ostream& f = std::cout ){
	// Define first-column width
	Int_t w = 7;
	
	// Header
	f << std::setw(w) << "KEY : " << std::setw(4) << "AMP" << "\t" << std::setw(4) << "MU" << "\t" << std::setw(4) << "SIG" << "\n";
	
	
	for ( Int_t i = 0; i < num_peaks + 1; i++ ){
		// Print first column
		if ( i == 0 ){ f << std::setw(w) << "BG : "; }
		else{ f << std::setw(w) << Form( "P%i : ", i-1 ); }
	
		// Print other columns
		for ( Int_t j = 0; j < 3; j++ ){
			if ( arr[i][j] != -1 ){ f << std::setw(4) << arr[i][j] << "\t" ; }
			else{ f << std::setw(4) << "X" << "\t"; }
		}
		
		// Finish the line
		f << "\n";
	}
	
	return;
}


// GENERATE TOTAL FIT STRINGS ------------------------------------------------------------------ //
TString GetFitString( Int_t num_peaks, Int_t bg_dim, Int_t** variable_type_arr, Int_t peak_num = 0 ){
	TString fit_str = "";
	
	// Add background to string
	fit_str += GetPolynomialString( 0, bg_dim );
	fit_str += " + ";
	Int_t var_ctr = 0;
	Int_t width_var = -1;
	Bool_t b_fill_arr = 1;
	Bool_t b_add_one_more = 0;
	if( variable_type_arr == NULL ){ b_fill_arr = 0; }
	
	// Add background to variable counter type
	for ( Int_t i = 0; i < 3; i++ ){
		if ( i < bg_dim + 1 ){ 
			if ( b_fill_arr ){ variable_type_arr[0][i] = var_ctr; }
			var_ctr++; 
		}
		else{ if ( b_fill_arr ){ variable_type_arr[0][i] = -1; } }
	}	

	// Add Gaussian peaks with some constrained widths
	for ( Int_t i = 0; i < num_peaks; i++ ){
		if ( peak_fix_widths[i + peak_num] == 0 ){
			// Don't fix the widths
			fit_str += GetGaussianString( var_ctr );
			for ( Int_t j = 0; j < 3; j++ ){
				if ( b_fill_arr ){ variable_type_arr[i+1][j] = var_ctr; }
				var_ctr++;
			}
		}
		else{
			// Fix the widths of the peaks
			// Check if the fixed width has been defined
			if ( width_var == -1 ){ width_var = var_ctr + 2; b_add_one_more = 1; }
			fit_str += GetGaussianStringFW( var_ctr, width_var );
			
			// Define the variable type array
			for ( Int_t j = 0; j < 2; j++ ){
				if ( b_fill_arr ){ variable_type_arr[i+1][j] = var_ctr; }
				var_ctr++;
			}
			
			// Add the fixed width to the array
			if ( b_fill_arr ){ variable_type_arr[i+1][2] = width_var; }
			if ( b_add_one_more == 1 ){ var_ctr++; b_add_one_more = 0; }
		}
		
		// Add a plus sign if more terms to be added
		if ( i != num_peaks - 1 ){ fit_str += " + "; }
	}
	return fit_str;
}

#endif
