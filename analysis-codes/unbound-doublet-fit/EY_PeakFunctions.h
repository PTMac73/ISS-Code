// EY_PeakFunctions.h
// Functions for calculating things with peaks
// ============================================================================================= //
// Patrick MacGregor
// Nuclear Physics Research Group
// Department of Physics and Astronomy
// The University of Manchester
// ============================================================================================= //
// These all require a TFitResultPtr, which assumes the following indices:
/*	* 0,1,...,n 			--> Background polynomial
	* n+1,n+2,n+3 			--> Peak 1
	* n+4,n+5,n+6			--> Peak 2
	...
	* n+3k+1,n+3k+2,n+3k+3	--> Peak k
*/

#ifndef EY_PEAK_FUNCTIONS_H_
#define EY_PEAK_FUNCTIONS_H_

#include <TMath.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TMatrixD.h>

#include <iostream>

const Int_t print_width = 12;
const Int_t print_precision = 10;

// ============================================================================================= //
// ============================================================================================= //
// ============================================================================================= //
// PEAK FUNCTIONS
struct Peak_t{
	Double_t amp;		// Amplitude
	Double_t mu;		// Mu
	Double_t sig;		// Sigma
	Int_t bg_dim;		// Dimension of bg
	Double_t* bg;		// Background (pointer to array)
	Double_t x;			// Position
	Double_t w;			// Width
	Double_t h;			// Height
	Double_t a;			// Area
	Double_t c;			// Centroid
	Double_t amp_err;
	Double_t mu_err;
	Double_t sig_err;
	Double_t* bg_err;
	Double_t x_err;
	Double_t w_err;
	Double_t h_err;
	Double_t a_err;
	Double_t c_err;
};


// --------------------------------------------------------------------------------------------- //
Double_t PeakWidth( Double_t sig ){
	return 2*sig*TMath::Sqrt( 2*TMath::Log(2) );
}


// --------------------------------------------------------------------------------------------- //
Double_t PeakHeight( Double_t amp, Double_t mu, Int_t bg_dim, Double_t *bg ){
	Double_t h = amp;
	for ( Int_t j = 0; j < bg_dim + 1; j++ ){
		h += -bg[j]*TMath::Power(mu,j);
	}
	return h;
}


// --------------------------------------------------------------------------------------------- //
Double_t PeakArea( Double_t amp, Double_t sig ){
	// Need to divide by the bin width in the units of the energy-axis (MeV)
	return amp*sig*TMath::Sqrt( TMath::TwoPi() )/( 0.02 );
}


// --------------------------------------------------------------------------------------------- //
Double_t PeakWidthErr( Double_t sig_err ){
	return 2*sig_err*TMath::Sqrt(2*TMath::Log(2));
}


// --------------------------------------------------------------------------------------------- //
Double_t PeakHeightErr( TFitResultPtr r, Int_t peak_index, Int_t bg_dim ){
	// bg_array_indices[bg_dim] = { 0, 1, ..., bg_dim }
	Double_t var = 0;

	// Define amplitude and mean more accessibly
	Double_t amp = r->Parameter( peak_index );
	Double_t mu = r->Parameter( peak_index + 1 );
	Double_t amp_err = r->ParError( peak_index );
	Double_t mu_err = r->ParError( peak_index + 1 );

	// Define covariance matrix
	TMatrixDSym cov = r->GetCovarianceMatrix();

	// Get the background values and errors
	Double_t bg[bg_dim+1], bg_err[bg_dim+1];
	for ( Int_t i = 0; i <= bg_dim; i++ ){
		bg[i] = r->Parameter( i );
		bg_err[i] = r->ParError( i );
	}

	// Calculate the derivative of h with mu:
	Double_t M = 0;
	for ( Int_t i = 1; i <= bg_dim; i++ ){
		M += i*TMath::Power( mu, i-1 )*bg[i];
	}

	// Add the non-summed terms
	var += TMath::Power( amp_err, 2 );
	var += TMath::Power( M*mu_err, 2 );
	var += -2*M*TMatrixDRow( cov, peak_index )( peak_index + 1 );

	// Add the summed terms
	for ( Int_t i = 0; i <= bg_dim; i++ ){
		var += TMath::Power( TMath::Power( mu, i )*bg_err[i], 2 );
		var += -2*TMath::Power( mu, i )*TMatrixDRow( cov, peak_index )( i );
		var += 2*TMath::Power( mu, i )*M*TMatrixDRow( cov, peak_index + 1 )( i );
		for ( Int_t k = i + 1; k <= bg_dim; k++ ){
			var += TMath::Power( mu, i+k )*TMatrixDRow( cov, i )( k );
		}
	}
	return TMath::Sqrt(var);
}


// --------------------------------------------------------------------------------------------- //
Double_t PeakAreaErr( TFitResultPtr r, Int_t peak_index, Int_t sig_index ){
	Double_t amp = r->Parameter( peak_index );
	Double_t sig = r->Parameter( sig_index );
	Double_t amp_err = r->ParError( peak_index );
	Double_t sig_err = r->ParError( sig_index );
	TMatrixDSym cov = r->GetCovarianceMatrix();

	return TMath::Sqrt( TMath::TwoPi()*( TMath::Power( amp*sig_err, 2 ) + TMath::Power( sig*amp_err, 2 ) + 2*amp*sig*TMatrixDRow( cov, peak_index )( sig_index ) ) )/0.02;
}


// --------------------------------------------------------------------------------------------- //
void CalculatePeakQuantities( TFitResultPtr r, Int_t bg_dim, Int_t peak_index, Int_t sig_index, Peak_t &p ){
	p.amp = r->Parameter(peak_index);
	p.mu = r->Parameter(peak_index + 1);
	p.sig = r->Parameter(sig_index);
	
	p.amp_err = r->ParError(peak_index);
	p.mu_err = r->ParError(peak_index + 1);
	p.sig_err = r->ParError(sig_index);

	p.bg_dim = bg_dim;
	Double_t* bg = new Double_t[bg_dim];
	Double_t* bg_err = new Double_t[bg_dim];
	for ( Int_t i = 0; i < bg_dim; i++ ){
		bg[i] = r->Parameter(i);
		bg_err[i] = r->ParError(i);
	}
	p.bg = bg;
	p.bg_err = bg_err;
	
	p.x = p.mu; p.x_err = p.mu_err;
	p.c = p.mu; p.c_err = p.mu_err;

	p.h = PeakHeight( p.amp, p.mu, p.bg_dim, p.bg );
	p.h_err = PeakHeightErr( r, peak_index, p.bg_dim );

	p.w = PeakWidth( p.sig );
	p.w_err = PeakWidthErr( p.sig_err );

	p.a = PeakArea( p.amp, p.sig );
	p.a_err = PeakAreaErr( r, peak_index, sig_index );
	
	return;
}


// --------------------------------------------------------------------------------------------- //
void PrintPeak( Peak_t &p, std::ostream& f = std::cout ){
	
	f << std::setprecision(print_precision) << 
	std::setw(print_width) << p.x << "\t" << 
	std::setw(print_width) << p.x_err << "\t" << 
	std::setw(print_width) << p.w << "\t" << 
	std::setw(print_width) << p.w_err << "\t" << 
	std::setw(print_width) << p.h << "\t" << 
	std::setw(print_width) << p.h_err << "\t" << 
	std::setw(print_width) << p.a << "\t" << 
	std::setw(print_width) << p.a_err << "\t" << 
	std::setw(print_width) << p.c << "\t" << 
	std::setw(print_width) << p.c_err << "\n";
	
	return;
}

void PrintPeakHeader( std::ostream& f = std::cout ){
	
	f <<  
	std::setw(print_width) << "POSITION" << 
	std::setw(print_width) << "" << "\t" << 
	std::setw(print_width) << "WIDTH" << "\t" << 
	std::setw(print_width) << "" << "\t" << 
	std::setw(print_width) << "HEIGHT" << "\t" << 
	std::setw(print_width) << "" << "\t" << 
	std::setw(print_width) << "AREA" << "\t" << 
	std::setw(print_width) << "" << "\t" << 
	std::setw(print_width) << "CENTROID" << "\t" << 
	std::setw(print_width) << "" << "\n";
	
	return;
}


// ============================================================================================= //
// ============================================================================================= //
// ============================================================================================= //
// REDUCED PEAK FUNCTIONS
struct RedPeak_t{
	Double_t* bg;	// Background array
	Double_t* bg_err;
	Int_t bg_dim;	// Background dimension
	Double_t lb;	// Lower bound for fitting
	Double_t ub;	// Upper bound for fitting
	Double_t a;		// Area
	Double_t c;		// Centroid
	Double_t a_err;
	Double_t c_err;
	
};


// --------------------------------------------------------------------------------------------- //
Double_t CalculateBackground( Double_t x, Int_t bg_dim, Double_t* bg_pars ){
	Double_t sum = 0;
	for ( Int_t i = 0; i < bg_dim + 1; i++ ){
		sum+= bg_pars[i]*TMath::Power(x,i);
	}
	return sum;
}


// --------------------------------------------------------------------------------------------- //
void CalculateRedPeakQuantities( TFitResultPtr r, Int_t bg_dim, Double_t lb, Double_t ub, TH1D* h, RedPeak_t &rp ){
	// Store the background stuff
	rp.bg_dim = bg_dim;
	rp.bg = new Double_t[bg_dim+1];
	rp.bg_err = new Double_t[bg_dim+1];
	
	for ( Int_t i = 0; i < bg_dim + 1; i++ ){
		rp.bg[i] = r->Parameter(i);
		rp.bg_err[i] = r->ParError(i);
	}
	
	rp.lb = lb;
	rp.ub = ub;
	
	// Calculate the area and centroid
	Double_t cent_num = 0, cent_den = 0;
	Double_t x;
	Double_t num_var = 0, den_var = 0;
	Double_t num_var_tot = 0, den_var_tot = 0;
	TMatrixDSym cov = r->GetCovarianceMatrix();
	
	// Loop over the bins and count the number of entries (inclusive)
	for ( Int_t i = h->FindBin(lb); i <= h->FindBin(ub); i++ ){
		
		// Find the bin position
		x = h->GetBinCenter(i);
		
		// Calculate the yield without the background
		cent_den += h->GetBinContent(i) - CalculateBackground( x, bg_dim, rp.bg );
		cent_num += x*( h->GetBinContent(i) - CalculateBackground( x, bg_dim ,rp.bg ) );
		
		// Caclulate the variance
		den_var = 0; num_var = 0;
		den_var += h->GetBinContent(i);
		for ( Int_t j = 0; j < bg_dim + 1; j++ ){
			den_var += TMath::Power( TMath::Power( x, j )*rp.bg_err[j], 2 );
			for ( Int_t k = j + 1; k < bg_dim + 1; k++ ){
				den_var += 2*TMath::Power( x, j + k )*TMatrixDRow( cov, j )( k );
			}
		}
		num_var = TMath::Power( x, 2 )*den_var;
		num_var_tot += num_var;
		den_var_tot += den_var;
	}
	
	// Assign the final values
	rp.a = cent_den;
	rp.c = cent_num/cent_den;
	rp.a_err = TMath::Sqrt( den_var_tot );
	rp.c_err = TMath::Sqrt( num_var_tot );
	return;
}

// --------------------------------------------------------------------------------------------- //
void PrintRedPeak( RedPeak_t &rp, std::ostream& f = std::cout ){
	TString empty_str = "";
	
	f << std::setprecision(print_precision) << 
	std::setw(print_width) << empty_str << "\t" << 
	std::setw(print_width) << empty_str << "\t" << 
	std::setw(print_width) << empty_str << "\t" << 
	std::setw(print_width) << empty_str << "\t" << 
	std::setw(print_width) << rp.lb << "\t" << 
	std::setw(print_width) << rp.ub << "\t" << 
	std::setw(print_width) << rp.a << "\t" << 
	std::setw(print_width) << rp.a_err << "\t" << 
	std::setw(print_width) << rp.c << "\t" << 
	std::setw(print_width) << rp.c_err << "\n";

	return;
}

// --------------------------------------------------------------------------------------------- //
void PrintFitPars( TFitResultPtr r, Int_t bg_dim, std::ostream& f = std::cout ){
	f << "[" << r->Chi2() << ",";			// Chi^2
	for ( Int_t i = 0; i < 2; i++ ){
		if ( i < bg_dim + 1 ){
			f << r->Parameter(i) << ",";	// Quadratic background parameters
		}
		else{
			f << "0.0,";
		}
	}
	f << "0.0,0.0,0.0]\n";					// The rest
	
	return;
}





























#endif
