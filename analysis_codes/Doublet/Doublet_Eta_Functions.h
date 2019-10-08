// Doublet_Eta_Functions.h
// Additional functions for calculating the etas
// ============================================================================================= //
// Patrick MacGregor
// Nuclear Physics Research Group
// School of Physics and Astronomy
// The University of Manchester
// ============================================================================================= //
#ifndef DOUBLET_ETA_FUNCTIONS_H_
#define DOUBLET_ETA_FUNCTIONS_H_

#include <TMatrixD.h>
#include <TMath.h>
#include <TGraph2D.h>
#include <TStyle.h>
#include <TStopwatch.h>

#include "Doublet_Globals.h"
/*
// Calculates the chi^2 value ------------------------------------------------------------------ //
Double_t CalculateChi2( TMatrixD &mu, TMatrixD &x ){
	TMatrixD Chi2(1,1);
	TMatrixD mu_t( mu.GetNcols(), mu.GetNrows() );
	TMatrixD x_t( x.GetNcols(), x.GetNrows() );
	mu_t.Transpose( mu );
	x_t.Transpose( x );
	Chi2 = mu_t*x_t*x*mu;
	return Chi2(0,0);
}
*/

// Calculates the chi^2 value ------------------------------------------------------------------ //
Double_t CalculateChi2( TMatrixD &mu, TMatrixD &x, TMatrixD &W ){
	TMatrixD Chi2(1,1);
	TMatrixD mu_t( mu.GetNcols(), mu.GetNrows() );
	TMatrixD x_t( x.GetNcols(), x.GetNrows() );
	mu_t.Transpose( mu );
	x_t.Transpose( x );
	Chi2 = mu_t*x_t*W*x*mu;
	return Chi2(0,0);
}

// Calculates the mu vector -------------------------------------------------------------------- //
TMatrixD CalculateMu( TMatrixD &beta_exp, TMatrixD &beta_0, TMatrixD &beta_1, Double_t eta_0, Double_t eta_1 ){
	TMatrixD mu ( beta_exp.GetNrows(), beta_exp.GetNcols() );
	mu = beta_exp - eta_0*beta_0 - eta_1*beta_1;
	return mu;
}

// Brute forces the solution and returns the best ---------------------------------------------- //
Double_t BruteForce( Double_t spacing, Double_t &eta_0, Double_t &eta_1, Double_t chi2, TMatrixD *matrix_array ){
	// Decompose matrix array
	TMatrixD x_matrix = matrix_array[0];
	TMatrixD error_matrix = matrix_array[1];
	TMatrixD beta_1_vector = matrix_array[2];
	TMatrixD beta_2_vector = matrix_array[3];
	TMatrixD beta_exp_vector = matrix_array[4];

	// Define initial quantities
	Double_t diff = 1.5*1e4*spacing;		// Defines the range over which you search
	Double_t lb_0 = 0.0, ub_0 = 0.0, lb_1 = 0.0, ub_1 = 0.0;

	if ( diff < eta_0 ){ lb_0 = eta_0 - diff; }	// Define the upper and lower bounds for each eta
	ub_0 = eta_0 + diff;
	if ( diff < eta_1 ){ lb_1 = eta_1 - diff; }
	ub_1 = eta_1 + diff;
	
	Double_t chi2_min = chi2, eta_0_min = eta_0, eta_1_min = eta_1;	// Define initial minima
	Int_t num_points_ctr = 0;										// Define a counter
	TMatrixD mu = CalculateMu( beta_exp_vector, beta_1_vector, beta_2_vector, eta_0, eta_1 );	// Define initial mu
	
	// Begin loop over each eta
	TStopwatch stopwatch;
	stopwatch.Start();
	Int_t divisor = (Int_t)( 4/(20*spacing*spacing) );


	for (Int_t i = 0; i < (Int_t)((ub_0 - lb_0)/spacing); i++ ){
		for (Int_t j = 0; j < (Int_t)((ub_1 - lb_1)/spacing); j++ ){
			
			// Calculate mu
			mu = CalculateMu( beta_exp_vector, beta_1_vector, beta_2_vector, lb_0 + i*spacing, lb_1 + j*spacing );
			
			// Calculate the chi^2
			chi2 = CalculateChi2( mu, x_matrix, error_matrix );

			// Iterate the counter for the number of points
			num_points_ctr++;

			if ( SWITCH_VERBOSE == 1 && num_points_ctr % divisor == 0){
				std::cout << std::setprecision(4) << num_points_ctr*spacing*spacing/4 << "% complete at " << spacing << " scale\n";
				std::cout << "Time elapsed: " << std::setprecision(4) << stopwatch.RealTime() << " s\n\n";
				stopwatch.Start(kFALSE);
			}
		
			// Store the new minimum
			if ( chi2 < chi2_min ){
				chi2_min = chi2;
				eta_0_min = lb_0 + i*spacing;
				eta_1_min = lb_1 + j*spacing;
			}
		}
	}
	stopwatch.Stop();
	eta_0 = eta_0_min;
	eta_1 = eta_1_min;
	
	// Print the final value
	if ( SWITCH_VERBOSE == 1 ){
		printf("%2.1e\t( %08.8f --> %08.8f )\t%08.8f\t( %08.8f --> %08.8f )\t%08.8f\t%8.4f\n", spacing, lb_0, ub_0, eta_0, lb_1, ub_1, eta_1, chi2_min );
	}
	// Return the final value
	return chi2_min;
}



// Calculates the best chi^2 using a hill-climbing algorithm ----------------------------------- //
Double_t HillClimb( Double_t &eta_0, Double_t &eta_1, TMatrixD *matrix_array ){
	// Decompose matrix array
	TMatrixD x_matrix = matrix_array[0];
	TMatrixD error_matrix = matrix_array[1];
	TMatrixD beta_1_vector = matrix_array[2];
	TMatrixD beta_2_vector = matrix_array[3];
	TMatrixD beta_exp_vector = matrix_array[4];
	
	// Define global constants for the hill-climb
	Double_t tolerance = 1e-10;			// Difference allowed between two chi^2 values before exiting
	Double_t step_tolerance = 1e-6;		// Step size minimum limit
	Double_t step_size = 0.1;			// Iniital step size to consider
	Int_t num_steps = 0;				// Counter for the number of steps taken

	// Define an initial mu and chi^2
	TMatrixD mu = CalculateMu( beta_exp_vector, beta_1_vector, beta_2_vector, eta_0, eta_1 );
	Double_t chi2 = CalculateChi2( mu, x_matrix, error_matrix );

	// Initialise other constants
	Double_t best_chi2 = chi2;		// Stores the minimum chi^2
	Double_t chi2_old = 0;			// Stores the old chi^2
	Int_t best_i = 0, best_j = 0;	// Stores the best step to take

	// Begin while loop
	while ( TMath::Abs( chi2 - chi2_old ) > tolerance || step_size > step_tolerance ){
		// Store the old chi2 values
		chi2_old = chi2;

		// Calculate the new chi2 values
		for ( Int_t i = -1; i < 2; i++ ){
			for ( Int_t j = -1; j < 2; j++ ){

				// Don't stay in the same place
				if ( i != 0 && j != 0 ){
					
					// Calculate mu and chi^2
					mu = CalculateMu( beta_exp_vector, beta_1_vector, beta_2_vector, eta_0 + i*step_size, eta_1 + j*step_size );
					chi2 = CalculateChi2( mu, x_matrix, error_matrix );

					// Check if it's a minimum, and that new eta's are greater than zero
					if ( chi2 < best_chi2 && eta_0 + i*step_size >= 0.0 && eta_1 + j*step_size >= 0.0 ){
						best_i = i;
						best_j = j;
						best_chi2 = chi2;

					}
				}
			}
		}
	
		// Store the best chi2
		chi2 = best_chi2;
		eta_0 = eta_0 + best_i*step_size;
		eta_1 = eta_1 + best_j*step_size;

		// Reduce the step size to home in
		if ( best_chi2 == chi2_old && step_size > step_tolerance ){
			step_size /= 10;
			//printf("%i\tSTEP SIZE NOW %e\n", num_steps, step_size );
		}
		
		// Iterate the number of steps
		if ( num_steps % 20000 == 0 && SWITCH_VERBOSE == 1){
			printf("%0000000i:\t%12.8f\t%12.8f\t%12.8f\n", num_steps, eta_0, eta_1, chi2 );
		}
		num_steps++;
	}
	if ( SWITCH_VERBOSE == 1 ){
		printf(">>>\t%8.4f\t%8.4f\t%8.4f\n", eta_0, eta_1, chi2 );
	}

	// BRUTE FORCE IT NOW IN RIGHT REGION
	// Now brute force the way there
	step_size = 1;
	chi2_old = chi2;
	while ( step_size > 1e-6 ){
		chi2 = BruteForce( step_size, eta_0, eta_1, chi2, matrix_array );
		if ( chi2 < chi2_old ){ chi2_old = chi2;}
		else{ chi2 = chi2_old; }
		step_size /= 10;
	}
	
	
	

	return chi2;
}


#endif
