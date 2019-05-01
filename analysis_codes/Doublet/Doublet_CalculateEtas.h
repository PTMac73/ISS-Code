// Doublet_CalculateEtas.h
// Calculates the etas based on the minimum chi^2 values being used
// ============================================================================================= //
// Patrick MacGregor
// Nuclear Physics Research Group
// School of Physics and Astronomy
// The University of Manchester
// ============================================================================================= //
#ifndef DOUBLET_CALCULATE_ETAS_H_
#define DOUBLET_CALCULATE_ETAS_H_

#include <TMatrixD.h>
#include <TMath.h>
#include <TGraph2D.h>
#include <TStyle.h>
#include "Doublet_Eta_Functions.h"
#include "Doublet_Globals.h"

Double_t CalculateEtas( Double_t *x, Double_t *error, TF1 *beta_exp, TF1 *beta_1, TF1 *beta_2, Double_t &eta_0, Double_t &eta_1){
	// Define preliminary quantities
	TMatrixD x_matrix(NUM_DATA_POINTS, NUM_DATA_POINTS);
	TMatrixD error_matrix(NUM_DATA_POINTS, NUM_DATA_POINTS);
	TMatrixD beta_1_vector(NUM_DATA_POINTS,1);
	TMatrixD beta_2_vector(NUM_DATA_POINTS,1);
	TMatrixD beta_exp_vector(NUM_DATA_POINTS,1);
	TMatrixD mu_var (NUM_DATA_POINTS, 1 );

	// Populate these matrices
	for ( Int_t i = 0; i < NUM_DATA_POINTS; i++ ){
		for ( Int_t j = 0; j < NUM_DATA_POINTS; j++ ){
			// Populate the x matrix
			x_matrix(i, j) = TMath::Power( x[i], j );

			// Populate the error matrix
			if ( i == j ){
				error_matrix(i, j) = 1.0/( error[i]*error[i] );
			}
			else{
				error_matrix(i, j) = 0;
			}
		}
		
		
		

		// Populate the beta matrices
		beta_1_vector(i, 0) = beta_1->GetParameter(i);
		beta_2_vector(i, 0) = beta_2->GetParameter(i);
		beta_exp_vector(i, 0) = beta_exp->GetParameter(i);
	}
	
	// Generate matrix_array
	TMatrixD matrix_array[5] = { x_matrix, error_matrix, beta_1_vector, beta_2_vector, beta_exp_vector };
	
	// CALCULATE MINIMUM CHI^2 VALUE
	if ( SWITCH_BRUTE_FORCE == 1 ){
		// Define initial variables
		Double_t step_size = 0.1;
		TMatrixD mu = CalculateMu( beta_exp_vector, beta_1_vector, beta_2_vector, eta_0, eta_1 );
		Double_t chi2 = CalculateChi2( mu, x_matrix, error_matrix );
		Double_t chi2_old = 1e10;

		// Zoom in on the minimum
		while ( step_size > 1e-8 ){
			chi2 = BruteForce( step_size, eta_0, eta_1, chi2, matrix_array );
			if ( chi2 < chi2_old ){ chi2_old = chi2;}
			else{ chi2 = chi2_old; }
			step_size /= 10;
		}
		// Print the best
		/*
		printf("\n\nBEST PARAMETERS:\n");
		printf("ETA-0\t\tETA-1\t\tCHI^2\n");
		printf("%8.8f\t%8.8f\t%8.8f\n", eta_0, eta_1, chi2 );
		*/

		return chi2;
	}
}

/*
		// Perform a hill climb
		// Lay out grid of starting positions in steps of 100
		Double_t initial_conditions[9][2] = {
			{0,0},
			{100,0},
			{200,0},
			{0,100},
			{100,100},
			{100,200},
			{200,0},
			{200,100},
			{200,200}
		};
		Double_t eta_array[9][3];
		Double_t eta_0, eta_1, chi2;
		Double_t best_eta_0, best_eta_1, best_chi2 = 0;

		// Loop over initial conditions
		for ( Int_t k = 0; k < 9; k++ ){

			// Define the starting eta's
			eta_0 = initial_conditions[k][0];
			eta_1 = initial_conditions[k][1];

			// Get the best chi^2 from the hill-climb
			printf("================================================================================\n");
			printf("STARTING FROM ( %f, %f )\n", eta_0, eta_1 );
			printf("================================================================================\n");
			chi2 = HillClimb( eta_0, eta_1, matrix_array );

			// Store the eta's
			eta_array[k][0] = eta_0;
			eta_array[k][1] = eta_1;
			eta_array[k][2] = chi2;
			
			// Print stuff
			printf(">>>\t%8.4f\t%8.4f\t%8.4f\n", eta_array[k][0], eta_array[k][1],  eta_array[k][2] );

			// Calculate the best chi^2
			if ( best_chi2 == 0 || chi2 < best_chi2 ){
				best_chi2 = chi2;
				best_eta_0 = eta_0;
				best_eta_1 = eta_1;
			}
		}
		// Print the best
		printf("\n\nBEST PARAMETERS:\n");
		printf("ETA-0\t\tETA-1\t\tCHI^2\n");
		printf("%8.4f\t%8.4f\t%8.4f\n", best_eta_0, best_eta_1, best_chi2 );
	
	*/


#endif
