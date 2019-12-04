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
		Double_t step_size = STEP_SIZE_I;
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
		return chi2;
	}
}


#endif
