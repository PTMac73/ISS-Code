// DoubletFitter.C
// Fits the ground state doublet with different l-values
// ============================================================================================= //
// Patrick MacGregor
// Nuclear Physics Research Group
// School of Physics and Astronomy
// The University of Manchester
// ============================================================================================= //

#include <TString.h>
#include <TF1.h>
#include <TMatrixD.h>

#include "Doublet_ImportData.h"
#include "Doublet_CreateFitFunc.h"
#include "Doublet_CalculateEtas.h"
#include "Doublet_Globals.h"


void DoubletFitter( TString in_file_loc = input_file_dir ){
	// Import the data for creating polynomial fits
	ImportData( in_file_loc );
	
	// Calculate the x range
	Double_t X_min = 10000.0; Double_t X_max = 0;
	FindMinMax( X, NUM_DATA_POINTS, X_min , X_max );

	// Generate the fit name
	TString fit_name = CreateFuncName(NUM_DATA_POINTS - 1);

	// Generate the fitting functions for experimental data
	TF1 *exp_fit = new TF1( "exp_fit", fit_name, GetRange( X_min, X_max, 1 ), GetRange( X_min, X_max, 2 ) );
	CreateFitFunc( X, Y, exp_fit );	

	// Generate the fitting functions for theoretical distributions
	TF1 *pt_fit[NUM_L];
	for ( Int_t i = 0; i < NUM_L; i++ ){
		pt_fit[i] = new TF1( Form( "pt_fit_%i", i ), fit_name, GetRange( X_min, X_max, 1 ), GetRange( X_min, X_max, 2 ) );
		CreateFitFunc( X, PT[i], pt_fit[i] );
	}
	
	// Now calculate the eta parameters for each set (eta[L1][L2][0 or 1 or 2 (for chi2)])
	Double_t Eta[NUM_L][NUM_L][3];

	// Now calculate the etas
	for ( Int_t i = 0; i < NUM_L; i++ ){
		for ( Int_t j = 0; j < NUM_L; j++ ){
			
			// Only fill if j < i so that you only test between two angular distributions once
			if ( i == 0 && j == 2 /*j < i*/ ){
				// Put in some initial estimates
				Eta[i][j][0] = 1.5;
				Eta[i][j][1] = 1.5;

				// Now calculate properly
				Eta[i][j][2] = CalculateEtas( X, E, exp_fit, pt_fit[i], pt_fit[j], Eta[i][j][0], Eta[i][j][1] );

				// Print the result
				if ( SWITCH_VERBOSE == 1 ){
					printf("%i --> %i:\t%8.8f\t%8.8f\t%8.8f\n", i, j, Eta[i][j][0], Eta[i][j][1], Eta[i][j][2] );
				}
			}
			else{
				// Fill eta with -1
				Eta[i][j][0] = -1.0;
				Eta[i][j][1] = -1.0;
				Eta[i][j][2] = -1.0;
			}
		}
	}
	
	// Print summary of useful information
		if ( SWITCH_VERBOSE == 1 ){
			printf("====================================================\n");
			printf(" L1 : L2 : Eta_1       : Eta_2       : Chi^2      \n");
			printf("----:----:-------------:-------------:--------------\n");
		}
		for ( Int_t i = 0; i < NUM_L; i++ ){
			for ( Int_t j = 0; j < NUM_L; j++ ){
				if ( i < j ){
					std::cout << std::setw(3) << i << " : " << std::setw(2) << j << " : " << std::setw(11) << std::fixed << std::setprecision(8) << Eta[i][j][0] << " : " << std::setw(11) << std::fixed << std::setprecision(8) << Eta[i][j][1] << " : " << std::setw(11) << std::fixed << std::setprecision(8) << Eta[i][j][2] << "\n";
					//printf("%i\t%i\t%8.8f\t%8.8f\t%8.8f\n", j, i, Eta[i][j][1], Eta[i][j][0], Eta[i][j][2] );
				}
			}
		} 
		if ( SWITCH_VERBOSE == 1 ){
			printf("====================================================\n");
		}
	

	// Loop over the different l-values and work out the "best chi^2"
	/*	
	Double_t F[3] = {1,3,5};
	Double_t F1[3] = {0,1,1};
	Double_t F2[3] = {1,0,2};
	Double_t X[3] = {-2,0,1};
	TF1 *F_fit = new TF1( "F_fit", "[0] + [1]*x + [2]*x^2", -5, 5 );
	TF1 *F1_fit = new TF1( "F1_fit", "[0] + [1]*x + [2]*x^2", -5, 5 );
	TF1 *F2_fit = new TF1( "F2_fit", "[0] + [1]*x + [2]*x^2", -5, 5 );

	for (Int_t i =0; i < 3; i++){
		F_fit->SetParameter(i, F[i] );
		F1_fit->SetParameter(i, F1[i] );
		F2_fit->SetParameter(i, F2[i] );
	}
	
	CalculateMinChi2( X, F_fit, F1_fit, F2_fit );
	
	// Return the best "chi^2" value with the angular momenta and the fitting functions and their
	// normalisations
	*/


}
































