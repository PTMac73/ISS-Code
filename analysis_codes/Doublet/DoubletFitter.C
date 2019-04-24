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

#include "Doublet_ImportData.h"
#include "Doublet_CreateFitFunc.h"
#include "Doublet_Globals.h"


void DoubletFitter(){
	// Import the data for creating polynomial fits
	ImportData();

	// Generate the fitting functions
	TF1 *exp_fit = new TF1();
	TF1 *pt_fit[NUM_L];

	CreateFitFunc( X, Y, exp_fit );
	/*
	for ( Int_t i = 0; i < NUM_L; i++ ){
		CreateFitFunc( X, PT[i], pt_fit[i] );
	}
	*/

	// Loop over the different l-values and work out the "best chi^2"

	// Return the best "chi^2" value with the angular momenta and the fitting functions and their
	// normalisations


}
