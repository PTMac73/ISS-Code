// Doublet_CreateFitFunc.h
// Creates a fitting function (Nth order polynomial) based on the X and Y values given
// ============================================================================================= //
// Patrick MacGregor
// Nuclear Physics Research Group
// School of Physics and Astronomy
// The University of Manchester
// ============================================================================================= //
#ifndef DOUBLET_CREATE_FIT_FUNC_H_
#define DOUBLET_CREATE_FIT_FUNC_H_

#include <TF1.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TMatrixD.h>
#include <TMath.h>
#include "Doublet_Globals.h"

// Create the polynomial name
TString CreateFuncName( Int_t N ){
	TString fit_name = "[0] + ";
	for ( Int_t i = 1; i <= N; i++ ){
		fit_name += Form( "[%i]*x^%i", i, i );
		if ( i != N ){
			fit_name += " + ";
		}
	}
	return fit_name;
}




void CreateFitFunc( Double_t *x, Double_t *y, TF1 *fit ){
	
	// Define the vector of Y and the matrix of X
	TMatrixD y_matrix(NUM_DATA_POINTS, 1);
	TMatrixD x_matrix(NUM_DATA_POINTS, NUM_DATA_POINTS);
	
	// Populate the matrices
	for ( Int_t i = 0; i < NUM_DATA_POINTS; i++ ){
		for ( Int_t j = 0; j < NUM_DATA_POINTS; j++ ){
			x_matrix(i, j) = TMath::Power( x[i], j );
		}
		// Populate the y_matrix
		y_matrix(i, 0) = y[i];
	}

	// Calculate the fit parameters
	TMatrixD fit_par = x_matrix.Invert()*y_matrix;

	// Construct the fit
	TString fit_name = CreateFuncName(NUM_DATA_POINTS - 1);
	fit = new TF1( "fit", fit_name );
	for ( Int_t i = 0; i < NUM_DATA_POINTS; i++ ){
		fit->SetParameter(i, fit_par(i, 0) );
	}

	// Draw the fit
	TGraph *graph = new TGraph( NUM_DATA_POINTS, x, y );
	TCanvas *c = new TCanvas();
	graph->Draw("AP");
	graph->SetMarkerStyle(3);
	fit->SetLineWidth(2);
	fit->Draw("same");
	fit_par.Print();

// Construct a TGraph
	/*TGraph *graph = new TGraph( NUM_DATA_POINTS, x, y );
	
	TCanvas *c = new TCanvas();
	graph->Draw("AP");

	// Fit the TGraph
	
	//std::cout << fit_name << std::endl;
	fit = new TF1( "fit", fit_name );
	graph->Fit( fit->GetName() );*/
	



}






#endif
