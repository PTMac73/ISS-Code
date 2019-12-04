// Doublet_CreateFitFunc.h
// Creates a fitting function (Nth order polynomial) based on the X and Y values given
// ============================================================================================= //
// Patrick MacGregor
// Nuclear Physics Research Group
// School of Physics and Astronomy
// The University of Manchester
// ============================================================================================= //
// CONTAINS
//  * CreateFuncName
//  * FindMinMax
//  * GetRange
//  * CreateFitFunc
// ============================================================================================= //

#ifndef DOUBLET_CREATE_FIT_FUNC_H_
#define DOUBLET_CREATE_FIT_FUNC_H_

#include <TF1.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TMatrixD.h>
#include <TMath.h>
#include <TAxis.h>
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

// Find the minimum and maximum of an array
void FindMinMax( Double_t *arr, Int_t size_of_array, Double_t &min, Double_t &max ){
	for ( Int_t i = 0; i < size_of_array; i++ ){
		if ( arr[i] < min ){ min = arr[i]; }
		if ( arr[i] > max ){ max = arr[i]; }
	}
}

// Calculate the range for the functions and for plotting graphs
Double_t GetRange( Double_t k_min, Double_t k_max, Int_t range_switch = 0 ){
	Double_t fraction = 0.33;
	Double_t result;
	switch ( range_switch ){
		// Return the lower bound
		case 1:
			result = k_min - fraction*(k_max - k_min);
			break;

		// Return the upper bound
		case 2:
			result = k_max + fraction*(k_max - k_min);
			break;

		// Return the range
		default:
			result = k_max - k_min;
			break;

	}
	return result;
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

	// Find the minimum and maximum values for x and y
	Double_t x_min = 10000.0; Double_t x_max = 0;
	Double_t y_min = 10000.0; Double_t y_max = 0;
	FindMinMax( x, NUM_DATA_POINTS, x_min , x_max );
	FindMinMax( y, NUM_DATA_POINTS, y_min, y_max );
	
	// Define the fit parameters
	for ( Int_t i = 0; i < NUM_DATA_POINTS; i++ ){
		fit->SetParameter(i, fit_par(i, 0) );
	}

	// Draw the fit
	if ( SWITCH_DRAW_FITS == 1 ){
		// Define a TCanvas
		TCanvas *c = new TCanvas();

		// Define a graph for the points, and format it
		TGraph *graph = new TGraph( NUM_DATA_POINTS, x, y );
		graph->SetMarkerStyle(20);
		graph->SetMarkerSize(1);

		// Format the fit
		fit->SetLineWidth(2);

		// Draw the fits
		fit->Draw("C");
		graph->Draw("P SAME");

		// Edit the view of the axes
		fit->GetXaxis()->SetRangeUser( GetRange( x_min, x_max, 1 ), GetRange( x_min, x_max, 2 ) );
		fit->GetYaxis()->SetRangeUser( GetRange( y_min, y_max, 1 ), GetRange( y_min, y_max, 2 ) );

		printf( "DRAWN A FIT\n" );

	}
}

#endif
