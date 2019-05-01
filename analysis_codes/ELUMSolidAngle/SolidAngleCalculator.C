// SolidAngleCalculator.C
// Calculates the solid angle for the luminosity detector, by solving a transcendental equation
// ============================================================================================= //
// Patrick MacGregor
// Nuclear Physics Research Group
// School of Physics and Astronomy
// The University of Manchester
// ============================================================================================= //
#include <TMath.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <iostream>
#include "conversion.h"

// Define variables
Double_t q = 1;						// Charge of deuteron [e]
Double_t B = 2.5;					// Magnetic Field [T]
Double_t mass_excess = 13.1357; 	// Mass excess of deuteron [MeV] (from NNDC)
Double_t A = 2;						// Mass number of deuteron [#]
Double_t mass = MassExcessToMass( 2, mass_excess );		// Mass of deuteron [u]
Double_t energy_per_u = 9.473;		// MeV per u	

Double_t z_elum = 168.4; 			// Distance from target to elum [mm]
Double_t z_elum_error = 0.7;		// Error on distance from target to elum [mm]
Double_t Al_distance = 10.0;		// Distance from elum to Al shield [mm]
Double_t radius_in = 48.00;			// Inner radius of elum [mm]
Double_t radius_out = 96.00;		// Outer radius of elum [mm]
Double_t velocity = TMath::Sqrt( 2*A*energy_per_u/( mass*AMU )  ); // Velocity of deuteron [1/c]

// Define radius function [ theta in rad, z in mm ]
Double_t CalculateRadius( Double_t theta, Double_t z ){
	Double_t amplitude = MToMM()*velocity*TMath::C()*mass*AMUToKg()*TMath::Sin( theta )/( q*EToCoulombs()*B );
	Double_t angle = q*EToCoulombs()*B*z/( MToMM()*velocity*TMath::C()*mass*AMUToKg()*TMath::Cos( theta ) );
	return amplitude*( 1.0 - TMath::Cos( angle ) );
}

// MAIN FUNCTION ------------------------------------------------------------------------------- //
void SolidAngleCalculator(){
	// Draw the function for a given theta
	Int_t size_of_array = 1685;
	Int_t size_of_theta = 180;
	Double_t theta[size_of_array*size_of_theta];
	Double_t z[size_of_array*size_of_theta];
	Double_t r[size_of_array*size_of_theta];

	for (Int_t i = 0; i < size_of_array; i++ ){
		for (Int_t j = 0; j < size_of_theta; j++ ){
			theta[size_of_theta*i + j] = 0.5*j;
			z[size_of_theta*i + j] = 0.1*i;
			r[size_of_theta*i + j] = CalculateRadius( theta[size_of_theta*i + j], z[size_of_theta*i + j] );
			if ( j % size_of_theta/10 == 0 && i % 20 == 0 ){
				printf("%f\t%f\t%f\n", z[i], theta[j], r[size_of_theta*i + j] );
			}
		}
	}
	gStyle->SetPalette( 55 );
	TGraph2D *graph2 = new TGraph2D( size_of_array*size_of_theta, z, theta, r );
	printf("Banana\n");
	TCanvas *c_graph2 = new TCanvas( "c_graph2", "GRAPH2", 1200, 900 );
	
	graph2->SetTitle("GRAPH; z; theta; r");
	graph2->Draw("surf");
	
	
	
}

