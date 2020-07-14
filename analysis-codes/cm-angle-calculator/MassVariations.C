// MassVariations.C
// Visualise the effect of changing the mass
// ============================================================================================= //
// Patrick MacGregor
// Nuclear Physics Research Group
// School of Physics and Astronomy
// The University of Manchester
// ============================================================================================= //

#include <TCanvas.h>
#include <TGraph.h>
#include <TLine.h>
#include <TMath.h>
#include <TText.h>

#include <iostream>

#include "CMAngleCalculator.h"
#include "MassVariations.h"

// Define the mass type
enum MassType_t { kM1, kM2, kM3, kM4 };

// Set the mass to a specific type
void SetMass( ReactionParameters& rp, MassType_t mt, Double_t mass ){
	switch(mt){
		case kM1 : rp.m1 = mass; break;
		case kM2 : rp.m2 = mass; break;
		case kM3 : rp.m3 = mass; break;
		case kM4 : rp.m4 = mass; break;
	}
}

// ---===---===---===---===---===---===---===---===---===---===---===---===---===---===---===--- //
// Calculates the CM Angle for a given mass
Double_t CMAngleCalculateMass( Double_t ex, Double_t z, MassType_t mt, Double_t mass ){
	// Define initial fixed quantities for the reaction
	ReactionParameters rp;
	FillRP(rp);
	SetMass(rp, mt, mass);

	// Now calculate derived quantities from fixed quantities
	KinematicsParameters kp;
	CalculateKP( rp, kp, ex );

	// Function-specific calculations
	Double_t p_para = rp.q*rp.B*z/( 2*TMath::Pi() );											// Parallel component of the momentum in MeV / c [LAB] -------> THIS IS A GUESS!!
	Double_t p_para_cm = p_para/kp.gamLAB_CM - kp.beta*kp.e3_cm;								// Parallel component of the momentum in MeV / c [CM]
	Double_t p_perp_cm = TMath::Sqrt( kp.e3_cm*kp.e3_cm - p_para_cm*p_para_cm - rp.m3*rp.m3 );	// Perpendicular component of the momentum in MeV / c [CM]
	

	// Calculate the function to be minimised and its first derivative
	Double_t fp = ( 2*p_perp_cm/( rp.q*rp.B ) )*TMath::Sin( rp.q*rp.B*z/( 2*kp.gamLAB_CM*( p_para_cm + kp.beta*kp.e3_cm ) ) ) - kp.rho;
	Double_t fd = ( ( -2*p_para_cm )/( rp.q*rp.B*p_perp_cm ) )*TMath::Sin( rp.q*rp.B*z/( 2*kp.gamLAB_CM*( p_para_cm + kp.beta*kp.e3_cm ) ) ) - ( ( p_perp_cm*z )/( kp.gamLAB_CM*TMath::Power( p_para_cm + kp.beta*kp.e3_cm , 2 ) ) )*TMath::Cos( rp.q*rp.B*z/( 2*kp.gamLAB_CM*( p_para_cm + kp.beta*kp.e3_cm ) ) );

	// Now perform Newton-Raphson calculation to get the angle in the CM frame
	Int_t n = 0;
	while ( TMath::Abs(fp) > 1e-5 && n < 100000){
		p_para_cm = p_para_cm - fp/fd;
		p_perp_cm = TMath::Sqrt( kp.e3_cm*kp.e3_cm - p_para_cm*p_para_cm - rp.m3*rp.m3 );

		fp = ( 2*p_perp_cm/( rp.q*rp.B ) )*TMath::Sin( rp.q*rp.B*z/( 2*kp.gamLAB_CM*( p_para_cm + kp.beta*kp.e3_cm ) ) ) - kp.rho;
		fd = ( ( -2*p_para_cm )/( rp.q*rp.B*p_perp_cm ) )*TMath::Sin( rp.q*rp.B*z/( 2*kp.gamLAB_CM*( p_para_cm + kp.beta*kp.e3_cm ) ) ) - ( ( p_perp_cm*z )/( kp.gamLAB_CM*TMath::Power( p_para_cm + kp.beta*kp.e3_cm , 2 ) ) )*TMath::Cos( rp.q*rp.B*z/( 2*kp.gamLAB_CM*( p_para_cm + kp.beta*kp.e3_cm ) ) );
		n++;
	}

	// Function minimised - now calculate the final quantities
	Double_t theta_cm = 180 - ( TMath::ACos( p_para_cm/kp.p3_cm )*180/TMath::Pi() );
	Double_t e3 = kp.gamLAB_CM*( kp.e3_cm + kp.beta*p_para_cm );
	Double_t p3 = TMath::Sqrt( e3*e3 - rp.m3*rp.m3 );
	Double_t p_perp = p_perp_cm;
	p_para = TMath::Sqrt( p3*p3 - p_perp*p_perp );

	Double_t theta = TMath::ACos( kp.gamLAB_CM*( p_para_cm + kp.beta*kp.e3_cm )/p3 )*180/TMath::Pi();

	return theta_cm;

}


// ---===---===---===---===---===---===---===---===---===---===---===---===---===---===---===--- //
// Batch plot the mass variations
void MassVariations(){
	// Declare variables
	MassType_t mass_type[NUM_MASS] = { kM1, kM2, kM4 };		// Define the different mass types
	
	const Int_t NUM_POINTS = 4000;						// Number of points in a given graph
	const Double_t STEP = 0.01;					// Step size in MeV
	
	TGraph* g[NUM_MASS];
	TGraph* gd[NUM_MASS];
	TCanvas* c[NUM_MASS];
	TCanvas* cd[NUM_MASS];
	
	const Int_t NUM_LINE = 13;
	TLine* l[NUM_LINE];
	TText* text[NUM_LINE];
	
	TLine* ld[NUM_LINE];
	TText* textd[NUM_LINE];
	
	
	// *LOOP* over mass (i)
	for ( Int_t i = 0; i < 3; i++ ){
	
		Double_t* theta_cm = new Double_t[NUM_POINTS];		// Define an array to hold the angles
		Double_t* theta_deriv = new Double_t[NUM_POINTS];		// Define an array to hold the angles
		Double_t* mass = new Double_t[NUM_POINTS];			// Define an array to hold the masses
	
		// Set mass and theta values to zero
		for ( Int_t j = 0; j < NUM_POINTS; j++ ){
			mass[j] = 0.0;
			theta_cm[j] = 0.0;
			theta_deriv[j] = 0.0;
		}
		
		
		// CALCULATE THE PARAMETERS	
		// *LOOP* over number of points (j)
		for ( Int_t j = 0; j < NUM_POINTS; j++ ){
			
			// Calculate the mass (pivot around base)
			mass[j] = mass_base[i] + STEP*(j - 0.5*NUM_POINTS);
			
			// Calculate theta
			Double_t temp_double = CMAngleCalculateMass( 0.0, -40.0, mass_type[i], mass[j] );
			if ( TMath::IsNaN(temp_double) ){ theta_cm[j] = 0.0; }
			else{ theta_cm[j] = temp_double; }
			
			// Calculate derivative of line
			if ( j > 2 ){
				theta_deriv[j] = (1/STEP)*( theta_cm[j] - theta_cm[j-1] );
			}
			else{
				theta_deriv[j] = 0.0;
			}
			
		}
		
		// Create TLine arrays
		for ( Int_t j = 0; j < NUM_LINE; j++ ){
			l[j] = new TLine( mass_base[i] - j*0.511, 0, mass_base[i] - j*0.511, 50 );
			l[j]->SetLineStyle(2);
			l[j]->SetLineWidth(1);
			text[j] = new TText( mass_base[i] - j*0.511, 52, Form( "%ie", j ) );
			text[j]->SetTextAlign(22);
			text[j]->SetTextFont(42);
			text[j]->SetTextSize(0.03);
			
			ld[j] = new TLine( mass_base[i] - j*0.511, -8, mass_base[i] - j*0.511, 4 );
			ld[j]->SetLineStyle(2);
			ld[j]->SetLineWidth(1);
			textd[j] = new TText( mass_base[i] - j*0.511, 5, Form( "%ie", j ) );
			textd[j]->SetTextAlign(22);
			textd[j]->SetTextFont(42);
			textd[j]->SetTextSize(0.03);
		}
		
		
		
		
		// GRAPH THE RESULTS
		// Theta_cm graph
		g[i] = new TGraph( NUM_POINTS, mass, theta_cm );
		CreateMassGraph( g[i], i );
		
		c[i] = new TCanvas( Form( "c%i", i + 1 + ( i == 2 ? 1 : 0 ) ), "CANVAS", C_WIDTH, C_HEIGHT );
		GlobSetCanvasMargins( c[i] );
		
		g[i]->Draw("APL");
		
		for ( Int_t j = 0; j < NUM_LINE; j++ ){
			if ( i == 0 && j < 10 ){ l[j]->Draw("SAME"); text[j]->Draw("SAME"); }
			if ( i == 1 && j < 2 ){ l[j]->Draw("SAME"); text[j]->Draw("SAME"); }
			if ( i == 2 && j < NUM_LINE ){ l[j]->Draw("SAME"); text[j]->Draw("SAME"); }
		}
		
		c[i]->Print( Form( "mass_variation_%i.svg", i + 1 + ( i == 2 ? 1 : 0 ) ) );
		
		// Derivative graph
		gd[i] = new TGraph( NUM_POINTS, mass, theta_deriv );
		CreateMassDerivGraph( gd[i], i );
		
		cd[i] = new TCanvas( Form( "cd%i", i + 1 + ( i == 2 ? 1 : 0 ) ), "CANVAS", C_WIDTH, C_HEIGHT );
		GlobSetCanvasMargins( cd[i], 0.15 );
		
		gd[i]->Draw("APL");
		
		for ( Int_t j = 0; j < NUM_LINE; j++ ){
			if ( i == 0 && j < 10 ){ ld[j]->Draw("SAME"); textd[j]->Draw("SAME"); }
			if ( i == 1 && j < 2 ){ ld[j]->Draw("SAME"); textd[j]->Draw("SAME"); }
			if ( i == 2 && j < NUM_LINE ){ ld[j]->Draw("SAME"); textd[j]->Draw("SAME"); }
		}
		
		cd[i]->Print( Form( "mass_variation_deriv_%i.svg", i + 1 + ( i == 2 ? 1 : 0 ) ) );
		
		// Empty the arrays
		delete[] mass;
		delete[] theta_cm;
		delete[] theta_deriv;
		
	}
	
	return;
}





































