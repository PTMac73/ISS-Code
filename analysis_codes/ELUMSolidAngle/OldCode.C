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
#include <TMultiGraph.h>
#include <TLine.h>
#include <TBox.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TString.h>
#include <TH1.h>
#include <iostream>
#include <vector>
#include "conversion.h"

// DEFINE GLOBALS
Bool_t SWITCH_DRAW_CANVAS = 1;
Bool_t SWITCH_PRINT_CANVAS = 0;


// Define radius function [ mm ] - > units are [ [s^-1], [mm], [1/c], [1/c] ]
Double_t CalculateRadius( Double_t cyclotron_freq, Double_t z, Double_t v_para, Double_t v_perp ){
	return ( v_perp*TMath::C()*MToMM()/cyclotron_freq )*( 1.0 - TMath::Cos( cyclotron_freq*z/( v_para*TMath::C()*MToMM() ) ) );
}

Double_t CalculateSolidAngle( Double_t theta1, Double_t theta2 ){
	return 2*TMath::Pi()*( TMath::Cos( theta1 ) - TMath::Cos( theta2 ) );
}

// MAIN FUNCTION =============================================================================== //
void SolidAngleCalculator(){
	// OPTIONS
	if ( SWITCH_DRAW_CANVAS == 0 ){
		  gROOT->SetBatch(kTRUE);
	}


	// CONSTANTS OF MOTION --------------------------------------------------------------------- //
	// M1 = 28Mg, M2 = 2H, M3 = 2H, M4 = 28Mg
	Double_t mass_excess[4] = { -15.018845, 13.13572176, 13.13572176, -15.018845};	// [MeV]
	Int_t mass_numbers[4] = { 28, 2, 2, 28 };										// Just a number
	Double_t mass[4];																// [AMU]
	for ( Int_t k = 0; k < 4; k++ ){
		mass[k] = MassExcessToMass( mass_numbers[k], mass_excess[k] );				
	}
	vector <Int_t> num_revolutions = {1};														// Number of turns
	Double_t ejectile_charge = 1;													// [e]
	Double_t b_field = 2.5;															// [T]
	Double_t beam_energy_per_nucleon = 9.473;										// MeV per u
	Double_t excitation_energy = 0.0;
	Double_t cyclotron_freq = ejectile_charge*EToCoulombs()*b_field/( mass[2]*AMUToKg() );	//[s^-1]

	// Define distances
	Double_t z_elum = 168.4; 			// Distance from target to elum [mm]
	Double_t z_elum_error = 0.7;		// Error on distance from target to elum [mm]
	Double_t Al_distance = 10.0;		// Distance from elum to Al shield [mm]
	Double_t Al_thickness = 10.0;		// Thickness of the Al shield [mm]
	Double_t radius_in = 48.00;			// Inner radius of elum [mm]
	Double_t radius_out = 96.00;		// Outer radius of elum [mm]


	// Calculate the range of angles (i.e. number of indices in array)
	Double_t angle_divisions = 0.01;
	Int_t lb_angle = -180;
	Int_t ub_angle = 180;
	const Int_t range_angles = TMath::Ceil( ( ub_angle - lb_angle )/angle_divisions ) + 1;

	// Constant kinetic terms (with angle)
	Double_t Q = ( ( mass[0] + mass[1] ) - ( mass[2] + mass[3] ) )*AMU;					// Nuclear Q-value [MeV]
	Double_t T1 = beam_energy_per_nucleon*mass_numbers[0];								// Beam energy [MeV]
	Double_t T_cm_i = ( mass[1]/( mass[0] + mass[1] ) )*T1;								// Initial CM energy [MeV]
	Double_t V_cm = TMath::Sqrt( 2*T1*mass[0]*AMU )/( AMU*( mass[0] + mass[1] ) );		// CM velocity [1/c]						
	Double_t T_cm_f = T_cm_i + Q - excitation_energy;									// Final energy in CM [MeV]
	Double_t T_cm_3 = mass[3]*T_cm_f/( mass[2] + mass[3] );								// Final energy of ejectile in CM [MeV]
	Double_t v3 = TMath::Sqrt( 2*T_cm_3/(mass[2]*AMU) );								// Final velocity of ejectile in CM [MeV]


	// CALCULATE ANGLE-DEPENDENT QUANTITIES ---------------------------------------------------- //
	// Define the quantites to be calculated in loop over angles
	Double_t phi_cm[range_angles];				// Angle of scattering (intuitive) [rad]
	Double_t theta_cm[range_angles];			// Angle of scattering in CM [rad]
	Double_t theta_lab[range_angles];			// Angle of scattering in LAB
	Double_t T_lab_3[range_angles];				// Lab energy of ejectile in CM [MeV]
	Double_t v_para[range_angles];				// Parallel velocity of ejectile in LAB [1/c]
	Double_t v_perp[range_angles];				// Perpendicular velocity of ejectile in LAB [1/c]
	Double_t r_shield[range_angles][2];			// Distance from beam axis at shield (both edges) [mm]
	Double_t r_elum[range_angles];				// Distance from beam axis at elum [mm]
	Int_t success_counter = 0;					// Counts the number of states that satisfy the requirements
	Int_t index_array[range_angles];			// Records the indices of all the successful runs


	// LOOP OVER THETA_CM
	for ( Int_t k = 0; k < range_angles; k++ ){
		theta_cm[k] = ( lb_angle + k*angle_divisions )*TMath::DegToRad();									// [rad]
		T_lab_3[k] = T_cm_3 + 0.5*mass[2]*AMU*V_cm*V_cm + mass[2]*AMU*v3*V_cm*TMath::Cos( theta_cm[k] );	// [MeV]
		v_para[k] = V_cm + v3*TMath::Cos( theta_cm[k] );													// [1/c]
		v_perp[k] = v3*TMath::Sin( theta_cm[k] );															// [1/c]
		
		// Calculate the distances from the array at specific values of z
		r_shield[k][0] = CalculateRadius( cyclotron_freq, z_elum - Al_distance - Al_thickness, v_para[k], v_perp[k] );		// [mm]
		r_shield[k][1] = CalculateRadius( cyclotron_freq, z_elum - Al_distance, v_para[k], v_perp[k] );						// [mm]
		r_elum[k] = CalculateRadius( cyclotron_freq, z_elum, v_para[k], v_perp[k] );										// [mm]

		// Calculate the lab angle
		theta_lab[k] = TMath::ASin( v3*TMath::Sin( theta_cm[k] )/TMath::Sqrt( v_para[k]*v_para[k] + v_perp[k]*v_perp[k] ) );
		
		// Check if it the particle satisfies the conditions
		if ( ( r_shield[k][0] > radius_out || r_shield[k][0] < radius_in ) && ( r_shield[k][1] > radius_out || r_shield[k][1] < radius_in ) && r_elum[k] > radius_in && r_elum[k] < radius_out ){
			success_counter++;
			index_array[k] = k;
		}
		else{
			index_array[k] = -1;
		}
	}
	
	
	// Now create an array to label all the successes
	const Int_t num_success = success_counter;
	Int_t index_success[num_success];
	Int_t k = 0;
	for ( Int_t i = 0; i < range_angles; i++ ){
		if ( index_array[i] != -1 ){
			index_success[k] = index_array[i];
			k++;
		}
	}
	
	// CALCULATE THE SUCCESSFUL TRAJECTORIES IN FULL
	// Define graph and set formatting options
	Double_t z_scale = 0.4;									// Size of step in z [mm]
	Int_t size_of_array = (Int_t)(z_elum/z_scale + 1);		// Number of steps in z in total
	Double_t z_success[size_of_array];						// Holds the z values (same for each trajectory) [mm]
	Double_t theta_success[num_success];					// Holds the angles of all the successful trajectories
	Double_t r_success[num_success][size_of_array];			// Holds the radii of all of the successful trajectories [mm]
	Double_t T_lab_3_success[num_success];					// Holds the lab energy of the ejectile [MeV]
	Double_t max_radius = 0;								// Work out the maximum distance from the beam axis [mm]
	Int_t num_turns_completed[num_success];					// Tracks the number of turns
	Double_t r_thresh = 5.0;								// Threshold to help count the number of turns  [mm]

	// Populate arrays for positions and radii
	for ( Int_t i = 0; i < size_of_array; i++ ){

		// Populate z in constant steps
		z_success[i] = z_scale*i;
		for ( Int_t j = 0; j < num_success; j++ ){
			if ( i == 0 ){
				// Populate theta
				theta_success[j] = theta_lab[index_success[j]];
				
				// Populate the lab energies
				T_lab_3_success[j] = T_lab_3[index_success[j]];
			}
		
			// Calculate the radius along the track
			r_success[j][i] = CalculateRadius( cyclotron_freq, z_success[i], v_para[index_success[j]], v_perp[index_success[j]] );
		
			// Discriminate based on the number of turns completed
			if ( i == 0 ){ num_turns_completed[j] = 0; }
			if ( i > 0 && r_success[j][i] > r_thresh && r_success[j][i-1] <= r_thresh ){
				num_turns_completed[j]++;
			}
		
			// Find the max radius
			if ( r_success[j][i] > max_radius ){ max_radius = r_success[j][i]; }
		}
	}
	
	// NOW CALCULATE THE SOLID ANGLE
	Double_t total_solid_angle = 0;
	
	for ( Int_t i = 0; i < num_success; i++ ){
		// Discriminate based on number of turns
		if ( num_turns_completed[i] < 100 ){
			// Discriminate on the energy
			if ( T_lab_3_success[i] < 4.0 && T_lab_3_success[i] > 3.0 ){
				total_solid_angle += CalculateSolidAngle( 0.5*( theta_lab[index_success[i]] + theta_lab[index_success[i]-1] ) , 0.5*( theta_lab[index_success[i]+1] + theta_lab[index_success[i]] ) );
				printf("%5.3f\u00b0\t%i\t%f MeV\t%f sr\n", theta_success[i]*TMath::RadToDeg(), num_turns_completed[i], T_lab_3_success[i], total_solid_angle );
			}
		}
	}
			

	// Define new graphs
	TGraph *g_rz[num_success];
	TMultiGraph *mg = new TMultiGraph();

	// Populate and format the new graphs
	Int_t allowed_colours[14] = {1, 632, 416, 600, 400, 920, 616, 432, 800, 820, 840, 860, 880, 900};
	for ( Int_t i = 0; i < num_success; i++ ){
		g_rz[i] = new TGraph( size_of_array, z_success, r_success[i] );
		g_rz[i]->SetLineWidth(1);
		g_rz[i]->SetLineColor( allowed_colours[i%14] );
		for ( Int_t j = 0; j < (Int_t)num_revolutions.size(); j++ ){
			if ( num_turns_completed[i] == num_revolutions[j] ){
				mg->Add( g_rz[i] );
			}
		}
	}

	// Define a canvas
	TCanvas *c_rz = new TCanvas( "c_rz", "CANVAS", 1200, 900 );

	// Draw the first TGraph with axes
	TH1F *c_frame = c_rz->DrawFrame(0, -2, 1.1*z_elum, 1.1*max_radius );
	c_frame->SetTitle( "Distance from beam axis (r) v.s. distance along beam axis (z); z / mm; r / mm" );

	// Draw the rest with just curves on the same axes
	mg->Draw("C");

	// Define lines for r1 and r2
	TBox *shield = new TBox( z_elum - Al_distance - Al_thickness, radius_in, z_elum - Al_distance, radius_out );
	TLine *elum = new TLine( z_elum, radius_in, z_elum, radius_out );

	// Format the lines
	shield->SetFillColor(kRed);
	elum->SetLineWidth(4);
	elum->SetLineColor(kBlue);

	// Draw the lines
	shield->Draw("same");
	elum->Draw("same");
	c_rz->SetFrameLineColor(0);
	c_rz->Modified(); c_rz->Update();
	
	// Print the canvas
	if ( SWITCH_PRINT_CANVAS == 1 ){
		c_rz->Print("/home/ptmac/Documents/07-CERN-ISS-Mg/analysis/PLOTS/MgTraj.pdf");
	}


	return;
}

