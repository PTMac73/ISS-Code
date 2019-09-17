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

// DEFINE GLOBAL SWITCHES
Bool_t SWITCH_DRAW_CANVAS = 0;
Bool_t SWITCH_PRINT_CANVAS = 0;
Bool_t SWITCH_VERBOSE = 0;
Int_t C_SCALE = 300;
Int_t C_WIDTH = 4*C_SCALE;
Int_t C_HEIGHT = 3*C_SCALE;

// Define radius function [ mm ] - > units are [ [s^-1], [mm], [1/c], [1/c] ]
Double_t CalculateRadius( Double_t cyclotron_freq, Double_t z, Double_t v_para, Double_t v_perp ){
	return ( v_perp*TMath::C()*MToMM()/cyclotron_freq )*( 1.0 - TMath::Cos( cyclotron_freq*z/( v_para*TMath::C()*MToMM() ) ) );
}

Double_t CalculateSolidAngle( Double_t theta1, Double_t theta2 ){
	return 2*TMath::Pi()*( TMath::Cos( theta1 ) - TMath::Cos( theta2 ) );
}



// MAIN FUNCTION =============================================================================== //
// z_elum is the distance from target to elum [mm]
Double_t SolidAngleCalculator( Double_t z_elum = 134.8 ){
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
	Double_t z_elum_error = 0.7;		// Error on distance from target to elum [mm]
	Double_t Al_distance = 13.42;		// Distance from elum to Al shield [mm]
	Double_t Al_thickness = 12.8;		// Thickness of the Al shield [mm]
	Double_t radius_in = 48.00;			// Inner radius of elum [mm]
	Double_t radius_out = 96.00;		// Outer radius of elum [mm]

	// Calculate the range of angles (i.e. number of indices in array)
	Double_t angle_divisions = 0.002;
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
	Double_t *theta_cm = new Double_t[range_angles];		// Angle of scattering in CM [rad]
	Double_t *theta_lab = new Double_t[range_angles];		// Angle of scattering in LAB
	Double_t *T_lab_3 = new Double_t[range_angles];			// Lab energy of ejectile in CM [MeV]
	Double_t *v_para = new Double_t[range_angles];			// Parallel velocity of ejectile in LAB [1/c]
	Double_t *v_perp = new Double_t[range_angles];			// Perpendicular velocity of ejectile in LAB [1/c]
	Double_t r_shield_0 = 0.0, r_shield_1 = 0.0;			// Distance from beam axis at shield (both edges) [mm]
	Double_t r_elum = 0.0;									// Distance from beam axis at elum [mm]
	vector <Int_t> index_temp;								// Records the indices of all the successful runs

	// LOOP OVER THETA_CM
	for ( Int_t k = 0; k < range_angles; k++ ){
		theta_cm[k] = ( lb_angle + k*angle_divisions )*TMath::DegToRad();									// [rad]
		T_lab_3[k] = T_cm_3 + 0.5*mass[2]*AMU*V_cm*V_cm + mass[2]*AMU*v3*V_cm*TMath::Cos( theta_cm[k] );	// [MeV]
		v_para[k] = V_cm + v3*TMath::Cos( theta_cm[k] );													// [1/c]
		v_perp[k] = v3*TMath::Sin( theta_cm[k] );															// [1/c]
		
		// Calculate the distances from the array at specific values of z
		r_shield_0 = CalculateRadius( cyclotron_freq, z_elum - Al_distance - Al_thickness, v_para[k], v_perp[k] );		// [mm]
		r_shield_1 = CalculateRadius( cyclotron_freq, z_elum - Al_distance, v_para[k], v_perp[k] );						// [mm]
		r_elum = CalculateRadius( cyclotron_freq, z_elum, v_para[k], v_perp[k] );										// [mm]

		// Calculate the lab angle
		theta_lab[k] = TMath::ASin( v3*TMath::Sin( theta_cm[k] )/TMath::Sqrt( v_para[k]*v_para[k] + v_perp[k]*v_perp[k] ) );
		
		// Check if it the particle satisfies the conditions of:
		/*	* The radius at the left and right-hand-side of the shield is > outer radius of the shield i.e. goes above the shield
			* The radius at the elum is > inner radius of elum and < outer radius i.e. it hits the detector
		*/	
		if ( r_shield_0 > radius_out && r_shield_1 > radius_out && r_elum > radius_in && r_elum < radius_out ){
			index_temp.push_back(k);
		}
	}

	// Found all the successful angles - now test to see how many turns they complete
	Double_t z_scale_1 = 1.0;									// Define a relatively large z over which to test things
	Int_t size_of_array_1 = (Int_t)(z_elum/z_scale_1 + 1);		// Define the size of an array with this z scale
	Double_t r_temp[size_of_array_1];							// Define a temporary radius
	vector < Int_t > index_success;								// Define a vector to store the successful indices
	vector < Int_t > num_turns_completed;						// Counter for the number of completed turns in a given trajectory
	Int_t num_turns_temp;										// Define a temporary counter for the number of turns
	Double_t r_thresh = 5.0;									// Threshold to help count the number of turns  [mm]

	// Loop over theta
	for ( Int_t i = 0; i < (Int_t)index_temp.size(); i++ ){
		// Reset the number of completed turns
		num_turns_temp = 0;

		// Loop over z
		for ( Int_t j = 0; j < size_of_array_1; j++ ){
			// Calculate the radius along the track
			r_temp[j] = CalculateRadius( cyclotron_freq, z_scale_1*j, v_para[index_temp[i]], v_perp[index_temp[i]] );
			
			// Calculate the number of turns
			if ( j > 0 && r_temp[j] > r_thresh && r_temp[j-1] <= r_thresh ){
				num_turns_temp++;
			}

		}
	
		// Store the index if it is good
		for ( Int_t k = 0; k < (Int_t)num_revolutions.size(); k++ ){
			if ( num_turns_temp == num_revolutions[k] ){
				index_success.push_back( index_temp[i] );
				num_turns_completed.push_back( num_turns_temp );
			}
		}
		
	}

	// CALCULATE THE SUCCESSFUL TRAJECTORIES IN FULL
	// Define graph and set formatting options
	Double_t z_scale = 0.1;										// Size of step in z [mm]
	Int_t size_of_array = (Int_t)(z_elum/z_scale + 1);			// Number of steps in z in total
	Double_t z_success[size_of_array];							// Holds the z values (same for each trajectory) [mm]
	Double_t theta_success[index_success.size()];				// Holds the angles of all the successful trajectories
	Double_t r_success[index_success.size()][size_of_array];	// Holds the radii of all of the successful trajectories [mm]
	Double_t T_lab_3_success[index_success.size()];				// Holds the lab energy of the ejectile [MeV]
	Double_t max_radius = 0;									// Work out the maximum distance from the beam axis [mm]
	Double_t total_solid_angle = 0;


	// Populate arrays for positions and radii - loop over theta
	if ( SWITCH_VERBOSE == 1 ){ printf( "NUM_REV \tTHETA_CM\tENERGY\t\tSOLID ANGLE\n" ); }
	for ( Int_t i = 0; i < (Int_t)index_success.size(); i++ ){
		// Populate theta
		theta_success[i] = theta_cm[index_success[i]];
				
		// Populate the lab energies
		T_lab_3_success[i] = T_lab_3[index_success[i]];

		// Calculate the solid angle (discriminating on the energy so that you get one peak)
		total_solid_angle += CalculateSolidAngle( 0.5*( theta_cm[index_success[i]] + theta_cm[index_success[i]-1] ) , 0.5*( theta_cm[index_success[i]+1] + theta_cm[index_success[i]] ) );
		if ( SWITCH_VERBOSE == 1 ){ printf("%8i\t%8.4f\u00b0\t%8f MeV\t%8f sr\n", num_turns_completed[i], theta_success[i]*TMath::RadToDeg(), T_lab_3_success[i], total_solid_angle ); }
		else{
			if ( i == 0 ){ std::cout << std::setw(8) << z_elum << "\t" << std::setw(7) << theta_success[i]*TMath::RadToDeg() << "\u00b0\t"; }
			if ( i == (Int_t)index_success.size() - 1 ){ std::cout << std::setw(7) << theta_success[i]*TMath::RadToDeg() << "\u00b0\t" << std::setw(12) << total_solid_angle << std::endl; }
		}
		

		// Loop over z
		for ( Int_t j = 0; j < size_of_array; j++ ){
			// Populate z in constant steps
			if ( i == 0 ){ z_success[j] = z_scale*j; }

			// Calculate the radius along the track
			r_success[i][j] = CalculateRadius( cyclotron_freq, z_success[j], v_para[index_success[i]], v_perp[index_success[i]] );
			
			// Find the max radius
			if ( r_success[i][j] > max_radius ){ max_radius = r_success[i][j]; }
		}
	}
	
	
	// Define new graphs
	TGraph *g_rz[index_success.size()];
	TMultiGraph *mg = new TMultiGraph();

	// Populate and format the new graphs
	Int_t allowed_colours[14] = {1, 632, 416, 600, 400, 920, 616, 432, 800, 820, 840, 860, 880, 900};
	for ( Int_t i = 0; i < (Int_t)index_success.size(); i++ ){
		g_rz[i] = new TGraph( size_of_array, z_success, r_success[i] );
		g_rz[i]->SetLineWidth(1);
		g_rz[i]->SetLineColor( allowed_colours[i%14] );
		mg->Add( g_rz[i] );
	}

	// Define a canvas
	TCanvas *c_rz = new TCanvas( "c_rz", "CANVAS", C_WIDTH, C_HEIGHT );

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
		c_rz->Print("/home/ptmac/Documents/07-CERN-ISS-Mg/analysis/PLOTS/MgTraj.png");
	}


	// Delete arrays that are not needed any more
	delete[] theta_lab;
	delete[] theta_cm;
	delete[] T_lab_3;
	delete[] v_para;
	delete[] v_perp;
	if ( SWITCH_DRAW_CANVAS == 0 ){
		delete c_rz;
		delete mg;
	}
	
	return total_solid_angle;
}



// BATCH FUNCTION ============================================================================== //
void ELUMSolidAngleCalculator(){
	// Calculate solid angles for different values of z_elum
	const Int_t num_z = 5;
	Double_t z_elum[num_z] = { 120.7, 124.7, 125.7, 126.7, 130.7 };
	Double_t solid_angle[num_z];
	if (SWITCH_VERBOSE == 0 ){ std::cout << "  Z (mm)\tmin(\u03b8CM)\tmax(\u03b8CM)\tSOLID_ANGLE" << std::endl; }
	for ( Int_t i = 0; i < num_z; i++ ){
		solid_angle[i] = SolidAngleCalculator( z_elum[i] );
	}

	return;
}





























