// ELUMSolidAngleCalculator.C
// Calculates the solid angle for the luminosity detector, by solving a transcendental equation
// ============================================================================================= //
// Patrick MacGregor
// Nuclear Physics Research Group
// School of Physics and Astronomy
// The University of Manchester
// ============================================================================================= //
#include <TBox.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TH1.h>
#include <THStack.h>
#include <TLine.h>
#include <TMath.h>
#include <TMultiGraph.h>
#include <TString.h>
#include <TStyle.h>

#include <iostream>

#include "ELUMSolidAngleCalculator.h"
#include "conversion.h"

// MAIN FUNCTION =============================================================================== //
// z_elum is the distance from target to elum [mm]
Double_t SolidAngleCalculator( Double_t z_elum = 125.7/*134.8*/ ){
	// OPTIONS
	if ( SWITCH_DRAW_CANVAS == 0 ){
		  gROOT->SetBatch(kTRUE);
	}


	// CONSTANTS OF MOTION --------------------------------------------------------------------- //
	// M1 = 28Mg, M2 = 2H, M3 = 2H, M4 = 28Mg
	// 28Mg(d,d)28Mg
	Double_t mass_excess[4] = { -15.018845, 13.13572176, 13.13572176, -15.018845};	// [MeV]
	Double_t charge_numbers[4] = { 12, 1, 1, 12 };									// [e]
	Int_t mass_numbers[4] = { 28, 2, 2, 28 };										// Just a number
	/*
	// 28Si(d,d)28Si
	Double_t mass_excess[4] = { -21.4927, 13.13572176, 13.13572176, -21.4927};	// [MeV]
	Double_t charge_numbers[4] = { 14, 1, 1, 14 };									// [e]
	Int_t mass_numbers[4] = { 28, 2, 2, 28 };										// Just a number
	*/
	Double_t mass[4];																// [AMU]
	for ( Int_t k = 0; k < 4; k++ ){
		mass[k] = MassExcessToMass( mass_numbers[k], mass_excess[k], charge_numbers[k] );				
	}

	Double_t ejectile_charge = 1;													// [e]
	Double_t b_field = 2.5;															// [T]
	Double_t beam_energy_per_nucleon = 9.473;										// MeV per u
	Double_t excitation_energy = 0.0;
	Double_t cyclotron_freq = ejectile_charge*EToCoulombs()*b_field/( mass[2]*AMUToKg() );	//[s^-1]

	// Define distances
	Double_t z_elum_error = 0.7;		// Error on distance from target to elum [mm]
	Double_t Al_distance = 13.42;		// Distance from elum to Al shield [mm] --> eLog 209
	Double_t Al_thickness = 12.8;		// Thickness of the Al shield [mm] --> eLog 209
	Double_t tube_thickness = 2.0;		// Thickness of the tube (unknown) [mm]
	Double_t radius_in = 48.00;			// Inner radius of elum [mm]
	Double_t radius_out = 96.00;		// Outer radius of elum [mm]

	// Define angle limits
	Double_t angle_limits[CHECK_NUM_REV][CHECK_TYPES][CHECK_MINMAX] = {
	{							// 1 REVOLUTION
		{ 100.0, 0.0 },			// Total angle limits
		{ 100.0, 0.0 },			// Below shield success angle limits
		{ 100.0, 0.0 },			// Above shield success angle limits
	},
	{							// 2 REVOLUTIONS
		{ 100.0, 0.0 },			// Total angle limits
		{ 100.0, 0.0 },			// Below shield success angle limits
		{ 100.0, 0.0 },			// Above shield success angle limits
	}
};

	// Calculate the range of angles (i.e. number of indices in array)
	Double_t angle_divisions = 0.001;
	const Double_t THETA_SPACING = 0.15;
	Int_t lb_angle = 0;
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
	
	MMM_t mmm_ejectile_energy;					// Calculate the mean ejectile energy
	mmm_ejectile_energy.min = 100000;
	mmm_ejectile_energy.max = 0;
	mmm_ejectile_energy.mean = 0;
	mmm_ejectile_energy.num = 0;

	// Define the quantites to be calculated in loop over angles
	Double_t *theta_cm = new Double_t[range_angles];		// Angle of scattering in CM [rad]
	Double_t *theta_lab = new Double_t[range_angles];		// Angle of scattering in LAB
	Double_t *T_lab_3 = new Double_t[range_angles];			// Lab energy of ejectile in CM [MeV]
	Double_t *v_para = new Double_t[range_angles];			// Parallel velocity of ejectile in LAB [1/c]
	Double_t *v_perp = new Double_t[range_angles];			// Perpendicular velocity of ejectile in LAB [1/c]
	Traj_t *traj_arr = new Traj_t[range_angles];			// Stores the trajectory information
	
	// Bools to work out their exit status
	Bool_t brun_red;
	Bool_t brun_orange;
	Bool_t brun_green;
	Bool_t in_shield_z = 0;
	Bool_t fail_flag = 0;									// Goes to 1 if you hit something
	
	// Define graph and set formatting options
	Double_t z_scale = 0.5;									// Size of step in z [mm]
	Int_t size_of_array = (Int_t)(z_elum/z_scale + 1);		// Number of steps in z in total
	Double_t z[size_of_array];								// Holds the z values [mm]
	Double_t r[size_of_array];								// Holds the r values [mm]
	Double_t max_radius = 0;								// Work out the maximum distance from the beam axis [mm]
	Double_t total_solid_angle = 0;
	Double_t fail_z, fail_r, fail_index;					// Define points where failure occurred
	TGraph *g_rz[range_angles];
	TMultiGraph *mg = new TMultiGraph();
	
	// Histogram for recording energies
	TH1D* h_elum_energy[7];
	TH1D* h_ptr;
	for ( Int_t i = 0; i < 7; i++ ){
		h_elum_energy[i] = new TH1D( Form( "h_elum_energy%i", i ), "ENERGY!", 1000, 0, 5 );
	}
	h_ptr = h_elum_energy[3];
	THStack* hs = new THStack( "hs", "" );


	// CALCULATE TRAJECTORIES ------------------------------------------------------------------ //
	// LOOP OVER THETA_CM
	for ( Int_t k = 0; k < range_angles; k++ ){
		// Calculate angle-specific quantities
		theta_cm[k] = ( lb_angle + k*angle_divisions )*TMath::DegToRad();									// [rad]
		T_lab_3[k] = T_cm_3 + 0.5*mass[2]*AMU*V_cm*V_cm + mass[2]*AMU*v3*V_cm*TMath::Cos( theta_cm[k] );	// [MeV]
		v_para[k] = V_cm + v3*TMath::Cos( theta_cm[k] );													// [1/c]
		v_perp[k] = v3*TMath::Sin( theta_cm[k] );															// [1/c]
		theta_lab[k] = TMath::ASin( v3*TMath::Sin( theta_cm[k] )/TMath::Sqrt( v_para[k]*v_para[k] + v_perp[k]*v_perp[k] ) );

		// Reset the number of completed turns
		//traj_arr[k].rev = 0;
		if ( v_para[k] == 0 ){
			traj_arr[k].rev = -1;
		}
		else{
			traj_arr[k].rev = TMath::Ceil( cyclotron_freq*z_elum*MMToM()/( TMath::TwoPi()*v_para[k]*TMath::C() ) ) - 1;
		}
		
		// Reset the flags and fail values
		fail_flag = 0;
		fail_z = 0.0;
		fail_r = 0.0;
		fail_index = 0.0;

		// LOOP OVER Z - calculate full trajectory and then fix later
		for ( Int_t j = 0; j < size_of_array; j++ ){
		
			// Populate z in constant steps
			z[j] = z_scale*j;
			
			// Calculate the radius along the track
			r[j] = CalculateRadius( cyclotron_freq, z[j], v_para[k], v_perp[k] );
		
			// Test conditions if it has not hit anything
			if ( fail_flag == 0 ){
			
				// Find the max radius
				if ( r[j] > max_radius ){ max_radius = r[j]; }
				
				// BAD TRAJECTORIES -- IMPOSSIBLE
				if ( TMath::IsNaN( r[j] ) ){
					fail_flag = 1;
					fail_index = j;
					fail_r = 0.0;
					fail_z = z[j];
				}
				
				// BAD TRAJECTORIES -- SHIELD
				// Check for bad trajectories from bottom, top, and left hitting shield
				if( z[j] >= z_elum - Al_distance - Al_thickness && z[j] <= z_elum - Al_distance ){
					if ( r[j] >= radius_in - tube_thickness && r[j-1] < radius_in - tube_thickness ){ // Approached from bottom
						Double_t m = ( r[j] - r[j-1] )/( z[j] - z[j-1] );
						fail_z = ( radius_in - tube_thickness - ( r[j] - m*z[j] ) )/m;
						fail_r = radius_in - tube_thickness;
						fail_index = j;
						fail_flag = 1;
					}
					else if ( r[j] <= radius_out && r[j-1] > radius_out ){ // Approached from top
						Double_t m = ( r[j] - r[j-1] )/( z[j] - z[j-1] );
						fail_z = ( radius_out - ( r[j] - m*z[j] ) )/m;
						fail_r = radius_out;
						fail_index = j;
						fail_flag = 1;
					}
					else if ( z[j-1] < z_elum - Al_distance - Al_thickness && r[j-1] >= radius_in - tube_thickness && r[j-1] <= radius_out ){// Approached from left - test at z_elum to see if it fails
						if ( CalculateRadius( cyclotron_freq, z_elum - Al_distance - Al_thickness, v_para[k], v_perp[k] ) >= radius_in - tube_thickness && CalculateRadius( cyclotron_freq, z_elum - Al_distance - Al_thickness, v_para[k], v_perp[k] ) <= radius_out ){
							fail_z = z_elum - Al_distance - Al_thickness;
							fail_r = CalculateRadius( cyclotron_freq, z_elum - Al_distance - Al_thickness, v_para[k], v_perp[k] );
							fail_index = j;
							fail_flag = 1;
						}
					}
				}
				// BAD TRAJECTORIES -- TUBE
				else if ( z[j] >= z_elum - Al_distance && z[j] <= z_elum ){
					if ( r[j] >= radius_in - tube_thickness && r[j-1] < radius_in - tube_thickness ){ // Approached from bottom
						Double_t m = ( r[j] - r[j-1] )/( z[j] - z[j-1] );
						fail_z = ( radius_in - tube_thickness - ( r[j] - m*z[j] ) )/m;
						fail_r = radius_in - tube_thickness;
						fail_index = j;
						fail_flag = 1;
					}
					else if ( r[j] <= radius_in && r[j-1] > radius_in ){ // Approached from top
						Double_t m = ( r[j] - r[j-1] )/( z[j] - z[j-1] );
						fail_z = ( radius_in - ( r[j] - m*z[j] ) )/m;
						fail_r = radius_in;
						fail_index = j;
						fail_flag = 1;
					}
				}
				
				// Calculate the number of turns (of full trajectory)
				/*if ( j > 1 && r[j] > r[j-1] && r[j-1] < r[j-2] ){
					traj_arr[k].rev++;
				}*/
				
				// Calculate if it goes above the shield or below it
				if ( CalculateRadius( cyclotron_freq, z_elum - Al_distance - Al_thickness, v_para[k], v_perp[k] ) < radius_out ){
					traj_arr[k].above_shield = 0;
				}
				else{ traj_arr[k].above_shield = 1; }
				
			}

		}

		// Fix trajectories if they failed
		if ( fail_flag == 1 ){
			for ( Int_t j = fail_index; j < size_of_array; j++ ){
				r[j] = fail_r;
				z[j] = fail_z;
			}
		}

		
		// Calculate the status of the particle
		// Calculate radius at end of array
		Double_t rtemp = CalculateRadius( cyclotron_freq, z_elum, v_para[k], v_perp[k] );

		// Particle hits shield
		if ( fail_flag == 1 ){
			// Particle would've hit elum
			if ( rtemp >= radius_in && rtemp <= radius_out ){
				traj_arr[k].status = 1;
			}
			// Particle misses elum
			else{
				traj_arr[k].status = 2;
			}
		}
		// Particle misses shield, but does it hit elum?
		else{
			if ( rtemp >= radius_in && rtemp <= radius_out ){
				traj_arr[k].status = 0;
			}
			// Particle misses elum
			else{
				traj_arr[k].status = 3;
			}
		}

		// Now check for angle limits etc.
		Double_t angle_thresh = 5;
		if ( traj_arr[k].rev >= 0 && traj_arr[k].rev <= 1  ){
			if ( traj_arr[k].status == 0 ){
				// Test above and below shield limits
				if ( theta_cm[k] < angle_limits[ traj_arr[k].rev ][ traj_arr[k].above_shield + 1 ][ 0 ] && ( TMath::Abs( theta_cm[k] - angle_limits[ traj_arr[k].rev ][ traj_arr[k].above_shield + 1 ][ 0 ] ) < angle_thresh*TMath::DegToRad() || angle_limits[ traj_arr[k].rev ][ traj_arr[k].above_shield + 1 ][ 0 ] > TMath::Pi() ) ){ 
					angle_limits[ traj_arr[k].rev ][ traj_arr[k].above_shield + 1 ][ 0 ] = theta_cm[k];
				}
				if ( theta_cm[k] > angle_limits[ traj_arr[k].rev ][ traj_arr[k].above_shield + 1 ][ 1 ] && ( TMath::Abs( theta_cm[k] - angle_limits[ traj_arr[k].rev ][ traj_arr[k].above_shield + 1 ][ 1 ] ) < angle_thresh*TMath::DegToRad() || angle_limits[ traj_arr[k].rev ][ traj_arr[k].above_shield + 1 ][ 1 ] < 0.01 ) ){
					angle_limits[ traj_arr[k].rev ][ traj_arr[k].above_shield + 1 ][ 1 ] = theta_cm[k];
				}
			}
			
			// Test regular maximum and minimum for all trajectories
			if ( theta_cm[k] < angle_limits[ traj_arr[k].rev ][ 0 ][ 0 ] && ( TMath::Abs( theta_cm[k] - angle_limits[ traj_arr[k].rev ][ 0 ][ 0 ] ) < angle_thresh*TMath::DegToRad() || angle_limits[ traj_arr[k].rev ][ 0 ][ 0 ] > TMath::Pi() ) ){
				angle_limits[ traj_arr[k].rev ][ 0 ][ 0 ] = theta_cm[k];
			}
			if ( theta_cm[k] > angle_limits[ traj_arr[k].rev ][ 0 ][ 1 ] && ( TMath::Abs( theta_cm[k] - angle_limits[ traj_arr[k].rev ][ 0 ][ 1 ] ) < angle_thresh*TMath::DegToRad() || angle_limits[ traj_arr[k].rev ][ 0 ][ 1 ] < 0.01 ) ){
				angle_limits[ traj_arr[k].rev ][ 0 ][ 1 ] = theta_cm[k];
			}
		}

		
		// Calculate the solid angle
		/*if ( traj_arr[k].status == 0 ){
			total_solid_angle += CalculateSolidAngle( 0.5*( theta_cm[k] + theta_cm[k - 1] ) , 0.5*( theta_cm[k + 1] + theta_cm[k] ) );
		}*/
		
		// Format graphs
		g_rz[k] = new TGraph( size_of_array, z, r );
		g_rz[k]->SetLineWidth(1);
		if ( traj_arr[k].status == 0 ){ g_rz[k]->SetLineColor( kGreen ); }
		else if ( traj_arr[k].status == 1 ){ g_rz[k]->SetLineColor( kOrange ); }
		else if ( traj_arr[k].status == 2 ){ g_rz[k]->SetLineColor( kRed ); }
		else if ( traj_arr[k].status == 3 ){ g_rz[k]->SetLineColor( kBlue ); }
		
		
		// Add to multigraph
		if ( ( k % (Int_t)(THETA_SPACING/angle_divisions) == 0 ) && traj_arr[k].rev == 0 ){ mg->Add( g_rz[k] ); }
		
		// Add results to histogram
		if ( traj_arr[k].status == 0 ){
			if ( traj_arr[k].rev == 0 && traj_arr[k].above_shield == 0 ){ h_elum_energy[0]->Fill( T_lab_3[k] ); }
			else if ( traj_arr[k].rev == 0 && traj_arr[k].above_shield == 1 ){ h_elum_energy[1]->Fill( T_lab_3[k] ); CalcRunningAverage( mmm_ejectile_energy, T_lab_3[k] ); }
			else if ( traj_arr[k].rev == 1 && traj_arr[k].above_shield == 0 ){ h_elum_energy[2]->Fill( T_lab_3[k] ); }
			else if ( traj_arr[k].rev == 1 && traj_arr[k].above_shield == 1 ){ h_elum_energy[3]->Fill( T_lab_3[k] ); }
			else if ( traj_arr[k].rev == 2 && traj_arr[k].above_shield == 0 ){ h_elum_energy[4]->Fill( T_lab_3[k] ); }
			else if ( traj_arr[k].rev == 2 && traj_arr[k].above_shield == 1 ){ h_elum_energy[5]->Fill( T_lab_3[k] ); }
			else{ h_elum_energy[6]->Fill( T_lab_3[k] ); }
		}
	}

	// Define a canvas
	TCanvas *c_rz = new TCanvas( "c_rz", "CANVAS", GetCanvasWidth(C_HEIGHT, 0, -0.05*max_radius, 1.1*z_elum, 1.1*max_radius ), C_HEIGHT );
	GlobSetCanvasMargins( c_rz );

	// Draw the frame with axes
	TH1F *c_frame = c_rz->DrawFrame(0, -0.05*max_radius, 1.1*z_elum, 1.1*max_radius );
	FormatFrame( c_frame, "z (mm)", "r (mm)" );
	c_rz->SetFrameLineColor(0);

	// Define shield and elum detector
	TBox *shield = new TBox( z_elum - Al_distance - Al_thickness, radius_in, z_elum - Al_distance, radius_out );
	TBox *tube = new TBox( z_elum - Al_distance - Al_thickness, radius_in - tube_thickness, z_elum + 0.5, radius_in );
	TBox *elum = new TBox( z_elum - 0.2, radius_in, z_elum + 0.5, radius_out );

	// Format the drawings
	shield->SetFillColor(kGray);
	tube->SetFillColor(kGray);
	elum->SetFillColor(pcb_green_i);

	// Draw the objects
	shield->Draw("same");
	tube->Draw("same");
	elum->Draw("same");

	// Draw the trajectories with just curves on the same axes
	mg->Draw("C");
	c_rz->Modified(); c_rz->Update();
	// Print the canvas
	if ( SWITCH_PRINT_CANVAS == 1 ){
		c_rz->Print("ELUMSolidAngle.png");
		c_rz->Print("ELUMSolidAngle.tex");
		c_rz->Print("ELUMSolidAngle.svg");
		c_rz->Print("ELUMSolidAngle.pdf");
	}

	// Print the angle limits
	/*for ( Int_t i = 0; i < CHECK_NUM_REV; i++ ){
		std::cout << "REV " << i + 1 << "\n";
		for ( Int_t j = 0; j < CHECK_TYPES; j++ ){
			std::cout << std::right << std::setw(8) << 180 - angle_limits[i][j][1]*TMath::RadToDeg() << " ---> " << std::setw(8) << 180 - angle_limits[i][j][0]*TMath::RadToDeg() << "\n";
		}
	}*/
	
	// Draw the histogram stack
	h_elum_energy[0]->SetFillColor( kBlack );
	h_elum_energy[1]->SetFillColor( kGreen );
	h_elum_energy[2]->SetFillColor( kBlack );
	h_elum_energy[3]->SetFillColor( kCyan );
	h_elum_energy[4]->SetFillColor( kBlack );
	h_elum_energy[5]->SetFillColor( kMagenta );
	h_elum_energy[6]->SetFillColor( kBlack );

	h_elum_energy[0]->SetLineColor( kBlack );
	h_elum_energy[1]->SetLineColor( kGreen );
	h_elum_energy[2]->SetLineColor( kBlack );
	h_elum_energy[3]->SetLineColor( kCyan );
	h_elum_energy[4]->SetLineColor( kBlack );
	h_elum_energy[5]->SetLineColor( kMagenta );
	h_elum_energy[6]->SetLineColor( kBlack );
	
	for ( Int_t i = 0; i < 7; i++ ){
		hs->Add( h_elum_energy[i] );
	}

	TCanvas* c_h = new TCanvas( "c_h", "CANVAS", 1200, 900 );
	GlobSetCanvasMargins( c_h );
	gStyle->SetOptStat();
	hs->Draw();
	c_h->Update();
	FormatFrame( hs, "T_{3} (MeV)", "Counts" );
	
	std::cout << "MEAN = " << mmm_ejectile_energy.mean << "\n";
	std::cout << "MIN  = " << mmm_ejectile_energy.min  << "\n";
	std::cout << "MAX  = " << mmm_ejectile_energy.max  << "\n";
	std::cout << "NUM  = " << mmm_ejectile_energy.num  << "\n";
	
	
	if ( SWITCH_PRINT_CANVAS == 1 ){
		c_h->Print("elum_energy_hist.png");
		c_h->Print("elum_energy_hist.tex");
		c_h->Print("elum_energy_hist.svg");
		c_h->Print("elum_energy_hist.pdf");
	}
	
	std::cout << z_elum << "\t" << angle_limits[0][2][0]*TMath::RadToDeg() << "\t" << angle_limits[0][2][1]*TMath::RadToDeg() << "\n";
	

	// Delete arrays that are not needed any more
	delete[] theta_lab;
	delete[] theta_cm;
	delete[] T_lab_3;
	delete[] v_para;
	delete[] v_perp;
	
	if ( SWITCH_DRAW_CANVAS == 0 ){
		delete c_rz;
		delete c_h;
		delete mg;
		for ( Int_t i = 0; i < 7; i++ ){
			if ( h_elum_energy[i]->IsOnHeap() ){ h_elum_energy[i]->Delete(); }
		}
	}
	
	return 0;
}



// BATCH FUNCTION ============================================================================== //
void ELUMSolidAngleCalculator(){
	// Calculate solid angles for different values of z_elum
	const Int_t num_z = 5;
	Double_t z_elum[num_z] = { 120.7, 124.7, 125.7, 126.7, 130.7 };
	Double_t solid_angle[num_z];
	if (SWITCH_VERBOSE == 0 ){ std::cout << "  Z (mm)\tmin(\u03b8CM)\tmax(\u03b8CM)" << std::endl; }
	for ( Int_t i = 2; i < 3; i++ ){
		solid_angle[i] = SolidAngleCalculator( z_elum[i] );
	}

	return;
}





























