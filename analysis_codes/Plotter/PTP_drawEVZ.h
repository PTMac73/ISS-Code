// PTP_drawEVZ.h
// Draw the E v.s. z plot
// ============================================================================================= //
// Patrick MacGregor
// Nuclear Physics Research Group
// School of Physics and Astronomy
// The University of Manchester
// ============================================================================================= //
#ifndef PTP_DRAW_EVZ_H_
#define PTP_DRAW_EVZ_H_

#include "PTPlotterINIT.h"
#include "PTF_GetPosNumber.h"
#include "PT_nuc_kinematics.h"
#include <TFile.h>
#include <TTree.h>
#include <TCutG.h>
#include <TString.h>
#include <TROOT.h>
#include <TSystem.h>
#include <THStack.h>
#include <TLine.h>
#include <TMath.h>
#include <TGraph.h>
#include <TH2.h>
#include <TCanvas.h>
#include <iostream>

struct plotterOptions;

// --------------------------------------------------------------------------------------------- //
// Draw EVZ plot
void DrawEVZ( TTree *t, plotterOptions &opt_s, Int_t pos_number ){
	gStyle->SetOptStat(kFALSE);
	
	// Print welcome message
	printDiv(); printf("PLOTTING THE E v.s. z PLOTS IN POSITION %i\n", pos_number ); printDiv();

	// Work out whether to iterate plots
	Int_t last_index;
	if ( SWITCH_ITERATE_PLOTS == 1 ){
		last_index = opt_s.numIter;
	}
	else{
		last_index = 1;
	}

	// Define variables
	TCanvas *c_evz[opt_s.numIter];
	TH2F *h_evz[opt_s.numIter];

	// Define axis limits
	Double_t x[4] = {-50,0,-10,9};
	const Int_t NUM_STATES = 9;		// Number of excited states in residual nucleus
	TGraph *g_EVZ[NUM_STATES];		// TGraphs for drawing excitation energy lines

	// CALCULATE THEORETICAL LINES (Non-relativistic)
	if ( DRAW_THEORETICAL_LINES_NR == 1 ){
		// Kinematic terms
		// M1 = 28Mg, M2 = 2H, M3 = 1H, M4 = 29Mg
		Double_t mass_excess[4] = { -15.018845, 13.13572176, 7.28897061, -10.602829};	// [MeV]
		Int_t mass_numbers[4] = { 28, 2, 1, 29 };										// Just a number
		Double_t mass[4];																// [AMU]
		for ( Int_t k = 0; k < 4; k++ ){
			mass[k] = MassExcessToMass( mass_numbers[k], mass_excess[k] );				
		}
		Int_t num_revolutions = 1;														// Number of turns
		Double_t ejectile_charge = 1;													// [e]
		Double_t b_field = 2.5;															// [T]
		Double_t beam_energy_per_nucleon = 9.473;										// MeV per u
		Double_t excitation_energy[NUM_STATES] = { 0, 1.054, 1.397, 2.225, 2.462, 2.937, 3.21, 3.937, 4.311 };
		Double_t rho_0 = ISSArrayRadius( -14.0, 14.0, 10.0 );							// Array radius in mm - this is estimated using width of Si strips

		// Calculate the range of angles (i.e. number of indices in array)
		Double_t angle_divisions = 0.2;
		Int_t lb_angle = 11;
		Int_t ub_angle = 50;
		const Int_t range_angles = TMath::Ceil( ( ub_angle - lb_angle )/angle_divisions ) + 1;

		// Constant kinetic terms
		Double_t Q = ( ( mass[0] + mass[1] ) - ( mass[2] + mass[3] ) )*AMU;					// Nuclear Q-value [MeV]
		Double_t T1 = beam_energy_per_nucleon*mass_numbers[0];								// Beam energy [MeV]
		Double_t T_cm_i = ( mass[1]/( mass[0] + mass[1] ) )*T1;								// Initial CM energy [MeV]
		Double_t V_cm = TMath::Sqrt( 2*T1*mass[0]*AMU )/( AMU*( mass[0] + mass[1] ) );		// CM velocity [1/c]
	
	
		// Define the quantites to be calculated in loop
		Double_t T_cm_f[NUM_STATES];							// Final energy in CM [MeV]
		Double_t v3[NUM_STATES];								// Final velocity of ejectile in CM [1/c]
		Double_t T_cm_3[NUM_STATES];							// Final energy of ejectile in CM [MeV]
		Double_t phi_cm[NUM_STATES][range_angles];				// Angle of scattering (intuitive) [rad]
		Double_t theta_cm[NUM_STATES][range_angles];			// Angle of scattering (true) [rad]
		Double_t T_lab_3[NUM_STATES][range_angles];				// Lab energy of ejectile in CM [MeV]
		Double_t v_para[NUM_STATES][range_angles];				// Parallel velocity of ejectile in CM
		Double_t v_perp[NUM_STATES][range_angles];				// Perpendicular velocity of ejectile in CM
		Double_t radius_circle[NUM_STATES][range_angles];		// Radius of circle traced out [mm]
		Double_t circle_path_length[NUM_STATES][range_angles];	// Length of circle traced out [mm]
		Double_t t_cyc[NUM_STATES][range_angles];				// Time to perform n revolutions (corrected)
		Double_t z[NUM_STATES][range_angles];					// z position [cm]
	
		// LOOP OVER NUMBER OF EXCITED STATES
		for ( Int_t j = 0; j < NUM_STATES; j++ ){
			// Calculate state-specific quantities
			T_cm_f[j] = T_cm_i + Q - excitation_energy[j];				// [MeV]
			T_cm_3[j] = mass[3]*T_cm_f[j]/( mass[2] + mass[3] );		// [MeV]
			v3[j] = TMath::Sqrt( 2*T_cm_3[j]/(mass[2]*AMU) );			// [1/c]
		
			// LOOP OVER THETA_CM
			for ( Int_t k = 0; k < range_angles; k++ ){
				// Calculate state and angle-specific quantities
				phi_cm[j][k] = ( lb_angle + k*angle_divisions )*TMath::DegToRad();
				theta_cm[j][k] = TMath::Pi() - phi_cm[j][k];
				T_lab_3[j][k] = T_cm_3[j] + 0.5*mass[2]*AMU*V_cm*V_cm + mass[2]*AMU*v3[j]*V_cm*TMath::Cos( theta_cm[j][k] );

				v_para[j][k] = V_cm + v3[j]*TMath::Cos( theta_cm[j][k] );	// [1/c]
				v_perp[j][k] = v3[j]*TMath::Sin( theta_cm[j][k] );		// [1/c]
			
				radius_circle[j][k] = MToMM()*AMUToKg( mass[2] )*( v_perp[j][k]*TMath::C() )/( ejectile_charge*ELECTRIC_CHARGE *b_field);
				circle_path_length[j][k] = 2*radius_circle[j][k]*( num_revolutions*TMath::Pi() -  TMath::ASin( rho_0/( 2*radius_circle[j][k] ) ) );
				t_cyc[j][k] = (circle_path_length[j][k]/MToMM() )/( v_perp[j][k]*TMath::C() );
				z[j][k] = v_para[j][k]*TMath::C()*t_cyc[j][k]*MToCM();
			}
		
			// Define graph and set formatting options
			g_EVZ[j] = new TGraph( range_angles, z[j], T_lab_3[j] );
			g_EVZ[j]->SetLineWidth(2);
		}
	}


	// Loop to create histograms
	TString plotTitle;
	for ( Int_t i = 0; i < last_index; i++ ){
		// Call the canvas
		c_evz[i] = new TCanvas( Form( "c_evz%i", i ), Form( "E v.s. z (%i/%i)", i+1, opt_s.numIter - 1 ), C_WIDTH, C_HEIGHT );
		
		
		// Define the histogram
		t->Draw( Form( "ecrr:z>>evz_%i(400, %f, %f, 900, %f, %f)", i, x[0], x[2], x[1], x[3] ), opt_s.cutList[i], "colz" );

		// Store the histogram
		h_evz[i] = (TH2F*)gDirectory->Get( Form( "evz_%i", i ) );

		// Format the histogram
		h_evz[i]->GetXaxis()->SetTitle("z / cm");
		h_evz[i]->GetYaxis()->SetTitle("E / MeV");
		
		// Produce title of histogram
		plotTitle = "E v.s. z with ";

		// Create the cut names
		if ( i == opt_s.numIter - 1 ){ plotTitle += opt_s.cutListDesc[i].Data(); }
		else{
			for ( Int_t j = i; j <= opt_s.numIter - 2; j++ ){
				plotTitle += opt_s.cutListDesc[j].Data();
				if ( j < opt_s.numIter - 2 ){ plotTitle += ", "; }
				if ( j == opt_s.numIter - 3 ){ plotTitle += "and "; }
			}
		}
		plotTitle +=" cut";
		if (i != opt_s.numIter - 2){ plotTitle +="s"; }
		printf("%i: %s\n", i, plotTitle.Data() );
		
		h_evz[i]->SetTitle( plotTitle.Data() );
		
		// DRAW SI STRIP DIVIDING LINES
		if ( DRAW_SI_STRIP_DIVIDERS == 1 ){
			// Define dividers
			TLine *l_divider[6][2];
			for ( Int_t j = 0; j < 6; j++ ){
				for ( Int_t k = 0; k < 2; k++ ){
					// Define positions for each line
					l_divider[j][k] = new TLine();
					l_divider[j][k]->SetX1( TMath::Power(-1.0,k)*2.5 - Z_OFF[pos_number-1] - Z_ARRAY_POS[j] );
					l_divider[j][k]->SetX2( l_divider[j][k]->GetX1() );
					l_divider[j][k]->SetY1( x[1] );
					l_divider[j][k]->SetY2( x[3] );

					// Define styles
					l_divider[j][k]->SetLineWidth(2);
					if ( j % 2 == 0 ){
						l_divider[j][k]->SetLineColor(kBlack);
					}
					else{
						l_divider[j][k]->SetLineColor(kBlue);
					}
					// Draw the line
					l_divider[j][k]->Draw("same");
				}
			}		
		}
		
		// Draw the theoretical lines
		if ( DRAW_THEORETICAL_LINES_NR == 1 ){
			for ( Int_t j = 0; j < NUM_STATES; j++ ){
				g_EVZ[j]->Draw("L");
			}
		}
		
		// Update the canvas
		c_evz[i]->Modified(); c_evz[i]->Update();

		// Print if desired
		// Print the canvas
		if ( SWITCH_PRINT_CANVAS == 1 ){
			c_evz[i]->Print( Form( "%sevz_%i%s", opt_s.printDir.Data(), i, PRINT_FORMAT.Data() ) );
		}
	}
	gStyle->SetOptStat(kTRUE);
	
}

#endif  // PTP_DRAW_EVZ_H_
