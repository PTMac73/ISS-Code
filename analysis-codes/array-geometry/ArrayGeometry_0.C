// ArrayGeometry.C
// Checks the excitation of each state and sees where it starts getting clipped
// ============================================================================================= //
// Patrick MacGregor
// Nuclear Physics Research Group
// School of Physics and Astronomy
// The University of Manchester
// ============================================================================================= //
#include "AG_physical_constants.h"
#include "AG_constants.h"
#include "AG_functions.h"
#include "AG_style.h"

#include <TCanvas.h>
#include <TGraph.h>
#include <TH1F.h>
#include <TBox.h>
#include <TEllipse.h>

#include <iostream>

void ArrayGeometry(){
	// Set style
	CreateStyle( ptm_style );
	ptm_style->cd();
	//gROOT->SetBatch(kTRUE);

	// DRAW STUFF ------------------------------------------------------------------------------ //
	TCanvas* c0 = new TCanvas( "c_array", "Array Trajectories", CANVAS_WIDTH, CANVAS_HEIGHT );
	
	// Draw frame
	TH1F* frame = c0->DrawFrame( X1, Y1, X2, Y2 );
	frame->SetTitle( ";z / cm; r / cm" );

	// Draw target
	TBox* b_target = new TBox( 0, -2, 0.5, 2 );
	b_target->SetFillColor(kBlack);
	b_target->Draw("SAME");

	// Draw PCB
	TBox* b_pcb = new TBox( - Si_centroids[POSITION][5] + LAST_CENTROID_TO_PCB - PCB_LENGTH, -PCB_WIDTH/2, - Si_centroids[POSITION][5] + LAST_CENTROID_TO_PCB, PCB_WIDTH/2 );
	b_pcb->SetFillColor(pcb_green_i);
	b_pcb->Draw("SAME");

	// Draw 4-jaw slits
	TBox* b_four_jaw = new TBox( b_pcb->GetX2(), -ARR_DIAM/2, b_pcb->GetX2() + SLIT_LENGTH, ARR_DIAM/2 );
	b_four_jaw->SetFillColor( fj_red_i );
	b_four_jaw->Draw("SAME");

	// Draw Si strips
	TBox* b_Si[6];
	for ( Int_t a = 0; a < 6; a++ ){
		b_Si[a] = new TBox( -SI_WIDTH/2 - Si_centroids[POSITION][a], -SI_HEIGHT/2, SI_WIDTH/2 - Si_centroids[POSITION][a], SI_HEIGHT/2 );
		b_Si[a]->SetFillColor(si_strip_i);
		b_Si[a]->Draw("SAME");
	}

	// Calculate the required z
	Double_t z_fj = b_pcb->GetX2() + SLIT_LENGTH;

	
	// CALCULATIONS ---------------------------------------------------------------------------- //
	// Declare variables
	Double_t* theta_cm = new Double_t[NUM_THETA];
	Double_t* r_fj = new Double_t[NUM_THETA];
	Double_t** r = new Double_t*[NUM_THETA];
	Double_t** z = new Double_t*[NUM_THETA];
	for ( Int_t i = 0; i < NUM_Z; i++ ){
		r[i] = new Double_t[NUM_Z];
		z[i] = new Double_t[NUM_Z];
	}
	
	// Declare doubles for angle-specific kinematic terms
	Double_t T_lab_3, v_para, v_perp, theta_var;

	// Declare flag to stop trajectory
	Bool_t stop_bool = 0;
	
	// Declare TGraph's to hold trajectory information
	TGraph *g_traj[NUM_THETA];
	
	// Loop over centre-of-mass angle (theta)
	for ( Int_t i = 0; i < NUM_THETA; i++ ){
		// Calculate angle-dependent kinematic terms
		theta_cm[i] = ( i*theta_spacing + theta_lb );
		T_lab_3 = T_cm_3 + 0.5*mass[2]*V_cm*V_cm + mass[2]*v3*V_cm*TMath::Cos( theta_cm[i]*TMath::DegToRad() );	// [MeV]
		theta_var = 180 - theta_cm[i];									// The true CM angle
		v_para = V_cm + v3*TMath::Cos( theta_var*TMath::DegToRad() );													// [c]
		v_perp = v3*TMath::Sin( theta_var*TMath::DegToRad() );														// [c]

		// Calculate the radius at the required z
		r_fj[i] = RadiusCalculator( v_para, v_perp, cyclotron_freq, z_fj );

		// Set stop bool
		stop_bool = 0;
		
		// Calculate the trajectory
		for ( Int_t j = 0; j < NUM_Z; j++ ){
			// Calculate trajectory if it hasn't stopped
			if ( stop_bool == 0 ){
				z[i][j] = -z_spacing*j;
				r[i][j] = RadiusCalculator( v_para, v_perp, cyclotron_freq, z[i][j] );
				
				// Check if it has stopped
				if ( r[i][j] < ARR_DIAM/2 && z[i][j] <= z_fj ){
					stop_bool = 1;
				}
			}
			else{
				r[i][j] = r[i][j-1];
				z[i][j] = z[i][j-1];
			}
			
		}
		
		// Create graph of trajectory
		g_traj[i] = new TGraph( NUM_Z, z[i], r[i] );

		// Format graph of trajectory
		if ( r_fj[i] < ARR_DIAM/2 ){
			g_traj[i]->SetLineColor(kRed);
		}
		else if ( r_fj[i] >= ARR_DIAM/2 && r_fj[i] < TMath::Sqrt(2)*ARR_DIAM/2 ){
			g_traj[i]->SetLineColor(kOrange);
		}
		else{
			g_traj[i]->SetLineColor(kGreen);
		}
		
		if ( i % 10 == 0 ){
			g_traj[i]->Draw("SAME");
			/*std::cout << theta_cm[i] << "\n";
			for ( Int_t j = 0; j < NUM_Z; j++ ){
				if ( j % 10 == 0 ){
					std::cout << "\t" << z[i][j] << "\t" << r[i][j] << "\n";
				}
			}*/
		}

	}

	// DRAW THE ARRAY HEAD ON
	TCanvas* c1 = new TCanvas( "c_clip", "Clipping on the array", CANVAS_HEIGHT, CANVAS_HEIGHT );
	c1->cd();
	c1->DrawFrame( -10, -10, 10, 10 );

	// Draw the four jaws and the inside of the array
	TBox* b_fj_head = new TBox( -ARR_DIAM/2, -ARR_DIAM/2, ARR_DIAM/2, ARR_DIAM/2 );
	b_fj_head->SetFillColor( fj_red_i );
	b_fj_head->Draw("SAME");

	TBox* b_fj_inside = new TBox( -ARR_IN_DIAM/2, -ARR_IN_DIAM/2, ARR_IN_DIAM/2, ARR_IN_DIAM/2 );
	b_fj_inside->SetFillColor( kBlack );
	b_fj_inside->Draw("SAME");

	// Draw circles of constant radius ( does not equal orbits! )
	TEllipse* r_fj_ellipse[NUM_THETA];
	Double_t* phi_frac = new Double_t[NUM_THETA];
	
	for ( Int_t i = 0; i < NUM_THETA; i++ ){
		// Define ellipse
		r_fj_ellipse[i] = new TEllipse( 0, 0, r_fj[i], r_fj[i] );

		// Set colour based on whether it hits or not
		if ( r_fj[i] >= TMath::Sqrt(2.0)*ARR_DIAM/2 ){
			r_fj_ellipse[i]->SetLineColor( kGreen );
			phi_frac[i] = 1;
		}
		else if ( r_fj[i] < TMath::Sqrt(2.0)*ARR_DIAM/2 && r_fj[i] >= ARR_DIAM/2 ){
			r_fj_ellipse[i]->SetLineColor( kOrange );
			phi_frac[i] = 4*TMath::ACos( ARR_DIAM/( 2*r_fj[i] ) )/TMath::Pi();
		}
		else{
			r_fj_ellipse[i]->SetLineColor( kRed );
			phi_frac[i] = 0;
		}

		r_fj_ellipse[i]->SetFillColorAlpha( kBlack, 0 );


		if ( i % 4 ){
			r_fj_ellipse[i]->Draw("SAME");
		}

	}


	// Draw the phi fraction
	TCanvas* c2 = new TCanvas( "c_phi_frac", "Fraction of phi obstructed", CANVAS_WIDTH, CANVAS_HEIGHT );
	c2->cd();
	TGraph* g_phi_frac = new TGraph( NUM_THETA, theta_cm, phi_frac );
	g_phi_frac->SetTitle("");
	g_phi_frac->GetXaxis()->SetTitle( "#theta_{cm} / ^{#circ}" );
	g_phi_frac->GetYaxis()->SetTitle( "f_{#phi}" );
	g_phi_frac->SetMarkerStyle(20);
	g_phi_frac->SetMarkerSize(0.5);
	g_phi_frac->Draw("ALP");


	// CALCULATE THE MINIMUM THETA WHERE CLIPPING DOES NOT OCCUR

	// Print the canvases
	c0->Print("array_side.pdf");
	c1->Print("array_head_on.pdf");
	c2->Print("phi_frac.pdf");
}

































