// ArrayGeometry.C
// Checks the excitation of each state and sees where it starts getting clipped
// ============================================================================================= //
// Patrick MacGregor
// Nuclear Physics Research Group
// Department of Physics and Astronomy
// The University of Manchester
// ============================================================================================= //
#include "AG_constants.h"
#include "AG_functions.h"
#include "AG_style.h"

#include <TCanvas.h>
#include <TGraph.h>
#include <TPad.h>
#include <TH1.h>
#include <TH2.h>
#include <TBox.h>
#include <TEllipse.h>
#include <TLine.h>
#include <TRandom.h>
#include <TFile.h>

#include <iostream>
#include <fstream>

TFile* outfile;

void ArrayGeometryEX( Double_t ex = EX, Double_t theta_lb = THETA_LB, Double_t theta_ub = THETA_UB ){
	// Batch mode
	if ( BATCH_MODE ){ gROOT->SetBatch(kTRUE); }

	// Define the plotting conditions
	const Double_t THETA_HEAD = 15.0;//0.5*(theta_ub + theta_lb);
	const Int_t NUM_THETA = (Int_t)( ( theta_ub - theta_lb )/theta_spacing ) + 1;
	const Int_t MOD_SIDE_TRAJECTORIES = (Int_t)( NUM_SIDE_ANGLE_INCREMENT/theta_spacing );
	const Int_t NUM_THETA_GRAPHS = (Int_t)TMath::Ceil( NUM_THETA/MOD_SIDE_TRAJECTORIES );
	Int_t side_traj_ctr = 0;


	if ( THETA_HEAD < theta_lb || THETA_HEAD > theta_ub ){
		std::cout << "Trajectories will not plot as " << THETA_HEAD << " is not within sampled angles!" << "\n";
	}


	// SETUP ----------------------------------------------------------------------------------- //
	// Declare arrays for storing things
	Double_t* theta_cm = new Double_t[NUM_THETA];
	Double_t* theta_lab = new Double_t[NUM_THETA];
	Double_t* theta_graph = new Double_t[NUM_THETA_GRAPHS];

	// Store trajectories in x,y,z for proton[z coordinate]
	Double_t* x3 = new Double_t[NUM_ZP];
	Double_t* y3 = new Double_t[NUM_ZP];
	Double_t* z3 = new Double_t[NUM_ZP];

	// Store trajectories in x,y,z for recoils[z coordinate]
	Double_t* x4 = new Double_t[NUM_ZR];
	Double_t* y4 = new Double_t[NUM_ZR];
	Double_t* z4 = new Double_t[NUM_ZR];
	Double_t* r4 = new Double_t[NUM_ZR];

	// Interpolating gradients and intercepts
	Double_t mx = 0.0;
	Double_t my = 0.0;
	Double_t cx = 0.0;
	Double_t cy = 0.0;

	// Calculate the number of successful proton hits [theta]
	Int_t* num_proton_hits = new Int_t[NUM_THETA];

	// Store the fraction of theta where there is a proton-recoil coincidence [theta]
	Double_t* theta_frac_proton_recoil = new Double_t[NUM_THETA];

	// Store solid angles [theta value][type: 1 = hit, 2 = hit array, miss detector, 3 = clip edge of array]
	Double_t** d_phi = new Double_t*[NUM_THETA];
	for ( Int_t i = 0; i < NUM_THETA; i++ ){
		d_phi[i] = new Double_t[3];
	}

	// Solid angle pointer
	Double_t solid_angle[NUM_THETA];

	// Status indicators: 0 = fail, 1 = success, 2 = questionable
	Int_t status_ejectile = 0;
	Int_t status_recoil = 0;
	Int_t success_green = (Int_t)kGreen;
	Int_t fail_red = (Int_t)kRed;

	// Open log file
	std::ofstream log_file;
	log_file.open( Form( "output-data/EX_%s-POS_%i.dat", DecimalDotToUnderscore( ex ).Data(), POSITION ), std::ofstream::out );

	// Check file is open
	if ( !log_file.is_open() ){
		std::cout << "LOG FILE FAILED TO OPEN!\n";
		std::exit(1);
	}

	// Open file to record successful angles
	std::ofstream angle_status_file;
	angle_status_file.open( Form( "output-data/EX_%s-POS_%i_angle_status.dat", DecimalDotToUnderscore( ex ).Data(), POSITION ), std::ofstream::out );

	// Check file is open
	if ( !angle_status_file.is_open() ){
		std::cout << "ANGLE STATUS FILE FAILED TO OPEN!\n";
		std::exit(1);
	}
	else{
		angle_status_file <<
			"    \u03b8CM (\u00b0)" <<
			"     EJECT." <<
			"     RECOIL" << "\n";
	}

	// Percentage progress holder
	Double_t prog = 0.0;

	// Boolean to test if the starting position is correct
	Bool_t STARTS_IN_VALID_SPOT = 0;




	// CALCULATIONS ---------------------------------------------------------------------------- //
	// Calculate angle-independent terms
	Double_t T_cm_f = T_cm_i + Q - ex;							// Final energy in CM [MeV]
	Double_t T_cm_3 = mass[3]*T_cm_f/( mass[2] + mass[3] );		// Final energy of ejectile in CM [MeV]
	Double_t v3 = TMath::Sqrt( 2*T_cm_3/mass[2] );				// Final velocity of ejectile in CM [c]
	Double_t v4 = -mass[2]*v3/mass[3];							// Final velocity of recoil in CM [c]

	// Declare doubles for angle-specific kinematic terms
	Double_t T_lab_3, u3_para, u3_perp, theta_var, phi, beam_offset_r, beam_offset_theta;
	Double_t v4_para, v4_perp;
	TRandom rand_phi, rand_beam_offset_r, rand_beam_offset_theta;
	rand_phi.SetSeed(0);
	rand_beam_offset_r.SetSeed(1);
	rand_beam_offset_theta.SetSeed(2);

	// Calculate contributions to solid angle (roughly)
	Double_t phi_cont = TMath::TwoPi()*1.0/NUM_EVENTS_PER_THETA;




	// CANVAS ITEMS ---------------------------------------------------------------------------- //
	// Declare TGraph's to hold trajectory information
	TGraph *g_traj_ejectile_side[NUM_THETA_GRAPHS];
	for ( Int_t i = 0; i < NUM_THETA_GRAPHS; i++ ){ g_traj_ejectile_side[i] = NULL; }
	TGraph *g_traj_ejectile_head[NUM_TRAJECTORIES];
	for ( Int_t i = 0; i < NUM_TRAJECTORIES; i++ ){ g_traj_ejectile_head[i] = NULL; }

	Bool_t* b_ejectile_success = new Bool_t[NUM_TRAJECTORIES];
	for ( Int_t i = 0; i < NUM_TRAJECTORIES; i++ ){
		b_ejectile_success[i] = 0;
	}

	TGraph *g_traj_recoil_side[NUM_THETA_GRAPHS];
	for ( Int_t i = 0; i < NUM_THETA_GRAPHS; i++ ){ g_traj_recoil_side[i] = NULL; }
	TGraph *g_traj_recoil_head[NUM_TRAJECTORIES];
	for ( Int_t i = 0; i < NUM_TRAJECTORIES; i++ ){ g_traj_recoil_head[i] = NULL; }

	// Histogram to look at distribution of phi
	TH1D *h_phi_dist = new TH1D( "h_phi_dist", "Phi distribution", 360, 0, 360 );
	TH2D *h_beam_width = new TH2D( "h_beam_width", "Beam width", 400, -3*BEAM_FWHM, 3*BEAM_FWHM, 400, -3*BEAM_FWHM, 3*BEAM_FWHM );

	// TGraph to look at theta acceptance on recoil detector
	TGraph* g_recoil_theta_frac = NULL;

	// Define a TEllipse to hold the start points of each trajectory
	//TEllipse *beam_start[NUM_EVENTS_PER_THETA];
	//for ( Int_t i = 0; i < NUM_EVENTS_PER_THETA; i++ ){ beam_start[i] = NULL; }


	// LOOP OVER CENTRE-OF-MASS ANGLE (THETA) (i) ---------------------------------------------- //
	for ( Int_t i = 0; i < NUM_THETA; i++ ){
		// Calculate angle-dependent kinematic terms
		theta_cm[i] = ( i*theta_spacing + theta_lb );					// The convenient CM angle (labelled eta in thesis)
		theta_var = 180 - theta_cm[i];									// The "proper" CM angle
		T_lab_3 = T_cm_3 + 0.5*mass[2]*V_cm*V_cm + mass[2]*v3*V_cm*TMath::Cos( theta_var*TMath::DegToRad() );		// [MeV]
		u3_para = v3*TMath::Cos( theta_var*TMath::DegToRad() ) + V_cm;												// LAB [c]
		u3_perp = v3*TMath::Sin( theta_var*TMath::DegToRad() );
		if ( u3_para <= 0 ){											// LAB [c]
			theta_lab[i] = TMath::ASin( u3_perp/TMath::Sqrt( u3_para*u3_para + u3_perp*u3_perp ) )*TMath::RadToDeg();	// [DEG]
		}
		else{
			theta_lab[i] = ( TMath::Pi() - TMath::ASin( u3_perp/TMath::Sqrt( u3_para*u3_para + u3_perp*u3_perp ) ) )*TMath::RadToDeg();	// [DEG]
		}

		// Calculate recoil quantities
		v4_para = v4*TMath::Cos( theta_var*TMath::DegToRad() ) + V_cm;												// LAB [c]
		v4_perp = v4*TMath::Sin( theta_var*TMath::DegToRad() );														// LAB [c]

		// Define x and y success bools, as well as a stop flag for phi
		Bool_t x_success, y_success, x_failure, y_failure, z_failure, on_z_strip;
		Bool_t x_success2, y_success2, x_failure2, y_failure2;
		Bool_t z_strip[6];
		Bool_t phi_stop = 0;
		Int_t traj_colour[NUM_EVENTS_PER_THETA + 1];
		for ( Int_t j = 0; j < NUM_EVENTS_PER_THETA + 1; j++ ){ traj_colour[j] = (Int_t)kRed; }

		// Initialise phi solid angle array parameters
		d_phi[i][0] = 0.0;
		d_phi[i][1] = 0.0;
		d_phi[i][2] = 0.0;

		// Initialise counters for how many successful hits there are
		num_proton_hits[i] = 0;
		theta_frac_proton_recoil[i] = 0.0;


		// LOOP OVER NUMBER OF EVENTS FOR A GIVEN THETA (j) ------------------------------------- //
		for ( Int_t j = 0; j < NUM_EVENTS_PER_THETA + 1; j++ ){
		// Calculate trajectories for a given theta
			STARTS_IN_VALID_SPOT = 1;
			if ( j == 0 ){
				// Choose non-random values for plotting purposes
				phi = CHOSEN_PHI;
				beam_offset_r = 0.0;
				beam_offset_theta = 0.0;
			}
			else{
				phi = rand_phi.Uniform(0,360);
				beam_offset_r = rand_beam_offset_r.Gaus(0, BEAM_FWHM*FWHMToSigma() );
				beam_offset_theta = rand_beam_offset_theta.Uniform(0,180);
			}			// This is initial phi angle

			phi_stop = 0;


			// LOOP OVER Z (PROTONS) (k) ------------------------------------------------------- //
			for ( Int_t k = 0; k < NUM_ZP; k++ ){
				// Calculate x,y,z
				if ( phi_stop == 0 ){
					if ( theta_lab[i] <= 90.0 ){
						z3[k] = -1.0*z_spacing*k;
					}
					else{
						z3[k] = z_spacing*k;
					}
					x3[k] = u3_perp*UC_CToMS()*UC_MToCM()*( TMath::Cos( phi*TMath::DegToRad() ) - TMath::Cos( phi*TMath::DegToRad() + cyclotron_freq*z3[k]/( u3_para*UC_CToMS()*UC_MToCM() ) ) )/cyclotron_freq + beam_offset_r*TMath::Cos( beam_offset_theta*TMath::DegToRad() ) + BEAM_SPOT_OFF_X;
					y3[k] = u3_perp*UC_CToMS()*UC_MToCM()*( -TMath::Sin( phi*TMath::DegToRad() ) + TMath::Sin( phi*TMath::DegToRad() + cyclotron_freq*z3[k]/( u3_para*UC_CToMS()*UC_MToCM() ) ) )/cyclotron_freq + beam_offset_r*TMath::Sin( beam_offset_theta*TMath::DegToRad() ) + BEAM_SPOT_OFF_Y;

					// Check if trajectory passed through the array
					if ( k == 0 ){
						if ( TMath::Abs( x3[k] ) > ARR_IN_DIAM/2 || TMath::Abs( y3[k] ) > ARR_IN_DIAM/2 ){
							// Check that beam offset is not too high
							if ( TMath::Abs( BEAM_SPOT_OFF_X ) > ARR_IN_DIAM/2 || TMath::Abs( BEAM_SPOT_OFF_Y ) > ARR_IN_DIAM/2 ){
								std::cerr << "Centre of beam spot does not pass through array ever! EXIT\n";
								std::exit(1);
							}

							// Outside defined region.
							STARTS_IN_VALID_SPOT = 0;
						}
						else{
							if ( j > 0 ){
								h_phi_dist->Fill(phi);
								h_beam_width->Fill( x3[k],y3[k] );
							}
						}
					}
				}
				else if ( phi_stop == 1 && k > 0 ){
					x3[k] = x3[k-1];
					y3[k] = y3[k-1];
					z3[k] = z3[k-1];
				}

				// Calculate success/failures if it's still moving but now in z-range of array
				if ( phi_stop == 0 && STARTS_IN_VALID_SPOT == 1 && z3[k] <= z_fj ){
					// Define bools to say it's hit the array, and define success and failure -- interpolate lines
					x_success = 0;	y_success = 0;
					x_failure = 0;	y_failure = 0;

					// Find when it is in the array
					if ( TMath::Abs( x3[k] ) <= ARR_DIAM/2 && TMath::Abs( y3[k] ) <= ARR_DIAM/2 && ( TMath::Abs( x3[k-1] ) > ARR_DIAM/2 || TMath::Abs( y3[k-1] ) > ARR_DIAM/2 ) ){
						// Now find interpolation parameters for x and y lines
						GetInterpolation( x3[k-1], y3[k-1], x3[k], y3[k], mx, cx );
						GetInterpolation( y3[k-1], x3[k-1], y3[k], x3[k], my, cy );

						// Define temporary x's and y's along boundaries
						Double_t tempy = mx*TMath::Sign( ARR_DIAM/2, x3[k-1] ) + cx;
						Double_t tempx = my*TMath::Sign( ARR_DIAM/2, y3[k-1] ) + cy;

						// Test whether points lie on a straight line in order
						Bool_t x_straight = ( ( tempx - x3[k-1] > 0 && x3[k] - tempx > 0 ) || ( tempx - x3[k-1] < 0 && x3[k] - tempx < 0 ) );
						Bool_t y_straight = ( ( tempy - y3[k-1] > 0 && y3[k] - tempy > 0 ) || ( tempy - y3[k-1] < 0 && y3[k] - tempy < 0 ) );

						// Test whether the point lies on an edge
						Bool_t x_on_edge = ( TMath::Abs( tempx ) <= ARR_DIAM/2 );
						Bool_t y_on_edge = ( TMath::Abs( tempy ) <= ARR_DIAM/2 );

						if ( x_straight && x_on_edge ){
							if ( TMath::Abs( tempx ) < SI_HEIGHT/2 ){
								x_success = 1;
								y3[k] = TMath::Sign( ARR_DIAM/2, y3[k-1] );
								x3[k] = tempx;
							}
							else if ( TMath::Abs( tempx ) >= SI_HEIGHT/2 && TMath::Abs( tempx ) <= ARR_DIAM/2 ){
								x_failure = 1;
								y3[k] = TMath::Sign( ARR_DIAM/2, y3[k-1] );
								x3[k] = tempx;
							}
						}
						else if ( y_straight && y_on_edge ){
							if ( TMath::Abs( tempy ) < SI_HEIGHT/2 ){
								y_success = 1;
								x3[k] = TMath::Sign( ARR_DIAM/2, x3[k-1] );
								y3[k] = tempy;
							}
							else if ( TMath::Abs( tempy ) >= SI_HEIGHT/2 && TMath::Abs( tempy ) <= ARR_DIAM/2 ){
								y_failure = 1;
								x3[k] = TMath::Sign( ARR_DIAM/2, x3[k-1] );
								y3[k] = tempy;
							}
						}
						else{
							std::cout << std::setprecision(6) << x3[k-1] << "    " << tempx << "    " << x3[k] << "    " << y3[k-1] << "    " << tempy << "    " << y3[k] << "\n";
						}

					}

					// Define z failure - z fails if on 4 jaws
					z_failure = ( z3[k] >= z_fj - 0.5*z_spacing - SLIT_LENGTH && TMath::Abs( x3[k] ) <= ARR_DIAM/2 && TMath::Abs( y3[k] ) <= ARR_DIAM/2 );


					/*
					if ( k > 0 && TMath::Abs( x3[k-1] ) > ARR_DIAM/2 && TMath::Abs( x3[k] ) <= ARR_DIAM/2 && TMath::Abs( y3[k-1] ) <= ARR_DIAM/2 ){
						// X is within array -- calculate gradient and intercept
						GetInterpolation( x3[k-1], y3[k-1], x3[k], y3[k], mx, cx );
						x3[k] = TMath::Sign( ARR_DIAM/2, x3[k] );
						y3[k] = mx*x3[k] + cx;
						x_success = ( TMath::Abs( y3[k] ) < SI_HEIGHT/2 );
						x_failure = ( TMath::Abs( y3[k] ) >= SI_HEIGHT/2 );
					}
					else if ( k > 0 && TMath::Abs( y3[k-1] ) > ARR_DIAM/2 &&  TMath::Abs( y3[k] ) <= ARR_DIAM/2 && TMath::Abs( x3[k-1] ) <= ARR_DIAM/2 ){
						// Y is within array -- calculate gradient and intercept
						GetInterpolation( y3[k-1], x3[k-1], y3[k], x3[k], my, cy );
						y3[k] = TMath::Sign( ARR_DIAM/2, y3[k] );
						x3[k] = my*y3[k] + cy;
						y_success = ( TMath::Abs( x3[k] ) < SI_HEIGHT/2 );
						y_failure = ( TMath::Abs( x3[k] ) >= SI_HEIGHT/2 );
					}
					*/






					/*x_success2 = ( TMath::Abs( x3[k] ) <= ARR_DIAM/2 && TMath::Abs( y3[k] ) < SI_HEIGHT/2 );
					y_success2 = ( TMath::Abs( y3[k] ) <= ARR_DIAM/2 && TMath::Abs( x3[k] ) < SI_HEIGHT/2 );
					x_failure2 = ( TMath::Abs( x3[k] ) <= ARR_DIAM/2 && TMath::Abs( y3[k] ) > SI_HEIGHT/2 && TMath::Abs( y3[k] ) <= ARR_DIAM/2 );
					y_failure2 = ( TMath::Abs( y3[k] ) <= ARR_DIAM/2 && TMath::Abs( x3[k] ) > SI_HEIGHT/2 && TMath::Abs( x3[k] ) <= ARR_DIAM/2 );*/


					// It hits the array
					if ( x_success || y_success ){
						phi_stop = 1;

						// See if z fails
						if ( z_failure ){
							traj_colour[j] = (Int_t)kBlue;
						}
						else{
							// Define z landing on the strip
							on_z_strip = 0;
							for ( Int_t a = 0; a < 6; a++ ){
								z_strip[a] = ( z3[k] >= -Si_centroids[POSITION][a] - SI_WIDTH/2 && z3[k] <= -Si_centroids[POSITION][a] + SI_WIDTH/2 );
								on_z_strip = on_z_strip || z_strip[a];
							}

							// Test to see if it lands on a strip
							if ( on_z_strip ){
								traj_colour[j] = (Int_t)kGreen;
								if ( j > 0 ){
									d_phi[i][0] += phi_cont;
									num_proton_hits[i]++;
								}
							}
							// Fails to land on strip, but on array
							else{
								traj_colour[j] = (Int_t)kRed;
							}
						}

					}
					// It hits the array, but not on a strip
					else if ( x_failure || y_failure ){
						phi_stop = 1;
						traj_colour[j] = (Int_t)kRed;
						if ( j > 0 ){
							d_phi[i][1] += phi_cont;
						}
					}
					else if ( z_failure ){
						phi_stop = 1;
						traj_colour[j] = (Int_t)kBlue;
						if ( j > 0 ){
							d_phi[i][2] += phi_cont;
						}

					}

				}

				//std::cout << "\n" << NUM_ZP << "--" << i << "--" << j << "--" << k << "\n";
				//log_file << i << "--" << j << "--" << k << "\n";


			}	// Loop over z (protons) (k)

			//std::cout << "\n" << i << "--" << j << "\n";
			//log_file << i << "--" << j << "\n";

			// LOOP OVER Z (RECOILS) (l) ------------------------------------------------------- //
			for ( Int_t l = 0; l < NUM_ZR; l++ ){
				// Calculate x,y,z (and phi) --> looking the other way so x is reversed! Also, y has to change so it spins the other way
				if ( l < NUM_ZR - 1 ){
					z4[l] = 1.0*z_spacing*l;
				}
				else{
					z4[l] = TARGET_RDT_DISTANCE;
				}

				// Reverse x because of the change of perspective, and then reverse back for plotting purposes (x_recoil = -x_protons)
				x4[l] = (-1.0)*(-1.0)*( v4_perp*UC_CToMS()*UC_MToCM()*( TMath::Cos( phi*TMath::DegToRad() ) - TMath::Cos( phi*TMath::DegToRad() + cyclotron_freq_recoil*z4[l]/( v4_para*UC_CToMS()*UC_MToCM() ) ) )/cyclotron_freq_recoil + beam_offset_r*TMath::Cos( beam_offset_theta*TMath::DegToRad() ) + BEAM_SPOT_OFF_X );
				y4[l] = (-1.0)*( v4_perp*UC_CToMS()*UC_MToCM()*( -TMath::Sin( phi*TMath::DegToRad() ) + TMath::Sin( phi*TMath::DegToRad() + cyclotron_freq_recoil*z4[l]/( v4_para*UC_CToMS()*UC_MToCM() ) ) )/cyclotron_freq_recoil + beam_offset_r*TMath::Sin( beam_offset_theta*TMath::DegToRad() ) + BEAM_SPOT_OFF_Y );
				r4[l] = TMath::Sqrt( x4[l]*x4[l] + y4[l]*y4[l] );
			}	// Loop over z (recoils) (l)


			// CHECK TRAJECTORIES FOR VALIDITY AND PLOT GRAPHS --------------------------------- //
			if ( STARTS_IN_VALID_SPOT == 1 ){
				// Draw successful trajectories in head-on plots
				if ( theta_cm[i] == THETA_HEAD && j > 0 & j < NUM_TRAJECTORIES + 1 ){
					// protons
					g_traj_ejectile_head[j - 1] = new TGraph( NUM_ZP, x3, y3 );
					g_traj_ejectile_head[j - 1]->SetLineColorAlpha( traj_colour[j], 0.3 );
					if ( traj_colour[j] == kGreen ){ b_ejectile_success[j-1] = 1; }


					//beam_start[j] = new TEllipse( x3[0], y3[0], 0.01);
					//beam_start[j]->SetFillColor( traj_colour[j] );

					// recoils
					g_traj_recoil_head[j - 1] = new TGraph( NUM_ZR, x4, y4 );
					g_traj_recoil_head[j - 1]->SetLineColor( kBlack );

				}

				// Test how many recoils worked
				Double_t rtemp = TMath::Sqrt( y4[NUM_ZR-1]*y4[NUM_ZR-1] + x4[NUM_ZR-1]*x4[NUM_ZR-1] );
				if ( rtemp >= RDT_RADIUS_TO_CLEAR  && rtemp <= RDT_SI_OUTER_RAD + 0.5*TMath::Sqrt(2)*RDT_DETECTOR_GAP ){
					theta_frac_proton_recoil[i] += 100.0/NUM_EVENTS_PER_THETA;
				}

				// Draw trajectories in theta for a random phi
				//if ( j == EVENT_NUMBER && i % TMath::Max( (Int_t)TMath::Floor( (Double_t)NUM_THETA/50.0 ), 1 ) == 0 ){

				if ( j == EVENT_NUMBER && i % MOD_SIDE_TRAJECTORIES == 0 ){
					// Ejectiles
					theta_graph[side_traj_ctr] = theta_cm[i];
					g_traj_ejectile_side[ side_traj_ctr ] = new TGraph( NUM_ZP, z3, y3 );
					g_traj_ejectile_side[ side_traj_ctr ]->SetLineColorAlpha( traj_colour[j], 0.3 );
					g_traj_ejectile_side[ side_traj_ctr ]->SetLineWidth(1);

					// recoils
					g_traj_recoil_side[side_traj_ctr] = new TGraph( NUM_ZR, z4, r4 );
					if ( rtemp >= RDT_SI_INNER_RAD + 0.5*TMath::Sqrt(2)*RDT_DETECTOR_GAP && rtemp <= RDT_SI_OUTER_RAD + 0.5*TMath::Sqrt(2)*RDT_DETECTOR_GAP ){
						g_traj_recoil_side[side_traj_ctr]->SetLineColor(kGreen);
					}
					else if ( rtemp >= RDT_PCB_INNER_RAD + 0.5*TMath::Sqrt(2)*RDT_DETECTOR_GAP && rtemp < RDT_SI_INNER_RAD + 0.5*TMath::Sqrt(2)*RDT_DETECTOR_GAP ){
						g_traj_recoil_side[side_traj_ctr]->SetLineColor(kOrange);
					}
					else{
						g_traj_recoil_side[side_traj_ctr]->SetLineColor(kRed);
					}

					side_traj_ctr++;
				}
				// Status file CALCULATIONS
				if ( j == EVENT_NUMBER ){
					if( traj_colour[j] == success_green ){ status_ejectile = 1; }
					else if( traj_colour[j] == fail_red ){ status_ejectile = 0; }
					else{ status_ejectile = 2; }

					if ( rtemp >= RDT_SI_INNER_RAD + 0.5*TMath::Sqrt(2)*RDT_DETECTOR_GAP && rtemp <= RDT_SI_OUTER_RAD + 0.5*TMath::Sqrt(2)*RDT_DETECTOR_GAP ){
						status_recoil = 1;
					}
					else if ( rtemp >= RDT_PCB_INNER_RAD + 0.5*TMath::Sqrt(2)*RDT_DETECTOR_GAP && rtemp < RDT_SI_INNER_RAD + 0.5*TMath::Sqrt(2)*RDT_DETECTOR_GAP ){
						status_recoil = 2;
					}
					else{
						status_recoil = 0;
					}
				}
			}
			else{
				j--;
				continue;
			}



		}	// Loop over events per theta (j)


		// Point to solid angle (for graphing purposes)
		solid_angle[i] = d_phi[i][0];
		// Print to file
		if ( i == 0 && log_file.is_open() ){
			log_file << theta_lb << "\u00b0 <= \u03b8 <= " << theta_ub << "\u00b0; \u03b8-SPACING = " << theta_spacing << "\u00b0; NUM EVENTS PER THETA = " << NUM_EVENTS_PER_THETA << "; z-SPACING = " << z_spacing << " cm\n";
			log_file <<
				"     THETA_CM(\u00b0)" <<
				"    THETA_LAB(\u00b0)" <<
				"     X_FINAL(cm)" <<
				"     Y_FINAL(cm)" <<
				"     Z_FINAL(cm)" <<
				"        SA_1(sr)" <<
				"        SA_2(sr)" <<
				"        SA_3(sr)" <<
				"    THETA_RDT(%)\n";
		}

		// Standard line to file
		log_file << std::fixed << std::setw(log_width) << std::setprecision( TMath::Abs( TMath::Floor( TMath::Log10( theta_lb ) ) ) + TMath::Abs( TMath::Floor( TMath::Log10( theta_spacing ) ) ) + 1 ) << theta_cm[i] <<
								  std::setw(log_width) << std::setprecision(log_prec) << theta_lab[i] <<
								  std::setw(log_width) << std::setprecision(log_prec) << x3[NUM_ZP-1] <<
								  std::setw(log_width) << std::setprecision(log_prec) << y3[NUM_ZP-1] <<
								  std::setw(log_width) << std::setprecision(log_prec) << z3[NUM_ZP-1] <<
								  std::setw(log_width) << std::setprecision(log_prec) << d_phi[i][0] <<
								  std::setw(log_width) << std::setprecision(log_prec) << d_phi[i][1] <<
								  std::setw(log_width) << std::setprecision(log_prec) << d_phi[i][2] <<
								  std::setw(log_width) << std::setprecision(log_prec) << theta_frac_proton_recoil[i] << "\n";

		// Write to status files
		angle_status_file << std::fixed << std::setw(11) << std::setprecision(8) << theta_cm[i] <<
										   std::setw(11) << std::setprecision(8) << status_ejectile <<
										   std::setw(11) << std::setprecision(8) << status_recoil << "\n";


		// Output percentage in terminal
		prog = 100.0*(Double_t)i/(Double_t)NUM_THETA;
		std::cout << std::fixed << std::setprecision(2) << prog << "% complete\r" << std::flush;
		//std::cout << std::fixed << std::setprecision(4) << theta_cm[i] << "\u00b0 --> " << 100.0*num_proton_hits[i]/NUM_EVENTS_PER_THETA << "% x " << 100.0*theta_frac_proton_recoil[i]/NUM_EVENTS_PER_THETA << "% = " << 100.0*num_proton_hits[i]*theta_frac_proton_recoil[i]/( NUM_EVENTS_PER_THETA*NUM_EVENTS_PER_THETA ) << "%\n";


		// Calculate the radius at the rdt for a given theta
		//std::cout << theta_cm[i] << "\t" << RDT_RADIUS_TO_CLEAR - TMath::Abs( v4_perp*UC_CToMS()*UC_MToCM()*TMath::Sqrt( 2.0*( 1.0 - TMath::Cos( cyclotron_freq_recoil*TARGET_RDT_DISTANCE/( v4_para*UC_CToMS()*UC_MToCM() ) ) ) )/cyclotron_freq_recoil ) << " cm\n";

	}	// Loop over theta (i)



	// CALCULATE THE AVERAGE SOLID ANGLE ------------------------------------------------------- //
	Double_t solid_angle_mav[NUM_THETA];
	Int_t half_window_size = 4; // Multiple by 2*i + 1 to get window width

	for ( Int_t i = 0; i < NUM_THETA; i++ ){
		if ( i == half_window_size ){
			Double_t sum = 0;
			for ( Int_t j = 0; j < 2*half_window_size + 1; j++ ){
				sum += solid_angle[i + j - half_window_size];
			}
			solid_angle_mav[i] = sum/(2*half_window_size + 1);
		}
		else if ( i < NUM_THETA - half_window_size && i > half_window_size ){
			solid_angle_mav[i] = solid_angle_mav[i-1] + ( solid_angle[i+half_window_size] - solid_angle[i-1-half_window_size] )/(2*half_window_size + 1);
		}
		else{
			solid_angle_mav[i] = 0.0;
		}
	}




	// ***************************************************************************************** //
	//                                    NOW PLOT EVERYTHING                                    //
	// ***************************************************************************************** //
	// Check the output root file is open
	Bool_t bFileOpen = outfile->IsOpen();

	// PROTONS SIDE ON ------------------------------------------------------------------------- //
	// Set style
	CreateStyle( ptm_style );
	ptm_style->cd();
	TCanvas* c_traj_ejectile_side = new TCanvas( Form( "c_traj_ejectile_side_%s", DecimalDotToUnderscore( ex ).Data() ), "EJECTILE SIDE", CANVAS_WIDTH, CANVAS_HEIGHT );

	// Draw frame
	TH1F* frame = c_traj_ejectile_side->DrawFrame( ej_side_X1, ej_side_Y1, ej_side_X2, ej_side_Y2 );
	FormatFrame( frame, "z (cm)", "r (cm)" );

	// Draw target
	TBox* b_target = new TBox( 0, -2, 0.5, 2 );
	b_target->SetFillColor(kBlack);
	b_target->Draw("SAME");

	// Draw Al bracket for array
	TBox* b_al_bracket = new TBox( - Si_centroids[POSITION][5] + LAST_CENTROID_TO_PCB - PCB_LENGTH - 1, -ARR_DIAM/2, - Si_centroids[POSITION][5] + LAST_CENTROID_TO_PCB, ARR_DIAM/2 );
	b_al_bracket->SetFillColor(kGray);
	b_al_bracket->Draw("SAME");

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

	// Draw trajectories
	for ( Int_t i = 0; i < NUM_THETA_GRAPHS; i++ ){
		//if ( g_traj_ejectile_side[i] != NULL && theta_graph[i] >= 0.0 && theta_graph[i] <= 44.0 ){
		if ( g_traj_ejectile_side[i] != NULL && theta_graph[i] >=0 && theta_graph[i] <= THETA_DT_LIM ){
			g_traj_ejectile_side[i]->Draw("SAME");
		}
	}

	if ( bFileOpen ){
		outfile->cd(); c_traj_ejectile_side->Write();
	}
	std::cout << "Plotted ejectile_side" << "\n";

	// PROTONS HEAD ON ------------------------------------------------------------------------- //
	TCanvas* c_traj_ejectile_head = new TCanvas( Form( "c_traj_ejectile_head_%s", DecimalDotToUnderscore( ex ).Data() ), "EJECTILE HEAD", CANVAS_HEIGHT, CANVAS_HEIGHT );
	c_traj_ejectile_head->cd();
	ptm_style->cd();
	Int_t traj_ejectile_radius = (Int_t)TMath::Ceil( 2.05*MaxRadius( v3, THETA_HEAD, cyclotron_freq ) );
	TH1F* frame1 = c_traj_ejectile_head->DrawFrame( -traj_ejectile_radius, -traj_ejectile_radius, traj_ejectile_radius, traj_ejectile_radius );
	FormatFrame( frame1, "x_{p} (cm)", "y (cm)" );

	// Draw the four jaws and the inside of the array
	TBox* b_fj_head = new TBox( -ARR_DIAM/2, -ARR_DIAM/2, ARR_DIAM/2, ARR_DIAM/2 );
	b_fj_head->SetFillColor( fj_red_i );
	b_fj_head->Draw("SAME");

	TBox* b_fj_inside = new TBox( -ARR_IN_DIAM/2, -ARR_IN_DIAM/2, ARR_IN_DIAM/2, ARR_IN_DIAM/2 );
	b_fj_inside->SetFillColor( kBlack );
	b_fj_inside->Draw("SAME");

	// Draw Si detectors
	TBox *b_si_headon[4];
	b_si_headon[0] = new TBox( -SI_HEIGHT/2, -ARR_DIAM/2, SI_HEIGHT/2, -0.9*ARR_DIAM/2 );
	b_si_headon[1] = new TBox( -SI_HEIGHT/2, ARR_DIAM/2, SI_HEIGHT/2, 0.9*ARR_DIAM/2 );
	b_si_headon[2] = new TBox( -ARR_DIAM/2, -SI_HEIGHT/2, -0.9*ARR_DIAM/2, SI_HEIGHT/2 );
	b_si_headon[3] = new TBox( ARR_DIAM/2, -SI_HEIGHT/2, 0.9*ARR_DIAM/2, SI_HEIGHT/2 );
	for ( Int_t j = 0; j < 4; j++ ){
		b_si_headon[j]->SetFillColor(kBlack);
		b_si_headon[j]->Draw("SAME");
	}
	// Draw trajectories
	Int_t ctr_traj_ejectile_head = 0;
	for ( Int_t j = 0; j < NUM_TRAJECTORIES; j++ ){
		if ( g_traj_ejectile_head[j] != NULL ){ g_traj_ejectile_head[j]->Draw("C SAME"); }
		//beam_start[j]->Draw("SAME");
	}
	if ( bFileOpen ){
		outfile->cd(); c_traj_ejectile_head->Write();
	}

	std::cout << "Plotted ejectile_head" << "\n";
	// DRAW THE PHI FRACTION ------------------------------------------------------------------- //
	TCanvas* c_phi_frac = new TCanvas( Form( "c_phi_frac_%s", DecimalDotToUnderscore( ex ).Data() ), "PHI FRACTION", CANVAS_WIDTH, CANVAS_HEIGHT );
	c_phi_frac->cd();
	ptm_style->cd();

	TH1F* frame2 = c_phi_frac->DrawFrame( TMath::Max( 0.0, theta_lb - theta_spacing ), 0, theta_ub + theta_spacing, 6.4 );
	FormatFrame( frame2, "#theta_{cm} (^{#circ})", "#Delta#phi (rad)" );

	TGraph* g_phi_frac = new TGraph( NUM_THETA, theta_cm, solid_angle );
	TGraph* g_phi_frac_mav = new TGraph( NUM_THETA, theta_cm, solid_angle_mav );

	g_phi_frac->SetMarkerStyle(20);
	g_phi_frac->SetMarkerSize(0.5);
	g_phi_frac->Draw("LP");

	g_phi_frac_mav->SetMarkerStyle(20);
	g_phi_frac_mav->SetMarkerSize(0.5);
	g_phi_frac_mav->SetLineColor(kBlue);
	g_phi_frac_mav->SetMarkerColor(kBlue);
	g_phi_frac_mav->Draw("LP SAME");

	TLine *l_solid_angle = new TLine( c_phi_frac->GetUxmin(), 8*TMath::ATan( 9.0/23.0 ), c_phi_frac->GetUxmax(), 8*TMath::ATan( 9.0/23.0 ) );
	l_solid_angle->SetLineStyle(2);
	l_solid_angle->SetLineColor(kRed);
	l_solid_angle->SetLineWidth(2);
	l_solid_angle->Draw("SAME");

	if ( bFileOpen ){
		outfile->cd(); c_phi_frac->Write();
	}

	std::cout << "Plotted phi_frac" << "\n";

	// DRAW THE RANDOM NUMBER DISTRIBUTIONS ---------------------------------------------------- //
	TCanvas* c_rand_numbers = new TCanvas( Form( "c_rand_numbers_%s", DecimalDotToUnderscore( ex ).Data() ), "RANDOM NUMBERS", CANVAS_WIDTH, 0.6465*CANVAS_HEIGHT);
	c_rand_numbers->Divide(2,1);

	c_rand_numbers->cd(1);
	ptm_style->cd();
	TPad* pad = (TPad*)c_rand_numbers->GetPad(1);
	pad->SetLeftMargin(0.13);
	pad->SetRightMargin(0.01);
	pad->SetBottomMargin(0.09);
	pad->SetTopMargin(0.02);
	FormatFrame( h_phi_dist, "#phi (#circ)", "Counts" );
	h_phi_dist->GetYaxis()->SetTitleOffset(2.0);
	h_phi_dist->GetYaxis()->SetRangeUser( 0, 1.05*h_phi_dist->GetMaximum() );
	h_phi_dist->Draw();

	c_rand_numbers->cd(2);
	ptm_style->cd();
	pad = (TPad*)c_rand_numbers->GetPad(2);
	pad->SetLeftMargin(0.09);
	pad->SetRightMargin(0.13);
	pad->SetBottomMargin(0.09);
	pad->SetTopMargin(0.02);

	FormatFrame( h_beam_width, "x (cm)", "y (cm)" );
	h_beam_width->Draw("colz");
	h_beam_width->GetYaxis()->SetTitleOffset(1.3);

	TLine* border[4];
	border[0] = new TLine( -ARR_IN_DIAM/2, -ARR_IN_DIAM/2, -ARR_IN_DIAM/2, ARR_IN_DIAM/2 );
	border[1] = new TLine( ARR_IN_DIAM/2, -ARR_IN_DIAM/2, ARR_IN_DIAM/2, ARR_IN_DIAM/2 );
	border[2] = new TLine( -ARR_IN_DIAM/2, -ARR_IN_DIAM/2, ARR_IN_DIAM/2, -ARR_IN_DIAM/2 );
	border[3] = new TLine( -ARR_IN_DIAM/2, ARR_IN_DIAM/2, ARR_IN_DIAM/2, ARR_IN_DIAM/2 );

	for ( Int_t i = 0; i < 4; i++ ){
		border[i]->SetLineWidth(1);
		border[i]->SetLineColor(fj_red_i);
		border[i]->Draw("SAME");
	}

	if ( bFileOpen ){
		outfile->cd(); c_rand_numbers->Write();
	}

	std::cout << "Plotted rand_numbers" << "\n";

	// RECOILS SIDE ON ------------------------------------------------------------------------- //
	// Set style
	CreateStyle( ptm_style );
	ptm_style->cd();
	TCanvas* c_traj_recoil_side = new TCanvas( Form( "c_traj_recoil_side_%s", DecimalDotToUnderscore( ex ).Data() ), "RECOIL SIDE", CANVAS_WIDTH, CANVAS_HEIGHT );

	// Draw frame
	TH1F* frame_rec_side = c_traj_recoil_side->DrawFrame( rec_side_X1, rec_side_Y1, rec_side_X2, rec_side_Y2 );
	FormatFrame( frame_rec_side, "z (cm)", "r (cm)" );

	// Draw target
	TBox* b_rec_target = new TBox( -2.0, -ARR_IN_DIAM/2, 0.0, ARR_IN_DIAM/2 );
	b_rec_target->SetFillColor(kBlack);
	b_rec_target->Draw("SAME");

	// Draw detector
	TBox* b_rec_Si_side[2];
	TBox* b_rec_pcb_side[2];

	for ( Int_t i = 0; i < 2; i++ ){
		b_rec_pcb_side[i] = new TBox( TARGET_RDT_DISTANCE, TMath::Power( -1, i )*( RDT_PCB_INNER_RAD + 0.5*RDT_DETECTOR_GAP ), TARGET_RDT_DISTANCE + 2.0, TMath::Power( -1, i )*( RDT_PCB_OUTER_RAD + 0.5*RDT_DETECTOR_GAP ) );
		b_rec_pcb_side[i]->SetFillColor( pcb_green_i );
		b_rec_pcb_side[i]->SetLineWidth( 0 );
		b_rec_pcb_side[i]->Draw("SAME");

		b_rec_Si_side[i] = new TBox( TARGET_RDT_DISTANCE, TMath::Power( -1, i )*( RDT_SI_INNER_RAD + 0.5*RDT_DETECTOR_GAP ), TARGET_RDT_DISTANCE + 1.0, TMath::Power( -1, i )*( RDT_SI_OUTER_RAD + 0.5*RDT_DETECTOR_GAP ) );
		b_rec_Si_side[i]->SetFillColor( si_strip_i );
		b_rec_Si_side[i]->SetLineWidth( 0 );
		b_rec_Si_side[i]->Draw("SAME");
	}

	// Draw trajectories
	Int_t ctr_traj_recoil_side = 0;
	for ( Int_t i = 0; i < NUM_THETA_GRAPHS; i++ ){
		if ( g_traj_recoil_side[i] != NULL && i % 4 == 0 ){
			g_traj_recoil_side[i]->Draw("SAME");
		}
	}

	if ( bFileOpen ){
		outfile->cd(); c_traj_recoil_side->Write();
	}

	std::cout << "Plotted recoil_side" << "\n";

	// RECOILS HEAD ON ------------------------------------------------------------------------- //
	TCanvas* c_traj_recoil_head = new TCanvas( Form( "c_traj_recoil_head_%s", DecimalDotToUnderscore( ex ).Data() ), "RECOIL HEAD", CANVAS_HEIGHT, CANVAS_HEIGHT );
	TPad* pad_traj_recoil_head = (TPad*)c_traj_recoil_head;
	c_traj_recoil_head->cd();
	ptm_style->cd();
	TH1F* frame_rec_head = c_traj_recoil_head->DrawFrame( -( 1.2*RDT_PCB_OUTER_RAD ), -( 1.2*RDT_PCB_OUTER_RAD ), ( 1.2*RDT_PCB_OUTER_RAD ), ( 1.2*RDT_PCB_OUTER_RAD ) );
	FormatFrame( frame_rec_head, "x_{rec} (cm)", "y (cm)" );

	// Draw the recoil detector
	TEllipse* rdt_pcb[4];
	TEllipse* rdt_Si[4];
	TEllipse* rdt_Si_inner_pcb[4];
	TEllipse* rdt_hole;

	for ( Int_t i = 0; i < 4; i++ ){
		// Draw large PCBs
		rdt_pcb[i] = new TEllipse(
			0.5*RDT_DETECTOR_GAP*TMath::Cos( RDT_ROTATION*TMath::DegToRad() + 0.25*TMath::Pi()*( 1 + 2*i ) ), 	// X1
			0.5*RDT_DETECTOR_GAP*TMath::Sin( RDT_ROTATION*TMath::DegToRad() + 0.25*TMath::Pi()*( 1 + 2*i ) ), 	// Y1
			RDT_PCB_OUTER_RAD, 														// R1
			RDT_PCB_OUTER_RAD, 														// R2
			i*90 + RDT_ROTATION, 													// Phi min
			(i+1)*90 + RDT_ROTATION 												// Phi max
		);
		rdt_pcb[i]->SetFillColor( pcb_green_i );
		rdt_pcb[i]->SetLineWidth(0);
		rdt_pcb[i]->Draw("SAME");



		// Draw Si detectors
		rdt_Si[i] = new TEllipse( 0.5*TMath::Sqrt(2)*RDT_DETECTOR_GAP*TMath::Cos( RDT_ROTATION*TMath::DegToRad() + 0.25*TMath::Pi()*( 1 + 2*i ) ), 0.5*TMath::Sqrt(2)*RDT_DETECTOR_GAP*TMath::Sin( RDT_ROTATION*TMath::DegToRad() + 0.25*TMath::Pi()*( 1 + 2*i ) ), RDT_SI_OUTER_RAD, RDT_SI_OUTER_RAD, i*90 + 0.8*( 90 - RDT_ANGULAR_COVERAGE ) + RDT_ROTATION, (i+1)*90 - 0.2*( 90 - RDT_ANGULAR_COVERAGE ) + RDT_ROTATION );
		rdt_Si[i]->SetFillColor( si_strip_i );
		rdt_Si[i]->SetLineWidth(0);
		rdt_Si[i]->Draw("SAME");

		// Draw inner PCBs
		rdt_Si_inner_pcb[i] = new TEllipse( 0.5*TMath::Sqrt(2)*RDT_DETECTOR_GAP*TMath::Cos( RDT_ROTATION*TMath::DegToRad() + 0.25*TMath::Pi()*( 1 + 2*i ) ), 0.5*TMath::Sqrt(2)*RDT_DETECTOR_GAP*TMath::Sin( RDT_ROTATION*TMath::DegToRad() + 0.25*TMath::Pi()*( 1 + 2*i ) ), RDT_SI_INNER_RAD, RDT_SI_INNER_RAD, i*90 + RDT_ROTATION, (i+1)*90 + RDT_ROTATION );
		rdt_Si_inner_pcb[i]->SetFillColor( pcb_green_i );
		rdt_Si_inner_pcb[i]->SetLineWidth(0);
		rdt_Si_inner_pcb[i]->Draw("SAME");
	}
	rdt_hole = new TEllipse( 0, 0, RDT_PCB_INNER_RAD + TMath::Sqrt(2)*0.5*RDT_DETECTOR_GAP, RDT_PCB_INNER_RAD + TMath::Sqrt(2)*0.5*RDT_DETECTOR_GAP, 0, 360 );
	rdt_hole->SetFillColor(kWhite);
	rdt_hole->SetLineWidth(0);
	rdt_hole->Draw("SAME");

	// Draw trajectories
	for ( Int_t j = 0; j < NUM_TRAJECTORIES; j++ ){
		if ( b_ejectile_success[j] == 1 ){
			// Calculate whether the trajectories hit the detector
			Double_t xtemp = g_traj_recoil_head[j]->GetX()[NUM_ZR-1];
			Double_t ytemp = g_traj_recoil_head[j]->GetY()[NUM_ZR-1];
			Double_t rtemp = TMath::Sqrt( xtemp*xtemp + ytemp*ytemp );

			if ( rtemp >= RDT_RADIUS_TO_CLEAR ){
				// Successful hit
				g_traj_recoil_head[j]->SetLineColorAlpha(kGreen, 0.5);
			}
			else if ( rtemp < RDT_RADIUS_TO_CLEAR && rtemp >= RDT_SI_INNER_RAD + 0.5*TMath::Sqrt(2)*RDT_DETECTOR_GAP ){
				// Hits PCB but not the Si
				g_traj_recoil_head[j]->SetLineColorAlpha(kOrange, 0.5);
			}
			else if ( rtemp < RDT_SI_INNER_RAD + 0.5*TMath::Sqrt(2)*RDT_DETECTOR_GAP ){
				// Misses both PCB and Si
				g_traj_recoil_head[j]->SetLineColorAlpha(kRed, 0.5);
			}
			else{
				// Does something else
				g_traj_recoil_head[j]->SetLineColorAlpha(kBlack,0.5);
				std::cout << xtemp << "\t" << ytemp << "\t" << rtemp << "\t" << "\n";
			}

			g_traj_recoil_head[j]->Draw("C SAME");
		}
	}

	if ( bFileOpen ){
		outfile->cd(); c_traj_recoil_head->Write();
	}

	std::cout << "Plotted recoil_head" << "\n";

	// FRACTION OF THETA ACCEPTED ON RECOIL DETECTORS ------------------------------------------ //
	TCanvas* c_recoil_theta_frac = new TCanvas( Form( "c_recoil_theta_frac_%s", DecimalDotToUnderscore( ex ).Data() ), "THETA FRACTION ON RDT", CANVAS_WIDTH, CANVAS_HEIGHT );
	c_recoil_theta_frac->cd();
	ptm_style->cd();

	TH1F* fr_recoil_theta_frac = c_recoil_theta_frac->DrawFrame( TMath::Max( 0.0, theta_lb - theta_spacing ), 0, theta_ub + theta_spacing, 105.0 );
	FormatFrame( fr_recoil_theta_frac, "#theta_{cm} (^{#circ})", "% hits on RDT" );

	g_recoil_theta_frac = new TGraph( NUM_THETA, theta_cm, theta_frac_proton_recoil );
	g_recoil_theta_frac->SetMarkerStyle(20);
	g_recoil_theta_frac->SetMarkerSize(0.4);
	g_recoil_theta_frac->GetXaxis()->CenterTitle();
	g_recoil_theta_frac->GetYaxis()->CenterTitle();
	g_recoil_theta_frac->SetMarkerColor( kRed );
	g_recoil_theta_frac->Draw("P");

	if ( bFileOpen ){
		outfile->cd(); c_recoil_theta_frac->Write();
	}


	// PRINT THE CANVASES AND FREE MEMORY------------------------------------------------------- //
	PrintCanvas( c_traj_ejectile_side, Form( "output-data/EX_%s-POS_%i-traj_ejectile_side", DecimalDotToUnderscore( ex ).Data(), POSITION ) );
	PrintCanvas( c_traj_ejectile_head, Form( "output-data/EX_%s-POS_%i-traj_ejectile_head", DecimalDotToUnderscore( ex ).Data(), POSITION ) );
	PrintCanvas( c_phi_frac,           Form( "output-data/EX_%s-POS_%i-phi_frac",           DecimalDotToUnderscore( ex ).Data(), POSITION ) );
	PrintCanvas( c_rand_numbers,       Form( "output-data/EX_%s-POS_%i-rand_numbers",       DecimalDotToUnderscore( ex ).Data(), POSITION ) );
	PrintCanvas( c_traj_recoil_side,   Form( "output-data/EX_%s-POS_%i-traj_recoil_side",   DecimalDotToUnderscore( ex ).Data(), POSITION ) );
	PrintCanvas( c_traj_recoil_head,   Form( "output-data/EX_%s-POS_%i-traj_recoil_head",   DecimalDotToUnderscore( ex ).Data(), POSITION ) );
	PrintCanvas( c_recoil_theta_frac,  Form( "output-data/EX_%s-POS_%i-recoil_theta_frac",  DecimalDotToUnderscore( ex ).Data(), POSITION ) );

	// Close the files
	log_file.close();
	angle_status_file.close();


	if ( !BATCH_MODE ){ gPad->WaitPrimitive("TPave"); }


	// Clear the memory
	// Arrays
	delete[] theta_cm;
	delete[] theta_graph;
	delete[] theta_lab;
	delete[] x3;
	delete[] y3;
	delete[] z3;
	delete[] x4;
	delete[] y4;
	delete[] z4;
	delete[] r4;
	delete[] num_proton_hits;
	delete[] theta_frac_proton_recoil;
	delete[] d_phi;

	std::cout << "Deleted arrays ... \n";

	// TCanvases
	if ( c_traj_ejectile_side->IsOnHeap() ){ delete c_traj_ejectile_side; }
	if ( c_traj_ejectile_head->IsOnHeap() ){ delete c_traj_ejectile_head; }
	if ( c_phi_frac->IsOnHeap() ){ delete c_phi_frac; }
	if ( c_rand_numbers->IsOnHeap() ){ delete c_rand_numbers; }
	if ( c_traj_recoil_side->IsOnHeap() ){ delete c_traj_recoil_side; }
	if ( c_traj_recoil_head->IsOnHeap() ){ delete c_traj_recoil_head; }
	if ( c_recoil_theta_frac->IsOnHeap() ){ delete c_recoil_theta_frac; }

	std::cout << "Deleted canvases ... \n";

	// TBoxes
	if ( b_target->IsOnHeap() ){ b_target->Delete(); }
	if ( b_pcb->IsOnHeap() ){ b_pcb->Delete(); }
	if ( b_four_jaw->IsOnHeap() ){ b_four_jaw->Delete(); }
	for ( Int_t i = 0; i < 6; i++ ){
		if ( b_Si[i]->IsOnHeap() ){ b_Si[i]->Delete(); }
	}
	if ( b_fj_head->IsOnHeap() ){ b_fj_head->Delete(); }
	if ( b_fj_inside->IsOnHeap() ){ b_fj_inside->Delete(); }
	for ( Int_t i = 0; i < 4; i++ ){
		if ( b_si_headon[i]->IsOnHeap() ){ b_si_headon[i]->Delete(); }
	}
	if ( b_rec_target->IsOnHeap() ){ b_rec_target->Delete(); }

	std::cout << "Deleted boxes ... \n";

	// TGraphs
	for ( Int_t i = 0; i < NUM_TRAJECTORIES; i++ ){
		if ( g_traj_ejectile_head[i] != NULL ){ if ( g_traj_ejectile_head[i]->IsOnHeap() ){ g_traj_ejectile_head[i]->Delete(); } }
		if ( g_traj_recoil_head[i] != NULL ){ if ( g_traj_recoil_head[i]->IsOnHeap() ){ g_traj_recoil_head[i]->Delete(); } }
	}
	for ( Int_t i = 0; i < NUM_THETA_GRAPHS; i++ ){
		if ( g_traj_ejectile_side[i] != NULL ){ if ( g_traj_ejectile_side[i]->IsOnHeap() ){ g_traj_ejectile_side[i]->Delete(); } }
		if ( g_traj_recoil_side[i] != NULL ){ if ( g_traj_recoil_side[i]->IsOnHeap() ){ g_traj_recoil_side[i]->Delete(); } }
	}
	if ( g_phi_frac->IsOnHeap() ){ g_phi_frac->Delete(); }
	if ( g_recoil_theta_frac->IsOnHeap() ){ g_recoil_theta_frac->Delete(); }

	std::cout << "Deleted graphs ... \n";

	// TH's

	if ( h_beam_width->IsOnHeap() ){ h_beam_width->Delete(); }
	if ( h_phi_dist->IsOnHeap() ){ h_phi_dist->Delete(); }
	//if ( frame->IsOnHeap() ){ frame->Delete(); }
	//if ( frame1->IsOnHeap() ){ frame1->Delete(); }
	//if ( frame2->IsOnHeap() ){ frame2->Delete(); }
	//if ( frame_rec_side->IsOnHeap() ){ frame_rec_side->Delete(); }
	//if ( frame_rec_head->IsOnHeap() ){ frame_rec_head->Delete(); }
	//if ( fr_recoil_theta_frac->IsOnHeap() ){ fr_recoil_theta_frac->Delete(); }

	std::cout << "Deleted histograms ... \n";

	// TLines
	if ( l_solid_angle->IsOnHeap() ){ l_solid_angle->Delete(); }
	for ( Int_t i = 0; i < 4; i++ ){
		if ( border[i]->IsOnHeap() ){ border[i]->Delete(); }
	}

	std::cout << "Deleted lines ... \n";

	// TEllipses
	//for ( Int_t i = 0; i < NUM_EVENTS_PER_THETA; i++ ){
		//if ( beam_start[i] != NULL ){ if ( beam_start[i]->IsOnHeap() ){ beam_start[i]->Delete(); } }
	//}


	// TPads
	//if ( pad->IsOnHeap() ){ pad->Delete(); }
	return;
}


// MAIN FUNCTION ------------------------------------------------------------------------------- //
void ArrayGeometry(){
	const Int_t NUM_STATES = 5;
	const Bool_t ALL = 0;
	/*Double_t STATES[NUM_STATES] = {
		0.00000,
		0.05460,
		1.09435,
		1.43465,
		2.24901,
		2.47816,
		2.87390,
		3.17944,
		3.85680,
		3.96590,
		//4.28395
		4.30000
	};*/

	Double_t theta_lb_arr[NUM_STATES] = { 16.3, 17.0, 17.6, 18.0, 19.8 };
	Double_t theta_ub_arr[NUM_STATES] = { 17.3, 18.0, 18.7, 19.0, 22.2 };

	Double_t STATES[NUM_STATES] = {
		0.06,
		1.45,
		2.50,
		3.20,
		6.15,
	};

	// Open ROOT file to write all the graphs etc
	outfile = new TFile( "output-data/array_geometry_files.root", "RECREATE" );


	// Plot all the things!
	if ( ALL == 1 ){
		for ( Int_t i = 0; i < NUM_STATES; i++ ){
			ArrayGeometryEX( STATES[i],  theta_lb_arr[i],  theta_ub_arr[i] );
			gROOT->ls();
		}
	}
	else{
		ArrayGeometryEX( 0.000 );
	}

	// Close the file
	if ( outfile->IsOpen() ){ outfile->Close(); }

}



/* *** TODO ***
 * Residuals for phi_frac plot
 * Error bars for phi_frac plot?
 * Save in ROOT files (in tree?) for ease of manipulation
 * Recoil solid angle stuff
*/
