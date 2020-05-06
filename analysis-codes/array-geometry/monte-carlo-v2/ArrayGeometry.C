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

void ArrayGeometryEX( Double_t ex = EX ){
	// Batch mode
	if ( BATCH_MODE ){ gROOT->SetBatch(kTRUE); }

	// Check the plotting conditions
	if ( THETA_HEAD < theta_lb || THETA_HEAD > theta_ub ){
		std::cout << "Trajectories will not plot as " << THETA_HEAD << " is not within sampled angles!" << "\n";
	}


	// SETUP ----------------------------------------------------------------------------------- //
	// Declare arrays for storing things
	Double_t* theta_cm = new Double_t[NUM_THETA];
	Double_t* theta_lab = new Double_t[NUM_THETA];

	// Store trajectories in x,y,z for proton[z coordinate]
	Double_t* x3 = new Double_t[NUM_ZP];
	Double_t* y3 = new Double_t[NUM_ZP];
	Double_t* z3 = new Double_t[NUM_ZP];

	// Store trajectories in x,y,z for recoils[z coordinate]
	Double_t* x4 = new Double_t[NUM_ZR];
	Double_t* y4 = new Double_t[NUM_ZR];
	Double_t* z4 = new Double_t[NUM_ZR];

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

	// Open log file
	std::ofstream log_file;
	log_file.open( Form( "output-data/EX_%s-POS_%i.dat", DecimalDotToUnderscore( ex ).Data(), POSITION ), std::ofstream::out );

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
	Double_t T_lab_3, v3_para, v3_perp, theta_var, phi, beam_offset_r, beam_offset_theta;
	Double_t v4_para, v4_perp;
	TRandom rand_phi, rand_beam_offset_r, rand_beam_offset_theta;
	rand_phi.SetSeed(0);
	rand_beam_offset_r.SetSeed(1);
	rand_beam_offset_theta.SetSeed(2);

	// Calculate contributions to solid angle (roughly)
	Double_t phi_cont = TMath::TwoPi()*1.0/NUM_EVENTS_PER_THETA;
	



	// CANVAS ITEMS ---------------------------------------------------------------------------- //
	// Declare TGraph's to hold trajectory information
	TGraph *g_traj_ejectile_side[NUM_THETA];
	for ( Int_t i = 0; i < NUM_THETA; i++ ){ g_traj_ejectile_side[i] = NULL; }
	TGraph *g_traj_ejectile_head[NUM_EVENTS_PER_THETA];
	for ( Int_t i = 0; i < NUM_EVENTS_PER_THETA; i++ ){ g_traj_ejectile_head[i] = NULL; }

	TGraph *g_traj_recoil_side[NUM_THETA];
	for ( Int_t i = 0; i < NUM_THETA; i++ ){ g_traj_recoil_side[i] = NULL; }
	TGraph *g_traj_recoil_head[NUM_EVENTS_PER_THETA];
	for ( Int_t i = 0; i < NUM_EVENTS_PER_THETA; i++ ){ g_traj_recoil_head[i] = NULL; }

	// Histogram to look at distribution of phi
	TH1D *h_phi_dist = new TH1D( "h_phi_dist", "Phi distribution", 360, 0, 360 );
	TH2D *h_beam_width = new TH2D( "h_beam_width", "Beam width", 400, -3*BEAM_FWHM, 3*BEAM_FWHM, 400, -3*BEAM_FWHM, 3*BEAM_FWHM );

	// TGraph to look at theta acceptance on recoil detector
	TGraph* g_recoil_theta_frac = NULL;

	// Define a TEllipse to hold the start points of each trajectory
	TEllipse *beam_start[NUM_EVENTS_PER_THETA];
	for ( Int_t i = 0; i < NUM_EVENTS_PER_THETA; i++ ){ beam_start[i] = NULL; }

	


	// LOOP OVER CENTRE-OF-MASS ANGLE (THETA) (i) ---------------------------------------------- //
	for ( Int_t i = 0; i < NUM_THETA; i++ ){
		// Calculate angle-dependent kinematic terms
		theta_cm[i] = ( i*theta_spacing + theta_lb );
		theta_var = 180 - theta_cm[i];									// The easier CM angle
		T_lab_3 = T_cm_3 + 0.5*mass[2]*V_cm*V_cm + mass[2]*v3*V_cm*TMath::Cos( theta_var*TMath::DegToRad() );		// [MeV]
		v3_para = v3*TMath::Cos( theta_var*TMath::DegToRad() ) + V_cm;												// LAB [c]
		v3_perp = v3*TMath::Sin( theta_var*TMath::DegToRad() );														// LAB [c]
		theta_lab[i] = TMath::ASin( v3_perp/TMath::Sqrt( v3_para*v3_para + v3_perp*v3_perp ) )*TMath::RadToDeg();	// [DEG]

		// Calculate recoil quantities
		v4_para = v4*TMath::Cos( theta_var*TMath::DegToRad() ) + V_cm;												// LAB [c]
		v4_perp = v4*TMath::Sin( theta_var*TMath::DegToRad() );														// LAB [c]

		// Define x and y success bools, as well as a stop flag for phi
		Bool_t x_success, y_success, x_failure, y_failure, z_failure;
		Bool_t phi_stop = 0;
		Int_t traj_colour[NUM_EVENTS_PER_THETA];
		for ( Int_t j = 0; j < NUM_EVENTS_PER_THETA; j++ ){ traj_colour[j] = (Int_t)kRed; }

		// Initialise phi solid angle array parameters
		d_phi[i][0] = 0.0;
		d_phi[i][1] = 0.0;
		d_phi[i][2] = 0.0;

		// Initialise counters for how many successful hits there are
		num_proton_hits[i] = 0;
		theta_frac_proton_recoil[i] = 0.0;
	



		// LOOP OVER NUMBER OF EVENTS FOR A GIVEN THETA (j) ------------------------------------- //
		for ( Int_t j = 0; j < NUM_EVENTS_PER_THETA; j++ ){
		// Calculate trajectories for a given theta
			STARTS_IN_VALID_SPOT = 1;
			phi = rand_phi.Uniform(0,360);			// This is initial phi angle
			beam_offset_r = rand_beam_offset_r.Gaus(0, BEAM_FWHM*FWHMToSigma() );
			beam_offset_theta = rand_beam_offset_theta.Uniform(0,180);
			phi_stop = 0;
			



			// LOOP OVER Z (PROTONS) (k) ------------------------------------------------------- //
			for ( Int_t k = 0; k < NUM_ZP; k++ ){
				// Calculate x,y,z
				if ( phi_stop == 0 ){
					z3[k] = -1.0*z_spacing*k;
					x3[k] = v3_perp*UC_CToMS()*UC_MToCM()*( TMath::Cos( phi*TMath::DegToRad() ) - TMath::Cos( phi*TMath::DegToRad() + cyclotron_freq*z3[k]/( v3_para*UC_CToMS()*UC_MToCM() ) ) )/cyclotron_freq + beam_offset_r*TMath::Cos( beam_offset_theta*TMath::DegToRad() ) + BEAM_SPOT_OFF_X;
					y3[k] = v3_perp*UC_CToMS()*UC_MToCM()*( -TMath::Sin( phi*TMath::DegToRad() ) + TMath::Sin( phi*TMath::DegToRad() + cyclotron_freq*z3[k]/( v3_para*UC_CToMS()*UC_MToCM() ) ) )/cyclotron_freq + beam_offset_r*TMath::Sin( beam_offset_theta*TMath::DegToRad() ) + BEAM_SPOT_OFF_Y;
					// Check if trajectory passed through the array
					if ( k == 0 ){
						if ( TMath::Abs( x3[k] ) > ARR_IN_DIAM/2 || TMath::Abs( y3[k] ) > ARR_IN_DIAM/2 ){
							// Check that beam offset is not too high
							if ( TMath::Abs( BEAM_SPOT_OFF_X ) > ARR_IN_DIAM/2 || TMath::Abs( BEAM_SPOT_OFF_Y ) > ARR_IN_DIAM/2 ){
								std::cerr << "Centre of beam spot does not pass through array ever! EXIT\n";
								std::exit(1);
							}

							// Outside defined region. Try again by decrementing j and doing loop again
							STARTS_IN_VALID_SPOT = 0;
						}
						else{
							h_phi_dist->Fill(phi);
							h_beam_width->Fill( x3[k],y3[k] );
						}
					}
				}
				else if ( phi_stop == 1 && k > 0 ){
					x3[k] = x3[k-1];
					y3[k] = y3[k-1];
					z3[k] = z3[k-1];
				}

				// Calculate success/failures
				if ( phi_stop == 0 && STARTS_IN_VALID_SPOT == 1 ){
					x_success = ( TMath::Abs( x3[k] ) <= ARR_DIAM/2 && TMath::Abs( y3[k] ) < SI_HEIGHT/2 && z3[k] <= z_fj );
					y_success = ( TMath::Abs( y3[k] ) <= ARR_DIAM/2 && TMath::Abs( x3[k] ) < SI_HEIGHT/2 && z3[k] <= z_fj );
					x_failure = ( TMath::Abs( x3[k] ) <= ARR_DIAM/2 && TMath::Abs( y3[k] ) > SI_HEIGHT/2 && TMath::Abs( y3[k] ) <= ARR_DIAM/2 && z3[k] <= z_fj );
					y_failure = ( TMath::Abs( y3[k] ) <= ARR_DIAM/2 && TMath::Abs( x3[k] ) > SI_HEIGHT/2 && TMath::Abs( x3[k] ) <= ARR_DIAM/2 && z3[k] <= z_fj );
					z_failure = ( TMath::Abs( x3[k] ) <= ARR_DIAM/2 && TMath::Abs( y3[k] ) <= ARR_DIAM/2 && z3[k] <= z_fj + 0.5*z_spacing && z3[k] >= z_fj - 0.5*z_spacing );
				

					if ( ( x_success || y_success ) ){
						phi_stop = 1;
						traj_colour[j] = (Int_t)kGreen;
						d_phi[i][0] += phi_cont;
						num_proton_hits[i]++;
					}
					else if ( ( x_failure || y_failure ) ){
						phi_stop = 1;
						traj_colour[j] = (Int_t)kOrange;
						d_phi[i][1] += phi_cont;
					}
					else if ( z_failure ){
						phi_stop = 1;
						traj_colour[j] = (Int_t)kRed;
						d_phi[i][2] += phi_cont;
					}

				}

			}	// Loop over z (protons) (k)




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
			}	// Loop over z (recoils) (l)


			// CHECK TRAJECTORIES FOR VALIDITY AND PLOT GRAPHS --------------------------------- //
			if ( STARTS_IN_VALID_SPOT == 1 ){
				// Draw successful trajectories in head-on plots
				if ( theta_cm[i] == THETA_HEAD ){
					// protons
					g_traj_ejectile_head[j] = new TGraph( NUM_ZP, x3, y3 );
					g_traj_ejectile_head[j]->SetLineColor( traj_colour[j] );

					beam_start[j] = new TEllipse( x3[0], y3[0], 0.01);
					beam_start[j]->SetFillColor( traj_colour[j] );

					// recoils
					g_traj_recoil_head[j] = new TGraph( NUM_ZR, x4, y4 );
					g_traj_recoil_head[j]->SetLineColor( kBlack );

					
				}
				
				// Test how many recoils worked
				Double_t rtemp = TMath::Sqrt( y4[NUM_ZR-1]*y4[NUM_ZR-1] + x4[NUM_ZR-1]*x4[NUM_ZR-1] );
				if ( rtemp >= RDT_RADIUS_TO_CLEAR  && rtemp <= RDT_SI_OUTER_RAD + 0.5*RDT_DETECTOR_GAP ){
					theta_frac_proton_recoil[i] += 100.0/NUM_EVENTS_PER_THETA;
				}

				// Draw trajectories in theta for a random phi
				if ( j == EVENT_NUMBER && i % TMath::Max( (Int_t)TMath::Floor( (Double_t)NUM_THETA/50.0 ), 1 ) == 0 ){
					// Ejectiles
					g_traj_ejectile_side[i] = new TGraph( NUM_ZP, z3, y3 );
					g_traj_ejectile_side[i]->SetLineColor( traj_colour[j] );

					// recoils
					g_traj_recoil_side[i] = new TGraph( NUM_ZR, z4, y4 );
					if ( rtemp >= RDT_SI_INNER_RAD + 0.5*RDT_DETECTOR_GAP && rtemp <= RDT_SI_OUTER_RAD + 0.5*RDT_DETECTOR_GAP ){
						g_traj_recoil_side[i]->SetLineColor(kGreen);
					}
					else if ( rtemp >= RDT_PCB_INNER_RAD + 0.5*RDT_DETECTOR_GAP && rtemp < RDT_SI_INNER_RAD + 0.5*RDT_DETECTOR_GAP ){
						g_traj_recoil_side[i]->SetLineColor(kOrange);
					}
					else{
						g_traj_recoil_side[i]->SetLineColor(kRed);
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
		if ( i == 0 ){
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




	/*********************************************************************************************/
	/*                                    NOW PLOT EVERYTHING                                    */
	/*********************************************************************************************/
	// Open TFile in which to save the canvases etc.
	TFile* f = new TFile("output-data/canvases.root", "RECREATE");


	// PROTONS SIDE ON ------------------------------------------------------------------------- //
	// Set style
	CreateStyle( ptm_style );
	ptm_style->cd();
	TCanvas* c_traj_ejectile_side = new TCanvas( "c_traj_ejectile_side", "EJECTILE SIDE", CANVAS_WIDTH, CANVAS_HEIGHT );
	
	// Draw frame
	TH1F* frame = c_traj_ejectile_side->DrawFrame( ej_side_X1, ej_side_Y1, ej_side_X2, ej_side_Y2 );
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

	// Draw trajectories
	Int_t ctr_traj_ejectile_side = 0;
	for ( Int_t i = 0; i < NUM_THETA; i++ ){
		if ( g_traj_ejectile_side[i] != NULL && ctr_traj_ejectile_side <= NUM_TRAJECTORIES ){
			g_traj_ejectile_side[i]->Draw("SAME");
			ctr_traj_ejectile_side++;
		}
	}
	

	// PROTONS HEAD ON ------------------------------------------------------------------------- //
	TCanvas* c_traj_ejectile_head = new TCanvas( "c_traj_ejectile_head", "EJECTILE HEAD", CANVAS_HEIGHT, CANVAS_HEIGHT );
	c_traj_ejectile_head->cd();
	ptm_style->cd();
	TH1F* frame1 = c_traj_ejectile_head->DrawFrame( -10, -10, 10, 10 );
	frame1->SetTitle("");
	frame1->GetXaxis()->SetTitle("x_{p} / cm");
	frame1->GetYaxis()->SetTitle("y / cm");

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
	for ( Int_t j = 0; j < NUM_EVENTS_PER_THETA; j++ ){
		if ( g_traj_ejectile_head[j] != NULL && ctr_traj_ejectile_head <= NUM_TRAJECTORIES ){
			g_traj_ejectile_head[j]->Draw("C SAME");
			beam_start[j]->Draw("SAME");
			ctr_traj_ejectile_head++;
		}
	}

	
	// DRAW THE PHI FRACTION ------------------------------------------------------------------- //
	TCanvas* c_phi_frac = new TCanvas( "c_phi_frac", "PHI FRACTION", CANVAS_WIDTH, CANVAS_HEIGHT );
	c_phi_frac->cd();	
	ptm_style->cd();

	TH1F* frame2 = c_phi_frac->DrawFrame( TMath::Max( 0.0, theta_lb - theta_spacing ), 0, theta_ub + theta_spacing, 6.5 );
	frame2->SetTitle("");
	frame2->GetXaxis()->SetTitle( "#theta_{cm} / ^{#circ}" );
	frame2->GetYaxis()->SetTitle( "#Delta#phi / rad" );

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


	// DRAW THE RANDOM NUMBER DISTRIBUTIONS ---------------------------------------------------- //
	TCanvas* c_rand_numbers = new TCanvas( "c4", "RANDOM NUMBERS", CANVAS_WIDTH, 0.6465*CANVAS_HEIGHT);
	c_rand_numbers->Divide(2,1);

	c_rand_numbers->cd(1);
	ptm_style->cd();
	TPad* pad = (TPad*)c_rand_numbers->GetPad(1);
	pad->SetLeftMargin(0.13);
	pad->SetRightMargin(0.01);
	pad->SetBottomMargin(0.09);
	pad->SetTopMargin(0.02);
	h_phi_dist->SetTitle("; #phi (#circ);Counts" );
	h_phi_dist->GetYaxis()->SetTitleOffset(2.0);
	h_phi_dist->GetYaxis()->SetRangeUser( 0, h_phi_dist->GetMaximum() + 2 );
	h_phi_dist->Draw();

	c_rand_numbers->cd(2);
	ptm_style->cd();
	pad = (TPad*)c_rand_numbers->GetPad(2);
	pad->SetLeftMargin(0.09);
	pad->SetRightMargin(0.13);
	pad->SetBottomMargin(0.09);
	pad->SetTopMargin(0.02);
	
	h_beam_width->SetTitle(";x / cm;y / cm");
	h_beam_width->Draw("colz");
	h_beam_width->GetYaxis()->SetTitleOffset(1.3);

	TLine* border[4];
	border[0] = new TLine( -ARR_IN_DIAM/2, -ARR_IN_DIAM/2, -ARR_IN_DIAM/2, ARR_IN_DIAM/2 );
	border[1] = new TLine( ARR_IN_DIAM/2, -ARR_IN_DIAM/2, ARR_IN_DIAM/2, ARR_IN_DIAM/2 );
	border[2] = new TLine( -ARR_IN_DIAM/2, -ARR_IN_DIAM/2, ARR_IN_DIAM/2, -ARR_IN_DIAM/2 );
	border[3] = new TLine( -ARR_IN_DIAM/2, ARR_IN_DIAM/2, ARR_IN_DIAM/2, ARR_IN_DIAM/2 );

	for ( Int_t i = 0; i < 4; i++ ){
		border[i]->SetLineWidth(1);
		border[i]->SetLineColor(kBlack);
		border[i]->Draw("SAME");
	}


	// RECOILS SIDE ON ------------------------------------------------------------------------- //
	// Set style
	CreateStyle( ptm_style );
	ptm_style->cd();
	TCanvas* c_traj_recoil_side = new TCanvas( "c_traj_recoil_side", "RECOIL SIDE", CANVAS_WIDTH, CANVAS_HEIGHT );
	
	// Draw frame
	TH1F* frame_rec_side = c_traj_recoil_side->DrawFrame( rec_side_X1, rec_side_Y1, rec_side_X2, rec_side_Y2 );
	frame_rec_side->SetTitle( ";z / cm; r / cm" );

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
	for ( Int_t i = 0; i < NUM_THETA; i++ ){
		if ( g_traj_recoil_side[i] != NULL && ctr_traj_recoil_side <= NUM_TRAJECTORIES ){
			g_traj_recoil_side[i]->Draw("SAME");
			ctr_traj_recoil_side++;
		}
	}


	// RECOILS HEAD ON ------------------------------------------------------------------------- //
	TCanvas* c_traj_recoil_head = new TCanvas( "c_traj_recoil_head", "RECOIL HEAD", CANVAS_HEIGHT, CANVAS_HEIGHT );
	TPad* pad_traj_recoil_head = (TPad*)c_traj_recoil_head;
	c_traj_recoil_head->cd();
	ptm_style->cd();
	TH1F* frame_rec_head = c_traj_recoil_head->DrawFrame( -( 1.2*RDT_PCB_OUTER_RAD ), -( 1.2*RDT_PCB_OUTER_RAD ), ( 1.2*RDT_PCB_OUTER_RAD ), ( 1.2*RDT_PCB_OUTER_RAD ) );
	frame_rec_head->SetTitle("");
	frame_rec_head->GetXaxis()->SetTitle("x_{rec} / cm");
	frame_rec_head->GetYaxis()->SetTitle("y / cm");

	// Draw the recoil detector
	TEllipse* rdt_pcb[4];
	TEllipse* rdt_Si[4];
	TEllipse* rdt_Si_inner_pcb[4];
	TEllipse* rdt_hole;

	for ( Int_t i = 0; i < 4; i++ ){
		// Draw large PCBs
		rdt_pcb[i] = new TEllipse( 0.5*RDT_DETECTOR_GAP*TMath::Power( -1.0, TMath::Floor( (i+1)/2 ) ), 0.5*RDT_DETECTOR_GAP*( i > 1 ? -1 : 1 ), RDT_PCB_OUTER_RAD, RDT_PCB_OUTER_RAD, i*90 + RDT_ROTATION, (i+1)*90 + RDT_ROTATION );
		rdt_pcb[i]->SetFillColor( pcb_green_i );
		rdt_pcb[i]->SetLineWidth(0);
		rdt_pcb[i]->Draw("SAME");

		// Draw Si detectors
		rdt_Si[i] = new TEllipse( 0.5*RDT_DETECTOR_GAP*TMath::Power( -1.0, TMath::Floor( (i+1)/2 ) ), 0.5*RDT_DETECTOR_GAP*( i > 1 ? -1 : 1 ), RDT_SI_OUTER_RAD, RDT_SI_OUTER_RAD, i*90 + 0.8*( 90 - RDT_ANGULAR_COVERAGE ) + RDT_ROTATION, (i+1)*90 - 0.2*( 90 - RDT_ANGULAR_COVERAGE ) + RDT_ROTATION );
		rdt_Si[i]->SetFillColor( si_strip_i );
		rdt_Si[i]->SetLineWidth(0);
		rdt_Si[i]->Draw("SAME");

		// Draw inner PCBs
		rdt_Si_inner_pcb[i] = new TEllipse( 0.5*RDT_DETECTOR_GAP*TMath::Power( -1.0, TMath::Floor( (i+1)/2 ) ), 0.5*RDT_DETECTOR_GAP*( i > 1 ? -1 : 1 ), RDT_SI_INNER_RAD, RDT_SI_INNER_RAD, i*90 + RDT_ROTATION, (i+1)*90 + RDT_ROTATION );
		rdt_Si_inner_pcb[i]->SetFillColor( pcb_green_i );
		rdt_Si_inner_pcb[i]->SetLineWidth(0);
		rdt_Si_inner_pcb[i]->Draw("SAME");
	}
	rdt_hole = new TEllipse( 0, 0, RDT_PCB_INNER_RAD, RDT_PCB_INNER_RAD, 0, 360 );
	rdt_hole->SetFillColor(kWhite);
	rdt_hole->SetLineWidth(0);
	rdt_hole->Draw("SAME");
	
	// Draw trajectories
	Int_t ctr_traj_recoil_head = 0;
	for ( Int_t j = 0; j < NUM_EVENTS_PER_THETA; j++ ){
		// Set the colour based on if the protons hit the array
		if ( g_traj_ejectile_head[j] != NULL && ctr_traj_recoil_head <= NUM_TRAJECTORIES ){
			if ( g_traj_ejectile_head[j]->GetLineColor() == kGreen ){
				g_traj_recoil_head[j]->SetLineColor(kGreen);
				
				// Calculate whether the trajectories hit the detector
				Double_t xtemp = g_traj_recoil_head[j]->GetX()[NUM_ZR-1];
				Double_t ytemp = g_traj_recoil_head[j]->GetY()[NUM_ZR-1];
				Double_t rtemp = TMath::Sqrt( xtemp*xtemp + ytemp*ytemp );

				if ( rtemp >= RDT_RADIUS_TO_CLEAR ){
					// Successful hit
					g_traj_recoil_head[j]->SetLineColor(kGreen);
				}
				else if ( rtemp < RDT_RADIUS_TO_CLEAR && rtemp >= RDT_SI_INNER_RAD + 0.5*RDT_DETECTOR_GAP ){
					// Hits PCB but not the Si
					g_traj_recoil_head[j]->SetLineColor(kOrange);
				}
				else if ( rtemp < RDT_SI_INNER_RAD + 0.5*RDT_DETECTOR_GAP ){
					// Misses both PCB and Si
					g_traj_recoil_head[j]->SetLineColor(kRed);
				}
				else{
					// Does something else
					g_traj_recoil_head[j]->SetLineColor(kBlack);
					std::cout << xtemp << "\t" << ytemp << "\t" << rtemp << "\t" << "\n";
				}
				
				g_traj_recoil_head[j]->Draw("C SAME");
				ctr_traj_recoil_head++;
			}
		}
	}
	

	// FRACTION OF THETA ACCEPTED ON RECOIL DETECTORS ------------------------------------------ //
	TCanvas* c_recoil_theta_frac = new TCanvas( "c_recoil_theta_frac", "THETA FRACTION ON RDT", CANVAS_WIDTH, CANVAS_HEIGHT );
	c_recoil_theta_frac->cd();	
	ptm_style->cd();

	TH1F* fr_recoil_theta_frac = c_recoil_theta_frac->DrawFrame( TMath::Max( 0.0, theta_lb - theta_spacing ), 0, theta_ub + theta_spacing, 105.0 );
	fr_recoil_theta_frac->SetTitle("");
	fr_recoil_theta_frac->GetXaxis()->SetTitle( "#theta_{cm} / ^{#circ}" );
	fr_recoil_theta_frac->GetYaxis()->SetTitle( "% hits on RDT" );

	g_recoil_theta_frac = new TGraph( NUM_THETA, theta_cm, theta_frac_proton_recoil );
	g_recoil_theta_frac->SetMarkerStyle(20);
	g_recoil_theta_frac->SetMarkerSize(0.2);
	g_recoil_theta_frac->GetXaxis()->CenterTitle();
	g_recoil_theta_frac->GetYaxis()->CenterTitle();
	g_recoil_theta_frac->SetMarkerColor( kRed );
	g_recoil_theta_frac->Draw("P");
	

	// PRINT THE CANVASES AND FREE MEMORY------------------------------------------------------- //
	c_traj_ejectile_side->Print( Form( "output-data/EX_%s-POS_%i-traj_ejectile_side.svg", DecimalDotToUnderscore( ex ).Data(), POSITION ) );
	c_traj_ejectile_head->Print( Form( "output-data/EX_%s-POS_%i-traj_ejectile_head.svg", DecimalDotToUnderscore( ex ).Data(), POSITION ) );
	c_phi_frac->Print( Form( "output-data/EX_%s-POS_%i-phi_frac.svg", DecimalDotToUnderscore( ex ).Data(), POSITION ) );
	c_rand_numbers->Print( Form( "output-data/EX_%s-POS_%i-rand_numbers.svg", DecimalDotToUnderscore( ex ).Data(), POSITION ) );
	c_traj_recoil_side->Print( Form( "output-data/EX_%s-POS_%i-traj_recoil_side.svg", DecimalDotToUnderscore( ex ).Data(), POSITION ) );
	c_traj_recoil_head->Print( Form( "output-data/EX_%s-POS_%i-traj_recoil_head.svg", DecimalDotToUnderscore( ex ).Data(), POSITION ) );
	c_recoil_theta_frac->Print( Form( "output-data/EX_%s-POS_%i-recoil_theta_frac.svg", DecimalDotToUnderscore( ex ).Data(), POSITION ) );

	// Write the canvases to the file
	c_traj_ejectile_side->Write();
	c_traj_ejectile_head->Write();
	c_phi_frac->Write();
	c_rand_numbers->Write();
	c_traj_recoil_side->Write();
	c_traj_recoil_head->Write();
	c_recoil_theta_frac->Write();


	
	// Close the files
	log_file.close();
	f->Close();


	if ( !BATCH_MODE ){ gPad->WaitPrimitive("TPave"); }
	std::cout << "Deleting variables...\n";

	// Clear the memory
	// Arrays
	delete[] theta_cm;
	delete[] theta_lab;
	delete[] x3;
	delete[] y3;
	delete[] z3;
	delete[] x4;
	delete[] y4;
	delete[] z4;
	delete[] num_proton_hits;
	delete[] theta_frac_proton_recoil;
	delete[] d_phi;

	// TCanvases
	if ( c_traj_ejectile_side->IsOnHeap() ){ delete c_traj_ejectile_side; }
	if ( c_traj_ejectile_head->IsOnHeap() ){ delete c_traj_ejectile_head; }
	if ( c_phi_frac->IsOnHeap() ){ delete c_phi_frac; }
	if ( c_rand_numbers->IsOnHeap() ){ delete c_rand_numbers; }
	if ( c_traj_recoil_side->IsOnHeap() ){ delete c_traj_recoil_side; }
	if ( c_traj_recoil_head->IsOnHeap() ){ delete c_traj_recoil_head; }
	if ( c_recoil_theta_frac->IsOnHeap() ){ delete c_recoil_theta_frac; }

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

	// TGraphs
	for ( Int_t i = 0; i < NUM_EVENTS_PER_THETA; i++ ){
		if ( g_traj_ejectile_head[i] != NULL ){ if ( g_traj_ejectile_head[i]->IsOnHeap() ){ g_traj_ejectile_head[i]->Delete(); } }
		if ( g_traj_recoil_head[i] != NULL ){ if ( g_traj_recoil_head[i]->IsOnHeap() ){ g_traj_recoil_head[i]->Delete(); } }
	}
	for ( Int_t i = 0; i < NUM_THETA; i++ ){
		if ( g_traj_ejectile_side[i] != NULL ){ if ( g_traj_ejectile_side[i]->IsOnHeap() ){ g_traj_ejectile_side[i]->Delete(); } }
		if ( g_traj_recoil_side[i] != NULL ){ if ( g_traj_recoil_side[i]->IsOnHeap() ){ g_traj_recoil_side[i]->Delete(); } }
	}
	if ( g_phi_frac->IsOnHeap() ){ g_phi_frac->Delete(); }
	if ( g_recoil_theta_frac->IsOnHeap() ){ g_recoil_theta_frac->Delete(); }
	

	// TH's
	
	if ( h_beam_width->IsOnHeap() ){ h_beam_width->Delete(); }
	if ( h_phi_dist->IsOnHeap() ){ h_phi_dist->Delete(); }
	/*if ( frame->IsOnHeap() ){ frame->Delete(); }
	if ( frame1->IsOnHeap() ){ frame1->Delete(); }
	if ( frame2->IsOnHeap() ){ frame2->Delete(); }
	if ( frame_rec_side->IsOnHeap() ){ frame_rec_side->Delete(); }
	if ( frame_rec_head->IsOnHeap() ){ frame_rec_head->Delete(); }
	if ( fr_recoil_theta_frac->IsOnHeap() ){ fr_recoil_theta_frac->Delete(); }*/


	// TLines
	if ( l_solid_angle->IsOnHeap() ){ l_solid_angle->Delete(); }
	for ( Int_t i = 0; i < 4; i++ ){
		if ( border[i]->IsOnHeap() ){ border[i]->Delete(); }
	}

	// TEllipses
	for ( Int_t i = 0; i < NUM_EVENTS_PER_THETA; i++ ){
		if ( beam_start[i] != NULL ){ if ( beam_start[i]->IsOnHeap() ){ beam_start[i]->Delete(); } }
	}


	// TPads
	//if ( pad->IsOnHeap() ){ pad->Delete(); }
}



// MAIN FUNCTION ------------------------------------------------------------------------------- //
void ArrayGeometry(){
	const Int_t NUM_STATES = 11;
	const Bool_t ALL = 0;
	Double_t STATES[NUM_STATES] = {
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
	};

	if ( ALL == 1 ){
		for ( Int_t i = 0; i < NUM_STATES; i++ ){
			ArrayGeometryEX( STATES[i] );
			gROOT->ls();
		}
	}
	else{
		ArrayGeometryEX( STATES[NUM_STATES-1] );
	}

}



/* *** TODO ***
 * Residuals for phi_frac plot
 * Error bars for phi_frac plot?
 * Save in ROOT files (in tree?) for ease of manipulation
 * Recoil solid angle stuff
*/



















