// AG_constants.C
// Constants for the file ArrayGeometry.C
// ============================================================================================= //
// Patrick MacGregor
// Nuclear Physics Research Group
// Department of Physics and Astronomy
// The University of Manchester
// ============================================================================================= //
#ifndef AG_CONSTANTS_H_
#define AG_CONSTANTS_H_

#include "AG_functions.h"
#include "AG_style.h"

#include <TColor.h>

// PLOTTING CONSTANTS -------------------------------------------------------------------------- //
const Int_t CANVAS_WIDTH = 1200;
const Int_t CANVAS_HEIGHT = 900;
const Bool_t BATCH_MODE = 1;

// Set colours
Int_t pcb_green_i = TColor::GetFreeColorIndex();
TColor *pcb_green = new TColor( pcb_green_i, 0.0, 111.0/255.0, 75.0/255.0 );

Int_t si_strip_i = TColor::GetFreeColorIndex();
TColor *si_strip = new TColor( si_strip_i, 141.0/255.0, 116.0/255.0, 71.0/255.0 );

Int_t fj_red_i = TColor::GetFreeColorIndex();
TColor *fj_red = new TColor( fj_red_i, 189.0/255.0, 124.0/255.0, 124.0/255.0 );

// EXPERIMENT-SPECIFIC CONSTANTS --------------------------------------------------------------- //
/* Particle labels:
	* 1 = 28Mg
	* 2 = d
	* 3 = p
	* 4 = 29Mg
*/

// Position
const Int_t POSITION = 1;

// Mass Excesses [MeV] (obtained from https://www.nndc.bnl.gov/wallet/wc8.html)
const Double_t Delta[4] = { -15.018, 13.136, 7.289, -10.60 };
const Int_t A[4] = { 28, 2, 1, 29 };

// Masses [MeV/c^2]
const Double_t mass[4] = {
	MassExcessToMass( Delta[0], A[0] ),
	MassExcessToMass( Delta[1], A[1] ),
	MassExcessToMass( Delta[2], A[2] ),
	MassExcessToMass( Delta[3], A[3] )
};

// Default excitation [MeV]
const Double_t EX = 0.0000;						// [MeV]

// Constants of reaction
Double_t ejectile_charge = 1;					// [e]
Double_t recoil_charge = 12;
Double_t b_field = 2.5;							// [T]
Double_t beam_energy_per_nucleon = 9.473;		// [MeV per u]
const Int_t POS = 1;

// Non-angle dependent quantities
Double_t Q = ( ( mass[0] + mass[1] ) - ( mass[2] + mass[3] ) );								// Nuclear Q-value [MeV]
Double_t T1 = beam_energy_per_nucleon*A[0];													// Beam energy [MeV]
Double_t T_cm_i = ( mass[1]/( mass[0] + mass[1] ) )*T1;										// Initial CM energy [MeV]
Double_t V_cm = TMath::Sqrt( 2*T1*mass[0] )/( mass[0] + mass[1] );							// CM velocity [1/c]
Double_t cyclotron_freq = ejectile_charge*UC_EToC()*b_field/( mass[2]/UC_KGToMEV() );		// [s^-1]
Double_t cyclotron_freq_recoil = recoil_charge*UC_EToC()*b_field/( mass[3]/UC_KGToMEV() );	// [s^-1]

// GEOMETRY ------------------------------------------------------------------------------------ //
// Geometry of the detector
// Position is first index, det # is second
Double_t Si_centroids[2][6] = {
	{ 45.366, 39.485, 33.609, 27.746, 21.910, 16.174 },
	{ 42.368, 36.487, 30.611, 24.748, 18.912, 13.176 }
};

const Double_t SI_WIDTH = 5.5;			// [cm]
const Double_t SI_HEIGHT = 0.9;			// [cm]

const Double_t PCB_WIDTH = 2.3;			// [cm]
const Double_t PCB_LENGTH = 38.8;		// [cm]

const Double_t ARR_DIAM = 2.3;			// [cm]
const Double_t ARR_IN_DIAM = 1.0;		// [cm]

const Double_t SLIT_LENGTH = 3.271;		// [cm]

const Double_t LAST_CENTROID_TO_PCB = 3.405;	// [cm]

const Double_t TARGET_RDT_DISTANCE = GetArrayRDTDistance(POS);	// [cm] --> Edith's report

// Calculate z at the four-jaws
const Double_t z_fj = - Si_centroids[POSITION][5] + LAST_CENTROID_TO_PCB + SLIT_LENGTH;

// Ejectile side-on plotting
Double_t ej_side_X1 = -50.0;														// Left hand edge
Double_t ej_side_X2 = 5.0;															// Right hand edge
Double_t ej_side_DY = (ej_side_X2-ej_side_X1)*(1 - marg_t - marg_b)*h/( (1 - marg_r - marg_l)*w );	// The resulting height of the plot
Double_t ej_side_Y1 = -0.5*ej_side_DY;														// Symmetrical spacing either side
Double_t ej_side_Y2 = 0.5*ej_side_DY;

// Ejectile side-on plotting
Double_t rec_side_X1 = -10.0;																				// Left hand edge
Double_t rec_side_X2 = TARGET_RDT_DISTANCE + 10.0;														// Right hand edge
Double_t rec_side_Y1 = -10.0;																// Symmetrical spacing either side
Double_t rec_side_Y2 = 10.0;

// Beam properties
const Double_t BEAM_FWHM = 0.2;					// [cm]
const Double_t BEAM_SPOT_OFF_X = 0.0;			// [cm] N.B. x is defined looking from the target towards the array
const Double_t BEAM_SPOT_OFF_Y = 0.0;			// [cm]

// Recoil detector properties (Micron Semiconductor QQQ1)
const Double_t RDT_SI_INNER_RAD = 0.9;			// [cm]		* Inner radius of silicon pad
const Double_t RDT_SI_OUTER_RAD = 5.0;			// [cm]		* Outer radius of silicon pad
const Double_t RDT_ANGULAR_COVERAGE = 82.0;		// [DEG]	* Angular coverage of silicon pad
const Double_t RDT_PCB_INNER_RAD = 0.525;		// [cm]		* Inner radius of PCB board
const Double_t RDT_PCB_OUTER_RAD = 5.15;		// [cm]		* Outer radius of PCB board
const Double_t RDT_ROTATION = 0.0;				// [DEG]	* Rotate the detector by this angle
const Double_t RDT_DETECTOR_GAP = 0.1;			// [cm]		* The gap between each detector from the others
const Double_t RDT_RADIAL_ALIGNMENT = 0.1;		// [cm]		* The radial tolerance that each trajectory must clear
const Double_t RDT_RADIUS_TO_CLEAR = RDT_SI_INNER_RAD + TMath::Sqrt(2)*0.5*RDT_DETECTOR_GAP + TMath::Sqrt(2)*0.5*RDT_RADIAL_ALIGNMENT;

// SAMPLING ------------------------------------------------------------------------------------ //
// CM Theta's to spam
Double_t theta_spacing = 0.01;	// [DEG]
Double_t THETA_LB= 15.00;		// [DEG]
Double_t THETA_UB = 19.00;		// [DEG]

// Number of events per theta
const Int_t NUM_EVENTS_PER_THETA = 10000;

// Proton z's to spam
Double_t z_spacing = 0.01;		// [cm]
const Int_t NUM_ZP = (Int_t)( ( 0 - ej_side_X1 )/z_spacing );
const Int_t NUM_THETA = (Int_t)( ( THETA_UB - THETA_LB )/theta_spacing ) + 1;

// Recoil z's to spam (add 1 so that the detector surface is definitely included)
const Int_t NUM_ZR = (Int_t)( ( TARGET_RDT_DISTANCE - 0 )/z_spacing ) + 1;

// Choose which angles to look at
const Double_t THETA_HEAD = 18.0;
const Double_t EVENT_NUMBER = 1;

// Log file precision
const Int_t log_prec = 8;
const Int_t log_width = 16;

// Number of trajectories to plot
const Int_t NUM_TRAJECTORIES = 200;




























#endif
