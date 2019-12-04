// AG_constants.C
// Constants for the file ArrayGeometry.C
// ============================================================================================= //
// Patrick MacGregor
// Nuclear Physics Research Group
// School of Physics and Astronomy
// The University of Manchester
// ============================================================================================= //
#ifndef AG_CONSTANTS_H_
#define AG_CONSTANTS_H_

#include "AG_physical_constants.h"
#include "AG_functions.h"
#include "AG_style.h"

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

// States [MeV]
const Double_t EX = 0.0000;						// [MeV]

Double_t ejectile_charge = 1;					// [e]
Double_t b_field = 2.5;							// [T]
Double_t beam_energy_per_nucleon = 9.473;		// [MeV per u]

// Non-angle dependent quantities
Double_t Q = ( ( mass[0] + mass[1] ) - ( mass[2] + mass[3] ) );							// Nuclear Q-value [MeV]
Double_t T1 = beam_energy_per_nucleon*A[0];												// Beam energy [MeV]
Double_t T_cm_i = ( mass[1]/( mass[0] + mass[1] ) )*T1;									// Initial CM energy [MeV]
Double_t V_cm = TMath::Sqrt( 2*T1*mass[0] )/( mass[0] + mass[1] );						// CM velocity [1/c]
Double_t cyclotron_freq = ejectile_charge*UC_EToC()*b_field/( mass[2]/UC_KGToMEV() );	// [s^-1]

// GEOMETRY ------------------------------------------------------------------------------------ //
// Theta's to spam
Double_t theta_spacing = 0.1;	// [DEG]
Double_t theta_lb = 11.0;		// [DEG]
Double_t theta_ub = 25.0;		// [DEG]

// Number of events per theta
const Int_t NUM_EVENTS_PER_THETA = 1000;

// z's to spam
Double_t z_spacing = 0.02;		// [cm]
const Int_t NUM_Z = (Int_t)( ( 0 - X1 )/z_spacing );
const Int_t NUM_THETA = (Int_t)( ( theta_ub - theta_lb )/theta_spacing ) + 1;

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


// Beam properties
const Double_t BEAM_FWHM = 0.2;					// [cm]
const Double_t BEAM_SPOT_OFF_X = ARR_IN_DIAM/2;	// [cm]
const Double_t BEAM_SPOT_OFF_Y = ARR_IN_DIAM/2;	// [cm]























#endif
