// CMAngleCalculator.h
// Contains useful functions for the CMAngleCalculator file
// ============================================================================================= //
// Patrick MacGregor
// Nuclear Physics Research Group
// School of Physics and Astronomy
// The University of Manchester
// ============================================================================================= //

#ifndef CM_ANGLE_CALCULATOR_H_
#define CM_ANGLE_CALCULATOR_H_

#include <TString.h>


// CONSTANTS ---===---===---===---===---===---===---===---===---===---===---===---===---===---== //
Double_t c = 299792458;					// Speed of light in m /s
Double_t u = 931.494;					// Unified atomic mass unit in MeV / c^2

// SWITCH
// 0 = d(28Mg,p)29Mg
// 1 = d(28Si,p)29Si
Int_t reaction = 0;
const Int_t NUM_REACTIONS = 2;

TString reaction_string_arr[2] = {
	"d(28Mg,p)29Mg",
	"d(28Si,p)29Si"
};


// GENERAL FUNCTIONS ---===---===---===---===---===---===---===---===---===---===---===---===--- //
// Print the current reaction
void PrintReaction(){
	for ( Int_t i = 0; i < NUM_REACTIONS; i++ ){
		std::cout << ( reaction == i ? "> " : "  " ) << Form( "%02d ", i ) << reaction_string_arr[i] << "\n";
	}
	return;
}

// Choose the reaction to use
void ChooseReaction(Int_t a){
	if ( a >= 0 && a < NUM_REACTIONS ){
		std::cout << "Reaction changed from " << reaction_string_arr[reaction] << " to " << reaction_string_arr[a] << "\n";
		reaction = a;
	}
	else{
		std::cout << "There are only " << NUM_REACTIONS << " available. Please choose a number between 0 and " << NUM_REACTIONS - 1 << ".\n";
		PrintReaction();
	}
	return;
}

// Calculates the average radius of the ISS array (be consistent with units!!)
/*/ The X-axis has its origin in the centre of the array (not necessarily the centre of the strip)
 *  X1 is the point furthest left of the origin
 *  X2 is the point furthest right of the origin
 *  height is the height of the array                                                            */
Double_t ISSArrayRadius( Double_t X1, Double_t X2, Double_t height ){
	Double_t A1 = X1/height;
	Double_t A2 = X2/height;
	return ( height*height/( 2*( X2 - X1 ) ) )*( A2*TMath::Sqrt( 1 + A2*A2 ) - A1*TMath::Sqrt( 1 + A1*A1 ) + TMath::Log( TMath::Abs( A2 + TMath::Sqrt( 1 + A2*A2 ) ) ) - TMath::Log( TMath::Abs( A1 + TMath::Sqrt( 1 + A1*A1 ) ) ) );
}

// Calculates nuclear mass from mass excess (and charge i.e. stripped electrons)
Double_t MassExcessToNuclearMass( Double_t Delta, Int_t A, Int_t q = 0 ){
	return A*u + Delta - q*0.511;
}

// Energies in MeV
Double_t BindEnToNuclearMass( Double_t bind_en_per_nucl, Int_t A, Int_t Z, Int_t q = 0 ){
	return Z*1.0078250322*u + (A-Z)*939.56542052 - A*bind_en_per_nucl - q*0.511;
}


// STRUCTS ---===---===---===---===---===---===---===---===---===---===---===---===---===---===- //
struct ReactionParameters{
	Double_t q;
	Double_t B;
	Double_t m1;
	Double_t m2;
	Double_t m3;
	Double_t m4;
	Double_t bepn;
};

struct KinematicsParameters{
	Double_t T1; 			// Kinetic energy of the beam [LAB] [ 9.473 MeV / u ]
	Double_t e1; 			// Total energy of beam particle [LAB] in MeV
	Double_t etot;			// Total energy [LAB] in MeV
	Double_t etot_cm ;		// Total energy [CM] (related to invariant mass) in MeV
	Double_t gamLAB_CM; 	// Gamma factor relating inertial frames (i.e. between LAB and CM)
	Double_t beta;			// Beta = v/c (v = velocity of CM frame)
	Double_t rho; 			// Radius of array in cm (approximate as circular - diameter of square is 23 mm).
	Double_t m4ex; 			// Invariant mass of the recoil nucleus
	Double_t e3_cm; 		// Energy of the outgoing proton [CM]
	Double_t p3_cm; 		// Momentum of the outgoing proton [CM]
};


// STRUCT FUNCTIONS ---===---===---===---===---===---===---===---===---===---===---===---===---= //
void Set28MgDP( ReactionParameters &rp ){
	rp.q = 1;												// Charge of outgoing proton in units of e
	rp.B = 2.5*c*(1e-8);									// Convert B-field to MeV / (cm*c*e)
	rp.m1 = BindEnToNuclearMass( 8.27241 , 28, 12,  12 ); 	// m1 = 28Mg in MeV / c^2		(beam)
	rp.m2 = BindEnToNuclearMass( 1.112283,  2,  1,  1 );	// m2 = d in MeV / c^2			(target)
	rp.m3 = 938.27208816;									// m3 = p in MeV / c^2			(ejectile)
	rp.m4 = BindEnToNuclearMass( 8.1132  , 29, 12, 12 );	// m4 = 29Mg in MeV / c^2		(recoil)
	rp.bepn = 9.473;										// bepn = beam energy per nucleon in MeV/u
	return;
}


void Set28SiDP( ReactionParameters &rp ){
	rp.q = 1;										// Charge of outgoing proton in units of e
	rp.B = 2.5*c*(1e-8);							// Convert B-field to MeV / (cm*c*e)
	rp.m1 = 0;//MassExcessToMass( -21.4927, 28,  9 ); 	// m1 = 28Mg in MeV / c^2		(beam)
	rp.m2 = 0;//MassExcessToMass(  13.1357,  2,  1 );	// m2 = d in MeV / c^2			(target)
	rp.m3 = 0;//MassExcessToMass(  7.28890,  1,  1 );	// m3 = p in MeV / c^2			(ejectile)
	rp.m4 = 0;//MassExcessToMass( -21.8950, 29, 14 );	// m4 = 29Mg in MeV / c^2		(recoil)
	rp.bepn = 9.473;								// bepn = beam energy per nucleon in MeV/u
	return;
}

// Fill structure with values
void FillRP( ReactionParameters &rp ){
	if ( reaction == 0 ){ Set28MgDP( rp ); }
	else if ( reaction == 1 ){ Set28SiDP( rp ); }
	else{ std::exit(1); }
	return;
}

// Calculate kinematic parameters
void CalculateKP( ReactionParameters &rp, KinematicsParameters &kp, Double_t ex ){
	kp.T1 = rp.bepn*rp.m1/u;																// Kinetic energy of the beam [LAB] [ 9.473 MeV / u ]
	kp.e1 = kp.T1 + rp.m1;																	// Total energy of beam particle [LAB] in MeV
	kp.etot = kp.e1 + rp.m2;																// Total energy [LAB] in MeV
	kp.etot_cm = TMath::Sqrt( rp.m1*rp.m1 + rp.m2*rp.m2 + 2*kp.e1*rp.m2 );					// Total energy [CM] (related to invariant mass) in MeV
	kp.gamLAB_CM = kp.etot/kp.etot_cm;														// Gamma factor relating inertial frames (i.e. between LAB and CM)
	kp.beta = TMath::Sqrt( 1 - TMath::Power( kp.gamLAB_CM, -2 ) );							// Beta = v/c (v = velocity of CM frame)
	kp.rho = ISSArrayRadius( -0.45, 0.45, 1.15 );											// Radius of array in cm (approximate as circular - diameter of square is 23 mm).
	kp.m4ex = rp.m4 + ex;																	// Invariant mass of the recoil nucleus
	kp.e3_cm = 0.5*( kp.etot_cm*kp.etot_cm + rp.m3*rp.m3 - kp.m4ex*kp.m4ex )/kp.etot_cm;	// Energy of the outgoing proton [CM]
	kp.p3_cm = TMath::Sqrt(kp.e3_cm*kp.e3_cm - rp.m3*rp.m3);								// Momentum of the outgoing proton [CM]
};


#endif

