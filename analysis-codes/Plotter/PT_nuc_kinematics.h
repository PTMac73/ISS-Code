// PT_nuc_kinematics.h
// Kinematics in the ISS and general nuclear physics formulae
// ============================================================================================= //
// Patrick MacGregor
// Nuclear Physics Research Group
// School of Physics and Astronomy
// The University of Manchester
// ============================================================================================= //
#ifndef PT_NUC_KINEMATICS_H_
#define PT_NUC_KINEMATICS_H_

#include <TMath.h>

// GLOBAL VARIABLES
Double_t AMU = 931.494;
Double_t ELECTRIC_CHARGE = 1.602E-19;

// Converts mass excess to mass in u
Double_t MassExcessToMass( Int_t A, Double_t mass_excess ){
	return (A*AMU + mass_excess)/AMU;
}

// Converts Joules to MeV
Double_t JoulesToMeV( Double_t energyJ ){
	return energyJ/( ELECTRIC_CHARGE*(1e6) );
}

// Converts AMU to kg
Double_t AMUToKg( Double_t mass ){
	return mass*AMU*ELECTRIC_CHARGE*(1e6)/( TMath::C()*TMath::C() );
}

Double_t MToMM(){
	return 1000.0;
}

Double_t MToCM(){
	return 100.0;
}


// Calculates the average radius of the ISS array (be consistent with units!!)
/*/ The X-axis has its origin in the centre of the array (not necessarily the centre of the strip)
 *  X1 is the point furthest left of the origin
 *  X2 is the point furthest right of the origin
 *  height is the height of the array                                                            */
Double_t ISSArrayRadius( Double_t X1, Double_t X2, Double_t height ){
	Double_t av_rad, term_1, term_2, term_3;
	term_1 = height*(X2-X1)/2;
	term_2 = ( TMath::Power(X2,3) - TMath::Power(X1,3) )/(2*height);
	term_3 = ( height*height/2.0 )*TMath::Log( TMath::Abs( ( height*height + X2*height + X2*X2 )/( height*height + X1*height + X1*X1 ) ) );
	av_rad = ( term_1 + term_2 + term_3 )/( X2 - X1 );
	return av_rad;
}


#endif
