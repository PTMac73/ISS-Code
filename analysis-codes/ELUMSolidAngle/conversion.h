// conversion.h
// Converts units
// ============================================================================================= //
// Patrick MacGregor
// Nuclear Physics Research Group
// School of Physics and Astronomy
// The University of Manchester
// ============================================================================================= //
#ifndef CONVERSION_H_
#define CONVERSION_H_

#include <TMath.h>

// GLOBAL VARIABLES
Double_t AMU = 931.494;
Double_t ELECTRIC_CHARGE = 1.602E-19;

// Converts Joules to MeV
Double_t JoulesToMeV( ){
	return 1.0/( ELECTRIC_CHARGE*(1e6) );
}
// Converts mass excess to mass in u
Double_t MassExcessToMass( Int_t A, Double_t mass_excess, Double_t q ){
	return (A*AMU + mass_excess - q*0.511)/AMU;
}
// Converts AMU to kg
Double_t AMUToKg( ){
	return AMU*ELECTRIC_CHARGE*(1e6)/( TMath::C()*TMath::C() );
}
// Converts metres to millimetres
Double_t MToMM(){
	return 1000.0;
}
// Converts millimetres to metres
Double_t MMToM(){
	return 0.001;
}
// Converts metres to centimetres
Double_t MToCM(){
	return 100.0;
}
// Converts units of e to Coulombs
Double_t EToCoulombs(){
	return ELECTRIC_CHARGE;
}

#endif

