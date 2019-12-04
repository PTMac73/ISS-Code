// PTF_GetPosNumber.h
// Initialises all the variables for the PTPlotter.C script
// ============================================================================================= //
// Patrick MacGregor
// Nuclear Physics Research Group
// School of Physics and Astronomy
// The University of Manchester
// ============================================================================================= //
#ifndef PTF_GET_POS_NUMBER_H_
#define PTF_GET_POS_NUMBER_H_

#include <TString.h>

// --------------------------------------------------------------------------------------------- //
// Gets the position number from a file name (character before the . in "XYZ-2.root")
Int_t GetPosNumber( TString s ){
	TString file_ending = ".root";
	Int_t pos_number;

	// Remove the file suffix
	s.Remove( s.Length() - file_ending.Length(), file_ending.Length() );
	
	// Remove all but the last character
	s.Remove( 0, s.Length() - 1 );

	// Calculate the number
	pos_number = s.Atoi();

	// Return it
	return pos_number;
}

#endif
