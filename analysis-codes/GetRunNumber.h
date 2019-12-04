// GetRunNumber.h
// Gets the run number from a root file named [PREFIX]123.root
// ============================================================================================= //
// Patrick MacGregor
// Nuclear Physics Research Group
// School of Physics and Astronomy
// The University of Manchester
// ============================================================================================= //
#ifndef GET_RUN_NUMBER_H_
#define GET_RUN_NUMBER_H_

#include <TString.h>

Int_t GetRunNumber( TString file_name ){
	// Strip the suffix
	TString file_ending = ".root";
	file_name.Remove( file_name.Length() - file_ending.Length(), file_ending.Length() );

	// Loop over the rest from the end and see if it's a number
	Int_t num_counter = 0;

	for ( Int_t i = file_name.Length() - 1; i >= 0; i-- ){
		if ( std::isdigit( file_name[i] ) != 0 ){
			num_counter++;
		}
		else{
			break;
		}

	}

	// Return the number	
	file_name.Remove( 0, file_name.Length() - num_counter );
	Int_t run_number = std::stoi( file_name.Data() );
	return run_number;
}


#endif
