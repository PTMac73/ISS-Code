// Doublet_ImportData.h
// Imports the data for the doublet
// Data format is:
// X [TAB] EXP [TAB] L=0 [TAB] L=1 [TAB] L=2 [TAB] L=3
// ============================================================================================= //
// Patrick MacGregor
// Nuclear Physics Research Group
// School of Physics and Astronomy
// The University of Manchester
// ============================================================================================= //
#ifndef DOUBLET_IMPORT_DATA_H_
#define DOUBLET_IMPORT_DATA_H_

#include <iostream>
#include <fstream>
#include <TString.h>
#include "Doublet_Globals.h"

void ImportData( TString in_file_dir ){
	// Open the file
	std::ifstream file_in;
	file_in.open( in_file_dir );
	
	// Check if the file is open
	if ( !file_in.is_open() ){
		std::cout << "File did not open." << std::endl;
		exit(1);
	}
	 
	// Store the data conveniently [ N.B. the L values are 0,1,2,2,3,4 ]
	Int_t i = 0;
	Double_t temp1;
	Double_t temp2;
	while ( true ){
		file_in >> X[i] >> Y[i] >> E[i] >> PT[0][i] >> PT[1][i] >> temp1 >> PT[2][i] >> PT[3][i] >> temp2;
		if ( !file_in.good() ){ break; }
		i++;
	}
	
	// Close the file
	file_in.close();
}

#endif
