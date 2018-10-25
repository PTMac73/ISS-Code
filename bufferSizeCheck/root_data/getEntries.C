// Writes the .chat files for the ISS data to check the best buffer size
// ============================================================================================== //
// Patrick MacGregor
// Nuclear Physics Research Group
// School of Physics and Astronomy
// The University of Manchester
// LAST EDITED: 24/10/18
// ============================================================================================== //
#include <vector>
#include "TFile.h"
#include "TTree.h"
#include "TString.h"

// GLOBAL VARIABLES
// bigbufsize obtained from the bufferWrite.py script
vector <Int_t> bigbufsize = {450, 5400, 10350, 15300, 20250, 25200, 30150, 35100, 40050, 45000};
TString prefix = "gen_run51-";
TString suffix = ".root";
TString fileName, sBBNum;
TFile *f;
TTree *t;
Int_t numEntries;

// MAIN FUNCTION
void getEntries(){
	// Loop over number of root files
	for ( Int_t i = 0; i < bigbufsize.size(); i++ ){
		// Make the file name
		sBBNum.Form( "%d", bigbufsize[i] );
		fileName = prefix + sBBNum + suffix;
		
		// Open the file
		f = new TFile( fileName.Data() );
		
		// Check the file has opened
		if ( f->IsOpen() ){
			// Get the number of entries in the TTree, and print next to the buffer size
			t = (TTree*)f->Get("gen_tree");
			
			numEntries = t->GetEntries();
			
			Printf( "%s & %i\\\\\n", sBBNum.Data(), numEntries );
			
			// Close the file
			f->Close();
			
		}
		else{
			continue;
		}
	}
}
