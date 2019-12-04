// PTCreateCutFile.C
// Combines the 4 cut files into 1 file that can be used in the PTMonitors.C script
// ============================================================================================= //
// Patrick MacGregor
// Nuclear Physics Research Group
// School of Physics and Astronomy
// The University of Manchester
// ============================================================================================= //
#include <TFile.h>
#include <TString.h>
#include <TObjArray.h>
#include <TCutG.h>

TString DIR_PREFIX = "../root_data/cut";

void PTCreateCutFile( TString outputName = "../working/rootcutFile.root" ){
	// Define a TObjArray to hold the cuts
	TObjArray *cutList = new TObjArray();
	
	// Loop over the four cuts and store the cut array
	for ( Int_t i = 0; i < 4; i++ ){
		TFile *cutFile = new TFile( Form( "%s%i.root", DIR_PREFIX.Data(), i ), "r" );
		TCutG *cut = (TCutG*)cutFile->Get( Form( "cut%i", i ) );
		cutList->Add(cut);
		if ( cutFile != NULL ) cutFile->Close();
	}

	// Write the four cuts to file
	TFile *outFile = new TFile( outputName.Data(), "RECREATE" );	
	cutList->Write("cutList", TObject::kSingleKey );	
	
	// Close the final cut file
	if ( outFile != NULL ) outFile->Close();
}
