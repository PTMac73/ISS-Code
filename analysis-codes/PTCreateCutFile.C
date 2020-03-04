// PTCreateCutFile.C
// Combines the 4 cut files into 1 file that can be used in the PTMonitors.C script
// ============================================================================================= //
// Patrick MacGregor
// Nuclear Physics Research Group
// Department of Physics and Astronomy
// The University of Manchester
// ============================================================================================= //
#include <TFile.h>
#include <TString.h>
#include <TObjArray.h>
#include <TCutG.h>

TString DIR_PREFIX = "../root_data/cut";

void PTCreateCutFile( TString outputName = "../working/custard.root" ){
	// Define a TObjArray to hold the cuts
	TObjArray *cut_list = new TObjArray();
	
	// Loop over the four cuts and store the cut array
	for ( Int_t i = 0; i < 4; i++ ){
		TFile *cut_file = new TFile( Form( "%s%i.root", DIR_PREFIX.Data(), i ), "r" );
		TCutG *cut = (TCutG*)cut_file->Get( Form( "cut%i", i ) );
		cut_list->Add(cut);
		if ( cut_file != NULL ) cutFile->Close();
	}

	// Write the four cuts to file
	TFile *out_file = new TFile( outputName.Data(), "RECREATE" );	
	cut_list->Write("cut_list", TObject::kSingleKey );	
	
	// Close the final cut file
	if ( out_file != NULL ) out_file->Close();
}
