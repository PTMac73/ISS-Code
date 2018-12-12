#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include <vector>

// GLOBAL VARIABLES
TChain *ch0, *ch1, *ch2, *chAlpha;
TFile *f0, *f1, *f2, *fAlpha;

// Run arrays (Array position - offset = 70 mm )??
vector <Int_t> runArrayAlpha = {0};
vector <Int_t> runArray0 = {50, 51, 52, 53, 54, 55}; // Array position = 119.765 => z_off = 3.00235
vector <Int_t> runArray1 = {56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 82, 83, 84, 85, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100}; // Array position = 164.98 mm => z_off = 9.498 cm
// Run 103 is the run with the 600 ug/cm2 target - use if you want to
vector <Int_t> runArray2 = {104, 105, 106, 107, 108, 110, 111, 112, 113, 115, 116, 117, 118}; // Array position = 135.0 mm => z_off = 6.5 cm

TString prefix = "../root_data/gen_run";
TString suffix = ".root";


// BEGIN FUNCTION THAT CHAINS THE RUNS TOGETHER
void chainRuns( vector <Int_t> runArray, TChain *ch, TFile *f, TString outName ){
	// Define the new tree
	ch = new TChain("gen_tree");
	
	// Delcare an iterator
	vector <Int_t>::iterator ptr;
	
	// Loop over the runs
	for ( ptr = runArray.begin(); ptr < runArray.end(); ptr++ ){
		
		// Generate the file name
		TString fileName;
		fileName += prefix; fileName += *ptr; fileName +=suffix;
		
		// Add the root files to the chain
		ch->AddFile( fileName.Data() );

	}
	
	// Open a file for writing
	f = new TFile( outName.Data(), "NEW" );
	
	// Check that it is open
	if ( f->IsOpen() ){
		Printf("File is open");
		ch->Write("gen_tree");
		f->Close();
		Printf("File is closed");
	}
}



void SharpChain(){
	// TChain the runs together
	//gDirectory->Delete("ch0");
	gDirectory->Delete("ch1");
	gDirectory->Delete("ch2");
	//gDirectory->Delete("chAlpha");
	
	//chainRuns( runArray0, ch0, f0, "../root_data/genPos0.root" );
	chainRuns( runArray1, ch1, f1, "../root_data/genPos1.root" );
	chainRuns( runArray2, ch2, f2, "../root_data/genPos2.root" );
	//chainRuns( runArrayAlpha, chAlpha, fAlpha, "../root_data/genAlpha.root" );
}
