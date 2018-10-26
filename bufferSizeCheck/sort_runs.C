// sort_runs.C
// Accompanies the shell script sort_runs.sh
// ============================================================================================== //
// Patrick MacGregor
// Nuclear Physics Research Group
// School of Physics and Astronomy
// The University of Manchester
// LAST EDITED: 24/10/18
// ============================================================================================== //
TString sortDir = "/home/ptmac/Documents/07-CERN-ISS-Mg/analysis/sort_codes";

void sort_runs(TString fileDir = "0"){
	if ( fileDir != "0" ){
		// Open the TFile
		TFile *f = new TFile( fileDir.Data() );
	
		// Get the TTree
		TTree *t = (TTree*)f->Get("gen_tree");

		// Process the TTree
		t->Process( Form( "%s/PTMonitors.C++", sortDir.Data() ) );

		// Close the TFile
		f->Close();
	}
}
