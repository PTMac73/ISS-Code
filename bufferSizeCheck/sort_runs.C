// sort_runs.C
// Accompanies the shell script sort_runs.sh
// ============================================================================================== //
// Patrick MacGregor
// Nuclear Physics Research Group
// School of Physics and Astronomy
// The University of Manchester
// LAST EDITED: 26/10/18
// ============================================================================================== //
TString sortDir = "/home/ptmac/Documents/07-CERN-ISS-Mg/analysis/sort_codes";

void sort_runs(TString fileDir = "0", Int_t mode = 0){
	if ( fileDir != "0" ){
		// 0 - PROCESS THE GEN_RUN FILES
		if ( mode == 0 ){
			// Open the TFile
			TFile *f = new TFile( fileDir.Data() );

			// Get the TTree
			TTree *t = (TTree*)f->Get("gen_tree");

			// Process the TTree
			t->Process( Form( "%s/PTMonitors.C++", sortDir.Data() ) );
	
			// List objects in memory
			//gDirectory->ls();

			// Close the TFile
			f->Close();
		}
		// 1 - PROCESS THE FIN_RUN FILES
		else if ( mode == 1 ){
			// Open the TFile
			TFile *f = new TFile( fileDir.Data() );

			// Get the two histograms
			TH2F *h0 = (TH2F*)f->Get("EVZ");
			TH1F *h1 = (TH1F*)f->Get("EXE");

			// Print to a table line
			Printf("%s\t%0.0f\t%0.0f\n", fileDir.Data(), h0->GetEntries(), h1->GetEntries() );

			// Close the TFile
			f->Close();
		}
		// DEFAULT - PRINT ERROR MESSAGE
		else{
			Printf("\nError. Mode is not an allowed value!\n\n");
		}
	}
}
