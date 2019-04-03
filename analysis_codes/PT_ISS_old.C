// Analyses all of the data at different array positions, and then spits out an E vs. z graph
#include "PTMonitors.h"
#include <TStyle.h>
#include <TCutG.h>
#include <TString.h>
#include <TMath.h>

// GLOBAL VARIABLES
// TTree to contain every single energy and position
TTree *full_tree;		// Final TTree
TBranch *bEnergy;		// Branches to store energy and position
TBranch *bPosition;
Float_t Z_final[24];	// Floats to hold energy and position
Float_t E_final[24];

TTree *t[70];			// TTree array to hold each TTree from the files
TList *treeList;		// TList to list the TTrees for merging


// Things pertinent for the run files
TString filePrefix = "/home/ptmac/Documents/07-CERN-ISS-Mg/analysis/root_data/gen_run";
TString fileSuffix = ".root";
const Int_t NexcludeRuns = 6;
Int_t excludeRuns[NexcludeRuns] = {81,86,101,102,109,114};
Bool_t skipRun;

// Offset of array position
// 6.5 is 
Float_t Z_off[2]={3.502,6.5};	// Distance to the physical end of array from the target (-20 to si start)
Float_t z_off[70];				// Holds the position of each run
Long64_t numEntries[70];		// Number of entries in each tree


// MAIN ACTION --------------------------------------------------------------------------------- //
void PT_ISS(){
	// Define the TTree
	full_tree = new TTree("full_tree", "Super TTree containing stuff");
	// Set the branch address for the TTree
	/*full_tree->SetBranchAddress("E", E_final, &bEnergy);
	full_tree->SetBranchAddress("Z", Z_final, &bPosition);*/
	
	// LOOP OVER THE FILES
	for (Int_t i = 50; i < 53/*119*/; i++){
		// Reset the skip run flag
		skipRun = 0;
	
		// Reject the run if in the excludeRuns array
		for (Int_t j = 0; j < NexcludeRuns; j++ ){
			if ( i == excludeRuns[j] ){
				skipRun = 1;
				break;
			}
		}
		
		// Now check whether to do the run
		if (skipRun == 1){
			continue;
		}
		else{
			// Define the file name
			TString sRunNum; sRunNum.Form("%d",i);
			TString fileName = filePrefix + sRunNum + fileSuffix;
			
			// Open the TFile
			TFile *f = new TFile(fileName);
			
			// Check if the file is open
			if ( f->IsOpen() ){
				Int_t runIndex = i-50;
				
				// Get the TTree and the number of entries
				t[runIndex] = (TTree*)f->Get("gen_tree");
				numEntries[runIndex] = t[runIndex]->GetEntries();
				
				// Define the array offset
				if ( i < 101 ){
					// Append the first array position
					z_off[runIndex] = Z_off[0];
				}
				else{
					// Append the first array position
					z_off[runIndex] = Z_off[1];
				}
				// Add the TTree to the TList
				treeList->Add(t[runIndex]);
				
				// Close the file
				f->Close();
			}
		}
		std::cout << numEntries << std::endl;
	}
	
	// Now run the TSelector!
	
	
	// Define the files to be analysed at a particular array position
	
	// Extract the relevant quantities for each run using a TSelector on each root file
	
	// Plot a super graph!
	// Make a new TStyle for the output graph
	TStyle *sharpyStyle = new TStyle("sharpyStyle","David Sharp's Style");
	sharpyStyle->SetOptStat(0);
	sharpyStyle->SetCanvasBorderMode(0);
	sharpyStyle->SetPadBorderMode(0);
	sharpyStyle->SetPadColor(0);
	sharpyStyle->SetCanvasColor(0);
	sharpyStyle->SetStatColor(0);
	sharpyStyle->SetTitleFillColor(0);
	sharpyStyle->SetTitleBorderSize(0);
	sharpyStyle->SetTitleAlign(23);
	sharpyStyle->SetTitleX(0.5);
	sharpyStyle->SetMarkerStyle(6);
	
	
}
