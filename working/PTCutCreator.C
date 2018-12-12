#include <TH2F.h>
#include <TFile.h>
#include <TChain.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TCutG.h>
#include <TString.h>
#include <TObjArray.h>

// Define global variables
TString inFileName = "../root_data/genPos1.root";
TString outFileName = "../working/ALL-MgCuts.root";

// MAIN FUNCTION
void PTCutCreator(){
	
	// Print message
	printf("================ Graphic Cut Creator for RDT ============== \n");
   
	// Get the file and the TChain
	TFile *f = new TFile( inFileName.Data() );
	TChain *chain = (TChain*)f->Get("gen_tree");

	// Print the list of files in the TChain
   	chain->GetListOfFiles()->Print();
   
	// Define some TStrings
	TString varX, varY, tag;
	
	// Style the graphs to be plotted
	gStyle->SetOptStat(00000);
	
	// Open a new canvas and ensure the toolbar is visible.
	TCanvas * cCutCreator = new TCanvas("cCutCreator", "RDT Cut Creator", 100, 100, 800, 800);
	if( !cCutCreator->GetShowToolBar() ) cCutCreator->ToggleToolBar();
	
   	// Open the out file
   	TFile * cutFile = new TFile( outFileName.Data(), "recreate" );
	cCutCreator->Update();
	
	// Define the cuts
	TCutG * cut = NULL;
	TObjArray * cutList = new TObjArray();
	TString expression[10];

	// Loop over the number of graphs to apply cuts to
	for (Int_t i = 0; i < 4; i++) {
		// Print a message
		printf("======== make a graphic cut on the plot, %d-th cut: \n", i );

		// Define the names of the histograms
		varX.Form("rdt[%d]",i+4); varY.Form("rdt[%d]",i);

		// Define the histogram plotting options - 2096 bins, from 0 to 8192
		expression[i].Form("%s:%s>>h(2096, 0, 8192, 2096, 0, 8192)", 
            varY.Data(),
            varX.Data());

		// Draw the histogram
		chain->Draw(expression[i], "", "col");
		
		// Update the canvas
		cCutCreator->Modified(); cCutCreator->Update();
		gPad->WaitPrimitive("CUTG");

		// Find the cut after it's been made
		cut = (TCutG*) gROOT->FindObject("CUTG");

		// Store the cut
		TString name; name.Form("cut%d", i);
		cut->SetName(name);
		cut->SetVarX(varX.Data());
		cut->SetVarY(varY.Data());
		cut->SetTitle(tag);
		cut->SetLineColor(1);
		cutList->Add(cut);

		// Print a message of congratulations
		printf(" cut-%d \n", i);
	}

	// Write the cut list to a file
	cutList->Write("cutList", TObject::kSingleKey);
	
	// Print exit message
	printf("====> saved %d cuts into %s\n", 4, outFileName.Data() );
}
