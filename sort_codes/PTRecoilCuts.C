// PTRecoilCuts.C
// Helps draw the recoil cuts for the ISS data, by continuously updating the pads when the cut has 
// been drawn
// ============================================================================================= //
// Patrick MacGregor
// Nuclear Physics Research Group
// School of Physics and Astronomy
// The University of Manchester
// LAST EDITED: 25/02/19
// ============================================================================================= //
#include <TTree.h>
#include <TFile.h>
#include <TString.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TH2.h>
#include <TCutG.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TROOT.h>
#include <iostream>

// GLOBAL VARIABLES
TString IN_FILE_NAME = "../root_data/fin.root";
TString OUT_FILE_PREFIX = "../root_data/";
Int_t NUM_CUTS = 100;


// --------------------------------------------------------------------------------------------- //
// ADDITIONAL FUNCTIONS
// Draw the E v.s. z plot
void DrawEVZ( TTree *t, TString cutName, Int_t n ){
	// Set up the histogram
	t->Draw( Form( "ecrr:z>>evz%i(401,-50,-10,901,0,9)", n ), cutName.Data(), "goff" );

	// Get the histogram
	TH2F *h = (TH2F*)gDirectory->Get( Form( "evz%i", n ) );

	// Draw appropriately
	h->GetXaxis()->SetTitle("z (cm)");
	h->GetYaxis()->SetTitle("E (MeV)");
	h->SetMarkerStyle(20);
	h->SetMarkerSize(0.5);

	h->Draw();

}


// Draw the excitation spectrum
void DrawEx( TTree *t, TString cutName, Int_t n ){
	// Set up the histogram
	t->Draw( Form( "Ex>>ex%i(451,-1,8)", n ), cutName.Data(), "goff" );

	// Get the histogram
	TH1F *h = (TH1F*)gDirectory->Get( Form( "ex%i", n ) );

	// Draw appropriately
	if ( n == 0 ){
		h->SetFillColor(2);
	}
	h->Draw("same");
}

// --------------------------------------------------------------------------------------------- //
// MAIN FUNCTION
void PTRecoilCuts( Int_t arrayNum = 0 ){

	// Clear the terminal
	gSystem->Exec("clear");
	printf( "========== RECOIL CUTS ON DETECTOR %i ==========\n", arrayNum );
	gStyle->SetOptStat("meni");

	// Get the TFile and the TTree
	TFile *f = new TFile ( IN_FILE_NAME.Data() );
	TTree *t = (TTree*)f->Get("fin_tree");

	// Draw the correct recoil array detector (with a toolbar for cutting)
	TCanvas *cRDT = new TCanvas( "cRDT", "Recoil Detector Plot", 1200, 900 );
	if ( !cRDT->GetShowToolBar() ){
		cRDT->ToggleToolBar();
	}
	t->Draw( Form( "rdt[%i]:rdt[%i]", arrayNum, arrayNum + 4 ) );

	// Define a counter for the number of cuts, and an array to store the cuts
	Int_t numCuts = 0;
	TCutG *cutArray[NUM_CUTS];
	TH2F *h_evz[NUM_CUTS];
	TH1F *h_ex[NUM_CUTS];
	
	// Define a diagnostic canvas to help see values
	TCanvas *cDiag = new TCanvas( "cDiag", "Diagnostic Plots", 935, 1049 );
	cDiag->Divide(1,2);
	
	// Define integral bounds
	Float_t lower, upper;
	lower = -0.18;
	upper = 0.18;
	
	// DO INITIAL START
	cDiag->cd(1); DrawEVZ( t, "", numCuts );
	cDiag->cd(2); DrawEx( t, "", numCuts );
	h_ex[numCuts] = (TH1F*)gDirectory->Get( Form( "ex%i", numCuts ) );
	h_evz[numCuts] = (TH2F*)gDirectory->Get( Form( "evz%i", numCuts ) );
			
	// Print the original integral
	TAxis *xaxis = h_ex[numCuts]->GetXaxis();
	Double_t integral = h_ex[numCuts]->Integral( xaxis->FindBin(lower), xaxis->FindBin(upper) );
	printf("Start: Integral from %f -> %f = %i\n", lower, upper, (Int_t)integral );

	// START THE LOOP
	while ( numCuts < NUM_CUTS ){
		// DRAW THE CUT
		printf("Draw a cut on the Recoil Detector Plot\n");
		cRDT->cd();
		TVirtualPad* padRDT = cRDT->GetPad(0);
		padRDT->WaitPrimitive("CUTG");

		// Cut drawn - get and store the cut
		cutArray[numCuts] = (TCutG*)gROOT->FindObject("CUTG");
		cutArray[numCuts]->SetName( Form( "tempCut%i", numCuts ) );
		
		// UPDATE THE LOOK OF THE CUT
		// Now redraw the canvas correctly, recolouring the cuts etc.	
		if ( numCuts > 2 ){
			cRDT->GetListOfPrimitives()->Remove(cutArray[numCuts-3]);
		}
		if ( numCuts > 1 ){
			cutArray[numCuts-2]->SetLineColor(4);
			cutArray[numCuts-2]->SetLineWidth(2);
		}
		if ( numCuts > 0 ){
			cutArray[numCuts-1]->SetLineColor(2);
			cutArray[numCuts-1]->SetLineWidth(2);
		}
		cutArray[numCuts]->SetLineColor(1);
		cutArray[numCuts]->SetLineWidth(2);
		
		// Update the canvas
		cRDT->Modified(); cRDT->Update();
		
		// Print a friendly message
		printf( "Cut #%i generated. Redrawing graphs...\n", numCuts );
		
		// EVZ PLOT UPDATE
		// Draw the new EVZ plot(and print helpful messages)
		TVirtualPad* padEVZ = cDiag->cd(1);
		DrawEVZ( t, cutArray[numCuts]->GetName(), numCuts+1 );
		h_evz[numCuts+1] = (TH2F*)gDirectory->Get( Form( "evz%i", numCuts+1 ) );
		
		// Remove previous plots and update the canvas
		cDiag->GetListOfPrimitives()->Remove(h_evz[numCuts]);
		padEVZ->Modified();padEVZ->Update();			
		printf("Updated E v.s. z plot\n");
		
		// EX PLOT UPDATE
		// Draw the new excitation spectrum on top of the old (and print helpful messages)
		TVirtualPad* padEx = cDiag->cd(2);
		padEx->Clear();
		h_ex[numCuts]->Draw();
		DrawEx( t, cutArray[numCuts]->GetName(), numCuts+1 );
		h_ex[numCuts+1] = (TH1F*)gDirectory->Get( Form( "ex%i", numCuts+1 ) );
		h_ex[numCuts]->SetFillColor(4);
		h_ex[numCuts+1]->SetFillColor(2);
		padEx->Modified(); padEx->Update();			
		printf("Updated excitation spectrum\n");
		
		// Calculate the integral from the new cut
		xaxis = h_ex[numCuts+1]->GetXaxis();
		integral = h_ex[numCuts+1]->Integral( xaxis->FindBin(lower), xaxis->FindBin(upper) );
		printf("Cut %i: Integral from %f -> %f = %i\n", numCuts, lower, upper, (Int_t)integral );

		// ASK IF THEY WOULD LIKE TO CONTINUE
		TString answer;
		printf( "Would you like to continue [Y/N]: ");
		
		// Check for any errors
		while ( !(std::cin>>answer) ){
			printf("Please enter a valid string\n");
			std::cin.clear();
			std::cin.ignore(10000, '\n');
		}
		// If answer is yes, continue chopping
		if ( answer.CompareTo("Y") == 0 || answer.CompareTo("Yes") == 0 || answer.CompareTo("y") == 0 || answer.CompareTo("yes") == 0 ){
			// Print a friendly message about continuing
			printf("Continuing to make data better ...\n");
			printf("================================================\n");
		}
		// If answer is no i.e. save the file and close
		else{
			// Define the ouput file name and write a congratulatory message
			TString outFileName = OUT_FILE_PREFIX + Form( "cut%i.root", arrayNum );
			printf("Congratulations! Cut %i is your cut that will be saved as \"cut%i\" in file %s\n", numCuts, arrayNum, outFileName.Data() );
			
			// Open a new TFile to store the cut
			TFile *outFile = new TFile( outFileName.Data(), "RECREATE" );
			
			// Rename the final cut to the required cut
			cutArray[numCuts]->SetName( Form( "cut%i", arrayNum ) );
			cutArray[numCuts]->SetTitle("");

			// Store the new cut
			cutArray[numCuts]->Write();

			// Close the file
			if ( outFile != NULL ) outFile->Close();

			// Break out of the loop			
			break;
		}
		
		// Iterate the number of cuts
		numCuts++;
	}
}

































