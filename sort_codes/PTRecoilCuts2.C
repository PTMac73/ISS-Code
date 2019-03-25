// PTRecoilCuts2.C
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
#include <TGraph.h>
#include <iostream>

// GLOBAL VARIABLES
TString IN_FILE_NAME = "../root_data/fin.root";
TString OUT_FILE_PREFIX = "../root_data/";
Int_t NUM_CUTS = 100;
Int_t NUM_BINS = 450;


// --------------------------------------------------------------------------------------------- //
// ADDITIONAL FUNCTIONS
// Draw the excitation spectrum

void DrawEx( TTree *t, TString cutName, Int_t n ){
	// Set up the histogram
	t->Draw( Form( "Ex>>ex%i(%i,-1,8)", n, NUM_BINS ), cutName.Data(), "goff" );

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
void PTRecoilCuts2( Int_t arrayNum = 0 ){

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
	t->Draw( Form( "rdt[%i]:rdt[%i]", arrayNum, arrayNum + 4 ), "Ex > -1 && Ex <= 8" );

	// Define a counter for the number of cuts, and an array to store the cuts
	Int_t numCuts = 0;
	TCutG *cutArray[NUM_CUTS];
	TH1F *h_ex[NUM_CUTS];
	TGraph *g_diff[NUM_CUTS];
	Double_t X[NUM_CUTS][NUM_BINS];
	Double_t Y[NUM_CUTS][NUM_BINS];
	
	// Define a diagnostic canvas to help see values
	TCanvas *cDiag = new TCanvas( "cDiag", "Diagnostic Plots", 935, 1049 );
	cDiag->Divide(1,2);
	
	// Define integral bounds
	Float_t lower, upper;
	lower = -0.18;
	upper = 0.18;
	
	// DO INITIAL START
	cDiag->cd(1); DrawEx( t, "", numCuts );
	h_ex[numCuts] = (TH1F*)gDirectory->Get( Form( "ex%i", numCuts ) );
			
	// Print the original integral
	TAxis *xaxis = h_ex[numCuts]->GetXaxis();
	Double_t integral = h_ex[numCuts]->Integral( xaxis->FindBin(lower), xaxis->FindBin(upper) );
	printf("Start: Int[%2.2f,%2.2f] = %i\n", lower, upper, (Int_t)integral );

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
		
		// EX PLOT UPDATE
		// Draw the new excitation spectrum on top of the old (and print helpful messages)
		TVirtualPad* padEx = cDiag->cd(1);
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
		printf(" cut%i: Int[%2.2f,%2.2f] = %i\n", numCuts, lower, upper, (Int_t)integral );

		// HIST DIFFERENCE PLOT UPDATE
		for ( Int_t i = 0; i < h_ex[numCuts]->GetNbinsX(); i++ ){
			X[numCuts][i] = h_ex[numCuts]->GetBinCenter(i);
			if ( i == 0 ){
				Y[numCuts][i] = 0;
			} else{
				Y[numCuts][i] = h_ex[numCuts+1]->GetBinContent(i) - h_ex[numCuts]->GetBinContent(i);
			}
		}
		g_diff[numCuts] = new TGraph( NUM_BINS, X[numCuts], Y[numCuts] );
		TVirtualPad* padDiff = cDiag->cd(2);
		padDiff->Clear();
		g_diff[numCuts]->Draw();
		g_diff[numCuts]->GetXaxis()->SetLimits(-1,8);
		padDiff->Modified(); padDiff->Update();
		printf("Histogram differences plotted\n");



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
			printf("Continuing to have fun :D ...\n");
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


/* #### TODO ####
 * Remove the E v.s. z plot
 * Add a histogram that calculates the difference between the two histograms.
 * Add the current excitation spectrum only
 * On the sole Ex spectrum, fit the first few states with some kind of linear background and display it.
 * Add the thetaCM condition and the timing cuts into play.
*/































