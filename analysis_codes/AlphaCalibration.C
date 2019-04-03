// AlphaCalibration.C
// Performs the alpha calibration for a given position
// ============================================================================================= //
// Patrick MacGregor
// Nuclear Physics Research Group
// School of Physics and Astronomy
// The University of Manchester
// ============================================================================================= //
#include <TFile.h>
#include <TTree.h>
#include <TH2.h>
#include <TDirectory.h>
#include <TCanvas.h>
#include <TStyle.h>


// GLOBAL VARIABLES
Bool_t SWITCH_DISPLAY_CANVAS = 1;
Bool_t SWITCH_PRINT_CANVAS = 0;


// MAIN FUNCTION
void AlphaCalibration( TFile *f ){
	
	// Get the TTree
	TTree *t = (TTree*)f->Get("fin_tree");
	
	// PLOT THE HISTOGRAM
	// Define bins
	Int_t hist_bins[2] = { 400, 900 };
	Double_t hist_limits_x[2] = { -50, -10 };
	Double_t hist_limits_y[2] = { 0, 9 };
	
	// Call a canvas
	TCanvas *c_EVZ = new TCanvas( "c_EVZ", "Canvas for plotting the alpha runs", 1200, 900 );
	
	// Draw the histogram
	t->Draw( Form( "ecrr:z>>h_EVZ(%i, %f, %f, %i, %f, %f)", hist_bins[0], hist_limits_x[0], hist_limits_x[1], hist_bins[1], hist_limits_y[0], hist_limits_y[1] ), "xcal[] > xcal_cuts[][0] && xcal[] < xcal_cuts[][1]", "goff" );
	
	// Retrieve the histogram
	TH2F* h_EVZ = (TH2F*)gDirectory->Get("h_EVZ");
	
	// FORMAT THE HISTOGRAM
	h_EVZ->SetTitle("Alpha Data with position cuts");
	h_EVZ->GetXaxis()->SetTitle("z / cm");
	h_EVZ->GetYaxis()->SetTitle("E / MeV");
	gStyle->SetOptStat("meni");
	h_EVZ->Draw("colz");
	gStyle->SetPalette(kRainBow);
	c_EVZ->Modified(); c_EVZ->Update();
	if ( SWITCH_PRINT_CANVAS == 1 ){
		c_EVZ->Print("../PLOTS/AlphaCuts.pdf");
	}
	
	
	// GET THE NUMBER OF COUNTS
	
	
	
	
	
	
	
}
