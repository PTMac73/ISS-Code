// PTFinFunc.C
// A group of functions for manipulating the fin.root files produced by the PTMonitors code

#include <TH2.h>
#include <TH1.h>
#include <TStyle.h>
#include <TCutG.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TString.h>
#include <TSystem.h>
#include <TObjArray.h>
#include <TMath.h>
#include <TFile.h>
#include "PTMonitors.h"

TCanvas *cRDT4;
TCanvas *cArray24;
TString blank = (TString)"";

// --------------------------------------------------------------------------------------------- //
// Draw all the recoil detector cuts on a single plot
void RDT_CUTS( TFile *f ){
	
	// Get the TTree
	TTree *t = (TTree*)f->Get("fin_tree");

	// Generate a new TCanvas
	cRDT4 = new TCanvas("c4", "Recoil detector cuts",1800,900);

	// Divide the TCanvas
	cRDT4->Divide(2,2);

	// Plot the graphs
	for (Int_t i = 0; i < 4; i++ ){
		cRDT4->cd(i+1);
		t->Draw( Form("rdt[%i]:rdt[%i]", i, i+4 ) );
		TCutG *cut = (TCutG*)f->Get( Form( "cut%i", i ) );
		cut->Draw("same");
	}
}
// --------------------------------------------------------------------------------------------- //
// Draw some function for each array detector
void arrayDetectorPlot( TFile *f, const char *str1, TString str2 = blank ){
	// Get the TTree and the cuts
	TTree *t = (TTree*)f->Get("fin_tree");
	TCutG *cut0 = (TCutG*)f->Get("cut0");
	TCutG *cut1 = (TCutG*)f->Get("cut1");
	TCutG *cut2 = (TCutG*)f->Get("cut2");
	TCutG *cut3 = (TCutG*)f->Get("cut3");

	// Now create the TCanvas
	cArray24 = new TCanvas("c24", "Individual array plot", 1920, 1080);
	
	// Divide the TCanvas
	cArray24->Divide(6,4);

	// Plot the graphs
	for (Int_t i = 0; i < 6; i++ ){
		for (Int_t j = 0; j < 4; j++ ){
			cArray24->cd( i + 6*j + 1);
			if ( strncmp( str2.Data() , "", 1 ) != 0 ){
				t->Draw( str1 , Form( "detID == %i && %s", i + 6*j, str2.Data() ) );
			}
			else{
				t->Draw( str1 , Form( "detID == %i", i + 6*j ) );
			}
		}
	}
}

// --------------------------------------------------------------------------------------------- //
// Plot E v.s. z from the fin.root file

// fin_tree->Draw("ecrr:z>>a(401,-50,-10,901,0,9)", "(cut0 || cut1 || cut2 || cut3) && (td_rdt_e > -30 && td_rdt_e < 30)");

// --------------------------------------------------------------------------------------------- //
// Plot timing cuts
void plotTimingCuts( TFile *f, Int_t arrayRow = 0 ){
	
	// Get the TTree
	TTree *t = (TTree*)f->Get("fin_tree");
	
	// Create 6 TCanvases (with 4 slots each) and histogram arrays
	TCanvas *c1;
	TH1F *a[4], *b[4];
	c1 = new TCanvas( Form( "c%i", arrayRow ), Form("Array Row %i", arrayRow), 1200, 900);
	c1->Divide(2,2);

	// Start populating the TCanvases
	for ( Int_t j = 0; j < 4; j++ ){
		// Change to the right pad
		c1->cd(j+1);

		// Draw using the TTree, but don't draw graphically
		t->Draw( Form( "td_rdt_e[%i]>>b%i(61,-30,30)", 6*j + arrayRow, 6*j + arrayRow ), Form( "td_rdt_e[%i] > %i && td_rdt_e[%i] < %i ", 6*j + arrayRow, -30, 6*j + arrayRow, 30 ), "goff" );
		t->Draw( Form( "td_rdt_e[%i]>>a%i(61,-30,30)", 6*j + arrayRow, 6*j + arrayRow ), Form( "td_rdt_e[%i] > td_rdt_e_cuts[%i][0] && td_rdt_e[%i] < td_rdt_e_cuts[%i][1] ", 6*j + arrayRow, 6*j + arrayRow, 6*j + arrayRow, 6*j + arrayRow ), "goff" );

		// Store the histogram
		a[j] = (TH1F*)gDirectory->Get( Form("a%i", 6*j + arrayRow ) );
		b[j] = (TH1F*)gDirectory->Get( Form("b%i", 6*j + arrayRow ) );

		// Draw the histogram
		a[j]->SetFillColor(5);
		b[j]->SetFillColor(2);
		b[j]->Draw();
		a[j]->Draw("same");
	}
	
}

void plotTimingCutsALL( TFile *f ){
	for (Int_t i = 0; i < 6; i++ ){
		plotTimingCuts( f, i );
	}
}

// --------------------------------------------------------------------------------------------- //
// Plot timing cuts
void plotXCAL( TFile *f ){
	
	// Get the TTree
	TTree *t = (TTree*)f->Get("fin_tree");

	// Get the cuts
	TCutG *cut0 = (TCutG*)f->Get("cut0");
	TCutG *cut1 = (TCutG*)f->Get("cut1");
	TCutG *cut2 = (TCutG*)f->Get("cut2");
	TCutG *cut3 = (TCutG*)f->Get("cut3");
	
	// Create 6 TCanvases (with 4 slots each) and histogram arrays
	TCanvas *c1;
	TH1F *a[24], *b[24];
	c1 = new TCanvas( "canvas", "xcal plots", 1350, 900);
	c1->Divide(6,4);

	// Start populating the TCanvases
	for ( Int_t i = 0; i < 6; i++ ){
		for ( Int_t j = 0; j < 4; j++ ){
			// Change to the right pad
			c1->cd(6*j + i + 1);

			// Draw using the TTree, but don't draw graphically
			t->Draw( Form( "xcal[%i]>>a%i(201,-0.5,1.5)", 6*j + i, 6*j + i ), "(cut0 || cut1 || cut2 || cut3) && xcal[] > xcal_cuts[][0] && xcal[] < xcal_cuts[][1]", "goff" );		
			t->Draw( Form( "xcal[%i]>>b%i(201,-0.5,1.5)", 6*j + i, 6*j + i ), "cut0 || cut1 || cut2 || cut3", "goff" );

			// Store the histogram
			a[6*j + i] = (TH1F*)gDirectory->Get( Form("a%i", 6*j + i ) );
			b[6*j + i] = (TH1F*)gDirectory->Get( Form("b%i", 6*j + i ) );

			// Draw the histogram
			b[6*j + i]->SetFillColor(2);
			b[6*j + i]->Draw();
			a[6*j + i]->SetFillColor(5);
			a[6*j + i]->SetTitle( Form( "xcal[%i]", 6*j + i ) );
			a[6*j + i]->Draw("SAME");
		}
	}
	c1->Print("/home/ptmac/Desktop/xcalCuts.pdf");
	
}

// --------------------------------------------------------------------------------------------- //
// Plot changes in the excitation spectrum from xcal
void xcalChanges( TFile *f, Int_t detID = 0 ){
	
	// Get the TTree and the cuts
	TTree *t = (TTree*)f->Get("fin_tree");
	TCutG *cut0 = (TCutG*)f->Get("cut0");
	TCutG *cut1 = (TCutG*)f->Get("cut1");
	TCutG *cut2 = (TCutG*)f->Get("cut2");
	TCutG *cut3 = (TCutG*)f->Get("cut3");
	
	// Define a series of variables
	const int N = 5;
	TH1F *u[N], *l[N];
	TCanvas *uc, *lc;
	Float_t xcal_lb[N] = { 0, 0.1, 0.2, 0.3, 0.4 };
	Float_t xcal_ub[N] = { 1, 0.98, 0.96, 0.94, 0.92 };

	// Set the style
	gStyle->SetOptStat("meni");
	
	// Draw and store the histograms
	for (Int_t i = 0; i < N; i++ ){
		// Draw
		t->Draw( Form( "Ex>>l%i(451, -1, 8 )", i ), Form( "( cut0 || cut1 || cut2 || cut3 )  &&  thetaCM > 11  &&  td_rdt_e[] > td_rdt_e_cuts[][0]  &&  td_rdt_e[] < td_rdt_e_cuts[][1]  &&  xcal[%i] >= %f  &&  xcal[%i] <= 1", detID, xcal_lb[i], detID ), "goff" );
		t->Draw( Form( "Ex>>u%i(451, -1, 8 )", i ), Form( "( cut0 || cut1 || cut2 || cut3 )  &&  thetaCM > 11  &&  td_rdt_e[] > td_rdt_e_cuts[][0]  &&  td_rdt_e[] < td_rdt_e_cuts[][1]  &&  xcal[%i] >= 0  &&  xcal[%i] <= %f", detID, detID, xcal_ub[i] ), "goff" );
		
		// Store
		l[i] = (TH1F*)gDirectory->Get( Form("l%i", i ) );
		u[i] = (TH1F*)gDirectory->Get( Form("u%i", i ) );

		// Set colours
		l[i]->SetFillColor(i+1);
		u[i]->SetFillColor(i+1);
	}

	// Now draw the histograms
	uc = new TCanvas( "uc", Form( "Upper (detID = %i)", detID ), 1200, 900 );
	for (Int_t i = 0; i < N; i++ ){
		u[i]->Draw("same");
	}

	lc = new TCanvas( "lc", Form( "Lower (detID = %i)", detID ), 1200, 900 );
	for (Int_t i = 0; i < N; i++ ){
		l[i]->Draw("same");
	}
}




















