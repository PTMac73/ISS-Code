// PT_ElumDetectorCounts.C
// Fit the deuteron scattering peak
// ============================================================================================= //
// Patrick MacGregor
// Nuclear Physics Research Group
// School of Physics and Astronomy
// The University of Manchester
// ============================================================================================= //
#include "PTF_GetPosNumber.h"
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TDirectory.h>
#include <TF1.h>
#include <TMath.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TString.h>
#include <iostream>

// DEFINE FUNCTIONS
void checkCounts( TTree *t );
Double_t bg( Double_t *x, Double_t *par );
Double_t gaussian( Double_t *x, Double_t *par );
Double_t fitf( Double_t *x, Double_t *par );

// MAIN FUNCTION
void PT_ElumDetectorCounts( TFile *f ){
	
	// Define the array-position number
	Int_t pos_num = GetPosNumber( f->GetName() );

	// Get the TTree
	TTree *t = (TTree*)f->Get("fin_tree");

	// Remove stats box
	gStyle->SetOptStat(0);

	// Check the number of counts
	//checkCounts(t);

	// Draw the full elum detector signal
	TCanvas *c_elum = new TCanvas( "c_elum", "ELUM detector", 1200, 900 );
	t->Draw("elum>>h_elum(2450, -800, 9000)");
	TH1F *h_elum = (TH1F*)gDirectory->Get("h_elum");

	// Change viewing limits and style histogram
	h_elum->GetXaxis()->SetRangeUser(-400,1600);
	h_elum->GetYaxis()->SetRangeUser(-20,400);
	h_elum->GetYaxis()->SetTitle("#");
	h_elum->GetXaxis()->SetTitle("Relative Energy");
	h_elum->SetTitle( Form( "(d,d) Scattering in ELUM (Pos %i)", pos_num ) );

	// Now define a fit - Gaussian with a quadratic background
	Double_t d_lb = 300;
	Double_t d_ub = 750;
	TF1 *f_fitting = new TF1("f_fitting", fitf, d_lb, d_ub, 6 );

	// Set the parameter limits
	/*   A */ f_fitting->SetParLimits(0, 0, 1000);
	/*  mu */ f_fitting->SetParLimits(1, 450, 600);
	/* sig */ f_fitting->SetParLimits(2, 0, 500);

	// Style the fit
	f_fitting->SetLineWidth(4);
	f_fitting->SetLineColor(kRed);

	// Fit the function
	h_elum->Fit( "f_fitting", "R" );

	// Draw the background and Gaussian separately
	TF1 *f_bg = new TF1("f_bg", bg, 0, 1400, 3);
	f_bg->SetLineColor(kBlue);
	f_bg->SetLineWidth(2);

	TF1 *f_gaussian = new TF1("f_gaussian", gaussian, 0, 1400, 3);
	f_gaussian->SetLineColor(kGreen+3);
	f_gaussian->SetLineWidth(2);

	// Get the fit parameters
	Double_t par[6], *e_par;
	f_fitting->GetParameters(par);
	f_gaussian->SetParameters(par);
	f_gaussian->Draw("same");
	
	f_bg->SetParameters(&par[3]);
	f_bg->Draw("same");

	// Add a legend
	TLegend *leg = new TLegend(0.65,0.65,0.9,0.9);
	leg->SetTextFont(42);
	leg->SetTextSize(0.04);
	leg->AddEntry(h_elum,"ELUM Signal","l");
	leg->AddEntry(f_bg,"Background fit","l");
	leg->AddEntry(f_gaussian,"Gaussian fit","l");
	leg->AddEntry(f_fitting,"Total Fit","l");
	leg->Draw();


	// Print window
	printf("===============================================================================================\n");
	printf(" GAUSSIAN FIT:\n");
	printf(" Amp. =\t%5.4e\t +/- \t%5.4e\n", par[0], f_fitting->GetParError(0) );
	printf(" Mu   =\t%5.4e\t +/- \t%5.4e\n", par[1], f_fitting->GetParError(1) );
	printf(" Sig. =\t%5.4e\t +/- \t%5.4e\n", par[2], f_fitting->GetParError(2) );
	printf(" Integral =\t %5.4e\n\n", f_gaussian->Integral(0,1000));

	printf(" BACKGROUND FIT (B2*x*x + B1*x + B0):\n");
	printf(" * B2 =\t%5.4e\t +/-\t%5.4e\n", par[3], f_fitting->GetParError(3) );
	printf(" * B1 =\t%5.4e\t +/-\t%5.4e\n", par[4], f_fitting->GetParError(4) );
	printf(" * B0 =\t%5.4e\t +/-\t%5.4e\n", par[5], f_fitting->GetParError(5) );
	printf("===============================================================================================\n");

	// Print the canvas
	c_elum->Print( Form( "../PLOTS/elum_plot%i.pdf", pos_num ) );


}
// --------------------------------------------------------------------------------------------- //
// FITTING FUNCTIONS
// Background
Double_t bg( Double_t *x, Double_t *par ){
	Double_t B2 = par[0];
	Double_t B1 = par[1];
	Double_t B0 = par[2];
	Double_t fitval = B2*x[0]*x[0] + B1*x[0] + B0;
	return fitval;
}

// Gaussian fit
Double_t gaussian( Double_t *x, Double_t *par ){
	Double_t A = par[0];
	Double_t mu = par[1];
	Double_t sig = par[2];
	Double_t exponent = 0;
	if (sig != 0 ){ exponent = - 0.5*TMath::Power( ( x[0] - mu )/sig, 2 ); }
	Double_t fitval = A*TMath::Exp( exponent );
	return fitval;
}

// Combination
Double_t fitf( Double_t *x, Double_t *par ){
	return gaussian(x, par) + bg(x, &par[3]);
}

// --------------------------------------------------------------------------------------------- //
// Checks the number of counts in each channel of the elum detector
void checkCounts( TTree *t ){
	// Variables
	TH1F *h_array[32];
	Double_t integral[32];
	Double_t sum = 0.0;

	// Loop over the number of array positions
	for ( Int_t i = 0; i < 33; i++ ){
		// Draw the histogram, store it, and get the number of counts.
		if (i == 32 ){
			t->Draw( Form( "elum>>hist%i", i ), "", "goff" );
		}
		else {
			t->Draw( Form( "elum[%i]>>hist%i", i, i ), "", "goff" );
		}
		
		h_array[i] = (TH1F*)gDirectory->Get( Form("hist%i", i ) );
		Int_t Nbins = h_array[i]->GetNbinsX();
		integral[i] = h_array[i]->Integral(0, Nbins);

		if (i == 32 ){
			printf("Sum --> %7.0f\n", sum);
			printf("Total --> %7.0f\n", integral[i]);
		}
		else {
			sum += integral[i];
			printf("%02i --> %7.0f \t (%7.0f)\n", i, integral[i], sum);
		}
	}
	
}
