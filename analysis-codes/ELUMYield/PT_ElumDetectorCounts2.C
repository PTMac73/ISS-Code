// PT_ElumDetectorCounts.C
// Fit the deuteron scattering peak
// ============================================================================================= //
// Patrick MacGregor
// Nuclear Physics Research Group
// School of Physics and Astronomy
// The University of Manchester
// ============================================================================================= //
#include "../Plotter/PTF_GetPosNumber.h"
#include "../WriteSPE.h"
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

// GLOBAL SWITCHES
Bool_t SWITCH_DRAW_CANVAS = 1;
Bool_t SWITCH_PRINT_CANVAS = 0;
Bool_t SWITCH_WRITE_SPECTRUM = 1;
TString spectrum_dir = "/home/ptmac/Documents/07-CERN-ISS-Mg/Mg-Analysis/ELUM";

// MAIN FUNCTION
void PT_ElumDetectorCounts2( TFile *f ){
	
	// Define the array-position number
	Int_t pos_num = GetPosNumber( f->GetName() );

	// Get the TTree
	TTree *t = (TTree*)f->Get("fin_tree");

	// Remove stats box
	gStyle->SetOptStat(0);

	// Check the number of counts
	//checkCounts(t);

	// Draw the full elum detector signal for each of the four quadrants
	TCanvas *c_elum = new TCanvas( "c_elum", "ELUM detector", 1200, 900 );
	c_elum->Divide(2,2);

	// Declare variables to be used in the loop
	const Int_t num_loops = 4;
	TH1F* h_elum[num_loops];
	Double_t d_lb = 300;
	Double_t d_ub = 750;
	TF1 *f_fitting[num_loops];
	TF1 *f_bg[num_loops];
	TF1 *f_gaussian[num_loops];
	Double_t par[6], *e_par;
	TLegend *leg[num_loops];
	
	// LOOP OVER THE FOUR QUADRANTS
	for (Int_t i = 0; i < num_loops; i++ ){
		// Change to the correct pad
		c_elum->cd(i+1);

		// Draw and store the histogram	
		t->Draw( Form( "elum[%i]>>h_elum_%i(2450, -800, 9000)", i + 4, i ) );
		h_elum[i] = (TH1F*)gDirectory->Get( Form( "h_elum_%i", i ) );

		// Now define a fit - Gaussian with a quadratic background
		f_fitting[i] = new TF1( Form("f_fitting_%i", i ), fitf, d_lb, d_ub, 6 );

		// Set the parameter limits
		/*   A */ f_fitting[i]->SetParLimits(0, 0, 1000);
		/*  mu */ f_fitting[i]->SetParLimits(1, 450, 600);
		/* sig */ f_fitting[i]->SetParLimits(2, 0, 500);

		// Help set the parameters to some sensible values
		f_fitting[i]->SetParameter(1, 520);
		f_fitting[i]->SetParameter(2, 25);

		// Style the fit
		f_fitting[i]->SetLineWidth(4);
		f_fitting[i]->SetLineColor(kRed);

		// Fit the function
		h_elum[i]->Fit( Form( "f_fitting_%i", i ), "R" );

		// Get the fit parameters
		f_fitting[i]->GetParameters(par);

		// Change viewing limits and style histogram
		h_elum[i]->GetXaxis()->SetRangeUser(-400,1600);
		h_elum[i]->GetYaxis()->SetRangeUser(-10,1.5*( par[0] + par[3]*par[1]*par[1] + par[4]*par[1] + par[5])  );
		h_elum[i]->GetYaxis()->SetTitle("#");
		h_elum[i]->GetXaxis()->SetTitle("Relative Energy");

		// Draw the background and Gaussian separately
		f_bg[i] = new TF1( Form( "f_bg_%i", i ), bg, d_lb, d_ub, 3);
		f_bg[i]->SetLineColor(kBlue);
		f_bg[i]->SetLineWidth(2);

		f_gaussian[i] = new TF1( Form( "f_gaussian_%i", i ), gaussian, d_lb, d_ub, 3);
		f_gaussian[i]->SetLineColor(kGreen+3);
		f_gaussian[i]->SetLineWidth(2);

		
		f_gaussian[i]->SetParameters(par);
		f_gaussian[i]->Draw("same");
		
		f_bg[i]->SetParameters(&par[3]);
		f_bg[i]->Draw("same");

		// Add a legend
		leg[i] = new TLegend(0.65,0.65,0.9,0.9);
		leg[i]->SetTextFont(42);
		leg[i]->SetTextSize(0.04);
		leg[i]->AddEntry(h_elum[i],"ELUM Signal","l");
		leg[i]->AddEntry(f_bg[i],"Background fit","l");
		leg[i]->AddEntry(f_gaussian[i],"Gaussian fit","l");
		leg[i]->AddEntry(f_fitting[i],"Total Fit","l");
		leg[i]->Draw();


		// Print window
		printf("===============================================================================================\n");
		printf(" CANVAS NUMBER %i -- GAUSSIAN FIT:\n", i);
		printf(" Amp. =\t%5.4e\t +/- \t%5.4e\n", par[0], f_fitting[i]->GetParError(0) );
		printf(" Mu   =\t%5.4e\t +/- \t%5.4e\n", par[1], f_fitting[i]->GetParError(1) );
		printf(" Sig. =\t%5.4e\t +/- \t%5.4e\n", par[2], f_fitting[i]->GetParError(2) );
		printf(" Integral =\t %5.4e\n\n", f_gaussian[i]->Integral(0,1000));

		printf(" BACKGROUND FIT (B2*x*x + B1*x + B0):\n");
		printf(" * B2 =\t%5.4e\t +/-\t%5.4e\n", par[3], f_fitting[i]->GetParError(3) );
		printf(" * B1 =\t%5.4e\t +/-\t%5.4e\n", par[4], f_fitting[i]->GetParError(4) );
		printf(" * B0 =\t%5.4e\t +/-\t%5.4e\n", par[5], f_fitting[i]->GetParError(5) );
		printf("===============================================================================================\n");

		// Write SPEs of the canvas
		if ( SWITCH_WRITE_SPECTRUM == 1 ){
			WriteSPE( h_elum[i]->GetName(), Form( "%s/elum_pos%i_%i", spectrum_dir.Data(), pos_num, i ) );
		}
	}
	// Print the canvas
	if ( SWITCH_PRINT_CANVAS == 1 ){
		c_elum->Print( Form( "/home/ptmac/Documents/07-CERN-ISS-Mg/analysis/PLOTS/ELUM_DD/elum_plot%i.pdf", pos_num ) );
	}

	
	
	return;
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
