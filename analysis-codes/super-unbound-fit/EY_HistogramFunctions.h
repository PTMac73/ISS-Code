// EY_HistogramFunctions.h
// Histogram functions for extracting yields
// ============================================================================================= //
// Patrick MacGregor
// Nuclear Physics Research Group
// Department of Physics and Astronomy
// The University of Manchester
// ============================================================================================= //

#ifndef EY_HISTOGRAM_FUNCTIONS_H_
#define EY_HISTOGRAM_FUNCTIONS_H_

#include <TFitResult.h>
#include <TH1.h>
#include <TMath.h>
#include <TStyle.h>

// --------------------------------------------------------------------------------------------- //
Double_t GetAmpEstimate( TH1D* h, Double_t lb, Double_t ub ){
	Int_t lb_bin = h->FindBin(lb);
	Int_t ub_bin = h->FindBin(ub);
	
	Double_t ae = 0.0;
	
	for ( Int_t i = lb_bin; i <= ub_bin; i++ ){
		if ( h->GetBinContent(i) > ae ){
			ae = h->GetBinContent(i);
		}
	}
	return ae;
}


// --------------------------------------------------------------------------------------------- //
void FormatHistogram( TH1D* h ){
	h->GetXaxis()->SetTitle( "Excitation Energy (MeV)" );
	h->GetXaxis()->CenterTitle();
	h->GetXaxis()->SetRangeUser(-0.5,7);
	
	h->GetYaxis()->SetTitle( "Counts per 20 keV" );
	h->GetYaxis()->CenterTitle();
	
	Double_t font_size = h->GetYaxis()->GetLabelSize();
	h->GetXaxis()->SetLabelSize(1.2*font_size);
	h->GetYaxis()->SetLabelSize(1.2*font_size);
	h->GetXaxis()->SetTitleSize(1.5*font_size);
	h->GetYaxis()->SetTitleSize(1.5*font_size);
	
	h->SetTitleFont(62, "xy");
	h->SetLabelFont(62, "xy");
	h->SetLineColor( kRed );
	gStyle->SetTitleFont(62,"t");

	return;
}


#endif
