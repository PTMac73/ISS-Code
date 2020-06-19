// AT_HistogramGlobals.h
// Global functions for AT_Histograms.h
// ============================================================================================= //
// Patrick MacGregor
// Nuclear Physics Research Group
// Department of Physics and Astronomy
// The University of Manchester
// ============================================================================================= //
#ifndef AT_HISTOGRAM_GLOBALS_H_
#define AT_HISTOGRAM_GLOBALS_H_

#include <TCanvas.h>
#include <TH1.h>
#include <TH2.h>
#include <TPaveText.h>
#include <TProfile.h>
#include <TString.h>
#include <iostream>


// --------------------------------------------------------------------------------------------- //
// GLOBAL FUNCTIONS
// Set canvas margins
void GlobSetCanvasMargins( TCanvas *c, Double_t l = 0.1, Double_t r = 0.02, Double_t t = 0.02, Double_t b = 0.1 ){
	TPad* pad = (TPad*)c;
	pad->SetLeftMargin( l );
	pad->SetRightMargin( r );
	pad->SetTopMargin( t );
	pad->SetBottomMargin( b );
	return;
}


// Set histogram fonts
void GlobSetHistFonts( TH1* h ){
	h->GetXaxis()->SetTitleFont(62);
	h->GetXaxis()->SetLabelFont(62);
	h->GetXaxis()->CenterTitle();
	h->GetYaxis()->SetTitleFont(62);
	h->GetYaxis()->SetLabelFont(62);
	h->GetYaxis()->CenterTitle();
	h->SetTitleFont(62);
	return;
}

void GlobSetHistFonts( TH2* h ){
	h->GetXaxis()->SetTitleFont(62);
	h->GetXaxis()->SetLabelFont(62);
	h->GetXaxis()->CenterTitle();
	h->GetYaxis()->SetTitleFont(62);
	h->GetYaxis()->SetLabelFont(62);
	h->GetYaxis()->CenterTitle();
	h->SetTitleFont(62);
	return;
}

void GlobSetHistFonts( TGraph* g ){
	g->GetXaxis()->SetTitleFont(62);
	g->GetXaxis()->SetLabelFont(62);
	g->GetXaxis()->CenterTitle();
	g->GetYaxis()->SetTitleFont(62);
	g->GetYaxis()->SetLabelFont(62);
	g->GetYaxis()->CenterTitle();
	return;
}


// Write error message about spe files
void ErrorSPE( TString str = "This mode" ){
	std::cout << "*** ERROR: " << str << " cannot have a valid .spe file" << "\n";
	return;
}


// Convert numbers into suitable strings for printing
TString DoubleToString( Double_t a ){
	TString b = Form( "%3.2f", a );
	if ( b.Contains('.') ){
		b.ReplaceAll( '.', '_' );
	}
	return b;
}


// Print functions
void PrintPDF( TCanvas* c, TString spec_name ){
	c->Print( ( spec_name + ".svg" ).Data() );
	return;
}

void PrintPNG( TCanvas* c, TString spec_name ){
	c->Print( ( spec_name + ".png" ).Data() );
	return;
}


void PrintAll( TCanvas* c, TString spec_name ){
	if ( PRINT_PDF ){ PrintPDF( c, spec_name ); }
	if ( PRINT_PNG ){ PrintPNG( c, spec_name ); }
	return;
}

// Set canvas title font
void SetCanvasTitleFont( TPad* pad ){
	pad->GetListOfPrimitives()->Print();
	std::cout << "===" << "\n";
	TPaveText* t = (TPaveText*)pad->GetPrimitive("title");
	if( t!= NULL ){ t->Print(); }
	/*TODO - Get me working!!!*/
	//t->SetTextFont();
	return;
}


// Create cherished spectra
void CreateExSpectrum( TH1F*& h, TString name ){
	h = new TH1F( name.Data(), name.Data(), 450, -1, 8 );
	h->SetTitle("");
	h->GetXaxis()->SetTitle("Excitation Energy (MeV)");
	h->GetYaxis()->SetTitle("Counts per 20 keV");	
	h->SetLineColor(kRed);	
	GlobSetHistFonts( h );
	return;
}

// Deifne lower and upper bounds
void CreateExSpectrum( TH1F*& h, TString name, Double_t lb, Double_t ub ){
	Int_t nbins = (Int_t)( 50*( ub - lb ) );
	h = new TH1F( name.Data(), name.Data(), nbins, lb, ub );
	h->SetTitle("");
	h->GetXaxis()->SetTitle("Excitation Energy (MeV)");
	h->GetYaxis()->SetTitle("Counts per 20 keV");	
	h->SetLineColor(kRed);	
	GlobSetHistFonts( h );
	return;
}



// Format generic 2D histograms
void GlobCreate2DHists( TH2F* h, TString x_label, TString y_label ){
	h->SetTitle("");
	h->GetXaxis()->SetTitle( x_label.Data() );
	h->GetYaxis()->SetTitle( y_label.Data() );
	h->SetMarkerStyle(20);
	h->SetMarkerSize(0.5);
	h->SetMarkerColor(kRed);	
	GlobSetHistFonts( h );
	return;
}

void CreateEVZSpectrum( TH2F*& h, TString name ){
	h = new TH2F( name.Data(), name.Data() , 400, -50, -10, 900, 0, 9 );
	GlobCreate2DHists( h, "z (cm)", "Energy (MeV)" );
	h->SetMarkerSize(0.2);
	return;
}



// Format generic profile plots
void GlobCreateProfile( TProfile* p, TString x_label, TString y_label ){
	p->SetTitle("");
	p->GetXaxis()->SetTitle( x_label.Data() );
	p->GetYaxis()->SetTitle( y_label.Data() );
	p->SetMarkerColor(kRed);
	p->SetLineColor(kRed);
	GlobSetHistFonts( p );
	return;
}


Double_t GetMeanBinPosition( TH1F* h, Double_t lb, Double_t ub, Double_t &err ){
	Double_t sum1 = 0;	// Weighted mean numerator
	Double_t sum2 = 0;	// Weighted mean denominator
	
	Double_t sumXXN = 0;	// X = bin position, N = number of counts
	Double_t sumXN = 0;
	Double_t sumN = 0;
	
	Int_t bub = h->FindBin(ub);
	Int_t blb = h->FindBin(lb);
	
	Int_t num_bins = bub - blb;
	Double_t w = h->GetBinWidth( blb );
	
	for ( Int_t i = 0; i < num_bins; i++ ){
		sum1 += h->GetBinContent( i + blb )*h->GetBinCenter( i + blb );
		sum2 += h->GetBinContent( i + blb );
		
		sumN += h->GetBinContent( i + blb );
		sumXN += h->GetBinContent( i + blb )*h->GetBinCenter( i + blb );
		sumXXN += h->GetBinContent( i + blb )*h->GetBinCenter( i + blb )*h->GetBinCenter( i + blb );
	}
	
	err = sumXXN/( sumN*sumN ) - sumXN*sumXN/( sumN*sumN*sumN ) + w*w/4.0;
	
	std::cout << sum1 << "\t" << sum2 << "\t" << sumN << "\t" << sumXN << "\t" << sumXXN << "\t" << w << "\t" << num_bins << "\t" << err << "\n";
	
	return sum1/sum2;
}



TString SideString( Int_t det ){
	if ( det >=  0 && det <=  5 ){ return "left"; }	
	else if ( det >=  6 && det <= 11 ){ return "bottom"; }
	else if ( det >= 12 && det <= 17 ){ return "right"; }
	else if ( det >= 18 && det <= 23 ){ return "top"; }
	else{ return ""; }
}






































#endif
