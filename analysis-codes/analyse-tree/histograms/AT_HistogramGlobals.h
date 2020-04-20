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
	h->GetYaxis()->SetTitleFont(62);
	h->GetYaxis()->SetLabelFont(62);
	h->SetTitleFont(62);
	return;
}

void GlobSetHistFonts( TH2* h ){
	h->GetXaxis()->SetTitleFont(62);
	h->GetXaxis()->SetLabelFont(62);
	h->GetYaxis()->SetTitleFont(62);
	h->GetYaxis()->SetLabelFont(62);
	h->SetTitleFont(62);
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
	h = new TH1F( name.Data(), "", 450, -1, 8 );
	h->SetTitle("");
	h->GetXaxis()->SetTitle("Excitation Energy (MeV)");
	h->GetYaxis()->SetTitle("Counts per 20 keV");	
	h->SetLineColor(kRed);	
	GlobSetHistFonts( h );
	return;
}

void CreateEVZSpectrum(){
	/*TODO*/
	return;
}





#endif
