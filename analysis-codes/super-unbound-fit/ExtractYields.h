// ExtractYields.h
// Contains peak information and useful data for ExtractYields.C
// ============================================================================================= //
// Patrick MacGregor
// Nuclear Physics Research Group
// Department of Physics and Astronomy
// The University of Manchester
// ============================================================================================= //

#ifndef EXTRACT_YIELDS_H_
#define EXTRACT_YIELDS_H_

#include <TFile.h>
#include <TMath.h>
#include <TPad.h>
#include <TString.h>

#include <iostream>

// CONSTANTS

// Location of ROOT files
TString root_file_dir = "pos%i_ex.root";
TString print_dir = "/home/ptmac/Documents/07-CERN-ISS-Mg/analysis/analysis-codes/super-unbound-fit";

// Things that probably won't change
const Int_t BG_DIM = 0; 		// (0 = const, 1 = linear etc.)
const Int_t NUM_PEAKS = 1;		// Number of peaks in the spectrum
const Double_t E_MIN = 5.38;	// Minimum energy on the spectrum
const Double_t E_MAX = 5.8;//6.17;	// Maximum energy on the spectrum
const Double_t PRE_MIN = 0.4;	// Minimum energy for the two strong states used to get widths
const Double_t PRE_MAX = 1.9;	// Maximum energy for the two strong states used to get widths

Double_t PEAK_WIDTH_ESTIMATE = 0.2;		// Range within which peak definitely falls
Double_t MAX_SF_ESTIMATE = 0.5;

const Double_t NEUTRON_SEP_ENERGY = 3.655;
const Double_t GS_DOUB_LB = -PEAK_WIDTH_ESTIMATE;	// Bounds for the ground state doublet yield evaluation
const Double_t GS_DOUB_UB = PEAK_WIDTH_ESTIMATE;

const Int_t NUM_SPECTRUM_SECTIONS = 2;

Int_t CANVAS_WIDTH = 1200;
Int_t CANVAS_HEIGHT = 900;

// SWITCHES
const Bool_t SW_BATCH_MODE = 0;
const Bool_t SW_DRAW_PRE_FIT = 1;
const Bool_t SW_HIST_TITLE = 0;


// Define information about the peaks
const Double_t pre_peak_energies[2] = {
	1.10,
	1.40
};

const Double_t peak_energies[NUM_PEAKS] = {
	5.55,
	//5.76,
	//6.07
};


// FUNCTIONS
// Generate title for histograms
TString GenerateHistTitle( Int_t pos, Int_t spectrum_type ){
	if ( spectrum_type == 2 ){
		return Form( "Position %i | All Rows", pos );
	}
	else{
		return Form( "Position %i | Row -- | States %s", pos, ( spectrum_type == 0 ? "1-6" : "7-9" ) );
	}
}

// Generates a name for a canvas
TString GenerateCanvasName( Int_t a ){
	return Form( "c_spec_pos%i", a );
}


// GENERATE NAMES FOR FILES
enum FileType_t { svg, dat, tex, pdf };

TString GenerateSuffix( FileType_t ft ){
	switch(ft){
		case svg  : return "svg";  break;
		case dat  : return "dat";  break;
		case tex  : return "tex";  break;
		case pdf  : return "pdf";  break;
		default   : return "confused"; break;
	}
}

TString GenerateFileName( FileType_t ft, Int_t pos, Int_t c_type = 0 ){
	// Define prefix if necessary
	TString prefix = "";
	if ( c_type == 1 ){ prefix = "pre_"; }
	
	// Define string relative to root file directory
	TString a = root_file_dir( 0, root_file_dir.First('.') );
	return Form( "%s/%s%s.%s",
		print_dir.Data(),
		prefix.Data(),
		Form( a, pos ), 
		GenerateSuffix(ft).Data()
	);
}

// Generates ROOT file name
TString GetROOTFileName( Int_t pos ){
	return Form( root_file_dir, pos  );
}


// Gets the number of histograms in a ROOT file
Int_t GetNumHists( TFile* f ){
	// Assumes file contains only histograms
	TList* l = (TList*)f->GetListOfKeys();
	Int_t num_hists = l->GetEntries();
	return num_hists;
}

// Gets the names of all the histograms
void GetHistNames( TList* l, TString* str_arr ){
	Int_t num_hists = l->GetEntries();
	for ( Int_t i = 0; i < num_hists; i++ ){
		str_arr[i] = (TString)l->At(i)->GetName();
	}
	return;
}

// Sets the margins for a canvas
void SetPadMargins( TPad* pad ){
	pad->SetLeftMargin(0.12);
	pad->SetRightMargin(0.01);
	if ( SW_HIST_TITLE == 1 ){ pad->SetTopMargin(0.08); }
	else{ pad->SetTopMargin(0.01); }
	pad->SetBottomMargin(0.12);
	return;
}


#endif
