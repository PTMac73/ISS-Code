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
#include <TPad.h>
#include <TString.h>

#include <iostream>

// CONSTANTS
const Int_t BG_DIM = 0; 		// (0 = const, 1 = linear etc.)
const Int_t NUM_PEAKS = 3;		// Number of peaks in the spectrum
const Double_t E_MIN = 3.6;		// Minimum energy on the spectrum
const Double_t E_MAX = 4.4;	// Maximum energy on the spectrum
const Double_t PRE_MIN = 0.4;	// Minimum energy for the two strong states used to get widths
const Double_t PRE_MAX = 1.9;	// Maximum energy for the two strong states used to get widths

Double_t PEAK_WIDTH_ESTIMATE = 0.1;		// Range within which peak definitely falls

const Double_t NEUTRON_SEP_ENERGY = 3.655;

const Double_t GS_DOUB_LB = -PEAK_WIDTH_ESTIMATE;	// Bounds for the ground state doublet yield evaluation
const Double_t GS_DOUB_UB = PEAK_WIDTH_ESTIMATE;

//const Int_t NUM_ROWS = 6;
const Int_t NUM_POSITIONS = 2;
const Int_t NUM_SPECTRUM_SECTIONS = 2;

Int_t CANVAS_WIDTH = 1200;
Int_t CANVAS_HEIGHT = 900;

// SWITCHES
const Bool_t SW_BATCH_MODE = 1;
const Bool_t SW_DRAW_PRE_FIT = 1;
const Bool_t SW_HIST_TITLE = 0;



// Location of ROOT files
TString root_file_dir = "pos%i_ex_best.root";
TString print_dir = "/home/ptmac/Documents/07-CERN-ISS-Mg/analysis/analysis-codes/unbound-doublet-fit";

Int_t num_peaks_in_spectrum[NUM_POSITIONS] = {
	NUM_PEAKS,
	NUM_PEAKS
};

// Define information about the peaks
const Double_t pre_peak_energies[2] = {
	1.10,
	1.40
};

const Double_t peak_energies[NUM_PEAKS] = {
	3.784951973,
	3.967225798,
	4.307444788
};

const Bool_t peak_fix_positions[NUM_PEAKS] = {
	1,	// 3.86 ---> UNBOUND
	1,	// 3.97 ---> UNBOUND
	1
};


// FUNCTIONS

TString GetROOTFileName( Int_t pos ){
	return Form( root_file_dir, pos  );
}


Int_t GetNumHists( TFile* f ){
	// Assumes file contains only histograms
	TList* l = (TList*)f->GetListOfKeys();
	Int_t num_hists = l->GetEntries();
	return num_hists;
}

void GetHistNames( TList* l, TString* str_arr ){
	Int_t num_hists = l->GetEntries();
	for ( Int_t i = 0; i < num_hists; i++ ){
		str_arr[i] = (TString)l->At(i)->GetName();
	}
	return;
}

/*TString GetHistName( Int_t row_num ){
	return Form( "h_spectrum0_%i_0", row_num );
}*/

TString GenerateHistTitle( Int_t pos, Int_t spectrum_type ){
	if ( spectrum_type == 2 ){
		return Form( "Position %i | All Rows", pos );
	}
	else{
		return Form( "Position %i | Row -- | States %s", pos, ( spectrum_type == 0 ? "1-6" : "7-9" ) );
	}
}

TString GenerateCanvasName( Int_t a ){
	return Form( "c_spec_pos%i", a );
}

TString GenerateCanvasPrintName(  Int_t pos, Int_t spectrum_type, Int_t c_type = 0 ){
	TString prefix = "";
	if ( c_type == 1 ){ prefix = "pre_"; }

	if ( spectrum_type == 2 ){
		TString a = root_file_dir( 0, root_file_dir.First('.') );
		return Form( "%s/%s.svg", print_dir.Data(), Form( a, pos ) );
	}
	else{
		return "FAIL.svg";
	}
}


TString GetDataFileName( Int_t pos, Int_t spectrum_type ){
	if ( spectrum_type == 2 ){
		TString a = root_file_dir( 0, root_file_dir.First('.') );
		return Form( "%s/%s.dat", print_dir.Data(), Form( a, pos ) );
	}
	else{
		return "FAIL.dat";
	}
}

void SetPadMargins( TPad* pad ){
	pad->SetLeftMargin(0.12);
	pad->SetRightMargin(0.01);
	if ( SW_HIST_TITLE == 1 ){ pad->SetTopMargin(0.08); }
	else{ pad->SetTopMargin(0.01); }
	pad->SetBottomMargin(0.12);
	return;
}

#endif
