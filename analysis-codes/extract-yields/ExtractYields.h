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
const Int_t NUM_PEAKS = 10;		// Number of peaks in the spectrum
const Double_t E_MIN = -1.0;	// Minimum energy on the spectrum
const Double_t E_MAX = 4.6;		// Maximum energy on the spectrum
const Double_t PRE_MIN = 0.4;	// Minimum energy for the two strong states used to get widths
const Double_t PRE_MAX = 1.9;	// Maximum energy for the two strong states used to get widths

Double_t PEAK_WIDTH_ESTIMATE = 0.3;		// Range within which peak definitely falls

const Double_t NEUTRON_SEP_ENERGY = 3.655;

const Double_t GS_DOUB_LB = -PEAK_WIDTH_ESTIMATE;	// Bounds for the ground state doublet yield evaluation
const Double_t GS_DOUB_UB = PEAK_WIDTH_ESTIMATE;

const Int_t NUM_ROWS = 6;
const Int_t NUM_POSITIONS = 2;
const Int_t NUM_SPECTRUM_SECTIONS = 2;

Int_t CANVAS_WIDTH = 1200;
Int_t CANVAS_HEIGHT = 900;

// SWITCHES
const Bool_t SW_BATCH_MODE = 1;
const Bool_t SW_DRAW_PRE_FIT = 1;
const Bool_t SW_HIST_TITLE = 0;



// Location of ROOT files
TString root_file_dir = "/home/ptmac/Documents/07-CERN-ISS-Mg/Mg-Analysis/SPE-Files/spe-files-004/root-files";
TString print_dir = "/home/ptmac/Documents/07-CERN-ISS-Mg/Mg-Analysis/SPE-Files/spe-files-004/yields";

Int_t num_peaks_in_spectrum[NUM_ROWS][NUM_POSITIONS] = {
	{1,3},
	{3,5},
	{5,NUM_PEAKS},
	{NUM_PEAKS,NUM_PEAKS},
	{NUM_PEAKS,NUM_PEAKS},
	{NUM_PEAKS,NUM_PEAKS}
};

// Define information about the peaks
const Double_t peak_energies[NUM_PEAKS] = {
	0.00,
	1.09462,
	1.4307,
	2.2659,
	2.4999,
	2.87,
	3.18,
	3.86,
	3.97,
	4.28,
};

const Bool_t peak_fix_widths[NUM_PEAKS] = {
	0,	// 0.0000 ---> DOUBLET
	1,	// 1.09
	1,	// 1.43
	1,	// 2.27
	1,	// 2.50
	1,	// 2.87
	1,	// 3.18
	0,	// 3.86 ---> UNBOUND
	0,	// 3.97 ---> UNBOUND
	0	// 4.28 ---> UNBOUND
};

const Bool_t peak_fix_positions[NUM_PEAKS] = {
	0,	// 0.0000 ---> DOUBLET 1/2
	0,	// 1.09
	0,	// 1.43
	0,	// 2.27
	0,	// 2.50
	0,	// 2.87
	0,	// 3.18
	0,	// 3.86 ---> UNBOUND
	0,	// 3.97 ---> UNBOUND
	0	// 4.28 ---> UNBOUND
};


// FUNCTIONS
/* spectrum_type = 0 --> states 1-6
                   1 --> states 7-9
                   2 --> full excitation spectrum
*/
TString CUTS_MODE_ID[2] = {
	"all",
	"sing"
};


TString GetROOTFileName( Int_t pos, Int_t spectrum_type ){
	if ( spectrum_type == 2 ){
		return Form( "%s/pos%i_all_full.root", root_file_dir.Data(), pos );
	}
	else{
		return Form( "%s/pos%i_all_row.root", root_file_dir.Data(), pos  );
	}
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

TString GenerateHistTitle( Int_t pos, Int_t row, Int_t spectrum_type ){
	if ( spectrum_type == 2 ){
		return Form( "Position %i | All Rows", pos );
	}
	else{
		return Form( "Position %i | Row %i | States %s", pos, row, ( spectrum_type == 0 ? "1-6" : "7-9" ) );
	}
}

TString GenerateCanvasName( Int_t row_num ){
	return Form( "c_spec_%i", row_num );
}

TString GenerateCanvasPrintName(  Int_t pos, Int_t row, Int_t spectrum_type, Int_t c_type = 0 ){
	TString prefix = "";
	if ( c_type == 1 ){ prefix = "pre_"; }

	if ( spectrum_type == 2 ){
		return Form( "%s/pos%i_all_rows_%sfit.svg", print_dir.Data(), pos, prefix.Data() );
	}
	else{
		return Form( "%s/pos%i_row%i_states%s_%sfit.svg", print_dir.Data(), pos, row, ( spectrum_type == 0 ? "1-6" : "7-9" ), prefix.Data() );
	}

}


TString GetDataFileName( Int_t pos, Int_t spectrum_type, Int_t row_num ){
	if ( spectrum_type == 2 ){
		return Form( "%s/pos%i_all_rows.dat", print_dir.Data(), pos );
	}
	else{
		return Form( "%s/pos%i_row%i_states%s.dat", print_dir.Data(), pos, row_num, ( spectrum_type == 0 ? "1-6" : "7-9" ) );
	}
}

void SetPadMargins( TPad* pad ){
	pad->SetLeftMargin(0.09);
	pad->SetRightMargin(0.01);
	if ( SW_HIST_TITLE == 1 ){ pad->SetTopMargin(0.08); }
	else{ pad->SetTopMargin(0.01); }
	pad->SetBottomMargin(0.09);
	return;
}

#endif
