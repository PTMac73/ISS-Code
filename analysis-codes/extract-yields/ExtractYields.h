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
TString root_file_dir = "pos%i_ex_corr.root";
TString print_dir = "/home/ptmac/Documents/07-CERN-ISS-Mg/analysis/analysis-codes/extract-yields/print";

// Things that probably won't change
const Int_t BG_DIM = 0; 		// (0 = const, 1 = linear etc.)

// REGION-SPECIFIC QUANTITIES
const Int_t NUM_PEAKS = 13;

// Min/max energy for the region
Double_t E_LIMITS[2] = {
	0.0,
	6.3
};

// Peak energies
/*const Double_t peak_energies[NUM_PEAKS] = {
	0.00000,
	1.09462,
	1.43070,
	2.26590,
	2.49990,
	2.87000,
	3.18000,

	3.854334927,
	3.989180979,
	4.309017061,

	5.570797291,
	5.759762758,
	5.989528736
};*/
const Double_t peak_energies[NUM_PEAKS] = {
	0.00000,
	1.09462,
	1.43070,
	2.26590,
	2.49990,
	2.87000,
	3.18000,

	3.90,
	4.04,
	4.4,

	5.629749246,
	5.825187047,
	6.056059684
};
/*const Double_t peak_energies[NUM_PEAKS] = {
	0.00000,
	1.09462,
	1.43070,
	2.26590,
	2.49990,
	2.87000,
	3.18000,
	3.874186853,
	4.011653487,
	4.313766313,
	5.6,
	5.72,
	6.0
	//5.671,
	//5.913,
	//6.101
};*/

Int_t peak_colours[NUM_PEAKS] = {
	kOrange,
	kRed,
	kRed,
	kRed,
	kRed,
	kRed,
	kRed,

	kBlue,
	kBlue,
	kBlue,

	kCyan,
	kCyan,
	kCyan
};

Double_t peak_label_pos_offset[NUM_PEAKS][2] = {
	{ 0.03, 25.0 },
	{ 0.00, 10.0 },
	{ 0.2, 10.0 },
	{ -0.05, 24.0 },
	{ 0.00, 26.0 },
	{ 0.00, 22.0 },
	{ 0.00, 25.0 },

	{ 0.00, 18.0 },
	{ 0.00, 20.0 },
	{ 0.00, 20.0 },

	{ 0.00, 25.0 },
	{ 0.00, 25.0 },
	{ 0.00, 15.0 }
};

TString peak_label_gs_str[5] = {
	"g.s.",
	"l=2",
	"+",
	"55",
	"l=0"
};

TString peak_label_ubd_str[4] = {
	"3929",
	"+",
	"4068",
	"l=1",
};

TString peak_label[NUM_PEAKS][2] = {
	{"",""},					// 00***
	{ "1092", "l=1",  },	// 01
	{ "1432", "l=3",  },	// 02
	{ "2270", "l=1",  },	// 03
	{ "2501", "l=2",  },	// 04
	{ "2900", "l=3",  },	// 05
	{ "3220", "l=2",  },	// 06
	{"",""},				// 07***
	{"",""},				// 08***
	{ "4360", "l=1",  },	// 09
	{ "", "5627",  },	// 10
	{ "", "5817",  },	// 11
	{ "", "6049",  },	// 12
};


// REGION BOOLEANS
Bool_t fix_widths[NUM_PEAKS] = {
	0,1,1,1,1,1,1,
	0,0,0,
	1,1,1
};
Bool_t fix_positions[NUM_PEAKS] = {
	0,0,0,0,0,0,0,
	0,0,0,
	0,0,0
};

// UNBOUND DOUBLET PARAMETERS
const Int_t FIRST_PEAK_J = 3;	// 1 or 3 for unbound doublet
const Double_t MAX_SF_WIDTH = 0.210;	// 210 keV max width due to sf
const Double_t sf_estimates[3] = {		// Estimates for current strength already in other peaks
	// 3/2- --> 0.364, 1/2- --> 0.809 [OLD]
	// 3/2- --> 0.650, 1/2- --> 0.827
	( FIRST_PEAK_J == 1 ? 0.827 : 0.650 ),	// SF of the first peak
	( FIRST_PEAK_J != 1 ? 0.827 : 0.650 ),	// SF of the second peak
	1.00	// SF of the third peak
};


// OTHER VARIABLES
const Double_t PRE_MIN = 0.4;	// Minimum energy for the two strong states used to get widths
const Double_t PRE_MAX = 1.9;	// Maximum energy for the two strong states used to get widths
Double_t PEAK_WIDTH_ESTIMATE = 0.1;		// Range within which peak definitely falls
Bool_t PRE_FIX_WID[2] = {1,1};
Bool_t PRE_FIX_POS[2] = {0,0};
const Double_t NEUTRON_SEP_ENERGY = 3.655;


//const Int_t NUM_ROWS = 6;
//const Int_t NUM_POSITIONS = 2;
//const Int_t NUM_SPECTRUM_SECTIONS = 2;

Int_t CANVAS_WIDTH = 1200;
Int_t CANVAS_HEIGHT = 900;

// SWITCHES
const Bool_t SW_BATCH_MODE = 1;
const Bool_t SW_DRAW_PRE_FIT = 1;
const Bool_t SW_HIST_TITLE = 0;


// Define information about the peaks
const Double_t pre_peak_energies[2] = {
	1.10,
	1.40
};



// FUNCTIONS
// Generate title for histograms
TString GenerateHistTitle( Int_t pos ){
	return Form( "Position %i | All Rows", pos );
}

// Generates a name for a canvas
TString GenerateCanvasName( Int_t a ){
	return Form( "c_spec_pos%i", a );
}

// Generate a string to denote the ordering of states
TString GenerateJString(){
	return Form( "%i%i", FIRST_PEAK_J, 4 - FIRST_PEAK_J );
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
	return Form( "%s/%s_%s%s.%s",
		print_dir.Data(),
		Form( a, pos ),
		prefix.Data(),
		GenerateJString().Data(),
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

// Calculates upper boundary for peak width
Double_t CalculatePeakUpperWidth( Int_t peak_num, Double_t w ){
	return TMath::Sqrt( w*w + TMath::Power( sf_estimates[peak_num]*MAX_SF_WIDTH, 2 ) );
}

#endif
