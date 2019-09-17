// CreateMgSpectra.C
// Creates a number of spectra for the Mg analysis
// ============================================================================================= //
// Patrick MacGregor
// Nuclear Physics Research Group
// School of Physics and Astronomy
// The University of Manchester
// ============================================================================================= //
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TCutG.h>
#include <TH1F.h>
#include <iostream>
#include "GetRunNumber.h"
#include "WriteSPE.h"

// GLOBAL CONSTANTS
TString print_dir = "/home/ptmac/Documents/07-CERN-ISS-Mg/Mg-Analysis/SPE-Files";
Int_t CANVAS_SCALE = 300;
Int_t CANVAS_WIDTH = 4*CANVAS_SCALE;
Int_t CANVAS_HEIGHT = 3*CANVAS_SCALE;

TString make_name( Int_t type, Int_t spectrum_type, Int_t row_number ){
	// TYPE 0 --> histogram
	// TYPE 1 --> canvas
	Char_t type_char;
	if ( type == 0 ){ type_char = 'h'; }
	else if ( type == 1 ){ type_char = 'c'; }
	else { std::cout << "Unknown type declared." << std::endl; exit(1); }

	TString name = TString::Format( "%c_spectrum%i-%i", type_char, spectrum_type, row_number );
	return name;
}

TString make_print_name( Int_t pos_number, TString LABEL , Int_t row_number, TString file_type ){
	TString name = TString::Format( "pos%i_row%i_%s%s", pos_number, row_number, LABEL.Data(), file_type.Data() );
	return name;
}

Int_t Get_NUM_DIFF_TYPES_SPECTRA( Int_t mode ){
	Int_t a = 0;

	if ( mode == 1 ){
		a = 2;
	}
	else if ( mode == 2 ){
		a = 1;
	}
	return a;
}

Bool_t use_detector_1[24][2] = {
	{ 1, 1 },  // DET #0
	{ 1, 1 },  // DET #1
	{ 1, 1 },  // DET #2
	{ 1, 1 },  // DET #3
	{ 1, 1 },  // DET #4
	{ 1, 1 },  // DET #5
	{ 1, 1 },  // DET #6
	{ 1, 1 },  // DET #7
	{ 1, 1 },  // DET #8
	{ 1, 1 },  // DET #9
	{ 1, 1 },  // DET #10
	{ 1, 1 },  // DET #11
	{ 1, 0 },  // DET #12
	{ 1, 0 },  // DET #13
	{ 1, 0 },  // DET #14
	{ 1, 0 },  // DET #15
	{ 1, 0 },  // DET #16
	{ 1, 0 },  // DET #17
	{ 1, 1 },  // DET #18
	{ 1, 1 },  // DET #19
	{ 1, 1 },  // DET #20
	{ 1, 1 },  // DET #21
	{ 1, 1 },  // DET #22
	{ 1, 1 }   // DET #23
};

Bool_t use_detector_2[24][1] = {
	{ 1 },  // DET #0
	{ 1 },  // DET #1
	{ 0 },  // DET #2
	{ 1 },  // DET #3
	{ 0 },  // DET #4
	{ 0 },  // DET #5
	{ 0 },  // DET #6
	{ 0 },  // DET #7
	{ 1 },  // DET #8
	{ 0 },  // DET #9
	{ 1 },  // DET #10
	{ 0 },  // DET #11
	{ 0 },  // DET #12
	{ 0 },  // DET #13
	{ 0 },  // DET #14
	{ 0 },  // DET #15
	{ 0 },  // DET #16
	{ 0 },  // DET #17
	{ 0 },  // DET #18
	{ 0 },  // DET #19
	{ 0 },  // DET #20
	{ 0 },  // DET #21
	{ 0 },  // DET #22
	{ 1 }   // DET #23
};

TString label_1[2] = {
	"states1-6",
	"states7-9"
};

TString label_2[1] = { "best_res" };


// CHOOSE MODE OF OPERATION
// * 1 = All detectors
// * 2 = Best detectors
const Int_t MODE = 2;


void CreateMgSpectra( TFile *f, Int_t mode = MODE ){
	// Set into batch mode
	gROOT->SetBatch(kTRUE);

	// Change the style
	gStyle->SetOptStat("meni");

	// Globals
	const Int_t NUM_DIFF_TYPES_SPECTRA = Get_NUM_DIFF_TYPES_SPECTRA( mode );
	const Int_t NUM_DETECTORS = 24;
	const Int_t NUM_ROWS = 6;
	const Int_t NUM_DETECTORS_PER_ROW = 4;
	TString label[NUM_DIFF_TYPES_SPECTRA];
	Bool_t use_detector[NUM_DETECTORS][NUM_DIFF_TYPES_SPECTRA];
	
	if ( mode == 1 ){
		std::copy( &use_detector_1[0][0], &use_detector_1[0][0] + NUM_DETECTORS*NUM_DIFF_TYPES_SPECTRA, &use_detector[0][0] );
		std::copy( &label_1[0], &label_1[0] + NUM_DIFF_TYPES_SPECTRA, &label[0] );
	}
	else if ( mode == 2 ){
		std::copy( &use_detector_2[0][0], &use_detector_2[0][0] + NUM_DETECTORS*NUM_DIFF_TYPES_SPECTRA, &use_detector[0][0] );
		std::copy( &label_2[0], &label_2[0] + NUM_DIFF_TYPES_SPECTRA, &label[0] );
	}
	else{
		std::cout << "MODE " << mode << " NOT DEFINED. FORCE EXIT" << std::endl;
		exit(1);
	}
	
	// Open the TTree
	TTree *t = (TTree*)f->Get("fin_tree");

	// Get the cuts
	TCutG *cut0 = (TCutG*)gDirectory->Get("cut0");
	TCutG *cut1 = (TCutG*)gDirectory->Get("cut1");
	TCutG *cut2 = (TCutG*)gDirectory->Get("cut2");
	TCutG *cut3 = (TCutG*)gDirectory->Get("cut3");

	// Get the position from the file name
	Int_t pos_number = GetRunNumber( f->GetName() );
	
	// Define cut string components
	TString cut_string;
	TString cut_string_i = "( cut0 || cut1 || cut2 || cut3 ) && ( ";
	TString cut_string_f = " ) && thetaCM > 11 && td_rdt_e[] >= td_rdt_e_cuts[][0] && td_rdt_e[] <= td_rdt_e_cuts[][1] && xcal[] >= xcal_cuts[][0] && xcal[] <= xcal_cuts[][1]";

	// Define canvas and histogram things
	TCanvas *c_spectrum[NUM_DIFF_TYPES_SPECTRA][NUM_ROWS];
	TH1F *h_spectrum[NUM_DIFF_TYPES_SPECTRA][NUM_ROWS];
	
	// Loop over rows, detectors per row, and different types of spectra
	for ( Int_t k = 0; k < NUM_DIFF_TYPES_SPECTRA; k++ ){
		
		// Define an output file to save each histogram
		TFile *f_out = new TFile( Form( "%s/pos%i_%s_hists.root", print_dir.Data(), pos_number, label[k].Data() ), "RECREATE" );

		for ( Int_t i = 0; i < NUM_ROWS; i++ ){

			// Make the cut string for each row
			cut_string = "";
			cut_string += cut_string_i;

			// See if spectrum will be filled
			Bool_t using_spectra = 0;

			for ( Int_t j = 0; j < NUM_DETECTORS_PER_ROW; j++ ){
				// Produce spectrum if correct
				if ( use_detector[i + NUM_ROWS*j][k] == 1 ){
					cut_string += Form( "detID == %i", i + NUM_ROWS*j );
					using_spectra = 1;
				}
				if ( use_detector[i + NUM_ROWS*(j+1)][k] == 1 && j != NUM_DETECTORS_PER_ROW - 1 && using_spectra == 1 ){
					cut_string += " || ";
				}
			}

			if ( using_spectra == 0 ){
				std::cout << "Row " << i << "not used." << std::endl;
			}
			else{
				cut_string += cut_string_f;
				std::cout << "CUT: " << cut_string << std::endl;
				
				// Make the spectra
				c_spectrum[k][i] = new TCanvas( make_name( 1, k, i ).Data(), "CANVAS", CANVAS_WIDTH, CANVAS_HEIGHT );
				t->Draw( Form( "Ex>>%s(450, -1, 8)", make_name( 0, k, i ).Data() ), cut_string );

				// Get the histogram
				h_spectrum[k][i] = (TH1F*)gDirectory->Get( make_name( 0, k, i ) );

				// Format the histogram
				h_spectrum[k][i]->SetTitle( Form( "Row %i for %s; Ex / MeV; #", i, label[k].Data() ) );

				// Print the canvas and write the spectrum
				c_spectrum[k][i]->Print( Form( "%s/%s", print_dir.Data(), make_print_name( pos_number, label[k], i, ".pdf").Data() ), "EmbedFonts" );
				WriteSPE( h_spectrum[k][i]->GetName(), Form( "%s/%s", print_dir.Data(), make_print_name(pos_number, label[k], i, "").Data() ) );
				std::cout << ">>> Completed position " << pos_number << ", row " << i << " for " << label[k] << std::endl << std::endl;
				
				// Write the histogram to file
				h_spectrum[k][i]->Write();
			}
		}
		
		// Close the output file
		f_out->Close();

	}
	
	return;
}


































