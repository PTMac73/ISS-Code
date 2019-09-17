// CreateStripSpectra.C
// Reads in one side of the array and creates a number of strip spectra for it
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
#include <iostream>
#include "GetRunNumber.h"
#include "WriteSPE.h"

// GLOBAL CONSTANTS
Double_t xcal_min = 0.11;
Double_t xcal_max = 0.84;
const Int_t num_strips_per_Si = 1;
const Int_t num_detectors_per_row = 6;
const Int_t row_number = 4;		// If this = 4, it does all of the rows
TString print_dir = "/home/ptmac/Documents/07-CERN-ISS-Mg/Mg-Analysis/DIAGNOSTICS";

Int_t CANVAS_SCALE = 300;
Int_t CANVAS_WIDTH = 4*CANVAS_SCALE;
Int_t CANVAS_HEIGHT = 3*CANVAS_SCALE;

/* 
>> Define a type for the reaction
     * 0 --> 28Mg(d,p)
     * 1 --> Alpha run
*/
const Bool_t TYPE = 0;

void CreateStripSpectra( TFile *f ){
	// Set into batch mode
	gROOT->SetBatch(kTRUE);

	// Change the style
	gStyle->SetOptStat("meni");

	// Define options (don't include run numbers)
	TString key_word;
	TString branch;
	TString bins;
	TString cut_string;
	if ( TYPE == 0 ){
		key_word = "position";
		branch = "Ex";
		bins = "(450,-1,8)";
		TCutG *cut0 = (TCutG*)gDirectory->Get("cut0");
		TCutG *cut1 = (TCutG*)gDirectory->Get("cut1");
		TCutG *cut2 = (TCutG*)gDirectory->Get("cut2");
		TCutG *cut3 = (TCutG*)gDirectory->Get("cut3");
		cut_string = Form("( cut0 || cut1 || cut2 || cut3 ) && thetaCM > 11 && td_rdt_e[] > td_rdt_e_cuts[][0] && td_rdt_e[] < td_rdt_e_cuts[][1] && ");
	}
	else if ( TYPE == 1 ){
		key_word = "alpha";
		branch = "e";
		bins = "(1000,0,2000)";
		cut_string = "";
	}

	// Check how many rows to loop over - define variables
	Int_t index0, index1;
	

	// Loop over all rows
	if ( row_number == 4 ){
		index0 = 0;
		index1 = row_number;
	}
	// Loop over 1 row
	else if ( row_number < 4  && row_number >= 0){
		index0 = row_number;
		index1 = row_number+1;
	}
	// Error somewhere
	else{
		std::cout << "Row number must be between 0 and 4 (inclusive)." << std::endl;
		exit(1);
	}

	// Open the TTree
	TTree *t = (TTree*)f->Get("fin_tree");

	// Get the run number
	Int_t run_number = GetRunNumber( (TString)f->GetName() );

	// Define the xcal conditions
	Double_t step_size = ( xcal_max - xcal_min )/num_strips_per_Si;

	// Populate an array holding the conditions
	Double_t xcal_values[num_strips_per_Si + 1];
	for ( Int_t i = 0; i < num_strips_per_Si + 1; i++ ){
		if ( i == 0 ){
			xcal_values[i] = xcal_min;
		}
		else if ( i == num_strips_per_Si ){
			xcal_values[i] = xcal_max;
		}
		else{
			xcal_values[i] = xcal_values[i-1] + step_size;
		}
	}
	
	// Create loop variables
	TH1F *h_strip[ ( index1 - index0 ) ][ num_detectors_per_row ][num_strips_per_Si ];
	TCanvas *c_strip[ ( index1 - index0 ) ][ num_detectors_per_row ][num_strips_per_Si ];

	// Loop over the number of rows
	for ( Int_t k = index0; k < index1; k++ ){

		// Loop over the number of detectors
		for ( Int_t i = 0; i < num_detectors_per_row; i++ ){

			// Loop over the number of mini-strips
			for ( Int_t j = 0; j < num_strips_per_Si; j++ ){
		
				// Create a canvas
				c_strip[k][i][j] = new TCanvas( Form( "c_strip%i_%i", num_detectors_per_row*k + i, j ), Form( "DET #%i (%i/%i)", num_detectors_per_row*k + i, j+1, num_strips_per_Si), CANVAS_WIDTH, CANVAS_HEIGHT );

				// Create a spectrum
				TString hist_name = Form( "h_strip%i_%i", num_detectors_per_row*k + i, j );
				t->Draw( Form( "%s>>%s%s", branch.Data(), hist_name.Data(), bins.Data() ), Form( "%sdetID == %i && xcal[] >= %f && xcal[] < %f", cut_string.Data(), num_detectors_per_row*k + i, xcal_values[j], xcal_values[j+1] ) );

				// Get the histogram
				h_strip[k][i][j] = (TH1F*)gDirectory->Get( hist_name );

				// Set the histogram title
				//h_strip[k][i][j]->SetTitle( Form( "Det #%i (%i/%i); Energy; # Counts", num_detectors_per_row*k + i, j + 1, num_strips_per_Si) );

				// Define a specific directory
				TString specific_dir = Form( "%s/%s%i" , print_dir.Data(), key_word.Data(), run_number );

				// Check if the print directory exists
				if ( gSystem->Exec( Form( "[ -d %s ]", specific_dir.Data() ) ) != 0 ){
					std::cout << std::endl << "Directory doesn't exist. Making directory..." << std::endl << std::endl;
					gSystem->Exec( Form( "mkdir -p %s", specific_dir.Data() ) );
				}

				// Print the spectrum to a file (and suppress print message)
				Int_t old_level = gErrorIgnoreLevel;
				gErrorIgnoreLevel = kWarning;
				TString file_name = Form( "%s%i_Det%i-%i", key_word.Data(), run_number, num_detectors_per_row*k + i, j );
				c_strip[k][i][j]->Print( Form( "%s/%s.pdf", specific_dir.Data(), file_name.Data() ) );
				gErrorIgnoreLevel = old_level;

				// Write the spectrum to a radware file
				WriteSPE( h_strip[k][i][j]->GetName(), Form( "%s/%s", specific_dir.Data(), file_name.Data() ) );
			}
		}
	}
	return;
}
