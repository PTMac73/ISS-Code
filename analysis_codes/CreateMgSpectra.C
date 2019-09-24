// CreateMgSpectra.C
// Creates a number of spectra for the Mg analysis - allows for creating strips within the detector
// as well.
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
#include "CreateMgSpectra.h"
#include "WriteSPE.h"

void CreateMgSpectra( TFile *f, Int_t mode = MODE ){
	// Set into batch mode
	gROOT->SetBatch(kTRUE);

	// Change the style
	gStyle->SetOptStat("meni");

	// Globals - decide what kind of operation is being run (best detectors, all detectors etc.)
	const Int_t NUM_DIFF_TYPES_SPECTRA = GetNumDiffTypesSpectra( mode );
	
	// Define objects to hold final information
	TString label[NUM_DIFF_TYPES_SPECTRA];
	Bool_t use_detector[NUM_DETECTORS][NUM_DIFF_TYPES_SPECTRA];
	

	// Copy the arrays into a global variable
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
	TString cut_string_m = " ) && thetaCM > 11 && td_rdt_e[] >= td_rdt_e_cuts[][0] && td_rdt_e[] <= td_rdt_e_cuts[][1] && ";

	// Define canvas and histogram things
	TCanvas *c_spectrum[NUM_DIFF_TYPES_SPECTRA][NUM_ROWS][NUM_STRIPS_PER_SI];
	TH1F *h_spectrum[NUM_DIFF_TYPES_SPECTRA][NUM_ROWS][NUM_STRIPS_PER_SI];
	
	// Loop over different types of spectra, rows, detectors per row, and number of strips
	for ( Int_t k = 0; k < NUM_DIFF_TYPES_SPECTRA; k++ ){
		
		// Define an output file to save each histogram
		TFile *f_out = new TFile( Form( "%s/pos%i_%s_hists.root", print_dir.Data(), pos_number, label[k].Data() ), "RECREATE" );

		// LOOP OVER THE NUMBER OF STRIPS
		for ( Int_t l = 0; l < NUM_STRIPS_PER_SI; l++ ){

			// LOOP OVER # ROWS
			for ( Int_t i = 0; i < NUM_ROWS; i++ ){

				// Make the cut string for each row
				cut_string = "";
				cut_string += cut_string_i;

				// See if spectrum will be filled
				Bool_t using_spectra = 0;

				// LOOP OVER # DETECTORS PER ROW
				for ( Int_t j = 0; j < NUM_DETECTORS_PER_ROW; j++ ){
					// Add detectors to the cuts
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
					cut_string += cut_string_m;

					// Now make the final part of the cut string
					cut_string += GenerateXCALCutString( NUM_STRIPS_PER_SI, l );
					std::cout << "CUT: " << cut_string << std::endl;
					
					// Make the spectra
					c_spectrum[k][i][l] = new TCanvas( MakeObjectName( 1, k, i, l ).Data(), "CANVAS", CANVAS_WIDTH, CANVAS_HEIGHT );
					t->Draw( Form( "Ex>>%s(450, -1, 8)", MakeObjectName( 0, k, i, l ).Data() ), cut_string );

					// Get the histogram
					h_spectrum[k][i][l] = (TH1F*)gDirectory->Get( MakeObjectName( 0, k, i, l ) );

					// Format the histogram
					h_spectrum[k][i][l]->SetTitle( Form( "Row %i, strip %i/%i for %s; Ex / MeV; #", i, l+1, NUM_STRIPS_PER_SI, label[k].Data() ) );

					// Print the canvas and write the spectrum
					c_spectrum[k][i][l]->Print( Form( "%s/%s", print_dir.Data(), MakePrintFileName( pos_number, label[k], i, l, NUM_STRIPS_PER_SI, ".pdf").Data() ), "EmbedFonts" );
					WriteSPE( h_spectrum[k][i][l]->GetName(), Form( "%s/%s", print_dir.Data(), MakePrintFileName(pos_number, label[k], i, l, NUM_STRIPS_PER_SI, "").Data() ) );
					std::cout << ">>> Completed position " << pos_number << ", row " << i << ", strip " << l+1 << "/" << NUM_STRIPS_PER_SI << " for " << label[k] << std::endl << std::endl;
					
					// Write the histogram to file
					h_spectrum[k][i][l]->Write();
				}
			}
		}

		// Close the output file
		f_out->Close();

	}
	
	return;
}


































