// CreateMgSpectra.C
// Creates a number of spectra for the Mg analysis - allows for creating strips within the detector
// as well.
// ============================================================================================= //
// Patrick MacGregor
// Nuclear Physics Research Group
// Department of Physics and Astronomy
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
#include "CreateMgSpectra_ProducePlot.h"
	/*TODO CHECK THAT THIS ALL WORKS*/
void CreateMgSpectra( TFile *f, Int_t det_mode = DETECTOR_MODE, Int_t spec_layout = SPECTRUM_LAYOUT, Int_t cuts_mode = CUTS_MODE ){
	// Set into batch mode
	gROOT->SetBatch(kTRUE);

	// Change the style
	gStyle->SetOptStat("meni");

	// Globals - decide what kind of operation is being run (best detectors, all detectors etc.)
	const Int_t NUM_SPECTRUM_VARIANTS = GetNumSpectrumVariants( det_mode );
	
	// Define objects to hold final information
	TString label[NUM_SPECTRUM_VARIANTS];
	Bool_t use_detector[NUM_DETECTORS][NUM_SPECTRUM_VARIANTS];
	
	// Copy the arrays into a global variable
	if ( det_mode == 0 ){
		std::copy( &use_detector_1[0][0], &use_detector_1[0][0] + NUM_DETECTORS*NUM_SPECTRUM_VARIANTS, &use_detector[0][0] );
		std::copy( &label_1[0], &label_1[0] + NUM_SPECTRUM_VARIANTS, &label[0] );
	}
	else if ( det_mode == 1 ){
		std::copy( &use_detector_2[0][0], &use_detector_2[0][0] + NUM_DETECTORS*NUM_SPECTRUM_VARIANTS, &use_detector[0][0] );
		std::copy( &label_2[0], &label_2[0] + NUM_SPECTRUM_VARIANTS, &label[0] );
	}
	else{
		std::cout << "DET_MODE " << det_mode << " NOT DEFINED. FORCE EXIT" << std::endl;
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
	TString cut_string_i = "( ";
	if ( cuts_mode == 0 ){ cut_string_i += "cut0 || cut1 || cut2 || cut3 ) && ( "; }
	TString cut_string_m = " ) && ";
	if ( cuts_mode == 0 ){ cut_string_m += "td_rdt_e[] >= td_rdt_e_cuts[][0] && td_rdt_e[] <= td_rdt_e_cuts[][1] && "; }

	// Define canvas and histogram things
	TCanvas *c_spectrum[NUM_SPECTRUM_VARIANTS][NUM_STRIPS_PER_SI][NUM_ROWS][NUM_DETECTORS_PER_ROW];
	TH1F *h_spectrum[NUM_SPECTRUM_VARIANTS][NUM_STRIPS_PER_SI][NUM_ROWS][NUM_DETECTORS_PER_ROW];

	// Define a bool to handle the cut strings for multiple detectors
	Bool_t detector_used = 0;
	
	// Define an output file to save each histogram
	TFile *f_out = new TFile( MakeRootFileName( pos_number, spec_layout, cuts_mode ), "RECREATE" );
	std::cout << "Beginning loops....." << "\n";
	// Loop over different types of spectra, rows, detectors per row, and number of strips ( here e.g. states1-6 )
	for ( Int_t k = 0; k < NUM_SPECTRUM_VARIANTS; k++ ){
		
		// LOOP OVER THE NUMBER OF STRIPS
		for ( Int_t l = 0; l < NUM_STRIPS_PER_SI; l++ ){

			// LOOP OVER # ROWS
			for ( Int_t i = 0; i < NUM_ROWS; i++ ){
				
				// LOOP OVER # DETECTORS PER ROW
				for ( Int_t j = 0; j < NUM_DETECTORS_PER_ROW; j++ ){
					// CONSTRUCT THE CUT STRINGS
					// Construct RBR cut string
					if ( spec_layout == 0 ){
						if ( j == 0 ){
							cut_string = cut_string_i;
							detector_used = 0;
						}
						if ( use_detector[i + NUM_ROWS*j][k] == 1 ){
								cut_string += Form( "detID == %i", i + NUM_ROWS*j );
								detector_used = 1;
						}
						if ( use_detector[i + NUM_ROWS*(j+1)][k] == 1 && j != NUM_DETECTORS_PER_ROW - 1 && detector_used == 1 ){
							cut_string += " || ";
						}
						
					}
					
					// Construct DBD cut strings and plot
					else if ( spec_layout == 1 ){
						// Construct DBD cut string and plot
						if ( use_detector[i + NUM_ROWS*j][k] == 1 ){
							cut_string = cut_string_i + Form( "detID == %i", i + NUM_ROWS*j ) + cut_string_m + GenerateXCALCutString( NUM_STRIPS_PER_SI, l );
							ProducePlot( c_spectrum[k][l][i][j], t, h_spectrum[k][l][i][j], cut_string, spec_layout, k, l, i + NUM_ROWS*j, pos_number, label[k], cuts_mode );
						}
					
					}
					
					// Construct full cut string
					else if ( spec_layout == 2 ){
						if ( i == 0 && j == 0 ){
							cut_string = cut_string_i;
							detector_used = 0;
						}
						if ( use_detector[i*NUM_DETECTORS_PER_ROW + j][k] == 1 ){
								cut_string += Form( "detID == %i", i + NUM_ROWS*j );
								detector_used = 1;
						}
						if ( use_detector[i*NUM_DETECTORS_PER_ROW + j + 1][k] == 1 && i*NUM_DETECTORS_PER_ROW + j + 1 < NUM_DETECTORS && detector_used == 1 ){
							cut_string += " || ";
						}
					}
					
				}
			
				if ( spec_layout == 0 ){
					// Plot the RBR detector from fully constructed cut string
					cut_string += cut_string_m;

					// Now make the final part of the cut string
					cut_string += GenerateXCALCutString( NUM_STRIPS_PER_SI, l ) + Form( " && thetaCM >= %f", THETA_CM_LIMIT[cuts_mode] );// + GenerateThetaCMCutString(9);
					//ProducePlot( c_spectrum[k][l][i][0], t, h_spectrum[k][l][i][0], cut_string, spec_layout, k, l, i, pos_number, label[k], cuts_mode );
					std::cout << cut_string << "\n";
				}
				
			} // LOOP OVER ROWS (i)
			
		} // LOOP OVER STRIPS (l)
		
	} // LOOP OVER VARIANTS (k)
	
	
	if ( spec_layout == 2 ){
		// Construct the excitation spectrum from everything
		cut_string += cut_string_m;
		cut_string += GenerateXCALCutString( NUM_STRIPS_PER_SI, 0 ) + Form( " && thetaCM >= %f", THETA_CM_LIMIT[cuts_mode] );
		std::cout << "Producing plot ..." << "\n";
		ProducePlot( c_spectrum[0][0][0][0], t, h_spectrum[0][0][0][0], cut_string, spec_layout, 0, 0, 0, pos_number, label[0], cuts_mode );

	}

	// Close the output file
	f_out->Close();
	
	return;
}

