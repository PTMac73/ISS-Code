// CreateMgSpectra_ProducePlot.C
// Produces a plot based on the information fed to it
// ============================================================================================= //
// Patrick MacGregor
// Nuclear Physics Research Group
// School of Physics and Astronomy
// The University of Manchester
// ============================================================================================= //
#ifndef CREATE_MG_SPECTRA_PRODUCE_PLOT_H_
#define CREATE_MG_SPECTRA_PRODUCE_PLOT_H_


#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TString.h>
#include "CreateMgSpectra.h"
#include "WriteSPE.h"


void ProducePlot( TCanvas *c, TTree *t, TH1F *h, TString cut_string, Int_t spectrum_type, Int_t num_spectrum_type, Int_t strip_number, Int_t id_number, Int_t pos_number, TString label ){
	// Make the spectra
	c = new TCanvas( MakeObjectName( 1, num_spectrum_type, id_number, strip_number ).Data(), "CANVAS", CANVAS_WIDTH, CANVAS_HEIGHT );
	t->Draw( Form( "Ex>>%s(450, -1, 8)", MakeObjectName( 0, num_spectrum_type, id_number, strip_number ).Data() ), cut_string );

	// Print the cut string
	std::cout << cut_string << "\n";

	// Get the histogram
	h = (TH1F*)gDirectory->Get( MakeObjectName( 0, num_spectrum_type, id_number, strip_number ) );

	// Format the histogram
	h->SetTitle( Form( "ID = %i, strip %i/%i for %s; Ex / MeV; #", id_number, strip_number + 1, NUM_STRIPS_PER_SI, label.Data() ) );

	// Print the canvas and write the spectrum
	c->Print( Form( "%s/%s", print_dir.Data(), MakePrintFileName( pos_number, label, id_number, strip_number, NUM_STRIPS_PER_SI, ".pdf", spectrum_type ).Data() ), "EmbedFonts" );

	WriteSPE( h->GetName(), Form( "%s/%s", print_dir.Data(), MakePrintFileName( pos_number, label, id_number, strip_number, NUM_STRIPS_PER_SI, "", spectrum_type ).Data() ) );

	std::cout << ">>> Completed position " << pos_number << ", " << ( spectrum_type == 1 ? "row" : "det" ) << id_number << ", strip " << strip_number + 1 << "/" << NUM_STRIPS_PER_SI << " for " << label << std::endl << std::endl;

	
	// Write the histogram to file
	h->Write();

}


#endif





































