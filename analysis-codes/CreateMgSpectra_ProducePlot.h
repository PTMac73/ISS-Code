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

#include <TCanvas.h>
#include <TFile.h>
#include <TText.h>
#include <TTree.h>
#include <TString.h>
#include "CreateMgSpectra.h"
#include "WriteSPE.h"


void ProducePlot( TCanvas *c, TTree *t, TH1F *h, TString cut_string, Int_t spec_layout, Int_t spec_variant, Int_t strip_number, Int_t id_number, Int_t pos_number, TString label, Int_t cuts_mode ){
	// Make the spectra
	c = new TCanvas( MakeObjectName( 1, spec_variant, id_number, strip_number, cuts_mode ).Data(), "CANVAS", CANVAS_WIDTH, CANVAS_HEIGHT );
	t->Draw( Form( "Ex_corrected>>%s(450, -1, 8)", MakeObjectName( 0, spec_variant, id_number, strip_number, cuts_mode ).Data() ), cut_string );

	// Print the cut string
	std::cout << cut_string << "\n";

	// Get the histogram
	h = (TH1F*)gDirectory->Get( MakeObjectName( 0, spec_variant, id_number, strip_number, cuts_mode ) );

	// Format the histogram
	h->SetTitle( cut_string.Data() );
	h->SetTitleFont(62);
	
	h->GetXaxis()->SetLabelFont(62);
	h->GetXaxis()->SetTitleFont(62);
	h->GetXaxis()->SetTitle("Excitation Energy (MeV)");
	h->GetXaxis()->CenterTitle();
	h->GetYaxis()->SetLabelFont(62);
	h->GetYaxis()->SetTitleFont(62);
	h->GetYaxis()->SetTitle("Counts per 20 keV");
	h->GetYaxis()->CenterTitle();
	
	TPad* pad = (TPad*)c->GetPad(0);
	pad->SetLeftMargin(0.08);
	pad->SetTopMargin(0.06);
	pad->SetBottomMargin(0.08);
	pad->SetRightMargin(0.02);
	

	// Print the canvas and write the spectrum
	c->Print( Form( "%s/%s", print_dir.Data(), MakePrintFileName( pos_number, label, id_number, strip_number, NUM_STRIPS_PER_SI, ".pdf", spec_layout, cuts_mode ).Data() ), "EmbedFonts" );

	WriteSPE( h->GetName(), Form( "%s/%s", print_dir.Data(), MakePrintFileName( pos_number, label, id_number, strip_number, NUM_STRIPS_PER_SI, "", spec_layout, cuts_mode ).Data() ) );

	std::cout << ">>> Completed position " << pos_number << ", " << ( spec_layout == 1 ? "row" : "det" ) << id_number << ", strip " << strip_number + 1 << "/" << NUM_STRIPS_PER_SI << " for " << label << std::endl << std::endl;

	
	// Write the histogram to file
	h->Write();

}


#endif





































