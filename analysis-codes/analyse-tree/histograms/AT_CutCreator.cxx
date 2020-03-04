// AT_CutCreator.cxx
// Feed in a histogram and get out a cut!
// ============================================================================================= //
// Patrick MacGregor
// Nuclear Physics Research Group
// Department of Physics and Astronomy
// The University of Manchester
// ============================================================================================= //

#ifndef AT_CUT_CREATOR_CXX
#define AT_CUT_CREATOR_CXX

#include <TCanvas.h>
#include <TCutG.h>
#include <TH2F.h>
#include <TFile.h>
#include <TObjArray.h>
#include <TROOT.h>
#include <TString.h>
#include <TStyle.h>

#include "../AT_Settings.h"


TCutG* CreateCut( TH2 *h, TCanvas* c, TString cut_name = "cut" ){
	TCutG * cut = NULL;
	
	// Draw the histogram and add a toolbar
	if( !c->GetShowToolBar() ){ c->ToggleToolBar(); }
	h->Draw();
	
	// Update the canvas
	c->Modified(); c->Update();
	gPad->WaitPrimitive("CUTG");

	// Find the cut after it's been made
	cut = (TCutG*) gROOT->FindObject("CUTG");
	cut->SetName(cut_name);
	cut->SetVarX( h->GetXaxis()->GetTitle() );
	cut->SetVarY( h->GetYaxis()->GetTitle() );
	cut->SetTitle("");
	cut->SetLineColor(1);
	
	std::cout << "Cut " << std::right << std::setw(10) << cut_name << " created.\n";
	return cut;
}

void WriteCutFile( TObjArray* arr ){
	TFile *f = new TFile( cut_dir.Data(), "RECREATE" );
	if ( f->IsOpen() ){
		arr->Write( arr->GetName(), TObject::kSingleKey );
		f->Close();
	}
	
	return;
}






#endif
