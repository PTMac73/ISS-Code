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
#include <TH1F.h>
#include <TH2F.h>
#include <TFile.h>
#include <TObjArray.h>
#include <TROOT.h>
#include <TString.h>
#include <TStyle.h>

#include "../AT_Settings.h"
#include "AT_HistogramGlobals.h"

TCanvas* c_test_ex;
TCanvas* c_test_evz;

// h is the rdt 2d spectrum
Int_t GetRDTNumber( TH2* h ){
	TString rdt_str = h->GetYaxis()->GetTitle();
	if ( ((TString)rdt_str[4]).IsDigit() ){
		return (Int_t)( rdt_str[4] - '0' );
	}
	else{
		return 0;
	}
}



// h is the rdt cut histogram
// c is the canvas defined for it
// cut_name is the name of the cut before outputting
// t is the fin_tree
// arr is a TOBjArray that contains additional things
TCutG* CreateCut( TH2 *h, TCanvas* c, TString cut_name = "cut", TTree* t = NULL, TObjArray *arr = NULL ){
	TCutG * cut0 = NULL; TCutG * cut1 = NULL; TCutG * cut2 = NULL;
	
	// Draw the histogram and add a toolbar
	if( !c->GetShowToolBar() ){ c->ToggleToolBar(); }
	h->Draw();
	
	Bool_t found_cut = 0;
	c_test_ex = new TCanvas( "c_test_ex", "TEST THE CUTS", C_WIDTH, C_HEIGHT );
	c_test_evz = new TCanvas( "c_test_evz", "TEST THE EVZ", C_WIDTH, C_HEIGHT );
	GlobSetCanvasMargins( c_test_ex );
	GlobSetCanvasMargins( c_test_evz );
	
	TH1F* h0 = NULL;
	TH1F* h1 = NULL;
	TH1F* h2 = NULL;
	TH1F* h_arr_ex; TH2F* h_arr_evz;
	TH2F* hevz0; TH2F* hevz1;
	Int_t cut_ctr = 0;
	
	h_arr_ex = (TH1F*)arr->FindObject( Form( "h_rdt_ex_mg_%i", GetRDTNumber(h) ) );
	c_test_ex->cd(); h_arr_ex->Draw();

	h_arr_evz = (TH2F*)arr->FindObject( Form( "h_rdt_evz_mg_%i", GetRDTNumber(h) ) );
	c_test_evz->cd(); h_arr_evz->Draw();
	
	while ( !found_cut ){
		// Update the rdt canvas and wait for a cut
		c->cd(); c->Modified(); c->Update();
		std::cout << "Draw a cut." << "\n";
		gPad->WaitPrimitive("CUTG");

		// Find the cut after it's been made
		cut0 = (TCutG*)gROOT->FindObject("CUTG");
		std::cout << std::left << std::setw(40) << std::setfill('.') << "Found cut";
		cut0->SetName("cut_test0");
		cut0->SetVarX( h->GetXaxis()->GetTitle() );
		cut0->SetVarY( h->GetYaxis()->GetTitle() );
		cut0->SetLineColor(kRed);
		cut0->SetLineWidth(2);
		c->Modified(); c->Update();
		
		// Draw the excitation spectrum and EVZ spectra
		c_test_ex->cd();
		std::cout << std::left << std::setw(40) << std::setfill('.') << "Drawing updated ex spectrum" << "\n";
		t->Draw("Ex>>h_test_ex0(450, -1, 8 )", "cut_test0", "goff" );
		h0 = (TH1F*)gDirectory->Get("h_test_ex0");
		h0->SetLineColor(kWhite);
		h0->SetLineWidth(0);
		h0->SetFillColor(kRed);
		h0->Draw("SAME");
		c_test_ex->Modified(); c_test_ex->Update();
		
		c_test_evz->cd();
		t->Draw("ecrr:z>>h_test_evz0(400, -50, -10, 450, 0, 9 )", "cut_test0", "goff" );
		std::cout << std::left << std::setw(40) << std::setfill('.') << "Drawing updated evz spectrum";
		hevz0 = (TH2F*)gDirectory->Get("h_test_evz0");
		hevz0->SetMarkerColor(kRed);
		hevz0->SetMarkerStyle(20);
		hevz0->SetMarkerSize(0.5);
		hevz0->Draw("SAME");
		c_test_evz->Modified(); c_test_evz->Update();
		std::cout << std::left << std::setw(40) << std::setfill('.') << "Updated all plots" << "\n";
		
		// Ask if they want to continue
		Bool_t valid_ans = 0;
		while ( !valid_ans ){
			std::cout << "Do you want to save this cut?[Y/N]" << "\n";
			TString ans;
			std::cin >> ans;
			
			if ( ans[0] == 'Y' || ans[0] == 'y' ){
				valid_ans = 1; found_cut = 1;
				std::cout << "Saving cut..." << "\n";
			}
			else if ( ans[0] == 'N' || ans[0] == 'n' ){
				std::cout << "Preparing for new cut..." << "\n";
				valid_ans = 1;
				
				// Reset ex spectrum for next time
				c_test_ex->Clear();
				h_arr_ex->Draw();
				
				if ( h1 != NULL ){
					h2 = (TH1F*)h1->Clone("h_test_ex2");
					h2->SetFillColor(kTeal);
					h2->SetTitle( Form("Cut %d", cut_ctr ) );
					h2->Draw("SAME");
				}
				
				if ( h0 != NULL ){
					h1 = (TH1F*)h0->Clone("h_test_ex1");
					h1->SetFillColor(kViolet);
					h1->Draw("SAME");
				}
				
				h0->Delete();
				c_test_ex->Modified(); c_test_ex->Update();
				cut_ctr++;
				
				
				// Reset rdt spectrum for next time
				c->cd();
				c->Clear();
				h->Draw();
				
				if ( cut1 != NULL ){
					cut2 = (TCutG*)cut1->Clone("cut_test2");
					cut2->SetLineColor(kTeal);
					cut2->Draw("SAME");
				}
				if ( cut0 != NULL ){
					cut1 = (TCutG*)cut0->Clone("cut_test1");
					cut1->SetLineColor(kViolet);
					cut1->SetLineWidth(2);
					cut1->Draw("SAME");
				}
				
				cut0->Delete();
				c->Modified(); c->Update();
				
				
				// Reset EVZ spectrum for next time
				hevz1 = (TH2F*)hevz0->Clone("h_test_evz1");
				hevz1->SetMarkerColor(kViolet);
				hevz1->SetTitle( Form("Cut %d", cut_ctr ) );
				c_test_evz->Clear();
				h_arr_evz->Draw();
				hevz1->Draw("SAME");
				c_test_evz->Modified(); c_test_evz->Update();
				hevz0->Delete();
			}
			else{
				std::cout << "Not a valid answer. Try again." << "\n";
			}
		}
	}
	
	cut0->SetName(cut_name);
	cut0->SetVarX( h->GetXaxis()->GetTitle() );
	cut0->SetVarY( h->GetYaxis()->GetTitle() );
	cut0->SetTitle("");
	cut0->SetLineColor(1);
	
	std::cout << "Cut " << std::right << std::setw(10) << cut_name << " created.\n";
	return cut0;
}

void WriteCutFile( TObjArray* arr ){
	TFile *f = new TFile( cut_dir_si.Data(), "RECREATE" );
	if ( f->IsOpen() ){
		arr->Write( arr->GetName(), TObject::kSingleKey );
		std::cout << "Cuttlefish created!" << "\n";
		f->Close();
	}
	else{
		std::cout << "Could not write to cut file!" << "\n";
	}
	
	return;
}






#endif
