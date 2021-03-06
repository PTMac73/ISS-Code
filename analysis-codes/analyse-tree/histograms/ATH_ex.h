// ATH_ex.h
// Excitation spectrum on a row-by-row/full basis
// ============================================================================================= //
// Patrick MacGregor
// Nuclear Physics Research Group
// Department of Physics and Astronomy
// The University of Manchester
// ============================================================================================= //

#include <TCanvas.h>
#include <TCutG.h>
#include <TFile.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <THStack.h>
#include <TString.h>
#include <TStyle.h>


#include <iostream>

#include "AT_HistogramGlobals.h"
#include "../AT_Settings.h"
#include "WriteSPE.h"


#ifndef ATH_EX_H_
#define ATH_EX_H_

bool b_get_widths = 1;			// Tells program to get FWHM for the first two excited states. Only works for DBD mode
void GetWidth(std::ostream &f = std::cout, TH1F* h = NULL, Int_t det_num = -1);	// Function call

// Switch for this is SW_EX
// Number of histograms = 24 (DBD) + 6 (RBR) + 1 (Full)
TH1F* h_ex_full = NULL;			// Full spectrum
TH1F* h_ex_full_corr = NULL;
TH1F* h_ex_full_best = NULL;	// Full spectrum (best)
TH1F* h_ex_rbr[6][2];			// RBR spectrum: [0] is full, [1] is best resolution detectors
TH1F* h_ex_dbd[24];				// DBD spectrum
TH1F* h_ex_full_evolution[5];	// Full spectrum evolution
TH1F* h_ex_before[4];
TH1F* h_ex_after[4];


Bool_t ex_print_opt[7] = {
	0,	// (0) full spectrum
	0,	// (1) RBR spectrum (all good dets)
	0,	// (2) RBR spectrum (best resolution dets)
	1,	// (3) DBD spectrum
	0,	// (4) spectrum evolution
	0,	// (5) full spectrum (best)
	0   // (6) full spectrum corrected with energy calibration
};

void HCreateEx(){
	// *LOOP* over detectors
	for ( Int_t i = 0; i < 24; i++ ){

		// Full spectrum
		if ( i == 0 && ALL_ROWS == 1 ){
			CreateExSpectrum( h_ex_full, "h_ex_full" );
			CreateExSpectrum( h_ex_full_best, "h_ex_full_best" );
			CreateExSpectrum( h_ex_full_corr, "h_ex_full_corr" );
		}

		// Evolution
		if ( i == 0 ){
			for ( Int_t j = 0; j < 5; j++ ){ CreateExSpectrum( h_ex_full_evolution[j], Form( "h_ex_full_evolution_%i", j ) ); }
		}

		// Row by row
		if ( ( i == ROW_NUMBER || ( ROW_NUMBER == -1 && i < 6 ) ) && ROW_BY_ROW == 1 ){
			CreateExSpectrum( h_ex_rbr[i][0], Form( "h_ex_rbr_%i", i ) );
			CreateExSpectrum( h_ex_rbr[i][1], Form( "h_ex_rbr_%i_best", i ) );
		}

		// Detector by detector
		if ( ( i == DET_NUMBER || ( DET_NUMBER == -1 ) ) && DET_BY_DET == 1 ){
			CreateExSpectrum( h_ex_dbd[i], Form( "h_ex_dbd_%i", i ) );

		}
	}
	return;
}


void HDrawEx(){
	// Define some local variables
	TCanvas* c_ex_full;
	TCanvas* c_ex_full_best;
	TCanvas* c_ex_full_corr;
	TCanvas* c_ex_rbr[6][2];
	TCanvas* c_ex_dbd[24];
	TCanvas* c_ex_dbd_comb;
	TCanvas* c_ex_evolution[4];

	if ( CANVAS_COMBINE == 1 ){
		c_ex_dbd_comb = new TCanvas( "c_ex_dbd_comb", "DBD Ex Spectrum", C_WIDTH, C_HEIGHT );
		c_ex_dbd_comb->Divide(6,4);
	}

	TString root_name = Form( "%s/pos%i_ex", print_dir.Data(), ARR_POSITION );
	TFile* f;
	TString spec_name;
	std::ofstream fwhm_data_file;
	fwhm_data_file.open( Form("%s/fwhm_data_file.txt", print_dir.Data()) );

	// Open root file if desired
	if ( SW_EX[1] == 1 && PRINT_ROOT == 1 ){
		f = new TFile( ( root_name + ".root" ).Data(), "RECREATE" );
	}


	// *LOOP* over detectors (i)
	for ( Int_t i = 0; i < 24; i++ ){

		// Full spectrum
		if ( i == 0 && ALL_ROWS == 1 ){
			// Plot spectrum
			c_ex_full = new TCanvas( "c_ex_full", "Full excitation spectrum",  C_WIDTH, C_HEIGHT );
			GlobSetCanvasMargins( c_ex_full );
			h_ex_full->Draw();

			c_ex_full_corr = new TCanvas( "c_ex_full_corr", "Full excitation spectrum (corrected)",  C_WIDTH, C_HEIGHT );
			GlobSetCanvasMargins( c_ex_full_corr );
			h_ex_full_corr->Draw();

			c_ex_full_best = new TCanvas( "c_ex_full_best", "Full excitation spectrum (best res)",  C_WIDTH, C_HEIGHT );
			GlobSetCanvasMargins( c_ex_full_best );
			h_ex_full_best->Draw();

			spec_name = Form( "%s/pos%i_ex_full", print_dir.Data(), ARR_POSITION );

			// Plot evolution spectrum
			for ( Int_t j = 0; j < 4; j++ ){
				c_ex_evolution[j] = new TCanvas( Form( "c_ex_evolution_%i", j ), Form( "EXCITATION SPECTRUM EVOLUTION | Case %i", j ), C_WIDTH, C_HEIGHT );
				GlobSetCanvasMargins( c_ex_evolution[j] );
				h_ex_before[j] = (TH1F*)h_ex_full_evolution[j]->Clone();
				h_ex_after[j] = (TH1F*)h_ex_full_evolution[j+1]->Clone();

				h_ex_before[j]->SetLineColor(kGray+2);
				h_ex_before[j]->SetFillColor(kGray+2);
				h_ex_after[j]->SetLineColor(kBlack);
				h_ex_after[j]->SetFillColor(kYellow);
				THStack *hs = new THStack("hs","");
				hs->Add(h_ex_before[j]);
				hs->Add(h_ex_after[j]);
				hs->Draw("nostack");
				GlobSetHistFonts(hs);
				hs->GetXaxis()->SetTitle("E_{x} (MeV)");
				hs->GetYaxis()->SetTitle("Counts per 20 keV");
			}


			// Print spectrum if desired
			if ( SW_EX[1] == 1 ){
				if( ex_print_opt[0] == 1 ){ PrintAll( c_ex_full, spec_name ); }
				if( ex_print_opt[5] == 1 ){ PrintAll( c_ex_full_best, spec_name + "_best" ); }
				if( ex_print_opt[6] == 1 ){ PrintAll( c_ex_full_corr, spec_name + "_corr" ); }
				if ( PRINT_ROOT == 1 ){
					f->cd();
					if( ex_print_opt[0] == 1 ){ h_ex_full->Write(); }
					if( ex_print_opt[5] == 1 ){ h_ex_full_best->Write(); }
					if( ex_print_opt[6] == 1 ){ h_ex_full_corr->Write(); }
				}

				for ( Int_t j = 0; j < 4; j++ ){
					if( ex_print_opt[4] == 1 ){
						PrintAll( c_ex_evolution[j], Form( "%s/pos%i_ex_evolution_%i", print_dir.Data(), ARR_POSITION, j ) );
						if ( PRINT_ROOT == 1 ){ f->cd(); h_ex_full_evolution[j]->Write(); }
					}
				}
			}

			// Write SPE file if desired
			if ( SW_EX[2] == 1 ){
				if( ex_print_opt[0] == 1 ){ WriteSPE( h_ex_full->GetName(), Form( "%s", spec_name.Data() ) ); }
				if( ex_print_opt[5] == 1 ){ WriteSPE( h_ex_full_best->GetName(), Form( "%s", ( spec_name + "_best" ).Data() ) ); }
			}

		}


		// Row by row
		if ( ( i == ROW_NUMBER || ( ROW_NUMBER == -1 && i < 6 ) ) && ROW_BY_ROW == 1 ){
			// Plot spectrum
			c_ex_rbr[i][0] = new TCanvas( Form( "c_ex_rbr_%i_LO", i ), Form( "Ex | Row %i | All rows", i ), C_WIDTH, C_HEIGHT );
			GlobSetCanvasMargins( c_ex_rbr[i][0] );
			h_ex_rbr[i][0]->Draw();

			c_ex_rbr[i][1] = new TCanvas( Form( "c_ex_rbr_%i_HI", i ), Form( "Ex | Row %i | Best dets", i ), C_WIDTH, C_HEIGHT );
			GlobSetCanvasMargins( c_ex_rbr[i][1] );
			h_ex_rbr[i][1]->Draw();

			spec_name = Form( "%s/pos%i_ex_rbr_%i", print_dir.Data(), ARR_POSITION, i );

			// Print spectrum if desired
			if ( SW_EX[1] == 1 ){
				if( ex_print_opt[1] == 1 ){ PrintAll( c_ex_rbr[i][0], spec_name + "" ); }
				if( ex_print_opt[2] == 1 ){ PrintAll( c_ex_rbr[i][1], spec_name + "_best" ); }
				if ( PRINT_ROOT == 1 ){
					f->cd();
					if( ex_print_opt[1] == 1 ){ h_ex_rbr[i][0]->Write(); }
					if( ex_print_opt[2] == 1 ){ h_ex_rbr[i][1]->Write(); }
				}
			}

			// Write SPE file if desired
			if ( SW_EX[2] == 1 ){
				if( ex_print_opt[1] == 1 ){ WriteSPE( h_ex_rbr[i][0]->GetName(), Form( "%s", ( spec_name + "" ).Data() ) ); }
				if( ex_print_opt[2] == 1 ){ WriteSPE( h_ex_rbr[i][1]->GetName(), Form( "%s", ( spec_name + "_best" ).Data() ) ); }
			}

		}


		// Detector by detector
		if ( ( i == DET_NUMBER || ( DET_NUMBER == -1 ) ) && DET_BY_DET == 1 ){
			spec_name = Form( "%s/pos%i_ex_dbd_%02i", print_dir.Data(), ARR_POSITION,  i );
			if ( CANVAS_COMBINE == 0 ){
				// Plot spectrum
				c_ex_dbd[i] = new TCanvas( Form( "c_ex_dbd_%i", i ), Form( "Ex | Det %i", i ), C_WIDTH, C_HEIGHT );
				GlobSetCanvasMargins( c_ex_dbd[i] );
				h_ex_dbd[i]->Draw();
				gStyle->SetTitleFont(62);

				// Print spectrum if desired
				if ( SW_EX[1] == 1 ){
					if( ex_print_opt[3] == 1 ){ PrintAll( c_ex_dbd[i], spec_name ); }
				}

				if ( b_get_widths == 1 && fwhm_data_file.is_open() ){
					GetWidth( fwhm_data_file, h_ex_dbd[i], i);
				}

			}
			else{
				c_ex_dbd_comb->cd(i+1);
				h_ex_dbd[i]->SetTitle( Form( "Det #%i", i ) );
				h_ex_dbd[i]->Draw();
				TPad* pad = (TPad*)c_ex_dbd_comb->GetPad(i+1);
				SetCanvasTitleFont( pad );

				// Print spectrum if desired
				if ( SW_EX[1] == 1 && i == 23 ){
					if( ex_print_opt[3] == 1 ){ PrintAll( c_ex_dbd_comb, Form( "%s/pos%i_ex_dbd_comb", print_dir.Data(), ARR_POSITION ) ); }
				}


			}

			// Write ROOT file if desired
			if ( PRINT_ROOT == 1 && SW_EX[1] == 1 && ex_print_opt[3] == 1 ){ f->cd(); h_ex_dbd[i]->Write(); }

			// Write SPE file if desired
			if ( SW_EX[2] == 1 && ex_print_opt[3] == 1 ){ WriteSPE( h_ex_dbd[i]->GetName(), Form( "%s", spec_name.Data() ) ); }
		}


	}
	// Close the root file
	if ( SW_EX[1] == 1 && PRINT_ROOT == 1 && f != NULL ){ if ( f->IsOpen() ){ f->Close(); } }

	// Close the text file
	if ( fwhm_data_file.is_open() ){ fwhm_data_file.close(); }
	return;
}


void GetWidth( std::ostream &f, TH1F* h, Int_t det_num ){
	// Define two Gaussians with same width and flat bg
	TF1* fit_func = new TF1( "fit_func", "[0]*exp(-0.5*((x-[1])/[2])**2) + [3]", -0.4, 0.4 );
	fit_func->SetNpx(2000);
	fit_func->SetLineColor(kBlue);

	// Define initial values
	//fit_func->SetParLimits(0,0,100);
	Int_t max_bin = 0;
	Double_t num = 0;	// Calculate weighted average position
	Double_t den = 0;
	for( Int_t i = 5; i < 47; i++ ){
		num += h->GetBinCenter(i)*h->GetBinContent(i);
		den += h->GetBinContent(i);
		if ( h->GetBinContent(i) > h->GetBinContent(max_bin) || max_bin == 0 ){
			max_bin = i;
		}
	}

	fit_func->SetParameter( 0, h->GetBinContent(max_bin) );
	fit_func->SetParLimits( 0, 0.5*h->GetBinContent(max_bin),1.5*h->GetBinContent(max_bin) );
	fit_func->SetParameter( 1, num/den );
	fit_func->SetParLimits( 2, 0,0.5 );
	fit_func->SetParameter( 2, 0.02 );
	fit_func->FixParameter( 3, 0 );

	// Define a Canvas and draw the histogram
	TCanvas* c_fit = new TCanvas( Form("c_fit_func%i", det_num), "CFITFUNC", C_WIDTH, C_HEIGHT );
	h->Draw();

	// Fit the function
	h->Fit(fit_func, "R");
	//fit_func->Draw("SAME");
	c_fit->Print( Form("%s/cfitfunc%i.png", print_dir.Data(), det_num) );
	Double_t sig = fit_func->GetParameter(2);

	f << "DET: " << std::setw(2) << det_num << " :\t" << std::setw(9) << sig << "\t" << "FWHM:\t" << GetFWHM(sig) << "\n";
	return;
}











#endif
