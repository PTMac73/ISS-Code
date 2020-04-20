#define AnalyseTree_cxx
// The class definition in AnalyseTree.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.


// The following methods are defined in this file:
//	 Begin():		  called every time a loop on the tree starts,
//						  a convenient place to create your histograms.
//	 SlaveBegin():	called after Begin(), when on PROOF called only on the
//						  slave servers.
//	 Process():		called for each event, in this function you decide what
//						  to read and fill your histograms.
//	 SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//						  called only on the slave servers.
//	 Terminate():	 called at the end of the loop on the tree,
//						  a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// root> T->Process("AnalyseTree.C")
// root> T->Process("AnalyseTree.C","some options")
// root> T->Process("AnalyseTree.C+")
//


#include "AnalyseTree.h"
#include "AT_Globals.h"
#include "AT_Histograms.h"
#include "AT_Settings.h"
#include <TCanvas.h>
#include <TCutG.h>
#include <TH1.h>
#include <TH2.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TStopwatch.h>
#include <TStyle.h>

// Progress report
TStopwatch stopwatch;
ULong64_t processed_entries = 0;
ULong64_t num_entries;
Double_t entry_frac = 0.1;



// BEGIN ANALYSIS
void AnalyseTree::Begin(TTree* t)
{
	// The Begin() function is called at the start of the query.
	// When running with PROOF Begin() is only called on the client.
	// The tree argument is deprecated (on PROOF 0 is passed).
	gStyle->SetOptStat(0);
	gStyle->SetTitleFont(62);
	
	if ( DISPLAY_CANVAS == 1 ){ gROOT->SetBatch(kFALSE); }
	else{ gROOT->SetBatch(kTRUE); }
	
	// Print summary of options
	PrintSummaryOfOptions();
	
	// Create histograms
	if ( SW_EX_COMPARE[0] == 1 ){ HCreateExCompare(); }
	if ( SW_RDT_CUTS[0] == 1 ){ HCreateRDTCuts(); }
	if ( SW_EVZ_COMPARE[0] == 1 ){ HCreateEVZCompare(); }
	if ( SW_EVZ[0] == 1 ){ HCreateEVZ(); }
	if ( SW_EVZ_SI[0] == 1 ){ HCreateEVZSi(); }
	if ( SW_EX_SI[0] == 1 ){ HCreateExSi(); }
	if ( SW_EX[0] == 1 ){ HCreateEx(); }
	if ( SW_XNXF[0] == 1 ){ HCreateXNXF(); }
	if ( SW_XCAL[0] == 1 ){ HCreateXCAL(); }
	
	// Get the number of entries
	num_entries = t->GetEntries();
	TString option = GetOption();
	
	// Get the cuts
	TFile *f = new TFile( cut_dir.Data() );
	if ( f->IsOpen() ){
		cut_list = (TObjArray*)f->FindObjectAny("cutList");
		f->Close();
	}
	else{
		std::cout << "NO MG CUTS" << "\n";
		std::exit(1);
	}

	// Get the Si cuts
	TFile* f_si = new TFile( cut_dir_si.Data() );
	if ( f_si->IsOpen() ){
		cut_list_si = (TObjArray*)f_si->FindObjectAny("cuttlefish");
		f_si->Close();
		found_si_cuts = 1;
	}
	else{
		std::cout << "NO SI CUTS FOUND. Will carry on without them." << "\n";
	}
	
	// Get the alpha cuts for XNXF stuff
	TFile* f_xnxf = new TFile( xnxfcut_dir.Data() );
	if ( f_xnxf->IsOpen() ){
		cut_list_xnxf = GetXNXFCutArray( f_xnxf );
		f_xnxf->Close();
	}
	else{
		std::cout << "NO XNXF CUTS FOUND. Will carry on without them." << "\n";
	}

	// Open the root file for writing
	if ( PRINT_ROOT == 1 ){ out_root_file = new TFile( Form( "%s/posXXX.root", print_dir.Data() ), "RECREATE" ); }

	// Start timing
	stopwatch.Start();
}

void AnalyseTree::SlaveBegin(TTree* t)
{
	// The SlaveBegin() function is called after the Begin() function.
	// When running with PROOF SlaveBegin() is called on each slave server.
	// The tree argument is deprecated (on PROOF 0 is passed).

	TString option = GetOption();
}

Bool_t AnalyseTree::Process(Long64_t entry)
{
	// The Process() function is called for each entry in the tree (or possibly
	// keyed object in the case of PROOF) to be processed. The entry argument
	// specifies which entry in the currently loaded tree is to be processed.
	// When processing keyed objects with PROOF, the object is already loaded
	// and is available via the fObject pointer.
	//
	// This function should contain the \"body\" of the analysis. It can contain
	// simple or elaborate selection criteria, run algorithms on the data
	// of the event and typically fill histograms.
	//
	// The processing can be stopped by calling Abort().
	//
	// Use fStatus to set the return value of TTree::Process().
	//
	// The return value is currently not used.
	
	// Count the entries and update the clock
	processed_entries++;
	if ( processed_entries > num_entries*entry_frac  ){
		std::cout << "Processed: " << std::setw(3) << 100*entry_frac << "%" << " in " << Form( "%6.1f", stopwatch.RealTime() ) << " s" << "\n";
		stopwatch.Start(kFALSE);
		entry_frac += 0.1;
	}

	// Get branches
	/*
	b_TAC->GetEntry(entry);               b_ELUM->GetEntry(entry);
	b_EZERO->GetEntry(entry);             b_EnergyTimestamp->GetEntry(entry);
	b_RDTTimestamp->GetEntry(entry);      b_TACTimestamp->GetEntry(entry);
	b_ELUMTimestamp->GetEntry(entry);     b_EZEROTimestamp->GetEntry(entry);
	b_X->GetEntry(entry);                 
	b_TD_RDT_ELUM->GetEntry(entry);       b_XOLD->GetEntry(entry);
	b_Ex_CORRECTED->GetEntry(entry);
	b_ECAL->GetEntry(entry);
	*/
	//std::cout << "Get branches" << "\n";
	b_Energy->GetEntry(entry);
	b_XF->GetEntry(entry);                b_XN->GetEntry(entry);
	b_XFCAL->GetEntry(entry);             b_XNCAL->GetEntry(entry);
	//std::cout << "Get more branches" << "\n";
	b_Z->GetEntry(entry);
	b_RDT->GetEntry(entry);
	b_XCAL->GetEntry(entry);              
	b_ECRR->GetEntry(entry);              b_TD_RDT_E->GetEntry(entry);
	b_Ex->GetEntry(entry);                
	b_ThetaCM->GetEntry(entry);           b_DetID->GetEntry(entry);
	b_TD_RDT_E_CUTS->GetEntry(entry);     b_XCAL_CUTS->GetEntry(entry);

	// Work out if it is inside the cut(s)
	is_in_rdt = 0; is_in_rdt_si = 0;
	for ( Int_t i = 0; i < 4; i++ ){
		if ( is_in_rdt == 0 ){
			TCutG* cut = (TCutG*)cut_list->At(i);
			is_in_rdt += cut->IsInside( rdt[i+4], rdt[i] );
		}
		
		if ( found_si_cuts && is_in_rdt_si == 0 ){
			TCutG* cut_si = (TCutG*)cut_list_si->At(i);
			is_in_rdt_si += cut_si->IsInside( rdt[i+4], rdt[i] );
		}
	}
	
	// Create some doubles for use in the loop
	Double_t XNCAL, XFCAL, XNcal, XFcal;
	
	// *LOOP* OVER DETECTORS IN THE ARRAY
	for ( Int_t i = 0; i < 24; i++ ){

		// Calculate cut booleans
		is_in_used_det = ( det_array[ i % 6 ][ (Int_t)TMath::Floor( i/6 ) ] == 1 );
		// is_in_rdt goes here
		is_in_td = ( ( td_rdt_e[i][0] > td_rdt_e_cuts[i][0] && td_rdt_e[i][0] < td_rdt_e_cuts[i][1] ) || \
		             ( td_rdt_e[i][1] > td_rdt_e_cuts[i][0] && td_rdt_e[i][1] < td_rdt_e_cuts[i][1] ) || \
		             ( td_rdt_e[i][2] > td_rdt_e_cuts[i][0] && td_rdt_e[i][2] < td_rdt_e_cuts[i][1] ) || \
		             ( td_rdt_e[i][3] > td_rdt_e_cuts[i][0] && td_rdt_e[i][3] < td_rdt_e_cuts[i][1] ) );
		is_in_theta_min = ( thetaCM[i] >= THETA_MIN );
		is_in_theta_range = ( thetaCM[i] >= THETA_LB && thetaCM[i] <= THETA_UB );
		is_in_xcal = ( xcal[i] >= xcal_cuts[i][0] && xcal[i] <= xcal_cuts[i][1] );
		is_in_xcal_mid = ( xcal[i] <= xcal_cuts[i][2] || xcal[i] >= xcal_cuts[i][3] );
		
		// XNXF cut boolean
		TCutG* cut_xnxf = (TCutG*)cut_list_xnxf->At(i);
		if ( cut_xnxf != NULL ){
			is_in_xnxf_cut = cut_xnxf->IsInside( xf[i], xn[i] );
		}
		else{
			is_in_xnxf_cut = 0;
		}
		
		// CREATE HISTOGRAMS ------------------------------------------------------------------- //
		// Order of cuts should be as written above
		XNCAL = xfxneCorr[i][1]*xnCorr[i]*xn[i] + xfxneCorr[i][0];
		XFCAL = xfxneCorr[i][1]*xf[i] + xfxneCorr[i][0];
		XNcal = xnCorr[i]*xn[i];
		XFcal = xf[i];

		// *HIST* XN-XF
		if ( SW_XNXF[0] == 1 && ( DET_NUMBER == i || DET_NUMBER == -1 ) ){ 
			h_xnxf[i]->Fill( xf[i], xn[i] );
			
			if ( !TMath::IsNaN( xf[i] ) && !TMath::IsNaN( xn[i] ) ){
				if ( XNCAL <  0.5*e[i] && XFCAL >= 0.5*e[i] ){ h_xnxf_colour[i][0]->Fill( xf[i], xn[i] ); }
				if ( XNCAL >= 0.5*e[i] && XFCAL >= 0.5*e[i] ){ h_xnxf_colour[i][1]->Fill( xf[i], xn[i] ); }
				if ( XNCAL <  0.5*e[i] && XFCAL <  0.5*e[i] ){ h_xnxf_colour[i][2]->Fill( xf[i], xn[i] ); }
				if ( XNCAL >= 0.5*e[i] && XFCAL <  0.5*e[i] ){ h_xnxf_colour[i][3]->Fill( xf[i], xn[i] ); }
			}
			if ( is_in_xnxf_cut ){ p_xnxf[i]->Fill( xf[i], xn[i] ); }
			
			h_xnE[i]->Fill( xn[i], e[i] );
			h_xfE[i]->Fill( xf[i], e[i] );
			
			if ( XNCAL <  0.5*e[i] ){ h_xnE_colour[i][0]->Fill( xn[i], e[i] ); }
			if ( XNCAL >= 0.5*e[i] ){ h_xnE_colour[i][1]->Fill( xn[i], e[i] ); }
			
			if ( XFCAL <  0.5*e[i] ){ h_xfE_colour[i][0]->Fill( xf[i], e[i] ); }
			if ( XFCAL >= 0.5*e[i] ){ h_xfE_colour[i][1]->Fill( xf[i], e[i] ); }

			if ( !TMath::IsNaN( xf[i] ) && !TMath::IsNaN( xn[i] ) && !TMath::IsNaN( e[i] ) ){
				if ( XNcal != 0.0 && XNcal/XFcal <= XNXF_FRAC ){ h_xnxfE_colour[i][0]->Fill( xnCorr[i]*xn[i] + xf[i], e[i] ); }
				else if ( XFcal != 0.0 && XNcal/XFcal >= ( 1 - XNXF_FRAC ) ){ h_xnxfE_colour[i][1]->Fill( xnCorr[i]*xn[i] + xf[i], e[i] ); }
				else{ 
					h_xnxfE_colour[i][2]->Fill( xnCorr[i]*xn[i] + xf[i], e[i] );
					h_xnxfE[i]->Fill( xnCorr[i]*xn[i] + xf[i], e[i] );
					p_xnxfE[i]->Fill( xnCorr[i]*xn[i] + xf[i], e[i] );
				}
			}

			
		}

		// *HIST* xcal no cuts
		if ( is_in_used_det && is_in_rdt && is_in_td && is_in_theta_min ){
			if ( SW_XCAL[0] == 1 ){ h_xcal[i]->Fill( xcal[i] ); }
		}
			
		// *HIST* xcal with cuts
		if ( is_in_used_det && is_in_rdt && is_in_td && is_in_theta_min && is_in_xcal ){ 
			if ( SW_XCAL[0] == 1 ){ h_xcal_cut[i]->Fill( xcal[i] );  }
		}
		
		// Do singles cuts
		if ( is_in_used_det && is_in_theta_min && is_in_xcal && is_in_xcal_mid ){
			// *HIST* E v.s. z plot - singles
			if ( SW_EVZ_COMPARE[0] == 1 ){ h_evz_compare[0]->Fill( z[i], ecrr[i] ); }
			if ( SW_EVZ_SI[0] == 1 ){ h_evz_si[0]->Fill( z[i], ecrr[i] ); }
			
			// Add additional angle cut
			if ( is_in_theta_range ){
				if ( SW_EVZ_COMPARE[0] == 1 ){ h_evz_compare[1]->Fill( z[i], ecrr[i] ); }
			}
			
			// Implement the desired row numbers only
			if ( i % 6 == ROW_NUMBER || ROW_NUMBER == -1 ){
				// *HIST* Compare excitation spectra 1
				if ( SW_EX_COMPARE[0] == 1 ){ h_ex_compare1[i % 6]->Fill( Ex[i] ); }
			}
		}
			
		// Do full cuts (Mg)
		if ( is_in_used_det && is_in_rdt && is_in_td && is_in_theta_min && is_in_xcal && is_in_xcal_mid ){
			
			// *HIST* Full E v.s. z (Mg)
			if ( SW_EVZ[0] == 1 ){ h_evz->Fill( z[i], ecrr[i] ); }
			if ( SW_EVZ_SI[0] == 1 ){ h_evz_si[1]->Fill( z[i], ecrr[i] ); }
			if ( SW_EVZ_COMPARE[0] == 1 ){ h_evz_compare[2]->Fill( z[i], ecrr[i] ); }
			if ( SW_RDT_CUTS[0] == 1 ){ h_rdt_evz_mg[ (Int_t)TMath::Floor( i/6 ) ]->Fill( z[i], ecrr[i] ); }
			
			// Add additional angle cut
			if ( is_in_theta_range ){
				if ( SW_EVZ_COMPARE[0] == 1 ){ h_evz_compare[3]->Fill( z[i], ecrr[i] ); }
			}
			
			
			// *HIST* Full excitation plot (Mg)
			if ( SW_RDT_CUTS[0] == 1 ){ h_rdt_ex_mg[ (Int_t)TMath::Floor( i/6 ) ]->Fill( Ex[i] ); }
			if ( SW_EX_SI[0] == 1 ){ h_ex_si[0]->Fill( Ex[i] ); }
			if ( SW_EX[0] == 1 ){
				if ( ALL_ROWS == 1 ){ h_ex_full->Fill( Ex[i] ); }
				if ( ROW_BY_ROW == 1 ){ h_ex_rbr[ i % 6 ]->Fill( Ex[i] ); }
				if ( DET_BY_DET == 1 ){ h_ex_dbd[i]->Fill( Ex[i] ); }
			}
			
			// Implement the desired row numbers only
			if ( i % 6 == ROW_NUMBER || ROW_NUMBER == -1 ){
				// *HIST* Compare excitation spectra 2
				if ( SW_EX_COMPARE[0] == 1 ){ h_ex_compare2[i % 6]->Fill( Ex[i] ); }
			}
		}


		// Do full cuts (Si)
		if ( is_in_used_det && found_si_cuts && is_in_rdt_si && is_in_td && is_in_theta_min && is_in_xcal && is_in_xcal_mid ){
			// *HIST* Full E v.s. z (Si)
			if ( SW_EVZ_SI[0] == 1 ){ h_evz_si[2]->Fill( z[i], ecrr[i] ); }
			
			// *HIST* Full Ex (Si)
			if ( SW_EX_SI[0] == 1 ){ h_ex_si[1]->Fill( Ex[i] ); }
		}
		
	} // *LOOP* over detectors


	// *LOOP* OVER RECOIL DETECTORS
	for ( Int_t i = 0; i < 4; i++ ){
	
		// *HIST* Recoil detectors
		if( SW_RDT_CUTS[0] == 1 ){
			h_rdt_cuts[i]->Fill( rdt[i+4], rdt[i] );
		}
	}
	
	return kTRUE;
}

void AnalyseTree::SlaveTerminate()
{
	// The SlaveTerminate() function is called after all entries or objects
	// have been processed. When running with PROOF SlaveTerminate() is called
	// on each slave server.

}

void AnalyseTree::Terminate()
{
	// The Terminate() function is the last function to be called during
	// a query. It always runs on the client, it can be used to present
	// the results graphically or save the results to file.
	
	// Draw stuff
	if ( SW_EX_COMPARE[0] == 1 ){ HDrawExCompare(); }
	if ( SW_RDT_CUTS[0] == 1 ){ HDrawRDTCuts( fChain ); }
	if ( SW_EVZ_COMPARE[0] == 1 ){ HDrawEVZCompare(); }
	if ( SW_EVZ[0] == 1 ){ HDrawEVZ(); }
	if ( SW_EVZ_SI[0] == 1 ){ HDrawEVZSi(); }
	if ( SW_EX_SI[0] == 1 ){ HDrawExSi(); }
	if ( SW_EX[0] == 1 ){ HDrawEx(); }
	if ( SW_XNXF[0] == 1 ){ HDrawXNXF(); }
	if ( SW_XCAL[0] == 1 ){ HDrawXCAL(); }
	
	if ( PRINT_ROOT == 1 ){ if ( out_root_file->IsOpen() ){ out_root_file->Close(); } }
	
	stopwatch.Start(kFALSE);
}
