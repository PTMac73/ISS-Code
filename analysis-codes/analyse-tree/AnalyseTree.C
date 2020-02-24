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


/* TODO
 * Header file to contain all of the histograms
 * Make a function to make new histograms and format them all
 * Header file to contain the settings and switches
 
 PLOTS
 * EVZ
 * RDT
 * Anything in Plotter



*/


// VARIABLES


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
	
	// Create histograms
	if ( SW_EX_COMPARE[0] == 1 ){ HCreateExCompare(); }
	if ( SW_RDT_CUTS[0] == 1 ){ HCreateRDTCuts(); }
	if ( SW_EVZ_COMPARE[0] == 1 ){ HCreateEVZCompare(); }
	if ( SW_EVZ[0] == 1 ){ HCreateEVZ(); }
	
	// Get the number of entries
	num_entries = t->GetEntries();
	TString option = GetOption();
	
	// Get the cuts
	TFile *f = new TFile( "/home/ptmac/Documents/07-CERN-ISS-Mg/analysis/working/ALL-MgCuts3.root" );
	if ( f->IsOpen() ){
		cut_list = (TObjArray*)f->FindObjectAny("cutList");
	}
	else{
		std::cout << "NO CUTS" << "\n";
		std::exit(1);
	}

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
	b_Energy->GetEntry(entry);            b_XF->GetEntry(entry);
	b_XN->GetEntry(entry);
	b_TAC->GetEntry(entry);               b_ELUM->GetEntry(entry);
	b_EZERO->GetEntry(entry);             b_EnergyTimestamp->GetEntry(entry);
	b_RDTTimestamp->GetEntry(entry);      b_TACTimestamp->GetEntry(entry);
	b_ELUMTimestamp->GetEntry(entry);     b_EZEROTimestamp->GetEntry(entry);
	b_X->GetEntry(entry);                 
	
	b_TD_RDT_ELUM->GetEntry(entry);       b_XOLD->GetEntry(entry);
	b_Ex_CORRECTED->GetEntry(entry);
	b_XFCAL->GetEntry(entry);             b_XNCAL->GetEntry(entry);
	b_ECAL->GetEntry(entry);
	*/
	b_Z->GetEntry(entry);
	b_RDT->GetEntry(entry);
	b_XCAL->GetEntry(entry);              
	b_ECRR->GetEntry(entry);              b_TD_RDT_E->GetEntry(entry);
	b_Ex->GetEntry(entry);                
	b_ThetaCM->GetEntry(entry);           b_DetID->GetEntry(entry);
	b_TD_RDT_E_CUTS->GetEntry(entry);     b_XCAL_CUTS->GetEntry(entry);


	// Work out if it is inside the cut(s)
	is_in_rdt = 0;
	for ( Int_t i = 0; i < 4; i++ ){
		if ( is_in_rdt == 0 ){
			TCutG* cut = (TCutG*)cut_list->At(i);
			is_in_rdt += cut->IsInside( rdt[i+4], rdt[i] );
		}
	}
	
	// *LOOP* OVER DETECTORS IN THE ARRAY
	for ( Int_t i = 0; i < 24; i++ ){
	
		// Calculate cut booleans
		is_in_used_det = ( det_array[ i % 6 ][ (Int_t)TMath::Floor( i/6 ) ] == 1 );
		is_in_td = ( ( td_rdt_e[i][0] > td_rdt_e_cuts[i][0] && td_rdt_e[i][0] < td_rdt_e_cuts[i][1] ) || \
		             ( td_rdt_e[i][1] > td_rdt_e_cuts[i][0] && td_rdt_e[i][1] < td_rdt_e_cuts[i][1] ) || \
		             ( td_rdt_e[i][2] > td_rdt_e_cuts[i][0] && td_rdt_e[i][2] < td_rdt_e_cuts[i][1] ) || \
		             ( td_rdt_e[i][3] > td_rdt_e_cuts[i][0] && td_rdt_e[i][3] < td_rdt_e_cuts[i][1] ) );
		is_in_xcal = ( xcal[i] >= xcal_cuts[i][0] && xcal[i] <= xcal_cuts[i][1] && ( xcal[i] <= xcal_cuts[i][2] || xcal[i] >= xcal_cuts[i][3] ) );
		is_in_theta_min = ( thetaCM[i] >= THETA_MIN );
		is_in_theta_range = ( thetaCM[i] >= THETA_LB && thetaCM[i] <= THETA_UB );
		
		// Do the required cuts for a singles spectrum, but not the upper angle
		if ( is_in_used_det && is_in_xcal && is_in_theta_min ){
			
			// *HIST* Compare E v.s. z plots no angle cut - singles 1
			if ( SW_EVZ_COMPARE[0] == 1 ){ h_evz_compare[0]->Fill( z[i], ecrr[i] ); }
			
			
			// Do full cuts
			if ( is_in_rdt && is_in_td ){
			
				// *HIST* Full E v.s. z
				if ( SW_EVZ[0] == 1 ){ h_evz->Fill( z[i], ecrr[i] ); }
				
				// *HIST* Compare E v.s. z plots no angle cut - clean 1
				if ( SW_EVZ_COMPARE[0] == 1 ){ h_evz_compare[2]->Fill( z[i], ecrr[i] ); }
				
			}	// If in the full cuts
			
			
			
			// Add lower and upper angular range cut
			if ( is_in_theta_range ){
			
				// *HIST* Compare E v.s. z plots no additional angle cut 2
				if ( SW_EVZ_COMPARE[0] == 1 ){ h_evz_compare[1]->Fill( z[i], ecrr[i] ); }
				
				// Add full cuts
				if ( is_in_rdt && is_in_td ){
				
					// *HIST* Compare E v.s. z plots no angle cut - clean 1
					if ( SW_EVZ_COMPARE[0] == 1 ){ h_evz_compare[3]->Fill( z[i], ecrr[i] ); }
				}
					
				// Implement the desired row numbers only
				if ( i % 6 == ROW_NUMBER || ROW_NUMBER == -1 ){
					// *HIST* Compare excitation spectra 1
					if ( SW_EX_COMPARE[0] == 1 ){ h_ex_compare1[i % 6]->Fill( Ex[i] ); }
						
					// Do the required cuts for a clean spectrum
					if ( is_in_rdt && is_in_td ){
						// *HIST* Compare excitation spectra 2
						if ( SW_EX_COMPARE[0] == 1 ){ h_ex_compare2[i % 6]->Fill( Ex[i] ); }
						
					} // If in the full spectrum
					
				}	// If in the right row
				
			}	// If in range of theta
			
		}	// If in the singles spectrum
		
	}	// Loop over detectors (i)
	
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
	if ( SW_RDT_CUTS[0] == 1 ){ HDrawRDTCuts(); }
	if ( SW_EVZ_COMPARE[0] == 1 ){ HDrawEVZCompare(); }
	if ( SW_EVZ[0] == 1 ){ HDrawEVZ(); }
	
	
	stopwatch.Start(kFALSE);
}
