#define PTMonitorsOLD_cxx

#include "PTMonitorsOLD.h"
#include <TH2.h>
#include <TH1.h>
#include <TStyle.h>
#include <TCutG.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TString.h>
#include <TSystem.h>
#include <TObjArray.h>
#include <TMath.h>
#include <TFile.h>

// SWITCHES FOR POST-PROCESSING
Bool_t qDrawGraphs = 1;
Bool_t qPrintGraphs = 0;
Bool_t qWriteData = 0;

// REDUCED TSelector CODE ---------------------------------------------------------------------- //
#define NUMPRINT 20 //>0
ULong64_t NUMSORT=100000000;
ULong64_t NumEntries = 0;
ULong64_t ProcessedEntries = 0;
Float_t Frac = 0.1; //Progress bar
TStopwatch StpWatch;

Int_t n=1;

TString cutName("cut1");
TCutG* cutG; //!
TObjArray * cutList;
TString cutTag;
Bool_t isCutFileOpen;
int numCut;
vector<int> countFromCut;

// Times in secondsf
Float_t timeZero=0;
Float_t timeCurrent=0;
Float_t timeRef=0;

Float_t x[24],z[24];
Float_t xcal[24],ecal[24],xfcal[24],xncal[24],ecrr[24],ezero[10];
Int_t tacA[24];
Float_t z_array_pos[6] = {35.868,29.987,24.111,18.248,12.412,6.676};//in cm

// z offset
Int_t OFF_POSITION = 1;
Float_t z_off;

Float_t xnCorr[24] = {0.907342,0.907342,0.976727,0.914866,1.021736,
		      0.887032,0.923250,0.953968,1.020180,0.918340,
		      0.983084,0.983084,0.997550,0.985319,0.959048,
		      1.008677,0.959601,1.066846,0.927771,0.985274,
		      0.921273,0.976498,1.062241,1.079507};
Float_t xfxneCorr[24][2] = {{29.091896,0.919262},{-0.744352,0.989133},{5.332432,1.046711},
			    {4.770114,1.073863},{-4.352881,0.901518},{-8.543459,0.995114},
			    {4.678705,1.015215},{3.955090,0.972769},{5.163730,0.998306},
			    {3.863314,  0.989275},{2.298429,  0.916884},{-17.435897,  0.897436},
			    {8.143049,  0.571533},{5.428828,  0.927071},{4.554876,  0.960028},
			    {4.423083,  0.967342},{1.436683,  1.026855},
			    {0.747782,  0.912706},
			    {6.048360, 0.914865},
			    {2.104460,  0.962689},
			    {1.011006,  1.034467},
			    {15.249334,  0.887257},
			    {14.071915,  1.095258},
			    {-2.256993,  0.896878}};
Float_t eCorr[24][2] = {{256.060637	,0.021569},
			{253.083810	,0.010404},
			{275.757609	,-0.012115},
			{266.830570	,0.028129},
			{247.134021	,0.013641},
			{244.161153	,0.002046},
			{263.857355	,0.042191},
			{250.108256	,-0.001003},
			{262.017938	,0.018393},
			{256.060637	,0.021569},
			{238.219726	,0.005357},
			{ 1.000000	,0.000000},
			{0	,0},
			{248.283604	,-0.026163},
			{242.321161	,-0.024002},
			{250.108256	,-0.001003},
			{262.017938	,-0.006414},
			{257.914882	,0.020954},
			{250.108256	,0.024985},
			{259.038694	,0.007406},
			{266.830570	,0.028129},
			{250.108256	,0.024985},
			{292.477670	,0.015062},
			{239.341772	-0.009266}};


Float_t exCorr[6] = { 938.272,  // mass of proton [MeV/c^2]
	                   1,        // charge of proton
	                   27954.0982, // cm frame total energy
	                   26996.5929, // mass of recoil
	                   0.132178, // beta to CM frame
	                   2.5}; // Bfield [T]
double a = 11.5 ; // perpendicular distance of detector to axis [mm]
double Ex, thetaCM;

double alpha, Et, beta, gamm, G, massB, mass; //variables for Ex calculation


Float_t tempTime=-1000;
Long64_t tempTimeLong=10001;

// HISTOGRAMS
TH2F* EVZ;			// Gated energy v.s. position
TH1F* EXE;			// Gated excitation spectrum
TH2F* EdE[4];		// Gated recoil detector E-dE plots
TH1F* TD_EBIS;		// Time difference on the EBIS-Energy time
TH1F* TD_Recoil;	// Time difference on the Energy-Recoil time
TH1F* EXE_Row[6];	// Gated excitation spectrum on the recoils.
TH2F* XN_XF[24];	// XN v.s. XF plots for each detector

// CANVASES
TCanvas *cEVZ, *cEXE;
TCanvas *cRecoilEdE[4];
TCanvas *cTD_EBIS, *cTD_Recoil;
TCanvas *cEXE_Row[6];
TCanvas *cXN_XF;

// OUTPUT FILE
TFile* outFile;

// CUTS FILE
TString cutFileDir = "/home/ptmac/Documents/07-CERN-ISS-Mg/analysis/working/ALL-MgCuts.root";

// NEW TREE STUFF
TTree* fin_tree;

// TSELECTOR BEGIN FUNCTION -------------------------------------------------------------------- //
void PTMonitorsOLD::Begin(TTree *tree){
	// Define offset (array position - offset position = 70mm???)
	if ( OFF_POSITION == 0 ){
		z_off = 4.9765;
	}
	else if ( OFF_POSITION == 1 ){
		z_off = 9.498;
	}
	else if ( OFF_POSITION == 2 ){
		z_off = 6.50;
	}
	Printf( "Z OFFSET = %f;\t ARRAY POSITION = %i", z_off, OFF_POSITION );


	TString option = GetOption();
	NumEntries = tree->GetEntries();

	//Get any cuts;
	TFile * fCut = new TFile( cutFileDir.Data() );			// open file
	isCutFileOpen = fCut->IsOpen(); 
	numCut = 0 ;
	if( isCutFileOpen ){
		cutList = (TObjArray *) fCut->FindObjectAny("cutList");
		numCut = cutList->GetEntries();
		printf("=========== found %d cutG in %s \n", numCut, fCut->GetName());
	
		cutG = new TCutG();
		for(int i = 0; i < numCut ; i++){
			printf(" cut name : %s , VarX: %s, VarY: %s, numPoints: %d \n",
				cutList->At(i)->GetName(),
		 		((TCutG*)cutList->At(i))->GetVarX(),
		 		((TCutG*)cutList->At(i))->GetVarY(),
		 		((TCutG*)cutList->At(i))->GetN()
	 		);
			countFromCut.push_back(0);
		}
	}

	
	alpha = 299.792458 * exCorr[5] * exCorr[1] / TMath::TwoPi() / 1000; //MeV/mm
	beta = exCorr[4];
	gamm = 1./TMath::Sqrt(1-beta*beta);
	G = alpha * gamm * beta * a;
	massB = exCorr[3];
	mass = exCorr[0];
	Et = exCorr[2];

	// SHARPY'S GRAPHS
	// Make a new TStyle
	TStyle *sharpyStyle = new TStyle("sharpyStyle","David Sharp's Style");
	sharpyStyle->SetOptStat(0);
	sharpyStyle->SetCanvasBorderMode(0);
	sharpyStyle->SetPadBorderMode(0);
	sharpyStyle->SetPadColor(0);
	sharpyStyle->SetCanvasColor(0);
	sharpyStyle->SetStatColor(0);
	sharpyStyle->SetTitleFillColor(0);
	sharpyStyle->SetTitleBorderSize(0);
	sharpyStyle->SetTitleAlign(23);
	sharpyStyle->SetTitleX(0.5);
	sharpyStyle->SetMarkerStyle(7);
	sharpyStyle->cd();
	
	// DEFINE HISTOGRAMS AND SET OPTIONS
	// Gated energy v.s. position
	EVZ = new TH2F("EVZ", "",700, -50, -5, 750 , 0 , 10);
	EVZ->GetXaxis()->SetTitle("z (cm)");
	EVZ->GetYaxis()->SetTitle("E (MeV)");
	
	// Gated excitation spectrum
	EXE = new TH1F("EXE", "", 400, -1, 8 );
	EXE->GetYaxis()->SetTitle("Counts");
	EXE->GetXaxis()->SetTitle("E (MeV)");
	EXE->SetFillColor(5);
	
	// Time difference on the EBIS-Energy time
	TD_EBIS = new TH1F("TD_EBIS", "", 100000, -1000, 1000);
	TD_EBIS->GetYaxis()->SetTitle("# counts");
	TD_EBIS->GetXaxis()->SetTitle("Time Difference / 10^{-8} s");
	TD_EBIS->SetFillColor(5);

	// Time difference on the Energy-Recoil time
	TD_Recoil = new TH1F("TD_Recoil", "", 100000, -1000, 1000);
	TD_Recoil->GetYaxis()->SetTitle("# counts");
	TD_Recoil->GetXaxis()->SetTitle("Time Difference / 10^{-8} s");
	TD_Recoil->SetFillColor(5);
	
	// Gated recoil detector E-dE plots	
	for ( Int_t ii = 0; ii < 4; ii++ ){
		EdE[ii] = new TH2F( Form("EdE%d",ii ), "", 1000, 0, 10000, 1000, 0, 4000 );
		EdE[ii]->SetTitle( Form( "Recoil %d", ii ) );
	}
	
	// Gated excitation spectrum on the recoils.
	for ( Int_t ii = 0; ii < 6; ii++ ){
		EXE_Row[ii] = new TH1F( Form( "EXE_Row%i", ii ), "", 450, -1, 8 );
		EXE_Row[ii]->GetYaxis()->SetTitle("Counts");
		EXE_Row[ii]->GetXaxis()->SetTitle("E (MeV)");
		EXE_Row[ii]->SetFillColor(5);
	}

	// XN v.s. XF plots for each detector
	for ( Int_t ii = 0; ii < 24; ii++ ){
		XN_XF[ii] = new TH2F( Form( "XN_XF: Row %i, Side %i", ii % 6, (int)TMath::Floor(ii/6) ), "", 5101, -100, 5000, 5101, -100, 5000 );
		XN_XF[ii]->GetYaxis()->SetTitle("XN");
		XN_XF[ii]->GetXaxis()->SetTitle("XF");
	}

	printf("======== number of cuts found : %d \n", numCut);
	StpWatch.Start();
}

// TSELECTOR SLAVEBEGIN FUNCTION --------------------------------------------------------------- //
void PTMonitorsOLD::SlaveBegin(TTree * /*tree*/){
	TString option = GetOption();
}

// TSELECTOR MAIN PROCESS ---------------------------------------------------------------------- //
Bool_t PTMonitorsOLD::Process(Long64_t entry){
	// Increment number of processed entries
	ProcessedEntries++;
	
	// Print out the progress of the sort
	if (ProcessedEntries<NUMSORT) {
		if (ProcessedEntries>NumEntries*Frac-1) {
			printf(" %3.0f%% (%llu/%llu k) processed in %6.1f seconds\n",
				Frac*100,ProcessedEntries/1000,NumEntries/1000,StpWatch.RealTime()
			);
			StpWatch.Start(kFALSE);
			Frac+=0.1;
		}
	
		// Get the entries from the defined TTree (populates each of the leaves for processing)
		b_Energy->GetEntry(entry);
		b_XF->GetEntry(entry);
		b_XN->GetEntry(entry);
		b_RDT->GetEntry(entry);
		b_TAC->GetEntry(entry);
		b_ELUM->GetEntry(entry);
		b_EZERO->GetEntry(entry);
		b_EnergyTimestamp->GetEntry(entry);
		b_RDTTimestamp->GetEntry(entry);
		b_TACTimestamp->GetEntry(entry);
		b_ELUMTimestamp->GetEntry(entry);
		b_EZEROTimestamp->GetEntry(entry);
		b_EBISTimestamp->GetEntry(entry);

		// DO CALCULATIONS
		/* ARRAY */
		for (Int_t i = 0; i < 24; i++) {
			// Calibrate each of the detectors
			xfcal[i] = xf[i]*xfxneCorr[i][1]+xfxneCorr[i][0];
			xncal[i] = xn[i]*xnCorr[i]*xfxneCorr[i][1]+xfxneCorr[i][0];
			ecal[i] = e[i]/eCorr[i][0]+eCorr[i][1];
			ecrr[i] = e[i]/eCorr[i][0]+eCorr[i][1];
		
			// Calculate the uncalibrated position on the strip
			if (xf[i]>0 || xn[i]>0 || !TMath::IsNaN(xf[i]) || !TMath::IsNaN(xn[i])) {
				x[i] = 0.5*((xf[i]-xn[i]) / (xf[i]+xn[i]))+0.5;
			}
		
			// Calculate the calibrated position on the strip
			if (xfcal[i]>0.5*e[i]) {
				xcal[i] = xfcal[i]/e[i];
			}else if (xncal[i]>=0.5*e[i]) {
				xcal[i] = 1.0 - xncal[i]/e[i];
			}
		  
			// Calculate the exact position on the z axis
			z[i] = 5.0*(xcal[i]-0.5) - z_off - z_array_pos[i%6];
			
			/* Fill the E-dE histograms if:
				* The position x (position on the strip) is between -1.1 and 1.1
				* The energy is greater than 100
				* One of xn or xf is greater than 0
			*/
			if ( x[i] > -1.1 && x[i] <1.1 && e[i] > 100 && ( xn[i] > 0 || xf[i] > 0 ) ){
				// Loop over the number of recoil detectors
	 			for ( Int_t ii = 0; ii < 4; ii++ ){
					EdE[ii]->Fill(rdt[ii+4],rdt[ii]);
				}
			}
		  
		} //Array loop
		/* TACs */
		for(Int_t i = 0; i < 4 ; i++){				// Loop over each side of array
			for(Int_t j = 0; j < 6; j++){			// Loop over each strip of side
		
				// Label the strip from 0 --> 23
				int detID = i*6+j;
		
				//======== Ex calculation by Ryan 
				double y = ecrr[detID] + mass; // to give the KE + mass of proton;
				double Z = alpha * gamm * beta * z[detID] * 10.;
				double H = TMath::Sqrt(TMath::Power(gamm * beta,2) * (y*y - mass * mass) ) ;

				// Calculate the angle
				if( TMath::Abs(Z) < H ) {
					// Use Newton's method to solve 0 ==  H * sin(phi) - G * tan(phi) - Z = f(phi) 
					double tolerance = 0.001;	// Desired precision
			 	 	double phi = 0; 			// Initial phi = 0 -> ensure the solution has f'(phi) > 0
					double nPhi = 0; 			// New phi

					int iter = 0;				// Number of iterations
						
					// Now calculate the angle
					do{
						phi = nPhi;
						nPhi = phi - (H * TMath::Sin(phi) - G * TMath::Tan(phi) - Z) / (H * TMath::Cos(phi) - G /TMath::Power( TMath::Cos(phi), 2));
						iter ++;
						if( iter > 10 || TMath::Abs(nPhi) > TMath::PiOver2()) break;
					} while( TMath::Abs(phi - nPhi ) > tolerance);
					phi = nPhi;

					// Check f'(phi) > 0
					double Df = H * TMath::Cos(phi) - G / TMath::Power( TMath::Cos(phi),2);
					if( Df > 0 && TMath::Abs(phi) < TMath::PiOver2()  ){
						double K = H * TMath::Sin(phi);
						double x = TMath::ACos( mass / ( y * gamm - K));
						double momt = mass * TMath::Tan(x); // momentum of particle b or B in CM frame
						double EB = TMath::Sqrt(mass*mass + Et*Et - 2*Et*TMath::Sqrt(momt * momt + mass * mass));
						Ex = EB - massB;
						
						double hahaha1 = gamm* TMath::Sqrt(mass * mass + momt * momt) - y;
						double hahaha2 = gamm* beta * momt;
						thetaCM = TMath::ACos(hahaha1/hahaha2) * TMath::RadToDeg();
				 
					}
					else{
						Ex = TMath::QuietNaN();
						thetaCM = TMath::QuietNaN();
					}	
				}
				else{
					Ex = TMath::QuietNaN();
					thetaCM = TMath::QuietNaN();
				}
				
				// Calculate the EBIS time - the array time and populate a histogram
				Long64_t X = (Long64_t)e_t[detID] - (Long64_t)ebis_t;
				TD_EBIS->Fill( X );
				

				// Calculate the recoil time stuff
				for (Int_t kk = 0; kk < 4; kk++ ){
					TD_Recoil->Fill( (Long64_t)(rdt_t[kk]-e_t[detID]) );
				}

				// Now look at cuts for gated spectra
				if( isCutFileOpen){
					for( int k = 0 ; k < numCut; k++ ){
						cutG = (TCutG *)cutList->At(k) ;
						if( cutG->IsInside(rdt[k+4], rdt[k]) ) { //CRH
							for (Int_t kk = 0; kk < 4; kk++) { 
								tacA[detID]= (int)(rdt_t[kk]-e_t[detID]);
								if(-30 < tacA[detID] && tacA[detID] < 30) {
									EVZ->Fill(z[detID],ecrr[detID]);
									EXE->Fill(Ex);
									EXE_Row[detID % 6]->Fill(Ex);
									XN_XF[detID]->Fill(xn[detID],xf[detID]);
								}
							}
						}
					}
				}
			} // Strip loop
		} // Side loop
	} // Processed entries
	return kTRUE;
}

// TSELECTOR SLAVE TERMINATE FUNCTION ---------------------------------------------------------- //
void PTMonitorsOLD::SlaveTerminate(){

}

// TSELECTOR TERMINATE FUNCTION ---------------------------------------------------------------- //
void PTMonitorsOLD::Terminate()
{
	// PLOT SHARPY'S GRAPHS AND WRITE TO FILE
	/*
	if ( qWriteData == 1 ){ outFile = new TFile("fin.root", "RECREATE"); }
	*/
	// E vs.z
	if ( qDrawGraphs == 1 ){
		cEVZ = new TCanvas("cEVZ","E v.s. z");
		EVZ->Draw("scat colz same");
		if ( qPrintGraphs == 1){ cEVZ->Print("EVZ.png"); }
	}
	if ( qWriteData == 1 ){ EVZ->Write(); }
	/*
	// Plot excitation energy
	if ( qDrawGraphs == 1 ){ 
		cEXE = new TCanvas( "cEXE","Excitation energy spectrum", 1080, 810 ); 
		EXE->Draw();
		if (qPrintGraphs == 1){ cEXE->Print("EXE.png"); }
	}
	
	if ( qWriteData == 1 ){ EXE->Write(); }


	// Plot the E-dE plots for making cuts
	for ( Int_t i = 0; i < 4; i++ ){
		if ( qDrawGraphs == 1 ){
			cRecoilEdE[i] = new TCanvas( Form("cRecoilEdE_%i", i), "E-dE plots", 540, 405 );
			cRecoilEdE[i]->ToggleToolBar();
			EdE[i]->Draw("box colz");
			if( isCutFileOpen ) {
				cutG = (TCutG *)cutList->At(i);
				cutG->Draw("same");
			}
		}
		if ( qWriteData == 1 ){ EdE[i]->Write( Form("EdE-%i",i) ); }
	}
	
	if ( qDrawGraphs == 1 ){
		cTD_EBIS = new TCanvas( "cTD_EBIS", "EBIS time difference with energy", 1080, 816 );
		TD_EBIS->Draw();
	}

	if ( qDrawGraphs == 1 ){
		cTD_Recoil = new TCanvas( "cTD_Recoil", "Recoil detector time difference with energy", 1080, 816 );
		TD_Recoil->Draw();
	}
	
	for ( Int_t i = 0; i < 6; i++ ){	
		if ( qDrawGraphs == 1 ){
			cEXE_Row[i] = new TCanvas( Form( "cEXE_Row%i", i ), Form( "Si detector E v.s. z: row %i", i ), 1080, 816 );
			EXE_Row[i]->Draw();
		}
		if ( qPrintGraphs == 1){ cEXE_Row[i]->Print( Form("cEXE_Row%i_Pos%i.png", i, OFF_POSITION ) ); }
		if ( qWriteData == 1 ){
			writespe( EXE_Row[i]->GetName(), Form( "Ex_Row%i_Pos%i", i, OFF_POSITION ) );
		}
	}
	
	
	if ( qDrawGraphs == 1 ){
		cXN_XF = new TCanvas( "cXN_XF", "XN v.s. XF plot", 1800, 900 );
		cXN_XF->Divide(6,4);
		for ( Int_t i = 0; i < 24; i++ ){
			cXN_XF->cd(i+1);
			XN_XF[i]->Draw();
		}
	}
	*/
	
	// Close the file
	if ( outFile != NULL ){ outFile->Close(); }
	
	
	// Print out some stuff
	if (ProcessedEntries>=NUMSORT){
		printf("Sorted only %llu\n",NUMSORT);
	}
	StpWatch.Start(kFALSE);
}
