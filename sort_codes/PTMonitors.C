#define PTMonitors_cxx

#include "PTMonitors.h"
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
Bool_t qDrawGraphs = 0;
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

//Float_t x[24],z[24];
//Float_t xcal[24],ecal[24],xfcal[24],xncal[24],ecrr[24];
Int_t tacA[24];
Float_t z_array_pos[6] = {35.868,29.987,24.111,18.248,12.412,6.676};//in cm

// z offset
Int_t OFF_POSITION = 1;
Bool_t ALPHA_RUN = 1;
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
TString cutFileDir = "/home/ptmac/Documents/07-CERN-ISS-Mg/analysis/working/ALL-MgCuts3.root";

// NEW TREE STUFF
TTree* fin_tree;

typedef struct {
	// Calculated quantities
	Float_t x[24];
	Float_t z[24];
	Float_t xcal[24];
	Float_t ecal[24];
	Float_t xfcal[24];
	Float_t xncal[24];
	Float_t ecrr[24];
	Float_t Ex[24];
	Float_t thetaCM[24];
	Int_t detID[24];
	int td_rdt_e[24][4];
	int td_rdt_elum[32][4];
	int td_e_ebis[24];
	TCutG* cut[100];
} FIN;

FIN fin;

// TSELECTOR BEGIN FUNCTION -------------------------------------------------------------------- //
void PTMonitors::Begin(TTree *tree){
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
	TD_EBIS = new TH1F("TD_EBIS", "", 10001, -5000, 5000);
	TD_EBIS->GetYaxis()->SetTitle("# counts");
	TD_EBIS->GetXaxis()->SetTitle("Time Difference / 10^{-8} s");
	TD_EBIS->SetFillColor(5);

	// Time difference on the Energy-Recoil time
	TD_Recoil = new TH1F("TD_Recoil", "", 2001, -1000, 1000);
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

	// NEW TTREE STUFF
	if ( ALPHA_RUN == 0 ){
		outFile = new TFile( Form( "fin%i.root", OFF_POSITION ), "RECREATE");
	}
	else{
		outFile = new TFile( "finAlpha2.root", "RECREATE");
	}
	
	fin_tree = new TTree( "fin_tree", "Tree containing everything" );
	fin_tree->Branch("e",e,"Energy[100]/F");
	fin_tree->Branch("e_t",e_t,"EnergyTimestamp[100]/l");
	fin_tree->Branch("xf",xf,"XF[100]/F");
	fin_tree->Branch("xf_t",xf_t,"XFTimestamp[100]/l");
	fin_tree->Branch("xn",xn,"XN[100]/F");
  	fin_tree->Branch("xn_t",xn_t,"XNTimestamp[100]/l"); 
	fin_tree->Branch("rdt",rdt,"RDT[100]/F");
	fin_tree->Branch("rdt_t",rdt_t,"RDTTimestamp[100]/l"); 
	fin_tree->Branch("tac",tac,"TAC[100]/F");
	fin_tree->Branch("tac_t",tac_t,"TACTimestamp[100]/l"); 
	fin_tree->Branch("elum",elum,"ELUM[32]/F");
	fin_tree->Branch("elum_t",elum_t,"ELUMTimestamp[32]/l");
	fin_tree->Branch("ezero",ezero,"EZERO[10]/F");
	fin_tree->Branch("ezero_t",ezero_t,"EZEROTimestamp[10]/l");
	//fin_tree->Branch("ebis_t",ebis_t,"EBISTimestamp/l"); 
	
	fin_tree->Branch("x",fin.x,"X[24]/F");
	fin_tree->Branch("z",fin.z,"Z[24]/F");
	fin_tree->Branch("xcal",fin.xcal,"XCalibrated[24]/F");
	fin_tree->Branch("ecal",fin.ecal,"ECalibrated[24]/F");
	fin_tree->Branch("xfcal",fin.xfcal,"XFCalibrated[24]/F");
	fin_tree->Branch("xncal",fin.xncal,"XNCalibrated[24]/F");
	fin_tree->Branch("ecrr",fin.ecrr,"ECalibrated[24]/F");
	fin_tree->Branch("td_rdt_e",fin.td_rdt_e,"RDT-E_TD[24][4]/I");
	fin_tree->Branch("td_rdt_elum",fin.td_rdt_elum,"RDT-ELUM_TD[32][4]/I");
	fin_tree->Branch("td_e_ebis",fin.td_e_ebis,"E-EBIS_TD[24]/I");
	fin_tree->Branch("Ex",fin.Ex,"Ex[24]/F");
	fin_tree->Branch("thetaCM",fin.thetaCM,"ThetaCM[24]/F");
	fin_tree->Branch("detID",fin.detID,"DetID[24]/I");
	fin_tree->Branch("td_rdt_e_cuts",td_rdt_e_cuts,"TD-RDT-E-CUTS[24][2]/I");
	fin_tree->Branch("xcal_cuts",xcal_cuts,"XCAL-E-CUTS[24][2]/F");

	printf("======== number of cuts found : %d \n", numCut);
	StpWatch.Start();
}

// TSELECTOR SLAVEBEGIN FUNCTION --------------------------------------------------------------- //
void PTMonitors::SlaveBegin(TTree * /*tree*/){
	TString option = GetOption();
}

// TSELECTOR MAIN PROCESS ---------------------------------------------------------------------- //
Bool_t PTMonitors::Process(Long64_t entry){
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
		
		// RESET ALL QUANTITIES TO NaN
		for ( Int_t i = 0; i < 32; i++ ){
			if ( i < 24 ){
				fin.x[i] = TMath::QuietNaN();
		 		fin.z[i] = TMath::QuietNaN();
		 		fin.xcal[i] = TMath::QuietNaN();
		 		fin.ecal[i] = TMath::QuietNaN();
		 		fin.xfcal[i] = TMath::QuietNaN();
		 		fin.xncal[i] = TMath::QuietNaN();
		 		fin.ecrr[i] = TMath::QuietNaN();
				fin.Ex[i] = TMath::QuietNaN();
		 		fin.thetaCM[i] = TMath::QuietNaN();
		 		fin.detID[i] = TMath::QuietNaN();
				fin.td_e_ebis[i] = TMath::QuietNaN();
			}
			for ( Int_t j = 0; j < 4; j++ ){
				if ( i < 24 ){
					fin.td_rdt_e[i][j] = TMath::QuietNaN();
				}				
				fin.td_rdt_elum[i][j] = TMath::QuietNaN();
			}
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
		/* RECOIL-ELUM */
		// Calculate the elum-recoil time, by first populating arrays with junk
		for ( Int_t i = 0; i < 32; i++ ){
			for ( Int_t j = 0; j < 4; j++ ){
				if ( rdt_t[j] != 0 && elum_t[i] != 0 ){
					fin.td_rdt_elum[i][j]= (int)(rdt_t[j]-elum_t[i]);
				}
				else {
					fin.td_rdt_elum[i][j] = 10000;
				}
			}
		}
			
		/* ARRAY */
		for (Int_t i = 0; i < 24; i++) {
			// Calibrate each of the detectors
			fin.xfcal[i] = xf[i]*xfxneCorr[i][1]+xfxneCorr[i][0];
			fin.xncal[i] = xn[i]*xnCorr[i]*xfxneCorr[i][1]+xfxneCorr[i][0];
			fin.ecal[i] = e[i]/eCorr[i][0]+eCorr[i][1];
			fin.ecrr[i] = e[i]/eCorr[i][0]+eCorr[i][1];
		
			// Calculate the uncalibrated position on the strip
			if (xf[i]>0 || xn[i]>0 || !TMath::IsNaN(xf[i]) || !TMath::IsNaN(xn[i])) {
				fin.x[i] = 0.5*((xf[i]-xn[i]) / (xf[i]+xn[i]))+0.5;
			}
		
			// Calculate the calibrated position on the strip
			if ( fin.xfcal[i] > 0.5*e[i] ) {
				fin.xcal[i] = fin.xfcal[i]/e[i];
			}else if ( fin.xncal[i] >= 0.5*e[i] ) {
				fin.xcal[i] = 1.0 - fin.xncal[i]/e[i];
			}
		  
			// Calculate the exact position on the z axis
			fin.z[i] = 5.0*( fin.xcal[i] - 0.5 ) - z_off - z_array_pos[i%6];
			
			/* Fill the E-dE histograms if:
				* The position x (position on the strip) is between -1.1 and 1.1
				* The energy is greater than 100
				* One of xn or xf is greater than 0
			*/
			if ( fin.x[i] > -1.1 && fin.x[i] <1.1 && e[i] > 100 && ( xn[i] > 0 || xf[i] > 0 ) ){
				// Loop over the number of recoil detectors
	 			for ( Int_t ii = 0; ii < 4; ii++ ){
					EdE[ii]->Fill( rdt[ii+4], rdt[ii] );
				}
			}
		  
		} //Array loop
		/* TACs */
		for(Int_t i = 0; i < 4 ; i++){				// Loop over each side of array
			for(Int_t j = 0; j < 6; j++){			// Loop over each strip of side
		
				// Label the strip from 0 --> 23
				Int_t index = i*6+j;
				fin.detID[index] = index;
		
				//======== Ex calculation by Ryan 
				double y = fin.ecrr[index] + mass; // to give the KE + mass of proton;
				double Z = alpha * gamm * beta * fin.z[index] * 10.;
				double H = TMath::Sqrt(TMath::Power(gamm * beta,2) * (y*y - mass * mass) ) ;

				// Calculate the angle
				if( TMath::Abs(Z) < H ) {
					// Use Newton's method to solve 0 ==  H * sin(phi) - G * tan(phi) - Z = f(phi) 
					double tolerance = 0.001;	// Desired precision
			 	 	double phi = 0; 			// Initial phi = 0 -> ensure the solution has f'(phi) > 0
					double nPhi = 0; 			// New phi

					int iter = 0;				// Number of iterations
						
					// Now calculate the angle using Newton-Raphson process
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
						// Found correct value of phi - now calculate everything else
						double K = H * TMath::Sin(phi);
						double x = TMath::ACos( mass / ( y * gamm - K));
						double momt = mass * TMath::Tan( x ); // momentum of particle b or B in CM frame
						double EB = TMath::Sqrt(mass*mass + Et*Et - 2*Et*TMath::Sqrt(momt * momt + mass * mass));
						fin.Ex[index] = EB - massB;
						
						double hahaha1 = gamm* TMath::Sqrt(mass * mass + momt * momt) - y;
						double hahaha2 = gamm* beta * momt;
						fin.thetaCM[index] = TMath::ACos(hahaha1/hahaha2) * TMath::RadToDeg();
				 
					}
					else{
						fin.Ex[index] = TMath::QuietNaN();
						fin.thetaCM[index] = TMath::QuietNaN();
					}	
				}
				else{
					fin.Ex[index] = TMath::QuietNaN();
					fin.thetaCM[index] = TMath::QuietNaN();
				}
				
				// Calculate the EBIS time - the array time and populate a histogram
				fin.td_e_ebis[index] = 10000;
				if ( ebis_t != 0 && e_t[index] != 0 ){
					fin.td_e_ebis[index] = (int)(e_t[index] - ebis_t);
				}
				TD_EBIS->Fill( fin.td_e_ebis[index] );
				

				// Calculate the recoil time stuff, by first populating arrays with junk
				for ( Int_t ii = 0; ii < 4; ii++ ){
						fin.td_rdt_e[index][ii] = 10000;
				}
				for ( Int_t kk = 0; kk < 4; kk++ ){
					if ( rdt_t[kk] != 0 && e_t[index] != 0 ){
						fin.td_rdt_e[index][kk]= (int)(rdt_t[kk]-e_t[index]);
					}
					TD_Recoil->Fill( fin.td_rdt_e[index][kk] );
				}

				// Now look at cuts for gated spectra
				if( isCutFileOpen){
					for( int k = 0 ; k < numCut; k++ ){
						fin.cut[k] = (TCutG *)cutList->At(k) ;
						if( fin.cut[k]->IsInside(rdt[k+4], rdt[k]) ) { //CRH
							for (Int_t kk = 0; kk < 4; kk++) { 
								if(-30 < fin.td_rdt_e[index][kk] && fin.td_rdt_e[index][kk] < 30) {
									EVZ->Fill( fin.z[index], fin.ecrr[index] );
									EXE->Fill(fin.Ex[index] );
									EXE_Row[index % 6]->Fill( fin.Ex[index] );
									XN_XF[index]->Fill( xn[index], xf[index] );
								}
							}
						}
					}
				}
			} // Strip loop
		} // Side loop

	// FILL THE NEW TTree BASED ON CALCULATIONS
	fin_tree->Fill();
	
	} // Processed entries
	return kTRUE;
}

// TSELECTOR SLAVE TERMINATE FUNCTION ---------------------------------------------------------- //
void PTMonitors::SlaveTerminate(){

}

// TSELECTOR TERMINATE FUNCTION ---------------------------------------------------------------- //
void PTMonitors::Terminate()
{
	// Write the cuts
	for ( int i = 0; i < 100; i++ ){
		if ( fin.cut[i] != NULL ){
			fin.cut[i]->Write();
		}
	}

	// Write the TTree
	fin_tree->Write();

	// Close the file
	if ( outFile != NULL ){ outFile->Close(); }
	
	// Print out some stuff
	if (ProcessedEntries>=NUMSORT){
		printf("Sorted only %llu\n",NUMSORT);
	}
	StpWatch.Start(kFALSE);
}
