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
Int_t numCut;
vector<Int_t> countFromCut;

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
Bool_t ALPHA_RUN = 0;
Float_t z_off;
//Float_t xcal_cuts[24][4];
Float_t thetaCM_lims[9];

// RYAN'S CORRECTION PARAMETERS
/*
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
*/
/* Ryan's code formatted better
Float_t xnCorr[24] = {
	0.907342,	// 00
	0.907342,	// 01
	0.976727,	// 02
	0.914866,	// 03
	1.021736,	// 04
	0.887032,	// 05
	0.923250,	// 06
	0.953968,	// 07
	1.020180,	// 08
	0.918340,	// 09
	0.983084,	// 10
	0.983084,	// 11
	0.997550,	// 12
	0.985319,	// 13
	0.959048,	// 14
	1.008677,	// 15
	0.959601,	// 16
	1.066846,	// 17
	0.927771,	// 18
	0.985274,	// 19
	0.921273,	// 20
	0.976498,	// 22
	1.062241,	// 22
	1.079507	// 23
};
Float_t xfxneCorr[24][2] = {
	{ 29.0918960, 0.919262 },	// 00
	{ -0.7443520, 0.989133 },	// 01
	{  5.3324320, 1.046711 },	// 02
	{  4.7701140, 1.073863 },	// 03
	{ -4.3528810, 0.901518 },	// 04
	{ -8.5434590, 0.995114 },	// 05
	{  4.6787050, 1.015215 },	// 06
	{  3.9550900, 0.972769 },	// 07
	{  5.1637300, 0.998306 },	// 08
	{  3.8633140, 0.989275 },	// 09
	{  2.2984290, 0.916884 },	// 10
	{ -17.435897, 0.897436 },	// 11
	{  8.1430490, 0.571533 },	// 12
	{  5.4288280, 0.927071 },	// 13
	{  4.5548760, 0.960028 },	// 14
	{  4.4230830, 0.967342 },	// 15
	{  1.4366830, 1.026855 },	// 16
	{  0.7477820, 0.912706 },	// 17
	{  6.0483600, 0.914865 },	// 18
	{  2.1044600, 0.962689 },	// 19
	{  1.0110060, 1.034467 },	// 20
	{ 15.2493340, 0.887257 },	// 21
	{ 14.0719150, 1.095258 },	// 22
	{ -2.2569930, 0.896878 }	// 23
};
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
			{ 1.000000	,0.000000},
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

*/
// PATRICK'S CORRECTION PARAMETERS
Float_t xcal_cuts[24][2] = {
	{ 0.01, 0.96 }, // 00
	{ 0.00, 0.96 }, // 01
	{ 0.04, 1.00 }, // 02
	{ 0.03, 0.97 }, // 03
	{ 0.03, 1.00 }, // 04
	{ 0.05, 1.00 }, // 05
	{ 0.13, 0.87 }, // 06
	{ 0.06, 0.89 }, // 07
	{ 0.00, 0.87 }, // 08
	{ 0.01, 0.99 }, // 09
	{ 0.05, 0.94 }, // 10
	{ 0.00, 1.00 }, // 11
	{ 0.00, 1.00 }, // 12
	{ 0.00, 0.85 }, // 13
	{ 0.08, 0.88 }, // 14
	{ 0.04, 0.89 }, // 15
	{ 0.00, 0.96 }, // 16
	{ 0.02, 1.00 }, // 17
	{ 0.05, 1.00 }, // 18
	{ 0.00, 0.96 }, // 19
	{ 0.04, 1.00 }, // 20
	{ 0.06, 1.00 }, // 21
	{ 0.06, 0.99 }, // 22
	{ 0.00, 0.95 }  // 23
};

Double_t xnCorr[24] = {
	0.966727,	// 00
	0.997927,	// 01
	0.989599,	// 02
	0.926229,	// 03
	1.047585,	// 04
	0.909146,	// 05
	0.929986,	// 06
	0.967946,	// 07
	1.024587,	// 08
	0.962987,	// 09
	1.004322,	// 10
	0.000000,	// 11
	0.990265,	// 12
	1.001459,	// 13
	1.030865,	// 14
	0.971375,	// 15
	1.030609,	// 16
	0.925257,	// 17
	0.947191,	// 18
	1.003948,	// 19
	0.922519,	// 20
	0.992650,	// 21
	1.048827,	// 22
	1.099749,	// 23
};

// 0 is intercept, 1 is gradient
Double_t xfxneCorr[24][2] = {
	{   3.375253, 0.887861 },	// 00
	{   9.919062, 0.918479 },	// 01
	{   0.271875, 1.039087 },	// 02
	{   3.018420, 1.063643 },	// 03
	{   0.694993, 0.882137 },	// 04
	{  12.472202, 0.972475 },	// 05
	{   4.401876, 1.008029 },	// 06
	{   2.290693, 0.964587 },	// 07
	{   0.896536, 0.986746 },	// 08
	{   0.066819, 0.970263 },	// 09
	{  -2.399563, 0.909507 },	// 10
	{   0.000000, 0.000000 },	// 11
	{   0.699133, 0.783315 },	// 12
	{  -0.418725, 0.920520 },	// 13
	{  -0.859874, 0.933119 },	// 14
	{   0.475548, 0.984137 },	// 15
	{  -2.388185, 0.991268 },	// 16
	{   8.690010, 0.966768 },	// 17
	{   0.462283, 0.910253 },	// 18
	{  20.404622, 0.924948 },	// 19
	{   0.399083, 1.036339 },	// 20
	{   1.132734, 0.888264 },	// 21
	{  27.645136, 1.051020 },	// 22
	{  -1.097706, 0.874151 },	// 23
};

Double_t eCorr[24][2] = {		
	{ 256.762578, 0.021931 },	// 00
	{ 254.152463, 0.017396 },	// 01
	{ 278.386749, 0.018180 },	// 02
	{ 269.216169, 0.022511 },	// 03
	{ 250.658532, 0.019702 },	// 04
	{ 248.773636, 0.050975 },	// 05
	{ 263.349381, 0.024779 },	// 06
	{ 253.917318, 0.033730 },	// 07
	{ 263.003792, 0.018787 },	// 08
	{ 256.833004, 0.022207 },	// 09
	{ 240.852271, 0.026121 },	// 10
	{ 000.000000, 0.000000 },	// 11
	{ 214.975175, 0.023908 },	// 12
	{ 250.237648, 0.009948 },	// 13
	{ 244.299282, 0.017374 },	// 14
	{ 253.116975, 0.026214 },	// 15
	{ 263.419502, 0.015851 },	// 16
	{ 258.624883, 0.037380 },	// 17
	{ 252.003156, 0.085604 },	// 18
	{ 261.440891, 0.063373 },	// 19
	{ 266.547243, 0.019006 },	// 20
	{ 251.115918, 0.060770 },	// 21
	{ 293.415905, 0.032550 },	// 22
	{ 240.751075, 0.041021 }	// 23
};		
	





Float_t exCorr[6] = { 938.272,  // mass of proton [MeV/c^2]
	                   1,        // charge of proton
	                   27954.0982, // cm frame total energy [correct]
	                   26996.5929, // mass of recoil [correct] (Fully stripped of electrons)
	                   0.132178, // beta to CM frame [correct]
	                   2.5}; // Bfield [T]
Float_t exCorr_si[6] = { 938.272,  // mass of proton [MeV/c^2]
	                   1,        // charge of proton
	                   27949.6742, // cm frame total energy
	                   26984.277, // mass of recoil
	                   0.132187, // beta to CM frame
	                   2.5}; // Bfield [T]
Double_t array_radius = ISSArrayRadius(-4.5,4.5,11.5); // perpendicular distance of detector to axis [mm]
//double Ex, thetaCM;

Double_t alpha, Et, beta, gamm, G, massB, mass; //variables for Ex calculation
Double_t alpha_si, Et_si, beta_si, gamm_si, G_si, massB_si, mass_si; //variables for Ex calculation

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
//TString cutFileDir = "/home/ptmac/Documents/07-CERN-ISS-Mg/analysis/working/ALL-MgCuts3.root";
TString cutFileDir = "/home/ptmac/Documents/07-CERN-ISS-Mg/analysis/working/mg_cuts.root";

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
	Float_t Ex_si[24];
	Float_t Ex_corrected[24];
	Float_t thetaCM[24];
	Int_t detID[24];
	int td_rdt_e[24][4];
	int td_rdt_elum[32][4];
	int td_e_ebis[24];
	TCutG* cut[100];
	Float_t xold[24];
} FIN;

FIN fin;

// TSELECTOR BEGIN FUNCTION -------------------------------------------------------------------- //
void PTMonitors::Begin(TTree *tree){
	// Define offset (array position - offset position = 70mm???)
	//OFF_POSITION = GetArrayPosition( tree );


	if ( OFF_POSITION == 0 ){
		z_off = 4.9765;
	}
	else if ( OFF_POSITION == 1 ){
		z_off = 9.498;
		//std::copy( &xcal_cuts1[0][0], &xcal_cuts1[0][0] + 24*4, &xcal_cuts[0][0] );
		//std::copy( &thetaCM_limsBOTH[0][0], &thetaCM_limsBOTH[0][9], &thetaCM_lims[9] );
	}
	else if ( OFF_POSITION == 2 ){
		z_off = 6.50;
		//std::copy( &xcal_cuts2[0][0], &xcal_cuts2[0][0] + 24*4, &xcal_cuts[0][0] );
		//std::copy( &thetaCM_limsBOTH[1][0], &thetaCM_limsBOTH[1][9], &thetaCM_lims[9] );
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
	G = alpha*gamm*beta*array_radius;
	massB = exCorr[3];
	mass = exCorr[0];
	Et = exCorr[2];
	
	alpha_si = 299.792458 * exCorr_si[5] * exCorr_si[1] / TMath::TwoPi() / 1000; //MeV/mm
	beta_si = exCorr_si[4];
	gamm_si = 1./TMath::Sqrt(1-beta_si*beta_si);
	G_si = alpha_si*gamm_si*beta_si*array_radius;
	massB_si = exCorr_si[3];
	mass_si = exCorr_si[0];
	Et_si = exCorr_si[2];

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
	std::cout << ConstructFinFileName( tree ) << "\n";
	outFile = new TFile( ConstructFinFileName( tree ), "RECREATE");

	fin_tree = new TTree( "fin_tree", "Tree containing everything" );
	fin_tree->Branch("e",e,"e[100]/F");
	fin_tree->Branch("e_t",e_t,"e_t[100]/l");
	fin_tree->Branch("xf",xf,"xf[100]/F");
	fin_tree->Branch("xf_t",xf_t,"xf_t[100]/l");
	fin_tree->Branch("xn",xn,"xn[100]/F");
  	fin_tree->Branch("xn_t",xn_t,"xn_t[100]/l");
	fin_tree->Branch("rdt",rdt,"rdt[100]/F");
	fin_tree->Branch("rdt_t",rdt_t,"rdt_t[100]/l");
	fin_tree->Branch("tac",tac,"tac[100]/F");
	fin_tree->Branch("tac_t",tac_t,"tac_t[100]/l");
	fin_tree->Branch("elum",elum,"elum[32]/F");
	fin_tree->Branch("elum_t",elum_t,"elum_t[32]/l");
	fin_tree->Branch("ezero",ezero,"ezero[10]/F");
	fin_tree->Branch("ezero_t",ezero_t,"ezero_t[10]/l");
	//fin_tree->Branch("ebis_t",ebis_t,"EBISTimestamp/l");

	fin_tree->Branch("x",fin.x,"x[24]/F");
	fin_tree->Branch("z",fin.z,"z[24]/F");
	fin_tree->Branch("xcal",fin.xcal,"xcal[24]/F");
	fin_tree->Branch("ecal",fin.ecal,"ecal[24]/F");
	fin_tree->Branch("xfcal",fin.xfcal,"xfcal[24]/F");
	fin_tree->Branch("xncal",fin.xncal,"xncal[24]/F");
	fin_tree->Branch("ecrr",fin.ecrr,"ecrr[24]/F");
	fin_tree->Branch("td_rdt_e",fin.td_rdt_e,"td_rdt_e[24][4]/I");
	fin_tree->Branch("td_rdt_elum",fin.td_rdt_elum,"td_rdt_elum[32][4]/I");
	fin_tree->Branch("td_e_ebis",fin.td_e_ebis,"td_e_ebis[24]/I");
	fin_tree->Branch("Ex",fin.Ex,"Ex[24]/F");
	fin_tree->Branch("Ex_si",fin.Ex_si,"Ex_si[24]/F");
	fin_tree->Branch("Ex_corrected",fin.Ex_corrected,"Ex_corrected[24]/F");
	fin_tree->Branch("thetaCM",fin.thetaCM,"thetaCM[24]/F");
	fin_tree->Branch("detID",fin.detID,"detID[24]/I");
	fin_tree->Branch("td_rdt_e_cuts",td_rdt_e_cuts,"td_rdt_e_cuts[24][2]/I");
	//fin_tree->Branch("xcal_cuts",xcal_cuts,"xcal_cuts[24][4]/F");
	fin_tree->Branch("xcal_cuts",xcal_cuts,"xcal_cuts[24][2]/F");
	fin_tree->Branch("xold",fin.xold,"xold[24]/F");
	fin_tree->Branch("thetaCM_lims", thetaCM_lims, "thetaCM_lims[9]/F");
	fin_tree->Branch("ex_lims", ex_lims, "ex_lims[10]/F");

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
				fin.Ex_si[i] = TMath::QuietNaN();
				fin.Ex_corrected[i] = TMath::QuietNaN();
		 		fin.thetaCM[i] = TMath::QuietNaN();
		 		fin.detID[i] = TMath::QuietNaN();
				fin.td_e_ebis[i] = TMath::QuietNaN();
				fin.xold[i] = TMath::QuietNaN();
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

			fin.xold[i] = 0.5*( ( fin.xfcal[i] - fin.xncal[i] )/e[i] + 1 );


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
						//fin.Ex_corrected[index] = excitation_energy_corr_pars[OFF_POSITION - 1][j][0]*(450.0/9.0)*( fin.Ex[index] + 1.0 ) + excitation_energy_corr_pars[OFF_POSITION - 1][j][1];
						fin.Ex_corrected[index] = ex_corr[0][0]*fin.Ex[index] + ex_corr[1][0];

						double hahaha1 = gamm* TMath::Sqrt(mass * mass + momt * momt) - y;
						double hahaha2 = gamm* beta * momt;
						fin.thetaCM[index] = TMath::ACos(hahaha1/hahaha2) * TMath::RadToDeg();

					}
					else{
						fin.Ex[index] = TMath::QuietNaN();
						fin.Ex_corrected[index] = TMath::QuietNaN();
						fin.thetaCM[index] = TMath::QuietNaN();
					}
				}
				else{
					fin.Ex[index] = TMath::QuietNaN();
					fin.Ex_corrected[index] = TMath::QuietNaN();
					fin.thetaCM[index] = TMath::QuietNaN();
				}
				
				// <> SI CALCULATION
				double y_si = fin.ecrr[index] + mass_si; // to give the KE + mass of proton;
				double Z_si = alpha * gamm_si * beta_si * fin.z[index] * 10.;
				double H_si = TMath::Sqrt(TMath::Power(gamm_si * beta_si,2) * (y_si*y_si - mass_si * mass_si) ) ;

				// Calculate the angle
				if( TMath::Abs(Z_si) < H_si ) {
					// Use Newton's method to solve 0 ==  H * sin(phi) - G * tan(phi) - Z = f(phi)
					double tolerance_si = 0.001;	// Desired precision
			 	 	double phi_si = 0; 			// Initial phi = 0 -> ensure the solution has f'(phi) > 0
					double nPhi_si = 0; 			// New phi

					int iter_si = 0;				// Number of iterations

					// Now calculate the angle using Newton-Raphson process
					do{
						phi_si = nPhi_si;
						nPhi_si = phi_si - (H_si * TMath::Sin(phi_si) - G_si * TMath::Tan(phi_si) - Z_si) / (H_si * TMath::Cos(phi_si) - G_si /TMath::Power( TMath::Cos(phi_si), 2));
						iter_si++;
						if( iter_si > 10 || TMath::Abs(nPhi_si) > TMath::PiOver2()) break;
					} while( TMath::Abs(phi_si - nPhi_si ) > tolerance_si);
					phi_si = nPhi_si;

					// Check f'(phi) > 0
					double Df_si = H_si * TMath::Cos(phi_si) - G_si / TMath::Power( TMath::Cos(phi_si),2);
					if( Df_si > 0 && TMath::Abs(phi_si) < TMath::PiOver2()  ){
						// Found correct value of phi - now calculate everything else
						double K_si = H_si * TMath::Sin(phi_si);
						double x_si = TMath::ACos( mass_si / ( y_si * gamm_si - K_si));
						double momt_si = mass_si * TMath::Tan( x_si ); // momentum of particle b or B in CM frame
						double EB_si = TMath::Sqrt(mass_si*mass_si + Et_si*Et_si - 2*Et_si*TMath::Sqrt(momt_si * momt_si + mass_si * mass_si));
						fin.Ex_si[index] = EB_si - massB_si;
						//fin.Ex_corrected[index] = excitation_energy_corr_pars[OFF_POSITION - 1][j][0]*(450.0/9.0)*( fin.Ex[index] + 1.0 ) + excitation_energy_corr_pars[OFF_POSITION - 1][j][1];

						double hahaha1_si = gamm_si* TMath::Sqrt(mass_si * mass_si + momt_si * momt_si) - y_si;
						double hahaha2_si = gamm_si* beta_si * momt_si;
						//fin.thetaCM[index] = TMath::ACos(hahaha1/hahaha2) * TMath::RadToDeg();

					}
					else{
						fin.Ex_si[index] = TMath::QuietNaN();
						//fin.Ex_corrected[index] = TMath::QuietNaN();
						//fin.thetaCM[index] = TMath::QuietNaN();
					}
				}
				else{
					fin.Ex_si[index] = TMath::QuietNaN();
					//fin.Ex_corrected[index] = TMath::QuietNaN();
					//fin.thetaCM[index] = TMath::QuietNaN();
				}

				// </> SI CALIBRATION


				// Calculate the EBIS time - the array time and populate a histogram
				fin.td_e_ebis[index] = 10000;
				if ( ebis_t != 0 && e_t[index] != 0 ){
					fin.td_e_ebis[index] = (int)(e_t[index] - ebis_t);
				}
				TD_EBIS->Fill( fin.td_e_ebis[index] );


				// Calculate the recoil time stuff, by populating arrays with junk if not satisfying requirements
				for ( Int_t kk = 0; kk < 4; kk++ ){
					if ( rdt_t[kk] > 0 && e_t[index] > 0 ){
						fin.td_rdt_e[index][kk]= (int)(rdt_t[kk]-e_t[index]);
					}
					else{
						fin.td_rdt_e[index][kk] = 10000;
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
