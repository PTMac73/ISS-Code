// Define PTMonitors_h
#ifndef PTMonitors_h
#define PTMonitors_h

// Include some stuff
#include "WriteSPE.h"
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TH1.h>
#include <TH2.h>
#include <TStopwatch.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TCutG.h>
#include <iostream>

// CUT ARRAYS
// Time difference between recoil detectors and array detectors
/*Int_t td_rdt_e_cuts[24][2] = {
	{-9, 6},
	{-11, 4},
	{-9, 3},
	{-14, 3},
	{-16, 4},
	{-20, 3},
	{-5, 6},
	{-10, 3},
	{-8, 4},
	{-13, 4},
	{-19, 4},
	{-200, 200},
	{-12, 5},
	{-7, 7},
	{-9, 7},
	{-8, 5},
	{-20, 5},
	{-20, 6},
	{-11, 6},
	{-15, 6},
	{-11, 8},
	{-18, 5},
	{-22, 6},
	{-21, 6} 
};*/

Int_t td_rdt_e_cuts[24][2] = {
	{  -7,  4 },	// 00
	{ -11,  4 },	// 01
	{  -9,  4 },	// 02
	{ -14,  3 },	// 03
	{ -16,  4 },	// 04
	{ -20,  3 },	// 05
	{  -5,  6 },	// 06
	{ -10,  3 },	// 07
	{  -10,  5 },	// 08
	{ -13,  4 },	// 09
	{ -19,  4 },	// 10
	{-200,200 },	// 11
	{ -12,  5 },	// 12
	{  -7,  7 },	// 13
	{  -9,  7 },	// 14
	{  -8,  5 },	// 15
	{ -20,  5 },	// 16
	{ -20,  6 },	// 17
	{ -11,  6 },	// 18
	{ -15,  6 },	// 19
	{ -11,  8 },	// 20
	{ -18,  5 },	// 21
	{ -22,  6 },	// 22
	{ -16,  6 } 	// 23
};

// Cuts on xcal for each of the detectors in the array for both positions
// Format is [ LB, UB, mid-LB, mid-UB ] for each detector
// POS 1
/*Float_t xcal_cuts1[24][4] = {
	{ 0.00, 0.96, 0.27, 0.50 }, 
	{ 0.00, 0.99, 0.49, 0.64 }, 
	{ 0.04, 0.98, 0.52, 0.66 }, 
	{ 0.00, 0.96, 0.48, 0.59 }, 
	{ 0.05, 1.00, 1.00, 0.00 }, 
	{ 0.07, 0.95, 0.48, 0.61 }, 
	{ 0.00, 1.00, 1.00, 0.00 }, 
	{ 0.00, 1.00, 1.00, 0.00 }, 
	{ 0.00, 0.89, 1.00, 0.00 }, 
	{ 0.04, 0.99, 0.48, 0.60 }, 
	{ 0.06, 1.00, 1.00, 0.00 }, 
	{ 0.01, 1.00, 1.00, 0.00 }, 
	{ 0.00, 1.00, 1.00, 0.00 }, 
	{ 0.00, 0.90, 1.00, 0.00 }, 
	{ 0.00, 0.90, 1.00, 0.00 }, 
	{ 0.00, 0.89, 0.43, 0.55 }, 
	{ 0.05, 0.89, 0.39, 0.58 }, 
	{ 0.11, 0.84, 1.00, 0.00 }, 
	{ 0.03, 0.90, 0.50, 0.64 }, 
	{ 0.10, 0.98, 0.58, 0.63 }, 
	{ 0.03, 1.00, 1.00, 0.00 }, 
	{ 0.03, 1.00, 1.00, 0.00 }, 
	{ 0.02, 1.00, 1.00, 0.00 }, 
	{ 0.00, 0.99, 1.00, 0.00 }
};

// POS 2
Float_t xcal_cuts2[24][4] = {
	{ 0.00, 0.96, 0.16, 0.50 }, 
	{ 0.00, 0.99, 0.46, 0.70 }, 
	{ 0.04, 0.98, 1.00, 0.00 }, 
	{ 0.00, 0.96, 0.47, 0.56 }, 
	{ 0.05, 1.00, 1.00, 0.00 }, 
	{ 0.07, 0.95, 1.00, 0.00 }, 
	{ 0.00, 1.00, 1.00, 0.00 }, 
	{ 0.00, 1.00, 1.00, 0.00 }, 
	{ 0.00, 0.89, 1.00, 0.00 }, 
	{ 0.04, 0.99, 0.55, 0.65 }, 
	{ 0.06, 1.00, 1.00, 0.00 }, 
	{ 0.01, 1.00, 1.00, 0.00 }, 
	{ 0.00, 1.00, 1.00, 0.00 }, 
	{ 0.00, 0.90, 1.00, 0.00 }, 
	{ 0.00, 0.90, 1.00, 0.00 }, 
	{ 0.00, 0.89, 1.00, 0.00 }, 
	{ 0.05, 0.89, 1.00, 0.00 }, 
	{ 0.11, 0.84, 1.00, 0.00 }, 
	{ 0.03, 0.90, 0.49, 0.64 }, 
	{ 0.10, 0.98, 1.00, 0.00 }, 
	{ 0.03, 1.00, 0.45, 0.62 }, 
	{ 0.03, 1.00, 0.48, 0.60 }, 
	{ 0.02, 1.00, 1.00, 0.00 }, 
	{ 0.00, 0.99, 1.00, 0.00 }
}; */

// Define the correction parameters based on [position][row #][gradient/intercept]
Double_t excitation_energy_corr_pars[2][6][2] = {
	{
		{ 0.01963, -0.92447 },
		{ 0.01888, -0.85955 },
		{ 0.02057, -1.06259 },
		{ 0.01988, -0.95548 },
		{ 0.01937, -0.87164 },
		{ 0.01944, -0.87309 }
	},
	{
		{ 0.01979, -0.93888 },
		{ 0.02025, -1.01499 },
		{ 0.01996, -0.99294 },
		{ 0.02007, -0.97508 },
		{ 0.01936, -0.86454 },
		{ 0.01931, -0.84685 }
	}
};

Double_t ex_corr[2][2] = {
	{ 1.0056800, 0.00411759 },
	{ 0.0127572, 0.00629672 }
};

// ThetaCM limits (based on eyeballing, with array size 2 for position)
Double_t thetaCM_limsBOTH[2][9] = {
	{ 15.0, 12.5, 15.5, 15.0, 13.5, 0.0, 13.0, 13.0, 13.0 },
	{ 19.0, 12.5, 16.0, 15.0, 13.0, 0.0, 13.5, 13.0, 10.5 }
};

// This needs to be continuous - excitation spectrum limits
Double_t ex_lims[10] = { -5, 0.5, 1.3, 1.9, 2.4, 2.7, 3.0, 3.5, 4.1, 8.0 };

// --------------------------------------------------------------------------------------------- //
TString ConstructFinFileName( TTree *t ){
	// Get the name of the file
	TString file_name;
	TString out_name = "";
	Bool_t alpha_run = 0;
	
	if ( t->GetCurrentFile() != NULL ){
		// Get the file name from the TTree
		file_name = t->GetCurrentFile()->GetName();
		
		if ( !file_name.Contains("/") ){
			// Manipulate the string - looks like gen_run##.root or genPos##.root
			if ( file_name != "" ){
				if ( file_name.Contains( "gen_run" ) ){
					file_name.Remove( 0, 7 );
				}
				else if ( file_name.Contains( "genPos" ) ){
					file_name.Remove( 0, 6 );
				}
				else if ( file_name.Contains( "genAlpha" ) ){
					file_name.Remove( 0, 8 );
					alpha_run = 1;
				}
				else{
					file_name = "XXX.root";
				}
				out_name = out_name + "fin" + ( alpha_run == 1 ? "Alpha" : "" ) + file_name;
			}
			else{
				out_name = "finERROR.root";
			}
		}
		else{
			out_name = "fin0.root";
		}
	}
	else{
		out_name = "fin0.root";
	}
	return out_name;
}
// --------------------------------------------------------------------------------------------- //
Int_t GetArrayPosition( TTree *t ){
	// Get the name of the file
	TString file_name;
	if ( t->GetCurrentFile() != NULL ){
		file_name = t->GetCurrentFile()->GetName();
	}
	Int_t pos = 1;	// Assume it's 1

	// Manipulate the string - looks like gen_run##.root or genPos##.root
	if ( file_name.Contains( "genPos" ) ){
		file_name.Remove( 0, 6 );
		file_name.Remove( file_name.Length() - 5, 5 );
		pos = file_name.Atoi();
	}
	else{
		std::cout << "Double check position! Assuming it's " << pos << " ...\n";
	}
	
	return pos;
}
// --------------------------------------------------------------------------------------------- //

// DEFINE PTMONITORS TSelector CLASS HERE ------------------------------------------------------ //
class PTMonitors : public TSelector {
public :
	TTree          *fChain;			// Pointer to the analyzed TTree or TChain

	// Declaration of leaf types
	Float_t         e[100];			// Energy detected in Si array
	ULong64_t       e_t[100];		// ^^ timestamp
	Float_t         xf[100];		// Far voltage in Si array
	ULong64_t       xf_t[100];		// ^^ timestamp
	Float_t         xn[100];		// Near voltage in Si array
	ULong64_t       xn_t[100];		// ^^ timestamp
	Float_t         rdt[100];		// Recoil detectors
	ULong64_t       rdt_t[100];		// ^^ timestamp
	Float_t         tac[100];		// ???
	ULong64_t       tac_t[100];		// ^^ timestamp
	Float_t         elum[32];		// Luminosity detectors
	ULong64_t       elum_t[32];		// ^^ timestamp
	Float_t         ezero[10];		// ???
	ULong64_t       ezero_t[10];	// ^^ timestamp
	ULong64_t		ebis_t;		// EBIS timestamp

	// List of branches to hold said leaves
	TBranch        *b_Energy;   //!
	TBranch        *b_EnergyTimestamp;   //!
	TBranch        *b_XF;   //!
	TBranch        *b_XFTimestamp;   //!
	TBranch        *b_XN;   //!
	TBranch        *b_XNTimestamp;   //!
	TBranch        *b_RDT;   //!
	TBranch        *b_RDTTimestamp;   //!
	TBranch        *b_TAC;   //!
	TBranch        *b_TACTimestamp;   //!
	TBranch        *b_ELUM;   //!
	TBranch        *b_ELUMTimestamp;   //!
	TBranch        *b_EZERO;   //!
	TBranch        *b_EZEROTimestamp;   //!
	TBranch		   *b_EBISTimestamp;

	// CLASS MEMBER FUNCTIONS
	PTMonitors(TTree * /*tree*/ =0) : fChain(0) { }		// Constructor
	virtual ~PTMonitors() { }							// Destructor
	virtual Int_t   Version() const { return 3; }		// Version of this class
	
	// Declare required TSelector Functions
	virtual void    Begin(TTree *tree);
	virtual void    SlaveBegin(TTree *tree);
	virtual void    Init(TTree *tree);
	virtual Bool_t  Notify();
	virtual Bool_t  Process(Long64_t entry);
	virtual void    SlaveTerminate();
	virtual void    Terminate();
	
	// Other member functions
	virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
	virtual void    SetOption(const char *option) { fOption = option; }
	virtual void    SetObject(TObject *obj) { fObject = obj; }
	virtual void    SetInputList(TList *input) { fInput = input; }
	virtual TList  *GetOutputList() const { return fOutput; }

	ClassDef(PTMonitors,0);
};

#endif

#ifdef PTMonitors_cxx
void PTMonitors::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
	if (!tree) return;
	fChain = tree;
	fChain->SetMakeClass(1);

	fChain->SetBranchAddress("e", e, &b_Energy);
	fChain->SetBranchAddress("e_t", e_t, &b_EnergyTimestamp);
	fChain->SetBranchAddress("xf", xf, &b_XF);
	fChain->SetBranchAddress("xf_t", xf_t, &b_XFTimestamp);
	fChain->SetBranchAddress("xn", xn, &b_XN);
	fChain->SetBranchAddress("xn_t", xn_t, &b_XNTimestamp);
	fChain->SetBranchAddress("rdt", rdt, &b_RDT);
	fChain->SetBranchAddress("rdt_t", rdt_t, &b_RDTTimestamp);
	fChain->SetBranchAddress("tac", tac, &b_TAC);
	fChain->SetBranchAddress("tac_t", tac_t, &b_TACTimestamp);
	fChain->SetBranchAddress("elum", elum, &b_ELUM);
	fChain->SetBranchAddress("elum_t", elum_t, &b_ELUMTimestamp);
	fChain->SetBranchAddress("ezero", ezero, &b_EZERO);
	fChain->SetBranchAddress("ezero_t", ezero_t, &b_EZEROTimestamp);
	fChain->SetBranchAddress("EBIS", &ebis_t, &b_EBISTimestamp);
}

Bool_t PTMonitors::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

// Calculates the average radius of the ISS array (be consistent with units!!)
/*/ The X-axis has its origin in the centre of the array (not necessarily the centre of the strip)
 *  X1 is the point furthest left of the origin
 *  X2 is the point furthest right of the origin
 *  height is the height of the array                                                            */
Double_t ISSArrayRadius( Double_t X1, Double_t X2, Double_t height ){
	Double_t A1 = X1/height;
	Double_t A2 = X2/height;
	return ( height*height/( 2*( X2 - X1 ) ) )*( A2*TMath::Sqrt( 1 + A2*A2 ) - A1*TMath::Sqrt( 1 + A1*A1 ) + TMath::Log( TMath::Abs( A2 + TMath::Sqrt( 1 + A2*A2 ) ) ) - TMath::Log( TMath::Abs( A1 + TMath::Sqrt( 1 + A1*A1 ) ) ) );
}

#endif // #ifdef PTMonitors_cxx
