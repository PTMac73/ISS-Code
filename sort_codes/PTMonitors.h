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
#include <TCutG.h>

// CUT ARRAYS
// Time difference between recoil detectors and array detectors
Int_t td_rdt_e_cuts[24][2] = 	{ {-9, 6},
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
	{-21, 6} };

// Cuts on xcal for each of the detectors in the array
Float_t xcal_cuts[24][2] = 	{ {0, 0.96},
	{0, 0.99},
	{0.04, 0.98},
	{0, 0.96},
	{0.05, 1},
	{0.07, 0.95},
	{0, 1},
	{0, 1},
	{0, 0.89},
	{0.04, 0.99},
	{0.06, 1},
	{0.01, 1},
	{0, 1},
	{0, 0.9},
	{0, 0.9},
	{0, 0.89},
	{0.05, 0.89},
	{0.11, 0.84},
	{0.03, 0.9},
	{0.1, 0.98},
	{0.03, 1},
	{0.03, 1},
	{0.02, 1},
	{0, 0.99} };

// Define the correction parameters based on [position][row #][gradient/intercept]
Double_t excitation_energy_corr_pars[2][6][2] = {
	{
		{ 0.02012, -1.00236 },
		{ 0.02037, -1.03572 },
		{ 0.01988, -0.97739 },
		{ 0.01962, -0.90813 },
		{ 0.01954, -0.89121 },
		{ 0.01949, -0.86839 }
	},
	{
		{ 0.02012, -1.00236 },
		{ 0.01987, -0.97342 },
		{ 0.02, -0.99215 },
		{ 0.01936, -0.87484 },
		{ 0.01945, -0.8742 },
		{ 0.01939, -0.86215 }
	}
};



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

#endif // #ifdef PTMonitors_cxx
