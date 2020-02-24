//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Feb 20 10:39:09 2020 by ROOT version 6.14/00
// from TTree fin_tree/Tree containing everything
// found on file: ../../root-data/fin1.root
//////////////////////////////////////////////////////////

#ifndef AnalyseTree_h
#define AnalyseTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

Bool_t det_array[6][4] = {
	{ 1, 1, 1, 1 },
	{ 1, 0, 0, 1 },
	{ 0, 1, 1, 1 },
	{ 1, 1, 0, 1 },
	{ 1, 1, 1, 1 },
	{ 1, 1, 0, 1 }
};

// Headers needed by this particular selector


class AnalyseTree : public TSelector {
public :
   TTree          *fChain = 0;   //!pointer to the analyzed TTree or TChain

	// Declare leaves
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
	// ----------------------------------------------------------------------------------------- //
	Float_t         x[24];
	Float_t         z[24];
	Float_t         xcal[24];
	Float_t         ecal[24];
	Float_t         xfcal[24];
	Float_t         xncal[24];
	Float_t         ecrr[24];
	Int_t           td_rdt_e[24][4];
	Float_t         Ex[24];
	Float_t         Ex_corrected[24];
	Float_t         thetaCM[24];
	Int_t           detID[24];
	Int_t           td_rdt_e_cuts[24][2];
	Float_t         xcal_cuts[24][4];
	Int_t           td_rdt_elum[32][4];
	Float_t         xold[24];
	
	// Branches to hold it all
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
	TBranch        *b_X;
	TBranch        *b_Z;
	TBranch        *b_XCAL;
	TBranch        *b_ECAL;
	TBranch        *b_XFCAL;
	TBranch        *b_XNCAL;
	TBranch        *b_ECRR;
	TBranch        *b_TD_RDT_E;
	TBranch        *b_Ex;
	TBranch        *b_Ex_CORRECTED;
	TBranch        *b_ThetaCM;
	TBranch        *b_DetID;
	TBranch        *b_TD_RDT_E_CUTS;
	TBranch        *b_XCAL_CUTS;
	TBranch        *b_TD_RDT_ELUM;
	TBranch        *b_XOLD;
	
	
   


   AnalyseTree(TTree *t = 0) { }
   virtual ~AnalyseTree() { }
   virtual Int_t   Version() const { return 2; }
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   virtual void    Init(TTree *tree);
   virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
   virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList  *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate();
   virtual void    Terminate();

   ClassDef(AnalyseTree,0);

};

#endif

#ifdef AnalyseTree_cxx
void AnalyseTree::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the reader is initialized.
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
	//fChain->SetBranchAddress("EBIS", &ebis_t, &b_EBISTimestamp);
	fChain->SetBranchAddress( "x", x, &b_X );
	fChain->SetBranchAddress( "z", z, &b_Z );
	fChain->SetBranchAddress( "xcal", xcal, &b_XCAL );
	fChain->SetBranchAddress( "ecal", ecal, &b_ECAL );
	fChain->SetBranchAddress( "xfcal", xfcal, &b_XFCAL );
	fChain->SetBranchAddress( "xncal", xncal, &b_XNCAL );
	fChain->SetBranchAddress( "ecrr", ecrr, &b_ECRR );
	fChain->SetBranchAddress( "td_rdt_e", td_rdt_e, &b_TD_RDT_E );
	fChain->SetBranchAddress( "Ex", Ex, &b_Ex );
	fChain->SetBranchAddress( "Ex_corrected", Ex_corrected, &b_Ex_CORRECTED );
	fChain->SetBranchAddress( "thetaCM", thetaCM, &b_ThetaCM );
	fChain->SetBranchAddress( "detID", detID, &b_DetID );
	fChain->SetBranchAddress( "td_rdt_e_cuts", td_rdt_e_cuts, &b_TD_RDT_E_CUTS );
	fChain->SetBranchAddress( "xcal_cuts", xcal_cuts, &b_XCAL_CUTS );
	fChain->SetBranchAddress( "td_rdt_elum", td_rdt_elum, &b_TD_RDT_ELUM );
	fChain->SetBranchAddress( "xold", xold, &b_XOLD );
}

Bool_t AnalyseTree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}


#endif // #ifdef AnalyseTree_cxx
