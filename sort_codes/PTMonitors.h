// Define PTMonitors_h
#ifndef PTMonitors_h
#define PTMonitors_h

// Include some stuff
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

// writespe function here - taken from gatemat2.cpp
void writespe(const Char_t *hisname, Char_t *spename, Char_t const*xy="X"){
	TH1		 *hist;
	Int_t	 i,j,NN,size;
	Int_t	 i1,i2;
	Int_t	 NBINS;
	Char_t	str[32];
	float *sp;
	NN = 4096;
	if ( strncmp( xy, "Y", 1 ) == 0 || strncmp( xy, "y", 1 ) == 0 ){ NN = 4096; }
	FILE *out;
	if ( !( sp = (float*) malloc( NN*sizeof( float ) ) ) ) {
		printf("\007	ERROR: Could not malloc data buffer.\n");
		exit(-1);
	}
	hist = (TH1*)gROOT->FindObject(hisname);
	if ( hist != NULL ){
		NBINS = hist->GetNbinsX();
		for( i = 1; i < NN + 1; i++ ){
			if ( i <= NBINS ){ 
				sp[i-1] = hist->GetBinContent(i);
			}
			else{
				sp[i-1] = 0.0;
			}
		}
		sprintf( str, "%s.spe", spename );
		out = fopen( str, "wb+" );
		i = 1;
		j = 24;
		fwrite( &j, 4, 1, out );
		fwrite( str, 8, 1, out );
		fwrite( &NN, 4, 1, out );
		fwrite( &i, 4, 1, out );
		fwrite( &i, 4, 1, out );
		fwrite( &i, 4, 1, out );
		fwrite( &j, 4, 1, out );
		size = sizeof(float)*NN;
		fwrite( &size, 4, 1, out );
		fwrite( sp, 4, NN, out );
		fwrite( &size, 4, 1, out );
		fclose( out );
		printf("wrote %i channels to %s\n", NN, str);
	} else {
		printf("spectrum %s not found\n", hisname);
	}
	free(sp);
	return;
}

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
