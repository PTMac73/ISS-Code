// PTP_DBD_classH24.h
// Class for drawing histograms
// ============================================================================================= //
// Patrick MacGregor
// Nuclear Physics Research Group
// School of Physics and Astronomy
// The University of Manchester
// ============================================================================================= //
#ifndef PTP_DBD_CLASS_HIST_24_H_
#define PTP_DBD_CLASS_HIST_24_H_

#include <vector>
#include <iostream>
#include <TString.h>

const Int_t NUM_PLOTS = 24;

class H24 {
	private:
		Int_t fDim;				// 1D or 2D histogram
		TString	fHistVars;		// E.g. xcal or ecrr:xcal				
		TString fHistName;		// E.g. h_xcal_%i
		TString fHistTitle;
		TString fCutString[2];	// 0 = fg, 1 = bg
		Bool_t fBG;
		Int_t* fNbins;
		Double_t* fLB;
		Double_t* fUB;
		TH1* fTHD_fg[NUM_PLOTS];
		TH1* fTHD_bg[NUM_PLOTS];
		TString fPrintFolder;
	
	public:
		// Constructor #1
		H24(){
			fDim = 1; 
			fHistVars =  "";
			fHistName = "";
			fHistTitle = "";
			fCutString[0] = "";
			fCutString[1] = "";
			fPrintFolder = "";
			fBG = 0;
			fNbins = new Int_t[1] { 100 };
			fLB = new Double_t[1] { 0.0 };
			fUB = new Double_t[1] { 1.0 };
			for ( Int_t i = 0; i < NUM_PLOTS; i++ ){
				fTHD_bg[i] = NULL;
				fTHD_fg[i] = NULL;
			}
		}

		// Constructor #2
		H24( Int_t hdim ){
			fDim =  hdim; 
			fHistVars = "";
			fHistName = "";
			fHistTitle = "";
			fCutString[0] = "";
			fCutString[1] = "";
			fPrintFolder = "";
			fBG = 0;
			fNbins = new Int_t[hdim];
			fLB = new Double_t[hdim];
			fUB = new Double_t[hdim];
			for ( Int_t i = 0; i < hdim; i++ ){
				fNbins[i] = 100;
				fLB[i] = 0.0;
				fUB[i] = 1.0;
			}
			for ( Int_t i = 0; i < NUM_PLOTS; i++ ){
				fTHD_bg[i] = NULL;
				fTHD_fg[i] = NULL;
			}
		}
	
		// Destructor
		~H24(){}

		// Setters
		void SetHistVars( TString h_vars ){ fHistVars = h_vars; }
		void SetHistName( TString h_name ){ fHistName = h_name; }
		void SetHistTitle( TString h_title ){ fHistTitle = h_title; }
		void SetBG( Bool_t h_bg ){ fBG = h_bg; }
		void SetPrintFolder( TString h_pf ){ fPrintFolder = h_pf; }

		void SetNbins( Int_t *hbins ){ fNbins = hbins; }
		void SetLB( Double_t* hlb ){ fLB = hlb; }
		void SetUB( Double_t* hub ){ fUB = hub; }

		void SetNbins( Int_t hbins, Int_t index ){ fNbins[index] = hbins; }
		void SetLB( Double_t hlb, Int_t index ){ fLB[index] = hlb; }
		void SetUB( Double_t hub, Int_t index ){ fUB[index] = hub; }
		void SetCutString( TString h_cut_string, Int_t index ){ fCutString[index] = h_cut_string; }

		void SetNbins( Int_t hbins ){ fNbins[0] = hbins; }
		void SetLB( Double_t hlb ){ fLB[0] = hlb; }
		void SetUB( Double_t hub ){ fUB[0] = hub; }
		void SetCutString( TString h_cut_string ){ fCutString[0] = h_cut_string; }

		void SetHBG( TH1* h_bg, Int_t index ){ fTHD_bg[index] = h_bg; }
		void SetHBG( TH1* h_bg ){ fTHD_bg[0] = h_bg; }
		void SetHFG( TH1* h_fg, Int_t index ){ fTHD_fg[index] = h_fg; }
		void SetHFG( TH1* h_fg ){ fTHD_fg[0] = h_fg; }

		// Getters
		Int_t GetDim(){ return fDim; }
		TString	GetHistVars(){ return fHistVars; }
		TString GetHistName(){ return fHistName; }
		TString GetHistName( Int_t index ){
			if ( index == 1 ){ return fHistName + "_bg"; }
			else{ return fHistName; }
		}
		TString GetPrintFolder(){ return fPrintFolder; }
		Bool_t GetBG(){ return fBG; }

		Int_t* GetNbins(){ return fNbins; }
		Double_t* GetLB(){ return fLB; }
		Double_t* GetUB(){ return fUB; }
		TString* GetCutString(){ return fCutString; }

		Int_t GetNbins( Int_t index ){ return ( ( index <= fDim ) ? fNbins[index] : -1 ); }
		Double_t GetLB( Int_t index ){ return ( ( index <= fDim ) ? fLB[index] : -1.0 ); }
		Double_t GetUB( Int_t index ){ return ( ( index <= fDim ) ? fUB[index] : -1.0 ); }
		TString GetCutString( Int_t index ){ return ( ( index <= fDim ) ? fCutString[index] : "" ); }

		TH1* GetHBG( Int_t index ){ return fTHD_bg[index]; }
		TH1* GetHFG( Int_t index ){ return fTHD_fg[index]; }

		// Print function
		void Print(){
			std::cout << "Dim      \t" << fDim << std::endl;
			std::cout << "Hist Name\t" << fHistName << std::endl;
			std::cout << "Hist Vars\t" << fHistVars << std::endl;
		}

		// Make hist string
		TString MakeHistStr( Int_t FGorBG){
			TString hist_str;
			hist_str +=Form( "%s>>", fHistVars.Data() );
			if ( FGorBG == 0 ){
				hist_str +=Form("%s(", fHistName.Data() );
			}
			else{
				hist_str +=Form("%s_bg(", fHistName.Data() );
			}
			for ( Int_t i = 0; i < fDim; i++ ){
				hist_str += Form( "%i, %f, %f", fNbins[i], fLB[i], fUB[i] );
				if ( i < fDim - 1 ){ hist_str += ", "; }
			}
			hist_str += ")";
			return hist_str;
		}

		TString MakeHistStr(){
			return MakeHistStr(0);
		}





};

#endif
