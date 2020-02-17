// AG_functions.C
// Constants for the file ArrayGeometry.C
// ============================================================================================= //
// Patrick MacGregor
// Nuclear Physics Research Group
// Department of Physics and Astronomy
// The University of Manchester
// ============================================================================================= //
#ifndef AG_FUNCTIONS_H_
#define AG_FUNCTIONS_H_

#include "AG_constants.h"

#include <TMath.h>
#include <TString.h>

// UNIT CONVERSION FUNCTIONS
const Double_t AMU = 931.494;		// MeV

Double_t UC_MToMM(){ return 1000.0; }
Double_t UC_MToCM(){ return 100.0; }
Double_t UC_KGToMEV(){ return TMath::C()*TMath::C()/( TMath::Qe()*1E6 ); }
Double_t UC_EToC(){ return TMath::Qe(); }
Double_t UC_CToMS(){ return TMath::C(); }

// MASS FUNCTIONS
Double_t MassExcessToMass( Double_t D, Int_t A ){ return A*AMU + D; } // MeV

/* Calculate the radius at a particular z [cm]
	* v_para, v_perp 	[c]
	* cyc_freq 			[s^-1]
	* z					[cm]
*/
Double_t RadiusCalculator( Double_t v_para, Double_t v_perp, Double_t cyc_freq, Double_t z ){
	Double_t r = v_perp*UC_CToMS()*UC_MToCM()/cyc_freq;		// [cm]
	return TMath::Abs( 2*r*TMath::Sin( cyc_freq*z/( 2*v_para*UC_CToMS()*UC_MToCM() ) ) );	// [cm]
}

Double_t FWHMToSigma(){ return 1.0/( 2*TMath::Sqrt(2*TMath::Log(2.0) ) ); }

Double_t GetArrayRDTDistance( Int_t POS ){
	//       ( Array-Target [Edith]          ) + ( Target-Recoil [Edith]         ) - Si  - centroid - ( offset )
	return ( ( 40.53 + 1713.4 - 717.5 + 0.75 ) + ( 168.3 + 419.5 - 32.225 - 9.10 ) - 25 - 66.76 - ( POS == 1 ? 94.98 : 65.0 ) )*0.1;
}

TString DecimalDotToUnderscore( Double_t num ){
	TString a = Form( "%5.4f", num );
	return a.ReplaceAll( '.', '_' );
}




#endif
