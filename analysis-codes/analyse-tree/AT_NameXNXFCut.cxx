// AT_NameXNXFCut.cxx
// Names cuts and stuff drawn by hand for the XN-XF cut stuff
// ============================================================================================= //
// Patrick MacGregor
// Nuclear Physics Research Group
// Department of Physics and Astronomy
// The University of Manchester
// ============================================================================================= //

#include <TCutG.h>
#include <TFile.h>


void WriteXNXFCut( TCutG* cut, Int_t i = 0, TFile* f = NULL ){
	cut->SetVarX( Form( "xf[%i]", i ) );
	cut->SetVarY( Form( "xn[%i]", i ) );
	cut->SetTitle("");
	cut->SetName( Form( "xnxf_cut_%i", i ) );
	
	if ( f != NULL ){
		f->cd();
		cut->Write();
	}
	
	return;
}
