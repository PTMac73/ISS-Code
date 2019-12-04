// PTP_DBD_SetType.h
// Set the type of histogram to be generated
// ============================================================================================= //
// Patrick MacGregor
// Nuclear Physics Research Group
// School of Physics and Astronomy
// The University of Manchester
// ============================================================================================= //
#ifndef PTP_DBD_SET_TYPE_H_
#define PTP_DBD_SET_TYPE_H_

#include "PTP_DBD_classH24.h"

// Make stuff based on type
/* TYPE DEFINITIONS
	* 0 --> xcal
	* 1 --> ecrr:xcal
	* 2 --> xn:xf
	* 3 --> xncal:xfcal
*/
H24 SetType( Int_t type ){
	// XCAL
	if ( type == 0 ){
		H24 h24(1);
		h24.SetHistVars( "xcal" );
		h24.SetHistName( "h_xcal_%i" );
		h24.SetHistTitle( "xcal plot #%i" );
		h24.SetCutString( "( cut0 || cut1 || cut2 || cut3 ) && thetaCM > 11 && td_rdt_e[] > td_rdt_e_cuts[][0] && td_rdt_e[] < td_rdt_e_cuts[][1] && xcal[] > xcal_cuts[][0] && xcal[] < xcal_cuts[][1] && detID == %i", 0);
		h24.SetCutString( "( cut0 || cut1 || cut2 || cut3 ) && thetaCM > 11 && td_rdt_e[] > td_rdt_e_cuts[][0] && td_rdt_e[] < td_rdt_e_cuts[][1] && detID == %i", 1);
		h24.SetNbins( 200 );
		h24.SetLB( -0.5 );
		h24.SetUB( 1.5  );
		h24.SetBG( true );
		h24.SetPrintFolder( "xcal" );
		return h24;
	}
	
	// ECRR:XCAL
	else if ( type == 1 ){
		H24 h24(2);
		h24.SetHistVars( "ecrr:xcal" );
		h24.SetHistName( "h_ecrr_xcal_%i" );
		h24.SetHistTitle( "ecrr/xcal plot #%i" );
		h24.SetCutString( "( cut0 || cut1 || cut2 || cut3 ) && thetaCM > 11 && td_rdt_e[] > td_rdt_e_cuts[][0] && td_rdt_e[] < td_rdt_e_cuts[][1] && xcal[] > xcal_cuts[][0] && xcal[] < xcal_cuts[][1] && detID == %i", 0);
		h24.SetCutString( "", 1);
		h24.SetNbins( new Int_t[2] { 200, 900 } );
		h24.SetLB( new Double_t[2] { -0.5, -1.0 }  );
		h24.SetUB( new Double_t[2] { 1.5, 8.0 }   );
		h24.SetBG( false );
		h24.SetPrintFolder( "e-xcal" );
		return h24;
	}

	// XN:XF
	else if ( type == 2 ){
		H24 h24(2);
		h24.SetHistVars( "xn:xf" );
		h24.SetHistName( "h_xn_xf_%i" );
		h24.SetHistTitle( "xn/xf plot #%i" );
		h24.SetCutString( "( cut0 || cut1 || cut2 || cut3 ) && thetaCM > 11 && td_rdt_e[] > td_rdt_e_cuts[][0] && td_rdt_e[] < td_rdt_e_cuts[][1] && xcal[] > xcal_cuts[][0] && xcal[] < xcal_cuts[][1] && detID == %i", 0);
		h24.SetCutString( "", 1);
		h24.SetNbins( new Int_t[2] { 1000, 1000 } );
		h24.SetLB( new Double_t[2] { 0.0, 0.0 }  );
		h24.SetUB( new Double_t[2] { 2000.0, 2000.0 }   );
		h24.SetBG( false );
		h24.SetPrintFolder( "xn-xf" );
		return h24;
	}
	// XNcal:XFcal
	else if ( type == 3 ){
		H24 h24(2);
		h24.SetHistVars( "xncal:xfcal" );
		h24.SetHistName( "h_xncal_xfcal_%i" );
		h24.SetHistTitle( "xncal/xfcal plot #%i" );
		h24.SetCutString( "( cut0 || cut1 || cut2 || cut3 ) && thetaCM > 11 && td_rdt_e[] > td_rdt_e_cuts[][0] && td_rdt_e[] < td_rdt_e_cuts[][1] && xcal[] > xcal_cuts[][0] && xcal[] < xcal_cuts[][1] && detID == %i", 0);
		h24.SetCutString( "", 1);
		h24.SetNbins( new Int_t[2] { 1000, 1000 } );
		h24.SetLB( new Double_t[2] { 0.0, 0.0 }  );
		h24.SetUB( new Double_t[2] { 2000.0, 2000.0 }   );
		h24.SetBG( false );
		h24.SetPrintFolder( "xncal-xfcal" );
		return h24;
	}
	// XN:XF - ALPHA
	else if ( type == 4 ){
		H24 h24(2);
		h24.SetHistVars( "xn:xf" );
		h24.SetHistName( "h_xn_xf_a_%i" );
		h24.SetHistTitle( "xn/xf plot #%i" );
		h24.SetCutString( "xcal[] > xcal_cuts[][0] && xcal[] < xcal_cuts[][1] && detID == %i", 0);
		h24.SetCutString( "", 1);
		h24.SetNbins( new Int_t[2] { 1000, 1000 } );
		h24.SetLB( new Double_t[2] { 0.0, 0.0 }  );
		h24.SetUB( new Double_t[2] { 2000.0, 2000.0 }   );
		h24.SetBG( false );
		h24.SetPrintFolder( "xn-xf" );
		return h24;
	}
	// XNcal:XFcal - ALPHA
	else if ( type == 5 ){
		H24 h24(2);
		h24.SetHistVars( "xncal:xfcal" );
		h24.SetHistName( "h_xncal_xfcal_a_%i" );
		h24.SetHistTitle( "xncal/xfcal plot #%i" );
		h24.SetCutString( "xcal[] > xcal_cuts[][0] && xcal[] < xcal_cuts[][1] && detID == %i", 0);
		h24.SetCutString( "", 1);
		h24.SetNbins( new Int_t[2] { 1000, 1000 } );
		h24.SetLB( new Double_t[2] { 0.0, 0.0 }  );
		h24.SetUB( new Double_t[2] { 2000.0, 2000.0 }   );
		h24.SetBG( false );
		h24.SetPrintFolder( "xncal-xfcal" );
		return h24;
	}


	else{
		std::cout << "Unknown type. Try again." << endl;
		exit(1);
	}
}




// FORMAT THE HISTOGRAMS
void FormatHistograms( Int_t type, TH1* h_fg, TH1* h_bg = NULL ){
	// GLOBAL SETTINGS
	h_fg->SetTitle("");
	if ( h_bg != NULL ){
		h_bg->SetTitle("");
	}


	// 0 - XCAL
	if ( type == 0 ){
		h_fg->SetFillColor(5);
		h_fg->GetXaxis()->SetTitle("xcal");
		h_fg->GetYaxis()->SetTitle("Counts");
		if ( h_bg != NULL ){
			h_bg->SetFillColor(2);
		}
	}
	// 1 - 
	else if ( type == 1 ){
		h_fg->SetOption("scat");
		h_fg->SetMarkerStyle(5);
		h_fg->SetMarkerSize(1.0);
		h_fg->GetXaxis()->SetTitle("xcal");
		h_fg->GetYaxis()->SetTitle("E (MeV)");
		h_fg->GetXaxis()->SetTitleSize(0.06);
		h_fg->GetYaxis()->SetTitleSize(0.06);
		h_fg->GetXaxis()->SetLabelSize(0.06);
		h_fg->GetYaxis()->SetLabelSize(0.06);
	}
	else if ( type == 2 || type == 4 ){
		h_fg->SetOption("scat");
		h_fg->SetMarkerStyle(5);
		h_fg->SetMarkerSize(1.0);
		h_fg->GetXaxis()->SetTitle("xf");
		h_fg->GetYaxis()->SetTitle("xn");
		h_fg->GetXaxis()->SetTitleSize(0.06);
		h_fg->GetYaxis()->SetTitleSize(0.06);
		h_fg->GetXaxis()->SetLabelSize(0.06);
		h_fg->GetYaxis()->SetLabelSize(0.06);
	}
	else if ( type == 3 || type == 5 ){
		h_fg->SetOption("scat");
		h_fg->SetMarkerStyle(5);
		h_fg->SetMarkerSize(1.0);
		h_fg->GetXaxis()->SetTitle("xfcal");
		h_fg->GetYaxis()->SetTitle("xncal");
		h_fg->GetXaxis()->SetTitleSize(0.06);
		h_fg->GetYaxis()->SetTitleSize(0.06);
		h_fg->GetXaxis()->SetLabelSize(0.06);
		h_fg->GetYaxis()->SetLabelSize(0.06);
	}
}



#endif
