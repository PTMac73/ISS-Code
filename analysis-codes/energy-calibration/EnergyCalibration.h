// EnergyCalibration.h
// Header file containing relevant variables and functions
// ============================================================================================= //
// Patrick MacGregor
// Nuclear Physics Research Group
// Department of Physics and Astronomy
// The University of Manchester
// ============================================================================================= //

#ifndef ENERGY_CALIBRATION_H_
#define ENERGY_CALIBRATION_H_

#include <TAxis.h>
#include <TCanvas.h>
#include <TFitResultPtr.h>
#include <TGraphErrors.h>
#include <TMatrixDSym.h>
#include <iostream>

// STRUCTS ---===---===---===---===---===---===---===---===---===---===---===---===---===---===- //
struct FitResult_t{
	Double_t m;
	Double_t c;
	TMatrixDSym* covmc;
	Bool_t empty_flag;
};



// CONSTANTS ---===---===---===---===---===---===---===---===---===---===---===---===---===---== //
// Numbers for arrays (data table attributes)
const Int_t NUM_POSI = 2;
const Int_t NUM_ROWS = 6;
const Int_t NUM_STAT = 6;
const Int_t ENERGY_INDEX = 2;
const Int_t POS_INDEX = 0;

// Canvas width
const Int_t C_WIDTH = 1600;
const Int_t C_HEIGHT = 900;

// Choose which states to use in the energy calibration
const Bool_t use_states[NUM_STAT] = {
		// E [MeV] (err ([keV])
	1,	// 1.09462 ( 0.21)
	1,	// 1.43070 ( 0.40)
	1,	// 2.26590 ( 0.70)
	1,	// 2.49990 ( 1.00)
	1,	// 3.22750 ( 0.60)
	1	// 4.28000 (40.00)
};

// SWITCHES ---===---===---===---===---===---===---===---===---===---===---===---===---===---=== //
Bool_t swDISP_CANV = 0;		// Display canvas (1) or run in batch mode (0)
Bool_t swVERBOSE = 0;		// Verbose output

// DATA TABLE (created in Mg-Analysis spreadsheet)
// Contains the x and y variables for each of the fits
// Order is [pos][row][state][col] 
// The four columns are position (x), position_err, energy (y), energy_err
Double_t data_table[NUM_POSI][NUM_ROWS][NUM_STAT][4] = {					
	{ 	{ 	{ -01.00, -1.0, 1.09462, 0.00021 }, 		
			{ -01.00, -1.0, 1.43070, 0.00040 }, 		
			{ -01.00, -1.0, 2.26590, 0.00070 }, 		
			{ -01.00, -1.0, 2.49990, 0.00100 }, 		
			{ -01.00, -1.0, 3.22750, 0.00060 }, 		
			{ -01.00, -1.0, 4.28000, 0.04000 }	},	
		{ 	{ 101.10, 0.30, 1.09462, 0.00021 }, 		
			{ 115.80, 0.60, 1.43070, 0.00040 }, 		
			{ -01.00, -1.0, 2.26590, 0.00070 }, 		
			{ -01.00, -1.0, 2.49990, 0.00100 }, 		
			{ -01.00, -1.0, 3.22750, 0.00060 }, 		
			{ -01.00, -1.0, 4.28000, 0.04000 }	},	
		{ 	{ 104.22, 0.20, 1.09462, 0.00021 }, 		
			{ 120.51, 0.11, 1.43070, 0.00040 }, 		
			{ 160.40, 0.70, 2.26590, 0.00070 }, 		
			{ 171.70, 0.30, 2.49990, 0.00100 }, 		
			{ -01.00, -1.0, 3.22750, 0.00060 }, 		
			{ -01.00, -1.0, 4.28000, 0.04000 }	},	
		{ 	{ 102.72, 0.23, 1.09462, 0.00021 }, 		
			{ 119.80, 0.13, 1.43070, 0.00040 }, 		
			{ 162.70, 0.60, 2.26590, 0.00070 }, 		
			{ 172.95, 0.22, 2.49990, 0.00100 }, 		
			{ 208.50, 0.30, 3.22750, 0.00060 }, 		
			{ -01.00, -1.0, 4.28000, 0.04000 }	},	
		{ 	{ 103.00, 0.30, 1.09462, 0.00021 }, 		
			{ 119.95, 0.18, 1.43070, 0.00040 }, 		
			{ 162.30, 0.60, 2.26590, 0.00070 }, 		
			{ 173.30, 0.40, 2.49990, 0.00100 }, 		
			{ 208.90, 0.50, 3.22750, 0.00060 }, 		
			{ 266.03, 0.19, 4.28000, 0.04000 }	},	
		{ 	{ 102.60, 0.30, 1.09462, 0.00021 }, 		
			{ 119.10, 0.30, 1.43070, 0.00040 }, 		
			{ 160.50, 0.80, 2.26590, 0.00070 }, 		
			{ 173.50, 0.60, 2.49990, 0.00100 }, 		
			{ 207.50, 0.70, 3.22750, 0.00060 }, 		
			{ 264.20, 0.21, 4.28000, 0.04000 }	}	}, 
	{ 	{ 	{ -01.00, -1.0, 1.09462, 0.00021 }, 		
			{ -01.00, -1.0, 1.43070, 0.00040 }, 		
			{ -01.00, -1.0, 2.26590, 0.00070 }, 		
			{ -01.00, -1.0, 2.49990, 0.00100 }, 		
			{ -01.00, -1.0, 3.22750, 0.00060 }, 		
			{ -01.00, -1.0, 4.28000, 0.04000 }	},	
		{ 	{ 103.20, 0.30, 1.09462, 0.00021 }, 		
			{ 118.80, 0.30, 1.43070, 0.00040 }, 		
			{ -01.00, -1.0, 2.26590, 0.00070 }, 		
			{ -01.00, -1.0, 2.49990, 0.00100 }, 		
			{ -01.00, -1.0, 3.22750, 0.00060 }, 		
			{ -01.00, -1.0, 4.28000, 0.04000 }	},	
		{ 	{ 103.60, 0.30, 1.09462, 0.00021 }, 		
			{ 120.84, 0.16, 1.43070, 0.00040 }, 		
			{ 163.50, 0.80, 2.26590, 0.00070 }, 		
			{ 173.80, 0.30, 2.49990, 0.00100 }, 		
			{ -01.00, -1.0, 3.22750, 0.00060 }, 		
			{ -01.00, -1.0, 4.28000, 0.04000 }	},	
		{ 	{ 103.20, 0.40, 1.09462, 0.00021 }, 		
			{ 120.00, 0.30, 1.43070, 0.00040 }, 		
			{ 160.00, 1.10, 2.26590, 0.00070 }, 		
			{ 172.90, 0.40, 2.49990, 0.00100 }, 		
			{ 208.50, 0.30, 3.22750, 0.00060 }, 		
			{ 264.60, 0.50, 4.28000, 0.04000 }	},	
		{ 	{ 102.80, 0.40, 1.09462, 0.00021 }, 		
			{ 119.70, 0.30, 1.43070, 0.00040 }, 		
			{ 160.70, 0.90, 2.26590, 0.00070 }, 		
			{ 172.50, 0.80, 2.49990, 0.00100 }, 		
			{ 208.90, 0.50, 3.22750, 0.00060 }, 		
			{ 266.30, 0.30, 4.28000, 0.04000 }	},	
		{ 	{ 101.70, 0.60, 1.09462, 0.00021 }, 		
			{ 119.20, 0.50, 1.43070, 0.00040 }, 		
			{ 159.70, 0.60, 2.26590, 0.00070 }, 		
			{ 171.80, 1.10, 2.49990, 0.00100 }, 		
			{ 207.50, 0.70, 3.22750, 0.00060 }, 		
			{ 264.10, 0.40, 4.28000, 0.04000 }	}	}, 
};					
		
/*Double_t data_table_pos1_full[NUM_POSI][NUM_STAT][4] = {		
	{	
		{ 102.91, 0.11, 1.09462, 0.00021 },
		{ 119.92, 0.08, 1.43070, 0.00040 },
		{ 161.60, 0.40, 2.26590, 0.00070 },
		{ 172.79, 0.17, 2.49990, 0.00100 },
		{ 208.40, 0.30, 3.22750, 0.00060 },
		{ 265.04, 0.14, 4.28000, 0.04000 }
	},{	
		{ 103.09, 0.17, 1.09462, 0.00021 },
		{ 120.02, 0.11, 1.43070, 0.00040 },
		{ 160.90, 0.50, 2.26590, 0.00070 },
		{ 173.16, 0.24, 2.49990, 0.00100 },
		{ 209.00, 0.50, 3.22750, 0.00060 },
		{ 265.18, 0.20, 4.28000, 0.04000 }
	}	
};		
*/
Double_t data_table_pos1_full[NUM_POSI][NUM_STAT][4] = {		
	{	
		{ 1.06784, 0.00252, 1.09462, 0.00021 }, 
		{ 1.40832, 0.00144, 1.43070, 0.00040 }, 
		{ 2.24114, 0.00606, 2.26590, 0.00070 }, 
		{ 2.46661, 0.00311, 2.49990, 0.00100 }, 
		{ 3.18367, 0.00491, 3.22750, 0.00060 }, 
		{ 4.30902, 0.00287, 4.28000, 0.04000 }, 
	},{	
		{ 1.07011, 0.00313, 1.09462, 0.00021 }, 
		{ 1.40890, 0.00219, 1.43070, 0.00040 }, 
		{ 2.23453, 0.00873, 2.26590, 0.00070 }, 
		{ 2.47067, 0.00385, 2.49990, 0.00100 }, 
		{ 3.17955, 0.00957, 3.22750, 0.00060 }, 
		{ 4.30902, 0.00000, 4.28000, 0.04000 }
	}	
};	

				


// FUNCTIONS ---===---===---===---===---===---===---===---===---===---===---===---===---===---== //
// Separates a subsection of the data table for a given graph
/*void CreateProbe( Double_t* pr[NUM_STAT], Int_t pos, Int_t row ){
	// Assigns the addresses of arrays in the data table to an array of pointers
	for ( Int_t i = 0; i < NUM_STAT; i++ ){
		pr[i] = data_table[pos][row][i];
	}
	return;
}


// Prints the probe
void PrintProbe( Double_t* pr[NUM_STAT] ){
	for ( Int_t i = 0; i < NUM_STAT; i++ ){
		for ( Int_t j = 0; j < 4; j++ ){
			std::cout << pr[i][j] << "\t";
		}
		std::cout << "\n";
	}
}*/


// Prints the arrays used for plotting the graphs
void PrintArray( Double_t* arr1, Double_t* arr2, Double_t* arr3, Double_t* arr4, Int_t num ){
	// Define the width of the field
	Int_t w = 8;
	
	// Print a header
	std::cout << std::setw(w) << "X" << std::setw(w) << "X_ERR" << std::setw(w) << "Y" << std::setw(w) << "Y_ERR" << "\n";
	
	// Print the values
	for ( Int_t i = 0; i < num; i++ ){
		std::cout << std::setw(w) << arr1[i] << std::setw(w) << arr2[i] << std::setw(w) << arr3[i] << std::setw(w) << arr4[i] << "\n";
	}
	std::cout << "\n";
}


// Calculate the number of usable states
Int_t NumUsableStates( Int_t pos, Int_t row ){
	Int_t sum = 0;
	
	// Sum up states which are being used in the calibration and contain actual data
	for ( Int_t i = 0; i < NUM_STAT; i++ ){
		sum += ( use_states[i] && ( data_table[pos][row][i][POS_INDEX] >= 0 ) );
	}
	return sum;
}



// GET FIT RESULTS FROM TFitResultPtr
void SetFitParameters( Double_t m, Double_t c, TMatrixDSym* covmc, FitResult_t& f ){
	f.m = m;
	f.c = c;
	f.covmc = (TMatrixDSym*)covmc->Clone();
	f.empty_flag = 0;
	return;
}


TString NegativeSpace( Double_t a ){
	if ( a < 0 ){ return ""; }
	else{ return " "; }
}


void Print2x2Matrix( TMatrixDSym* a ){
	Int_t w = 8;
	std::cout << "[ " << NegativeSpace( (*a)(0,0) ) << std::setw(w) << (*a)(0,0) << "\t" << NegativeSpace( (*a)(1,0) ) << std::setw(w) << (*a)(1,0) << " ]\n";
	std::cout << "[ " << NegativeSpace( (*a)(1,0) ) << std::setw(w) << (*a)(0,1) << "\t" << NegativeSpace( (*a)(1,1) ) << std::setw(w) << (*a)(1,1) << " ]\n";
	return;
}


void PrintFitParameters( FitResult_t &t ){
	Int_t w = 8;
	std::cout << "M = " << std::setw(w) << t.m << " \u00b1 " << std::setw(w) << TMath::Sqrt( (*t.covmc)(1,1) ) << "\n";
	std::cout << "C = " << std::setw(w) << t.c << " \u00b1 " << std::setw(w) << TMath::Sqrt( (*t.covmc)(0,0) ) << "\n";
	std::cout << "Covariance Matrix:\n";
	Print2x2Matrix( t.covmc );
	std::cout << "\n";
	return;
}


void CalculateNewFitParameters( FitResult_t& src, FitResult_t& dest ){
	TMatrixDSym covMC(2);
	TMatrixDSym covmc = (*src.covmc);
	
	// Calculate gradient and intercept
	dest.m = 1/src.m;
	dest.c = -src.c/src.m;
	
	// Variance on intercept
	covMC(0,0) = TMath::Power( dest.c, 2 )*covmc(1,1) + TMath::Power( dest.m, 2 )*covmc(0,0) + 2*TMath::Power( dest.m, 2 )*dest.c*covmc(1,0);
	
	// Variance on gradient
	covMC(1,1) = TMath::Power( dest.m, 4 )*covmc(1,1);
	
	// Covariance term
	covMC(1,0) = TMath::Power( dest.m, 3 )*( covmc(1,0) + dest.c );
	covMC(0,1) = covMC(1,0);
	
	dest.covmc = (TMatrixDSym*)covMC.Clone();
	
	return;
}

void SetEmptyFit( FitResult_t &f ){
	f.m = 0;
	f.c = 0;
	TMatrixDSym a(2);
	f.covmc = (TMatrixDSym*)a.Clone();
	f.empty_flag = 1;
	return;
}




// FORMAT STUFF ---===---===---===---===---===---===---===---===---===---===---===---===---===-- //
// Format LOBF
void FormatLine( TF1* f ){
	f->SetLineWidth(1);
	f->SetLineColor(kRed);
	return;
}


// Set graph fonts
void GlobSetFonts( TGraphErrors* g ){
	g->GetXaxis()->SetTitleFont(62);
	g->GetXaxis()->SetLabelFont(62);
	g->GetXaxis()->CenterTitle();
	g->GetYaxis()->SetTitleFont(62);
	g->GetYaxis()->SetLabelFont(62);
	g->GetYaxis()->CenterTitle();
	return;
}


// Format the graph using a convenient function
void FormatGraph( TGraphErrors* g ){
	g->GetYaxis()->SetTitle("Energy (MeV)");
	g->GetXaxis()->SetTitle("gf3 Position");
	g->SetMarkerStyle(20);
	g->SetMarkerSize(0.1);
	GlobSetFonts(g);
	return;
}


// Set the margins of the canvas
void GlobSetCanvasMargins( TCanvas *c, Double_t l = 0.1, Double_t r = 0.02, Double_t t = 0.02, Double_t b = 0.1 ){
	TPad* pad = (TPad*)c;
	pad->SetLeftMargin( l );
	pad->SetRightMargin( r );
	pad->SetTopMargin( t );
	pad->SetBottomMargin( b );
	return;
}

void PrintSummaryTable(){
	std::cout << "USING STATES:" << "\n";
	std::cout << "    E (MeV)\n";
	for ( Int_t i = 0; i < NUM_STAT; i++ ){
		std::cout << " " << ( use_states[i] == 1 ? "X" : " " ) << "  " << data_table[0][0][i][ENERGY_INDEX] << "\n";
	}
	std::cout << "\n\n";
	return;
}


Int_t GetCanvasNumber( Int_t pos, Int_t row ){
	return 2*row + pos + 1;
}
	






























#endif
