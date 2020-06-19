// EnergyCalibration.C
// Calibrates the energies based on gf3 positions and ENSDF data
// ============================================================================================= //
// Patrick MacGregor
// Nuclear Physics Research Group
// Department of Physics and Astronomy
// The University of Manchester
// ============================================================================================= //

#include <iostream>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TF1.h>
#include <TMatrixDSym.h>
#include <TROOT.h>
#include "EnergyCalibration.h"

// Define global fit parameters
FitResult_t fit1[NUM_POSI][NUM_ROWS];	// For calibrating positions to energy
//FitResult_t fit2[NUM_POSI][NUM_ROWS];	// For calibrating energies to position


void EnergyCalibration(){

	// Decide whether to print graphs
	if ( swDISP_CANV == 1 ){ gROOT->SetBatch(kFALSE); }
	else{ gROOT->SetBatch(kTRUE); }

	// Print summary table
	PrintSummaryTable();


	// Declare some variables
	TGraphErrors *g_ecal[NUM_POSI][NUM_ROWS];	// Graphs to hold the calibrations of position to energy
	TCanvas *c_ecal;							// Canvas to hold the calibrations of position to energy
	Int_t num_points = 0;						// Variable holding the number of points for each calibration
	TF1* sline = new TF1( "sline", "[0] + [1]*x", 0, 300 );	// Generic straight line fit to be used for each calibration
	FormatLine(sline);
	
	
	// Set up canvas for drawing
	c_ecal = new TCanvas( "c_ecal", "Calibrate position to energy fits", C_WIDTH, C_HEIGHT );
	c_ecal->SetCanvasSize(2*C_WIDTH,2*C_HEIGHT);
	GlobSetCanvasMargins( c_ecal );
	c_ecal->Divide(4,3);
	

	// *LOOP* OVER POSITIONS (i)
	for ( Int_t i = 0; i < NUM_POSI; i++ ){
		
		// *LOOP* OVER ROWS (j)
		for ( Int_t j = 0; j < NUM_ROWS; j++ ){
			
			// Calculate the number of usable points
			num_points = NumUsableStates(i,j);
			
			// Print line divider
			if ( swVERBOSE ){ std::cout << std::setfill('-') << std::setw(80) << "" << std::setfill(' ') << "\n"; }
			
			// If not enough points to form
			if ( num_points < 2 ){
				// Don't carry out rest if not plottable
				std::cout << "Pos " << i << ", Row " << j << " is not plottable. Continuing." << "\n";
				SetEmptyFit( fit1[i][j] );
				//SetEmptyFit( fit2[i][j] );
				continue;
			}
			
			if ( swVERBOSE ){ std::cout << "Plotting Pos " << i << ", Row " << j << "\n"; }
			
			// Define arrays to hold the variables from the data table
			Double_t* x = new Double_t[num_points];
			Double_t* x_err = new Double_t[num_points];
			Double_t* y = new Double_t[num_points];
			Double_t* y_err = new Double_t[num_points];
			
			// Populate arrays
			Int_t l = 0;	// Fixes which index of the array it enters
			for ( Int_t k = 0; k < NUM_STAT; k++ ){
				if ( use_states[k] && data_table[i][j][k][POS_INDEX] >= 0 ){
					x[k - l] = data_table[i][j][k][0];
					x_err[k - l] = data_table[i][j][k][1];
					y[k - l] = data_table[i][j][k][2];
					y_err[k - l] = data_table[i][j][k][3];
				}
				else{
					l++;
				}
			}
			
			// Print arrays
			if ( swVERBOSE ){
				std::cout << "PRINT ARRAYS: Pos " << i << ", Row " << j << ", # = " << num_points <<  "\n";
				PrintArray( x, x_err, y, y_err, num_points );
			}

			
			// Create TGraphs and TCanvases and format them
			
			g_ecal[i][j] = new TGraphErrors( num_points, x, y, x_err, y_err );
			FormatGraph( g_ecal[i][j] );
			g_ecal[i][j]->SetTitle( Form( "Position %i | Row %i", i + 1, j ) );
			
			c_ecal->cd( GetCanvasNumber(i,j) );
			g_ecal[i][j]->Draw("AP");
			
			
			// Get straight line parameters
			TFitResultPtr s = g_ecal[i][j]->Fit(sline, "SQ");
			if ( swVERBOSE ){ std::cout << "\u03c7^2 = " << s->Chi2() << "\n"; }
			TMatrixDSym* covmc = (TMatrixDSym*)s->GetCovarianceMatrix().Clone();
			SetFitParameters( sline->GetParameter(1), sline->GetParameter(0), covmc, fit1[i][j] );
			if ( swVERBOSE ){ PrintFitParameters( fit1[i][j] ); }
			
			// Calculate the reversed parameters
			/*CalculateNewFitParameters( fit1[i][j], fit2[i][j] );
			
			if ( swVERBOSE ){
				std::cout << "NEW PARAMETERS\n";
				PrintFitParameters( fit2[i][j] );
			}*/
			

			// Delete the arrays
			delete[] x;
			delete[] x_err;
			delete[] y;
			delete[] y_err;
			

		}// *LOOP* OVER ROWS (j)
		
	}// *LOOP* OVER POSITIONS (i)
	
	
	// CALCULATE THE AVERAGE FIT (FOR EMPTY STATES)
	FitResult_t average_fit;
	SetEmptyFit( average_fit );
	Int_t num_full_fits = 0;
	
	// Get the average and the max errors
	for ( Int_t i = 0; i < NUM_POSI; i++ ){
		for ( Int_t j = 0; j < NUM_ROWS; j++ ){
			if ( fit1[i][j].empty_flag == 0 ){
				average_fit.m += fit1[i][j].m;
				average_fit.c += fit1[i][j].c;
				num_full_fits++;
				
				if ( (*fit1[i][j].covmc)(0,0) >= (*average_fit.covmc)(0,0) ){
					(*average_fit.covmc)(0,0) = (*fit1[i][j].covmc)(0,0);
				}
				if ( (*fit1[i][j].covmc)(1,1) >= (*average_fit.covmc)(1,1) ){
					(*average_fit.covmc)(1,1) = (*fit1[i][j].covmc)(1,1);
				}
				if ( (*fit1[i][j].covmc)(1,0) >= (*average_fit.covmc)(1,0) || (*average_fit.covmc)(1,0) == 0  ){
					(*average_fit.covmc)(1,0) = (*fit1[i][j].covmc)(1,0);
					(*average_fit.covmc)(0,1) = (*fit1[i][j].covmc)(0,1);
				}
			}
		}
	}
	
	average_fit.m = average_fit.m/num_full_fits;
	average_fit.c = average_fit.c/num_full_fits;
	
	
	
	// Print the fit parameters to a table based on spreadsheet
	// Table 1
	std::cout << "\n\n>>> TABLE 1 <<<" << "\n\n";
	Int_t w = 12;
	TString sp = "    ";
	std::cout << std::left << std::setw(w) << "> M <" << sp << std::setw(w) << "> C <" << sp << std::setw(w) << "> M <" << sp << std::setw(w) << "> C <" << "\n";
	
	for ( Int_t j = 0; j < NUM_ROWS; j++ ){
		std::cout 
			<< std::setw(w) << ( fit1[0][j].empty_flag == 0 ? fit1[0][j].m : average_fit.m ) << sp 
			<< std::setw(w) << ( fit1[0][j].empty_flag == 0 ? fit1[0][j].c : average_fit.c ) << sp 
			<< std::setw(w) << ( fit1[1][j].empty_flag == 0 ? fit1[1][j].m : average_fit.m ) << sp 
			<< std::setw(w) << ( fit1[1][j].empty_flag == 0 ? fit1[1][j].c : average_fit.c ) << "\n";
	}
	
	// Table 2
	std::cout << "\n\n>>> TABLE 2 <<<" << "\n\n";
	std::cout << std::left << std::setw(w) << "> \u03c3(M) <" << sp << std::setw(w) << "> \u03c3(C) <" << sp << std::setw(w) << "> cov(M,C) <" << sp << std::setw(w) << "> \u03c3(M) <" << sp << std::setw(w) << "> \u03c3(C) <" << sp << std::setw(w) << "> cov(M,C) <" << "\n";
	
	for ( Int_t j = 0; j < NUM_ROWS; j++ ){
		std::cout 
			<< std::setw(w) << TMath::Sqrt( ( fit1[0][j].empty_flag == 0 ? (*fit1[0][j].covmc)(1,1)  : (*average_fit.covmc)(1,1) ) ) << sp
			<< std::setw(w) << TMath::Sqrt( ( fit1[0][j].empty_flag == 0 ? (*fit1[0][j].covmc)(0,0)  : (*average_fit.covmc)(0,0) ) ) << sp
			<< std::setw(w) << ( fit1[0][j].empty_flag == 0 ? (*fit1[0][j].covmc)(1,0)  : (*average_fit.covmc)(1,0) ) << sp
			<< std::setw(w) << TMath::Sqrt( ( fit1[1][j].empty_flag == 0 ? (*fit1[1][j].covmc)(1,1)  : (*average_fit.covmc)(1,1) ) ) << sp
			<< std::setw(w) << TMath::Sqrt( ( fit1[1][j].empty_flag == 0 ? (*fit1[1][j].covmc)(0,0)  : (*average_fit.covmc)(0,0) ) ) << sp
			<< std::setw(w) << ( fit1[1][j].empty_flag == 0 ? (*fit1[1][j].covmc)(1,0)  : (*average_fit.covmc)(1,0) ) << sp << "\n";
	}
	
	// Print the canvas
	c_ecal->Print("energy_cal_fits.pdf");
	

	std::cout << "Deleted arrays from memory\n";

	return;
}
