// EnergyCalibration.C
// Calibrates the energies based on ENSDF data and the fit data produced from Ryan's code in 
// PTMonitors.C
// ============================================================================================= //
// Patrick MacGregor
// Nuclear Physics Research Group
// Department of Physics and Astronomy
// The University of Manchester
// ============================================================================================= //

#include <TCanvas.h>
#include <TFitResultPtr.h>
#include <TGraphErrors.h>

#include "EnergyCalibration.h"



void EnergyCalibration(){
	// Prepare canvas
	TCanvas* c_graph = new TCanvas( "c_graph", "CANVAS", C_WIDTH, C_HEIGHT );
	TPad* pad = (TPad*)c_graph->GetPad(0);
	pad->SetLeftMargin(0.08);
	pad->SetRightMargin(0.01);
	pad->SetTopMargin(0.01);
	pad->SetBottomMargin(0.08);

	// Declare arrays
	Double_t* x = new Double_t[NUM_POINTS];
	Double_t* x_err = new Double_t[NUM_POINTS];
	Double_t* y = new Double_t[NUM_POINTS];
	Double_t* y_err = new Double_t[NUM_POINTS];
	
	for ( Int_t i = 0; i < NUM_POINTS; i++ ){
		x[i] = monitors_data[i][0];
		x_err[i] = monitors_data[i][1];
		y[i] = ensdf_data[i][0];
		y_err[i] = ensdf_data[i][1];
	}
	
	TGraphErrors* g = new TGraphErrors( NUM_POINTS, x, y, x_err, y_err );
	g->SetTitle("");
	g->GetXaxis()->SetTitle("PTMonitors.C");
	g->GetXaxis()->SetLabelFont(62);
	g->GetXaxis()->SetTitleFont(62);
	g->GetYaxis()->SetTitle("ENDSF");
	g->GetYaxis()->SetLabelFont(62);
	g->GetYaxis()->SetTitleFont(62);
	g->Draw("AP");
		
	TF1* fit_func = new TF1( "fit_func", "[0]*x + [1]", 0, 5 );
	fit_func->SetLineWidth(1);
	TFitResultPtr r = g->Fit( fit_func, "SQ" );
	

	std::cout << "m\t" << std::setprecision(6) << std::setw(12) << r->Parameter(0) << "\t" << r->ParError(0) << "\n";
	std::cout << "c\t" << std::setprecision(6) <<  std::setw(12) << r->Parameter(1) << "\t" << r->ParError(1) << "\n";
	
	c_graph->Print("energy_calibration.pdf");
	

	return;
}
