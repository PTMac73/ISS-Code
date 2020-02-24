// BeamContamination.C
// Calculates the beam contamination for a given position using the zero degrees detector
// ============================================================================================= //
// Patrick MacGregor
// Nuclear Physics Research Group
// Department of Physics and Astronomy
// The University of Manchester
// ============================================================================================= //
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <iostream>


void BeamContamination( TFile *f ){
	// Get the TTree
	TTree *t = (TTree*)f->Get("fin_tree");
	
	// Plot the spectrum
	gStyle->SetOptStat(0);
	TCanvas *c_bc = new TCanvas( "c_bc", "BEAM CONTAMINATION", 1200, 900 );
	t->Draw("ezero[0]:ezero[1]>>h_bc(500, 0, 2000, 500, 0, 2000)","","goff");
	TH2F* h_bc = (TH2F*)gDirectory->Get("h_bc");
	h_bc->SetTitle("E#DeltaE detector at 0#circ; #DeltaE; E");
	h_bc->Draw("colz");
	
	// Define points for cuts
	TCutG *cutMG = new TCutG("cutMG");
	cutMG->SetPoint( 0, 993.74, 865.561 );
	cutMG->SetPoint( 1, 27.5459, 1274.6 );
	cutMG->SetPoint( 2, 25.4591, 1077.23 );
	cutMG->SetPoint( 3, 1020.87, 659.611 );
	cutMG->SetPoint( 4, 1008.35, 856.979 );
	cutMG->SetPoint( 5, 993.74, 865.561 );

	TCutG *cutSI = new TCutG("cutSI");
	cutSI->SetPoint( 0, 680.718, 1062.93 );
	cutSI->SetPoint( 1, 92.237, 1340.39 );
	cutSI->SetPoint( 2, 83.8898, 1469.11 );
	cutSI->SetPoint( 3, 680.718, 1174.49 );
	cutSI->SetPoint( 4, 680.718, 1080.09 );
	cutSI->SetPoint( 5, 680.718, 1062.93 );
	
	// Format the cuts
	cutMG->SetLineColor(kRed);
	cutMG->SetLineWidth(2);
	cutSI->SetLineColor(kBlue);
	cutSI->SetLineWidth(2);
	
	cutMG->Draw("same");
	cutSI->Draw("same");

	// Print info about the integrals
	std::cout << std::setw(14) << "28Mg TOTAL: " << cutMG->IntegralHist( h_bc ) << std::endl;
	std::cout << std::setw(14) << "28Si TOTAL: " << cutSI->IntegralHist( h_bc ) << std::endl;
	std::cout << std::setw(14) << "28Mg + 28Si: " << cutSI->IntegralHist( h_bc ) + cutMG->IntegralHist( h_bc ) << std::endl;
	std::cout << std::setw(14) << "TOTAL: " << h_bc->Integral() << std::endl;
	
	// Print the graphs
	c_bc->Print("/home/ptmac/Documents/07-CERN-ISS-Mg/analysis/PLOTS/BeamContamination.pdf");
	
	
}
/*
	// Count the numbers
	Int_t counter[10];
	Float_t ezero[10];
	t->SetBranchAddress("ezero", &ezero );
	
	for ( Int_t i = 0; i < t->GetEntries(); i++ ){
		t->GetEntry(i);
		for ( Int_t j = 0; j < 10; j++ ){
			if ( i == 0 ){
				counter[j] = 0;
				printf("%f\t", ezero[j]);
			}
			if ( ! std::isnan(ezero[j]) ){
				counter[j] += 1;
			}
		}
	}
	printf("\n");
	std::cout << counter[0] << "\t" << counter[1] << std::endl;
	std::cout << counter[2] << "\t" << counter[3] << std::endl;
	std::cout << counter[4] << "\t" << counter[5] << std::endl;
	std::cout << counter[6] << "\t" << counter[7] << std::endl;
	std::cout << counter[8] << "\t" << counter[9] << std::endl;
}
*/
