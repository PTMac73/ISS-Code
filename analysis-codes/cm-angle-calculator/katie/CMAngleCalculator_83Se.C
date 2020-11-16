// Calculates the CM angle for the reaction 28Mg(d,p)29Mg in inverse kinematics given an array position and an excitation energy
#include <iostream>
#include <fstream>
#include <TMath.h>
#include <TString.h>
#include <TGraph.h>
#include <TGraph2D.h>

//TString iFileDir = "/home/ptmac/Documents/07-CERN-ISS-Mg/Mg-Analysis/excitationEnergyList.txt";

TString iFileDir = "/home/katie/Documents/h073_83Kr/working/runs164_220ExFitting83Se/Ex_fits/new_finalsTHETA>9/z_edges2.txt";

// MAIN FUNCTION
// ex in MeV, z in cm
// In order for this to work, quote:
//	* Masses, momenta, and energies in terms of MeV and factors of c
//	* Velocities in terms of c
//	* Charges in terms of e
//	* Distances in terms of cm 
Double_t CMAngleCalculator(double ex, double z){
	// Define initial fixed quantities for the reaction
	int n = 0;								// Number of iterations
	double c = 299792458;					// Speed of light in m /s
	double q = 1;							// Charge of outgoing proton in units of e
	double u = 931.494;						// Unified atomic mass unit in MeV / c^2
	double B = 2.85*c*(1e-8);				// Convert B-field to MeV / (cm*c*e)
	double m1 = 81.916*u;					// m1 = 82Se in MeV / c^2		    (beam)  https://pubchem.ncbi.nlm.nih.gov/#query=82Se
	double m2 = 2.01410177785*u;			// m2 = d in MeV / c^2			(target)
	double m3 = 1.00782503207*u;			// m3 = p in MeV / c^2			(ejectile)
	double m4 = 82.919*u;					// m4 = 83Se in MeV / c^2		(recoil)	https://pubchem.ncbi.nlm.nih.gov/#query=83Se

	// Now calculate derived quantities from fixed quantities
	double T1 = 10.012*m1/u;									// Kinetic energy of the beam [LAB] [ 9.473 MeV / u ]
	double e1 = T1 + m1;										// Total energy of beam particle [LAB] in MeV
	double etot = e1 + m2;										// Total energy [LAB] in MeV
	double etot_cm = TMath::Sqrt( m1*m1 + m2*m2 + 2*e1*m2 );	// Total energy [CM] (related to invariant mass) in MeV
	double gamLab_CM = etot/etot_cm;									// Gamma factor relating inertial frames (i.e. between LAB and CM)
	double beta = TMath::Sqrt( 1 - TMath::Power( gamLab_CM, -2 ) );		// Beta = v/c (v = velocity of CM frame)
	double rho = 1.15;													// Radius of array in cm (approximate as circular - diameter of square is 23 mm).

	double m4ex = m4 + ex;															// Invariant mass of the recoil nucleus
	double e3_cm = 0.5*( etot_cm*etot_cm + m3*m3 - m4ex*m4ex )/etot_cm;				// Energy of the outgoing proton [CM]
	double p3_cm = TMath::Sqrt(e3_cm*e3_cm - m3*m3);								// Momentum of the outgoing proton [CM]

	double p_para = q*B*z/( 2*TMath::Pi() );										// Parallel component of the momentum in MeV / c [LAB] -------> THIS IS A GUESS!!
	double p_para_cm = p_para/gamLab_CM - beta*e3_cm;								// Parallel component of the momentum in MeV / c [CM]
	double p_perp_cm = TMath::Sqrt( e3_cm*e3_cm - p_para_cm*p_para_cm - m3*m3 );	// Perpendicular component of the momentum in MeV / c [CM]

	// Calculate the function to be minimised and its first derivative
	double fp = ( 2*p_perp_cm/( q*B ) )*TMath::Sin( q*B*z/( 2*gamLab_CM*( p_para_cm + beta*e3_cm ) ) ) - rho;
	double fd = ( ( -2*p_para_cm )/( q*B*p_perp_cm ) )*TMath::Sin( q*B*z/( 2*gamLab_CM*( p_para_cm + beta*e3_cm ) ) ) - ( ( p_perp_cm*z )/( gamLab_CM*TMath::Power( p_para_cm + beta*e3_cm , 2 ) ) )*TMath::Cos( q*B*z/( 2*gamLab_CM*( p_para_cm + beta*e3_cm ) ) );

	// Now perform Newton-Raphson calculation to get the angle in the CM frame
	while ( TMath::Abs(fp) > 1e-5 && n < 100000){
		p_para_cm = p_para_cm - fp/fd;
		p_perp_cm = TMath::Sqrt( e3_cm*e3_cm - p_para_cm*p_para_cm - m3*m3 );

		fp = ( 2*p_perp_cm/( q*B ) )*TMath::Sin( q*B*z/( 2*gamLab_CM*( p_para_cm + beta*e3_cm ) ) ) - rho;
		fd = ( ( -2*p_para_cm )/( q*B*p_perp_cm ) )*TMath::Sin( q*B*z/( 2*gamLab_CM*( p_para_cm + beta*e3_cm ) ) ) - ( ( p_perp_cm*z )/( gamLab_CM*TMath::Power( p_para_cm + beta*e3_cm , 2 ) ) )*TMath::Cos( q*B*z/( 2*gamLab_CM*( p_para_cm + beta*e3_cm ) ) );
		n++;
	}

	// Function minimised - now calculate the final quantities
	double theta_cm = 180 - ( TMath::ACos( p_para_cm/p3_cm )*180/TMath::Pi() );
	double e3 = gamLab_CM*( e3_cm + beta*p_para_cm );
	double p3 = TMath::Sqrt( e3*e3 - m3*m3 );
	double p_perp = p_perp_cm;
	p_para = TMath::Sqrt( p3*p3 - p_perp*p_perp );
	double theta = TMath::ACos( gamLab_CM*( p_para_cm + beta*e3_cm )/p3 )*180/TMath::Pi();
	std::cout << ex << "\t" << z << "\t" << theta_cm << std::endl;
	//std::cout << n << "\t" << ex << "\t" << z << "\t" << theta_cm << "\t" << theta << "\t" << e3-m3 << "\t" << std::endl;
	//std::cout << p3 << "\t" << p_para << "\t" << p_perp << "\t" << p3_cm << "\t" << p_para_cm << "\t" << p_perp_cm << "\t" << p3_cm*TMath::Sin(theta_cm) << "\t" << p3*TMath::Sin(theta) << std::endl;

	return theta_cm;
}

/*************************************************************************************************/
void CMAngleBatch( TString iFileDir ){
	// Open the file for reading	
	ifstream iFile;
	iFile.open( iFileDir.Data() );

	// Check that it opened
	if ( iFile.is_open() ){
		// Print success message
		std::cout << "File opened successfully!" << std::endl;

		// Read the data
		double ex, z;
		std::cout << "EX\tz\tTH_CM" << std::endl;
		//std::cout << "n\tEx\tz\tTH_CM\tTH\tE3-m3" << std::endl;
		//std::cout << "P3\tP3_PARA\tP3_PERP\tP3CM\tP3_PACM\tP3_PECM\tP_PerpCalcCM\tP_PerpCalc" << std::endl;
		while ( !iFile.eof() ){
			// Store the values
			iFile >> ex >> z;
			
			// Calculate the angles
			CMAngleCalculator( ex, z );
		}

		// Close the file
		iFile.close();
		std::cout << "File closed (I think...)!" << std::endl;
	}
	else{
		// Print failure message
		std::cout << "File did not open. Process terminated." << std::endl;
	}
}
/*************************************************************************************************/
/*
void CMAngleCalculatorRyan(){

	//======== Ex calculation by Ryan 
	double y = ecrr[detID] + mass; // to give the KE + mass of proton;
	double Z = alpha * gamm * beta * z[detID] * 10.;
	double H = TMath::Sqrt(TMath::Power(gamm * beta,2) * (y*y - mass * mass) ) ;

	// Calculate the angle
	if( TMath::Abs(Z) < H ) {
		// Use Newton's method to solve 0 ==  H * sin(phi) - G * tan(phi) - Z = f(phi) 
		double tolerance = 0.001;	// Desired precision
 	 	double phi = 0; 			// Initial phi = 0 -> ensure the solution has f'(phi) > 0
		double nPhi = 0; 			// New phi

		int iter = 0;				// Number of iterations
			
		// Now calculate the angle
		do{
			phi = nPhi;
			nPhi = phi - (H * TMath::Sin(phi) - G * TMath::Tan(phi) - Z) / (H * TMath::Cos(phi) - G /TMath::Power( TMath::Cos(phi), 2));
			iter ++;
			if( iter > 10 || TMath::Abs(nPhi) > TMath::PiOver2()) break;
		} while( TMath::Abs(phi - nPhi ) > tolerance);
		phi = nPhi;

		// Check f'(phi) > 0
		double Df = H * TMath::Cos(phi) - G / TMath::Power( TMath::Cos(phi),2);
		if( Df > 0 && TMath::Abs(phi) < TMath::PiOver2()  ){
			double K = H * TMath::Sin(phi);
			double x = TMath::ACos( mass / ( y * gamm - K));
			double momt = mass * TMath::Tan(x); // momentum of particle b or B in CM frame
			double EB = TMath::Sqrt(mass*mass + Et*Et - 2*Et*TMath::Sqrt(momt * momt + mass * mass));
			Ex = EB - massB;
			
			double hahaha1 = gamm* TMath::Sqrt(mass * mass + momt * momt) - y;
			double hahaha2 = gamm* beta * momt;
			thetaCM = TMath::ACos(hahaha1/hahaha2) * TMath::RadToDeg();
	 
		}
		else{
			Ex = TMath::QuietNaN();
			thetaCM = TMath::QuietNaN();
		}	
	}
	else{
		Ex = TMath::QuietNaN();
		thetaCM = TMath::QuietNaN();
	}

}
*/

void CMAngleGraph2D( ){
	// Define the z's
	Double_t z_spacing = 0.05;
	Double_t ex_spacing = 0.05;
	Double_t z_lb = -45.0;
	Double_t z_ub = -10.0;
	Double_t ex_lb = 0.0;
	Double_t ex_ub = 4.4;
	
	Int_t num_z = (Int_t)((z_ub - z_lb)/z_spacing);
	Int_t num_ex = (Int_t)((ex_ub - ex_lb)/ex_spacing);

	Double_t *z = new Double_t[ num_z*num_ex ];
	Double_t *ex = new Double_t[ num_z*num_ex ];
	Double_t *theta = new Double_t[ num_z*num_ex ];

	TGraph2D *g = new TGraph2D();
	gStyle->SetPalette(1);
	g->SetTitle("SURFACE; Ex (MeV); z (cm); #theta_{cm} (#circ)");
	
	
	for ( Int_t i = 0; i < num_z; i++ ){
		for ( Int_t j = 0; j < num_ex; j++ ){
			//ex[num_ex*i + j] = j*ex_spacing + ex_lb;
			//z[num_ex*i + j] = i*z_spacing + z_lb;
			/*theta[num_ex*i + j] = */CMAngleCalculator( ex[j], z[i]);
			//g->SetPoint(num_ex*i + j, ex[num_ex*i + j], z[num_ex*i + j], theta[num_ex*i + j] );
		}
	}
	//TCanvas *c = new TCanvas("c","CANVAS", 1200, 900 );
	//g->Draw("surf1");


	delete[] z;
	delete[] theta;
	return;

}


void CMAngleGraph1D( ){
	// Define the z's
	Double_t z_spacing = 0.005;
	Double_t z_lb = -50.0;
	Double_t z_ub = -10.0;
	Double_t ex = 0.582;
	
	Int_t num_z = (Int_t)((z_ub - z_lb)/z_spacing);

	Double_t *z = new Double_t[ num_z ];
	Double_t *theta = new Double_t[ num_z ];

	TGraph *g = new TGraph();
	gStyle->SetPalette(1);
	g->SetTitle( Form( "GRAPH: Ex = %2.1f; z (cm); #theta_{cm} (#circ)", ex ) );
	
	
	for ( Int_t i = 0; i < num_z; i++ ){
		z[i] = i*z_spacing + z_lb;
		theta[i] = CMAngleCalculator( ex, z[i]);
		g->SetPoint(i, z[i], theta[i] );
	}
	TCanvas *c = new TCanvas("c","CANVAS", 1200, 900 );
	g->Draw("AL");


	delete[] z;
	delete[] theta;
	return;

}








