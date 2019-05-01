// Calculates the CM angle for the reaction 28Mg(d,p)29Mg in inverse kinematics given an array position and an excitation energy
#include <iostream>
#include <fstream>
#include "TMath.h"
#include "TString.h"

//TString iFileDir = "/home/ptmac/Documents/07-CERN-ISS-Mg/Mg-Analysis/excitationEnergyList.txt";


// MAIN FUNCTION
// ex in MeV, z in cm
// In order for this to work, quote:
//	* Masses, momenta, and energies in terms of MeV and factors of c
//	* Velocities in terms of c
//	* Charges in terms of e
//	* Distances in terms of cm 
void CMAngleCalculator(double ex, double z){
	// Define initial fixed quantities for the reaction
	int n = 0;								// Number of iterations
	double c = 299792458;					// Speed of light in m /s
	double q = 1;							// Charge of outgoing proton in units of e
	double u = 931.494;						// Unified atomic mass unit in MeV / c^2
	double B = 2.5*c*(1e-8);				// Convert B-field to MeV / (cm*c*e)
	double m1 = 27.9838768*u;				// m1 = 28Mg in MeV / c^2		(beam)
	double m2 = 2.01410177785*u;			// m2 = d in MeV / c^2			(target)
	double m3 = 1.00782503207*u;			// m3 = p in MeV / c^2			(ejectile)
	double m4 = 28.988600*u;				// m4 = 29Mg in MeV / c^2		(recoil)

	// Now calculate derived quantities from fixed quantities
	double T1 = 9.473*m1/u;										// Kinetic energy of the beam [LAB] [ 9.473 MeV / u ]
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













