// Calculates the CM angle for the reaction 28Mg(d,p)29Mg in inverse kinematics given an array position and an excitation energy
#include <iostream>
#include <fstream>
#include "TMath.h"
#include "TString.h"

TString iFileDir = "/home/ptmac/Documents/07-CERN-ISS-Mg/Mg-Analysis/excitationEnergyList.txt";


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
	double T1 = 9.473*m1/u;										// Kinetic energy of the beam [LAB] [ 10 MeV / u ] 265.16
	double e1 = T1 + m1;										// Total energy of beam particle [LAB] in MeV
	double etot = e1 + m2;										// Total energy [LAB] in MeV
	double etot_cm = TMath::Sqrt( m1*m1 + m2*m2 + 2*e1*m2 );	// Total energy [CM] (related to invariant mass) in MeV
	double gam = etot/etot_cm;									// Gamma factor
	double beta = TMath::Sqrt( 1 - TMath::Power( gam, -2 ) );	// Beta = v/c
	double rho = 1.15;											// Radius of array in cm (approximate as circular - diameter of square is 23 mm).

	double m4ex = m4 + ex;													// Invariant mass of the recoil nucleus
	double e3_cm = 0.5*( etot_cm*etot_cm + m3*m3 - m4ex*m4ex )/etot_cm;		// Energy of the outgoing proton [CM]
	double p3_cm = TMath::Sqrt(e3_cm*e3_cm - m3*m3);						// Momentum of the outgoing proton [CM]
	double p_para = q*B*z/( 2*TMath::Pi()*gam );							// Parallel component of the momentum in MeV / c [LAB]
	double p_para_cm = p_para/gam - beta*e3_cm;								// Parallel component of the momentum in MeV / c [CM]
	double p_perp_cm = TMath::Sqrt( e3_cm*e3_cm - p_para_cm*p_para_cm - m3*m3 );		// Perpendicular component of the momentum in MeV / c [CM]

	// Calculate the function to be minimised and its first derivative
	double fp = 2*TMath::Sqrt( p_perp_cm )/( q*B )*TMath::Sin( q*B*z/( 2*gam*( p_para_cm + beta*e3_cm ) ) ) - rho;
	double fd=-(2/q*B)*(TMath::Sqrt(p_perp_cm)*q*B*z*TMath::Cos(q*B*z/(2*gam*(p_para_cm+beta*e3_cm)))/(2*gam*(p_para_cm+beta*e3_cm)*(p_para_cm+beta*e3_cm))-(p_para_cm/TMath::Sqrt(p_perp_cm)*TMath::Sin(q*B*z/(2*gam*(p_para_cm+beta*e3_cm)))));
	
	//double fp = ( ( 2*p_perp_cm )/( q*B ) )*TMath::Sin( ( q*B*z )/( 2*( p_para_cm + beta*e3_cm ) ) ) - rho;






	// Now perform Newton-Raphson calculation to get the angle in the CM frame
	while (TMath::Abs(fp)>1e-5){
			p_para_cm = p_para_cm - fp/fd;
			p_perp_cm = TMath::Sqrt( e3_cm*e3_cm - p_para_cm*p_para_cm - m3*m3 );
			//p_perp_cm = ( e3_cm*e3_cm - p_para_cm*p_para_cm - m3*m3 );

			fp=2*TMath::Sqrt(p_perp_cm)/(q*B)*TMath::Sin(q*B*z/(2*gam*(p_para_cm+beta*e3_cm)))-rho;
			fd=-(2/q*B)*(TMath::Sqrt(p_perp_cm)*q*B*z*TMath::Cos(q*B*z/(2*gam*(p_para_cm+beta*e3_cm)))/(2*gam*(p_para_cm+beta*e3_cm)*(p_para_cm+beta*e3_cm))-(p_para_cm/TMath::Sqrt(p_perp_cm)*TMath::Sin(q*B*z/(2*gam*(p_para_cm+beta*e3_cm)))));
		n++;
	}

	// Function minimised - now calculate the final quantities
	double theta_cm = 180 - ( TMath::ACos( p_para_cm/p3_cm )*180/TMath::Pi() );
	double e3 = gam*( e3_cm + beta*p_para_cm );
	double p3 = TMath::Sqrt( e3*e3 - m3*m3 );
	double theta = TMath::ACos( gam*( p_para_cm + beta*e3_cm )/p3 )*180/TMath::Pi();
	std::cout << n << "\t" << ex << "\t" << z << "\t" << theta_cm << "\t" << theta << "\t" << e3-m3 << "\t" << std::endl;
		
}

/*************************************************************************************************/
void CMAngleBatch(){
	// Open the file for reading	
	ifstream iFile;
	iFile.open( iFileDir.Data() );

	// Check that it opened
	if ( iFile.is_open() ){
		// Print success message
		std::cout << "File opened successfully!" << std::endl;

		// Read the data
		double ex, z;
		std::cout << "n\tEx\tz\tTH_CM\tTH\tE3-m3" << std::endl;
		while ( !iFile.eof() ){
			// Store the values
			iFile >> ex >> z;

			// Calculate the angles
			CMAngleCalculator( ex, z );
		}

		// Close the file
		iFile.close();
	}
	else{
		// Print failure message
		std::cout << "File did not open. Process terminated." << std::endl;
	}
}















