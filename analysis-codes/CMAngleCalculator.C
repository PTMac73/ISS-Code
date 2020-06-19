// Calculates the CM angle for the reaction 28Mg(d,p)29Mg in inverse kinematics given an array 
// position and an excitation energy
// ============================================================================================= //
// Patrick MacGregor
// Nuclear Physics Research Group
// School of Physics and Astronomy
// The University of Manchester
// ============================================================================================= //
#include <iostream>
#include <fstream>
#include <vector>
#include <TMath.h>
#include <TString.h>
#include <TGraph.h>
#include <TGraph2D.h>

// CONSTANTS
Double_t c = 299792458;					// Speed of light in m /s
Double_t u = 931.494;					// Unified atomic mass unit in MeV / c^2

// Calculates the average radius of the ISS array (be consistent with units!!)
/*/ The X-axis has its origin in the centre of the array (not necessarily the centre of the strip)
 *  X1 is the point furthest left of the origin
 *  X2 is the point furthest right of the origin
 *  height is the height of the array                                                            */
Double_t ISSArrayRadius( Double_t X1, Double_t X2, Double_t height ){
	Double_t A1 = X1/height;
	Double_t A2 = X2/height;
	return ( height*height/( 2*( X2 - X1 ) ) )*( A2*TMath::Sqrt( 1 + A2*A2 ) - A1*TMath::Sqrt( 1 + A1*A1 ) + TMath::Log( TMath::Abs( A2 + TMath::Sqrt( 1 + A2*A2 ) ) ) - TMath::Log( TMath::Abs( A1 + TMath::Sqrt( 1 + A1*A1 ) ) ) );
}

// VARIABLES
struct ReactionParameters{
	Double_t q;
	Double_t B;
	Double_t m1;
	Double_t m2;
	Double_t m3;
	Double_t m4;
	Double_t bepn;
};

Double_t MassExcessToMass( Double_t Delta, Int_t A, Int_t q = 0 ){
	return A*u + Delta - q*0.511;
}

// FUNCTIONS TO DEFINE VARIABLES
void Set28MgDP( ReactionParameters &rp ){
	rp.q = 1;										// Charge of outgoing proton in units of e
	rp.B = 2.5*c*(1e-8);							// Convert B-field to MeV / (cm*c*e)
	rp.m1 = MassExcessToMass( -15.0188, 28,  9 ); 	// m1 = 28Mg in MeV / c^2		(beam)
	rp.m2 = MassExcessToMass(  13.1357,  2,  1 );	// m2 = d in MeV / c^2			(target)
	rp.m3 = MassExcessToMass(  7.28890,  1,  1 );	// m3 = p in MeV / c^2			(ejectile)
	rp.m4 = MassExcessToMass( -10.6028, 29, 12 );	// m4 = 29Mg in MeV / c^2		(recoil)
	rp.bepn = 9.473;								// bepn = beam energy per nucleon in MeV/u
	return;
}


void Set28SiDP( ReactionParameters &rp ){
	rp.q = 1;										// Charge of outgoing proton in units of e
	rp.B = 2.5*c*(1e-8);							// Convert B-field to MeV / (cm*c*e)
	rp.m1 = MassExcessToMass( -21.4927, 28,  9 ); 	// m1 = 28Mg in MeV / c^2		(beam)
	rp.m2 = MassExcessToMass(  13.1357,  2,  1 );	// m2 = d in MeV / c^2			(target)
	rp.m3 = MassExcessToMass(  7.28890,  1,  1 );	// m3 = p in MeV / c^2			(ejectile)
	rp.m4 = MassExcessToMass( -21.8950, 29, 14 );	// m4 = 29Mg in MeV / c^2		(recoil)
	rp.bepn = 9.473;								// bepn = beam energy per nucleon in MeV/u
	return;
}


// SWITCH
// 0 = d(28Mg,p)29Mg
// 1 = d(28Si,p)29Si
Int_t reaction = 0;
const Int_t NUM_REACTIONS = 2;

TString reaction_string_arr[2] = {
	"d(28Mg,p)29Mg",
	"d(28Si,p)29Si"
};

void PrintReaction(){
	for ( Int_t i = 0; i < NUM_REACTIONS; i++ ){
		std::cout << ( reaction == i ? "> " : "  " ) << Form( "%02d ", i ) << reaction_string_arr[i] << "\n";
	}
	return;
}

void ChooseReaction(Int_t a){
	if ( a >= 0 && a < NUM_REACTIONS ){
		std::cout << "Reaction changed from " << reaction_string_arr[reaction] << " to " << reaction_string_arr[a] << "\n";
		reaction = a;
	}
	else{
		std::cout << "There are only " << NUM_REACTIONS << " available. Please choose a number between 0 and " << NUM_REACTIONS - 1 << ".\n";
		PrintReaction();
	}
	return;
}

// Fill structure with values
void FillRP( ReactionParameters &rp ){
	if ( reaction == 0 ){ Set28MgDP( rp ); }
	else if ( reaction == 1 ){ Set28SiDP( rp ); }
	else{ std::exit(1); }
	return;
}


struct KinematicsParameters{
	Double_t T1; 			// Kinetic energy of the beam [LAB] [ 9.473 MeV / u ]
	Double_t e1; 			// Total energy of beam particle [LAB] in MeV
	Double_t etot;			// Total energy [LAB] in MeV
	Double_t etot_cm ;		// Total energy [CM] (related to invariant mass) in MeV
	Double_t gamLAB_CM; 	// Gamma factor relating inertial frames (i.e. between LAB and CM)
	Double_t beta;			// Beta = v/c (v = velocity of CM frame)
	Double_t rho; 			// Radius of array in cm (approximate as circular - diameter of square is 23 mm).
	Double_t m4ex; 			// Invariant mass of the recoil nucleus
	Double_t e3_cm; 		// Energy of the outgoing proton [CM]
	Double_t p3_cm; 		// Momentum of the outgoing proton [CM]
};

void CalculateKP( ReactionParameters &rp, KinematicsParameters &kp, Double_t ex ){
	kp.T1 = rp.bepn*rp.m1/u;																// Kinetic energy of the beam [LAB] [ 9.473 MeV / u ]
	kp.e1 = kp.T1 + rp.m1;																	// Total energy of beam particle [LAB] in MeV
	kp.etot = kp.e1 + rp.m2;																// Total energy [LAB] in MeV
	kp.etot_cm = TMath::Sqrt( rp.m1*rp.m1 + rp.m2*rp.m2 + 2*kp.e1*rp.m2 );					// Total energy [CM] (related to invariant mass) in MeV
	kp.gamLAB_CM = kp.etot/kp.etot_cm;														// Gamma factor relating inertial frames (i.e. between LAB and CM)
	kp.beta = TMath::Sqrt( 1 - TMath::Power( kp.gamLAB_CM, -2 ) );							// Beta = v/c (v = velocity of CM frame)
	kp.rho = ISSArrayRadius( -0.45, 0.45, 1.15 );											// Radius of array in cm (approximate as circular - diameter of square is 23 mm).
	kp.m4ex = rp.m4 + ex;																	// Invariant mass of the recoil nucleus
	kp.e3_cm = 0.5*( kp.etot_cm*kp.etot_cm + rp.m3*rp.m3 - kp.m4ex*kp.m4ex )/kp.etot_cm;	// Energy of the outgoing proton [CM]
	kp.p3_cm = TMath::Sqrt(kp.e3_cm*kp.e3_cm - rp.m3*rp.m3);								// Momentum of the outgoing proton [CM]
	return;
};






// MAIN FUNCTION
// ex in MeV, z in cm
// In order for this to work, quote:
//	* Masses, momenta, and energies in terms of MeV and factors of c
//	* Velocities in terms of c
//	* Charges in terms of e
//	* Distances in terms of cm 
Double_t CMAngleCalculator(Double_t ex, Double_t z, Bool_t print = 1 ){
	// Define initial fixed quantities for the reaction
	ReactionParameters rp;
	FillRP(rp);

	// Now calculate derived quantities from fixed quantities
	KinematicsParameters kp;
	CalculateKP( rp, kp, ex );

	// Function-specific calculations
	Double_t p_para = rp.q*rp.B*z/( 2*TMath::Pi() );											// Parallel component of the momentum in MeV / c [LAB] -------> THIS IS A GUESS!!
	Double_t p_para_cm = p_para/kp.gamLAB_CM - kp.beta*kp.e3_cm;								// Parallel component of the momentum in MeV / c [CM]
	Double_t p_perp_cm = TMath::Sqrt( kp.e3_cm*kp.e3_cm - p_para_cm*p_para_cm - rp.m3*rp.m3 );	// Perpendicular component of the momentum in MeV / c [CM]

	// Calculate the function to be minimised and its first derivative
	Double_t fp = ( 2*p_perp_cm/( rp.q*rp.B ) )*TMath::Sin( rp.q*rp.B*z/( 2*kp.gamLAB_CM*( p_para_cm + kp.beta*kp.e3_cm ) ) ) - kp.rho;
	Double_t fd = ( ( -2*p_para_cm )/( rp.q*rp.B*p_perp_cm ) )*TMath::Sin( rp.q*rp.B*z/( 2*kp.gamLAB_CM*( p_para_cm + kp.beta*kp.e3_cm ) ) ) - ( ( p_perp_cm*z )/( kp.gamLAB_CM*TMath::Power( p_para_cm + kp.beta*kp.e3_cm , 2 ) ) )*TMath::Cos( rp.q*rp.B*z/( 2*kp.gamLAB_CM*( p_para_cm + kp.beta*kp.e3_cm ) ) );

	// Now perform Newton-Raphson calculation to get the angle in the CM frame
	Int_t n = 0;
	while ( TMath::Abs(fp) > 1e-5 && n < 100000){
		p_para_cm = p_para_cm - fp/fd;
		p_perp_cm = TMath::Sqrt( kp.e3_cm*kp.e3_cm - p_para_cm*p_para_cm - rp.m3*rp.m3 );

		fp = ( 2*p_perp_cm/( rp.q*rp.B ) )*TMath::Sin( rp.q*rp.B*z/( 2*kp.gamLAB_CM*( p_para_cm + kp.beta*kp.e3_cm ) ) ) - kp.rho;
		fd = ( ( -2*p_para_cm )/( rp.q*rp.B*p_perp_cm ) )*TMath::Sin( rp.q*rp.B*z/( 2*kp.gamLAB_CM*( p_para_cm + kp.beta*kp.e3_cm ) ) ) - ( ( p_perp_cm*z )/( kp.gamLAB_CM*TMath::Power( p_para_cm + kp.beta*kp.e3_cm , 2 ) ) )*TMath::Cos( rp.q*rp.B*z/( 2*kp.gamLAB_CM*( p_para_cm + kp.beta*kp.e3_cm ) ) );
		n++;
	}

	// Function minimised - now calculate the final quantities
	Double_t theta_cm = 180 - ( TMath::ACos( p_para_cm/kp.p3_cm )*180/TMath::Pi() );
	Double_t e3 = kp.gamLAB_CM*( kp.e3_cm + kp.beta*p_para_cm );
	Double_t p3 = TMath::Sqrt( e3*e3 - rp.m3*rp.m3 );
	Double_t p_perp = p_perp_cm;
	p_para = TMath::Sqrt( p3*p3 - p_perp*p_perp );
	Double_t theta = TMath::ACos( kp.gamLAB_CM*( p_para_cm + kp.beta*kp.e3_cm )/p3 )*180/TMath::Pi();
	if ( print == 1 ){ std::cout << std::left << std::setprecision(12) << std::setw(12) << ex << "\t" << std::setprecision(12) << std::setw(12) << z << "\t" << std::setprecision(12) << std::setw(12) << theta_cm << std::endl; }
	//std::cout << n << "\t" << ex << "\t" << z << "\t" << theta_cm << "\t" << theta << "\t" << e3-m3 << "\t" << std::endl;
	//std::cout << p3 << "\t" << p_para << "\t" << p_perp << "\t" << p3_cm << "\t" << p_para_cm << "\t" << p_perp_cm << "\t" << p3_cm*TMath::Sin(theta_cm) << "\t" << p3*TMath::Sin(theta) << std::endl;

	return theta_cm;
}

/*************************************************************************************************/
void CMAngleBatch( TString in_file_Dir, Bool_t b_draw = 0 ){
	// Open the file for reading	
	ifstream in_file;
	in_file.open( in_file_Dir.Data() );

	// Check that it opened
	if ( in_file.is_open() ){
		// Print success message
		std::cout << "File opened successfully!" << std::endl;

		// Read the data
		Double_t ex, z, theta;
		std::string line;
		TString ex_str, z_str;
		Int_t num_lines = 0;
		std::cout << std::left << std::setw(12) << "EX" << "\t" << std::setw(12) << "z" << "\t" << std::setw(12) << "TH_CM" << std::endl;
		
		while( getline( in_file, line ) ){
			TObjArray* obj_arr = ( (TString)line ).Tokenize("\t");
			
			if ( obj_arr->GetEntries() >= 2 ){
				ex_str = ( (TObjString*)( obj_arr->At(0) ) )->String();
				z_str = ( (TObjString*)( obj_arr->At(1) ) )->String();
				
				if ( ex_str.IsFloat() && z_str.IsFloat() ){
					ex = ex_str.Atof();
					z = z_str.Atof();
					theta = CMAngleCalculator( ex, z );
				}
				else{
					std::cout << std::left << std::setw(12) << "--" << "\t" << std::setw(12) << "--" << "\t" << std::setw(12) << "--" << "\n";
				}
				
				num_lines++;
			}
		}
		
		
		
/*
		while ( !in_file.eof() ){
			// Store the values
			in_file >> ex >> z;

			num_lines++;
			
			// Calculate the angles
			theta = CMAngleCalculator( ex, z );

		}
*/
		// Close the file
		in_file.close();
		std::cout << "File closed (I think...)!" << std::endl;

		// Draw it if desired
		if ( b_draw == 1 ){
			const Int_t arr_size = num_lines;
			Double_t* arr_z = new Double_t[arr_size];
			Double_t* arr_th = new Double_t[arr_size];
			in_file.open( in_file_Dir.Data() );
			Int_t i = 0;
			while ( !in_file.eof() ){
				in_file >> ex >> arr_z[i];
				arr_th[i] = CMAngleCalculator( ex, arr_z[i], 0 );
				i++;
			}
			in_file.close();
	

			TGraph *g = new TGraph( arr_size, arr_z, arr_th);
			TCanvas *c = new TCanvas( "canvas", "canvas", 1200, 900 );
			g->Draw();


			delete[] arr_z;
			delete[] arr_th;
		}

	}
	else{
		// Print failure message
		std::cout << "File did not open. Process terminated." << std::endl;
	}
}

/*************************************************************************************************/
// Try and find an angle for a given theta and excitation
Double_t CMAngleFindZ( Double_t ex, Double_t theta_cm ){
	// Define initial fixed quantities for the reaction
	ReactionParameters rp;
	FillRP(rp);

	// Now calculate derived quantities from fixed quantities
	KinematicsParameters kp;
	CalculateKP( rp, kp, ex );

	Double_t theta_cm_true = 180 - theta_cm;
	Double_t p_para_cm = kp.p3_cm*TMath::Cos( theta_cm_true*TMath::DegToRad() );	// Parallel component of the momentum in MeV / c [CM]
	Double_t p_perp = kp.p3_cm*TMath::Sin( theta_cm_true*TMath::DegToRad() );		// Perpendicular component of the momentum in MeV / c [CM/LAB]
	Double_t p_para_lab = kp.gamLAB_CM*(p_para_cm + kp.beta*kp.e3_cm);				// Parallel component of the momentum in MeV / c [LAB]
	
	Double_t r = TMath::Abs( p_perp/(rp.q*rp.B) );						// Radius or orbit / cm [LAB]
	Double_t lc = 2*r*( TMath::Pi() - TMath::ASin( 0.5*kp.rho/r ) );	// Distance travelled perp to z / cm [LAB]
	Double_t z = p_para_lab*lc/p_perp;									// z / cm [LAB] --> This is a guess!

	// Newton-Raphson time!
	Int_t n = 0;
	Double_t fp = 2*r*TMath::Sin( p_perp*z/( 2*r*p_para_lab ) ) - kp.rho;
	Double_t fd = p_perp*TMath::Cos( p_perp*z/( 2*r*p_para_lab ) )/p_para_lab;
	while ( n < 10000 && TMath::Abs( fp ) > 1e-5 ){
		z = z - fp/fd;
		fp = 2*r*TMath::Sin( p_perp*z/( 2*r*p_para_lab ) ) - kp.rho;
		fd = p_perp*TMath::Cos( p_perp*z/( 2*r*p_para_lab ) )/p_para_lab;
		n++;
	}

	return z;
}

/*************************************************************************************************/
void CMAngleFindMinZ( Double_t theta ){
	const Int_t NUM_STATES = 10;
	Double_t STATES[NUM_STATES] = {
		0.00000,
		0.05460,
		1.09860,
		1.43186,
		2.25419,
		2.48316,
		2.88716,
		3.18338,
		3.91978,
		4.28368
	};

	Int_t n;
	Double_t step_size;

	for ( Int_t i = 0; i < NUM_STATES; i++ ){
		n = 0;
		Double_t step_size = 0.1;
		Double_t z = CMAngleFindZ( STATES[i], theta );
		Double_t test;
		Double_t min, min_z;

		// Get initial estimate of angle
		while ( step_size > 1e-8 ){
			min = 5;
			min_z = z;
			for ( Int_t j = -10; j <= 10; j++ ){
				test = CMAngleCalculator( STATES[i], z + j*step_size, 0 );
				if ( TMath::Abs( theta - test ) < min ){
					min = TMath::Abs( theta - test );
					min_z = z + j*step_size;
				}
			}
			z = min_z;
			step_size /= 10.0;
		}
		std::cout << STATES[i] << "\t" << z << "\t" << std::setprecision(12) << CMAngleCalculator( STATES[i], z, 0 ) << "\n";
	}

}


/*************************************************************************************************/
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

/*************************************************************************************************/
void CMAngleGraph1D( ){
	// Define the z's
	Double_t z_spacing = 0.005;
	Double_t z_lb = -50.0;
	Double_t z_ub = -10.0;
	Double_t ex = 0.0;
	
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









