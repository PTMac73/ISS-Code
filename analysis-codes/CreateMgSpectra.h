// CreateMgSpectra.h
// Header file containing lots of useful quantities and functions for creating Mg spectra (split or
// not)
// ============================================================================================= //
// Patrick MacGregor
// Nuclear Physics Research Group
// Department of Physics and Astronomy
// The University of Manchester
// ============================================================================================= //

#ifndef CREATE_MG_SPECTRA_H_
#define CREATE_MG_SPECTRA_H_

// GLOBAL CONSTANTS ---------------------------------------------------------------------------- //
TString print_dir = "/home/ptmac/Documents/07-CERN-ISS-Mg/Mg-Analysis/SPE-Files";
Int_t CANVAS_SCALE = 300;
Int_t CANVAS_WIDTH = 4*CANVAS_SCALE;
Int_t CANVAS_HEIGHT = 3*CANVAS_SCALE;

// Define properties of the system
const Int_t NUM_ROWS = 6;
const Int_t NUM_DETECTORS_PER_ROW = 4;
const Int_t NUM_DETECTORS = NUM_ROWS*NUM_DETECTORS_PER_ROW;
const Int_t NUM_STRIPS_PER_SI = 1;

// Cut value
const Double_t THETA_CM_LIMIT[2] = {
	16.6278,
	11.0000
};

// CHOOSE WHICH DETECTORS TO USE --------------------------------------------------------------- //
// * 0 = All detectors
// * 1 = Best detectors
const Int_t DETECTOR_MODE = 0;

// CHOOSE HOW TO LAYOUT SPECTRA ---------------------------------------------------------------- //
// * 0 = Row-by-row
// * 1 = Detector-by-detector
// * 2 = Full excitation spectrum
const Int_t SPECTRUM_LAYOUT = 0;

// CHOOSE WHICH CUTS TO USE -------------------------------------------------------------------- //
// * 0 = Full cuts
// * 1 = Cuts for singles spectra
const Int_t CUTS_MODE = 1;

// This determines the detectors to use if all the detectors are being used. The split is for
// between states 1-6 and 7-9.
Bool_t use_detector_1[24][2] = {
	{ 1, 1 },  // DET #0
	{ 1, 1 },  // DET #1
	{ 0, 0 },  // DET #2
	{ 1, 1 },  // DET #3
	{ 1, 1 },  // DET #4
	{ 1, 1 },  // DET #5
	{ 1, 1 },  // DET #6
	{ 0, 0 },  // DET #7
	{ 0, 0 },  // DET #8
	{ 1, 1 },  // DET #9
	{ 1, 1 },  // DET #10
	{ 1, 1 },  // DET #11
	{ 1, 0 },  // DET #12
	{ 0, 0 },  // DET #13
	{ 1, 0 },  // DET #14
	{ 0, 0 },  // DET #15
	{ 1, 0 },  // DET #16
	{ 0, 0 },  // DET #17
	{ 1, 1 },  // DET #18
	{ 1, 1 },  // DET #19
	{ 1, 1 },  // DET #20
	{ 1, 1 },  // DET #21
	{ 1, 1 },  // DET #22
	{ 1, 1 }   // DET #23
};

// This determines the detectors to use if the best resolution detectors are being used.
Bool_t use_detector_2[24][1] = {
	{ 1 },  // DET #0
	{ 1 },  // DET #1
	{ 0 },  // DET #2
	{ 1 },  // DET #3
	{ 0 },  // DET #4
	{ 0 },  // DET #5
	{ 0 },  // DET #6
	{ 0 },  // DET #7
	{ 1 },  // DET #8
	{ 0 },  // DET #9
	{ 1 },  // DET #10
	{ 0 },  // DET #11
	{ 0 },  // DET #12
	{ 0 },  // DET #13
	{ 0 },  // DET #14
	{ 0 },  // DET #15
	{ 0 },  // DET #16
	{ 0 },  // DET #17
	{ 0 },  // DET #18
	{ 0 },  // DET #19
	{ 0 },  // DET #20
	{ 0 },  // DET #21
	{ 0 },  // DET #22
	{ 1 }   // DET #23
};


// This labels the mode of operation
TString label_1[2] = {
	"states1-6",
	"states7-9"
};

TString label_2[1] = { "best_res" };

// ID string array for different spectrum types
TString SPEC_TYPE_ID[4] = {
	"",			// 0
	"row",		// 1
	"det",		// 2
	"full",		// 3
};

TString CUTS_MODE_ID[2] = {
	"all",
	"sing"
};



// MODE OF OPERATION FUNCTIONS + VARIABLES ----------------------------------------------------- //
// This tells you how large the use_detector arrays are based on the mode
Int_t GetNumSpectrumVariants( Int_t det_mode ){
	Int_t a = 0;
	if ( det_mode == 0 ){ a = 2; }
	else if ( det_mode == 1 ){ a = 1; }
	return a;
}


// FUNCTIONS ----------------------------------------------------------------------------------- //
// Makes the names for an object depending on whether it's a histogram or a canvas
// N.B. ID number can be detID or row number
TString MakeObjectName( Int_t type, Int_t spec_layout, Int_t id_number, Int_t strip_number, Int_t cuts_mode ){
	// TYPE 0 --> histogram
	// TYPE 1 --> canvas
	Char_t type_char;
	if ( type == 0 ){ type_char = 'h'; }
	else if ( type == 1 ){ type_char = 'c'; }
	else { std::cout << "Unknown type declared." << std::endl; exit(1); }

	TString name = TString::Format( "%c_spec%i_%i_%i_%i", type_char, spec_layout, id_number, strip_number, cuts_mode );
	return name;
}



// Makes the file names for all of the printed spectra
TString MakePrintFileName( Int_t pos_number, TString LABEL , Int_t id_number, Int_t strip_number, Int_t num_strips, TString file_type, Int_t spec_layout, Int_t cuts_mode ){
	TString name;
	TString det_identifier = SPEC_TYPE_ID[spec_layout+1];
	TString cut_identifier = CUTS_MODE_ID[cuts_mode];

	if ( strip_number == 0 && num_strips == 1 ){
		// Full spectrum
		if ( spec_layout == 3 ){
			name = TString::Format( "pos%i_%s_%s-cuts%s", pos_number, det_identifier.Data(), cut_identifier.Data(), file_type.Data() );
		}
		// Rows/Dets
		else{
			name = TString::Format( "pos%i_%s%i_%s_%s-cuts%s", pos_number, det_identifier.Data(), id_number, LABEL.Data(), cut_identifier.Data(), file_type.Data() );
		}
	}
	else{
		name = TString::Format( "pos%i_%s%i_strip%i-of-%i_%s_%s-cuts%s", pos_number, det_identifier.Data(), id_number, strip_number+1, num_strips, LABEL.Data(), cut_identifier.Data(), file_type.Data() );
	}
	return name;
}


TString MakeRootFileName( Int_t pos_number, Int_t spec_layout, Int_t cuts_mode ){
	TString det_identifier = SPEC_TYPE_ID[spec_layout+1];
	TString cut_identifier = CUTS_MODE_ID[cuts_mode];
	TString name = TString::Format( "%s/pos%i_%s_%s.root", print_dir.Data(), pos_number, cut_identifier.Data(), det_identifier.Data() );
	return name;
	
}


// Function to spit out the thetaCM cuts **OBSOLETE**
TString GenerateThetaCMCutString( Int_t NUM_STATES = 9 ){
	TString out_str = "( ";
	for ( Int_t i = 0; i < NUM_STATES; i++ ){
		out_str += Form( "( Ex > ex_lims[%i] && Ex <= ex_lims[%i] && thetaCM >= thetaCM_lims[%i] )", i, i+1, i );
		if ( i < NUM_STATES - 1 ){
			out_str += " || ";
		}
	}
	out_str += " )";
	return out_str;
}



// Write a function to spit out the xcal cuts
TString GenerateXCALCutString( Int_t NUM_STRIPS_PER_SI, Int_t strip_number ){
	// Define a string
	TString out_str;

	// Calculate what the out string is
	if ( strip_number == 0 && NUM_STRIPS_PER_SI == 1){
		// No strips, just the regular thing - add the excluded middle as well
		out_str = "xcal[] >= xcal_cuts[][0] && xcal[] <= xcal_cuts[][1] && ( xcal[] <= xcal_cuts[][2] || xcal[] >= xcal_cuts[][3] )";
	}
	else if ( strip_number < NUM_STRIPS_PER_SI ){
		out_str = Form( "xcal[] >= ( xcal_cuts[][0] + ( %i.0*( xcal_cuts[][1] - xcal_cuts[][0]) )/%i.0 ) &&  xcal[] <= ( xcal_cuts[][0] + ( %i.0*( xcal_cuts[][1] - xcal_cuts[][0]) )/%i.0 ) && ( xcal[] <= xcal_cuts[][2] || xcal[] >= xcal_cuts[][3] )", strip_number, NUM_STRIPS_PER_SI, strip_number + 1, NUM_STRIPS_PER_SI );
	}
	else{
		std::cout << "Warning. Strip number exceeds total number of strips. Examine code closely!.\n";
		exit(1);
	}

	// Return the string
	return out_str;
}









#endif
