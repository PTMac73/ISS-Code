// AG_physical_constants.C
// Physical constants for the file ArrayGeometry.C
// ============================================================================================= //
// Patrick MacGregor
// Nuclear Physics Research Group
// School of Physics and Astronomy
// The University of Manchester
// ============================================================================================= //
#ifndef AG_PHYSICAL_CONSTANTS_H_
#define AG_PHYSICAL_CONSTANTS_H_

#include <TColor.h>

// PHYSICAL CONSTANTS
const Double_t AMU = 931.494;		// MeV

// CONSTANTS ABOUT THE SCRIPT
const Int_t CANVAS_WIDTH = 1200;
const Int_t CANVAS_HEIGHT = 900;
const Bool_t BATCH_MODE = 0;

// Set colours
Int_t pcb_green_i = TColor::GetFreeColorIndex();
TColor *pcb_green = new TColor( pcb_green_i, 0.0, 111.0/255.0, 75.0/255.0 );

Int_t si_strip_i = TColor::GetFreeColorIndex();
TColor *si_strip = new TColor( si_strip_i, 141.0/255.0, 116.0/255.0, 71.0/255.0 );

Int_t fj_red_i = TColor::GetFreeColorIndex();
TColor *fj_red = new TColor( fj_red_i, 189.0/255.0, 124.0/255.0, 124.0/255.0 );


#endif
