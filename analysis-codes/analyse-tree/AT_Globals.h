// AT_Globals.h
// Global functions for the AnalyseTree.C script
// ============================================================================================= //
// Patrick MacGregor
// Nuclear Physics Research Group
// Department of Physics and Astronomy
// The University of Manchester
// ============================================================================================= //

#ifndef AT_GLOBALS_H_
#define AT_GLOBALS_H_

#include "AT_Settings.h"


Int_t col_spacing[4] = { 20, 10, 10, 10 };

// N.B. Do width and then the column of that width

void PrintHorzDiv(){
	std::cout << std::left << std::setfill('-') << \
		std::setw(col_spacing[0])  << "+"  << \
		std::setw(col_spacing[1])  << "+"  << \
		std::setw(col_spacing[2])  << "+"  << \
		std::setw(col_spacing[3])  << "+"  << "+\n";
}

TString AppendColDiv( TString a ){
	return ( "| " + a );
}

TString BoolToStr( Bool_t a ){
	if ( a == 0 ){ return "-"; }
	else { return "x"; };
}

void PrintColumn( TString col1, TString col2, TString col3, TString col4 ){
	col1 = AppendColDiv(col1);
	col2 = AppendColDiv(col2);
	col3 = AppendColDiv(col3);
	col4 = AppendColDiv(col4);


	std::cout << std::left << std::setfill(' ') << \
		std::setw( col_spacing[0] ) << col1 << \
		std::setw( col_spacing[1] ) << col2 << \
		std::setw( col_spacing[2] ) << col3 << \
		std::setw( col_spacing[3] ) << col4 << "|\n";

}

void PrintSingleOption( Bool_t opt, TString text ){
	std::cout << " * " << std::left << std::setfill('.') << std::setw(46) << text << " " << ( opt == 1 ? "Y" : "." ) << "\n";
}

void PrintSummaryOfOptions(){
	// PRINT HEADER
	PrintHorzDiv();
	PrintColumn( "NAME", "ON", "PRINT", "SPE" );
	PrintHorzDiv();
	
	// PRINT CONTENTS
	PrintColumn( "EVZ"        , BoolToStr( SW_EVZ[0]         ) , BoolToStr( SW_EVZ[1]         ), BoolToStr( SW_EVZ[2]         ) );
	PrintColumn( "EVZ_compare", BoolToStr( SW_EVZ_COMPARE[0] ) , BoolToStr( SW_EVZ_COMPARE[1] ), BoolToStr( SW_EVZ_COMPARE[2] ) );
	PrintColumn( "EVZ_Si"     , BoolToStr( SW_EVZ_SI[0]      ) , BoolToStr( SW_EVZ_SI[1]      ), BoolToStr( SW_EVZ_SI[2]      ) );
	PrintColumn( "Ex"         , BoolToStr( SW_EX[0]          ) , BoolToStr( SW_EX[1]          ), BoolToStr( SW_EX[2]          ) );
	PrintColumn( "Ex_compare" , BoolToStr( SW_EX_COMPARE[0]  ) , BoolToStr( SW_EX_COMPARE[1]  ), BoolToStr( SW_EX_COMPARE[2]  ) );
	PrintColumn( "Ex_Si"      , BoolToStr( SW_EX_SI[0]       ) , BoolToStr( SW_EX_SI[1]       ), BoolToStr( SW_EX_SI[2]       ) );
	PrintColumn( "RDT_cuts"   , BoolToStr( SW_RDT_CUTS[0]    ) , BoolToStr( SW_RDT_CUTS[1]    ), BoolToStr( SW_RDT_CUTS[2]    ) );
	PrintColumn( "XCAL"       , BoolToStr( SW_XCAL[0]        ) , BoolToStr( SW_XCAL[1]        ), BoolToStr( SW_XCAL[2]        ) );
	PrintColumn( "XN:XF"      , BoolToStr( SW_XNXF[0]        ) , BoolToStr( SW_XNXF[1]        ), BoolToStr( SW_XNXF[2]        ) );
	
	// PRINT DIVIDER
	PrintHorzDiv();
	PrintSingleOption( PRINT_PDF, "Printing pdfs?" );
	PrintSingleOption( PRINT_PNG, "Printing pngs?" );
	PrintSingleOption( PRINT_ROOT, "Producing ROOT files?" );
	PrintSingleOption( CANVAS_COMBINE, "Combining canvases?" );
	PrintSingleOption( DRAW_NEW_CUTS, "Drawing new cuts?" );
	PrintSingleOption( ALL_ROWS, "Drawing full Ex spectrum?" );
	PrintSingleOption( ROW_BY_ROW, "Drawing RBR Ex spectrum?" );
	PrintSingleOption( DET_BY_DET, "Drawing DBD Ex spectrum?" );


	
	// RESET COUT STREAM
	std::cout << std::right << std::setfill(' ');
		
		

	return;
}








#endif
