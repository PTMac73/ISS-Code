// writespe function here - adapted from gatemat2.cpp
// ============================================================================================= //
// Patrick MacGregor
// Nuclear Physics Research Group
// Department of Physics and Astronomy
// The University of Manchester
// ============================================================================================= //
#ifndef WRITESPE_h_
#define WRITESPE_h_

#include <TH1.h>
#include <TROOT.h>

void WriteSPE(const Char_t *hisname, Char_t *spename, Char_t const*xy="X"){
	TH1		 *hist;
	Int_t	 i,j,NN,size;
	Int_t	 i1,i2;
	Int_t	 NBINS;
	Char_t	str[256];
	float *sp;
	NN = 4096;
	if ( strncmp( xy, "Y", 1 ) == 0 || strncmp( xy, "y", 1 ) == 0 ){ NN = 4096; }
	FILE *out;
	if ( !( sp = (float*) malloc( NN*sizeof( float ) ) ) ) {
		printf("\007	ERROR: Could not malloc data buffer.\n");
		exit(-1);
	}
	hist = (TH1*)gROOT->FindObject(hisname);
	
	if ( hist != NULL ){
		
		NBINS = hist->GetNbinsX();
		for( i = 1; i < NN + 1; i++ ){
			if ( i <= NBINS ){ 
				sp[i-1] = hist->GetBinContent(i);
			}
			else{
				sp[i-1] = 0.0;
			}
		}

		sprintf( str, "%s.spe", spename );
		
		out = fopen( str, "wb" );
		i = 1;
		j = 24;
		fwrite( &j, 4, 1, out );
		fwrite( str, 8, 1, out );
		fwrite( &NN, 4, 1, out );
		fwrite( &i, 4, 1, out );
		fwrite( &i, 4, 1, out );
		fwrite( &i, 4, 1, out );
		fwrite( &j, 4, 1, out );
		size = sizeof(float)*NN;
		fwrite( &size, 4, 1, out );
		fwrite( sp, 4, NN, out );
		fwrite( &size, 4, 1, out );
		fclose( out );
		
		printf("Wrote %i channels to %s\n", NN, str);
	}
	else {
		printf("Spectrum %s not found\n", hisname);
	}
	//free(sp);
	return;
}


#endif
