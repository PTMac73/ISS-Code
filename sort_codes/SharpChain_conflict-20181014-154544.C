#include "TChain.h"
#include "TFile.h"
#include "TTree.h"

void SharpChain(Int_t which = 0){
	// TChain the runs together
	gDirectory->Delete("ch");
	TChain *ch = new TChain("gen_tree");
	TFile *f;
	
	// Check that "which" has been assigned
	if (which == 0){
		ch->AddFile("../root_data/gen_run50.root");
		ch->AddFile("../root_data/gen_run51.root");
		ch->AddFile("../root_data/gen_run52.root");
		ch->AddFile("../root_data/gen_run53.root");
		ch->AddFile("../root_data/gen_run54.root");
		ch->AddFile("../root_data/gen_run55.root");
		ch->AddFile("../root_data/gen_run56.root");
		ch->AddFile("../root_data/gen_run57.root");
		ch->AddFile("../root_data/gen_run58.root");
		ch->AddFile("../root_data/gen_run59.root");
		ch->AddFile("../root_data/gen_run60.root");
		ch->AddFile("../root_data/gen_run61.root");
		ch->AddFile("../root_data/gen_run62.root");
		ch->AddFile("../root_data/gen_run63.root");
		ch->AddFile("../root_data/gen_run64.root");
		ch->AddFile("../root_data/gen_run65.root");
		ch->AddFile("../root_data/gen_run66.root");
		ch->AddFile("../root_data/gen_run67.root");
		ch->AddFile("../root_data/gen_run68.root");
		ch->AddFile("../root_data/gen_run69.root");
		ch->AddFile("../root_data/gen_run70.root");
		ch->AddFile("../root_data/gen_run71.root");
		ch->AddFile("../root_data/gen_run72.root");
		ch->AddFile("../root_data/gen_run73.root");
		ch->AddFile("../root_data/gen_run74.root");
		ch->AddFile("../root_data/gen_run75.root");
		ch->AddFile("../root_data/gen_run76.root");
		ch->AddFile("../root_data/gen_run77.root");
		ch->AddFile("../root_data/gen_run78.root");
		ch->AddFile("../root_data/gen_run79.root");
		ch->AddFile("../root_data/gen_run80.root");
		ch->AddFile("../root_data/gen_run82.root");
		ch->AddFile("../root_data/gen_run83.root");
		ch->AddFile("../root_data/gen_run84.root");
		ch->AddFile("../root_data/gen_run85.root");
		ch->AddFile("../root_data/gen_run87.root");
		ch->AddFile("../root_data/gen_run88.root");
		ch->AddFile("../root_data/gen_run89.root");
		ch->AddFile("../root_data/gen_run90.root");
		ch->AddFile("../root_data/gen_run91.root");
		ch->AddFile("../root_data/gen_run92.root");
		ch->AddFile("../root_data/gen_run93.root");
		ch->AddFile("../root_data/gen_run94.root");
		ch->AddFile("../root_data/gen_run95.root");
		ch->AddFile("../root_data/gen_run96.root");
		ch->AddFile("../root_data/gen_run97.root");
		ch->AddFile("../root_data/gen_run98.root");
		ch->AddFile("../root_data/gen_run99.root");
		ch->AddFile("../root_data/gen_run100.root");
		
		// Write to a file
		f = new TFile("../root_data/gen_runPos1.root","NEW");
	}
	else{
		ch->AddFile("../root_data/gen_run103.root");
		ch->AddFile("../root_data/gen_run104.root");
		ch->AddFile("../root_data/gen_run105.root");
		ch->AddFile("../root_data/gen_run106.root");
		ch->AddFile("../root_data/gen_run107.root");
		ch->AddFile("../root_data/gen_run108.root");
		ch->AddFile("../root_data/gen_run110.root");
		ch->AddFile("../root_data/gen_run111.root");
		ch->AddFile("../root_data/gen_run112.root");
		ch->AddFile("../root_data/gen_run113.root");
		ch->AddFile("../root_data/gen_run115.root");
		ch->AddFile("../root_data/gen_run116.root");
		ch->AddFile("../root_data/gen_run117.root");
		ch->AddFile("../root_data/gen_run118.root");
		// Write to a file
		f = new TFile("../root_data/gen_runPos2.root","NEW");
	}
	
	// Write to a file
	if ( f->IsOpen() ){
		printf("File is open\n");
		ch->Write("gen_tree");
		f->Close();
	}
	
}
