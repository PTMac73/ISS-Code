TString dir="/home/ptmac/Documents/07-CERN-ISS-Mg/analysis/bufferSizeCheck";
TString expName = "iss000";

void process_run(Int_t RUNNUM=5, Int_t BBNUM=450){
	TString name;
	if ( RUNNUM > 9 & RUNNUM < 100){
		name.Form("%s/root_data/run%02d-%d.root", dir.Data(), RUNNUM, BBNUM);
	}
	else if ( RUNNUM > 100 ){
		name.Form("%s/root_data/run%03d-%d.root", dir.Data(), RUNNUM, BBNUM);
	}
	TFile f(name);
	TTree *t1 = (TTree*)f.Get("tree");

	TString processCmd;
	processCmd.Form("%s/GeneralSort.C+", dir.Data());

	t1->Process(processCmd);
	f.Close();
}

