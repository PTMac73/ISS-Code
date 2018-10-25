TString dir="/home/ptmac/Documents/07-CERN-ISS-Mg";
TString expName = "iss000";

void process_run(Int_t RUNNUM=5, Int_t SORTNUM=0){
  if (SORTNUM==0) {
    TString name;
    if ( RUNNUM > 9 & RUNNUM < 100){
    	name.Form("%s/analysis/root_data/run%02d.root", dir.Data(), RUNNUM);
    }
    else if ( RUNNUM > 100 ){
    	name.Form("%s/analysis/root_data/run%03d.root", dir.Data(), RUNNUM);
    }
    TFile f(name);
    TTree *t1 = (TTree*)f.Get("tree");

    TString processCmd;
    processCmd.Form("%s/analysis/sort_codes/GeneralSort.C+", dir.Data(),expName.Data());

    t1->Process(processCmd);
    f.Close();
  }

  else if (SORTNUM==1) {
    TString name("gen.root");
    TFile ff(name);
    TTree *t2 = (TTree*)ff.Get("gen_tree");
    t2->Process("../codes/Monitors.C++");
  }
}
