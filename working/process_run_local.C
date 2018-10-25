TString expName = "iss631";

void process_run_local(Int_t RUNNUM=5, Int_t SORTNUM=0)
{
  if (SORTNUM==0) {
    TString name;
    name.Form("/Users/heliosdigios/experiments/%s/root_data/run", expName.Data());
    name+=RUNNUM;
    name+=".root";
    TFile f(name);
    TTree *t1 = (TTree*)f.Get("tree");

    TString processCmd;
    processCmd.Form("/Users/heliosdigios/experiments/%s/analysis/sort_codes/GeneralSort.C++", expName.Data());

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
