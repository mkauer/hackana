// This macro will extract the number of generated MC,
// the number of kept events from slimming, the activity,
// and the runtime. It then computes the efficiency and
// the number of expected events.

// This doesn't work for rootfiles that used inplimentation
// of victor's background model or vera's radon model because
// the original number of kept events gets erased and replaced
// with the number of expected events.

// Version 09.09.18


void mcstats(char *rootfile){
  ///////////////////////////
  string name="etot";
  ///////////////////////////
  TFile *file=new TFile(rootfile,"READ");
  TDirectory *d;
  TH1D *hist;
  char ananame[256];
  //char *ananame;
  Long64_t bmcevga;
  Float_t activity;
  Float_t time;
  TTree *tree=(TTree*)file->Get("ana10");
  tree->SetBranchAddress("nameBranch",&ananame);
  tree->SetBranchAddress("mcevgaBranch",&bmcevga);
  tree->SetBranchAddress("activityBranch",&activity);
  tree->SetBranchAddress("timeBranch",&time);
  Int_t nana=tree->GetEntries();
  
  printf("\n\n");
  printf("%-25s  %13s  %8s  %10s  %10s  %8s  %8s  \n","  NAME","GENERATED","KEPT ","EFF.(%)","ACT. ","TIME ","EXP. ");
  printf("%-25s  %13s  %8s  %10s  %10s  %8s  %8s  \n","=========================","===========","========","==========","==========","========","========");
    
  Float_t mytime;
  Float_t eff;
  Float_t kept;
  Float_t expect;
  Float_t mcevga;
  string temp;
  for(Int_t i=0;i<nana;i++){
    mytime=kept=eff=expect=mcevga=0;
    tree->GetEntry(i);
    temp=(string(ananame)+name).c_str();
    d=(TDirectory*)file->GetDirectory(ananame);
    hist=(TH1D*)d->Get(temp.c_str());
    kept=hist->GetEntries();
    mcevga=bmcevga;
    //if(kept && mcevga) eff=mcevga/kept;
    if(kept && mcevga) eff=kept/mcevga*100;
    mytime=time/60/60/24;
    //if(eff) expect=(activity*time)/eff;
    expect=eff*activity*time;
    
    printf("%-25s  %13.0f  %8.0f  %10.1e  %10.1e  %8.1f  %8.1f  \n",ananame,mcevga,kept,eff,activity,mytime,expect);
    
  }
  printf("\n\n");
  delete file;
}

