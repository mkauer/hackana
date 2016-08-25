// Implementation of ana10 classes and subclasses
#include "ana.hpp"

#include <string>
#include <sstream>
#include <iomanip>
#include "../../utils/h10.h"
#include "TH2.h"
#include "TH1.h"
#include "TLegend.h"
#include "TFractionFitter.h"
#include "TMath.h"
#include "TLimit.h"
#include "TLimitDataSource.h"
#include "TConfidenceLevel.h"

  ana10::ana10(){
    //Init additional tree to store selected events parameters
    events= new TTree((name+"events").c_str(),("Selected events for "+name).c_str());
    events->Branch("Esum",&Esum_ana10,"Esum/D");
    events->Branch("Ediff",&Ediff_ana10,"Ediff/D");
    events->Branch("Cos",&Cos_ana10,"Cos/D");
    events->Branch("w",&w_ana10,"w/D");

    //initialize runlist period and run quality to use
    P1=true;
    P2a=false;
    P2b=false;

    set0=true;
    set1=false;
    set2=false;
    set3=false;
    set4=false;
    extended=false;

}
  ana10::ana10(const char *prefix){
    string pref=prefix;
    name=pref;

  //check if to use individual swactivity for this background component
    if(stest(pref,"sw")){  
      if(stest(pref,"bi210")) swactivity_use = true;
      if(stest(pref,"bi214")) swactivity_use = true;
      if(stest(pref,"pb214")) swactivity_use = true;
    }
    if(stest(pref,"ex") || stest(pref,"pm") || stest(pref,"sc") || stest(pref,"ss")){
      exactivity_use = true;
    }

    exactivity_ini = false;
    noextruemc_warned = false;
    swactivity_ini = false;


    //Init additional tree to store selected events parameters
    events= new TTree((name+"events").c_str(),("Selected events for "+name).c_str());
    events->Branch("Esum",&Esum_ana10,"Esum/D");
    events->Branch("Ediff",&Ediff_ana10,"Ediff/D");
    events->Branch("Cos",&Cos_ana10,"Cos/D");
    events->Branch("w",&w_ana10,"w/D");
    blankh1 = new TH1D((name+"blank").c_str(),"dummy",10,0.,1.);
    blankh2 = new TH2D((name+"blank2").c_str(),"dummy",10,0.,1.,10,0.,1.);
}

Bool_t ana10::Addh1(const char *hnamec, const char * title,Int_t nch, Float_t bmin, Float_t bmax){
  string hname=hnamec;
  //check that there is no such a histogram already
  if(d1map.count(hname)>0) return false;

  //create a hystogram in ROOT memory
  d1map[hname]= new TH1D((name+hname).c_str(),title,nch,bmin,bmax);
  //d1map[hname]->Sumw2();
  return true;
}
Bool_t ana10::Addh1(const char *hnamec, TH1D * source){
  string hname=hnamec;
  //check that there is no such a histogram already
  if(d1map.count(hname)>0) return false;

  //create a hystogram in ROOT memory
  d1map[hname]= (TH1D *)source->Clone((name+hname).c_str());
  d1map[hname]->Sumw2();
  return true;
}
Bool_t ana10::Addh2(const char *hnamec, const char * title,Int_t nch, Float_t bmin, Float_t bmax,Int_t nchy, Float_t ymin, Float_t ymax){
  string hname=hnamec;
  //check that there is no such a histogram already
  if(d2map.count(hname)>0) return false;

  //create a hystogram in ROOT memory
  d2map[hname]= new TH2D((name+hname).c_str(),title,nch,bmin,bmax,nchy,ymin,ymax);
  d2map[hname]->Sumw2();
  return true;
}
TH1D * ana10::Geth1(const char *hnamec){
  string hname=hnamec;
  if(d1map.count(hname)>0){
    return d1map[hname];
  }else{
    cout<<"Geth1 warning: requested histogram "<<hname<<" has not been booked!"<<endl;
    return blankh1;
  }
}
TH2D * ana10::Geth2(const char *hnamec){
  string hname=hnamec;
  if(d2map.count(hname)>0){
    return d2map[hname];
  }else{
    cout<<"Geth2 warning: requested histogram "<<hname<<" has not been booked!"<<endl;
    return blankh2;
  }
}
Bool_t ana10::Checkh1(const char * hnamec){
  string hname=hnamec;
  return (d1map[hname]!=NULL);  
}
Float_t ana10::GetEventWeight(){
  Float_t w;
  w=1;
  if(run<0 && date<100000){
    //all MC should be weighted as (runlist time)/NEVGA
    // date = runlist time in hours, convert to seconds first
    // time = NEVGA
    w*=getActivityWeight();
    Float_t nevganorm;
    if(time>0){
      nevganorm=time;
    }else{
      nevganorm = -1000.*time;
    }
    w*=((Float_t)date)*3600/nevganorm;
    return w;
  }
  if(run<0 && date>100000){
    //all MC should be weighted as (runlist time)/NEVGA
    // date = runlist time in hours, convert to seconds first
    // time = NEVGA/1000
    w=1.0;
    return w;
  }
  return 1.0;
}

void ana10::AddRunTime(Int_t runnumber,Float_t runtime){
  // Adds time of the  run  runnumber into array runtimes
  // All runs from analysis runlist should be here, as it used for MC normalizatiopn later
  
  if(runtimes.count(runnumber)>0) return;
  runtimes[runnumber]=runtime;
  return;
}

Float_t ana10::MCTime(){
  Float_t t=0;
  Int_t runmin=1000000;
  Int_t runmax=0;
  Int_t iitime;
  typedef Timesmap::iterator titerator;
  for(titerator i=runtimes.begin();i!=runtimes.end();i++){
    if(i->first > runmax) runmax= i->first;
    if(i->first < runmin) runmin = i->first;
  }

  for(Int_t i=runmin;i<runmax;i++){
      if(RunList(i,iitime)){
	t+=iitime;
      }
  }

  return t;
}

Float_t ana10::ExpTime(){
  Float_t t=0;
  Int_t runmin=1000000;
  Int_t runmax=0;
  Int_t iitime;
  typedef Timesmap::iterator titerator;
  for(titerator i=runtimes.begin();i!=runtimes.end();i++){
    t+=i->second;
  }

  return t;
}


Float_t ana10:: CalculateTime(){
  Int_t i;
  Float_t t=0;
  Int_t iitime;
  for(i=1000;i<10000;i++){
    if(RunList(i,iitime)){
      t+=iitime;
    }
  }
  return t;
}

Float_t ana10:: getRnActivityRunList(){

  //this function calculates mean Rn activity for current runlist
  Int_t i;
  Float_t rn=0;
  Float_t t=0;
  Int_t iitime;
  for(i=1000;i<RUNMAX;i++){
    if(RunList(i,iitime)){
      rn+=rrn_act[i]*iitime;
      t+=iitime;
    }
  }
  rn = rn / t;
  return rn;
}

ana10::~ana10(){
  //Destruction is tricky since map has only pointers,
  //and when map deleted only pointers will be deleted, not histos themselves!
  typedef TH1Dmap::iterator d1mapiterator;
  for(d1mapiterator i=d1map.begin();i!=d1map.end();i++){ // loop through the map
    delete i->second;//delete hystogram from heap
  }
  typedef TH2Dmap::iterator d2mapiterator;
  for(d2mapiterator i=d2map.begin();i!=d2map.end();i++){ // loop through the map
    delete i->second;//delete hystogram from heap
  }
  delete blankh1;
  delete blankh2;
  delete events;
  delete bEsum;
  delete bEdiff;
  delete bCos;
  delete bw;

}
void ana10::Write(TFile *f){
  f->cd();
  //create sub dir for all elemets
  f->mkdir(name.c_str());
  f->cd(name.c_str());
  //loop through all histograms and svae them into the file
  typedef TH1Dmap::iterator d1mapiterator;
  for(d1mapiterator i=d1map.begin();i!=d1map.end();i++){ // loop through the map
    if(i->second!=NULL){
      i->second->Write();
    }else{
      cout<<"ana10 Write() error, NULL histo in the 1D list!? "<<endl;
    }
  }
  typedef TH2Dmap::iterator d2mapiterator;
  for(d2mapiterator i=d2map.begin();i!=d2map.end();i++){ // loop through the map
    if(i->second!=NULL){
      i->second->Write();
    }else{
      cout<<"ana10 Write() error, NULL histo in the 2D list!? "<<endl;
    }
  }
  //save selected events tree
  events->Write();


}
void ana10::Read(TFile *f){
  const char *cname=name.c_str();
  TDirectory *dir;
  f->GetObject(cname,dir);
  TList *d1list=dir->GetListOfKeys();
  TIter next(d1list);
  TKey * temp;
    while(temp = (TKey *)next()){
      string clname;
      clname=temp->GetClassName();
      if(clname=="TH1D"){
	string hname;
	hname=temp->GetName();
	Int_t gn=name.size();
	hname.erase(0,gn);
	d1map[hname]=(TH1D *)dir->Get(temp->GetName());
      }
      if(clname=="TH2D"){
	string hname;
	hname=temp->GetName();
	Int_t gn=name.size();
	hname.erase(0,gn);
	d2map[hname]=(TH2D *)dir->Get(temp->GetName());
      }
      if(clname=="TTree"){
	if(events) delete events;
	events = (TTree *)dir->Get(temp->GetName());
	bEsum = events->GetBranch("Esum");
	bEdiff = events->GetBranch("Ediff");
	bCos = events->GetBranch("Cos");
	bw = events->GetBranch("w");
	bEsum->SetAddress(&Esum_ana10);
	bEdiff->SetAddress(&Ediff_ana10);
	bCos->SetAddress(&Cos_ana10);
	bw->SetAddress(&w_ana10);
      }
    }
}
void ana10::Read(TFile *f, string prefix ){
  // use this read, if prefix in the file not coincide with this->name
  const char *cname=prefix.c_str();
  TDirectory *dir;
  f->GetObject(cname,dir);
  TList *d1list=dir->GetListOfKeys();
  TIter next(d1list);
  TKey * temp;
    while(temp = (TKey *)next()){
      string clname;
      clname=temp->GetClassName();
      if(clname=="TH1D"){
	string hname;
	hname=temp->GetName();
	Int_t gn=prefix.size();
	hname.erase(0,gn);
	d1map[hname]=(TH1D *)dir->Get(temp->GetName())->Clone((name+hname).c_str());
      }
      if(clname=="TH2D"){
	string hname;
	hname=temp->GetName();
	Int_t gn=prefix.size();
	hname.erase(0,gn);
	d2map[hname]=(TH2D *)dir->Get(temp->GetName())->Clone((name+hname).c_str());
      }
      if(clname=="TTree"){
	if(events) delete events;
	string hname;
	hname=temp->GetName();
	events = (TTree *)dir->Get(temp->GetName())->Clone((name+hname).c_str());
	bEsum = events->GetBranch("Esum");
	bEdiff = events->GetBranch("Ediff");
	bCos = events->GetBranch("Cos");
	bw = events->GetBranch("w");
	bEsum->SetAddress(&Esum_ana10);
	bEdiff->SetAddress(&Ediff_ana10);
	bCos->SetAddress(&Cos_ana10);
	bw->SetAddress(&w_ana10);
      }
    }
}
void ana10::Getname(char *cname){
  strcpy(cname,name.c_str());
  return;
}
const char * ana10::Getname(){
  return name.c_str();
}

void ana10::Geth1list(TObjArray *a){
  //loop through all histograms and save them into the file
  typedef TH1Dmap::iterator d1mapiterator;
  for(d1mapiterator i=d1map.begin();i!=d1map.end();i++){ // loop through the map
    a->AddLast(i->second);
  }

}
void ana10::Geth2list(TObjArray *a){
  //loop through all histograms and save them into the file
  typedef TH2Dmap::iterator d2mapiterator;
  for(d2mapiterator i=d2map.begin();i!=d2map.end();i++){ // loop through the map
    a->AddLast(i->second);
  }

}
bana10::bana10(const char *pref,Long64_t nmc,Float_t a):ana10(pref){


    mcevga=nmc;
    activity=a;
    eactivity=0;
    norm=a;
    enorm=0;
    parent=NULL;
    dependent=false;
    scaled=false;
    cout<<"Init "<<pref<<" "<<nmc<<" "<<a<<endl;

}
bana10::~bana10(){
  delete Esum_set;
  delete Hist_set;
  delete Esum_pdf;
  delete E2D_pdf;
  delete d1pdf;
}


Int_t bana10::ConstructPdf(RooRealVar& roo,Bool_t key)
{

  if(key){
    RooRealVar w("w","w",0.,100000.);
    RooArgSet tmp(roo,w);
    Esum_set = new RooDataSet((name + "esumset").c_str(),"Esum data set",events,tmp,0,"w");
    Esum_pdf = new RooKeysPdf((name + "_esumpdf").c_str(),"Esum pdf",roo,* (RooDataSet *)Esum_set);
    d1pdf=Esum_pdf;
  }else{
    //prepare hist
    Int_t nentries = events->GetEntries();
    Float_t *points = new Float_t[nentries];
    if(nentries){
      for(Int_t i =0;i<nentries;i++){
	events->GetEntry(i);
	points[i]=Esum_ana10;
      }
      //sort all the points
      qsort(points,nentries,sizeof(Float_t),compare_doubles);
      Float_t *xbins = new Float_t[1000];
      Int_t i,j;
      xbins[0]=roo.getMin();
      j=1;
      xbins[j]=xbins[0];
      Int_t ncount=0;
      Float_t step = 0.012;
      Int_t full=2;
      for(i=0;i<nentries;i++){
	if(points[i]>xbins[0]){
	  if(points[i]>xbins[j] && ncount<full){
	    while (points[i]>xbins[j]) xbins[j]+=step;
	  }
	  if(points[i]>xbins[j] && ncount>=full){
	    j++;
	    xbins[j]=xbins[j-1]+step;
	    ncount=0;
	  }
	  if(points[i]<xbins[j]) ncount++;
	}
      }
      j++;
      xbins[j]=roo.getMax();
      TH1F *dist = new TH1F((name + "esumhist").c_str(),"Esum data set",j,xbins);// init with variable bin size
      for(i =0;i<nentries;i++){
	events->GetEntry(i);
	dist->Fill(Esum_ana10,w_ana10);
      }
      Esum_set = new RooDataHist((name + "esumset").c_str(),"Esum data set",RooArgSet(roo),dist);
      d1pdf = new RooHistPdf((name + "_esumpdf").c_str(),"Esum pdf",roo,* (RooDataHist *)Esum_set);
    }else{
      Esum_set = new RooDataHist((name + "esumset").c_str(),"Esum data set",RooArgSet(roo));
      d1pdf = new RooHistPdf((name + "_esumpdf").c_str(),"Esum pdf",roo,* (RooDataHist *)Esum_set);
    }
  }
  
  return 1;
}

Int_t bana10::ConstructSet(RooRealVar& roo,Bool_t key)
{

  if(key){
    roo.Print();
    RooRealVar w("w","w",0.,100000.);
    RooArgSet tmp(roo,w);
    Esum_set = new RooDataSet((name + "esumset").c_str(),"Esum data set",events,tmp,0,"w");
    Esum_set->Print();
    d1pdf=Esum_pdf;
  }else{
    //prepare hist
    Int_t nentries = events->GetEntries();
    Float_t *points = new Float_t[nentries];
    if(nentries){
      for(Int_t i =0;i<nentries;i++){
	events->GetEntry(i);
	points[i]=Esum_ana10;
      }
      //sort all the points
      qsort(points,nentries,sizeof(Float_t),compare_doubles);
      Float_t *xbins = new Float_t[1000];
      Int_t i,j;
      xbins[0]=roo.getMin();
      j=1;
      xbins[j]=xbins[0];
      Int_t ncount=0;
      Int_t full=2;
      for(i=0;i<nentries;i++){
	if(points[i]>xbins[0]){
	  if(points[i]>xbins[j] && ncount<full){
	    while (points[i]>xbins[j]) xbins[j]+=0.01;
	  }
	  if(points[i]>xbins[j] && ncount>=full){
	    j++;
	    xbins[j]=xbins[j-1]+0.01;
	    ncount=0;
	  }
	  if(points[i]<xbins[j]) ncount++;
	}
      }
      j++;
      xbins[j]=roo.getMax();
      TH1F *dist = new TH1F((name + "esumhist").c_str(),"Esum data set",j,xbins);// init with variable bin size
      for(i =0;i<nentries;i++){
	events->GetEntry(i);
	if(Esum_ana10>roo.getMin())dist->Fill(Esum_ana10,w_ana10);
      }
      Esum_set = new RooDataHist((name + "esumset").c_str(),"Esum data set",RooArgSet(roo),dist);
    }else{
      Esum_set = new RooDataHist((name + "esumset").c_str(),"Esum data set",RooArgSet(roo));
    }
  }
  
  return 1;
}

Int_t bana10::Construct2DEPdf(RooRealVar& roosum,RooRealVar& roodiff,Bool_t key){

  RooRealVar w("w","w",0.,100000.);
  RooArgSet tmp(roosum,roodiff,w);
  Esum_set = new RooDataSet((name + "esumset").c_str(),"Esum data set",events,tmp,0,"w");

  if(key){
    E2D_pdf = new Roo2DKeysPdf((name + "_esumpdf").c_str(),"Esum pdf",roosum,roodiff,*(RooDataSet *)Esum_set);
  }else{
    RooArgSet tmp2(roosum,roodiff);
    Hist_set = new RooDataHist((name + "ehistset").c_str(),"Esum data hist",tmp2,*Esum_set);
    E2D_pdf = new  RooHistPdf((name + "_esumpdf").c_str(),"Esum pdf",tmp2,*Hist_set);  
  }
  return 1;
}



void bana10::setnorm(Float_t t){
  //set norm for objects created from slim ROOT
  //if t!=0 => old normalisation scheme (time/NEVGA)
  //need to rescale histograms for using activity as a direct normalisation

  //SNSW root files 
  //no scaling
  if(run==-1){
    cout<<"Skip MC scaling"<<endl;
    return;

  }

  if(t==0){
    //    totaltime=ExpTime();
    totaltime=CalculateTime();
  }else{
    totaltime=t;
  }
  Float_t scalef;

  if(date>100000){
    scalef=totaltime/mcevga;
  }else{
    if(MCTime()){
      scalef=totaltime/MCTime();
    }else{
      scalef=0;
    }
  }

  if(mcevga && !scaled){
    scaled=true;
    TObjArray dummy,dummy2;
    Geth1list(&dummy);
    for(Int_t i=0;i<dummy.GetEntries();i++){
      ((TH1D*)dummy[i])->Scale(scalef);
    }
        Geth2list(&dummy2);
    for(Int_t i=0;i<dummy2.GetEntries();i++){
      ((TH2D*)dummy2[i])->Scale(scalef);
    }
  }
  if(mcevga){
      norm=activity;
      enorm=eactivity;
  }else{
    norm=1;
    enorm=0;
  }

  //  cout << "setnorm "<<name<<"  t="<<t<<" totaltime="<<totaltime<<endl;
  //cout<<"     mcevga = "<<mcevga<<" mctime="<<MCTime()<<endl;
  //cout<<" activity  = "<<activity<<" eactivity = "<<eactivity<<endl;
}

void bana10::scale(Float_t s){
  cout<<"Scale "<<name<<" factor "<<s<<endl;
  TObjArray dummy,dummy2;
  Geth1list(&dummy);
  for(Int_t i=0;i<dummy.GetEntries();i++){
    ((TH1D*)dummy[i])->Scale(s);
  }
  Geth2list(&dummy2);
  for(Int_t i=0;i<dummy2.GetEntries();i++){
    ((TH2D*)dummy2[i])->Scale(s);
  }
}


void bana10::setnormR(Float_t t){
  // set norm for files reread from saved ROOT histograms
  if(t==0){
    totaltime=CalculateTime();
  }else{
    totaltime=t;
  }    
  if(mcevga){
      norm=activity;
      enorm=eactivity;
    }
}



void bana10::setnormalization(Float_t v,Float_t ev){
  norm=v;
  enorm=ev;
  activity=norm;
  eactivity=enorm;
}

Double_t bana10::GetExpected(const char *name, Int_t lx,Int_t ux){
  //return expected number of events for 1D histogram
  Double_t ni;
  if(lx==0 && ux==0){
    ni = Geth1(name)->Integral();
  }else{
    ni = Geth1(name)->Integral(lx,ux);
  }
  ni = ni*norm;
  return ni;
}
Double_t bana10::GetErrorExpected(const char *name, Int_t lx, Int_t ux){
  Double_t ni;
  Double_t dni;
  TH1D * hd1 = Geth1(name);
  if(lx==0 && ux==0){
    ni = Geth1(name)->Integral();
  }else{
    ni = Geth1(name)->Integral(lx,ux);
  }

  if(lx==0 && ux==0){
    lx = 0;
    ux = hd1->GetNbinsX()+1;
  }
  dni=0;
  for(Int_t i=lx;i<=ux;i++){
    dni+=pow(hd1->GetBinError(i),2);
  }
  dni = dni*norm*norm + pow(ni*enorm,2);
  //  cout<<" GetErrorExpected "<<this->Getname()<<" enorm="<<enorm<<" ni="<<ni<<" dni="<<dni<<endl;
  dni = sqrt(dni);
  return dni;
}

int bana10::Print(){
  typedef TH1Dmap::iterator d1mapiterator;
  d1mapiterator i=d1map.begin(); // take any histogram
  if(name=="data" || name=="DATA" || name=="Phase1" || name=="Phase2" ){
    anaout<<"Experimental time = "<<totaltime/3600/24<<" days, "<<i->second->Integral()<<" events"<<endl;
  }else{       // MC background
    if(enorm){ // results of Fit
      anaout<<name<<"   A = "<<activity<<"("<< eactivity<<")   eff = "<<i->second->Integral()/totaltime<<"   Nexp = "<<i->second->Integral()*norm<<endl;
    }else{     // known background
      anaout<<name<<"   A = "<<activity<<"   eff = "<<i->second->Integral()/totaltime<<"   Nexp = "<<i->second->Integral()*norm<<endl;
    }
  }
  return 0;
}

int bana10::Print(const char *hname){
  if(name=="data" || name=="DATA" || name=="Phase1" || name=="Phase2" ){
    anaout<<"Experimental time = "<<totaltime/3600/24<<" days, "<<Geth1(hname)->Integral()<<" events"<<endl;
  }else{        // MC background
    if(enorm){  // results of Fit
      anaout<<name<<"   A =  "<<activity<<" +/- "<< eactivity<<"   eff = "<<GetExpected(hname)/totaltime/activity<<"   Nexp = "<<GetExpected(hname)<<" +/- "<<GetErrorExpected(hname)<<endl;
    }else{      // known background
      anaout<<name<<"   A =  "<<activity<<"   eff = "<<GetExpected(hname)/totaltime/activity<<"   Nexp = "<<GetExpected(hname)<<" +/- "<<GetErrorExpected(hname)<<endl;
    }
  }
  return 0;
}

void bana10::Addrelative(bana10 * child, Float_t factor){
  relative[nrelative]=child;
  relfactor[nrelative]=factor;
  child->activity=activity*factor;
  child->dependent=true;
  child->parent=this;
  nrelative++;
}



/// Saving-Reading data from the files routines


Int_t ReadSavedRoot(const char * fname,const char * name,Float_t activity,bana10 *&ana){
//read just 1 component
   TFile *oldfile=new TFile(fname);
   char ananame[256];
   Long64_t bmcevga;
   Float_t time;
   TTree *maintree=(TTree *)oldfile->Get("ana10");
   maintree->SetBranchAddress("nameBranch",&ananame);
   maintree->SetBranchAddress("mcevgaBranch",&bmcevga);
   maintree->SetBranchAddress("timeBranch",&time);
   Int_t nana=maintree->GetEntries();
   for(Int_t i=0;i<nana;i++){
     maintree->GetEntry(i);
     if(string(ananame)==name){
       if(string(ananame)=="data"){
	 activity=time;
	 if(ana!=NULL) {
	   cout<<"Error! Two experimantal data ntuples are not allowed"<<endl;
	   return 1;
	 }
       }
       ana = new bana10(ananame,bmcevga,activity);
       ana->scaled=true;
       ana->setnorm(time);
       ana->Read(oldfile);
     }
   }
 return 0;
}

Int_t ReadSavedRootCh(const char * fname,const char * name,string nchnl,Float_t activity,bana10 *&ana, Float_t forcemcevga, Float_t scalebkg){
  //read just 1 component
  //put it in channel nchnl
  string sfname = string("root://")+string(fname);
  TFile *oldfile= new TFile(fname);
  Int_t chs=nchnl.size();
  string rname=string(name).erase(0,chs);
  char ananame[256];
  Long64_t bmcevga;
  Float_t time;
  Float_t saved_activity;
  TTree *maintree=(TTree *)oldfile->Get("ana10");
  maintree->SetBranchAddress("nameBranch",&ananame);
  maintree->SetBranchAddress("mcevgaBranch",&bmcevga);
  maintree->SetBranchAddress("activityBranch",&saved_activity);
  maintree->SetBranchAddress("timeBranch",&time);
  Int_t nana=maintree->GetEntries();
  for(Int_t i=0;i<nana;i++){
    maintree->GetEntry(i);
    if(string(ananame)==rname){
      if(string(ananame)=="data"){
	activity=time;
	if(ana!=NULL) {
	  cout<<"Error! Two experimantal data ntuples are not allowed"<<endl;
	  return 1;
	}
      }
      Float_t scalef=1;
      if(forcemcevga > 0){
	scalef = ((Float_t)bmcevga)/forcemcevga;
	bmcevga=(Long64_t)forcemcevga;
      }
      
      if(activity < 0 ){
	activity = saved_activity;
	cout<<"For the component "<<name<<", contral.dat has negative A. Use A saved in ntuple:  "<<scientific<<setprecision(2)<<activity<<endl;
      }
      
      if(scalebkg >= 0){
	activity = activity * scalebkg;
	cout<<"For the component "<<name<<", contral.dat has 'S' option to scale by:  "<<fixed<<scalebkg<<endl;
	cout<<"For the component "<<name<<", the new activity is now set to:  "<<scientific<<setprecision(2)<<activity<<endl;
      }
      
      ana = new bana10(name,bmcevga,activity);
      ana->scaled=true;
      ana->setnorm(time);
      ana->Read(oldfile,rname);
      ana->scale(scalef);
    }
  }
  //   oldfile->Close();
  //   delete oldfile;
  return 0;
}


Int_t ReadSavedRoot(const char * fname,bana10 *&data,bana10 **bgr,Int_t &nbgr,Float_t &exptime){

   TFile *oldfile=new TFile(fname);
   char ananame[256];
   Long64_t bmcevga;
   Float_t activity;
   Float_t time;
   TTree *maintree=(TTree *)oldfile->Get("ana10");
   maintree->SetBranchAddress("nameBranch",&ananame);
   maintree->SetBranchAddress("mcevgaBranch",&bmcevga);
   maintree->SetBranchAddress("activityBranch",&activity);
   maintree->SetBranchAddress("timeBranch",&time);
   Int_t nana=maintree->GetEntries();
   for(Int_t i=0;i<nana;i++){
     maintree->GetEntry(i);
     if(string(ananame)=="data"){
       if(data!=NULL) {
	 cout<<"Double record for data, skip it"<<endl;
       }else{
	 data = new bana10(ananame,bmcevga,activity);
	 data->scaled=true;
	 data->setnorm(time);
	 data->Read(oldfile);
	 exptime=time;
       }
     }else{
       Bool_t exist=false;
       for(Int_t j=0;j<nbgr;j++){
	 if(string(ananame)==bgr[j]->Getname()){
	   cout<<"Double record for "<<bgr[j]->Getname()<<" skip it"<<endl;
	   exist=true;
	 }
       }
       if(!exist){
	 bgr[nbgr] = new bana10(ananame,bmcevga,activity);
	 bgr[nbgr]->scaled=true;
	 bgr[nbgr]->setnorm(time);
	 bgr[nbgr]->Read(oldfile);
	 nbgr++;
       }
     }
   }
   return 0;
}
Int_t ReadSavedRoot(const char * fname,bana10 *&data,bana10 **bgr,Int_t &nbgr,bana10 ** sig, Int_t nsig,Float_t &exptime){
  TFile *oldfile=new TFile(fname);
  char ananame[256];
  Long64_t bmcevga;
  Float_t activity;
  Float_t time;
  TTree *maintree=(TTree *)oldfile->Get("ana10");
  maintree->SetBranchAddress("nameBranch",&ananame);
  maintree->SetBranchAddress("mcevgaBranch",&bmcevga);
  maintree->SetBranchAddress("activityBranch",&activity);
  maintree->SetBranchAddress("timeBranch",&time);
  Int_t nana=maintree->GetEntries();
  for(Int_t i=0;i<nana;i++){
     maintree->GetEntry(i);
     if(string(ananame)=="data"){
       if(data!=NULL) {
	 cout<<"Double record for data, skip it"<<endl;
       }else{
	 data = new bana10(ananame,bmcevga,activity);
	 data->scaled=true;
	 data->setnorm(time);
	 data->Read(oldfile);
	 exptime=time;
       }
     }else{
       Bool_t exist=false;
       for(Int_t j=0;j<nbgr;j++){
	 if(string(ananame)==bgr[j]->Getname()){
	   cout<<"Double record for "<<bgr[j]->Getname()<<" skip it"<<endl;
	   exist=true;
	 }
       for(Int_t j=0;j<nsig;j++){
	 if(string(ananame)==sig[j]->Getname()){
	   cout<<"Double record for "<<sig[j]->Getname()<<" skip it"<<endl;
	   exist=true;
	 }
       }
      }
       if(!exist){
	 bgr[nbgr] = new bana10(ananame,bmcevga,activity);
	 bgr[nbgr]->scaled=true;
	 bgr[nbgr]->setnorm(time);
	 bgr[nbgr]->Read(oldfile);
	 nbgr++;
       }
     }
   }
 return 0;
}

Int_t ReadSlim(const char * fname,const char *name,Float_t activity,Long64_t mcevga,bana10 *&ana,Int_t version){
 	TChain *chain = new TChain("h10");
	Int_t runnumber;
	Int_t timet;
	chain->SetBranchAddress("run",&runnumber);
	chain->SetBranchAddress("time",&timet);
	chain->Add(fname);	      
	Int_t nentries = Int_t(chain->GetEntries());
	//cout<<"!!!!!!!!!!  In "<<fname<<" in "<<name<<" found "<<nentries<<" entries. !!!!!!!!!!"<<endl;
	chain->GetEntry(nentries-1);
	if(timet!=0) {
	  Long64_t rootmcevga;
	  if(timet<0){
	    rootmcevga=-(Long64_t)1000*timet;
	  }else{
	    rootmcevga=timet;
	  }
	  cout<<"NEVGA="<<rootmcevga<<" found in the root file "<<endl;
	  cout<<"Provided value "<<mcevga<<endl;
	  if(mcevga==0){
	    mcevga=rootmcevga;
	    cout<<"Use mcevga="<<mcevga<<endl;
	  }
	}
	ana = new bana10(name,mcevga,activity);
	if(ana->swactivity_use) cout<<"Will try to use individual layer/sector activity"<<endl;
	if(ana->exactivity_use) cout<<"Will try to use EXBG MODEL F activity from exbg_modelf.dat"<<endl;
	ana->nemosv=version;
	ana->Init(chain);
	ana->Loop(0);
	return 0;
}
Int_t SaveRoot(const char *fname,bana10 *data,bana10** bgr,Int_t nbgr,Bool_t recreate){
  TFile * newfile;
  if(recreate){
    newfile = new TFile(fname,"RECREATE");
  }else{
    newfile = new TFile(fname,"UPDATE");    
  }
  TTree * maintree;
  char ananame[256];
  Long64_t bmcevga;
  Float_t activity;
  Float_t time;

  if(recreate){
    maintree=new TTree("ana10","Analysisi results tree");
    maintree->Branch("nameBranch",&ananame,"name/C");
    maintree->Branch("mcevgaBranch",&bmcevga,"mcevga/L");
    maintree->Branch("activityBranch",&activity,"activity/F");
    maintree->Branch("timeBranch",&time,"time/F");
  }else{
    maintree= (TTree * )newfile->Get("ana10");
    maintree->SetBranchAddress("nameBranch",&ananame);
    maintree->SetBranchAddress("mcevgaBranch",&bmcevga);
    maintree->SetBranchAddress("activityBranch",&activity);
    maintree->SetBranchAddress("timeBranch",&time);
  }
      
      for(Int_t i=0;i<nbgr;i++){
	bgr[i]->Getname(ananame);
	bmcevga=bgr[i]->mcevga;
	activity=bgr[i]->activity;
	time=bgr[i]->totaltime;
	bgr[i]->Write(newfile);
	maintree->Fill();
      }
      if(data){
	data->Getname(ananame);
	bmcevga=data->mcevga;
	activity=data->activity;
	time=data->totaltime;
	data->Write(newfile);
	maintree->Fill();
      }
      newfile->cd();
      maintree->Write();
      newfile->Write();
      newfile->Close();
}

TH1 * GetCumulativeDown(TH1 * src,Int_t nx)
{
  string name = src->GetName();
  name = "cum_"+name;
  TH1 * dst = (TH1 *)src->Clone(name.c_str());
  dst->Reset();
  if(nx > src->GetNbinsX()) nx = src->GetNbinsX();
  Double_t bin=0;
  Double_t error=0;
  for(Int_t i=nx;i>0;i--){
    bin+=src->GetBinContent(i);
    error+=pow(src->GetBinError(i),2);
    dst->SetBinContent(i,bin);
    dst->SetBinError(i,sqrt(error));
  }
  return dst;
}
TH1 * GetCumulativeDown(TH1 * src, Double_t bup){
  for(Int_t i = 0;i<src->GetNbinsX();i++){
    if(bup<src->GetBinLowEdge(i)){
      return GetCumulativeDown(src,i-1);
    }
  }
  return GetCumulativeDown(src,src->GetNbinsX());
}
Bool_t stest(string s1, string s2){
  return (s1.find(s2) != string::npos);
}

