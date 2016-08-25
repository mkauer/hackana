#define h10_cxx
#define h10_cut

#include "ana.hpp"
#include "TH2.h"
#include "TStyle.h"
#include "TNtuple.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TApplication.h"
#include "TRint.h"
#include "TGApplication.h"
#include "TLine.h"
#include "TPaveLabel.h"

#include <fstream>
#include <sstream>
#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <string>

#include "TROOT.h"

//patch against old gcc 3.2.3 error
char* operator+( std::streampos&, char* );

//extra functions
void Usage();

//TROOT root("app","N3 analysis program");


using namespace std;
ofstream anaout;
ofstream erlog;
//Pmstatus *pmstatus;

Int_t rtime[RUNMAX];
Float_t rrn_act[RUNMAX];
Int_t rstatus[RUNMAX];
Int_t run_date[RUNMAX];


int main(int argc, char *argv[]){
  //argv[0] is name of programme 
  std::cout<<"\n\t ROOT Analysis (RootAna) Program \n";
  std::cout<<  "\t Author V. Vasiliev 01.2006 \n";
  std::cout<<  "\t ./ana --help for options \n\n";
  
  
  Int_t timetot=0;
  Bool_t mc=false;
  Bool_t filelist=false;
  Bool_t rootinput=true;
  TTree* tree=0;
  
  const Int_t MAXCHNLS=10;
  const Int_t MAXCOMP=400;
  Int_t nbgr[MAXCHNLS];
  Float_t exptime[MAXCHNLS];
  Int_t nchnls=1;// number of channels
  string namechnl[MAXCHNLS];
  bana10 *bgr[MAXCHNLS][MAXCOMP];
  bana10 *data[MAXCHNLS];
  bana10 *sig[MAXCHNLS][MAXCOMP];
  Int_t nsig[MAXCHNLS];
  bana10 *lim[MAXCHNLS][MAXCOMP];
  Int_t nlim[MAXCHNLS];

  //init Pmstatus cache
  // pmstatus = new Pmstatus();

  //init part
  for(Int_t i=0;i<MAXCHNLS;i++){
    nbgr[i]=0;
    nsig[i]=0;
    exptime[i]=0;
    data[i]=NULL;
    nlim[i]=0;
  }
  
  
  //use UCL mirror server
  //nemo3::DbMgr::Init("nemodb.hep.ucl.ac.uk",3306);   // old ucl server
  //nemo3::DbMgr::Init("nemodb1.hep.ucl.ac.uk",3306);  // new ucl server
  //nemo3::DbMgr::Init("ccnemodb.in2p3.fr",13306);     // in2p3 server
  // nemo3::DbMgr::Init("db1.hep.ucl.ac.uk",3306);         // newer ucl server
  
  
  //LOOP root ntuples, read parameters from dat file provided
  char *control;
  char *rooty;
  Int_t numFiles;
  ifstream in0, in1, in2;
  int batch=0;
  for(int i=1;i<=argc;i++){
    if(string(argv[i])=="-h"||string(argv[i])=="help"||string(argv[i])=="-help"||string(argv[i])=="--help"){
      Usage();
      return 1;
    }
    if(string(argv[i])=="-b"){
      batch=1;
      i++;
    }
    if(argv[i] && ! argv[i+1]){
      control=argv[i];
      string stemp=string(argv[i]);
      int len=strlen(stemp.c_str())-1;
      while(stemp[len] != '.') len--;
      stemp.erase(len,stemp.length());
      stemp=stemp+".root";
      rooty=new char[stemp.length()+1];
      strcpy(rooty,stemp.c_str());
      i++;
    }else if(argv[i] && argv[i+1] && !(string(argv[i+1])=="-b")){
      control=argv[i];
      rooty=argv[i+1];
      i=i+2;
    }else if(argv[i] && argv[i+1] && string(argv[i+1])=="-b"){
      batch=1;
      control=argv[i];
      string stemp=string(argv[i]);
      int len=strlen(stemp.c_str())-1;
      while(stemp[len] != '.') len--;
      stemp.erase(len,stemp.length());
      stemp=stemp+".root";
      rooty=new char[stemp.length()+1];
      strcpy(rooty,stemp.c_str());
      i=i+2;
    }else{
      control="control.dat";
      rooty="default.root";
    }
  }
  
  std::cout<<""<<std::endl;
  std::cout<<"Using control file  ==> "<<control<<std::endl;
  std::cout<<"Using rootfile file ==> "<<rooty<<std::endl;
  std::cout<<"Using batch-mode    ==> "<<batch<<std::endl;
  std::cout<<""<<std::endl;
  //return 0;
  
  in0.open(control);
  if(in0.is_open()){
  char filenames[MAXCOMP][256];
  string bname[MAXCOMP];
  Long64_t mcevga[MAXCOMP];
  Float_t bactivity[MAXCOMP];
  Float_t eactivity[MAXCOMP];
  Int_t runnumber[MAXCOMP];
  Int_t datenumber[MAXCOMP];
  string tag[MAXCOMP];
  string relative[MAXCOMP];
  Float_t rfactor[MAXCOMP];
  string nchannel[MAXCOMP];
  Int_t channel[MAXCOMP];
  Float_t scalebkg[MAXCOMP];
  
  erlog.open("errorlog.dat");
  anaout.open("ana10.dat");

  Int_t i=0;
  Bool_t endf=false;
  nchnls=1;
  namechnl[0]="";
  while(!endf){
    filenames[i][0]=0;
    bname[i]="";
    bactivity[i]=0;
    mcevga[i]=0;
    tag[i]="";
    relative[i]="";
    rfactor[i]=0;
    eactivity[i]=0;
    channel[i]=0;
    nchannel[i]="";
    scalebkg[i]=-1;
    in0 >> filenames[i]>>bname[i]>>bactivity[i]>>mcevga[i]>>tag[i];
    channel[i]=0;
    eactivity[i]=0;
    while(in0.peek()!=10 && in0.peek()!=EOF){
      in0.get();
      if(in0.peek()=='D'){ // Dependancy
	in0.get();
	in0>>relative[i]>>rfactor[i];
      }
      if(in0.peek()=='E'){ // Error
	in0.get();
	in0>>eactivity[i];
      }
      if(in0.peek()=='S'){ // Scale
	in0.get();
	in0>>scalebkg[i];
      }
      if(in0.peek()=='C'){ // Channel
	in0.get();
	in0>>channel[i]>>nchannel[i];
	namechnl[channel[i]]=nchannel[i];
	if(channel[i]+1>nchnls) nchnls=channel[i]+1;
      }      
    }
    string s(filenames[i]);
    if(s.size()==0) endf=true;
    if(filenames[i][0]!='#') {
      cout<<i<<" "<<filenames[i]<<" "<<bname[i]<<" "<<bactivity[i]<<" "
	  <<mcevga[i]<<" "<<tag[i]<<" "<<relative[i]<<" "<<rfactor[i]
	  <<" "<<scalebkg[i]<<" "<<channel[i]<<" "<<nchannel[i]<<endl;
      i++;
    }else{
      cout<<"skipping --> "<<filenames[i]<<endl;
      filenames[i][0]=0;
    }
  } // end of while(!endf)
  in0.close();
  
  numFiles=i-1;
  std::cout<<numFiles<<" files to process"<<std::endl;
  for (int i=0;i<numFiles;i++){  
    Int_t ich=channel[i];
    Int_t nemosv=7;
    std::cout<<"Adding: "<<filenames[i]<<" channel "<<namechnl[ich]<<std::endl;
    bname[i]=namechnl[ich]+bname[i];
    if(tag[i]=="L" || tag[i]=="L7" || tag[i]=="L8"){
      if(tag[i]=="L8") nemosv=8;
      ReadSlim(filenames[i],bname[i].c_str(),bactivity[i],mcevga[i],lim[ich][nlim[ich]],nemosv);
      nlim[ich]++;
    }
    if(tag[i]=="LR"){
      ReadSavedRootCh(filenames[i],bname[i].c_str(),namechnl[ich],bactivity[i],lim[ich][nlim[ich]]);
      nlim[ich]++;
    }
    if(tag[i]=="S" || tag[i]=="S7" || tag[i]=="S8"){
      if(tag[i]=="S8") nemosv=8;
      ReadSlim(filenames[i],bname[i].c_str(),bactivity[i],mcevga[i],sig[ich][nsig[ich]],nemosv);
      nsig[ich]++;
    }
    if(tag[i]=="SR"){
      ReadSavedRootCh(filenames[i],bname[i].c_str(),namechnl[ich],bactivity[i],sig[ich][nsig[ich]],mcevga[i]);
      nsig[ich]++;
    }
    if(tag[i]=="B"|| tag[i]=="B7" || tag[i]=="B8"){
      if(tag[i]=="B8") nemosv=8;
      ReadSlim(filenames[i],bname[i].c_str(),bactivity[i],mcevga[i],bgr[ich][nbgr[ich]],nemosv);
      bgr[ich][nbgr[ich]]->eactivity=eactivity[i];
      nbgr[ich]++;
    }
    if(tag[i]=="BR"){
      //ReadSavedRootCh(filenames[i],bname[i].c_str(),namechnl[ich],bactivity[i],bgr[ich][nbgr[ich]]);
      ReadSavedRootCh(filenames[i],bname[i].c_str(),namechnl[ich],bactivity[i],bgr[ich][nbgr[ich]],0,scalebkg[i]);
      //if non-zero mcevga in dat control file, override old mcevga, saved in root
      if(mcevga[i]!=0) bgr[ich][nbgr[ich]]->mcevga=mcevga[i];
      bgr[ich][nbgr[ich]]->eactivity=eactivity[i];
      nbgr[ich]++;
    }
    if(tag[i]=="D"|| tag[i]=="D7" || tag[i]=="D8"){
      std::cout<<"Find data slim "<<bname[i]<<std::endl;
      if(tag[i]=="D8") nemosv=8;
      ReadSlim(filenames[i],bname[i].c_str(),bactivity[i],mcevga[i],data[ich],nemosv);
      exptime[ich]=bactivity[i];
      data[ich]->mcevga=0;
      data[ich]->setnorm(exptime[ich]);

    }
    if(tag[i]=="DR"){
      std::cout<<"Find data slim "<<bname[i]<<std::endl;
      ReadSavedRootCh(filenames[i],bname[i].c_str(),namechnl[ich],bactivity[i],data[ich]);
      data[ich]->mcevga=0;
      exptime[ich]=data[ich]->totaltime;
    }
    //Forbid to read all contents of ROOT file. Temporarily, because not clear how to be in the case of multiple channels. Demand user explicitly tell which channels and components to load from root file.
    //    if(tag[i]=="R"){
    //      ReadSavedRoot(filenames[i],data[ich],bgr[ich],nbgr,sig,nsig,exptime);
    //    }
  }

  cout<<"Start normalization"<<endl;
  for(Int_t k=0;k<nchnls;k++){
    for(Int_t i=0;i<nbgr[k];i++){
      bgr[k][i]->setnorm(exptime[k]);
    }
    for(Int_t i=0;i<nsig[k];i++){
      sig[k][i]->setnorm(exptime[k]);
    }
    for(Int_t i=0;i<nlim[k];i++){
      lim[k][i]->setnorm(exptime[k]);
    }
  }
  cout<<"End normalization"<<endl;
  //add dependencies between components for fit and drawing purposes
  for(Int_t i=0;i<numFiles;i++){
    bana10 *parent,*child;

    child=NULL;
    parent=NULL;
    Int_t ich = channel[i];
    if(relative[i].size()!=0){
      //look for references
      for(Int_t j=0;j<nbgr[ich];j++){
	if(bgr[ich][j]->Getname()==bname[i]) child=bgr[ich][j];
	if(bgr[ich][j]->Getname()==namechnl[ich]+relative[i]) parent=bgr[ich][j];
      }
      for(Int_t j=0;j<nsig[ich];j++){
	if(sig[ich][j]->Getname()==bname[i]) child=sig[ich][j];
	if(sig[ich][j]->Getname()==namechnl[ich]+relative[i]) parent=sig[ich][j];
      }
      //if references found, add as relative
      if(child && parent){
	parent->Addrelative(child,rfactor[i]);
      }else{
	if(!child) cout<<bname[i]<<" not found, no relation with parent "<<relative[i]<<" established."<<endl;
	if(!parent) cout<<relative[i]<<" not found, no relation with child "<<bname[i]<<" established."<<endl;
      }
    }
  }
    
  //saving into a file
  SaveRoot(rooty,data[0],bgr[0],nbgr[0]);
  if(nsig[0]>0){
    SaveRoot(rooty,NULL,sig[0],nsig[0],false);
  }
  if(nlim[0]>0){
    SaveRoot(rooty,NULL,lim[0],nlim[0],false);
  }
  
  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  //***  USER PLOTTING, FITTING AND OTHER STUFF DOING CODE HERE*******************
  //******************************************************************************
  
  TApplication theApp("App", &argc, argv);
  ROOT_settings();
  Urun(data,bgr,sig,lim,nbgr,nsig,nlim,nchnls,namechnl);
  if(!batch) theApp.Run();
  //theApp.Terminate(0);
    
  }else{
    std::cout<<"\n\nCould not find file ==> "<<control<<std::endl;
    Usage();
    return 1;
  }
  
  anaout.close();
  erlog.close();
  
  return 0;
}


void Usage(){
  std::cout<<"\n============================================================================="<<std::endl;
  std::cout<<"Usage: ./ana <control.dat> <outfile.root> "<<std::endl;
  std::cout<<"<control.dat> = contains information about root files to load, format:"<<std::endl;
  std::cout<<"    [filename] [component name] [activity Bq] [MC evs] [tags]"<<std::endl;
  std::cout<<"       [filename] = filename with full path"<<std::endl;
  std::cout<<"       [component name] = what it is (bi214pmt, 2b2n etc....), should be unique"<<std::endl;
  std::cout<<"       [activity Bq] = component activity (for known background)/time for exp. data (s)"<<std::endl;
  std::cout<<"       [MC evs] = number of generated events in initial MC sample"<<std::endl;
  std::cout<<"       [tags] = describes what is inside"<<std::endl;
  std::cout<<"            D = data from experimant"<<std::endl;
  std::cout<<"            B = background MC"<<std::endl;
  std::cout<<"            R = reuse histrograms from previous rootfile"<<std::endl; 
  std::cout<<"            S = signal to fit"<<std::endl;
  std::cout<<"            L = signal to put limit on"<<std::endl;
  std::cout<<"<outfile.root> -- save all histograms here"<<std::endl;
  std::cout<<"    If left blank, name will be the same as control.dat name but with .root"<<std::endl;
  std::cout<<"=============================================================================\n"<<std::endl;
}

