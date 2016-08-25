#define h10_cxx
#define h10_cut

#include "h10.h"
#include "TH2.h"
#include "TStyle.h"
#include "TNtuple.h"
#include "TCanvas.h"
#include "TMath.h"

#include <fstream>
#include <sstream>
#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <string>

 Long64_t nevgamc[100];
 Int_t mctype[100];
 Int_t filemctype[15000];
 Int_t ntypes=0;
 Long64_t date_orig;

//Pmstatus *pmstatus;

using namespace std;
int main(int argc, char* argv[]){
  Long64_t timetot=0;
  Long64_t longdummy=0;
  Bool_t mc=false;


 Bool_t filelist=true;
 TTree* tree=0;



  //init Pmstatus cache
  //pmstatus = new Pmstatus();

  
  //use UCL mirror server
  //nemo3::DbMgr::Init("nemodb.hep.ucl.ac.uk",3306);  // old ucl server
  //nemo3::DbMgr::Init("nemodb1.hep.ucl.ac.uk",3306);   // new ucl server
  // nemo3::DbMgr::Init("ccnemodb.in2p3.fr",13306);
  //  nemo3::DbMgr::Init("zion",3306);

 std::cout<<"Slim program. Preselects events."<<std::endl;
 std::cout<<"Author V. Vasiliev 01.2006"<<std::endl;
 std::cout<<"Usage slim2 [-c  <control file>] [-r <root file>] [-n <event number>] -new/old"<<std::endl;
 cout<< "-c option -- control file with reconstructed root ntuples to preselect from, control.dat by default"<<std::endl;
 cout<< "-r option -- root file to save results, rooty.root by default"<<std::endl;
 cout<< "-n option -- start from the event n, 1 by default"<<std::endl;
 cout<< "-old/new switch -- old ntuples, perform estimation of MCEVGA, not write to the output ntuple; new ntuples, MCEVGA stored in the ntuples, sum and copy it to the results ntuple."<<std::endl;

 //Analyse command line
 char controlfile[256]="control.dat";
 char rootfile[256]="rooty.root";
 Long64_t nstart=0;

 Bool_t fnew=false;
 for(Int_t i=1;i<argc;i++){
   Int_t inc=0;
   if(string(argv[i])=="-c"){
       strcpy(controlfile,argv[i+1]);
       inc=1;
   }
   if(string(argv[i])=="-r"){
       strcpy(rootfile,argv[i+1]);
       inc=1;
   }
   if(string(argv[i])=="-n"){
     stringstream oss;
     oss<<argv[i+1];
     oss>> nstart;
     inc=1;
   }
   if(string(argv[i])=="-new"){
       fnew=true;
   }
   i+=inc;
 }

 Int_t treecount = 0;
 if (filelist)  
   {  
      char control[256];
      Int_t numFiles;
      ifstream in0, in1, in2;

      in0.open(controlfile);
      char filenames[15000][256];
      Int_t runnumber[15000], datenumber[15000];
      

      Int_t i=0;
      while(in0 >> filenames[i]){
	i++;
      }
      numFiles=i;
      timetot=0;
      std::cout<<numFiles<<" files to process"<<std::endl;
      if (tree==0)
        {  TChain *chain = new TChain("h10");
	   Int_t runnumber;
	   Int_t timet;
	   Int_t date;
	   chain->SetBranchAddress("run",&runnumber);
	   chain->SetBranchAddress("time",&timet);
	   chain->SetBranchAddress("date",&date);
           for (int i=0;i<numFiles;i++) 	     
	     { 
	       std::cout<<" going to add "<<filenames[i]<<std::endl;
	       Int_t errcode = chain->AddFile(filenames[i],-1);	      
	       if(errcode == 1){
		
		 std::cout<<"File is added "<<filenames[i]<<std::endl;
		 Int_t nentries = Int_t(chain->GetEntries());
		 chain->GetEntry(nentries-1);
		 if(i==0) date_orig = date;
		 fnew=(runnumber<0);
		 if(timet==0 || runnumber<0){		  
		   // mc files
		   mc=true;
		   std::cout<<"MC file, run "<<runnumber;
		   if(timet==0){		
		     std::cout<<", no NEVGA info available. Estimate later."<<std::endl;
		   }else{
		     //detect mc runlist type
		     longdummy=timet;
		     if(longdummy<0) longdummy*=-1000;
		     Bool_t exist=false;
		     for(Int_t i1=0;i1<ntypes;i1++){
		       if(abs(runnumber-mctype[i1])<500){
			 nevgamc[i1]+=longdummy;
			 filemctype[treecount]=mctype[i1];
			 exist = true;
		       }
		     }
		     if(!exist){
		       mctype[ntypes]=runnumber;
		       nevgamc[ntypes]=longdummy;
		       filemctype[treecount]=runnumber;
		       ntypes++;
		     }
		     timetot+=longdummy;
		     std::cout<<", generated "<<timet<<" events type "<<filemctype[i]<<std::endl;
		   }
		 }else{
		   //data files
		   mc=false;
		   longdummy=timet;
		   if(longdummy<0) longdummy*=-1000;
		   std::cout<<"DATA file, run "<<runnumber;
		   std::cout<<" time "<<longdummy<<" s"<<std::endl;
		   timetot+=longdummy;
		 }
		 treecount++;
	       }	     
	     }
	   tree=chain;
	}
   }

  h10 george(tree);
  if(!fnew && mc ){
    george.timetot=0;
  }else{
    george.timetot=-timetot/1000;
  }
  strcpy(george.rfile,rootfile);

  if(mc){
      if(fnew){
       std::cout<<"Total MC events in the slim file "<<timetot/1000<<" E3"<<std::endl;
       cout<<"Saved in the slim tree!"<<endl;
      }
  }else{
    std::cout<<"Total raw data time in the slim file "<<timetot/1000<<"E3 s"<<std::endl;
  }
  george.Loop(nstart);

  if(mc){
      std::cout<<"Estimated MC events in the slim "<<george.timetot<<"E3."<<std::endl;
  }

}


void h10::Loop(Long64_t eventnumber)
{
   Float_t c=3e8;
   char control[256];
   Int_t numFiles;
   ifstream in0;
   Int_t mcevents=0;
   Int_t mcevc=0;
   

   TFile *newfile = new TFile(rfile,"RECREATE");

   if (fChain == 0) return;

   fChain->SetCacheSize(0);
   
   // Events passage control, 0-19 preselect, 20-99 slim, 100+ ana10
   Int_t results[100]={0};
   Long64_t nbytes = 0;
   
   TTree* newh10 = (TTree*)fChain->GetTree()->CloneTree(0);
   Long64_t counter = 0;
   Long64_t nentries = fChain->GetEntries();
   Long64_t nb = -1;
   Long64_t treenum=-1;
   std::cout<<"Number of entries in the chain "<<nentries<<std::endl;

       for (Long64_t jentry=eventnumber; jentry<nentries;jentry++) {
       Bool_t preselect = false;
       Long64_t ientry = LoadTree(jentry);
       Long64_t i,j,k;
       if (ientry < 0) break;
       if (treenum != fChain->GetTreeNumber()) 
         { treenum = fChain->GetTreeNumber();
           fChain->CopyAddresses(newh10);
         }

        nb = fChain->GetEntry(jentry);
        nbytes += nb;
	//results[0]++;
	//check that the date is correct
	/**
	   if(date!=date_orig and run<0){
	   std::cout<<"WARNING!! date = "<< date<<" for the file "<<fChain->GetDirectory()->GetName();
	   date=date_orig;
	   std::cout<<" new date = "<<date<<std::endl;
	   }
	**/
	//If non-zero time was found in the trees during loading -- save it!
	for(Int_t i1=0;i1<ntypes;i1++){
	  if(mctype[i1]==filemctype[treenum] and run<0){
	    time=-(Int_t)(nevgamc[i1]/1000);
	  }
	}
	if(timetot!=0 and run>0){
	  time=timetot;
	}
	// And do estimation for MC
	  if(Myievent<mcevc) {
	    Double_t divider = TMath::Power(10,round(TMath::Log10(mcevc)-2));
	    mcevents+=Int_t((mcevc/divider + 1.)*divider/1000);
	    
	  }
	  mcevc=Myievent;
	
	  if(jentry%10000==0) std::cout<<"Processing event "<<jentry<<" file "<<fChain->GetDirectory()->GetName()<<std::endl;
	  if(jentry%10000==0) newh10->AutoSave();
	  
	  //Bool_t preselect=preselect_1tk(jentry,results);
	  //Bool_t preselect=preselect_1e1g(jentry,results);
	  //Bool_t preselect=preselect_1e(jentry,results);
	  //if(!preselect) continue;  
	  
	  if(Cut(results)==1){
	    counter++;
	    newh10->Fill();
	  }
	  
       }
       mcevents+=mcevc/1000;
       timetot=mcevents;
       std::cout<<"number of events selected:  "<<counter<<std::endl;
       std::cout<<" now writing "<<std::endl;
       newh10->Write();
   
   newfile->Write();
   newfile->Close(); 
   for(Int_t i=0;i<99;i++)
     { 
       if(results[i]) std::cout<<" results for cut "<<i<<" "<<results[i]<<std::endl;
     } 
}

