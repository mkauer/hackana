#include "ana.hpp"
#include "TH2.h"
#include "TStyle.h"
#include "TNtuple.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TApplication.h"
#include "TGApplication.h"
#include "TLine.h"
#include "TPaveLabel.h"

#include <fstream>
#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <string>


void ana10::Loop(int eventnumber)
{
  // LOOP for analysis events selection. Author V. Vasiliev. 01.2006
  // This function creates histogram, loop through the chain and fills histogram for further displaying/analysis
  //Should be provided by user

  //technical variables 


  Int_t results[20]={0};// Events passage control, 0-10 slim, 11-20 ana10
  
   Int_t nbytes = 0;
   Int_t counter = 0;
   Int_t nentries = Int_t(fChain->GetEntries());
   Int_t nb = -1;
   int treenum=-1;

   //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   //
   //    USER DEFINITION FOR EVENT SELECTION PARAMETERS
   //
   //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   //------- PUT YOUR CUTS HERE --------------------------------
   Float_t Etot_min=0.2;


   //Set desired run period
   // see ana10/utils.hpp for the defenition of the constants
   SetPeriod(ana::P1,true);
   SetPeriod(ana::P2a,false);
   SetPeriod(ana::P2b,false);

   SetRuns(ana::GOOD_RUN,false);
   SetRuns(ana::OK_RUN,false);
   SetRuns(ana::NOAIR_RUN,false);
   SetRuns(ana::AFTER_EC_RUN,false);
   SetRuns(ana::HV_RUN,false);
   SetRuns(ana::ALL_BUT_BAD_RUN,false);
   SetRuns(ana::STANDARD_RUN,true);

   //--------End of CUTS----------------------------------------

   //-------- BOOK YOUR HISTOGRAMS HERE-------------------------
   Addh1("etot","Total energy",100,0,4.2);
   //---------End of BOOK---------------------------------------

   //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   //
   //            END OF USER SECTION
   //
   //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   if (fChain == 0) return;
   std::cout<<"Number of entries in the chain "<<nentries<<std::endl;

   //main selection LOOP, fills booked hystos
   for (Int_t jentry=eventnumber; jentry<nentries;jentry++) {

       Bool_t preselect = false;
       Int_t ientry = LoadTree(jentry);
       Int_t nb = fChain->GetEntry(jentry);

       Int_t i,j,k;
       if (ientry < 0) break;
	results[0]++;

	if(jentry%50000==0) std::cout<<jentry<<" events passed "<<std::endl;

        //It is highly recommended to call correct SLIM selection function (same as was used in slim2 program), to ensure that even unslimed root files anaysed properly with every the same cuts applied. If you don't want to call slim Cut(), provide dummy vertion of this function, since its presence is obligatory.

	if(!Cut(results)) continue;


	//check that run is in the runlist 
	//save selected runs duration for normalization of MC later.	
	Int_t iitime;
	Bool_t goodrun=RunList(abs(run),iitime);
	if(!goodrun) continue;
	AddRunTime(abs(run),iitime);


	//Check PMT status. Ignore events w/o LTC correction
	Bool_t goodpms=true;
	for(Int_t k=0;k<Nsc;k++){
	  Int_t pmstat=pmstatus->get(run,Sc[k][1],Sc[k][2],Sc[k][3],Sc[k][4]);
	  goodpms = goodpms && GoodPMStatus(pmstat);
       	}
	if(!goodpms) continue;

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//
//             START OF USER DEFINED SECTION
//
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

//!!!!!!!!!!!!!ENTER YOUR SELECTION CODE HERE!!!!!!!!!!!!!!!!

	Float_t Etot=0;
	for(Int_t i=0;i<Nsc;i++){
	  Etot+=Sc[i][8]*1000;
	}
	if(Etot<Etot_min) continue;
	results[11]++;



        // attach weight to each MC event. 
	Float_t w;
	if(run<0){
	  w= GetEventWeight();// MC
	}else{
	  w=1.; //REAL
	}

//!!!!!!!!!!!!!!!!!FILL HISTOGRAMS HERE!!!!!!!!!!!!!!!!!!!!!!

	//Fill histograms for sucessfull events
	Geth1("etot")->Fill(Etot,w);
	


//!!!!!!!!!!!!!END OF USER SECTION!!!!!!!!!!!!!!!!!!!!!!!!!!
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	counter++;
   }
   std::cout<<"number of events:  "<<counter<<std::endl;
   for(Int_t i=0;i<20;i++)
     { 
       if(results[i]) std::cout<<" results for cut "<<i<<" "<<results[i]<<std::endl;
     } 
}
