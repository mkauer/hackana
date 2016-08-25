// Calculate event weight to adjust MC to individual activities per sector/gg layer
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


Float_t ana10::getActivityWeight(){
  //cout<<"getActivityWeight called, use? "<<swactivity_use<<endl;
  if(swactivity_use){
    return getswWeight(); 
  }else if(exactivity_use){
    return getexWeight();
  }
  else{
    return 1;
  }
}

Float_t ana10::getswWeight(){
  if(!swactivity_use) return 1;
  if(!swactivity_ini) initswActivity();
  Int_t sec, layer;
  if(swactivity_ini && swactivity_use){
    sec = FindTrueSector(0);
    layer = FindWireLayer(0);
    //    if(name == "swbi210") cout<<"swbi210 sector "<<sec<<" layer "<<layer<<" weight="<<swweight[sec][layer]<<endl;
    return swweight[sec][layer];
  }else{
    return 1;
  }
}

void ana10::initswActivity(){
  //cout<<"initswActivity called"<<endl;
  //cout<<"trueMC_flag = "<<trueMC_flag<<endl;
  if(swactivity_ini) return; // already initialised

  if(!swactivity_use){
    swactivity_ini = true;
    return;
  }
  
  if(!trueMC_flag){
    cout<<"No true MC information for the component "<<name<<endl;
    cout<<"Individual activity per sector/gg layer cant be used!"<<endl;
    cout<<"Will use mean activity "<<activity<<" from control.dat for the component "<<name<<endl;
    swactivity_use = false;
    swactivity_ini=true;
    return;
  }

  //read individual activities from the Vera's files
  Int_t io,layer,sector;
  Float_t a, ea;
  string file;
  if(stest(name,"bi210") && stest(name,"sw")){
    file="swbi210_activities.dat";
    cout<<"read "<<name<<" activities from the file "<<file<<endl;
    ifstream in0;
    in0.open(file.c_str());
    if(!in0){
      cout<<"File "<<file<<" with activities data cant be found"<<endl;
      cout<<"Will use mean activity "<<activity<<" from control.dat for the component "<<name<<endl;
      swactivity_use=false;
      swactivity_ini = true;
      return;
    }
    while(in0){
      in0 >> io >> layer>>sector>>a>>ea;
      swweight[sector][layer + io*10] = a/1000.;
      if (a==0){
	swweight[sector][layer + io*10] = swweight[sector][5 + io*10]/1000.;
      }
    }
    // Vera did not measured activities of the layers 0,6,7,8
    for(int s = 0;s<20;s++){
      if(swweight[s][0]==0) swweight[s][0] = swweight[s][1];
      if(swweight[s][10]==0) swweight[s][10] = swweight[s][11];
      if(swweight[s][6]==0) swweight[s][6] = swweight[s][5];
      if(swweight[s][7]==0) swweight[s][7] = swweight[s][5];
      if(swweight[s][8]==0) swweight[s][8] = swweight[s][5];
      if(swweight[s][16]==0) swweight[s][16] = swweight[s][15];
      if(swweight[s][17]==0) swweight[s][17] = swweight[s][15];
      if(swweight[s][18]==0) swweight[s][18] = swweight[s][15];
    }

    
    swactivity_ini = true;
  }
  if(stest(name,"bi214") && stest(name, "sw")){
    file="swrn_activities.dat";
    cout<<"read "<<name<<" activities from the file "<<file<<endl;
    ifstream in0;
    in0.open(file.c_str());
    if(!in0){
      cout<<"File "<<file<<" with activities data cant be found"<<endl;
      cout<<"Will use mean activity "<<activity<<" from control.dat for the component "<<name<<endl;
      swactivity_use=false;
      swactivity_ini = true;
      return;
    }
    while(in0){
      in0 >> io >> layer >> sector >>a>>ea;
      if(sector > 0) {
	swweight[sector-1][layer + io*10] = a/1000.;
	if (a==0){
	  swweight[sector-1][layer + io*10] = swweight[sector-1][5 + io*10]/1000.;
	}
      }
    }
    swactivity_ini = true;
  }
  if(stest(name,"pb214") && stest(name,"sw")){
    file="swrn_activities.dat";
    cout<<"read "<<name<<" activities from the file "<<file<<endl;
    ifstream in0;
    in0.open(file.c_str());
    if(!in0){
      cout<<"File "<<file<<" with activities data cant be found"<<endl;
      cout<<"Will use mean activity "<<activity<<" from control.dat for the component "<<name<<endl;
      swactivity_use=false;
      swactivity_ini = true;
      return;
    }
    while(in0){
      in0 >> io >> layer>> sector >>a>>ea;
      if(sector > 0){
	swweight[sector-1][layer + io*10] = a/1000.;
	if (a==0){
	  swweight[sector-1][layer + io*10] = swweight[sector-1][5 + io*10]/1000.;
	}
      }
    }
    swactivity_ini = true;
  }

  Float_t mactivity = getswMeanActivity();
  Float_t dbactivity, correction=1;;
  if(stest(name,"bi214") || stest(name,"pb214")){
    dbactivity = getRnActivityRunList();
    cout << "Mean Rn activity for current RunList is "<<getRnActivityRunList()<<" mBq "<<endl;

    //Vera's measurements were done for the 2b2n Phase 1 and Phase 2 runlists.
    //If different runlist is used, and also in the case of 2007 data correction is needed
    //Currently simply scale total activity according to Rn activity stored in the DB

    if(dbactivity > 800){
      correction = dbactivity / 1200.6;// during 2b2n Phase1 A = 1200.6
    }else{
      correction = dbactivity / 196.96;// during 2b2n Phase2 A = 196.96
    }
    if(activity == 0){
      cout<<"Reset component "<<name<<" activity to "<<dbactivity<<endl;;
      activity = dbactivity;
      mactivity = dbactivity;
    }

  }else{
    if(activity == 0){
      cout<<"Reset component "<<name<<" activity to "<<mactivity<<endl;;
      activity = mactivity;
    }
  }

  if(mactivity !=activity){
    cout<<"Warning! Component "<<name<<" activity in control.dat is"<<activity<<endl;
    cout<<"Mean activity from "<<file<<" is "<<mactivity<<endl;
    cout<<"Put this or 0 into control.dat"<<endl;
  }  
  

  for(Int_t s=0;s<20;s++){
    for(Int_t l=0;l<18;l++){
      swweight[s][l]*=(correction / mactivity);
    }
  }

  return;
}

Float_t ana10::getswMeanActivity(){
  Float_t a = 0;
  Float_t n = 0;
  for(Int_t s=0;s<20;s++){
    for(Int_t l=0;l<18;l++){
      a+=swweight[s][l];
      if(swweight[s][l]>0) n++;
    }
  }
  return a/n;
}
