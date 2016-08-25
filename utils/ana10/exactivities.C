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


Float_t ana10::getexWeight(){
  if(!exactivity_use) return 1;
  if(!exactivity_ini) initexActivity();
  Int_t pos;// internal, outer or petal
  if(exactivity_ini && exactivity_use){
    if(trueMC_flag) {
      pos = FindEXBGposition(0);
    }else{
      pos = 0;// use 0 component if true vtx is not preserved
    }
    return exweight[pos];
  }else{
    return 1;
  }
}

void ana10::initexActivity(){
  //cout<<"initexActivity called"<<endl;
  //cout<<"trueMC_flag = "<<trueMC_flag<<endl;
  if(exactivity_ini) return; // already initialised

  if(!exactivity_use){
    exactivity_ini = true;
    return;
  }
  
  if(!trueMC_flag && !noextruemc_warned){
    cout<<"No true MC information for the component "<<name<<endl;
    cout<<"Individual activity for I/O/Petal can't be used "<<endl;
    noextruemc_warned = true;
  }

  //read individual activities from the Victor F model

  exactivity_ini = true;
  exactivity_use = false; // do not use this mechanism if component is not listed in exbg_modelf.dat
 
  Float_t ai,ao,ap;
  char mname[12];
  string file;
  if(stest(name,"ex") || stest(name,"exbg") || stest(name,"scin") || stest(name,"ss") || stest(name,"sc") || stest(name,"pm")){
    file = "exbg_modelf.dat";
    cout<<"\t LOOKING FOR CONTROL FILE ==> "<<file<<endl;
    cout<<"read "<<name<<" activities from the file "<<file<<endl;
    ifstream in0;
    in0.open(file.c_str());
    if(!(in0.is_open())) cout<<"\t COULD NOT OPEN CONTROL FILE ==> "<<file<<endl;
    while(in0){
      in0 >> mname >> ai >> ao >> ap;
      //cout << "Component read is "<<mname<<" "<<ai<<" "<<ao<<" "<<ap<<endl;
      if(stest(mname,"exbg2") && (stest(name,"ex2") || stest(name,"exbg2"))){
	//exbg2 Bi214
	if(stest(name,"bi") && stest(mname,"bi214")) {
	  exweight[0]=ai;
	  exweight[1]=ao;
	  exweight[2]=ap;
	  exactivity_use=true;
	  exactivity_ini = true;
	}
	//exbg2 Pb214
	if(stest(name,"pb") && stest(mname,"pb214")) {
	  exweight[0]=ai;
	  exweight[1]=ao;
	  exweight[2]=ap;
	  exactivity_use=true;
	  exactivity_ini = true;
	}
	//exbg2 Tl208
	if(stest(name,"tl") && stest(mname,"tl208")) {
	  exweight[0]=ai;
	  exweight[1]=ao;
	  exweight[2]=ap;
	  exactivity_use=true;
	  exactivity_ini = true;
	}
      }
      if(stest(mname,"exbg4") && (stest(name,"ex4") || stest(name,"exbg4"))){
	//exbg4 Bi214
	if(stest(name,"bi") && stest(mname,"bi214")) {
	  exweight[0]=ai;
	  exweight[1]=ao;
	  exweight[2]=ap;
	  exactivity_use=true;
	  exactivity_ini = true;
	}
	//exbg4 Ac228
	if(stest(name,"ac") && stest(mname,"ac228")) {
	  exweight[0]=ai;
	  exweight[1]=ao;
	  exweight[2]=ap;
	  exactivity_use=true;
	  exactivity_ini = true;
	}
	//exbg4 Tl208
	if(stest(name,"tl") && stest(mname,"tl208")) {
	  exweight[0]=ai;
	  exweight[1]=ao;
	  exweight[2]=ap;
	  exactivity_use=true;
	  exactivity_ini = true;
	}
      }
      if(stest(mname,"exbg11") && (stest(name,"ex11") || stest(name,"exbg11") || stest(name,"pm"))){
	//exbg11 Bi214
	if(stest(name,"bi") && !stest(name,"bi212") && stest(mname,"bi214")) {
	  exweight[0]=ai;
	  exweight[1]=ao;
	  exweight[2]=ap;
	  exactivity_use=true;
	  exactivity_ini = true;
	}
	//exbg11 Pb214
	if(stest(name,"pb") && stest(mname,"pb214")) {
	  exweight[0]=ai;
	  exweight[1]=ao;
	  exweight[2]=ap;
	  exactivity_use=true;
	  exactivity_ini = true;
	}
	//exbg11 Tl208
	if(stest(name,"tl") && stest(mname,"tl208")) {
	  exweight[0]=ai;
	  exweight[1]=ao;
	  exweight[2]=ap;
	  exactivity_use=true;
	  exactivity_ini = true;
	}
	//exbg11 Ac228
	if(stest(name,"ac") && stest(mname,"ac228")) {
	  exweight[0]=ai;
	  exweight[1]=ao;
	  exweight[2]=ap;
	  exactivity_use=true;
	  exactivity_ini = true;
	}
	//exbg11 K40
	if(stest(name,"k") && stest(mname,"k40")) {
	  exweight[0]=ai;
	  exweight[1]=ao;
	  exweight[2]=ap;
	  exactivity_use=true;
	  exactivity_ini = true;
	}
      }
      if(stest(mname,"exbg7") && (stest(name,"ex7") || stest(name,"exbg7"))){
	//exbg7 Bi214
	if(stest(name,"bi") && stest(mname,"bi214")) {
	  exweight[0]=ai;
	  exweight[1]=ao;
	  exweight[2]=ap;
	  exactivity_use=true;
	  exactivity_ini = true;
	}
	//exbg7 Pb214
	if(stest(name,"pb") && stest(mname,"pb214")) {
	  exweight[0]=ai;
	  exweight[1]=ao;
	  exweight[2]=ap;
	  exactivity_use=true;
	  exactivity_ini = true;
	}
	//exbg7 Tl208
	if(stest(name,"tl") && stest(mname,"tl208")) {
	  exweight[0]=ai;
	  exweight[1]=ao;
	  exweight[2]=ap;
	  exactivity_use=true;
	  exactivity_ini = true;
	}
	//exbg7 Ac228
	if(stest(name,"ac") && stest(mname,"ac228")) {
	  exweight[0]=ai;
	  exweight[1]=ao;
	  exweight[2]=ap;
	  exactivity_use=true;
	  exactivity_ini = true;
	}
	//exbg7 K40
	if(stest(name,"k") && stest(mname,"k40")) {
	  exweight[0]=ai;
	  exweight[1]=ao;
	  exweight[2]=ap;
	  exactivity_use=true;
	  exactivity_ini = true;
	}
	//exbg7 Co60
	if(stest(name,"co") && stest(mname,"co60")) {
	  exweight[0]=ai;
	  exweight[1]=ao;
	  exweight[2]=ap;
	  exactivity_use=true;
	  exactivity_ini = true;
	}
	//exbg7 pa234m
	if(stest(name,"pa") && stest(mname,"pa234m")) {
	  exweight[0]=ai;
	  exweight[1]=ao;
	  exweight[2]=ap;
	  exactivity_use=true;
	  exactivity_ini = true;
	}
      }
      if(stest(mname,"exbg6") && (stest(name,"ex6") || stest(name,"exbg6"))){
	//exbg6 Co60
	if(stest(name,"co") && stest(mname,"co60")) {
	  exweight[0]=ai;
	  exweight[1]=ao;
	  exweight[2]=ap;
	  exactivity_use=true;
	  exactivity_ini = true;
	}
      }
      if(stest(mname,"exbg17") && (stest(name,"ex17") || stest(name,"exbg17"))){
	//exbg17 Co60
	if(stest(name,"co") && stest(mname,"co60")) {
	  exweight[0]=ai;
	  exweight[1]=ao;
	  exweight[2]=ap;
	  exactivity_use=true;
	  exactivity_ini = true;
	}
      }
      if(stest(mname,"exbg9") && (stest(name,"ex9") || stest(name,"exbg9"))){
	//exbg9 Co60
	if(stest(name,"co") && stest(mname,"co60")) {
	  exweight[0]=ai;
	  exweight[1]=ao;
	  exweight[2]=ap;
	  exactivity_use=true;
	  exactivity_ini = true;
	}
      }
      if(stest(mname,"exbg8") && (stest(name,"ex8") || stest(name,"exbg8"))){
	//exbg8 Co60
	if(stest(name,"co") && stest(mname,"co60")) {
	  exweight[0]=ai;
	  exweight[1]=ao;
	  exweight[2]=ap;
	  exactivity_use=true;
	  exactivity_ini = true;
	}
	//exbg8 Pa234m
	if(stest(name,"pa") && stest(mname,"pa234m")) {
	  exweight[0]=ai;
	  exweight[1]=ao;
	  exweight[2]=ap;
	  exactivity_use=true;
	  exactivity_ini = true;
	}
      }
      if(stest(mname,"sci") && (stest(name,"sc") || stest(name,"sci"))){
	//sci K40
	if(stest(name,"k") && stest(mname,"k40")) {
	  exweight[0]=ai;
	  exweight[1]=ao;
	  exweight[2]=ap;
	  exactivity_use=true;
	  exactivity_ini = true;
	}
      }
      if(stest(mname,"ss") && stest(name,"ss")){
	//ss bi210
	if(stest(name,"bi210") && stest(mname,"bi210")) {
	  exweight[0]=ai;
	  exweight[1]=ao;
	  exweight[2]=ap;
	  exactivity_use=true;
	  exactivity_ini = true;
	}
	//ss bi214
	if(stest(name,"bi214") && stest(mname,"bi214")) {
	  exweight[0]=ai;
	  exweight[1]=ao;
	  exweight[2]=ap;
	  exactivity_use=true;
	  exactivity_ini = true;
	}
	//ss pb214
	if(stest(name,"pb") && stest(mname,"pb214")) {
	  exweight[0]=ai;
	  exweight[1]=ao;
	  exweight[2]=ap;
	  exactivity_use=true;
	  exactivity_ini = true;
	}
      }
    }
    Float_t mactivity = getexMeanActivity();
    if(activity == 0){
      cout<<"Reset component's "<<name<<" activity to "<<mactivity<<endl;;
      activity = mactivity;
    }
    if(mactivity !=activity){
      cout<<"Warning! Component "<<name<<" activity in control.dat is "<<activity<<endl;
      cout<<"Mean activity from "<<file<<" is "<<mactivity<<endl;
      cout<<"Put this or 0 into control.dat"<<endl;
    }  
    for(Int_t s=0;s<3;s++){
      if (mactivity != 0) exweight[s]/=mactivity;
      //cout<<"Set weight for component "<<name<<" iop="<<s<<" W="<<exweight[s]<<endl;
    }
    in0.close();
  }
  return;
}

Float_t ana10::getexMeanActivity(){
  Float_t a = 0;
  Float_t n = 0;
  for(Int_t s=0;s<3;s++){
      a+=exweight[s];
      if(exweight[s]>0) n++;
  }
  if(n==0) return 0;
  return a/n;
}
