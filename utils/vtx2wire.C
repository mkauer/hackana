#include "h10.h"
#include <fstream>
#include <iostream>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>


//Find that which layer of the tracker the event is coming from 
//N.Fatemi-Ghomi 
Int_t h10::FindWireLayer(Int_t trk){
Int_t layer=100;
Float_t Rtrue=sqrt(Xvntu[trk]*Xvntu[trk]+Yvntu[trk]*Yvntu[trk]);  
// Float_t Rtrue=X_foil[0];
//internal ----------------------------------------------------
//layers 0,1,2,3
 if(Rtrue<153.6 && Rtrue>150.7) layer=0;
 if(Rtrue<150.7 && Rtrue>147.9) layer=1;
 if(Rtrue<147.9 && Rtrue>145.1) layer=2;
 if(Rtrue<145.1 && Rtrue>142.2) layer=3;
 //layers 4,5
 if(Rtrue<128.6 && Rtrue>125.7) layer=4;
 if(Rtrue<125.7 && Rtrue>122.8) layer=5;
 //layers 6,7,8
 if(Rtrue<109.2 && Rtrue>106.3) layer=6;
 if(Rtrue<106.3 && Rtrue>103.5) layer=7;
 if(Rtrue<103.5 && Rtrue >100.6)layer=8;
 //extenal-------------------------------------------------------
 //layers 10,11,12,13
 if(Rtrue<159.3 && Rtrue>156.4) layer=10;
 if(Rtrue<162.1 && Rtrue>159.3) layer=11;
 if(Rtrue<164.9 && Rtrue>162.1) layer=12;
 if(Rtrue<167.8 && Rtrue>164.9) layer=13;
 //layers 14,15
 if(Rtrue<184.3 && Rtrue>181.4) layer=14;
 if(Rtrue<187.2 && Rtrue>184.3) layer=15;
 //layer 16,17,18
 if(Rtrue<203.7 & Rtrue>200.8) layer=16;
 if(Rtrue<206.5 && Rtrue>203.7) layer=17;
 if(Rtrue<209.304 && Rtrue>206.5) layer=18; 
 if(layer==100){
   std::cout<<"layer "<<layer<<" R "<<Rtrue<<std::endl;
   layer=0;
 }
				
 return layer;

}

Int_t h10::FindTrueSector(Int_t trk){
  // this function calculates sector of a true MC vertex
  // author V. Vasiliev 01/2008

  Float_t theta = atan2(Yvntu[trk],Xvntu[trk]);
  theta = theta/3.1415926*10;
  if(theta > 0){
    return (Int_t) theta;
  }else{
    theta +=20;
    return (Int_t) theta;
  }
}
