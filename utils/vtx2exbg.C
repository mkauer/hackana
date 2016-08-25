#include "h10.h"
#include <fstream>
#include <iostream>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>


//For external background determine whether the vertex is in the
// internal tower, pos = 0
// external wall, pos = 1
// petals , pos = 2

Int_t h10::FindEXBGposition(Int_t trk){
  
  if(!trueMC_flag) return 0;// no true vertex stored

  Int_t pos = 100;
  
  Float_t Rtrue=sqrt(Xvntu[trk]*Xvntu[trk]+Yvntu[trk]*Yvntu[trk]); 
  Float_t Ztrue = Zvntu[trk];


  //internal ----------------------------------------------------

  if(Rtrue < 155) pos = 0;
  //extenal-------------------------------------------------------

  if(Rtrue > 155) pos = 1;

  //petal
  if(fabs(Ztrue)>125) pos = 2;

 if(pos==100){
   std::cout<<"pos "<<pos<<" R "<<Rtrue<<" Z "<<Ztrue<<std::endl;
   pos=0;
 }
				
 return pos;

}
