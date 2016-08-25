////////////////////////////////////////////////////////////
//
//  NEMO source information, and some useful routines     //
//  V.Vasiliev. V1.0 11.2004
////////////////////////////////////////////////////////////

#ifndef foils_h
#define foils_h

#include "TMath.h"
#include <stdio.h>
#include <iostream>
#include "h10.h"
#include <fstream>

// Names of sources and corresponding nuclide by source number 
char Snuclide[11][16]={"mo","mo","se","cd","nd","ca","zr","te","ten","cu","calib tube"};
char Sfoil[11][16]={"mo metal","mo composit","se composit","cd metal","nd","ca","zr","te-130 enriched","te natural","copper","calib tube"};
 
// Source mass [g] and isotope  mass by source number.
Float_t Smass[11]={2545.17,5578.29,1127.07,491.18,56.68,18.52,26.51,756.38,893.62,620.8,0.1};
Float_t Amass[11]={100,100,82,116,150,48,96,130,130,50,50};
Float_t Snmass[11]={2479.,4435.,934.6,405.,37.,6.99,9.4,453.9,166.34,620.8,0.1};


Float_t A2HalfLife(Float_t A, Int_t* isrc, Int_t nsrc){
  // This function converts activity into half0life for given source
  // Int_t -- array af sources numbers according tp Sfoil array
  // Returns T1/2 in years

  Float_t mtot=0;
  Int_t i;
  for(i=0;i<nsrc;i++){
    mtot+=Snmass[isrc[i]]/Amass[isrc[i]];
  }
  Float_t T=mtot/A*1.3232*1e16;
  return T;
}

Float_t h10::Calcsec(Float_t xvert, Float_t yvert)
{
  //this function returns vertex sector position from 0 to 20.   

  
  Float_t xa=xvert/TMath::Sqrt(xvert*xvert+yvert*yvert);
  xa=TMath::ACos(xa);
  if(yvert<0){
    xa=2.*TMath::Pi()-xa;
  }
  xa=20.*(xa/2./TMath::Pi());
  return xa;
}
Int_t h10::Decodesrc(Float_t xvert,Float_t yvert, Float_t zvert)
{
  // this function calculates source number by the vertex position
  //source codes: 0-- mo metal
  // 1 -- mo composit
  // 2 -- se composit
  // 3 -- cd
  // 4 -- nd
  // 5 -- ca
  // 6 -- Zr
  // 7 -- Te-130
  // 8 -- Te natural
  // 9 -- Cu foil
  // 10 -- Calibration tube (for calibration source) in the beginning of each sector
  // 99 -- no intersection with the source where found in the tracking.
  // 100 -- error, source not found
  // Use provided arrays for convinience:
  //                     Snuclide -- source nuclide name
  //                     Sfoil -- source  type name
  //                     Smass -- source foil mass
  //                     Snmass -- mass of isotope in the source

  if(xvert==0&&yvert==0&&zvert==0) return 99;

  Float_t xsecv=Calcsec(xvert,yvert);
  //check calibration tube
  if(xsecv-int(xsecv)<0.07){
    return 10;
  }
        //mo metal 
	if((xsecv>=1. &&xsecv <=1.736)||(xsecv>2.&&xsecv<5.)
	   ||(xsecv>=5.&&xsecv<=5.333)){ 
	  return 0;
	}
	//mo composit
        if((xsecv>=1.736&&xsecv<=2)||(xsecv>=5.333&&xsecv<=5.736)
	   ||(xsecv>=10.&&xsecv<=17.)){
	  return 1;
	}
	//Se 82 foils
  	if((xsecv>=7.&&xsecv<=8.0)||(xsecv>=8.&&xsecv<=8.2&&zvert<=-38.6))
	  {
	    //new se foils
	    return 2;
	  }
  	if((xsecv>=6.&&xsecv<=7.)||(xsecv>=8.&&xsecv<=8.2&&zvert>=-38.6)
	   ||(xsecv>=8.2&&xsecv<=8.33)){
	  //old se from NEMO 2
	  return 2;
	}
	//Cd 116 foils
	if(xsecv>=18.&&xsecv<=19.){
	  return 3;
	}
	//Nd 150 foils
	if(xsecv>=5.7371&&xsecv<=5.8706){
	  return 4;
	}
	// Ca-48
	if(xsecv>=5.8717&&xsecv<=6.0009){
	  if(zvert<=-5.&&zvert>=-65.){
	  return 5;
	  }
	}
	// Zr-96
	if(xsecv>=5.8717&&xsecv<=6.0009){ 
	if(zvert<=115.&&zvert>=2.){
	  return 6;
	}}
	//Te-130
	if((xsecv>=9.&&xsecv<=10.)||(xsecv>=17.&&xsecv<=18.)){ 
	  return 7;
        }
	//Te nat
	if((xsecv>=19.&&xsecv<=20.)||(xsecv>=8.33&&xsecv<=9.)){
	  return 8;
	}
	//CU
	if(xsecv>=0.&&xsecv<=1.){
	  return 9;
	}
	//		std::cout<<"Error decoding source (100) SEC:";
	//std::cout<<xsecv<<" X:Y:Z"<<xvert<<" "<<yvert<<" "<<zvert<<std::endl;
	return 100;
}


#endif
