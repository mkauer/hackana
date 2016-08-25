#include "h10.h"
#include <fstream>
#include <iostream>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>


void h10::testmodule()
{ 
  std::cout<<" hello world "<<std::endl;
}

void h10::assign()
{ 
  for (Int_t i=0; i<Nbr_tks; i++) 
    {Int_t assign = 0;  
     for (Int_t j=0; j<Nsc; j++)
      { Float_t deltax = (X_scintil[i]-Sc[j][9]);
        Float_t deltay = (Y_scintil[i]-Sc[j][10]);
        Float_t deltaz = (Z_scintil[i]-Sc[j][11]);
        Float_t delta = sqrt(deltax*deltax+deltay*deltay+deltaz*deltaz);
        if (delta<10)
          { Ind_scintil[i]=j;
            assign++;
          }
      }
     if(assign>1)std::cout<<" error: more than one track wants this scint"<<std::endl;
    }
}



Float_t   h10::TOF2b_int(Float_t &tf, Float_t &sf,Int_t tk1, Int_t tk2)
{
  //Calculates TOF prpobability for internal hypothesis for two electrons
  // Retruns Chi2 prob.
  Float_t prob=1.E-30;
  // Theoretical TOF for internal hypothesis
  Float_t tofb[2];
  Float_t dtofb[2];
  //Measured TOF for internal hypothesis
  Float_t tofbm[2];
  Float_t rme=0.511;
  Float_t vlight=3.E10, dlen;
  Int_t tk[2];
  tk[0]=tk1;
  tk[1]=tk2;

  for(Int_t j=0;j<2;j++){
    Int_t isc=Ind_scintil[tk[j]]-1;
    Float_t eb=Sc[isc][8]*1000.;
    Float_t ve=sqrt(eb*eb+2.*rme*eb)/(rme+eb)*vlight;
    Float_t dstrack=gethelixl(tk[j],dlen);
    Float_t fwhm;
    Float_t sige;
    tofb[j]=dstrack/ve*1.E9;
    tofbm[j]=getsctimebb(isc,sige);
    dtofb[j]=sige;

    if(Sc[isc][2]==1){
      fwhm=0.14;
    }else{
      fwhm=0.17;
    }
    Float_t se1 = fwhm*sqrt(eb)/2.354; 
    sige = tofb[j]*rme*rme/eb/(eb+rme)/(eb+2*rme)*se1;
    dtofb[j]=sqrt(sige*sige+dtofb[j]*dtofb[j]);
  }

  Float_t dtt=tofbm[1]-tofbm[0]-(tofb[1]-tofb[0]);
  Float_t chiint2=dtt*dtt/(dtofb[1]*dtofb[1]+dtofb[0]*dtofb[0]);
  prob = TMath::Prob(chiint2,1); 
  tf=(tofbm[0]-tofb[0])/(dtofb[0]*dtofb[0])+(tofbm[1]-tofb[1])/(dtofb[1]*dtofb[1]);
  sf=sqrt(1./(dtofb[0]*dtofb[0])+1/(dtofb[1]*dtofb[1]));
  tf=tf/(sf*sf);
  sf=1./sf;
  //  std::cout<<"Run "<<run<<" event "<<Myievent<<std::endl;
  // std::cout<<"T2-T1 "<<tofbm[1]-tofbm[0]<<" TOF1,2 "<<tofb[0]<<","<<tofb[1]<<" S1,2 "<<dtofb[0]<<","<<dtofb[1]<<std::endl;
  //std::cout<<"dtt="<<dtt<<" chi2="<<chiint2<<std::endl;
  return prob;  
}

Float_t   h10::TOF2b_ext(Int_t tk1, Int_t tk2)
{
  //Calculates TOF prpobability for external hypothesis for two electrons
  // Retruns Chi2 prob.
  Float_t prob=1.E-30;
  // Theoretical TOF for external hypothesis
  Float_t tofb[2];
  Float_t dtofb[2];
  //Measured times 
  Float_t tofbm[2];
  Float_t rme=0.511;
  Float_t vlight=3.E10;
  Float_t fwhm, dlen;
  Float_t eb,sige;

  Int_t tk[2];
  tk[0]=tk1;
  tk[1]=tk2;

  //Estimate energy losses in the foil
    Float_t cos_dir1=(Cos_dir[tk[0]][0]*Xvert+Cos_dir[tk[0]][1]*Yvert)/sqrt(Xvert*Xvert+Yvert*Yvert);
    Float_t cos_dir2=(Cos_dir[tk[1]][0]*Xvert+Cos_dir[tk[1]][1]*Yvert)/sqrt(Xvert*Xvert+Yvert*Yvert);
    Float_t floss=0.037/fabs(cos_dir1)+0.037/fabs(cos_dir2);
  

  //1st hypothesis 0->1
  for(Int_t j=0;j<2;j++){
    Int_t isc=Ind_scintil[tk[1]]-1;
    eb=Sc[isc][8]*1000.;
    //    if(j==0) {std::cout<<" eb before "<<eb<<std::endl;eb+=floss; std::cout<<" adding floss"<<floss<<" eb after "<<eb<<std::endl;}

    Float_t ve=sqrt(eb*eb+2.*rme*eb)/(rme+eb)*vlight;
    Float_t dstrack=gethelixl(tk[j],dlen);
    tofb[j]=dstrack/ve*1.E9;
    tofbm[j]=getsctimebb(Ind_scintil[tk[j]]-1,sige);
    dtofb[j]=sige;

  }
  if(Sc[Ind_scintil[tk[1]]-1][2]==1){
     fwhm=0.14;
  }else{
     fwhm=0.17;
  }
  Float_t se1 = fwhm*sqrt(eb)/2.354; 
  sige = (tofb[0]+tofb[1])*rme*rme/eb/(eb+rme)/(eb+2*rme)*se1;

  Float_t dtt=tofbm[1]-tofbm[0]-(tofb[1]+tofb[0]);
  //  std::cout<<" Vlad1 "<<tofbm[1]<<" "<<tofbm[0]<<" "<<tofb[1]<<" "<<tofb[0]<<" energy "<<eb<<" floss "<<floss<<std::endl;
  Float_t chiint2=dtt*dtt/(dtofb[1]*dtofb[1]+dtofb[0]*dtofb[0]+sige*sige);
  //    std::cout<<"TOFextmy dtt1 "<<dtt<<std::endl;
  Float_t prob1 = TMath::Prob(chiint2,1); 


  //2nd hypothesis 1->0
  for(Int_t j=0;j<2;j++){
    Int_t isc=Ind_scintil[tk[0]]-1;
    eb=Sc[isc][8]*1000.;
    //    if(j==1) {std::cout<<" eb before "<<eb<<std::endl;eb+=floss; std::cout<<" adding floss"<<floss<<" eb after "<<eb<<std::endl;}
    Float_t ve=sqrt(eb*eb+2.*rme*eb)/(rme+eb)*vlight;
    Float_t dstrack=gethelixl(tk[j],dlen);
    tofb[j]=dstrack/ve*1.E9;
    tofbm[j]=getsctime(Ind_scintil[tk[j]]-1,sige);

    dtofb[j]=sige;
  }

  if(Sc[Ind_scintil[tk[0]]-1][2]==1){
    fwhm=0.14;
  }else{
    fwhm=0.17;
  }
  se1 = fwhm*sqrt(eb)/2.354; 
  sige = (tofb[0]+tofb[1])*rme*rme/eb/(eb+rme)/(eb+2*rme)*se1;

  dtt=tofbm[0]-tofbm[1]-(tofb[1]+tofb[0]);
  //std::cout<<" Vlad2 "<<tofbm[1]<<" "<<tofbm[0]<<" "<<tofb[1]<<" "<<tofb[0]<<" energy "<<eb<<" floss "<<floss<<std::endl;
  //    std::cout<<"TOFextmy dtt2 "<<dtt<<std::endl;
  chiint2=dtt*dtt/(dtofb[1]*dtofb[1]+dtofb[0]*dtofb[0]+sige*sige);
  Float_t prob2 = TMath::Prob(chiint2,1); 
  if(prob1<prob2) prob=prob2;
  if(prob2<prob1) prob=prob1;

  return prob;
}


Float_t h10::getsctime(Int_t isc,Float_t &sigma)
{
  // This function returns time of PMT signal and its error (sigma)



  return getsctimebb(isc,sigma);

  // Previous version was needed to treat beta and gamma with different decay constant
  // this problem was solved in recent time calibrations produced (started form nemos8)
  // so this code become obsolete now.
  /**
  Float_t t=0.;
  Float_t fwhm;
  if(isc>=Nsc) return 0.;
    t=-Sc[isc][6]*0.053;
    if(Sc[isc][2]==1){
      fwhm=0.14;
    }else{
      fwhm=0.17;
    }
    if(run>1000){
      if(Sc[isc][2]==0) fwhm=grndm.Gaus(0.173,0.014);
      if(Sc[isc][2]==1) fwhm=grndm.Gaus(0.142,0.015);
      if(Sc[isc][2]>1){
	if(Sc[isc][3]<3) fwhm=grndm.Gaus(0.163,0.022);
	if(Sc[isc][3]>2) fwhm=grndm.Gaus(0.146,0.024);
      }
    }
  Float_t phen=(Sc[isc][8]*1000.)*(2.3548/fwhm)*(2.3548/fwhm);
  //-- sigma from transit time fluctuation and TDC convertion    
  
 
  Float_t tau_sc=6.;
  Bool_t ass=false;
  Float_t fwhm_tr=1.;
  sigma = (tau_sc*tau_sc+(fwhm_tr/2.3548)*(fwhm_tr/2.3548))/phen;
  sigma =sqrt(sigma);

  return t;
  **/
}

Float_t h10::getsctimebb(Int_t isc,Float_t &sigma)
{
  // This function returns time of PMT signal and its error (sigma)



  Float_t t=0.;
  Float_t fwhm;



  Bool_t raw=false;
  
  if(run > 0)  raw = true;// RAW/MC data switch

  if(isc>=Nsc) return 0.;

  t=-Sc[isc][6]*SC_TDC_CH;

  if(!raw){
      //use resolution from the nemos v. 8 simulations
    if(Sc[isc][2]==0) fwhm=0.171;
    if(Sc[isc][2]==1) {
      if(Sc[isc][3]==0 || Sc[isc][3]==2)   fwhm=0.139;
      if(Sc[isc][3]==1)   fwhm=0.140;
    }
    if(Sc[isc][2]>1){
      if(Sc[isc][3]==0) fwhm=0.158;
      if(Sc[isc][3]==1) fwhm=0.158;
      if(Sc[isc][3]==2) fwhm=0.169;
      if(Sc[isc][3]==3) fwhm=0.148;
    }
  }
  //if(raw){
      // use real resolution measured during calibration
      //fwhm = pmstatus->get_pm_resol(run,(Int_t)Sc[isc][1],(Int_t)Sc[isc][2],(Int_t)Sc[isc][3],(Int_t)Sc[isc][4]);
  //}
		

  
  // Use Vera's value for experimental electron tau 
  Float_t tau_sc;
  //For MC with nemos v 7 tau=6 ns
  if(nemosv==7) tau_sc=6.;
  // In new nemos tau = 4.3 ns
  if(nemosv==8) tau_sc=4.3;
  // If tau_sc saved in the root file (new reconstructions)
  // use it instead 
  if(run<0 && tau_sc_save>1 && tau_sc_save<10) tau_sc=tau_sc_save;
  

  Bool_t ass=false;
  if(raw){ //RAW data
    for(Int_t i=0;i<Nbr_tks;i++){
      ass=ass||(isc==Ind_scintil[i]-1);
    }

    if(ass){
      tau_sc=4.3;//electron track, use tau=4.3ns
    }else{
      tau_sc=4.1;
    }
  }
  //number of photoelectrons
  Float_t phen=(Sc[isc][8]*1000.)*(2.3548/fwhm)*(2.3548/fwhm);



  //PM gain variation for the 
  Float_t VG = 0.0025;

  //-- sigma from transit time fluctuation and TDC convertion    
  Float_t fwhm_tr=1.;
  sigma = (tau_sc*tau_sc+(fwhm_tr/2.3548)*(fwhm_tr/2.3548))/phen/(1 + VG);

  // add additional component due to time calibration for the runs 
  // in the year 2007 ~ 50ps
  if(run > 5469){
    //if(run > 7800){
    sigma = sqrt(sigma*sigma + 0.05*0.05);
  }

  //try to correct TOF distribution for 2MeV events:
  //if(raw){
    //get TDC fluctuation from LTC measurement in ns
    //Float_t tdc_fluct = pmstatus->get_pm_ltce(run,(Int_t)Sc[isc][1],(Int_t)Sc[isc][2],(Int_t)Sc[isc][3],(Int_t)Sc[isc][4]); 
    //sigma = (sigma + tdc_fluct*tdc_fluct);
  //}


  //add additional systematic error
  if(raw){
    sigma = sigma+0.13*0.13;
  }

  sigma =sqrt(sigma);
  return t;
}

Float_t h10::gethypoth(Int_t track, Int_t object, Int_t flag, Int_t ieflag, Float_t* tf)
{  
   Float_t c = 3e8;
   Float_t me = 0.511;
   Float_t i, scdec=4.3;


// What we are trying to do is to see whether two objects
// are consistent with coming from the same place on the foil (internal)
// or whether
// they are consistent with coming from outside, and unrelated to each 
// other or 
// external AND related to each other.

// First, what are the objects? First object must be a track
// Second object is a track if flag=0, and is a gamma if flag=1

// What hypthesis are we trying to check? ieflag=0, internal
// ieflag=1, external.

  //first lets deal with the electron and do the internal hypothesis
  //X_Scintil is projected position in the scintillator
  //[0] is the TRACK index, not the Scint index
  //xfoil is from tracking,the start position of the hypothesis

   Float_t dstrack, sigt1, sigl1, fwhm, ee, beta, dbeta, T1, T2;
   Float_t gamma, tfoil1, tfoil2;
   Int_t track2;
   T1=-Sc[Ind_scintil[track]-1][6]*0.053; //time of e scint hit

   dstrack=gethelixl(track,sigl1);
   Float_t time1=timeofflight(sigt1, tfoil1, track, Ind_scintil[track]-1,0);
   //   std::cout<<" sigt1 "<<sigt1<<" sigl1 "<<sigl1<<std::endl;
   Float_t sectore = Sc[Ind_scintil[track]-1][1];
   Float_t time2, sigt2, sigl2, T, pvec[3];

   //OK, so now lets compare the track with something:
   Float_t delsgam;
   if (flag==1) //a gamma
     { if (object!=Ind_scintil[track]-1)
         { Float_t dsgamma = getstraightlength(object,X_foil[track],Y_foil[track],Z_foil[track],sigl2,pvec); //straight from foil
           T2 = -(Sc[object][6])*0.053;
           if (Sc[object][2]==0)
             { fwhm = 0.17;
             }
           else
             { fwhm = 0.14;
             }
           Float_t npe=Sc[object][8]*1000.*(2.3548*2.3548)/(fwhm*fwhm);
           sigt2 = sqrt(((scdec*scdec)+1./(2.3548*2.3548))/npe);
           time2 = dsgamma/c*1e7; //straight flight, no curvey
           tfoil2 = T2-time2;
           Float_t sectorg = Sc[object][1];
         } 		 
     }
   if (flag==0) //the second object is a track 
     { 
       track2=object;
       Int_t scinthit = Ind_scintil[track2]-1;
       dstrack=gethelixl(track2,sigl2);
       time2=timeofflight(sigt2, tfoil2, track2, scinthit, 0);
       //       std::cout<<" a track sigt2 "<<sigt2<<" sigl2 "<<sigl2<<std::endl;
       T2=-Sc[scinthit][6]*0.053; //time of e scint hit
     }
   Float_t sigt=(sigt1*sigt1+sigt2*sigt2+sigl1*sigl1+sigl2*sigl2);


//         Electron travel time (curvey) is timens.
//         Gamma travel time is time2
//         Second track travel time is time2
//         Time of electron scint hit is T0
//         Time of gamma scint hit is T

// This is the internal hypothesis
   if (ieflag==0)
     {  //Float_t chi=(time2-timens)-(T-T0);
        Float_t chi = tfoil1-tfoil2;
        sigt+=(0.0015*0.0015);
        Float_t chisq = chi*chi/sigt;

	//	std::cout<<" internal chisq "<<chisq<<" sigt "<<sigt<<std::endl;
        Float_t prob = TMath::Prob(chisq,1); 

	//get the weighted mean of the time at the foil 
	//        Float_t tfoil1 = T0 - timens;
        //Float_t tfoil2 = T - time2;
        Float_t dfoil1 = sigt1*sigt1+sigl1*sigl1;
        Float_t dfoil2 = sigt2*sigt2+sigl2*sigl2;
	//	std::cout<<" dfoil1 "<<dfoil1<<" dfoil2 "<<dfoil2<<std::endl;
        tf[1] = 1./sqrt(1./dfoil1+1./dfoil2);
        tf[0] = (tfoil1/dfoil1+tfoil2/dfoil2)/(1./dfoil1+1./dfoil2);
	//	std::cout<<" tf[1] "<<tf[1]<<" tf[0] "<<tf[0]<<std::endl;

        return prob;
     }
// This is the external hypothesis
   if (ieflag==1)
     {  tf[0] = 0;
        Float_t chi=(T1-time1-time2-T2); //ABS
        Int_t scinthit = Ind_scintil[track]-1;
        Float_t dstrack=gethelixl(track2,sigl2);
        time2=timeofflight(sigt2, tfoil2, track2, scinthit, 0);
        dstrack=gethelixl(track,sigl1);
        time1=timeofflight(sigt1, tfoil1, track, scinthit, 0);

        Float_t sigt=(sigt1*sigt1+sigt2*sigt2+sigl1*sigl1+sigl2*sigl2);
        sigt+=(0.0015*0.0015);
        Float_t chisq = chi*chi/sigt;
        Float_t prob = TMath::Prob(chisq,1); 
	//	std::cout<<" me1 "<<T2<<" "<<T1<<" "<<time2<<" "<<time1<<" energy "<<Sc[scinthit][8]*1000.<<std::endl;
        return prob;
     }
   //other external hypothesis
   if (ieflag==2)
     {  tf[0] = 0;
        Int_t scinthit = Ind_scintil[track2]-1;        
        Float_t chi=(T2-time1-time2-T1); //ABS
        Float_t dstrack=gethelixl(track2,sigl2);
        time2=timeofflight(sigt2, tfoil2, track2,scinthit, 0);
        dstrack=gethelixl(track,sigl1);
        time1=timeofflight(sigt1, tfoil1, track,scinthit, 0);

        Float_t sigt=(sigt1*sigt1+sigt2*sigt2+sigl1*sigl1+sigl2*sigl2);
        sigt+=(0.0015*0.0015);
        Float_t chisq = chi*chi/sigt;
        Float_t prob = TMath::Prob(chisq,1); 
	//	std::cout<<" me1 "<<T2<<" "<<T1<<" "<<time2<<" "<<time1<<" energy "<<Sc[scinthit][8]*1000.<<std::endl;
        return prob;
     }
      
}   

Float_t h10::getscatterhyp(Int_t scint1, Int_t scint2, Float_t& chi)
{  
   Float_t sigt, chisq, chi1, chi2;
   Float_t sigl3, sigt3, T1, T2, tof, prob1, prob2;

   // ieflag = internal hypothesis
   // ieflag = external hypothesis

   dgammatogamma(scint1, scint2, tof, sigl3, sigt3); 
   //   std::cout<<" siglen "<<sigl3<<" sigtime "<<sigt3<<std::endl;
   T1 = -(Sc[scint1][6])*0.053;
   T2 = -(Sc[scint2][6])*0.053;
   sigt=(sigt3*sigt3+sigl3*sigl3);


   chi=(T2-T1-tof);
   chisq = chi*chi/sigt;
   chi1 = chi/sqrt(sigt);
   prob1 = TMath::Prob(chisq,1);

   chi=(T1-T2-tof);
   chisq = chi*chi/sigt;
   chi2 = chi/sqrt(sigt);
   prob2 = TMath::Prob(chisq,1);

   if(prob1<prob2)
     { chi = chi2;
       return prob2;
     }
   else
     { chi = chi1;
       return prob1;
     }
}   



Float_t h10::getsinglehyp(Float_t tfoil, Float_t dt0, Float_t vertx, Float_t verty, Float_t vertz, Float_t sigl1, Int_t object, Int_t flag, Int_t ieflag, Float_t& chi,Float_t energy)
{  

   Float_t c = 3e8;
   Float_t me = 0.511;
   Float_t scdec = 4.3;
   Float_t i, fwhm;

   // inputs to this module are time of first object at foil, error on time,
   // time of second object, length of second objects travel, objects index
   // flag about what the second object is, and internal of external hyp flag
   // Corrected, clustered, energy of both objects, e1=gamma, e2=original track

   // flag = 0 object is track index
   // flag = 1 object is gamma index

   // ieflag = internal hypothesis
   // ieflag = external hypothesis

   Float_t time2, sigt2, sigl2, T2, tfoil2, pvec[3];
   if (flag==0) //object is a track
     { T2 = -(Sc[object][6])*0.053;
       Float_t dsegment = gethelixl(object,sigl2);
       time2 = timeofflight(sigt2, tfoil2, object, Ind_scintil[object]-1,0);
       Float_t sectore = Sc[Ind_scintil[object]-1][1];
       tfoil2 = T2-time2;
     }
   if (flag==1) //object is a gamma
     { T2 = -(Sc[object][6])*0.053;
     Float_t dsgamma = getstraightlength(object,vertx,verty,vertz,sigl2,pvec);
       if (Sc[object][2]==0)
         { fwhm = 0.17;
          }
       else
          { fwhm = 0.14;

          }
       Float_t npe=energy*(2.3548*2.3548)/(fwhm*fwhm);
       sigt2 = sqrt(((scdec*scdec)+1./(2.3548*2.3548))/npe);
       time2 = dsgamma/c*1e7; //straight flight, no curvey
       tfoil2=T2-time2;
       Float_t sectorg = Sc[object][1];
 }

// This is the internal hypothesis
   if (ieflag==0)
     {  chi=(tfoil2)-(tfoil);
        Float_t sigt=(dt0*dt0+sigt2*sigt2+sigl2*sigl2+sigl1*sigl1);
        sigt+=(0.0015*0.0015);
	//       	std::cout<<" chi "<<chi<<" sigt1 "<<dt0<<" sigt2 "<<sigt2<<" sigl1 "<<sigl1<<" sigl2 "<<sigl2<<" sigt "<<sigt<<std::endl;
        Float_t chisq = chi*chi/sigt;
        chi = chi/sqrt(sigt);
        Float_t prob = TMath::Prob(chisq,1);
        return prob;
     }
// This is the external hypothesis
   if (ieflag==1)
     {  chi=tfoil-tfoil2-time2;
        Float_t sigt=(dt0*dt0+sigt2*sigt2+sigl2*sigl2+sigl1*sigl1);
        sigt+=(0.0015*0.0015);
        Float_t chisq = chi*chi/sigt;
        Float_t prob = TMath::Prob(chisq,1); 
	//	if(prob<0.1)std::cout<<" ext "<<" timens "<<timens<<" T "<<T<<" t0 "<<t0<<" sigt "<<sigt<<" chisq "<<chisq<<std::endl; 
        chi = chi/sqrt(sigt);
        return prob;
     }
}   

Float_t h10::gethelixl(Int_t tr, Float_t &dlen)
{


      Float_t c = 3e8;
      Float_t me = 0.511;
      Float_t fwhm;
      Float_t beta, dbeta;
      Float_t vsc[2];
      Float_t vvr[2];
      vsc[0]=X_scintil[tr]-Xc[tr];
      vsc[1]=Y_scintil[tr]-Yc[tr];
	vvr[0]=X_foil[tr]-Xc[tr];
	vvr[1]=Y_foil[tr]-Yc[tr];
	Float_t rvsc=sqrt(vsc[0]*vsc[0]+vsc[1]*vsc[1]);
	Float_t rvvr=sqrt(vvr[0]*vvr[0]+vvr[1]*vvr[1]);
	Float_t prod=vsc[0]*vvr[0]+vsc[1]*vvr[1];
	Float_t	cos_check=prod/rvsc/rvvr;
	//Float_t	cos_check=prod/pow(Radc[tr],2);
	if(cos_check>1.){
	  //	  std::cout<<"TRLEN WARNING cos_check="<<cos_check<<std::endl;
	  cos_check=1.	;	   
	} 
	if(cos_check<-1.){
	  //	  std::cout<<"TRLEN WARNING cos_check="<<cos_check<<std::endl;
	  cos_check=-1.	;	   
	}
	Float_t dtetha = acos(cos_check);
	Float_t dstrack=dtetha*sqrt(Radc[tr]*Radc[tr]+Hc[tr]*Hc[tr]);

        Float_t ee=Sc[Ind_scintil[tr]-1][8]*1000.0; //energy of scint
        ee = ee + me; //total energy
        beta = sqrt(1.-1./pow(ee/me,2));

        if(Sc[Ind_scintil[tr]-1][2]==0)
          { fwhm = 0.17;
          }
        else
          { fwhm = 0.14;
          }
        dbeta = 1/beta*me*me/ee*fwhm/2.;
        dlen = dstrack*dbeta/(beta*beta*c)*1e7;// in ns
	return dstrack;
}
Float_t h10::getcurveylength(Int_t trackindex, Int_t scintindex, Float_t& dlen)
  {
      Float_t c = 3e8;
      Float_t me = 0.511;
      Float_t fwhm;
      Float_t beta, dbeta;

      Float_t dsegment=sqrt(pow(X_scintil[trackindex]-X_foil[trackindex],2)+pow(Y_scintil[trackindex]-Y_foil[trackindex],2));

      Float_t term1 = dsegment/Radc[trackindex]/2.;
      Float_t term2 = asin(term1);
      Float_t dstrack = sqrt(pow(2.*term2*Radc[trackindex],2)+pow(Z_scintil[trackindex]-Z_foil[trackindex],2));   //curvy length
      //      std::cout<<" Radc[trackindex] "<<Radc[trackindex]<<" trackindex "<<trackindex<<std::endl;
       Float_t ee=Sc[scintindex][8]*1000.0; //energy of scint

       //       ee = sqrt(pow(ee,2)+pow(me,2));
       ee = ee + me; //total energy
       beta = sqrt(1.-1./pow(ee/me,2));

       if(Sc[scintindex][2]==0)
         { fwhm = 0.17;
         }
       else
         { fwhm = 0.14;
         }
       dbeta = 1/beta*me*me/ee*fwhm/2.;
       //       std::cout<<" beta"<<beta<<" fwhm "<<fwhm<<" ee "<<ee<<" dstrack "<<dstrack<<std::endl;
       dlen = dstrack*dbeta/(beta*beta*c)*1e7;// in ns
       //      std::cout<<" getcurveylen: debta "<<dbeta<<" dlen "<<dlen<<std::cout;
       return dstrack;
  }



Float_t h10::timeofflight(Float_t& sigmat, Float_t& tfoil, Int_t object, Int_t scintobject, Int_t flag)
{  Float_t c = 3e8;
   Float_t me = 0.511;
   Float_t sigl1;

   //only works for tracks: you dont know the starting position for the gammas 
   if (flag==0) //track
      { Float_t ee=Sc[scintobject][8]*1000.0; //energy of scint hit
      Float_t sT0;
        Float_t T0=getsctime(scintobject,sT0); //time of e scint hit
        ee = ee + me;
        Float_t gamma = ee/me;
        Float_t beta = sqrt(1.-1./pow(gamma,2));

        Float_t fwhm;
        if (Sc[scintobject][2]==1)
          { fwhm = 0.14;
          }
        else
          { fwhm = 0.17;
          }

        Float_t dbeta = 1/beta*me*me/ee*fwhm/2; // velocity uncertainty
        Float_t dstrack = gethelixl(object, sigl1); //curvey length
        Float_t timens = dstrack/(c*beta)*1e7; //curvy track going at beta, time
        sigmat = sqrt(sT0*sT0+sigl1*sigl1+timens*timens*dbeta*dbeta);
        tfoil = T0 - timens; //time at foil
        return timens;
      }
}
Float_t h10::getfoiltime(Float_t& sigmat, Float_t& sigl1, Int_t object)
{  
  // This function returns time, when an electron[object] was emitted from the foil. 

   Float_t c = 3e8;
   Float_t me = 0.511;
   Float_t ee=Sc[Ind_scintil[object]-1][8]*1000.0; //energy of scint hit
   Float_t T0=-Sc[Ind_scintil[object]-1][6]*0.053; //time of e scint hit
   ee = ee + me;
   Float_t gamma = ee/me;
   Float_t beta = sqrt(1.-1./pow(gamma,2));

   Float_t fwhm;
    if (Sc[Ind_scintil[object]-1][2]==0)
       { fwhm = 0.17;        
       }
    else
       { fwhm = 0.14;
       }

    //    Float_t dbeta = 1/beta*me*me/ee*fwhm/2;
    Float_t npe=Sc[Ind_scintil[object]-1][8]*1000.*(2.3548*2.3548)/(fwhm*fwhm);
    Float_t dstrack = getcurveylength(object, Ind_scintil[object]-1,sigl1); //curvey length

    sigmat = sqrt(((6.0*6.0)+1./(2.3548*2.3548))/npe);
    Float_t timens = dstrack/(c*beta)*1e7; //curvy track going at beta, time
    Float_t tfoil = T0 - timens; //time at foil
    return tfoil;
 
}

Float_t h10::getstraightlength(Int_t scint,Float_t vertx, Float_t verty, Float_t vertz, Float_t& sigl, Float_t* pvec)
{  Float_t c = 3e8;
   Float_t dsgamma = sqrt(pow(Sc[scint][9]-vertx,2)+pow(Sc[scint][10]-verty,2)+pow(Sc[scint][11]-vertz,2));
   Float_t xg=Sc[scint][9];
   Float_t yg=Sc[scint][10];
   Float_t zg=Sc[scint][11];

   pvec[0] = (Sc[scint][9]-vertx)/dsgamma;
   pvec[1] = (Sc[scint][10]-verty)/dsgamma;
   pvec[2] = (Sc[scint][11]-vertz)/dsgamma;

//    std::cout<<" xg "<<xg<<" yg "<<yg<<" zg "<<zg<<std::endl;
//    std::cout<<" xv "<<vertx<<" yv "<<verty<<" zv "<<vertz<<std::endl;
//    std::cout<<" pgx "<<pvec[0]<<" pgy "<<pvec[1]<<" pgz "<<pvec[2]<<std::endl;

   dsgamma = dsgamma+19.0; //add on depth into scintillator on length
   Float_t delsgam = sqrt(15.0*15.0+dsgamma*dsgamma)-dsgamma;
   delsgam = sqrt(delsgam*delsgam+19.0*19.0); //add on depth into scintillator on error
   sigl = delsgam/c*1e7; //in ns
   return dsgamma;
}


void h10::dgammatogamma(Int_t scint1, Int_t scint2, Float_t& tof, Float_t& dl, Float_t& dt)
{  Float_t c = 3e8;
   Float_t ds = 0.0; //extra #cm into scint block the gamma travelled on av
   Float_t fwhm1, fwhm2;
   //correct for approximate position inside scintillator
   Float_t theta1 = TMath::ATan(Sc[scint1][9]/Sc[scint1][10]);
   Float_t theta2 = TMath::ATan(Sc[scint2][9]/Sc[scint2][10]);
   Float_t yprime1 = ds*TMath::Sin(theta1);
   Float_t yprime2 = ds*TMath::Sin(theta2);
   Float_t xprime1 = ds*TMath::Cos(theta1);
   Float_t xprime2 = ds*TMath::Cos(theta2);
   Float_t zprime1 = ds*TMath::Tan(theta1);
   Float_t zprime2 = ds*TMath::Tan(theta2);


   Float_t x1 = Sc[scint1][9];
   Float_t y1 = Sc[scint1][10];
   Float_t z1 = Sc[scint1][11];
   Float_t x2 = Sc[scint2][9];
   Float_t y2 = Sc[scint2][10];
   Float_t z2 = Sc[scint2][11];


   x1 = x1 + xprime1;
   x2 = x2 + xprime2;
   y1 = y1 + yprime1;
   y2 = y2 + yprime2;
   z1 = z1 + zprime1;
   z2 = z2 + zprime2;

   Float_t dsgamma = sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2));
   Float_t delsgam = 10.0;
   tof = dsgamma/c*1e7; //in ns
   dl = delsgam/c*1e7; //in ns

   Float_t e1=Sc[scint1][8]*1000.0; //energy of scint hit
   Float_t T1=-Sc[scint1][6]*0.053; //time of e scint hit
   Float_t e2=Sc[scint2][8]*1000.0; //energy of scint hit
   Float_t T2=-Sc[scint2][6]*0.053; //time of e scint hit

   Float_t fwhm;
    if (Sc[scint1][2]==0)
       { fwhm1 = 0.17;        
       }
    else
       { fwhm1 = 0.14;
       }

    if (Sc[scint2][2]==0)
       { fwhm2 = 0.17;        
       }
    else
       { fwhm2 = 0.14;
       }


    Float_t npe=Sc[scint1][8]*1000.*(2.3548*2.3548)/(fwhm1*fwhm1);
    Float_t sigmat1 = sqrt(((6.0*6.0)+1./(2.3548*2.3548))/npe);
    npe=Sc[scint2][8]*1000.*(2.3548*2.3548)/(fwhm2*fwhm2);
    Float_t sigmat2 = sqrt(((6.0*6.0)+1./(2.3548*2.3548))/npe);
    dt = sqrt(sigmat1*sigmat1+sigmat2*sigmat2);

    // return dsgamma;
}
Bool_t h10::checkfirstlayer(Int_t trackindex, Int_t iflag)
{    Bool_t foil=false;
     Float_t smin=100.;
     for (Int_t j=0;j<Nbr_pts[trackindex];j++)
        { 
          Int_t ggindex=Ind_points[trackindex][j]-1;
          Int_t layer=int(Gg[ggindex][3]+0.5);
          Int_t cellno=int(Gg[ggindex][4]+0.5);
          Int_t io=int(Gg[ggindex][2]+0.5);
          Int_t sec=int(Gg[ggindex][1]+0.5);

          Float_t ggx=Gg[ggindex][9];
          Float_t ggy=Gg[ggindex][10];
          Float_t ggz=Gg[ggindex][11];
          Float_t radius=sqrt(pow(ggx,2)+pow(ggy,2));
          Float_t times=Gg[ggindex][6];
 	  Float_t timeq=Gg[ggindex][5];
          Float_t tfast=(307-timeq)*0.02;
          Float_t tslow=(35430-times)*0.02;
          foil=(foil||layer==0||layer==1);
          if (layer==0||layer==1)
            { Float_t xf = ggx;
 	      Float_t yf = ggy;
              Float_t s = sqrt(pow(X_foil[trackindex]-ggx,2)+pow(Y_foil[trackindex]-ggy,2)+pow(Z_foil[trackindex]-ggz,2));
              if(s<smin)smin=s;
             }
        }
      if(iflag=1)foil = foil&&(smin<20);

      return foil;
}
Bool_t h10::checklayer(Int_t trackindex, Int_t row)
{   
  // This subroutine checks if track[trackindex] has a hit in specified layer
  // Author V.Vasiliev

     Bool_t foil=false;
     for (Int_t j=0;j<Nbr_pts[trackindex];j++)
        { 
          Int_t ggindex=Ind_points[trackindex][j]-1;
          Int_t layer=int(Gg[ggindex][3]);
          foil=(foil||layer==row);
        }
      return foil;
}
Int_t h10::getsector(Int_t trackindex)
{    Bool_t foil=false;
     Int_t sec;
     for (Int_t j=0;j<Nbr_pts[trackindex];j++)
        { Int_t ggindex=Ind_points[trackindex][j]-1;
          Int_t layer=int(Gg[ggindex][3]+0.5);
          if (layer==0||layer==1)
            { sec=int(Gg[ggindex][1]+0.5);
             }
        }
      return sec;
}
Bool_t h10::amIalone(Int_t scintNum1, Int_t scintNum2, Int_t flag)
{ 
  //flag=0, ignore diagonal neighbours
  //flag=1,count diagonal neighours as associated

  Int_t block, layer, area, sector;
  Int_t blockn, layern, arean, sectorn;
  Bool_t alone=true;

  block = int(Sc[scintNum1][4]);
  layer = int(Sc[scintNum1][3]);
  area = int(Sc[scintNum1][2]);
  sector = int(Sc[scintNum1][1]);  

  blockn = int(Sc[scintNum2][4]);
  layern = int(Sc[scintNum2][3]);
  arean = int(Sc[scintNum2][2]);
  sectorn = int(Sc[scintNum2][1]);

  if(area!=arean)
    { alone = true;
      return alone;
    }  
  if(flag==0&&(area==3||area==2))
    { alone = true;
      return alone;
    }
  //give the block a unique code

  Int_t blockx, blocky, blockxplus,blockxminus,blockyplus, blockyminus;
  Int_t newlayer;
  if(area==1) //external
    { blockx = sector*3+layer;
      blockxplus=blockx+1;
      blockxminus=blockx-1;
      blocky = block;
      if(blockx==0)blockxminus = 59;
      if(blockx==59)blockxplus=0;
    }
  if(area==0) //internal
    { newlayer = layer;
      if(newlayer==2)newlayer=1; //get rid of stupid missing layer
      blockx = sector*2+newlayer;
      blocky = block;
      blockxplus=blockx+1;
      blockxminus=blockx-1;
      if(blockx==0)blockxminus = 39;
      if(blockx==39)blockxplus=0;
    }

  if(area==2) //bottom petal
    { blockx = sector*3+block;
      blockxplus=blockx+1;
      blockxminus=blockx-1;
      blocky = layer;
      if(blockx==0)blockxminus = 59;
      if(blockx==59)blockxplus=0;
    }

  if(area==3) //top petal
    { blockx = sector*3+block;
      blockxplus=blockx+1;
      blockxminus=blockx-1;
      blocky = layer;
      if(blockx==0)blockxminus = 59;
      if(blockx==59)blockxplus=0;
    }

  blockyplus=blocky+1;
  blockyminus=blocky-1;

  Int_t nblockx=999, nblocky=999;
  Bool_t bingo = false;

  // these are the second scintillators unique numbers.

  if (area==0&&int(Sc[scintNum2][2])==area)
    { newlayer = layern;
      if(newlayer==2)newlayer=1; //get rid of stupid missing layer
      nblockx = sectorn*2+newlayer;
      nblocky = blockn;
    }
  if(area==1&&int(Sc[scintNum2][2])==area)
    { nblockx = sectorn*3+layern;
      nblocky = blockn;
    }
  if(area==2&&int(Sc[scintNum2][2])==area)
    { nblockx = sectorn*3+blockn;
      nblocky = layern;
    }
  if(area==3&&int(Sc[scintNum2][2])==area)
    { nblockx = sectorn*3+blockn;
      nblocky = layern;
    }

  bingo=bingo||(nblocky==blocky&&(nblockx==blockxplus));
  bingo=bingo||(nblocky==blocky&&(nblockx==blockxminus));
  bingo=bingo||(nblockx==blockx&&(nblocky==blockyplus));
  bingo=bingo||(nblockx==blockx&&(nblocky==blockyminus));
  bingo=bingo||((nblocky==blockyplus)&&(nblockx==blockxplus));
  bingo=bingo||((nblocky==blockyplus)&&(nblockx==blockxminus));
  bingo=bingo||((nblocky==blockyminus)&&(nblockx==blockxplus));
  bingo=bingo||((nblocky==blockyminus)&&(nblockx==blockxminus));

  //      std::cout<<" nblocky "<<nblocky<<" blocky "<<blocky<<" nblockx "<<nblockx<<" blockxplus "<<blockxplus<<" blockxminus "<<blockxminus<<std::endl;

  //    std::cout<<" nblockx "<<nblockx<<" blockx "<<blocky<<" nblocky "<<nblocky<<" blockyplus "<<blockyplus<<" blockyminus "<<blockyminus<<std::endl;

       if(bingo&&(area==arean))
         { 
           alone=false;
         }
       if(!bingo)
         { alone=true;
         }

       return alone;     
      			     

 }


void h10::Init_scintillator(scintillator* sci, Bool_t* Isocluster, Int_t iflag)
{
  //Looks for scintillator clusters
  //iflag = 0,1 group by side/side+corner 
  //iflag=2 do not group scintillators in cluster
        Int_t jflag;
        Bool_t used[250];
        for (Int_t i=0; i<250; i++){ used[i]=false; }
        Int_t Icluster[250]={0}, nisol=0;
	sci->icluster=0;
        for (Int_t i=0; i<Nsc; i++){  
	  if(used[i])continue;

	  sci->ncluster[sci->icluster].Esum=0;
	  Isocluster[i] = true;
	  sci->ncluster[sci->icluster].index[0]=i;
	  sci->ncluster[sci->icluster].ngroup=1;
	  
	  if(iflag<2)jflag=iflag;
	  if(iflag==2)jflag=0;  

	  for (Int_t j=0;j<Nsc; j++){ 
	    if (i!=j){ 
	      Bool_t solo = amIalone(i,j,jflag);
	      Isocluster[i] = Isocluster[i] && solo;
	      if(!solo&&iflag<2&&!used[j]){
		Int_t n=sci->ncluster[sci->icluster].ngroup;
		sci->ncluster[sci->icluster].index[n]=j;
		sci->ncluster[sci->icluster].ngroup++;
		used[j] = true;
	      }
	    }
	  }
	  
	  //Check if cluster is assigned to a track
	  sci->ncluster[sci->icluster].assigned=false;
	  for (Int_t k=0; k<Nbr_tks; k++)
	    for(Int_t j=0;j<sci->ncluster[sci->icluster].ngroup;j++){
	      if(Ind_scintil[k]-1==sci->ncluster[sci->icluster].index[j])sci->ncluster[sci->icluster].assigned=true;
	    }  
	  
	  sci->ncluster[sci->icluster].Iseed = i;
	  sci->ncluster[sci->icluster].Esum = 0;
	  Float_t Eseed = 0;
	  for (Int_t k=0; k<sci->ncluster[sci->icluster].ngroup; k++){ 
	    Float_t Es = Sc[sci->ncluster[sci->icluster].index[k]][8]*1000.;
	    sci->ncluster[sci->icluster].Esum+=Es;
	    if (Es > Eseed){
	      sci->ncluster[sci->icluster].Iseed=sci->ncluster[sci->icluster].index[k];
	      Eseed = Es;
	    }
	  }
	  sci->icluster++;
	  used[i]=true;
	}
	/*	std::cout<<"Clusters "<<sci->icluster<<std::endl;
	 *	for(Int_t i=0;i<sci->icluster;i++){
	 * std::cout<<" ass["<<i<<"]= "<<sci->ncluster[i].assigned;
	 * std::cout<<" CL["<<i<<"] esum="<<sci->ncluster[i].Esum<<" seed "<<sci->ncluster[i].Iseed<<" ## ";
	 * for(Int_t k=0;k<sci->ncluster[i].ngroup;k++){
	 *   std::cout<<sci->ncluster[i].index[k]<<"("<<Sc[sci->ncluster[i].index[k]][8]*1000<<") ";
	 * }
	 *	std::cout<<std::endl;
	 *}
	 std::cout<<Ind_scintil[0]-1<<std::endl; */
}
void h10::addclusters(scintillator *sci, Int_t cluster1, Int_t cluster2)
{       

        std::cout<<" addcluster: sci.icluster "<<sci->icluster<<std::endl;
        Int_t idominant, ideferential;
        sci->icluster = sci->icluster-1;
	//reduce number of clusters by one
        if(cluster1<cluster2){idominant=cluster1;ideferential=cluster2;}
        if(cluster2<cluster1){idominant=cluster2;ideferential=cluster1;}

// 	if(cluster1<cluster2)std::cout<<" 1 < 2 "<<std::endl;
// 	if(cluster2<cluster1)std::cout<<" 2 < 1 "<<std::endl;
 
	//make lowest index the dominant one: add the energies
        sci->ncluster[idominant].Esum+=sci->ncluster[ideferential].Esum;
        for(Int_t i=0; i<=sci->ncluster[ideferential].ngroup; i++)
	   { sci->ncluster[idominant].ngroup++;
	     sci->ncluster[idominant].index[sci->ncluster[idominant].ngroup]=sci->ncluster[ideferential].index[i];
	   }
        sci->ncluster[idominant].Iseed = 0;
        Float_t Eseed = 0;
        sci->ncluster[idominant].Esum=0;
	for (Int_t k=0; k<=sci->ncluster[idominant].ngroup; k++)
           { Float_t Es = Sc[sci->ncluster[idominant].index[k]][8]*1000.;
	     if (Es > Eseed)
               { sci->ncluster[idominant].Esum+=Es;
                 sci->ncluster[idominant].Iseed=sci->ncluster[idominant].index[k];  
                 Eseed = Es;
               }

	  }
        std::cout<<" addcluster: sci.icluster new "<<sci->icluster<<std::endl;
	for (Int_t k=0; k<=sci->icluster; k++)
           { Float_t Es = sci->ncluster[k].Esum;
      	     std::cout<<" k "<<k<<" Esum "<<Es<<" Iseed "<<sci->ncluster[k].Iseed<<std::endl;
      	     std::cout<<" ngroup "<<sci->ncluster[k].ngroup<<std::endl;
           }
        for (Int_t k=ideferential; k<sci->icluster; k++)
	  { sci->ncluster[k].Esum = sci->ncluster[k+1].Esum;
	    sci->ncluster[k].Iseed = sci->ncluster[k+1].Iseed;
	    sci->ncluster[k].ngroup = sci->ncluster[k+1].ngroup;
            for (Int_t j=0; j< sci->ncluster[k].ngroup; j++)
	      { sci->ncluster[k].index[j] = sci->ncluster[k+1].index[j];
              }
	    std::cout<<" new k " <<k<<" Esum "<<sci->ncluster[k].Esum<<" Iseed "<<sci->ncluster[k].Iseed<<std::endl;
          }
}
Int_t h10::getalphas(Int_t* alparray)
{	//sort out the assigned hits from the spare ones.
        //loop through ALL hits and find the slow ones
        //then loop through assigned hits and see what is left

        int alp=0;
        Int_t nslow=0;
        for (Int_t j=0; j<Ngg; j++)
          {  Bool_t bingo=false;
             Float_t times=Gg[j][6];
             Float_t timeq=Gg[j][5];

             if (times>0)
               { Float_t tslow = Gg[j][12]/1000.;
                 if (tslow >55&& tslow< 700)
 	           { Int_t georgie = int(times)-int(times/4096)*4096;
 	             if (georgie==int(timeq))
                       { nslow++;
                         alp++;
	                 alparray[alp]=j;
                       }
                   }
 	       }
           }
	return alp;
}

Float_t h10::find_cos_eg(Int_t gcluster, Int_t track, scintillator  sci){
  
  // This function finds  cos of angle between electron and photon
	Bool_t Isocluster[250];

  Int_t sciobj=sci.ncluster[gcluster].Iseed;

  Float_t dif[3];
  dif[0]=-X_foil[track]+Sc[sciobj][9];

  dif[1]=-Y_foil[track]+Sc[sciobj][10];
  dif[2]=-Z_foil[track]+Sc[sciobj][11];
  //normalize the vector
  Float_t rdif=sqrt(dif[0]*dif[0]+dif[1]*dif[1]+dif[2]*dif[2]);
  
  Float_t prod=Cos_dir[track][0]*dif[0]+Cos_dir[track][1]*dif[1]+Cos_dir[track][2]*dif[2];
  Float_t coseg=prod/rdif;
  
  return coseg;
}


Bool_t CheckRunStatus(Int_t rs,Int_t s){
  //Checks if runstatus rs contains status s returns true or false
  
  if(s<10) return rs==s;
  while(s>=10){
    s=s/10;
    rs=rs/10;
    if(s<10){
      return (s==rs%10);
    }
  }
  return false;
}

