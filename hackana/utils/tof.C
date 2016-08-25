#include "h10.h"
#include <fstream>
#include <iostream>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>


Int_t h10::twogamma(Float_t X, Float_t Y,Float_t Z,Float_t tfoil, Float_t dfoil, Float_t &Eg1,Float_t &Eg2){
  //This function checks if there are two gamma real gamma inside the event
  //with 2 tracks. Supposed that gammas are emitted from (X,Y,Z) at tfoil+-dfoil
  //returns 1 if success, Eg1 and Eg2 -- total gammas energy MeV

  	  //Extended gamma cuts
  Float_t Egth=0.150;//Threshold for main gammas scintillators. Should be high enough for reliable time calculations

	  Float_t p2gintmax=0.;
	  Int_t isc1,isc2;
	  Bool_t twogamma=false;
	  //Check that event topology is correct
	  if(Nbr_tks!=2) return 0;
	  if(Nsc<4) return 0;
	  

	  //Look for internal gammas above threshold Egth
	  for(Int_t ii=0;ii<Nsc;ii++){
	    if(ii==Ind_scintil[0]-1 || ii==Ind_scintil[1]-1) continue;
	    if(Sc[ii][8]*1E3<Egth) continue;
	    for(Int_t jj=ii+1;jj<Nsc;jj++){
	      if(jj==Ind_scintil[0]-1 || jj==Ind_scintil[1]-1) continue;
	      if(Sc[jj][8]*1E3<Egth) continue;
	      Float_t p2gint=TOF2g_int(ii,jj,X,Y,Z,tfoil,dfoil);
	      if(p2gint>p2gintmax){
		p2gintmax=p2gint;
		isc1=ii;
		isc2=jj;
	      }
	    }
	  }
	  if(p2gintmax==0) return 0;
	  Eg1=Sc[isc1][8]*1000;
	  Eg2=Sc[isc2][8]*1000;
	  twogamma=(p2gintmax>0.05);
	  //	  std::cout<<"Internal gamma prob: "<<p2gintmax<<" ISC: "<<isc1<<" "<<isc2<<std::endl;
	  //Look if they are not one rescattered
	  Float_t p2gscat1=TOF2g_scat(isc1,isc2,X,Y,Z,tfoil,dfoil);
	  Float_t p2gscat2=TOF2g_scat(isc2,isc1,X,Y,Z,tfoil,dfoil);
	  Float_t p2gscat;
	  if(p2gscat1<p2gscat2) p2gscat=p2gscat2;
	  if(p2gscat1>p2gscat2) p2gscat=p2gscat1;
	  //	  std::cout<<"Scatter gamma prob: "<<p2gscat<<std::endl;
	  //	  twogamma=twogamma&&(p2gscat>p2gintmax);
	  //	  if(twogamma) std::cout<<"BINGO"<<std::endl;
	  
	  //Look if other gammas are rescattered from these two 
	  for(Int_t ii=0;ii<Nsc;ii++){
	    if(ii==Ind_scintil[0]-1 || ii==Ind_scintil[1]-1) continue;
	    if(ii==isc1 || ii==isc2) continue;
	    Float_t pgscat1=TOF2g_scat(isc1,ii,X,Y,Z,tfoil,dfoil);
	    Float_t pgscat2=TOF2g_scat(isc2,ii,X,Y,Z,tfoil,dfoil);
	    Float_t pgscat;
	    if(pgscat1<pgscat2) {
	      pgscat=pgscat2;
	      Eg2+=Sc[ii][8]*1000.;
	    }
	    if(pgscat1>pgscat2) {
	      pgscat=pgscat1;
	      Eg1+=Sc[ii][8]*1000.;
	    }
	    twogamma=twogamma&&(pgscat>0.01);
	  }
	  if(twogamma) return 1;
	  return 0;
}
Float_t   h10::TOFbg_ext(Int_t iscg, Float_t &tf, Float_t &sf)
{
  // Calculates TOF prpobability for external hypothesis for electron-gamma events
  // Only topology gamma->electron is considered here
  // Retruns Chi2 prob.
  Float_t prob=1.E-30;
  // Theoretical TOF for internal hypothesis
  Float_t tofb[2];
  Float_t dtofb[2];
  //Measured TOF for internal hypothesis
  Float_t tofbm[2];
  Float_t rme=0.511;
  Float_t vlight=3.E10, dlen;
  Float_t fwhm;
  Float_t sige;
  
  //electron tof, measured time and error
  if(Nbr_tks!=1) return prob;
  for(Int_t j=0;j<Nbr_tks;j++){
    Int_t isc=Ind_scintil[j]-1;
    Float_t eb=Sc[isc][8]*1000.;
    Float_t ve=sqrt(eb*eb+2.*rme*eb)/(rme+eb)*vlight;
    Float_t dstrack=gethelixl(j,dlen);
    tofb[j]=dstrack/ve*1.E9;
    tofbm[j]=getsctime(isc,sige);
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

  //gamma tof, measured time and error
  Float_t X,Y,Z;
  X=X_foil[0];
  Y=Y_foil[0];
  Z=Z_foil[0];
  tofbm[1]=getsctime(iscg,sige);
  tofb[1]=-getstraightlengthv(iscg,X,Y,Z,dlen)/30.;
  dtofb[1]=sqrt(sige*sige+dlen*dlen);
    

  //calculate chi2 and probability
  Float_t dtt=tofbm[1]-tofbm[0]-(tofb[1]-tofb[0]);
  Float_t chiint2=dtt*dtt/(dtofb[1]*dtofb[1]+dtofb[0]*dtofb[0]);
  prob = TMath::Prob(chiint2,1); 
  tf=(tofbm[0]-tofb[0])/(dtofb[0]*dtofb[0])+(tofbm[1]-tofb[1])/(dtofb[1]*dtofb[1]);
  sf=sqrt(1./(dtofb[0]*dtofb[0])+1/(dtofb[1]*dtofb[1]));
  tf=tf/(sf*sf);
  sf=1./sf;
  //std::cout<<"T1-T2 "<<dtt<<" TOFe,g "<<tofb[0]<<","<<tofb[1]<<" Se,g "<<dtofb[0]<<","<<dtofb[1]<<std::endl;
  return prob;  
}

Float_t   h10::TOFbg_int(Int_t iscg, Float_t &tf, Float_t &sf)
{
  //Calculates TOF prpobability for internal hypothesis for electron-gamma events
  // Retruns Chi2 prob.
  Float_t prob=1.E-30;
  // Theoretical TOF for internal hypothesis
  Float_t tofb[2];
  Float_t dtofb[2];
  //Measured TOF for internal hypothesis
  Float_t tofbm[2];
  Float_t rme=0.511;
  Float_t vlight=3.E10, dlen;
  Float_t fwhm;
  Float_t sige;


  //electron tof, measured time and error
  if(Nbr_tks!=1) return prob;
  for(Int_t j=0;j<Nbr_tks;j++){
    Int_t isc=Ind_scintil[j]-1;
    Float_t eb=Sc[isc][8]*1000.;
    Float_t ve=sqrt(eb*eb+2.*rme*eb)/(rme+eb)*vlight;
    Float_t dstrack=gethelixl(j,dlen);
    tofb[j]=dstrack/ve*1.E9;
    tofbm[j]=getsctime(isc,sige);
    dtofb[j]=sige;

    fwhm = getscfwhm_real(isc);

    Float_t se1 = fwhm*sqrt(eb)/2.354; 
    sige = tofb[j]*rme*rme/eb/(eb+rme)/(eb+2*rme)*se1;
    dtofb[j]=sqrt(sige*sige+dtofb[j]*dtofb[j]);
  }

  //gamma tof, measured time and error

  Float_t X,Y,Z;
  X=X_foil[0];
  Y=Y_foil[0];
  Z=Z_foil[0];
  tofbm[1]=getsctime(iscg,sige);
  tofb[1]=(getstraightlengthv(iscg,X,Y,Z,dlen))/30.;
  dtofb[1]=sqrt(sige*sige+dlen*dlen);
    

  //calculate chi2 and probability
  Float_t dtt=tofbm[1]-tofbm[0]-(tofb[1]-tofb[0]);
  Float_t chiint2=dtt*dtt/(dtofb[1]*dtofb[1]+dtofb[0]*dtofb[0]);
  prob = TMath::Prob(chiint2,1); 
  tf=(tofbm[0]-tofb[0])/(dtofb[0]*dtofb[0])+(tofbm[1]-tofb[1])/(dtofb[1]*dtofb[1]);
  sf=sqrt(1./(dtofb[0]*dtofb[0])+1/(dtofb[1]*dtofb[1]));
  tf=tf/(sf*sf);
  sf=1./sf;
  //td::cout<<Ind_scintil[0]-1<<" "<<iscg<<" T1-T2 "<<tofbm[1]-tofbm[0]<<" TOFe-g "<<tofb[0]-tofb[1]<<" Se,g "<<dtofb[0]<<","<<dtofb[1]<<" "<<prob<<std::end;l
  return prob;  
}

//************************************************************************************


Float_t   h10::TOFbg_scat(Int_t iscg)
{

  //*******************************************************************************
  //
  //  Calculates TOF prpobability for hypothesis of scattered gamma quanta:
  //   g -> Scg -> Sce or g -> Sce ->Scg 
  //  
  //  Retruns Chi2 prob.
  //
  //*******************************************************************************
  Float_t prob=1.E-30;

  Float_t rme=0.511;
  Float_t vlight=3.E10, dlen;
  Float_t fwhm;
  Float_t sigma1, sigma2;

  Int_t isce = Ind_scintil[0]-1;

  Float_t tof_cal, tof_meas;
  Float_t dtof_cal, dtof_meas;

  Float_t chi21,chi22;

  // 1st case g -> Sce -> Scg
  // gamma is always quicker then electron, so do not check the winner
  Float_t X,Y,Z;
  X = Sc[isce][9];
  Y = Sc[isce][10];
  Z = Sc[isce][11];

  tof_cal=(getstraightlengthv(iscg,X,Y,Z,dlen))/30.;
  dtof_cal = dlen;

  tof_meas = getsctime(iscg,sigma1)-getsctime(isce,sigma2);
  dtof_meas = sqrt(sigma1*sigma1 + sigma2*sigma2);
  
  chi21 = pow(tof_cal - tof_meas,2)/(dtof_meas*dtof_meas + dtof_cal*dtof_cal);

  // 2nd case g ->Scg ->Sce
  // electron can be first so need to check the winner, but not for now
  X = Sc[iscg][9];
  Y = Sc[iscg][10];
  Z = Sc[iscg][11];

  tof_cal=(getstraightlengthv(isce,X,Y,Z,dlen))/30.;
  dtof_cal = dlen;

  tof_meas = getsctime(isce,sigma1)-getsctime(iscg,sigma2);
  dtof_meas = sqrt(sigma1*sigma1 + sigma2*sigma2);
  
  chi22 = pow(tof_cal - tof_meas,2)/(dtof_meas*dtof_meas + dtof_cal*dtof_cal);

  if(chi21 < chi22){
    prob = TMath::Prob(chi21,1);
  }else{
    prob = TMath::Prob(chi22,1);
  }
  return prob;  
}

//************************************************************************************

Float_t h10::TOF1g_int(Int_t isc1,Float_t X,Float_t Y,Float_t Z,Float_t tfoil,Float_t dtfoil){
  // This function calculated probability that 1 gamma was emmited from 
  //one point (X,Y,Z) at the time tfoil+-dtfoil

  Float_t sig[2];
  Float_t sct[2];
  Float_t glen[2];
  Float_t sglen[2];
  Float_t tsig;
  sct[0]=getsctime(isc1,sig[0]);
  glen[0]=getstraightlengthv(isc1,X,Y,Z,sglen[0])/30.;

  tsig=(sig[0]*sig[0]+sglen[0]*sglen[0]+dtfoil*dtfoil);
  Float_t chi2=pow(sct[0]-(tfoil+glen[0]),2)/tsig;
  Float_t prob=TMath::Prob(chi2,1);
  //  std::cout<<"TOFg1_int output, dt= "<<(sct[0]-(tfoil+glen[0]))<<" prob="<<prob<<std::endl;
  //  std::cout<<" sigsc = "<<sig[0]<<" sigl="<<sglen[0]<<" sigfoil="<<dtfoil<<" sigtot="<<tsig<<std::endl;
  return prob;
}

Float_t h10::TOF1g_ext(Int_t isc1,Float_t X,Float_t Y,Float_t Z,Float_t tfoil,Float_t dtfoil){
  // This function calculated probability that 1 gamma was emmited from 
  //scintillator and made scattering in the point (X,Y,Z) at the time tfoil+-dtfoil

  Float_t sig[2];
  Float_t sct[2];
  Float_t glen[2];
  Float_t sglen[2];
  sct[0]=getsctime(isc1,sig[0]);
  glen[0]=getstraightlengthv(isc1,X,Y,Z,sglen[0])/30.;

  sglen[0]=(sig[0]*sig[0]+sglen[0]*sglen[0]+dtfoil*dtfoil);

  Float_t chi2=pow(sct[0]-(tfoil-glen[0]),2)/sglen[0];
  Float_t prob=TMath::Prob(chi2,1);
  return prob;
}

Int_t h10::twogammav2(Float_t X, Float_t Y,Float_t Z,Float_t tfoil, Float_t dfoil, Float_t &Egamma1,Float_t &Egamma2,Float_t &critery){
  //This function checks if there are two real gamma inside the event
  //with 2 tracks. Supposed that gammas are emitted from (X,Y,Z) at tfoil+-dfoil
  //returns 1 if success, Eg1 and Eg2 -- total gammas energy MeV
  //Another algorithm

	  Float_t pmax=0.;
	  Float_t Eg1,Eg2;
	  Int_t isc1,isc2;
	  Bool_t twogamma=false;
	  Bool_t thresh1,thresh2, hth1,hth2;
	  Float_t Egth=0.00;//Threshold for main gammas scintillators. Should be high enough for reliable time calculations
	  Float_t Egthh=0.550;//High threshold for gammas scintillators.(Compton edge)
	  //Check that event topology is correct
	  critery=1.E30;
	  if(Nbr_tks!=2) return 0;
	  if(Nsc<4) return 0;
	  Int_t nmax=0;
	  Int_t scint[15];
	  for(Int_t i=0;i<Nsc;i++){
	    if(i!=Ind_scintil[0]-1 && i!=Ind_scintil[1]-1){
	      scint[nmax]=i;
	      nmax++;
	    }
	    if(nmax>14) return 0;
	  }

	  Int_t gamma1[15],ng1;
	  Int_t gamma2[15],ng2;
	  Int_t nc=Int_t((pow(2,nmax))-0.5);
	  //	  std::cout<<"Gamma sc "<<nmax<<" "<<nc<<std::endl;
	  for(Int_t i=1;i<nc;i++){
	    // Use binary presentation for looping over all possible 
	    // PMT combination.
	    //eg for nmax=3 we have 001 010 011 100 101 110 
            // Place is PMT number (1,2,3), bit -- to wich gamma it corresponds
	    //(1st or 2nd).
	    // Look for the most probable combination and order.
	    ng1=0;
	    ng2=0;
	    Eg1=0;
	    Eg2=0;
	    thresh1=false;
	    thresh2=false;
	    hth1=true;
	    hth2=true;
	    for(Int_t j=0;j<nmax;j++){
	      if(1&i>>j){
		gamma1[ng1]=scint[j];
		ng1++;
		thresh1=thresh1||(Sc[scint[j]][8]*1E3>Egth);
		hth1=hth1&&(Sc[scint[j]][8]*1E3<Egthh);
		Eg1+=Sc[scint[j]][8]*1E3;
	      }else{
		gamma2[ng2]=scint[j];
		ng2++;
		thresh2=thresh2||(Sc[scint[j]][8]*1E3>Egth);
		hth2=hth2&&(Sc[scint[j]][8]*1E3<Egthh);
		Eg2+=Sc[scint[j]][8]*1E3;
	      }
	    }
	      Egamma1=Eg1;
	      Egamma2=Eg2;
	    if(ng1==0 || ng2==0) continue;
	    //if(!thresh1 || !thresh2 || ! hth1 || !hth2) continue;
	    Float_t chi21=chi2_tof1g(gamma1,ng1,X,Y,Z,tfoil,dfoil);
	    Float_t chi22=chi2_tof1g(gamma2,ng2,X,Y,Z,tfoil,dfoil);
	    Float_t chi2=chi22+chi21;
	    Float_t prob=TMath::Prob(chi2,ng1+ng2);
	    if(chi2<critery){
	      pmax=prob;
	      critery=chi2;
	    }

	  }
	  critery=critery/(ng1+ng2);
	  if(pmax>0.000005) twogamma=true;
	  std::cout<<"Event "<< Myievent<<" prob "<<pmax<<" chi2 "<<critery<<" E: "<<Egamma1<<" "<<Egamma2<<" "<<nmax<<std::endl; 	  
	    if(twogamma) return 1;
	  return 0;
}


Int_t h10::onegammav2(Float_t X, Float_t Y,Float_t Z,Float_t tfoil, Float_t dfoil, Float_t &Egamma1,Float_t &critery){
  //This function checks if there is a real gamma inside the event
  //with 2 tracks. Supposed that gamma is emitted from (X,Y,Z) at tfoil+-dfoil
  //returns 1 if success, Eg1  -- total gamma energy MeV
  //Another algorithm

	  Float_t pmax=0.;
	  Float_t Eg1,Eg2;
	  Int_t isc1,isc2;
	  Bool_t onegamma=false;
	  //Check that event topology is correct
	  critery=1.E30;
	  if(Nbr_tks!=2) return 0;
	  if(Nsc<3) return 0;
	  Int_t nmax=0;
	  Int_t scint[15];
	  for(Int_t i=0;i<Nsc;i++){
	    if(i!=Ind_scintil[0]-1 && i!=Ind_scintil[1]-1){
	      scint[nmax]=i;
	      nmax++;
	    }
	    if(nmax>14) return 0;
	  }

	  Int_t gamma1[15],ng1;
	  Int_t gamma2[15],ng2;
	  Int_t nc=Int_t((pow(2,nmax))-0.5);
	  //	  std::cout<<"Gamma sc "<<nmax<<" "<<nc<<std::endl;

	    ng1=0;
	    ng2=0;
	    Eg1=0;
	    Eg2=0;
	    Int_t i=nc;
	    for(Int_t j=0;j<nmax;j++){
	      if(1&i>>j){
		gamma1[ng1]=scint[j];
		ng1++;
		Eg1+=Sc[scint[j]][8]*1E3;
	      }else{
		gamma2[ng2]=scint[j];
		ng2++;
		Eg2+=Sc[scint[j]][8]*1E3;
	      }
	    }  
	    Float_t chi21=chi2_tof1g(gamma1,ng1,X,Y,Z,tfoil,dfoil);
	    Float_t chi2=chi21;
	    Float_t prob=TMath::Prob(chi2,ng1+ng2);
	      Egamma1=Eg1;
	      pmax=prob;
	      critery=chi2/(ng1+ng2);

	      //critery=pmax;
	  if(pmax>0.000005) onegamma=true;
	  std::cout<<"One gamma prob "<<pmax<<" E: "<<Egamma1<<" "<<nmax<<std::endl; 	  
	    if(onegamma) return 1;
	  return 0;
}


Float_t h10::TOF2g_int(Int_t isc1,Int_t isc2,Float_t X,Float_t Y,Float_t Z,Float_t tfoil,Float_t dtfoil){
  // This function calculated probability that 2 gammas where emmited from 
  //one point (X,Y,Z) at the time tfoil+-dtfoil

  if(isc1==isc2) return 0;
  Float_t sig[2];
  Float_t sct[2];
  Float_t glen[2];
  Float_t sglen[2];
  sct[0]=getsctime(isc1,sig[0]);
  sct[1]=getsctime(isc2,sig[1]);

  glen[0]=getstraightlengthv(isc1,X,Y,Z,sglen[0])/30.;
  glen[1]=getstraightlengthv(isc2,X,Y,Z,sglen[1])/30.;
  //  std::cout<<"Gtof int "<<glen[0]<<" "<<glen[1]<<std::endl;
  //  std::cout<<" sig[0] "<<sig[0]<<" sglen[0] "<<sglen[0]<<" dtfoil "<<dtfoil<<std::endl;
  //  std::cout<<" sig[1] "<<sig[1]<<" sglen[1] "<<sglen[1]<<" dtfoil "<<dtfoil<<std::endl;

  sglen[0]=(sig[0]*sig[0]+sglen[0]*sglen[0]+dtfoil*dtfoil);
  sglen[1]=(sig[1]*sig[1]+sglen[1]*sglen[1]+dtfoil*dtfoil);
  Float_t chi2=pow(sct[0]-(tfoil+glen[0]),2)/sglen[0] + pow(sct[1]-(tfoil+glen[1]),2)/sglen[1];
  Float_t prob=TMath::Prob(chi2,2);
  //  std::cout<<" Prob "<<prob<<" chi2 "<<chi2<<std::endl;
  return prob;
}

Float_t h10::TOF2g_ext(Int_t isc1,Int_t isc2,Float_t X,Float_t Y,Float_t Z,Float_t tfoil,Float_t dtfoil){
  // This function calculated probability that 2 gammas in the event are just one external one.

  if(isc1==isc2) return 0;
  Float_t sig[2];
  Float_t sct[2];
  Float_t glen[2];
  Float_t sglen[2];
  sct[0]=getsctime(isc1,sig[0]);
  sct[1]=getsctime(isc2,sig[1]);
  glen[0]=getstraightlengthv(isc1,X,Y,Z,sglen[0])/30.;
  glen[1]=getstraightlengthv(isc2,X,Y,Z,sglen[1])/30.;
  //  std::cout<<"Gtof int "<<glen[0]<<" "<<glen[1]<<std::endl;

  sglen[0]=(sig[0]*sig[0]+sglen[0]*sglen[0]+dtfoil*dtfoil);
  sglen[1]=(sig[1]*sig[1]+sglen[1]*sglen[1]+dtfoil*dtfoil);
  //1st hyp isc1->isc2
  Float_t chi2_h1=pow(sct[0]-(tfoil-glen[0]),2)/sglen[0] + pow(sct[1]-(tfoil+glen[1]),2)/sglen[1];
  //2nd hyp isc2->isc1
  Float_t chi2_h2=pow(sct[0]-(tfoil+glen[0]),2)/sglen[0] + pow(sct[1]-(tfoil-glen[1]),2)/sglen[1];
  Float_t chi2;
  if(chi2_h1<chi2_h2) chi2=chi2_h1;
  if(chi2_h2<chi2_h1) chi2=chi2_h2;

  
  Float_t prob=TMath::Prob(chi2,2);
  return prob;
}

Float_t h10::TOF2g_scat(Int_t isc1,Int_t isc2,Float_t X,Float_t Y,Float_t Z,Float_t tfoil,Float_t dtfoil){
  // this function calculates probability that there is one gamma emitted from
  //X,Y,Z at tfoil+-dtfoil and fired scintillator isc1 and than isc2

  if(isc1==isc2) return 0;
  Float_t sig[2];
  Float_t sct[2];
  Float_t glen[2];
  Float_t sglen[2];
  sct[0]=getsctime(isc1,sig[0]);
  sct[1]=getsctime(isc2,sig[1]);
  glen[0]=getstraightlengthv(isc1,X,Y,Z,sglen[0])/30.;
  Float_t X1=Sc[isc1][9];
  Float_t Y1=Sc[isc1][10];
  Float_t Z1=Sc[isc1][11];
  glen[1]=glen[0]+getstraightlengthv(isc2,X1,Y1,Z1,sglen[1])/30.;
  //  std::cout<<"Gtof scat "<<glen[0]<<" "<<glen[1]<<std::endl;
  sglen[1]=sqrt(pow(sglen[0],2)+pow(2*sglen[1],2));
  sglen[0]=(sig[0]*sig[0]+sglen[0]*sglen[0]+dtfoil*dtfoil);
  sglen[1]=(sig[1]*sig[1]+sglen[1]*sglen[1]+dtfoil*dtfoil);
  
  Float_t chi2=pow(sct[0]-(tfoil+glen[0]),2)/sglen[0] + pow(sct[1]-(tfoil+glen[1]),2)/sglen[1];
  Float_t prob=TMath::Prob(chi2,2);
  return prob;
}

Float_t h10::chi2_tof1g(Int_t *scint, Int_t nscint,Float_t x,Float_t y,Float_t z, Float_t tfoil, Float_t dfoil)
{
  // This function checks if PMTs from scint array 
  // are hitted by one gamma emitted from (x,y,z) at tfoil+-dfoil
  // and rescattered inside NEMO.
  // It investigates all possible orderings of PMTs in scint array, 
  // and chooses the best chi2 value.
  // Returns chi2 value for nscint degree of freedom for further calculations.

  //nscint max=15
  Int_t order[15];
  Int_t scintcopy[15];
  Bool_t chkorder[15];
  Float_t chi2min=-1;
  if(nscint>6) return 15000.;

  for(Int_t i=0;i<15;i++){
    order[i]=0;
  }

  Int_t add=0;
  while(add==0){
    
    //scan all possible orders via loop
    //eg nscint=3:  000 100 200 010 110 210 020 120 ...
    add=1;
    //    std::cout<<"While ";
    for(Int_t i=0;i<nscint;i++){
      order[i]+=add;
      if(order[i]==nscint){
	order[i]=0;
	add=1;
      }else{
	add=0;
      }
      //    std::cout<<order[i];
    }
    //    std::cout<<std::endl;

    // check that there is no repetions, eg 123 is good order, 113 is not
    for(Int_t i=0;i<nscint;i++){
      chkorder[i]=false;
    }
    Bool_t goodorder=true;
    for(Int_t i=0;i<nscint;i++){
      if(chkorder[order[i]]){ 
	goodorder=false;
	continue;
      }else{
	chkorder[order[i]]=true;
      }	       
      scintcopy[i]=scint[order[i]];
    }
    //If order is correct calculate chi2
    if(goodorder){
      Float_t chi2=chi2_tof1g0(scintcopy,nscint,x,y,z,tfoil,dfoil);
      if(chi2<chi2min || chi2min<0) {
	chi2min =chi2;
      }
    }
  }
  return chi2min;
}

Float_t h10::chi2_tof1g0(Int_t *scint,Int_t nscint, Float_t x,Float_t y,Float_t z, Float_t tfoil, Float_t dfoil)
{
  // This function checks if PMTs from scint array, ordered as it is, 
  // are hitted by one gamma emitted from (x,y,z) at tfoil+-dfoil
  // and rescattered inside NEMO.
  // Returns chi2 value for nscint degree of freedom for further calculations.
  Float_t lx=x;
  Float_t ly=y;
  Float_t lz=z;
  Float_t chi2=0;
  Float_t path=0;
  Float_t dpath=0;
  Float_t dpathtmp;
  Float_t time,dtime;

  for(Int_t i=0;i<nscint;i++){
    Int_t isc=scint[i];
    path+=getstraightlengthv(isc,lx,ly,lz,dpathtmp);
    dpath+=dpathtmp*dpathtmp;
    time=getsctime(isc,dtime);
    chi2+=pow(time-(tfoil+path/30.),2)/(dpath+dtime*dtime+dfoil*dfoil);

    lx =Sc[isc][9];
    ly =Sc[isc][10];
    lz =Sc[isc][11];  
  }
  return chi2;
}

Float_t h10::getstraightlengthv(Int_t scint, Float_t X,Float_t Y,Float_t Z, Float_t& sigl){  
  Float_t c = 3e8;
  Float_t DEPTH,DELTAG;

  if(run > 1000) {
    //VK parameters to tune gamma internal prob.
    // see VK Feb 2007 talk in Orsay
    //real data
    DEPTH = 10.;
    DELTAG = 13;
  } else{
    //mc data
    DEPTH = 10.9;
    DELTAG = 8.5;
  }
  Float_t dsgamma;

  //  dsgamma = sqrt(pow(Sc[scint][9]-X,2)+pow(Sc[scint][10]-Y,2)+pow(Sc[scint][11]-Z,2));
  
 
  Float_t nx,ny;
  Float_t vx,vy;
  Float_t xx,yy,zz;
  xx = Sc[scint][9];
  yy = Sc[scint][10];
  zz = Sc[scint][11];
  Float_t theta = atan2(xx,yy);
    
  if(Sc[scint][2]==0){
    xx = xx - DEPTH * cos(theta);
    yy = yy - DEPTH * sin(theta);
  }
  if(Sc[scint][2]==1){
    xx = xx + DEPTH * cos(theta);
    yy = yy + DEPTH * sin(theta);
  }
  if(Sc[scint][2]==2){
    zz = zz - DEPTH;
  }
  if(Sc[scint][2]==3){
    zz = zz + DEPTH;
  }
  dsgamma = sqrt(pow(xx - X,2) +pow(yy - Y,2) +pow(zz - Z,2));
  sigl = DELTAG / c*1E7;
  
  /**  
  nx = -Sc[scint][10]/sqrt(pow(Sc[scint][9],2)+pow(Sc[scint][10],2));   
  ny = -Sc[scint][9]/sqrt(pow(Sc[scint][9],2)+pow(Sc[scint][10],2));
  vx = (Sc[scint][9]-X)/dsgamma;   
  vy = (Sc[scint][10]-Y)/dsgamma;
  Float_t cosa = (nx*vx+ny*vy);
   
  Float_t dxysc=7.5; // half dimention of mean face scintillator
  Float_t dzsc=15;// half dimension of scintillator depth
  Float_t delsgam =(2.*dsgamma*cosa*dxysc-pow(dxysc,2))/2./dsgamma;
  //delsgam = sqrt(15.0*15.0+19.0*19.0)/2.; //add on depth into scintillator on error
  delsgam=sqrt(delsgam*delsgam + dzsc*dzsc);
  sigl = delsgam/c*1e7; //in ns
  **/
  return dsgamma;
}

void  h10::TOF_bsigle(Float_t &tf, Float_t &sf,Int_t tk1)
{
  //This function assumes that the electron was emmited from the foil
  //Based on that time on the foil tfoil+-dtfoil is calculated.

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

  Int_t isc=Ind_scintil[tk[0]]-1;
  //electron energy
  Float_t eb=Sc[isc][8]*1000.;
  //the electron speed
  Float_t ve=sqrt(eb*eb+2.*rme*eb)/(rme+eb)*vlight;
  // the electron path
  Float_t dstrack=gethelixl(tk[0],dlen);
  Float_t fwhm;
  Float_t sige;
  //the electron time of flight
  tofb[0]=dstrack/ve*1.E9;
  // the electron measured arrival time
  tofbm[0]=getsctimebb(isc,sige);
  // error on time measurement in PM 
  dtofb[0]=sige;

  if(Sc[isc][2]==1){
    fwhm=0.14;
  }else{
    fwhm=0.17;
  }
  Float_t se1 = fwhm*sqrt(eb)/2.354; 
  //error due to electron speed uncertainty
  sige = tofb[0]*rme*rme/eb/(eb+rme)/(eb+2*rme)*se1;

  dtofb[0]=sqrt(sige*sige+dtofb[0]*dtofb[0]);

  tf=(tofbm[0]-tofb[0]);
  sf=dtofb[0];
  tf=tf;
  return;  
}


Float_t   h10::TOF2b_gas_int(Float_t &tf, Float_t &sf,Float_t *vtx,Int_t tk1, Int_t tk2)
{
  //Calculates TOF prpobability for internal hypothesis for two electrons
  // shen vertex is in the gas
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
    Float_t dstrack=gethelixl(tk[j],vtx,dlen);
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
  //  std::cout<<"T1-T2 "<<dtt<<" TOF1,2 "<<tofb[0]<<","<<tofb[1]<<" S1,2 "<<dtofb[0]<<","<<dtofb[1]<<std::endl;
  return prob;  
}

Float_t   h10::TOF2b_gas_ext(Float_t *vtx,Int_t tk1, Int_t tk2)
{
  //Calculates TOF prpobability for external hypothesis for two electrons
  // when event vertex is in the gas volume
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

  //No energy losses in the gas
    Float_t floss=0.0;
  

  //1st hypothesis 0->1
  for(Int_t j=0;j<2;j++){
    Int_t isc=Ind_scintil[tk[1]]-1;
    eb=Sc[isc][8]*1000.;
    //    if(j==0) {std::cout<<" eb before "<<eb<<std::endl;eb+=floss; std::cout<<" adding floss"<<floss<<" eb after "<<eb<<std::endl;}

    Float_t ve=sqrt(eb*eb+2.*rme*eb)/(rme+eb)*vlight;
    Float_t dstrack=gethelixl(tk[j],vtx,dlen);
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
  Float_t chiint2=dtt*dtt/(dtofb[1]*dtofb[1]+dtofb[0]*dtofb[0]+sige*sige);
  Float_t prob1 = TMath::Prob(chiint2,1); 


  //2nd hypothesis 1->0
  for(Int_t j=0;j<2;j++){
    Int_t isc=Ind_scintil[tk[0]]-1;
    eb=Sc[isc][8]*1000.;
    Float_t ve=sqrt(eb*eb+2.*rme*eb)/(rme+eb)*vlight;
    Float_t dstrack=gethelixl(tk[j],vtx,dlen);
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
  chiint2=dtt*dtt/(dtofb[1]*dtofb[1]+dtofb[0]*dtofb[0]+sige*sige);
  Float_t prob2 = TMath::Prob(chiint2,1); 
  if(prob1<prob2) prob=prob2;
  if(prob2<prob1) prob=prob1;

  return prob;
}


Float_t h10::gethelixl(Int_t tr, Float_t* vtx,Float_t &dlen)
{

  //Find electron track length strating from a middle point, vtx in the 
  //tracking volume. Not from the foil

      Float_t c = 3e8;
      Float_t me = 0.511;
      Float_t fwhm;
      Float_t beta, dbeta;
      Float_t vsc[2];
      Float_t vvr[2];

      vsc[0]=X_scintil[tr]-Xc[tr];
      vsc[1]=Y_scintil[tr]-Yc[tr];

      vvr[0]=vtx[0]-Xc[tr];
      vvr[1]=vtx[1]-Yc[tr];

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

//***************************************************************
//
// calculate approx sc resolution for electron speed uncertainty
//
//***************************************************************
Float_t h10::getscfwhm_real(Int_t isc){
  Float_t fwhm;
  if(Sc[isc][2]==0) fwhm=0.173;
  if(Sc[isc][2]==1) fwhm=0.142;
  if(Sc[isc][2]>1){
    if(Sc[isc][3]<3) fwhm=0.163;
    if(Sc[isc][3]>2) fwhm=0.146;
  }
  return fwhm;
}

