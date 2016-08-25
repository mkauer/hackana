#include "h10.h"
#include "TH2.h"
#include "TStyle.h"
#include "TNtuple.h"
#include "TCanvas.h"
#include "TMath.h"

#include <fstream>
#include <iostream>
#include <cmath>
#include <stdlib.h>

Bool_t h10::Preselect(int eventnumber, Int_t* results)
{

   Float_t sigl2;
   Float_t Eelectron, Ebiggam, ee, dee, dtimens, dbeta, sigt1, sigt2, sigl1;
   Float_t c=3e8;
   Bool_t foil, secsig, secback, twotrack=false;
   Bool_t twotrack0, twotrack1, twotrack2, twotrack3, twotrack4, twotrack5;
   Bool_t twotrack6, twotrack7, twotrack8, twotrack9, twotrack10, twotrack11;
   Float_t dstrack, chisq, chi, beta, timens, T, T0;
   Float_t dsegment;
   Float_t dsgamma, tdgamma[10];
   Float_t probcut = 0.1;
   Float_t delseccut = 8;
   Float_t eelecmin = 0.3;


   Bool_t  mc=false;
   Float_t biggamlow = 1.7;
   Float_t lin1[10];
   Float_t gprob[10], lprob;
   Int_t alphas = 2, alparray[1000];
   Int_t tot[4];
   Int_t nbytes = 0;
   Float_t xf, yf, Escint, E1, E2;

   //two tracks only
   if(Nsc>50)return twotrack;
   twotrack0=(Nbr_tks==2);
   if(!twotrack0)return twotrack;
   twotrack = twotrack0;

   //both assigned to scint hits
   //if(twotrack)results[0]++;
   //twotrack1 = (Ind_scintil[0]!=0&&Ind_scintil[1]!=0);
   //twotrack1 = twotrack1&&Ind_scintil[0]!=Ind_scintil[1];
   //twotrack = twotrack&&twotrack1;
   //if(twotrack)results[1]++;

   Float_t dvert, dvertz, dvertxy, vertx, verty, vertz;

   vertx=(X_foil[0]+X_foil[1])/2;
   verty=(Y_foil[0]+Y_foil[1])/2;
   vertz=(Z_foil[0]+Z_foil[1])/2;

   dvertxy = sqrt(pow(X_foil[0]-X_foil[1],2)+pow(Y_foil[0]-Y_foil[1],2));
   dvertz=sqrt(pow(Z_foil[0]-Z_foil[1],2));

   Bool_t disqualify=true;
   if (sqrt(X_foil[0]*X_foil[0]+Y_foil[0]*Y_foil[0])>150&&(sqrt(X_foil[1]*X_foil[1]+Y_foil[1]*Y_foil[1])>150))disqualify=false;

   //ensure tracks come from the same vertex
   //   twotrack2=(dvertxy<6.0&&dvertz<8.0&&fabs(vertz)<125.0&&!disqualify);
   twotrack2=(dvertxy<4.0&&dvertz<4.0&&fabs(vertz)<125.0&&!disqualify);
   twotrack = twotrack&&twotrack2;
   if (twotrack)results[2]++;

   //both tracks must have hits in either layer0 or layer1
   //Bool_t foil1 = checkfirstlayer(0,0);
   //Bool_t foil2 = checkfirstlayer(1,0);
   //twotrack3 = foil1&&foil2;
   //twotrack = twotrack&&twotrack3;
   //if (twotrack)results[3]++;

   //cluster up the neighbouring hits (amIalone(i,j,0)=nearest neighbours,
   //                                  amIalone(i,j,1)=diagonal neighbours
   Int_t ncluster=0;
   Float_t Ecluster[50]={0};
   Int_t Icluster[50]={0}, nisol=0;
   Bool_t Isocluster[50];
   for (Int_t i=0; i<Nsc; i++)
      {  Isocluster[i] = true;
         for (Int_t j=0;j<Nsc; j++)
            { if (i!=j)
	        { 
	          Bool_t solo = amIalone(i,j,0);
                  Isocluster[i] = Isocluster[i] && solo;
                 }
             }
         if (Sc[i][8]*1000. > 0.0)
	    { Ecluster[ncluster] = Sc[i][8]*1000.;
              Icluster[ncluster] = i;
              ncluster++;
            } 
        }

   //are the scint hits associated with the tracks isolated?
   twotrack4 = (Isocluster[Ind_scintil[0]-1]&&Isocluster[Ind_scintil[1]-1]);
   twotrack = twotrack&&twotrack4;

   //scintillator deposits for tracks > 0.1MeV
   //E1 = Sc[Ind_scintil[0]-1][8]*1000.;
   //E2 = Sc[Ind_scintil[1]-1][8]*1000.;
   //twotrack5 = (E1>0.1&&E2>0.1);
   //twotrack = twotrack&&twotrack5;
   //if (twotrack)results[5]++;

   //are both tracks electrons?
   twotrack6 = (Qc[0]<0&&Qc[1]<0);  
   twotrack = twotrack&&twotrack6;
   if (twotrack)results[6]++;

	//	Int_t ifoil = Decodesrc(vertx,verty,vertz);
        //twotrack7 = (ifoil==0||ifoil==1);
        //twotrack = twotrack&&twotrack7;
        //if(twotrack)results[7]++;


   Float_t tfoil, dfoil, dummy; //put this here temporarily
   Float_t pull, tf[2];
   //check internal hypothesis for two tracks
   Float_t prob2tk = gethypoth(0,1,0,0,tf);
   //check external hypothesis for two tracks
   Float_t prob2tkext1 = gethypoth(0,1,0,1,tf);
   Float_t prob2tkext2 = gethypoth(0,1,0,2,tf);
   //hypothesis consistent with internal and not external
   twotrack8=(prob2tk>0.01&&prob2tkext1<0.001&&prob2tkext2<0.001);
   twotrack = twotrack&&twotrack8;
   if(twotrack)results[8]++;

   //track lengths > 30cm
   //Float_t ds;
   //Float_t dlen1 = gethelixl(0,ds);
   //Float_t dlen2 = gethelixl(1,ds);
   //twotrack9=(dlen1>30&&dlen2>30);
   //twotrack = twotrack&&twotrack9;
   //if(twotrack)results[9]++;

   //is there at least one extra scintillator deposit
   twotrack10=(ncluster>2);
   twotrack = twotrack&&twotrack10;
   if (twotrack)results[10]++;

   return twotrack;
}
