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

Bool_t h10::preselect_1tk(int eventnumber, Int_t* results)
{

   Float_t sigl2;
   Float_t Eelectron, Ebiggam, ee, dee, dtimens, dbeta, sigt1, sigt2, sigl1;
   Float_t c=3e8;
   Bool_t foil, secsig, secback, onetrack=false;
   Bool_t onetrack0, onetrack1, onetrack2, onetrack3, onetrack4, onetrack5;
   Bool_t onetrack6, onetrack7, onetrack8, onetrack9, onetrack10, onetrack11;
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

   //Int_t ientry = LoadTree(eventnumber);
   //if (ientry < 0) return onetrack;
   //Int_t nb = fChain->GetEntry(eventnumber);

   //one track only
   if(Nsc>50)return onetrack;
   onetrack0=(Nbr_tks==1);
   if(!onetrack0)return onetrack;
   onetrack = onetrack0;

   //both assigned to scint hits
   if(onetrack)results[0]++;
   onetrack1 = (Ind_scintil[0]!=0);
   onetrack = onetrack&&onetrack1;
   if(onetrack)results[1]++;

   Float_t dvert, dvertz, dvertxy, vertx, verty, vertz;

   vertx=(X_foil[0]);
   verty=(Y_foil[0]);
   vertz=(Z_foil[0]);

   Bool_t disqualify=true;
   //   if (sqrt(X_foil[0]*X_foil[0]+Y_foil[0]*Y_foil[0])>150&&(sqrt(X_foil[1]*X_foil[1]+Y_foil[1]*Y_foil[1])>150))disqualify=false;

   //ensure tracks come from the same vertex
   //onetrack2=(dvertxy<6.0&&dvertz<8.0&&fabs(vertz)<125.0&&!disqualify);
   //onetrack2=(dvertxy<4.0&&dvertz<4.0&&fabs(vertz)<125.0&&!disqualify);
   //onetrack = onetrack&&onetrack2;
   if (onetrack)results[2]++;

   //both tracks must have hits in either layer0 or layer1
   Bool_t foil1 = checklayer(0,0)||checklayer(0,1)||checklayer(0,2)||checklayer(0,3);
   onetrack3 = foil1;
   onetrack = onetrack&&onetrack3;
   if(onetrack) results[3]++;

   //cluster up the neighbouring hits (amIalone(i,j,0)=side neighbours,
   //                                  amIalone(i,j,1)=diagonal (corner) neighbours
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
   onetrack4 = Isocluster[Ind_scintil[0]-1];
   onetrack = onetrack&&onetrack4;
   if(onetrack) results[4]++;
   //scintillator deposits for tracks > 0.1MeV
   E1 = Sc[Ind_scintil[0]-1][8]*1000.;
   onetrack5 = (E1>0.1);
   onetrack = onetrack&&onetrack5;
   if(onetrack) results[5]++;

   //are both tracks electrons?
   onetrack6 = (Qc[0]<0);  
   onetrack = onetrack&&onetrack6;
   if (onetrack)results[6]++;

	//Int_t ifoil = Decodesrc(vertx,verty,vertz);
        //onetrack7 = (ifoil==0||ifoil==1);
        //onetrack = onetrack&&onetrack7;
        //if(onetrack)results[7]++;


   Float_t tfoil, dfoil, dummy; //put this here temporarily
   Float_t pull, tf[2];

   //track lengths > 30cm
   Float_t ds;
   Float_t dlen1 = gethelixl(0,ds);
   onetrack9=(dlen1>30);
   onetrack = onetrack&&onetrack9;
   if(onetrack) results[8]++;
   //is there at least one extra scintillator deposit
   onetrack10=(ncluster>1);
   onetrack = onetrack&&onetrack10;
   if(onetrack) results[9]++;
   return onetrack;
}

Bool_t h10::preselect_1e1g(int eventnumber, Int_t* results)
{
  //preselects events for 1e1g analysis
  // demand hit in 1st 2 layers of GG
  // Electron energy >200 keV 
  // Gamma energy > 200 keV
  // At list one gamma scintillator >125 keV

   Float_t c=3e8;
   Bool_t foil, onetrack=false;
   Bool_t onetrack0, onetrack1, onetrack2, onetrack3, onetrack4, onetrack5;
   Bool_t onetrack6, onetrack7, onetrack8, onetrack9, onetrack10, onetrack11;
   Float_t ebmin=0.2;
   Float_t egmin=0.2;
   Float_t eseedmin=0.125;
   Float_t etotmin=0.4;

   Bool_t  mc=false;
   Float_t E1, E2;

   //Int_t ientry = LoadTree(eventnumber);
   //if (ientry < 0) return onetrack;
   //Int_t nb = fChain->GetEntry(eventnumber);

   //one track only
   if(Nsc>50)return onetrack;
   onetrack0=(Nbr_tks==1);
   if(!onetrack0)return onetrack;
   onetrack = onetrack0;
   if(onetrack)results[0]++;

   //track assigned to scint hits
   onetrack1 = (Ind_scintil[0]!=0);
   onetrack = onetrack&&onetrack1;
   if(!onetrack) return onetrack;
   if(onetrack)results[1]++;

   Float_t dvert, dvertz, dvertxy, vertx, verty, vertz;

   vertx=(X_foil[0]);
   verty=(Y_foil[0]);
   vertz=(Z_foil[0]);
   //check if vertex is in the foil
   onetrack2=((vertx*vertx+verty*verty)>100)&&fabs(vertz)<125;
   onetrack=onetrack2&&onetrack;
   if(!onetrack) return onetrack;
   if (onetrack)results[2]++;

   //track must have hits in either layer0 or layer1
   Bool_t foil1 = checklayer(0,0)||checklayer(0,1);
   onetrack3 = foil1;
   onetrack = onetrack&&onetrack3;
   if(!onetrack) return onetrack;
   if(onetrack) results[3]++;

   //cluster up the neighbouring hits (amIalone(i,j,0)=side neighbours,
   //                                  amIalone(i,j,1)=diagonal (corner) neighbours
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
   onetrack4 = Isocluster[Ind_scintil[0]-1];
   onetrack = onetrack&&onetrack4;
   if(!onetrack) return onetrack;
   if(onetrack) results[4]++;

   //scintillator deposits for tracks > ebmin
   E1 = Sc[Ind_scintil[0]-1][8]*1000.;
   onetrack5 = (E1>ebmin);
   onetrack = onetrack&&onetrack5;
   if(!onetrack) return onetrack;
   if(onetrack) results[5]++;

   //gamma deposit= total-electron deposit > egmin
   //and theres one seed > 0.125 keV
   E2=0;
   Bool_t threshold=false;
   for(Int_t i=0;i<Nsc;i++){
     E2+=Sc[i][8]*1000;
     threshold=threshold ||((Sc[i][8]*1000>eseedmin)&&(Ind_scintil[0]-1!=i));
   }
   onetrack5=(E2-E1>egmin)&&threshold;
   onetrack=onetrack&&onetrack5;
   if(!onetrack) return onetrack;
   if(onetrack) results[5]++;


   //are both tracks electrons?
   onetrack6 = (Qc[0]<0);  
   onetrack = onetrack&&onetrack6;
   if(!onetrack) return onetrack;
   if (onetrack)results[6]++;

   //track lengths > 30cm
   Float_t ds;
   Float_t dlen1 = gethelixl(0,ds);
   onetrack9=(dlen1>30);
   onetrack = onetrack&&onetrack9;
   if(onetrack) results[8]++;
   //is there at least one extra scintillator deposit
   onetrack10=(ncluster>1);
   onetrack = onetrack&&onetrack10;
   if(onetrack) results[9]++;
   return onetrack;
}
Bool_t h10::preselect_1e(int eventnumber, Int_t* results)
{

   Float_t sigl2;
   Float_t Eelectron, Ebiggam, ee, dee, dtimens, dbeta, sigt1, sigt2, sigl1;
   Float_t c=3e8;
   Bool_t foil, secsig, secback, onetrack=false;
   Bool_t onetrack0, onetrack1, onetrack2, onetrack3, onetrack4, onetrack5;
   Bool_t onetrack6, onetrack7, onetrack8, onetrack9, onetrack10, onetrack11;
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

   //Int_t ientry = LoadTree(eventnumber);
   //if (ientry < 0) return onetrack;
   //Int_t nb = fChain->GetEntry(eventnumber);

   //one track only
   if(Nsc>50)return onetrack;
   onetrack0=(Nbr_tks==1);
   if(!onetrack0)return onetrack;
   onetrack = onetrack0;

   //both assigned to scint hits
   if(onetrack)results[0]++;
   onetrack1 = (Ind_scintil[0]!=0);
   onetrack = onetrack&&onetrack1;
   if(onetrack)results[1]++;

   Float_t dvert, dvertz, dvertxy, vertx, verty, vertz;

   vertx=(X_foil[0]);
   verty=(Y_foil[0]);
   vertz=(Z_foil[0]);

   Bool_t disqualify=true;
   //   if (sqrt(X_foil[0]*X_foil[0]+Y_foil[0]*Y_foil[0])>150&&(sqrt(X_foil[1]*X_foil[1]+Y_foil[1]*Y_foil[1])>150))disqualify=false;

   //ensure tracks come from the same vertex
   //onetrack2=(dvertxy<6.0&&dvertz<8.0&&fabs(vertz)<125.0&&!disqualify);
   //onetrack2=(dvertxy<4.0&&dvertz<4.0&&fabs(vertz)<125.0&&!disqualify);
   //onetrack = onetrack&&onetrack2;
   if (onetrack)results[2]++;

   //both tracks must have hits in either layer0 or layer1
   Bool_t foil1 = checklayer(0,0)||checklayer(0,1)||checklayer(0,2)||checklayer(0,3);
   onetrack3 = foil1;
   onetrack = onetrack&&onetrack3;
   if(onetrack) results[3]++;

   //cluster up the neighbouring hits (amIalone(i,j,0)=side neighbours,
   //                                  amIalone(i,j,1)=diagonal (corner) neighbours
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
   onetrack4 = Isocluster[Ind_scintil[0]-1];
   onetrack = onetrack&&onetrack4;
   if(onetrack) results[4]++;
   //scintillator deposits for tracks > 0.1MeV
   E1 = Sc[Ind_scintil[0]-1][8]*1000.;
   onetrack5 = (E1>0.1);
   onetrack = onetrack&&onetrack5;
   if(onetrack) results[5]++;

   //are both tracks electrons?
   onetrack6 = (Qc[0]<0);  
   onetrack = onetrack&&onetrack6;
   if (onetrack)results[6]++;

   Float_t tfoil, dfoil, dummy; //put this here temporarily
   Float_t pull, tf[2];

   //track lengths > 30cm
   Float_t ds;
   Float_t dlen1 = gethelixl(0,ds);
   onetrack7=(dlen1>30);
   onetrack = onetrack&&onetrack7;
   if(onetrack) results[7]++;
   return onetrack;
}
