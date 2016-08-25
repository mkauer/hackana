#define h10_cxx
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

using namespace std;

class flap
{ public:
  Int_t once;
  Int_t twice;
};

class stupid
{ public:
  flap george;
  Float_t johnny;
  Float_t mary(int i, int j)
  { return i*j;
  }
};

int main(int argc, char* argv[])
{ //argv[0] is name of programme 

 Bool_t filelist=true;
 TTree* tree=0;

 stupid bloody;
 std::cout<<" inside slim"<<std::endl;

 if (filelist)  
   {  std::cout<<" here "<<std::endl;
      char control[256];
      Int_t numFiles;
      ifstream in0, in1, in2;
      in0.open("control.dat");
      in0 >> numFiles >> control;
      std::cout<<" control "<<control<<" numfile "<<numFiles<<std::endl;
      in0.close();
      in0.open(control);      
      char filenames[200][256];
      Int_t runnumber[200], datenumber[200];

      for (int i=0;i<numFiles;i++)
         {  in0 >> filenames[i];
         }

      if (tree==0)
        {  TChain *chain = new TChain("h10");
	   Int_t run;
           for (int i=0;i<numFiles;i++) 	     
              { 
                chain->Add(filenames[i]);	      
	        std::cout<<filenames[i]<<std::endl;
                std::cout<<"Adding: "<<filenames[i]<<std::endl;
                Int_t nentries = Int_t(chain->GetEntries());
              }
           tree=chain;
        }

   }

  h10 george(tree);
  if(argc>1) //if main has an argument, its an event number
  { george.Loop(atoi(argv[1]));
    george.testmodule();
  }
  else
    { 
      std::cout<<" calling george"<<std::endl;
      george.Loop(0);
      std::cout<<"george"<<std::endl;
    }
}


void h10::Loop(Long64_t eventnumber)
{
   stupid bleeding;
   scintillator sci;
   TFile *newfile = new TFile("rooty.root","RECREATE");

   Int_t NGAM=2;
   Float_t sigl2;
   Float_t Eelectron, Ebiggam, ee, dee, dtimens, dbeta, sigt1, sigt2, sigl1;
   Float_t c=3e8;
   Bool_t foil, secsig, secback, twotrack;
   Bool_t twotrack0, twotrack1, twotrack2, twotrack3, twotrack4, twotrack5;
   Bool_t twotrack6, twotrack7, twotrack8, twotrack9, twotrack10, twotrack11;
   Bool_t mo100x, mo100x0, mo100x1, mo100x2, mo100x3, mo100x4, mo100x5;
   Bool_t mo100x6, mo100x7, mo100x8, mo100x9, mo100x10, mo100x11, mo100x12;

   Float_t dstrack, chisq, chi, beta, timens, T, T0;
   Float_t dsegment;
   Float_t dsgamma, tdgamma[10];
   Float_t probcut = 0.1;
   Float_t delseccut = 8;
   Float_t Egmin=0.15;
   Float_t Eemin=0.15;
   Float_t Eeemaxcut=2.0;
   Float_t Eeemincut=0.2;
   Float_t Egmax=0.59;
   Float_t Eemax=0.55;
   Float_t Egamsummin=0.35;
   Float_t Egamsummax=1.0;

   TH1F *deltata = new TH1F("deltata","deltata",200,-20,20);
   TH1F *deltata2 = new TH1F("deltata2","deltata2",200,-20,20);
   TH1F *chisq1 = new TH1F("chisq1","chisq1",100,0,10);
   TH1F *chisq2 = new TH1F("chisq2","chisq2",100,0,10);
   TH2F *egamvsee = new TH2F("egamvsee","egamvsee",100,0,5,100,0,5);
   TH2F *e1vse2 = new TH2F("e1vse2","e1vse2",100,0,2,100,0,2);
   TH2F *e1vse2a = new TH2F("e1vse2a","e1vse2a",100,0,2,100,0,2);
   TH2F *pscatvsdelt = new TH2F("pscatvsdelt","pscatvsdelt",100,0,1,100,-12,12);
   TH2F *scintxy = new TH2F("scintxy","scintxy",100,-300,300,100,-300,300);
   TH2F *scintrz1 = new TH2F("scintrz1","scintrz1",100,0,300,100,-300,300);
   TH2F *scintrz2 = new TH2F("scintrz2","scintrz2",100,0,300,100,-300,300);
   TH1F *probby1 = new TH1F("probby1","probby1",200,0,1);
   TH1F *lprobint = new TH1F("lprobint","lprobint",100,0,20);
   TH1F *lprobext = new TH1F("lprobext","lprobext",100,0,20);
   TH1F *lprobscat = new TH1F("lprobscat","lprobscat",100,0,20);
   TH1F *logprob1 = new TH1F("logprob1","logprob1",100,0,20);
   TH1F *logprob2 = new TH1F("logprob2","logprob2",100,0,20);
   TH1F *sector2 = new TH1F("sector2","sector2",21,0,20);
   TH1F *sector5 = new TH1F("sector5","sector5",21,0,20);
   TH1F *logprob3 = new TH1F("logprob3","logprob3",100,0,20);
   TH1F *logprob4 = new TH1F("logprob4","logprob4",100,0,20);
   TH2F *loglog = new TH2F("loglog","loglog",100,0,20,100,0,20);
   TH1F *probby3 = new TH1F("probby3","probby3",200,0,1);
   TH1F *probby2 = new TH1F("probby2","probby2",200,0,1);
   TH1F *deltat2 = new TH1F("deltat2","deltat2",200,-10,10);
   TH1F *deltat1 = new TH1F("deltat1","deltat1",200,-10,10);
   TH1F *deltat3 = new TH1F("deltat3","deltat3",200,-10,10);
   TH1F *deltat4 = new TH1F("deltat4","deltat4",200,-10,10);
   TH1F *probby4 = new TH1F("probby4","probby4",200,0,1);
   TH1F *Etime = new TH1F("Etime","Etime",100,0,100000);
   TH1F *Etotal = new TH1F("Etotal","Etotal",100,0.05,5);
   TH1F *Emini = new TH1F("Emini","Emini",100,0.0,2);
   TH1F *Eee = new TH1F("Eee","Eee",100,0.0,3);
   TH1F *Eg = new TH1F("Eg","Eg",100,0.0,3);
   TH1F *GamSum = new TH1F("GamSum","GamSum",100,0.0,5);
   TH1F *ngam = new TH1F("ngam","ngam",10,0.0,10);
   TH1F *Ngammas = new TH1F("Ngammas","Ngammas",10,0.0,10);
   TH1F *Ngammas1 = new TH1F("Ngammas1","Ngammas1",10,0.0,10);
   TH1F *Ngammas2 = new TH1F("Ngammas2","Ngammas2",10,0.0,10);
   TH1F *Nscint = new TH1F("Nscint","Nscint",10,0.0,10);
   TH1F *Nscint2 = new TH1F("Nscint2","Nscint2",10,0.0,10);
   TH1F *Nalphas = new TH1F("Nalphas","Nalphas",10,0.0,10);
   TH1F *cost = new TH1F("cost","cost",100,-1.0,1.0);
   TH1F *cuts0 = new TH1F("cuts0","cuts0",100,0.0,3);
   TH1F *cuts1 = new TH1F("cuts1","cuts1",100,0.0,3);
   TH1F *cuts2 = new TH1F("cuts2","cuts2",100,0.0,3);
   TH1F *cuts3 = new TH1F("cuts3","cuts3",100,0.0,3);
   TH1F *cuts4 = new TH1F("cuts4","cuts4",100,0.0,3);
   TH1F *cuts5 = new TH1F("cuts5","cuts5",100,0.0,3);
   TH1F *cuts6 = new TH1F("cuts6","cuts6",100,0.0,3);
   TH1F *cuts7 = new TH1F("cuts7","cuts7",100,0.0,3);
   TH1F *cuts8 = new TH1F("cuts8","cuts8",100,0.0,5);
   TH1F *cuts9 = new TH1F("cuts9","cuts9",100,0.0,5);
   TH1F *cuts10 = new TH1F("cuts10","cuts10",100,0.0,5);
   TH1F *Emaxi = new TH1F("Emaxi","Emaxi",100,0.0,3);
   TH2F *lprobscatve = new TH2F("lprobscatve","lprobscatve",100,0,200,100,0,2);
   TH1F *probscat = new TH1F("probscat","probscat",100,0,1);
   TH1F *pullscat = new TH1F("pullscat","pullscat",100,-10,10);
   TH1F *biggam = new TH1F("biggam","biggam",250,0.2,4.0);
   TH2F *probivse = new TH2F("probivse","probivse",100,0,30,100,0,30);
   TH1F *biggam1 = new TH1F("biggam1","biggam1",250,0.2,4.0);
   TH1F *bigelec = new TH1F("bigelec","bigelec",250,0.2,4.0);
   TH1F *biggam2 = new TH1F("biggam2","biggam2",250,0.2,4.0);
   TH1F *biggam3 = new TH1F("biggam3","biggam3",250,0.2,4.0);
   TH1F *result3 = new TH1F("result3","result3",20,-0.5,19.5);
   TH1F *result1 = new TH1F("result1","result1",20,-0.5,19.5);
   TH1F *result2 = new TH1F("result2","result2",20,-0.5,19.5);
   TH1F *result4 = new TH1F("result4","result4",20,-0.5,19.5);
   TH1F *result5 = new TH1F("result5","result5",20,-0.5,19.5);
   TH1F *Electron = new TH1F("Electron","Electron",250,0.02,3);
   TH1F *Electronadc = new TH1F("Electronadc","Electronadc",250,0.0,1000);
   TH1F *charge = new TH1F("charge","charge",5,-2,2);
   TH1F *helix = new TH1F("helix","helix",100,0.,1.0);
   TH1F *distancex = new TH1F("distancex","distancex",100,0,50);
   TH1F *distancey = new TH1F("distancey","distancey",100,0,50);
   TH1F *distances = new TH1F("distances","distances",100,0,50);
   TH1F *deltaxyz = new TH1F("deltaxyz","deltaxyz",100,0,50);
   TH1F *vertex = new TH1F("vertex","vertex",100,0,50);
   TH1F *vertexy = new TH1F("vertexy","vertexy",100,0,50);
   TH1F *vertexz = new TH1F("vertexz","vertexz",100,0,50);
   TH2F *erat1 = new TH2F("erat1","erat1",100,0.0,3.0,100,0.0,3.0);
   TH2F *erat1a = new TH2F("erat1a","erat1a",100,0.0,3.0,100,0.0,3.0);
   TH1F *timeslow1 = new TH1F("timeslow1","timeslow1",100,0.0,800);
   TH1F *timeslow = new TH1F("timeslow","timeslow",100,0.0,800);
   

   char control[256];
   char rootfile[256];

   Int_t numFiles;

   Bool_t  mc=false;
   if (fChain == 0) return;

   Float_t biggamlow = 1.7;
   Float_t lin1[10], tfoil, dfoil;
   Float_t gprob[10], lprob;
   Int_t alphas = 2, alparray[1000];
   Int_t tot[4];
   Int_t nbytes = 0;
   Float_t xf, yf, Escint, E1, E2, Eg1, Eg2;

   //   Int_t results[20]={0}; this is equivalent to the following
   Int_t *results = new Int_t [20];
   for(Int_t i=0;i<11;i++) {
      results[i]=0;
   }
   Int_t results2[20]={0}; 


   Int_t nentries = Int_t(fChain->GetEntries());
   Int_t nb = -1;
   int treenum=-1;
   std::cout<<" number of entries "<<nentries<<std::endl;

   for (Int_t jentry=0; jentry<nentries;jentry++){
        Int_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        Int_t nb = fChain->GetEntry(jentry);
        Nscint->Fill(Nsc);
        Int_t alp = getalphas(alparray);
        Int_t nafgg=NAfasthits();
        checkdelayed();
        Float_t atime1,atime2;
        Int_t alpha1=checksingledhit(atime1);
        Int_t alpha2=checkgroupdhit(atime2);
	if(alpha1>0)timeslow1->Fill(atime1/1000.);
	if(alpha2>0)timeslow1->Fill(atime2/1000.);

	if (!Preselect(jentry, results))continue;
        Nscint2->Fill(Nsc);

        Int_t scintil1 = Ind_scintil[0]-1;
        Int_t scintil2 = Ind_scintil[1]-1;

        Float_t dvert, dvertz, dvertxy, vertx, verty, vertz, dvertxyz;

        vertx=(X_foil[0]+X_foil[1])/2;
        verty=(Y_foil[0]+Y_foil[1])/2;
        vertz=(Z_foil[0]+Z_foil[1])/2;

        Int_t sector = int(Calcsec(vertx, verty));

        dvertxy = sqrt(pow(X_foil[0]-X_foil[1],2)+pow(Y_foil[0]-Y_foil[1],2));
        dvertxyz = sqrt(pow(X_foil[0]-X_foil[1],2)+pow(Y_foil[0]-Y_foil[1],2)+pow(Z_foil[0]-Z_foil[1],2));
        dvertz=sqrt(pow(Z_foil[0]-Z_foil[1],2));
        Float_t prob2tk = gethypoth(0,1,0,0,&tfoil);

// 	mo100x0=(dvertxy<4.0&&dvertz<4.0&&abs(vertz)<125.0);
        mo100x0 = true;
        mo100x = mo100x0;

        E1 = Sc[scintil1][8]*1000.;
        E2 = Sc[scintil2][8]*1000.;
        Float_t Etwo = E1 + E2;
        mo100x = mo100x&&(E1>Eemin&&E1<Eemax)&&(E2>Eemin&&E2<Eemax);
        if (mo100x)results2[0]++;

//find dot product of the two tracks
        Float_t costhe=0;
        for (Int_t i=0; i<3; i++)
	   { costhe += Cos_dir[0][i]*Cos_dir[1][i];
           }

        if(mo100x)
          { cuts0->Fill(Etwo);
            cost->Fill(costhe);
          }

        if (mo100x)results2[1]++;

        Bool_t Isocluster[50];
	Init_scintillator(&sci,Isocluster,0);
 
	mo100x2 = (Isocluster[scintil1]&&Isocluster[scintil2]);
        mo100x = mo100x&&mo100x2;
        if(!mo100x)continue;
        if(mo100x)results2[2]++;

        mo100x3=(sci.icluster>2);
        mo100x = mo100x&&mo100x3;
        if (mo100x)results2[3]++;
        
	Int_t ifoil = Decodesrc(vertx,verty,vertz);
        
        Float_t area1 = int(Sc[scintil1][2]);
        Float_t area2 = int(Sc[scintil2][2]);
        

        Eee->Fill(E1);
        Eee->Fill(E2);
        Bool_t wall1=(area1!=2&&area1!=3);
        Bool_t wall2=(area2!=2&&area2!=3);
        mo100x=mo100x&&wall1&&wall2;
        if(mo100x)
          { cuts1->Fill(Etwo);
          }

	sector2->Fill(sector);
        Float_t Egamsum = 0;
        Int_t gam=0; 
        Int_t iclusterold = sci.icluster;
        Int_t candgam[20]={20*0};

        if (mo100x)
          { for (Int_t is=0; is<iclusterold; is++)
 	       { if (sci.ncluster[is].Esum>Egmin&&sci.ncluster[is].Esum<Egmax&&!sci.ncluster[is].assigned)
 		   { Float_t chi, sigl1, delt;
                     Int_t seedy = sci.ncluster[is].Iseed;
 		     deltat1->Fill(chi);
 		     chisq1->Fill(chi*chi);
		     for(Int_t j=0; j<sci.icluster; j++)
//could this be a scattered gamma which needs adding to cluster?
//only if its in the wall scintillators
		       { if(!sci.ncluster[j].assigned&&(Sc[sci.ncluster[j].Iseed][2]==0||Sc[sci.ncluster[j].Iseed][2]==1))
                           { Int_t seedz = sci.ncluster[j].Iseed;
//make sure they arnt the same cluster
                             if (seedz!=seedy)
			       { 
//get their relative scatter hypothysis
                                 Float_t pscat=getscatterhyp(seedz,seedy,delt);
			         deltat2->Fill(delt);
			         if(pscat>0.005)
                                   { addclusters(&sci,is,j);
			           }
                                 Float_t Lpscat1 = -TMath::Log10(pscat);
                                 lprobscat->Fill(Lpscat1);
                                 pscatvsdelt->Fill(pscat,delt);
			       }
			   }
		       }
		   }
	      }
	  }

        if(mo100x)Ngammas->Fill(sci.icluster);
	//get the probs for the clusters
        Int_t gcount=0;
        if (mo100x)
          { for (Int_t is=0; is<sci.icluster; is++)
 	       { if (!sci.ncluster[is].assigned)
 		   { Int_t seedy = sci.ncluster[is].Iseed;
                     sci.ncluster[is].probint = getsinglehyp(tfoil,dfoil,vertx,verty,vertz,dvertxyz/2.,seedy,1,0,chi,sci.ncluster[is].Esum);
                     sci.ncluster[is].probext = getsinglehyp(tfoil,dfoil,vertx,verty,vertz,dvertxyz/2.,seedy,1,1,chi,sci.ncluster[is].Esum);
                     Float_t Lpint = -TMath::Log10(sci.ncluster[is].probint);
                     Float_t Lpext = -TMath::Log10(sci.ncluster[is].probext);
                     lprobint->Fill(Lpint);
                     lprobext->Fill(Lpext);
                     loglog->Fill(Lpint,Lpext);
                     gcount++;
                     if(gcount==1)Eg1=sci.ncluster[is].Esum;
                     if(gcount==2)Eg2=sci.ncluster[is].Esum;
                   }
               }
	  }
        if(mo100x)Ngammas1->Fill(sci.icluster);

        mo100x = mo100x&&(sci.icluster==(2+NGAM));
	if(!mo100x)continue;
        Bool_t inthyp=true;
        Bool_t exthyp=true;
	for(Int_t i=0;i<sci.icluster;i++)
	  { if (!sci.ncluster[i].assigned)
              { inthyp = inthyp&&(sci.ncluster[i].probint>0.04);
	        exthyp = exthyp&&(sci.ncluster[i].probext>=0.0);
//	        exthyp = exthyp&&(sci.ncluster[i].probext<0.04);
              }
          }
        mo100x=mo100x&&inthyp&&exthyp;
        if(!mo100x)continue;
        if(mo100x)Ngammas2->Fill(sci.icluster);
        if (mo100x)results2[4]++;
	Float_t Egam=0;
        Int_t igam=0;
        Int_t jgam=0;
//check that there are the correct number of gamma clusters with >Egmax which are not assigned to the track which have a good internal prob

        for(Int_t k=0; k<sci.icluster;k++)
	  { if(!sci.ncluster[k].assigned)
	      { Egam+=sci.ncluster[k].Esum;
	        Eg->Fill(sci.ncluster[k].Esum);
                jgam++;
                if (sci.ncluster[k].Esum<Egmax)igam++;
	      }
	  }

        erat1->Fill(Etwo,Egam);
	GamSum->Fill(Egam);
	//        mo100x5=(Etwo<Eeemaxcut&&Egam>0.25*Etwo&&igam==2&&jgam==2);
	//        mo100x5=(Etwo<Eeemaxcut&&igam==2&&jgam==2);

        mo100x5=(igam==2&&jgam==2);
	mo100x4=(ifoil==0||ifoil==1);
        mo100x=mo100x&&mo100x5;
        if (mo100x)results2[5]++;

        mo100x6=!(alpha1>0||alpha2>0);
	//        if(alpha1>0||alpha2>0||nafgg>0)

        if (mo100x&&mo100x6)
          { cuts3->Fill(Etwo);
   	    cuts6->Fill(Egam);
            cuts9->Fill(Etwo+Egam);
            e1vse2->Fill(E1,E2);
            result4->Fill(ifoil);
	  }

        if (mo100x)results2[6]++;
	//        mo100x=mo100x&&(E2<=(-E1+1.2));
	//	mo100x7=(Egam<Egamsummax&&Egam>Egamsummin);
        mo100x7 = abs(E1-E2)<0.25&&abs(Eg1-Eg2)<0.25;
   //        mo100x7=mo100x7&&costhe<-0.3;

        if (mo100x&&mo100x7&&mo100x6)
          { cuts4->Fill(Etwo);
   	    cuts7->Fill(Egam);
            result2->Fill(sector);
            result5->Fill(ifoil);
            cuts10->Fill(Etwo+Egam);
            results2[7]++;
	  }


        if (mo100x&&mo100x7&&!mo100x6)
          { result3->Fill(sector);
            cuts2->Fill(Etwo);
	    cuts5->Fill(Egam);
	    if(alpha1>0)timeslow->Fill(atime1/1000.);
	    if(alpha2>0)timeslow->Fill(atime2/1000.);
            result1->Fill(ifoil);
            cuts8->Fill(Etwo+Egam);
            erat1a->Fill(Etwo,Egam);
            e1vse2a->Fill(E1,E2);
          }
   }


   newfile->Write();
   newfile->Close(); 
   for(Int_t i=0;i<11;i++)
     { std::cout<<" results for cut "<<i<<" "<<results[i]<<std::endl;
     } 
   for(Int_t i=0;i<11;i++)
     { std::cout<<" results2 for cut "<<i<<" "<<results2[i]<<std::endl;
     } 
}

