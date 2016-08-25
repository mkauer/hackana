#ifndef ANA10_H
#define ANA10_H


//includes for RooFit
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooKeysPdf.h"
#include "RooHistPdf.h"
#include "Roo2DKeysPdf.h"
#include "RooPlot.h"
#include "RooFitResult.h"
#include "RooAddPdf.h"
#include "RooExtendPdf.h"
#include "RooNLLVar.h"
#include "RooMCStudy.h"


#include <string>
// Include database interface header
//#include "N3Db.h"
//#include "N3DbMgr.h"

#include "../h10.h"
#include "../helen.hpp"
#include "util.hpp"

#include "TFoam.h"
#include "TH2.h"
#include "TH1.h"
#include "TMap.h"
#include "TTree.h"
#include "TObjArray.h"
#include "TKey.h"
#include "TString.h"
#include "TRegexp.h"
#include "TPRegexp.h"

#include <map>
#include <iostream>
#include <fstream>
#include <vector>
//#include <vector.h>

//stream for analysis output
//defined in ana_main

extern ofstream anaout;
extern ofstream erlog;  // Matt Kauer added this in!!
//extern Pmstatus *pmstatus;


//NEMO DB interface class, initialised in main()
//extern nemo3::N3Db nemodb;



using namespace std;

typedef map<string,TH1D *> TH1Dmap;
typedef map<string,TH2D *> TH2Dmap;

typedef map<Int_t,Float_t> Timesmap;

  const Int_t MAXCHNLS=10;
  const Int_t MAXCOMP=400;


class ana10 : public h10{
public:
  //Hystogrames created are stored in map containers!
  //Each ana10 instance created with its own unique name aname.
  //In the Loop cycle, one can create TH1D or TH2D hystorgma via call to 
  //Addh1(hname...) Addh22, providing hyst name hname, and parameters.
  //Created hystos have name=aname+hname, so it is always unique in ROOT memory.
  //Hystogram hname is automatically created in every ana10 instance, 
  //and user provided access to it  through Geth1(hname), Geth2(hname).
  //Thus it is doesn't metter how many background components are in the analysis,
  //and there is no need to invent new names, since its done automatically.
  TH1Dmap d1map;//1Dim hystos
  TH2Dmap d2map;//2Dim hystos
  Timesmap runtimes;

  TH1D* blankh1;
  TH2D* blankh2;

  string name;

  //Selected events parameters for RooFit 
  TTree *events;
  Double_t Esum_ana10;
  Double_t Ediff_ana10;
  Double_t Cos_ana10;
  Double_t w_ana10;
  TBranch *bEsum,*bEdiff,*bCos,*bw;

  Long64_t mcevga;
  Float_t MCscalefactor;// = mcevga / N decays in real experiment 
  Float_t activity;
  Float_t eactivity;
  Float_t norm;
  Float_t enorm;
  Float_t totaltime;

  //References for all hystogramm instances created 
  //----- no up to the moment

//Manager for calculating pmstatuses


  ana10();
  ana10(const char *prefix);
  ~ana10();
  void Loop(Int_t);
  Bool_t Addh1(const char *,const char *,Int_t,Float_t,Float_t);
  Bool_t Addh1(const char *hnamec, TH1D * source);
  TH1D * Geth1(const char *);
  Bool_t Checkh1(const char * hname);
  Bool_t Addh2(const char *,const char *,Int_t,Float_t,Float_t,Int_t,Float_t, Float_t);
  TH2D * Geth2(const char *);
  void Write(TFile *);
  void Read(TFile *);
  void Read(TFile *,string);
  void Getname(char *);
  const char * Getname();
  void Geth1list(TObjArray *);
  void Geth2list(TObjArray *);
  Float_t GetEventWeight();

  //This stuff is responsible for runlist selection  
  Bool_t RunList(Int_t,Int_t &);
  void SetPeriod(ana::period,Bool_t);
  void SetRuns(ana::run,Bool_t);
  void AddRunTime(Int_t, Float_t);
  Float_t CalculateTime();
  Float_t MCTime();
  Float_t ExpTime();
  
  Bool_t swactivity_use;//use individual activities per sector/gg plane
  Bool_t swactivity_ini;// activities initialised
  Float_t swweight[20][19];// activities per sector/gg plane

  Bool_t exactivity_use;//use individual activities per sector/gg plane
  Bool_t exactivity_ini;// activities initialised
  Bool_t noextruemc_warned;// warning that no true MC info available already issued
  Float_t exweight[3];// activities per sector/gg plane

  //work with individual activities
  Float_t getActivityWeight();

  Float_t getswWeight();
  void initswActivity();
  Float_t getswMeanActivity();

  Float_t getexWeight();
  void initexActivity();
  Float_t getexMeanActivity();

  Float_t getRnActivityRunList();
  
  Int_t GetRunDate(int);
  Int_t GetYear(int);
  Int_t GetMonth(int);
  Int_t GetDay(int);
  Int_t Days_Nemo_Era(int);
  
  
protected:
  // Define period to be used in the analysis
  Bool_t P1;// Phase 1 period from Feb 2003 to Oct 2004
  Bool_t P2a; // Phase 2 period a) from Oct 2004 to 22 may 2006
  Bool_t P2b; // Phase 2 period b) from 22 may 2006 to 31 Dec 2007
  Bool_t PBA; // BA test list from Yuriy
  Bool_t P3; // new 2007 data
  Bool_t P1_0nu; // Phase 1 for 0nu analysis
  Bool_t P2_0nu; // Phase 2 for 0nu analysis

  //Define Runs to be used in the analysis
  Bool_t set0,set1,set2,set3,set4,extended,standard,bb0nu;
  

};


class bana10: public ana10{
public:

  //related backgrounds (belong to one chain), e.g. Ac228 for Tl208,
  //and e.g. related[i]->activity=this->activity*relfactor[i]
  bana10 * relative[10];
  Float_t  relfactor[10];
  Int_t nrelative;
  Bool_t dependent;
  Bool_t scaled;
  bana10 * parent;


  RooAbsData *Esum_set;
  RooDataHist *Hist_set;
  RooKeysPdf *Esum_pdf;
  RooAbsPdf *E2D_pdf;
  RooAbsPdf *d1pdf;



  bana10(const char *pref,Long64_t nmc,Float_t a);
  ~bana10();
  void setnorm(Float_t t);// set norm factor for t sec.
  void setnormR(Float_t t);// set norm factor for t sec.
  void setnormalization(Float_t,Float_t);
  int Print(const char *hname);
  int Print();// Prints analysis details to std::cout
  void Addrelative(bana10 *,Float_t );
  Double_t GetExpected(const char *,Int_t lx=0,Int_t ux=0);
  Double_t GetErrorExpected(const char *,Int_t lx=0 ,Int_t ux=0);
  void scale(Float_t s);

  //RooFit usage for ana10
  Int_t ConstructPdf(RooRealVar&,Bool_t key=true);
  Int_t ConstructSet(RooRealVar&,Bool_t key=true);
  Int_t Construct2DEPdf(RooRealVar&,RooRealVar&,Bool_t key=true);


};

void Drawall_d1(const char * hname, ana10 *data,bana10 **bgr, Int_t nbgr);

//void Drawall_d1(const char * hname, ana10 *data,bana10 **bgr, Int_t nbgr,bana10 **,Int_t,Bool_t kLegend = true,Bool_t kSum=false,Bool_t cumlative=false,Double_t blow=0., Double_t bup=0.,Int_t rebin=1);

void Drawsn_d1(const char * hname, ana10 *data,bana10 **bgr, Int_t nbgr,bana10 **,Int_t,Bool_t kLegend = true,Bool_t kSum=false,Bool_t cumlative=false,Double_t blow=0., Double_t bup=0.,Int_t rebin=1);

//void Drawconf_d1(const char * hname, ana10 *data,bana10 **bgr, Int_t nbgr,bana10 **,Int_t,Bool_t kLegend = true,Bool_t kSum=false,Bool_t cumlative=false,Double_t blow=0., Double_t bup=0.,Int_t rebin=1,Bool_t draw_on=false,bana10 **lim=NULL, Int_t nlim=0);
void Drawconf_d1(const char * hname, ana10 *data,bana10 **bgr, Int_t nbgr,bana10 **,Int_t,Bool_t kLegend = true,Bool_t kSum=false,Bool_t cumlative=false,Double_t blow=0., Double_t bup=0.,Int_t rebin=1,Int_t fillstyle=1001,Bool_t draw_on=false,bana10 **lim=NULL, Int_t nlim=0);
void Drawconf_d1_Ch(const char * hname, bana10 *data[],bana10 *bgr[][MAXCOMP], Int_t nbgr[], bana10 *sig[][MAXCOMP], Int_t nsig[],string namechnl[],Int_t nchnls,Bool_t kLegend=true, Bool_t kSum=false,Bool_t cumulative=false, Double_t blow=0,Double_t bup=0.,Int_t rebin=1,Int_t fillstyle=1001,Bool_t draw_0n=false,bana10 **lim=NULL,Int_t nlim=0);

void Drawall_d1(const char * hname,ana10 *data,bana10 **bgr,Int_t nbgr,bana10 **,Int_t,Int_t kLegend=1,Bool_t kSum=false,Bool_t cumlative=false,Double_t blow=0.,Double_t bup=0.,Int_t rebin=1,Int_t fillstyle=1001,Int_t printnum=10,Double_t yscale=1.0);
void Drawall_d1_Ch(const char * hname, bana10 *data[],bana10 *bgr[][MAXCOMP], Int_t nbgr[], bana10 *sig[][MAXCOMP], Int_t nsig[],string namechnl[],Int_t nchnls,Bool_t kLegend=true, Bool_t kSum=false,Bool_t cumulative=false, Double_t blow=0.,Double_t bup=0.,Int_t rebin=1,Int_t fillstyle=1001);

//void Drawconf_d1(const char * hname, ana10 *data,bana10 **bgr, Int_t nbgr,bana10 **,Int_t,Bool_t kLegend = true,Bool_t kSum=false,Bool_t cumlative=false,Double_t blow=0., Double_t bup=0.,Int_t rebin=1);
void Drawdiff_d1(const char * hname, ana10 *data,bana10 **bgr, Int_t nbgr, bana10 ** sig, Int_t nsig,Bool_t kLegend=false);
void Drawall_d2x(const char * hname, ana10 *data,bana10 **bgr, Int_t nbgr, bana10 ** sig, Int_t nsig,Bool_t kLegend=false,Int_t rebin=1);
void Drawall_d2y(const char * hname, ana10 *data,bana10 **bgr, Int_t nbgr, bana10 ** sig, Int_t nsig,Bool_t kLegend=false,Int_t rebin =1);

void Draw_Ch_merge(const char * hname, bana10 *data[],bana10 *bgr[][MAXCOMP], Int_t nbgr[], bana10 *sig[][MAXCOMP], Int_t nsig[],string namechnl[],Int_t nchnls,bana10 *& newdata,bana10* newbgr[MAXCOMP],Int_t& newnbgr, bana10 * newsig[MAXCOMP], Int_t& newnsig);


// added by matt kauer
///////////////////////////////////////////////////////
void DrawBreakdown(const char *hname,ana10 *data,bana10 **bgr,Int_t nbgr,bana10 **sig,Int_t nsig,Int_t *hdraw,vector<string> foilbkgs,vector<string> vetobkgs,Double_t xmin=0,Double_t xmax=0,Int_t rebin=1,string s_type="",Double_t yscale=1,string isotope="");
void DrawBkgConf(const char *hname,ana10 *data,bana10 **bgr,Int_t nbgr,bana10 **sig,Int_t nsig,Double_t xmin=0,Double_t xmax=0,Int_t rebin=1,string s_type="",Double_t yscale=1,Int_t notused=0,Int_t order=0,string isotope="");
void DrawBkgAdd(const char *hname,ana10 *data,bana10 **bgr,Int_t nbgr,bana10 **sig,Int_t nsig,Double_t xmin=0,Double_t xmax=0,Int_t rebin=1,string s_type="",string isotope="");
void DrawBkgSub(const char *hname,ana10 *data,bana10 **bgr,Int_t nbgr,bana10 **sig,Int_t nsig,Double_t xmin=0,Double_t xmax=0,Int_t rebin=1,string s_type="",string isotope="");
void DrawLim(const char *hname,ana10 *data,bana10 **bgr,Int_t nbgr,bana10 **sig,Int_t nsig,bana10 **lim,Int_t nlim,Double_t xmin=0,Double_t xmax=0,Int_t rebin=1,Double_t limEvents=-1,string s_type="",string isotope="",string decaymode="");
string printHalfStatSyst(Float_t halflife, Float_t staterr, Float_t systerr, Int_t precision);
void createCollieHists(const char *hname,ana10 *data,bana10 **bgr,Int_t nbgr,bana10 **sig,Int_t nsig,bana10 **lim,Int_t nlim,Int_t rebin=1,Double_t limEvents=-1,string decaymode="");
///////////////////////////////////////////////////////

Int_t FitSignal(const char * name, bana10 *data,bana10 **bgr,Int_t nbgr,bana10 ** sig,Int_t nsig,Int_t,Int_t);
Int_t FitSignalCh(const char * name, bana10 *data[],bana10 *bgr[][MAXCOMP],Int_t nbgr[],bana10 * sig[][MAXCOMP],Int_t nsig[],Int_t nchnls,Int_t lx, Int_t ux);
Int_t FitSignal(const char * name, bana10 *data,bana10 **bgr,Int_t nbgr,bana10 ** sig,Int_t nsig, Float_t blow, Float_t bup, Int_t rebin=1);
Int_t LFitSignal(const char * name, bana10 *data,bana10 **bgr,Int_t nbgr,bana10 ** sig,Int_t nsig, Float_t blow, Float_t bup, Int_t rebin=1, Double_t scalebkgs=1.0);
Int_t LFitSignalCh(const char * name, bana10 *data[],bana10 *bgr[][MAXCOMP],Int_t nbgr[],bana10 * sig[][MAXCOMP],Int_t nsig[], Int_t nchnls, Float_t blow, Float_t bup, Int_t rebin);


//Int_t HelenLimit(const char * name, bana10 *data,bana10 **bgr,Int_t nbgr,bana10 ** sig,Int_t nsig,bana10 ** lim, Int_t nlim,Int_t lx, Int_t ux);
//Int_t HelenLimitf(const char * name, bana10 *data,bana10 **bgr,Int_t nbgr,bana10 ** sig,Int_t nsig,bana10 ** lim, Int_t nlim,Float_t lx, Float_t ux);
Double_t HelenLimit(const char * name, bana10 *data,bana10 **bgr,Int_t nbgr,bana10 ** sig,Int_t nsig,bana10 ** lim,Int_t nlim,Int_t lx,Int_t ux);
Double_t HelenLimitf(const char * name, bana10 *data,bana10 **bgr,Int_t nbgr,bana10 ** sig,Int_t nsig,bana10 ** lim,Int_t nlim,Float_t lx,Float_t ux);
//Int_t MCLimit(const char * name, bana10 *data,bana10 **bgr,Int_t nbgr,bana10 ** sig,Int_t nsig,bana10 ** lim, Int_t nlim,Int_t rebin=1,Double_t es=0.,Double_t eb=0.);
Double_t MCLimit(const char * name, bana10 *data,bana10 **bgr,Int_t nbgr,bana10 ** sig,Int_t nsig,bana10 ** lim, Int_t nlim,Int_t rebin=1,Double_t es=0.,Double_t eb=0.);
Int_t MCLimit2d(const char * name, bana10 *data,bana10 **bgr,Int_t nbgr,bana10 ** sig,Int_t nsig,bana10 ** lim, Int_t nlim,Int_t rebinx=1,Int_t rebiny=1,Double_t es=0.,Double_t eb=0);
Int_t MCLimitch(const char * name, bana10 *data[],bana10 *bgr[][MAXCOMP],Int_t nbgr[],bana10 * sig[][MAXCOMP],Int_t nsig[],bana10 * lim[][MAXCOMP], Int_t nlim[],Int_t nchnls,Int_t rebin,Double_t es=0.,Double_t eb=0);
Int_t MCLimit2dch(const char * name, bana10 *data[],bana10 *bgr[][MAXCOMP],Int_t nbgr[],bana10 * sig[][MAXCOMP],Int_t nsig[],bana10 * lim[][MAXCOMP], Int_t nlim[],Int_t nchnls,Int_t rebinx=1, Int_t rebiny=1,Double_t es=0.,Double_t eb=0);
Int_t MajoronLimit(const char * name, bana10 *data,bana10 **bgr,Int_t nbgr,bana10 ** sig,Int_t nsig,bana10 ** lim, Int_t nlim,Int_t rebin =1,Double_t es=0.,Double_t eb=0);


Int_t LExclusion(const char * name, bana10 **bgr,Int_t nbgr,bana10 ** sig,Int_t nsig, bana10 **lim, Int_t nlim, Float_t blow=0, Float_t bup=0, Int_t rebin=1);
Int_t LExclusion_pessimistic(const char * name, bana10 **bgr,Int_t nbgr,bana10 ** sig,Int_t nsig, bana10 **lim, Int_t nlim, Float_t blow=0, Float_t bup=0, Int_t rebin=1,Int_t SN=1);



Int_t ReadSavedRoot(const char * fname,bana10 *&data,bana10 ** bgr,Int_t &nbgr,Float_t &exptime);//read all components saved
Int_t ReadSavedRoot(const char * fname,bana10 *&data,bana10 ** bgr,Int_t &nbgr,bana10 **, Int_t,Float_t &exptime);//read all components saved
Int_t ReadSavedRoot(const char * fname,const char * name,Float_t activity,bana10 *&ana);
//Int_t ReadSavedRootCh(const char * fname,const char * name,string nchnl,Float_t activity,bana10 *&ana, Float_t forcemcevga=0);
Int_t ReadSavedRootCh(const char * fname,const char * name,string nchnl,Float_t activity,bana10 *&ana, Float_t forcemcevga=0,Float_t scalebkg=-1);
Int_t ReadSlim(const char * fname,const char *name,Float_t activity,Long64_t mcevga,bana10 *&ana,Int_t version=7);
Int_t SaveRoot(const char *fname,bana10 *data,bana10** bgr,Int_t nbgr,Bool_t recreate=true);

//Main user provided subroutine to display/print analysis results.
Int_t Urun(bana10 *data[],bana10 *bgr[][MAXCOMP],bana10 * sig[][MAXCOMP],bana10 * lim[][MAXCOMP],Int_t nbgr[],Int_t nsig[],Int_t nlim[],Int_t nchnls, string namechnl[]);
//Int_t Urun(bana10 *data[],bana10 *bgr[][MAXCOMP],bana10 * sig[][MAXCOMP],Int_t nbgr[],Int_t nsig[],Int_t nchnls, string namechnl[]);

void ROOT_settings();

string printenumber(Float_t number, Float_t error, Int_t precision);


TH1 * GetCumulativeDown(TH1 * src,Int_t nx);
TH1 * GetCumulativeDown(TH1 * src,Double_t);
Int_t CopyShape1D(TH1 *src, TH1 *tgt);

//********************** Likelihood limit service functions*******************
Double_t X(RooAbsPdf &sbpdf, RooAbsPdf &bpdf,RooDataSet &data);
Double_t LLHCL(RooAbsPdf *bpdf,RooAbsPdf *sbpdf,RooDataSet *expdata,const Double_t precision=0.02);
Double_t LLHlimit(TList* bpdfs, TList* bnorms,RooAbsPdf *spdf,RooDataSet *expdata,const Double_t CL=0.9, const Double_t eps=0.1,const Double_t guesslimit=100);
Double_t RooLLHLimit(const char * name, bana10 *data,bana10 **bgr,Int_t nbgr,bana10 ** sig,Int_t nsig,bana10 ** lim, Int_t nlim,Double_t xmin=2.551);


int compare_doubles(const void *a,const void * b);

Bool_t stest(string s1, string s2);


#endif // ANA10_H
