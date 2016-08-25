#ifndef PMSTATUS_H
#define PMSTATUS_H


//#include "N3Db.h"
//#include "N3DbMgr.h"
#include <vector>
#include <string>
#include <TROOT.h>
#include <TH1.h>
#include <TCanvas.h>
#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>



#define RUNMIN 1800
#define RUNMAX 7000
#define PMMAX  1940


//define new PM statuses for ana10 analysis
#define BAD_ALPHA -10
#define BAD_ALPHA_ERR -100
#define BAD_RESOL -1000
#define BAD_RESOL_ERR -10000



Bool_t GoodPMStatus(Int_t);
Bool_t CheckStatus(Int_t , Int_t);
Bool_t Check_LEC_jump(Int_t );

class Pmstatus { 
 public:
  Pmstatus();
  ~Pmstatus();
  Int_t get(Int_t run,Int_t pmnum);
  Int_t get(Int_t run,Int_t sec, Int_t io,Int_t fcll, Int_t blk);
  void draw_calib_stat();

  Float_t get_pm_resol(Int_t run, Int_t pmnum);
  Float_t get_pm_resol(Int_t run, Int_t sec, Int_t io,Int_t fcll, Int_t blk);
  Float_t get_pm_resol_err(Int_t run, Int_t pmnum);
  Float_t get_pm_resol_err(Int_t run, Int_t sec, Int_t io,Int_t fcll, Int_t blk);

  Float_t get_pm_alpha(Int_t run, Int_t pmnum);
  Float_t get_pm_alpha(Int_t run, Int_t sec, Int_t io,Int_t fcll, Int_t blk);
  Float_t get_pm_alpha_err(Int_t run, Int_t pmnum);
  Float_t get_pm_alpha_err(Int_t run, Int_t sec, Int_t io,Int_t fcll, Int_t blk);

  Float_t get_pm_ltce(Int_t run, Int_t pmnum);
  Float_t get_pm_ltce(Int_t run, Int_t sec, Int_t io,Int_t fcll, Int_t blk);


 protected:
  Int_t statusmap[RUNMAX-RUNMIN][PMMAX+1];
  Float_t tse[RUNMAX-RUNMIN][PMMAX+1];
  Float_t ltce[RUNMAX-RUNMIN][PMMAX+1];
  Float_t pm_resol[RUNMAX-RUNMIN][PMMAX+1];
  Float_t pm_resol_err[RUNMAX-RUNMIN][PMMAX+1];
  Float_t pm_alpha[RUNMAX-RUNMIN][PMMAX+1];
  Float_t pm_alpha_err[RUNMAX-RUNMIN][PMMAX+1];
  Int_t runcache[RUNMAX-RUNMIN];

  //some statistic on calibration parameters
  
  TH1D * alpha,*alpha_err,*resol, *resol_err, *tdcf; 

  FILE *fp;
  
};

#endif  //PMSTATUS_H
