//////////////////////////////////////////////////////////
//   This class has been automatically generated 
//     (Thu Jun 17 10:28:30 2004 by ROOT version3.02/07)
//   from TTree h10/MYNTUPLE
//   found on file: 2285-2290.root
//////////////////////////////////////////////////////////


#ifndef __h10_h
#define __h10_h
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TMath.h>   // added for ROOT compatibility
#include <TRandom3.h>
#include "scintillator.h"
//#include "pmstatus.h"

#include <iostream>

#define SC_TDC_CH 0.053 // scintillator TDC channle in ns

//extern Pmstatus *pmstatus;


extern char Snuclide[][16];
extern char Sfoil[][16];

Float_t A2HalfLife(Float_t A, Int_t* isrc, Int_t nsrc);


class h10 {
   public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain
   Long_t timetot;
   char  rfile[250];
   TRandom3 grndm;
   Int_t nemosv;



//Declaration of leaves types

   Int_t           Nsc;
   Int_t           Ngg;
   Float_t         Sc[2000][12];   //[Nsc]
   /**--
*--   SCINTILLATOR DATA:
*--   NSC         : NUMBER OF SCINTILLATORS HIT (ISC RUNS FROM 1 TO NSC)
*--   SC(1,ISC)   : PM STATUS
*--   SC(2,ISC)   : SECTOR NUMBER (0-19 AS IN KDIGI)
*--   SC(3,ISC)   : IOBT FLAG (0- 3 AS IN KDIGI)
*--   SC(4,ISC)   : FCLL    (0- 3 AS IN KDIGI)
*--   SC(5,ISC)   : BLOC NUMBER (0-16 AS IN KDIGI)
*--   SC(6,ISC)   : ADC CONTENT OF THE SCINTILLATOR HIT
*--   SC(7,ISC)   : TDC CONTENT OF THE SCINTILLATOR HIT
*--   SC(8,ISC)   : FLAG    SET TO 1 IF SIGNAL ABOVE HIGH THRESHOLD
*--   SC(9,ISC)   : ENERGY DEPOSIT (GeV)
*--   SC(10,ISC)  : X COORDINATE OF THE FRONT FACE OF THE SCINTILLATOR HIT
*--                 (M.R.S.)
*--   SC(11,ISC)  : Y COORDINATE OF THE FRONT FACE OF THE SCINTILLATOR HIT
*--                 (M.R.S.)
*--   SC(12,ISC)  : Z COORDINATE OF THE FRONT FACE OF THE SCINTILLATOR HIT
*--                 (M.R.S.)
**/

   Float_t         Gg[3000][15];   //[Ngg]

   /**--   GEIGER CELL DATA:
*--   NGG         : NUMBER OF CELLS HIT (IGG RUNS FROM 1 TO NGG)
*--   GG(1,IGG)   : WSTATUS
*--   GG(2,IGG)   : SECTOR NUMBER (0-19 AS IN KDIGI)
*--   GG(3,IGG)   : IO FLAG (0- 1 AS IN KDIGI)
*--   GG(4,IGG)   : ABCD..I FLAG  (0- 8 AS IN KDIGI)
*--   GG(5,IGG)   : CELL NUMBER
*--   GG(6,IGG)   : FAST TDC CONTENT
*--   GG(7,IGG)   : SLOW TDC CONTENT
*--   GG(8,IGG)   : TDC CONTENT FOR CATHODE 1 (BOTTOM)
*--   GG(9,IGG)   : TDC CONTENT FOR CATHODE 1 (TOP)
*--   GG(10,IGG)  : X COORDINATE OF THE WIRE HIT (M.R.S.)
*--   GG(11,IGG)  : Y COORDINATE OF THE WIRE HIT (M.R.S.)
*--   GG(12,IGG)  : Z COORDINATE OF THE HIT (M.R.S.) -
*--                 = 0. IF THE Z COORDINATE IS NOT  MEASURED
*--   GG(13,IGG)  : ELECTRON DRIFT DISTANCE FOR FAST HITS -
*--                 TIME (ns) FOR DELAYED HITS
*--   GG(14,IGG)  : ERROR ON THE Z COORDINATE OF THE HIT (M.R.S.) -
*--                 >= 5000. IF THE Z COORDINATE IS NOT OR POORLY MEASURED
*--   GG(15,IGG)  : RECIPROCAL OF THE WEIGHT USED FOR CIRCLE FIT -
*--                 FIXED TO 1.cm FOR DELAYED HITS
**/

   Int_t           Nbr_tks;
   Char_t          Nbr_pts[10];   //[Nbr_tks]
   UChar_t         Ind_points[10][200];   //[Nbr_tks]
   Float_t         Xc[10];   //[Nbr_tks]
   Float_t         Yc[10];   //[Nbr_tks]
   Float_t         Zc[10];   //[Nbr_tks]
   Float_t         Radc[10];   //[Nbr_tks]
   Float_t         Hc[10];   //[Nbr_tks]
   Float_t         Ec[10];   //[Nbr_tks]
   Float_t         Dec[10];   //[Nbr_tks]
   Float_t         Qc[10];   //[Nbr_tks]
   Float_t         Prob_radc[10];   //[Nbr_tks]
   Float_t         Prob_hc[10];   //[Nbr_tks]
   Float_t         Prob_helix[10];   //[Nbr_tks]
   Float_t         X_foil[10];   //[Nbr_tks]
   Float_t         Y_foil[10];   //[Nbr_tks]
   Float_t         Z_foil[10];   //[Nbr_tks]
   Float_t         Cos_dir[10][3];   //[Nbr_tks]
   Int_t           Myievent;
   Int_t           Mynbr_tks;
   Float_t         X_scintil[10];   //[Mynbr_tks]
   Float_t         Y_scintil[10];   //[Mynbr_tks]
   Float_t         Z_scintil[10];   //[Mynbr_tks]
   Char_t          Ind_scintil[10];   //[Mynbr_tks] = IND_SCINT
   Char_t          Myimpact[10];   //[Mynbr_tks] = IMPACT
   Float_t         Xvert;// = X_VERT
   Float_t         Yvert;
   Float_t         Zvert;
   /**--
*--   NBR_TKS         : NUMBER OF RECONSTRUCTED TRACKS
*--   NBR_PTS(IT)     : NUMBER OF ASSOCIATED HITS ON TRACK IT
*--   IND_POINTS(NP,IT) IGG VALUE FOR HIT NP ON TRACK IT
*--   XC(IT)          : X COORDINATE OF THE CENTER OF CURVATURE OF TRACK IT
*--   YC(IT)          : Y COORDINATE OF THE CENTER OF CURVATURE OF TRACK IT
*--   ZC(IT)          : Z COORDINATE OF THE CENTER OF CURVATURE OF TRACK IT
*--   RADC(IT)        : RADIUS OF CURVATURE OF TRACK IT
*--   HC(IT)          : HELIX PITCH OF TRACK IT
*--   EC(IT)          : ENERGY OF TRACK IT CALCULATED FROM THE CURVATURE
*--   QC(IT)          : CHARGE OF TRACK IT CALCULATED FROM THE CURVATURE
*--   PROB_RADC(IT)   : PROBABILITY OF CIRCLE FIT FOR TRACK IT
*--                                               IN THE HORIZONTAL PLANE
*--   PROB_HC(IT)     : PROBABILITY OF LINE FIT FOR TRACK IT
*--                                               IN THE R*THETA PLANE
*--   PROB_HELIX(IT)  : PROBABILITY OF THE HELIX FIT FOR TRACK IT
*--   X_FOIL(IT)      : X COORDINATE OF THE INTERSECTION OF TRACK IT
*--                                     WITH THE FOIL CYLINDER (FROM TRACKING)
*--   Y_FOIL(IT)      : Y COORDINATE OF THE INTERSECTION OF TRACK IT
*--                                     WITH THE FOIL CYLINDER (FROM TRACKING)
*--   Z_FOIL(IT)      : Z COORDINATE OF THE INTERSECTION OF TRACK IT
*--                                     WITH THE FOIL CYLINDER (FROM TRACKING)
*--   COS_DIR(3,IT)   : DIRECTION COSINES OF TRACK IT
*--                     (FROM TRACKING OR FROM TRACK FIT AND
*--                                                 INTERSECTION WITH THE FOIL)
*--   IND_SCINT(IT)   : ISC VALUE OF THE SCINTILLATOR ASSOCIATED TO TRACK IT
*--   IMPACT(IT)      : IMPACT REGION OF THE TRACK ONTO THE SCINTILLATOR BLOCK
*--   X_VERT          : X COORDINATE OF THE EVENT VERTEX
*--   Y_VERT          : Y COORDINATE OF THE EVENT VERTEX
*--   Z_VERT          : Z COORDINATE OF THE EVENT VERTEX
**/
   
   Int_t             run;  // event RUN number for data, for MC = -run which conditions used for detector simulation
   Int_t            date;  // event date
   Int_t            time;  // duration of run in s RAW data, number MC generated for MC
   Float_t   tau_sc_save;  // scintillator decay constant
   
 /** KINEMATICS SECTION
       True MC information about particles generated. Stored only if reconstructed with special flag for nemor. Majority of reco files do not have it due to filesize issue.

c     GENERATED DATA:
c     NVNTU  NUMBER OF VERTICES
c     XVNTU  X COORDINATE OF THE VERTEX
c     YVNTU  Y COORDINATE OF THE VERTEX
c     ZVNTU  Z COORDINATE OF THE VERTEX
c     TOFVNTU  TIME OF CREATION OF VERTEX
c     NTNTU  NUMBER OF TRACKS
c     PXNTU  PX OF THE TRACK
c     PYNTU  PY OF THE TRACK
c     PZNTU  PZ OF THE TRACK
c     TOFTNTU TIME OF CREATION OF TRACK
c     IVNTU  VERTEX OF ATTACHMENT
c     IPNTU  GEANT PARTICLE TYPE

 **/

   Bool_t          trueMC_flag;// indicates if true MC present in the file
   Int_t           Nvntu;
   Float_t         Xvntu[7];   //[Nvntu]
   Float_t         Yvntu[7];   //[Nvntu]
   Float_t         Zvntu[7];   //[Nvntu]
   Float_t         Tofvntu[7]; //[Nvntu]
   Int_t           Ntntu;
   Float_t         Pxntu[31];  //[Ntntu]
   Float_t         Pyntu[31];  //[Ntntu]
   Float_t         Pzntu[31];  //[Ntntu]
   Float_t         Toftntu[31];//[Ntntu]
   UChar_t         Ivntu[31];  //[Ntntu]
   UChar_t         Ipntu[31];  //[Ntntu]


 //-------------------------------------------

   //Useful variables
   Int_t Ndgg;// Number of delayed hits
   Int_t dgg[1000];// List of delayed hits found after checkdelayed() call

//List of branches
   TBranch        *b_tau_sc_save;   //!
   TBranch        *b_Nsc;   //!
   TBranch        *b_Ngg;   //!
   TBranch        *b_Sc;   //!
   TBranch        *b_Gg;   //!
   TBranch        *b_Nbr_tks;   //!
   TBranch        *b_Nbr_pts;   //!
   TBranch        *b_Ind_points;   //!
   TBranch        *b_Xc;   //!
   TBranch        *b_Yc;   //!
   TBranch        *b_Zc;   //!
   TBranch        *b_Radc;   //!
   TBranch        *b_Hc;   //!
   TBranch        *b_Ec;   //!
   TBranch        *b_Dec;   //!
   TBranch        *b_Qc;   //!
   TBranch        *b_Prob_radc;   //!
   TBranch        *b_Prob_hc;   //!
   TBranch        *b_Prob_helix;   //!
   TBranch        *b_X_foil;   //!
   TBranch        *b_Y_foil;   //!
   TBranch        *b_Z_foil;   //!
   TBranch        *b_Cos_dir;   //!
   TBranch        *b_Myievent;   //!
   TBranch        *b_Mynbr_tks;   //!
   TBranch        *b_X_scintil;   //!
   TBranch        *b_Y_scintil;   //!
   TBranch        *b_Z_scintil;   //!
   TBranch        *b_Ind_scintil;   //!
   TBranch        *b_Myimpact;   //!
   TBranch        *b_Xxvert;   //!
   TBranch        *b_Yyvert;   //!
   TBranch        *b_Zzvert;   //!
   TBranch        *b_run;
   TBranch        *b_date;
   TBranch        *b_time;
   TBranch        *b_Nvntu;   //!
   TBranch        *b_Xvntu;   //!
   TBranch        *b_Yvntu;   //!
   TBranch        *b_Zvntu;   //!
   TBranch        *b_Tofvntu;   //!
   TBranch        *b_Ntntu;   //!
   TBranch        *b_Pxntu;   //!
   TBranch        *b_Pyntu;   //!
   TBranch        *b_Pzntu;   //!
   TBranch        *b_Toftntu;   //!
   TBranch        *b_Ivntu;   //!
   TBranch        *b_Ipntu;   //!

   h10(TTree *tree=0);
   ~h10();

// Functions for alpha particle analisis and Radon rejection.

   Int_t NAfasthits();
   Int_t NA_fast_hits(Float_t max_dst);
   Int_t fastGG_other_side(Float_t max_dst);

   Int_t mywstatus(Int_t igg);
   Int_t checkdelayed();
   Int_t checksingledhit(Float_t &atime);
   Int_t checksingledhit(Float_t &atime,Float_t &dist);
   Int_t checksingledhit2(Float_t &atime,Float_t &dist);
   Int_t checkggclose(Int_t igg1, Int_t igg2,Float_t dxymax,Float_t dzmax);
   Int_t checkgroupdhit(Float_t &atime);
   Int_t checkgroupdhit(Float_t &atime, Float_t &adtime,Float_t& alength);
   Int_t checkgroupdhit_vera(Float_t &atime, Float_t &adtime,Float_t& alength);
  


   // TOF functions
   Float_t getsctime(Int_t isc,Float_t &sigma);
   Float_t getsctimebb(Int_t isc,Float_t &sigma);
   void   TOF_bsigle(Float_t &tf, Float_t &sf,Int_t tk1);
   Float_t TOFbg_ext(Int_t iscg, Float_t &tf, Float_t &sf);
   Float_t TOFbg_int(Int_t iscg, Float_t &tf, Float_t &sf);
   Float_t TOFbg_scat(Int_t iscg);
   Float_t TOF2b_int(Float_t &tf, Float_t &sf, Int_t tk1=0,Int_t tk2=1);
   Float_t TOF2b_ext(Int_t tk1=0,Int_t tk2=1);

   Float_t TOF2b_gas_int(Float_t &tf, Float_t &sf, Float_t *vtx,Int_t tk1=0,Int_t tk2=1);
   Float_t TOF2b_gas_ext(Float_t *vtx,Int_t tk1=0,Int_t tk2=1);

   Float_t TOF2g_int(Int_t isc1,Int_t isc2,Float_t X,Float_t Y,Float_t Z,Float_t tfoil,Float_t dtfoil);
   Float_t TOF1g_int(Int_t isc1, Float_t X,Float_t Y,Float_t Z,Float_t tfoil,Float_t dtfoil);
   Float_t TOF1g_ext(Int_t isc1, Float_t X,Float_t Y,Float_t Z,Float_t tfoil,Float_t dtfoil);
   Float_t TOF2g_scat(Int_t isc1,Int_t isc2,Float_t X,Float_t Y,Float_t Z,Float_t tfoil,Float_t dtfoil);
   Float_t TOF2g_ext(Int_t isc1,Int_t isc2,Float_t X,Float_t Y,Float_t Z,Float_t tfoil,Float_t dtfoil);
   Int_t twogamma(Float_t X, Float_t Y,Float_t Z,Float_t tfoil, Float_t dfoil, Float_t &Eg1,Float_t &Eg2);

   Int_t twogammav2(Float_t X, Float_t Y,Float_t Z,Float_t tfoil, Float_t dfoil, Float_t &Egamma1,Float_t &Egamma2,Float_t &);
   Int_t onegammav2(Float_t X, Float_t Y,Float_t Z,Float_t tfoil, Float_t dfoil, Float_t &Egamma1,Float_t &critery);
   Float_t chi2_tof1g(Int_t *scint, Int_t nscint,Float_t x,Float_t y,Float_t z, Float_t tfoil, Float_t dfoil);
   Float_t chi2_tof1g0(Int_t *scint,Int_t nscint, Float_t x,Float_t y,Float_t z, Float_t tfoil, Float_t dfoil);



   Float_t gethelixl(Int_t tr,Float_t &dlen);
   Float_t gethelixl(Int_t tr,Float_t *vtx,Float_t &dlen);

   Int_t  Cut(Int_t entry);
   Int_t  Cut(Int_t *results);
   Int_t  GetEntry(Int_t entry);
   Int_t  LoadTree(Int_t entry);
   void   Init(TTree *tree);
   void   Loop(Long64_t entry);
   void   testmodule();
   Float_t gethypoth(Int_t track, Int_t scint, Int_t flag, Int_t ieflag, Float_t* tf);
   Float_t getcurveylength(Int_t trackindex, Int_t scintindex, Float_t& dlen);
   Bool_t Notify();
   Float_t timeofflight(Float_t& sigmat,  Float_t& tfoil, Int_t object, Int_t scintobject, Int_t flag);

   Float_t getsinglehyp(Float_t t0, Float_t dt0, Float_t vertx, Float_t verty, Float_t vertz, Float_t sigl1, Int_t object, Int_t flag, Int_t ieflag, Float_t&, Float_t ecluster); 
   Float_t getstraightlengthv(Int_t scintind, Float_t X,Float_t Y,Float_t Z, Float_t& dlen);
   Float_t getstraightlength(Int_t scintind, Float_t X,Float_t Y,Float_t Z, Float_t& dlen, Float_t* pvec);
   Bool_t checkfirstlayer(Int_t trackindex, Int_t iflag);
   Bool_t checklayer(Int_t trackindex, Int_t row);

   Float_t getfoiltime(Float_t& sigmat, Float_t& len, Int_t object);
   Bool_t amIalone(Int_t scintNum1, Int_t scintNum2, Int_t flag);
   void dgammatogamma(Int_t scint1, Int_t scint2, Float_t& tof1, Float_t& dl, Float_t& dt);
   Int_t getsector(Int_t trackindex);
   Float_t getscatterhyp(Int_t scint1, Int_t scint2, Float_t& chi);
   void   Show(Int_t entry = -1);
   void assign();
   Float_t Calcsec(Float_t xvert, Float_t yvert);
   Int_t Decodesrc(Float_t xvert,Float_t yvert, Float_t zvert);
   void testme(Int_t eventnumber);
   Bool_t Preselect(Int_t eventnumber, Int_t* john);
   Bool_t preselect_1tk(Int_t eventnumber, Int_t* john);   
   Bool_t preselect_1e1g(Int_t eventnumber, Int_t* john);
   Bool_t preselect_1e(int eventnumber, Int_t* results);
   void Init_scintillator(scintillator* sci, Bool_t* isol, Int_t iflag);
   void addclusters(scintillator* sci, Int_t cluster1, Int_t cluster2);
   Int_t getalphas(Int_t* alparray);
   Bool_t Preselect_ntk(int eventnumber, Int_t* results, scintillator* sci, Int_t ntracks);
   Float_t find_cos_eg(Int_t gcluster,Int_t track,scintillator  sci);

   Float_t helix_cylinder_intersection(Int_t trk,Float_t t0, Float_t radius);
   Int_t vtx_volume(Int_t trk1, Int_t trk2, Float_t *vtx, Float_t& dvx, Float_t& dvz,Bool_t);

   Int_t FindWireLayer(Int_t trk);
   Int_t FindTrueSector(Int_t trk);
   Int_t FindEXBGposition(Int_t trk);
   Float_t getscfwhm_real(Int_t isc);

};

Bool_t CheckRunStatus(Int_t , Int_t);

#endif
#ifdef h10_cxx
#ifndef h10_cxx_done
#define h10_cxx_done

h10::h10(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.

  nemosv=7;

  grndm=TRandom3(0);

      //keep this
      Init(tree);
}

h10::~h10()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t h10::GetEntry(Int_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Int_t h10::LoadTree(Int_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Int_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->IsA() != TChain::Class()) return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void h10::Init(TTree *tree)
{
//   Set branch addresses
   if (tree == 0){
     //generate new tree
     fChain = new TTree("h10","blank h10 created");
     fChain->Branch("tau_sc_save",&Nsc,"tau_sc_save/F");
     fChain->Branch("Nsc",&Nsc,"Nsc/I");
     fChain->Branch("Ngg",&Ngg,"Ngg/I");
     fChain->Branch("Sc",Sc,"Sc[Nsc][12]/F");
     fChain->Branch("Gg",Gg,"Gg[Ngg][15]/F");
     fChain->Branch("Nbr_tks",&Nbr_tks,"Nbr_tks/I");
     fChain->Branch("Nbr_pts",Nbr_pts,"Nbr_pts[Nbr_tks]/B");
     fChain->Branch("Ind_points",Ind_points,"Ind_points[Nbr_tks][200]/b");
     fChain->Branch("Xc",Xc,"Xc[Nbr_tks]/F");
     fChain->Branch("Yc",Yc,"Yc[Nbr_tks]/F");
     fChain->Branch("Zc",Zc,"Zc[Nbr_tks]/F");
     fChain->Branch("Radc",Radc,"Radc[Nbr_tks]/F");
     fChain->Branch("Hc",Hc,"Hc[Nbr_tks]/F");
     fChain->Branch("Ec",Ec,"Ec[Nbr_tks]/F");
     fChain->Branch("Dec",Dec,"Dec[Nbr_tks]/F");
     fChain->Branch("Qc",Qc,"Qc[Nbr_tks]/F");
     fChain->Branch("Prob_radc",Prob_radc,"Prob_radc[Nbr_tks]/F");
     fChain->Branch("Prob_hc",Prob_hc,"Prob_hc[Nbr_tks]/F");
     fChain->Branch("Prob_helix",Prob_helix,"Prob_helix[Nbr_tks]/F");
     fChain->Branch("X_foil",X_foil,"X_foil[Nbr_tks]/F");
     fChain->Branch("Y_foil",Y_foil,"Y_foil[Nbr_tks]/F");
     fChain->Branch("Z_foil",Z_foil,"Z_foil[Nbr_tks]/F");
     fChain->Branch("Cos_dir",Cos_dir,"Cos_dir[Nbr_tks][3]/F");
     fChain->Branch("Myievent",&Myievent,"Myievent/I");
     fChain->Branch("Mynbr_tks",&Mynbr_tks,"Mynbr_tks/I");
     fChain->Branch("X_scintil",X_scintil,"X_scintil[Nbr_tks]/F");
     fChain->Branch("Y_scintil",Y_scintil,"Y_scintil[Nbr_tks]/F");
     fChain->Branch("Z_scintil",Z_scintil,"Z_scintil[Nbr_tks]/F");
     fChain->Branch("Ind_scintil",Ind_scintil,"Ind_scintil[Nbr_tks]/B");
     fChain->Branch("Myimpact",Myimpact,"Myimpact[Nbr_tks]/B");
     fChain->Branch("Xvert",&Xvert,"Xvert/F");
     fChain->Branch("Yvert",&Yvert,"Yvert/F");
     fChain->Branch("Zvert",&Zvert,"Zvert/F");
     fChain->Branch("run",&run,"run/I");
     fChain->Branch("date",&date,"date/I");
     fChain->Branch("time",&time,"time/I");
     return;
   } 
   fChain    = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("tau_sc_save",&tau_sc_save);
   fChain->SetBranchAddress("Nsc",&Nsc);
   fChain->SetBranchAddress("Ngg",&Ngg);
   fChain->SetBranchAddress("Sc",Sc);
   fChain->SetBranchAddress("Gg",Gg);
   fChain->SetBranchAddress("Nbr_tks",&Nbr_tks);
   fChain->SetBranchAddress("Nbr_pts",Nbr_pts);
   fChain->SetBranchAddress("Ind_points",Ind_points);
   fChain->SetBranchAddress("Xc",Xc);
   fChain->SetBranchAddress("Yc",Yc);
   fChain->SetBranchAddress("Zc",Zc);
   fChain->SetBranchAddress("Radc",Radc);
   fChain->SetBranchAddress("Hc",Hc);
   fChain->SetBranchAddress("Ec",Ec);
   fChain->SetBranchAddress("Dec",Dec);
   fChain->SetBranchAddress("Qc",Qc);
   fChain->SetBranchAddress("Prob_radc",Prob_radc);
   fChain->SetBranchAddress("Prob_hc",Prob_hc);
   fChain->SetBranchAddress("Prob_helix",Prob_helix);
   fChain->SetBranchAddress("X_foil",X_foil);
   fChain->SetBranchAddress("Y_foil",Y_foil);
   fChain->SetBranchAddress("Z_foil",Z_foil);
   fChain->SetBranchAddress("Cos_dir",Cos_dir);
   fChain->SetBranchAddress("Myievent",&Myievent);
   fChain->SetBranchAddress("Mynbr_tks",&Mynbr_tks);
   fChain->SetBranchAddress("X_scintil",X_scintil);
   fChain->SetBranchAddress("Y_scintil",Y_scintil);
   fChain->SetBranchAddress("Z_scintil",Z_scintil);
   fChain->SetBranchAddress("Ind_scintil",Ind_scintil);
   fChain->SetBranchAddress("Myimpact",Myimpact);
   fChain->SetBranchAddress("Xvert",&Xvert);
   fChain->SetBranchAddress("Yvert",&Yvert);
   fChain->SetBranchAddress("Zvert",&Zvert);
   fChain->SetBranchAddress("run",&run);
   fChain->SetBranchAddress("date",&date);
   fChain->SetBranchAddress("time",&time);
   
   //check if true MC information is there
   if(fChain->GetBranch("Nvntu")){
     fChain->SetBranchAddress("Nvntu",&Nvntu);
     fChain->SetBranchAddress("Xvntu",Xvntu);
     fChain->SetBranchAddress("Yvntu",Yvntu);
     fChain->SetBranchAddress("Zvntu",Zvntu);
     fChain->SetBranchAddress("Tofvntu",Tofvntu);
     fChain->SetBranchAddress("Ntntu",&Ntntu);
     fChain->SetBranchAddress("Pxntu",Pxntu);
     fChain->SetBranchAddress("Pyntu",Pyntu);
     fChain->SetBranchAddress("Pzntu",Pzntu);
     fChain->SetBranchAddress("Toftntu",Toftntu);
     fChain->SetBranchAddress("Ivntu",Ivntu);
     fChain->SetBranchAddress("Ipntu",Ipntu);
     trueMC_flag = true;
   }else{
     trueMC_flag = false;
   }
   

   Notify();
}

Bool_t h10::Notify()
{
   // Called when loading a new file.
   // Get branch pointers.
   b_tau_sc_save = fChain->GetBranch("tau_sc_save");
   b_Nsc = fChain->GetBranch("Nsc");
   b_Ngg = fChain->GetBranch("Ngg");
   b_Sc = fChain->GetBranch("Sc");
   b_Gg = fChain->GetBranch("Gg");
   b_Nbr_tks = fChain->GetBranch("Nbr_tks");
   b_Nbr_pts = fChain->GetBranch("Nbr_pts");
   b_Ind_points = fChain->GetBranch("Ind_points");
   b_Xc = fChain->GetBranch("Xc");
   b_Yc = fChain->GetBranch("Yc");
   b_Zc = fChain->GetBranch("Zc");
   b_Radc = fChain->GetBranch("Radc");
   b_Hc = fChain->GetBranch("Hc");
   b_Ec = fChain->GetBranch("Ec");
   b_Dec = fChain->GetBranch("Dec");
   b_Qc = fChain->GetBranch("Qc");
   b_Prob_radc = fChain->GetBranch("Prob_radc");
   b_Prob_hc = fChain->GetBranch("Prob_hc");
   b_Prob_helix = fChain->GetBranch("Prob_helix");
   b_X_foil = fChain->GetBranch("X_foil");
   b_Y_foil = fChain->GetBranch("Y_foil");
   b_Z_foil = fChain->GetBranch("Z_foil");
   b_Cos_dir = fChain->GetBranch("Cos_dir");
   b_Myievent = fChain->GetBranch("Myievent");
   b_Mynbr_tks = fChain->GetBranch("Mynbr_tks");
   b_X_scintil = fChain->GetBranch("X_scintil");
   b_Y_scintil = fChain->GetBranch("Y_scintil");
   b_Z_scintil = fChain->GetBranch("Z_scintil");
   b_Ind_scintil = fChain->GetBranch("Ind_scintil");
   b_Myimpact = fChain->GetBranch("Myimpact");
   b_Xxvert = fChain->GetBranch("Xxvert");
   b_Yyvert = fChain->GetBranch("Yyvert");
   b_Zzvert = fChain->GetBranch("Zzvert");
   b_run = fChain->GetBranch("run");
   b_date = fChain->GetBranch("date");
   b_time= fChain->GetBranch("time");
 
   b_Nvntu = fChain->GetBranch("Nvntu");
   b_Xvntu = fChain->GetBranch("Xvntu");
   b_Yvntu = fChain->GetBranch("Yvntu");
   b_Zvntu = fChain->GetBranch("Zvntu");
   b_Tofvntu = fChain->GetBranch("Tofvntu");
   b_Ntntu = fChain->GetBranch("Ntntu");
   b_Pxntu = fChain->GetBranch("Pxntu");
   b_Pyntu = fChain->GetBranch("Pyntu");
   b_Pzntu = fChain->GetBranch("Pzntu");
   b_Toftntu = fChain->GetBranch("Toftntu");
   b_Ivntu = fChain->GetBranch("Ivntu");
   b_Ipntu = fChain->GetBranch("Ipntu");
   

  return kTRUE;
}

void h10::Show(Int_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}

#ifndef h10_cut
#define h10_cut
Int_t h10::Cut(Int_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
// standard implementation if not defined earlier
   return 1;
}
Int_t h10::Cut(Int_t * results){
  return 1;
}
#endif //h10_cut



#endif // #ifndef h10_cxx_done
#endif // #ifdef h10_cxx



