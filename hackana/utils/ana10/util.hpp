#ifndef A10UTIL_H
#define A10UTIL_H

//#include "../pmstatus.h"

#define RUNMIN 1800
#define RUNMAX 7000
#define PMMAX  1940


extern Int_t rtime[RUNMAX];
extern Int_t rstatus[RUNMAX];
extern Int_t run_date[RUNMAX];
extern Float_t rrn_act[RUNMAX];



namespace ana{
  // Possible run periods:
  // Phase 1 period from Feb 2003 to Oct 2004
  // Phase 2 period a) from Oct 2004 to 22 may 2006
  // Phase 2 period b) from 22 may 2006 to 31 Dec 2007
  // PBA run list for blind analysis test from Yuriy 
  enum period{P1,P2a,P2b,PBA,P2007,P3,P1_0nu,P2_0nu};


  // Possible Run sets:
	// GOOD_RUN = Status 1
	// OK_RUN = Possibly good runs, exclude statuses 1,2,5,4,1000,6,100000,10000,20000,100000000,1000000000
	// NOAIR_RUN = Ventillation off (status 6 & 100000)
	// AFTER_EC_RUN = After absolute calibration (4 & 1000)
	// HV_RUN = Runs with HV problems, statuses 5 ,10000 ,20000 ,100000000, 1000000000
	// ALL_BUT_BAD_RUN = all except statuses 0, 6, 7, 100000
        //STANDARD_RUN = runs from starndard runlist, decided at Feb'07 Orsay analysis meeting. See meeting minutes.

  enum run{GOOD_RUN,OK_RUN,NOAIR_RUN,AFTER_EC_RUN,HV_RUN,ALL_BUT_BAD_RUN,STANDARD_RUN,BB0NU_RUN};
}


#endif
