// Implementation of ana10 classes and subclasses
#include "ana.hpp"

#include <string>
#include <sstream>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>


Bool_t ana10::RunList(int runnumber, Int_t &timerun)
{
  Bool_t bset0,bset1,bset2,bset3,bset4,b5,bstandard,bbb0nu;
  int iival,iival1;
  int dateival;
  timerun=0;
  Int_t r = abs(runnumber);

  Int_t barunlist[]={4657,4658,4659,4660,4661,4663,4664,4665,4667,4668,4669,4670,4671,4672,4673,4674,4675,4676,4677,4678,4679,4680,4681,4682,4683,4684,4685,4686,4687,4688,4689,4690,4691,4692,4693,4694,4697,4698,4699,4700,4701,4702,4703,4704,4706,4707,4708,4709,4710,4711,4712,4713,4714,4715,4716,4717,4718,4719,4720,4721,4722,4723,4744,4745,4746,4747,4749,4751,4752,4753,4754,4755,4756,4757,4758,4759,4760,4761,4762,4763,4764,4765,4766,4767,4768,4769,4770,4772,4773,4774,4775,4776,4777,4778,4779,4780,4781,4782,4783,4784,4786,4787,4788,4789,4790};

  
  Bool_t good_runnum=false;

  if(r<1865) return false;
  if(r>6556) return false;

  //Phase 1 data
  if(P1_0nu && (r<3408 || (r>3553 && r<3566)) ) good_runnum=true;
  //if(P1_0nu && (r<3396 || (r>3553 && r<3566)) ) good_runnum=true;
  if(P1 && (r<3396)) good_runnum=true;
  
  
  //Phase2 data
  if(P2a && (r>3395 && r<4859)) good_runnum=true;
  if(P2b && (r>4858 && r<5469)) good_runnum=true;
  
  if(P2_0nu && (r>3407 && r<5469)) good_runnum=true;
  
  //Phase3 , 2007, data
  if(P3 && (r>5468 && r<6556)) good_runnum=true;
  
  //yuri's BA test period
  if(PBA){
    good_runnum=true;
    if(r<4656) good_runnum=false;
    if(r>4791) good_runnum=false;
  }
  
  if(!good_runnum) return false;
  

  if(PBA){
    iival =0;
    for(Int_t j=0;j<106;j++){
      if(r==barunlist[j]) iival=1;
    }
    if(iival==0) return false;
  }	


  if(rstatus[r] == 0){
    FILE *fp=freopen("db.dump","a",stdout);
 
    //Bool_t rini=nemo3::DbMgr::Db().AnalyzeFile(r,"betabeta");
    float Act;
    int Status;
    /*
    if(rini){
      nemo3::DbMgr::Db().Get("RUNSTORE","Status",&iival);
      nemo3::DbMgr::Db().Get("RUNSTORE","TotalTime",&iival1);
      nemo3::DbMgr::Db().Get("PD","DateStart",&dateival);
      cout<<"DB test, run "<<r<<" date "<<dateival<<endl;;

      nemo3::DbMgr::Db().Get("RN","Status",&Status);
      if (Status==0){
	nemo3::DbMgr::Db().Get("RN","Act",&Act);
      }
      if(Status!=0 || Act < 1){
	cout<<"Warning! Run "<<r<<" has no RN activity stored in DB"<<endl;
	Act = 213;// Use mean acticity durin P2_0nu period as a temporar solution
      }
      
      fclose(fp);
      //restore stdout
      fp=freopen("/dev/tty","w",stdout);
      if(PBA){
	for(Int_t j=0;j<106;j++){
	  if(r==barunlist[j]) iival=1;
	}
      }	
      if(Act>400 && (r>3395)) iival+=200000; // apply for Phase II runs only

      //play with activity cut
      //   if (Act > 450) iival=0;
      }else{
      rstatus[r]=-2;
      timerun=0;
      return false;
    }*/
    rstatus[r]=iival+1;
    rtime[r]=iival1;
    rrn_act[r]=Act;
    run_date[r]=dateival;
  }else{
    iival = rstatus[r]-1;
    iival1=rtime[r];
  }

  timerun=iival1;
  //set0 good runs with status 1
  bset0=set0 & (iival==1);

  //set1 of possibly good runs, exclude statuses 1,2,5,4,1000,6,100000,10000,20000,100000000,1000000000
  bset1=set1& !CheckRunStatus(iival,0);
  bset1=bset1& !CheckRunStatus(iival,1);
  bset1=bset1& !CheckRunStatus(iival,2);
  bset1=bset1& !CheckRunStatus(iival,5);
  bset1=bset1& !CheckRunStatus(iival,6);
  bset1=bset1& !CheckRunStatus(iival,100000);
  bset1=bset1& !CheckRunStatus(iival,4);
  bset1=bset1& !CheckRunStatus(iival,1000);
  bset1=bset1& !CheckRunStatus(iival,10000);
  bset1=bset1& !CheckRunStatus(iival,20000);
  bset1=bset1& !CheckRunStatus(iival,100000000);
  bset1=bset1& !CheckRunStatus(iival,1000000000);
  
  //set2 runs with ventillation off (6 & 100000)
  bset2 = set2 & (CheckRunStatus(iival,6) | CheckRunStatus(iival,100000) | CheckRunStatus(iival,300000));
  
  //set3 runs after absolute calibration (4 & 1000)
  bset3 = set3 & (CheckRunStatus(iival,4) | CheckRunStatus(iival,1000));
  
  //set4 runs with HV problems
  bset4 = set4 & (CheckRunStatus(iival,5) | CheckRunStatus(iival,10000) | CheckRunStatus(iival,20000) | CheckRunStatus(iival,100000000) | CheckRunStatus(iival,1000000000) );
  
  //standard set of runs, proposed by VT for common analysis


  bstandard= standard & !CheckRunStatus(iival,0);
  bstandard= bstandard & !CheckRunStatus(iival,100);
  bstandard= bstandard & !CheckRunStatus(iival,100000000);
  bstandard= bstandard & !CheckRunStatus(iival,1000000000);
  bstandard= bstandard & !CheckRunStatus(iival,3);
  bstandard= bstandard & !CheckRunStatus(iival,7);
  bstandard= bstandard & !CheckRunStatus(iival,8);
  bstandard= bstandard & !CheckRunStatus(iival,100000);
  bstandard= bstandard & !CheckRunStatus(iival,200000);//runs with hi Rn > 400 mBq in PII
  bstandard= bstandard & !CheckRunStatus(iival,300000);//runs with hi Rn > 400 mBq in PII
  bstandard= bstandard & !CheckRunStatus(iival,6);
  bstandard= bstandard & !CheckRunStatus(iival,2000000);
  bstandard= bstandard & !CheckRunStatus(iival,3000000);

  //runs for bb0nu analysis

  bbb0nu = bb0nu & (iival==1 || iival==2 || iival==10 || iival==200002 || iival==200010 || iival==200001);

  if(!extended){
    return bset0 | bset1 | bset2 | bset3 | bset4 | bstandard | bbb0nu;
  }else{
    if(CheckRunStatus(iival,0)) return false;
    if(CheckRunStatus(iival,6)) return false;
    if(CheckRunStatus(iival,7)) return false;
    if(CheckRunStatus(iival,100000)) return false;
    return true;
  }
}

void ana10::SetPeriod(ana::period pr,Bool_t flag){
  
  if(pr==ana::P1)     P1=flag;
  if(pr==ana::P2a)    P2a=flag;
  if(pr==ana::P2b)    P2b=flag;
  if(pr==ana::PBA)    PBA=flag;
  if(pr==ana::P2007)  P3=flag;
  if(pr==ana::P3)     P3=flag;
  if(pr==ana::P1_0nu) P1_0nu=flag;
  if(pr==ana::P2_0nu) P2_0nu=flag;
  

}

void ana10::SetRuns(ana::run rn,Bool_t flag){
  if(rn==ana::GOOD_RUN)      set0 = flag;
  if(rn==ana::OK_RUN)        set1 = flag;
  if(rn==ana::NOAIR_RUN)     set2 = flag;
  if(rn==ana::AFTER_EC_RUN)  set3 = flag;
  if(rn==ana::HV_RUN)        set4 = flag;
  if(rn==ana::ALL_BUT_BAD_RUN) extended = flag;
  if(rn==ana::STANDARD_RUN)  standard = flag;
  if(rn==ana::BB0NU_RUN)     bb0nu = flag;
}

Int_t ana10::GetRunDate(int r){

 Int_t dummy;
 //have we retrived DB info about  this run already?
 if(rstatus[r] == 0){
   //if not do so
   RunList(r,dummy);
 }  
 return run_date[r];
}

Int_t ana10::GetYear(int date){
  return (date/10000);

}
Int_t ana10::GetMonth(int date){
  return (date%10000)/100;
}

Int_t ana10::GetDay(int date){
  return (date%100);
}

//number of days since NEMO era start 01.01.2003
Int_t ana10::Days_Nemo_Era(int date){
  float days;
  float m_days[] = {31,28.25,31,30,31,30,31,31,30,31,30,31};
  days = (GetYear(date) - 2003)*365.25;
  for(int i =0 ;i < GetMonth(date)-1;i++){
    days += m_days[i];
  }
  days+=GetDay(date);
  return int(days);
}
