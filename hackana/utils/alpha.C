/////////////////////////////////////////////////////////////////////
//                                                                 //
//  Routines to search for alpha-particle presense in the event    //
//  V.Vasiliev. V1.0 11.2004                                       //
/////////////////////////////////////////////////////////////////////


#include "TMath.h"
#include <stdio.h>
#include <iostream>
#include "h10.h"
#include <fstream>




// Functions for alpha particle analisis and Radon rejection.

Int_t h10::NAfasthits()
{
  //returns 0 if everything is OK
  // returns >=2 if there are too many NA fast hits

  if(fastGG_other_side(15)==1) return 3;
  return NA_fast_hits(15);

}

Int_t h10::NA_fast_hits(Float_t max_dst)
{
  //calculates number of NA GG fast hits close to the vertex]
  // (distance < max_dst) 
  //if there are too much of them returns n -- group of n NA hits close to vertex
  Int_t nagg[3000];

  Int_t Nnagg=0;
  for(Int_t i=0;i<Ngg;i++){
    if(mywstatus(i)!=1) continue;
    Int_t flag=0;
    for(Int_t j=0;j<Nbr_tks;j++){
      for(Int_t g=0;g<Nbr_pts[j];g++){
	if(i==Ind_points[j][g]-1) flag=1;
      }
    }
    if(flag==0){
      nagg[Nnagg]=i;
      Nnagg++;
    }
  }

  //Calculate total number NA fast hits close to vertex
  Int_t counter[2]={0,0}; 
  for(Int_t i=0;i<Nnagg;i++){
      Int_t igg=nagg[i];
      Float_t dst=(Gg[igg][9]-Xvert)*(Gg[igg][9]-Xvert);
      dst += (Gg[igg][10]-Yvert)*(Gg[igg][10]-Yvert);
      dst=sqrt(dst);
      if(dst<max_dst){
        counter[int(Gg[igg][2]+0.5)]++;
      }
  }
  Int_t cmax=counter[0];
  if(cmax<counter[1]) cmax=counter[1];
  if(cmax>2) return cmax;
  return 0;
}
Int_t h10::fastGG_other_side(Float_t max_dst)
{
  // if two tracks are in the same half of the detector
  // checks if there are any fast GG hits on the other side of the foil
  // close to the vertex (distance < max_dst)

  Int_t nagg[3000];

  Int_t Nnagg=0;
  for(Int_t i=0;i<Ngg;i++){
    if(mywstatus(i)!=1) continue;
    Int_t flag=0;
    for(Int_t j=0;j<Nbr_tks;j++){
      for(Int_t g=0;g<Nbr_pts[j];g++){
	if(i==Ind_points[j][g]-1) flag=1;
      }
    }
    if(flag==0){
      nagg[Nnagg]=i;
      Nnagg++;
    }
  }
  //Check if tracks are on the same side
  Int_t ioevent=int(Gg[Ind_points[0][0]-1][2]+0.5);
  for(int k=0;k<Nbr_tks;k++){
    if(ioevent!=Gg[Ind_points[k][0]-1][2]){
      ioevent=-1;
    }
  }
  if(ioevent>=0){
    for(Int_t i=0;i<Nnagg;i++){
      if(Gg[nagg[i]][2]!=ioevent){
      //Calculate VTX distance
      Int_t igg=nagg[i];
      Float_t dst=(Gg[igg][9]-Xvert)*(Gg[igg][9]-Xvert);
      dst += (Gg[igg][10]-Yvert)*(Gg[igg][10]-Yvert);
      dst=sqrt(dst);
      if(dst<max_dst){
	        return 1;
      }
      }
    }
  }
  return 0;
}

Int_t h10::mywstatus(Int_t igg)
{
  //This function calculates gg hit status:
  // 1 -- Good fast hit
  // 10 -- Good delayed hit
  // 0 -- failed to calculate status
  Int_t ifast,islow;
  Int_t mystatus=0;
  islow=int(Gg[igg][6]+0.5);
  ifast=int(Gg[igg][5]+0.5);
  //  std::cout<<" mywstatus "<<islow<<" "<<ifast<<std::endl;
  if(run>1865){ // REAL data
    if(islow==0 && ifast > 220 && ifast<307){
      //   hit in time
      return 1;
    }
    if(islow-4096*(islow/4096)!= ifast &&  ifast>220 && ifast<307){
      //   hit in time, electronic problems 
      return 1;
    }
    if(islow-4096*(islow/4096)== ifast){
      //   delayed hit 
      return 10;
    }
  }else{ //MTCA data
    if(islow==0){
      return 1;
    }else{
      return 10;
    }
  }
    // else it is noise
  return 0;
}

Int_t h10::checkdelayed()
{
  // Search for delayed hits in the event
  // Produces list dgg with  Ndgg records in it
  Ndgg=0;
  for(Int_t i=0;i<Ngg;i++){
    if(mywstatus(i)==10){
      dgg[Ndgg]=i;
      Ndgg++;
    }
  }
  return 1;
}

Int_t h10::checksingledhit(Float_t &atime){
  Float_t dist;
  return checksingledhit(atime,dist);
}
 Int_t h10::checksingledhit(Float_t &atime, Float_t &dist)
   {
     //Check if there is a single delayed hit close to the vertex of the event?
     // returns number of such a hits and there mean  time
     Float_t dxymax=25.;
     Float_t dzmax=30.;
     Float_t timemin=100000.; // safe time, beyond wich there is no refiring

     Float_t dtime;
     Int_t counter=0;
     Float_t dstmin=100;
     atime=0.;
     for(Int_t i=0;i<Ndgg;i++){
       Float_t dstxy=(Gg[dgg[i]][9]-Xvert)*(Gg[dgg[i]][9]-Xvert);
       dstxy+=(Gg[dgg[i]][10]-Yvert)*(Gg[dgg[i]][10]-Yvert);
       dstxy=sqrt(dstxy);
       Float_t dstz=fabs(Gg[dgg[i]][11]-Zvert);

       dtime=Gg[dgg[i]][12];
       if(dstxy<dxymax && dstz<dzmax && dtime>timemin){
	 counter++;
	 atime+=dtime;
	 if(dstxy<dstmin){
	   dstmin=dstxy;
	 }
       }   
     }
     if(counter !=0) atime=atime/counter;
     dist=dstmin;
     if(counter==1) return counter;
     return 0;
   }
 Int_t h10::checksingledhit2(Float_t &atime, Float_t &dist)
   {
     //Check if there is a single delayed hit close to the vertex of the event?
     // returns number of such a hits and there mean  time

     //dist -- distance between closest delayed hit and electron track start

     Float_t dxymax=25.;
     Float_t dzmax=30.;
     Float_t timemin=100000.; // safe time, beyond wich there is no refiring

     Float_t dtime;
     Int_t counter=0;
     Float_t dstmin=100;
     atime=0.;
     //find electron track start
     Float_t bx,by,bz,blr;
     blr=8;
     for(Int_t i=0;i<Nbr_tks;i++)
       for(Int_t j=0;j<Nbr_pts[i];j++){
	 Int_t igg=Ind_points[i][j]-1;
	 if(Gg[igg][3]<blr){
	   blr = Gg[igg][3];
	   bx=Gg[igg][9];
	   by=Gg[igg][10];
	   bz=Gg[igg][11];
	 }
     }

     for(Int_t i=0;i<Ndgg;i++){
       Float_t dstxy=(Gg[dgg[i]][9]-Xvert)*(Gg[dgg[i]][9]-Xvert);
       dstxy+=(Gg[dgg[i]][10]-Yvert)*(Gg[dgg[i]][10]-Yvert);
       dstxy=sqrt(dstxy);
       Float_t dstz=fabs(Gg[dgg[i]][11]-Zvert);

       dtime=Gg[dgg[i]][12];
       if(dstxy<dxymax && dstz<dzmax && dtime>timemin){
	 counter++;
	 atime+=dtime;
	 Float_t bdstxy=(Gg[dgg[i]][9]-bx)*(Gg[dgg[i]][9]-bx);
	 bdstxy+=(Gg[dgg[i]][10]-by)*(Gg[dgg[i]][10]-by);
	 bdstxy=sqrt(bdstxy);
	 if(bdstxy<dstmin){
	   dstmin=bdstxy;
	 }
       }   
     }
     if(counter !=0) atime=atime/counter;
     dist=dstmin;
     if(counter==1) return counter;
     return 0;
   }

 Int_t h10::checkggclose(Int_t igg1, Int_t igg2,Float_t dxymax,Float_t dzmax){
   Float_t distxy=(Gg[igg1][9]-Gg[igg2][9])*(Gg[igg1][9]-Gg[igg2][9]);
   distxy+=(Gg[igg1][10]-Gg[igg2][10])*(Gg[igg1][10]-Gg[igg2][10]);
   distxy=sqrt(distxy);

   Float_t distz=fabs(Gg[igg1][11]-Gg[igg2][11]);
   Int_t ret=0;
   if(distxy < dxymax && distz < dzmax && Gg[igg1][2]==Gg[igg2][2]) ret=1;
   return ret;
 }
 Int_t h10::checkgroupdhit(Float_t &atime){
   Float_t dist;
   Float_t dtmax;
   Int_t result= checkgroupdhit(atime, dtmax,dist);
   if(dist>36) result=0; // alpha track cannot be longer then 35 cm
   return result;
 }

 Int_t h10::checkgroupdhit(Float_t &atime, Float_t &adtime, Float_t &alength){
   //Check if there is a group of aligned int time delayed hits, close to any track
   Float_t dtmax=4000.;// Maximal interval between aligned hits 2000 ns.
   Float_t tdhitmin=20000.;//Minimal time for delayed hit, to reject refiring
   Int_t dhitsmin=2; // Minimal Number of minimal dhits in the group demanded
   Float_t gdxymax=25.;
   Float_t gdzmax=30.;


   Int_t ndgroup=0;
   Int_t ndg[100];
   Int_t dgroup[100][100];
   

   //Search for groups of aligned hits
   for(Int_t i=0;i<Ndgg;i++){
     Int_t igg=dgg[i];
     Float_t t1=Gg[igg][12];
     if(t1>tdhitmin){
       Int_t flag=0;
       for(Int_t j=0;j<ndgroup;j++){
	 for(Int_t k=0;k<ndg[j];k++){
	   Float_t dtt=fabs(t1-Gg[dgroup[j][k]][12]);
	   if(dtt<dtmax && flag==0 && Gg[dgroup[j][k]][12] > tdhitmin ){
	     dgroup[j][ndg[j]]=igg;
	     ndg[j]++;
	     flag=1;
	   }
	 }
       }
       if(flag==0){
	 ndg[ndgroup]=1;
	 dgroup[ndgroup][0]=igg;
	 ndgroup++;
       }
     }
   }
   //Now look of any group is close to any track
     Int_t flag=0;
     atime=0.;
     Float_t gtime;
     Int_t idgg,igg;
     Int_t counter=0;
     Int_t gsize;
     for(Int_t i=0;i<ndgroup;i++){
       flag=0;
       gtime=0;
       Float_t dtmin=1000000;
       Float_t dtmax=0;
       if(ndg[i]>=dhitsmin){

	 // Claculate 2 lengths: using precise hits and hits with big Z error separetly
	 Float_t length=0;
	 Int_t rhit=-1;
	 Int_t lhit=-1;

	 Float_t u_length=0;
	 Float_t u_length_err = 0;
	 Int_t u_rhit=dgroup[i][0];
	 Int_t u_lhit=u_rhit;

	 for(Int_t gg=0;gg<ndg[i];gg++){
	   idgg=dgroup[i][gg];
	   gtime+=Gg[idgg][12];
	   
 

	   //calculating delayed hit length as distance between edge hits
	   //for good Z resolution first
	   if (Gg[idgg][13]<2.){
	     if(rhit==-1){
	       rhit=idgg;
	       lhit=idgg;
	     }
	     Float_t dist1=sqrt(pow(Gg[idgg][9]-Gg[rhit][9],2)+pow(Gg[idgg][10]-Gg[rhit][10],2)+pow(Gg[idgg][11]-Gg[rhit][11],2));
	     Float_t dist2=sqrt(pow(Gg[idgg][9]-Gg[lhit][9],2)+pow(Gg[idgg][10]-Gg[lhit][10],2)+pow(Gg[idgg][11]-Gg[lhit][11],2));

	     if(Gg[idgg][12]>dtmax) dtmax = Gg[idgg][12];
	     if(Gg[idgg][12]<dtmin) dtmin = Gg[idgg][12];
	     if(dist1>length){
	       lhit = idgg;
	       length=dist1;
	     }
	     if(dist2>length){
	       rhit=idgg;
	       length=dist2;
	     }
	   }
	   //now using all GG hits in the group
	   Float_t dist1=sqrt(pow(Gg[idgg][9]-Gg[u_rhit][9],2)+pow(Gg[idgg][10]-Gg[u_rhit][10],2)+pow(Gg[idgg][11]-Gg[u_rhit][11],2));
	   Float_t dist2=sqrt(pow(Gg[idgg][9]-Gg[u_lhit][9],2)+pow(Gg[idgg][10]-Gg[u_lhit][10],2)+pow(Gg[idgg][11]-Gg[u_lhit][11],2));
	   
	   if(Gg[idgg][12]>dtmax) dtmax = Gg[idgg][12];
	   if(Gg[idgg][12]<dtmin) dtmin = Gg[idgg][12];
	   if(dist1>u_length){
	     u_lhit = idgg;
	     u_length=dist1;
	   }
	   if(dist2>u_length){
	     u_rhit=idgg;
	     u_length=dist2;
	   }
	   

	   //       //Check if this group is close to any track
	   //       for(Int_t j=0;j<Nbr_tks;j++){
	   // for(Int_t k=0;k<Nbr_pts[j];k++){
	   //   igg=Ind_points[j][k]-1;
	   //     if(checkggclose(igg,idgg,gdxymax,gdzmax)==1) flag=1;
	   // }

	   //Check if this group is close to vertex
	   Float_t dstxy=pow(Gg[idgg][9]-Xvert,2);
	   dstxy+=pow(Gg[idgg][10]-Yvert,2);
	   dstxy=sqrt(dstxy);
	   Float_t dstz=fabs(Gg[idgg][11]-Zvert);
	   if(dstxy<gdxymax && dstz<gdzmax) flag=1;
	   
	 }


	 //Compare lengths at this stage
	 u_length_err = fabs(Gg[u_rhit][11]-Gg[u_lhit][11])/u_length*sqrt(pow(Gg[u_rhit][13],2)+pow(Gg[u_lhit][13],2));
	 if(length < (u_length - u_length_err) ) length = u_length - u_length_err;


	 gtime=gtime/ndg[i];
	 if(flag!=0 && gtime>tdhitmin) {
	   atime=gtime;
	   counter++;
	   adtime = (dtmax - dtmin)/1000.;
	   alength=length;
	   gsize=ndg[i];
	 }
       }
     }
     if(atime!=0 and counter==1) {
       return gsize;
     }
     return 0;
 }
 Int_t h10::checkgroupdhit_vera(Float_t &atime, Float_t &adtime, Float_t &alength){
   //Check if there is a group of aligned int time delayed hits, close to any track
   Float_t dtmax=1600.;// Maximal interval between aligned hits 2000 ns.
   Float_t tdhitmin=20000.;//Minimal time for delayed hit, to reject refiring
   Int_t dhitsmin=2; // Minimal Number of dhits in the group demanded
   Float_t gdxymax=10.;
   Float_t gdzmax=15.;


   Int_t ndgroup=0;
   Int_t ndg[100];
   Int_t dgroup[100][100];
   

   //Search for groups of aligned hits
   for(Int_t i=0;i<Ndgg;i++){
     Int_t igg=dgg[i];
     Float_t t1=Gg[igg][12];
     if(t1>tdhitmin){
       Bool_t assigned = false;
       for(Int_t j=0;j<ndgroup;j++){
	 Bool_t flag=true;
	 for(Int_t k=0;k<ndg[j];k++){
	   Float_t dtt=fabs(t1-Gg[dgroup[j][k]][12]);
	   flag = flag & (dtt<dtmax); 
	 }
	 if(flag){
	   dgroup[j][ndg[j]]=igg;
	   ndg[j]++;
	   assigned = true;
	 }
       }
       if(!assigned){
	 ndg[ndgroup]=1;
	 dgroup[ndgroup][0]=igg;
	 ndgroup++;
       }
     }
   }
   //Now look if any group is close to any track
     Int_t flag=0;
     atime=0.;
     Float_t gtime;

     Int_t idgg,igg;
     Int_t gsize=0;
     Int_t counter=0;
     for(Int_t i=0;i<ndgroup;i++){
       flag=0;
       gtime=0;
       Float_t dtmin=1000000;
       Float_t dtmax=0;
       if(ndg[i]>=dhitsmin){
	 Float_t dst_wire=50;
	 Float_t dst_foil=50;
	 Float_t length=0;
	 Int_t rhit=dgroup[i][0];
	 Int_t lhit=rhit;
	 for(Int_t gg=0;gg<ndg[i];gg++){
	   idgg=dgroup[i][gg];
	   gtime+=Gg[idgg][12];

	   //calculating delayed hit length as distance between edge hits
	   Float_t dist1=sqrt(pow(Gg[idgg][9]-Gg[rhit][9],2)+pow(Gg[idgg][10]-Gg[rhit][10],2)+pow(Gg[idgg][11]-Gg[rhit][11],2));
	   Float_t dist2=sqrt(pow(Gg[idgg][9]-Gg[lhit][9],2)+pow(Gg[idgg][10]-Gg[lhit][10],2)+pow(Gg[idgg][11]-Gg[lhit][11],2));
	   if(Gg[idgg][12]>dtmax) dtmax = Gg[idgg][12];
	   if(Gg[idgg][12]<dtmin) dtmin = Gg[idgg][12];
	   Int_t tmpgg=idgg;
	   if(dist1>length){
	     Int_t tmp=lhit;
	     lhit = tmpgg;
	     tmpgg=lhit;
	     length=dist1;
	   }
	   if(dist2>length){
	     Int_t tmp=rhit;
	     rhit=tmpgg;
	     tmpgg=tmp;
	     length=dist2;
	   }
	 }

       //Check if this group is close to vertex
       idgg=rhit;
       Float_t dstxy=pow(Gg[idgg][9]-Xvert,2);
       dstxy+=pow(Gg[idgg][10]-Yvert,2);
       dstxy=sqrt(dstxy);
       Float_t dstz=fabs(Gg[idgg][11]-Zvert);
       if(dstxy<gdxymax && dstz<gdzmax){
	 flag=1;
	 if(dst_foil>sqrt(dstxy*dstxy+dstz*dstz)) dst_foil=sqrt(dstxy*dstxy+dstz*dstz) ;
       }
       //Check if this group is close to the fast track
       for(Int_t it=0;it<Nbr_tks;it++){
	 Int_t fgg=Ind_points[it][0]-1;
	 Float_t fx = Gg[fgg][9];
	 Float_t fy = Gg[fgg][10];
	 Float_t fz = Gg[fgg][11];
	 dstxy=pow(Gg[idgg][9]-fx,2);
	 dstxy+=pow(Gg[idgg][10]-fy,2);
	 dstxy=sqrt(dstxy);
	 dstz=fabs(Gg[idgg][11]-fz);
	 if(dstxy<gdxymax && dstz<gdzmax){
	   flag=1;
	   if(dst_wire>sqrt(dstxy*dstxy+dstz*dstz)) dst_wire=sqrt(dstxy*dstxy+dstz*dstz) ;
	   
	 }
	 if(Gg[idgg][2]!=Gg[fgg][2]) flag = 0;

       idgg=lhit;
       dstxy=pow(Gg[idgg][9]-Xvert,2);
       dstxy+=pow(Gg[idgg][10]-Yvert,2);
       dstxy=sqrt(dstxy);
       dstz=fabs(Gg[idgg][11]-Zvert);
       if(dstxy<gdxymax && dstz<gdzmax){
	 flag=1;
	 if(dst_foil>sqrt(dstxy*dstxy+dstz*dstz)) dst_foil=sqrt(dstxy*dstxy+dstz*dstz) ;
       }
       //Check if this group is close to the fast track
       for(Int_t it=0;it<Nbr_tks;it++){
	 Int_t fgg=Ind_points[it][0]-1;
	 Float_t fx = Gg[fgg][9];
	 Float_t fy = Gg[fgg][10];
	 Float_t fz = Gg[fgg][11];
	 dstxy=pow(Gg[idgg][9]-fx,2);
	 dstxy+=pow(Gg[idgg][10]-fy,2);
	 dstxy=sqrt(dstxy);
	 dstz=fabs(Gg[idgg][11]-fz);
	 if(dstxy<gdxymax && dstz<gdzmax){
	   flag=1;
	   if(dst_wire>sqrt(dstxy*dstxy+dstz*dstz)) dst_wire=sqrt(dstxy*dstxy+dstz*dstz) ;
	 }	   
       
	 if(Gg[idgg][2]!=Gg[fgg][2]) flag = 0;
       }
       }
	 gtime=gtime/ndg[i];
	 if(flag!=0 ) {
	   atime=gtime;
	   counter++;
	   adtime = (dtmax - dtmin)/1000.;
	   alength=length;
	   gsize=ndg[i] ;
	   return gsize; //return as soon as first suitable group is found.
	 }
       }
     }
     
     if(counter == 1) return gsize;

     return 0;
 }
