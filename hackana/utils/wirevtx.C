#include "h10.h"
#include <fstream>
#include <iostream>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

Int_t h10::vtx_volume(Int_t trk1, Int_t trk2, Float_t *vtx, Float_t& dvx, Float_t& dvz, Bool_t finescan=false){
  /** calculates vertex between the two tracks in the volume (not on the foil)

  **/

  Int_t trk[2] = {trk1,trk2};
  

  //Use directional cosines for t0 calculations;

  Float_t t[2];
  for(Int_t i = 0 ;i<2;i++){
    // std::cout<<"X,Y,R = "<<Xc[trk[i]]<<" "<<Yc[trk[i]]<<" "<<Radc[trk[i]]<<std::endl;
   t[i] = acos(-Qc[trk[i]]*Cos_dir[trk[i]][1]/sqrt(pow(Cos_dir[trk[i]][0],2) + pow(Cos_dir[trk[i]][1],2)));
    if(-Qc[trk[i]]*Cos_dir[trk[i]][0]>0){
      t[i] *= -1;
    }
    if(fabs(Zc[trk[i]]+Hc[trk[i]]*t[i]-Z_foil[trk[i]])>0.2){
      Float_t t2=(Z_foil[trk[i]]-Zc[trk[i]])/Hc[trk[i]];
      Int_t  k = (int)( (t[i]-t2)/2/3.14159265);
      t[i] += k*2*3.14159265;
    }
  }

  //find the minimum now

  Float_t vdist[200][4];// radius, t1,t2, & delta V

  Float_t mindist = 100000;
  Int_t minpos = 0;
  Int_t count=0;
  Float_t Rstart, Rend, Rstep;
  if(finescan){
    Rstart = 149;
    Rend=161;
    Rstep=0.25;
  }else{
    Rstart=102;
    Rend=208;
    Rstep=1.;
  }
  for(Float_t r=Rstart;r<Rend;r+=Rstep){
    Float_t t1 = helix_cylinder_intersection(trk[0],t[0],r);
    Float_t t2 = helix_cylinder_intersection(trk[1],t[1],r);
    
    Float_t dist = pow((Xc[trk[0]]+Radc[trk[0]]*cos(t1) - Xc[trk[1]] - Radc[trk[1]]*cos(t2)),2);
    dist += pow((Yc[trk[0]]+Radc[trk[0]]*sin(t1) - Yc[trk[1]] - Radc[trk[1]]*sin(t2)),2);
    dist += pow((Zc[trk[0]]+Hc[trk[0]]*t1 - Zc[trk[1]] - Hc[trk[1]]*t2),2);

    dist = sqrt(dist);
    vdist[count][0] = r;
    vdist[count][1] = t1;
    vdist[count][2] = t2;
    vdist[count][3] = dist;

    if(dist < mindist){
      mindist = dist;
      minpos = count;
    }
    count++;
  }

  //Find the vertex
  vtx[0] =(Xc[trk[0]]+Radc[trk[0]]*cos(vdist[minpos][1]) + Xc[trk[1]] + Radc[trk[1]]*cos(vdist[minpos][2]))/2;
  vtx[1] =(Yc[trk[0]]+Radc[trk[0]]*sin(vdist[minpos][1]) + Yc[trk[1]] + Radc[trk[1]]*sin(vdist[minpos][2]))/2;
  vtx[2] =(Zc[trk[0]]+Hc[trk[0]]*vdist[minpos][1] + Zc[trk[1]] + Hc[trk[1]]*vdist[minpos][2])/2;

  dvx =  pow( (Xc[trk[0]]+Radc[trk[0]]*cos(vdist[minpos][1]) - Xc[trk[1]] - Radc[trk[1]]*cos(vdist[minpos][2])),2);
  dvx +=  pow( (Yc[trk[0]]+Radc[trk[0]]*sin(vdist[minpos][1]) - Yc[trk[1]] - Radc[trk[1]]*sin(vdist[minpos][2])),2);
  dvx = sqrt(dvx);
  dvz =  fabs(Zc[trk[0]]+Hc[trk[0]]*vdist[minpos][1] - Zc[trk[1]] - Hc[trk[1]]*vdist[minpos][2]);

  //  std::cout<<"Volume vertex: minpos = "<<minpos<<" r="<<vdist[minpos][0]<<" t1= "<<vdist[minpos][1]<<" t2="<<vdist[minpos][2]<<std::endl;
  //   std::cout<<"Volume vertex  "<<vtx[0]<<" "<<vtx[1]<<" "<<vtx[2]<<" distnaces X:"<< dvx<<" Z:"<<dvz<<std::endl;
  return 1;

}

Float_t h10::helix_cylinder_intersection(Int_t trk,Float_t t0, Float_t radius){
  /** Calculates intersection between the track trk and the cylinder of
      radius radius. Pics up the solution closest to the t0
  **/

  Float_t R2=pow(Radc[trk],2);
  Float_t R4 = R2*R2;
  Float_t X2=pow(Xc[trk],2);
  Float_t X4 = X2*X2;
  Float_t Y2=pow(Yc[trk],2);
  Float_t Y4 = Y2*Y2;
  Float_t radius2=radius*radius;
  Float_t radius4=radius2*radius2;

  Float_t a5 = 1/(R2*(X2+Y2));
  Float_t a1 = a5*R2*Y2*radius2;
  Float_t a2 = a5*R4*Y2;
  Float_t a3 = a5*R2*Y2*X2;
  Float_t a4 = a5*R2*Y4;
  Float_t a6 = a5*Radc[trk]*Yc[trk];

  Float_t D = -radius4*R2*X2+2*radius2*R4*X2-R2*R4*X2+
    2*radius2*R2*X4+ 2*R4*X4- R2*X2*X4+
    2*radius2*R2*X2*Y2+ 2*R4*X2*Y2- 2*R2*X4*Y2- R2*X2*Y4;

  //Four different solutions (with +- 2Pi accuracy) possible
  Float_t t[4];
  t[0] = -acos((radius2-R2-X2-Y2-a1+a2+a3+a4-a6*sqrt(D))/2/Radc[trk]/Xc[trk]);
  t[1] = +acos((radius2-R2-X2-Y2-a1+a2+a3+a4-a6*sqrt(D))/2/Radc[trk]/Xc[trk]);
  t[2] = -acos((radius2-R2-X2-Y2-a1+a2+a3+a4+a6*sqrt(D))/2/Radc[trk]/Xc[trk]);
  t[3] = +acos((radius2-R2-X2-Y2-a1+a2+a3+a4+a6*sqrt(D))/2/Radc[trk]/Xc[trk]);

  //find the closest to the vertex on the foil, t0
  Float_t mindist = 10.;
  Int_t solution =0;
  for(Int_t i=0;i<4;i++){
    Int_t k=0;
    Float_t dist = t[i]-t0;
    if(fabs(dist)>6.28){
      k = (int)(dist / 2 /3.1415926538);
      dist = dist - k * 2 * 3.141592638;
    }
    if (fabs(dist)< mindist){
      solution = i;
      t[i] = t[i] -  k * 2 * 3.141592638;
      mindist = fabs(dist);
    }
  }
  /**
  std::cout<<"X,Y,R,Q = "<<Xc[trk]<<" "<<Yc[trk]<<" "<<Radc[trk]<<" "<<Qc[trk]<<" event="<<Myievent<<std::endl;
  std::cout<<"a1-6, D, radius "<<a1<<" "<<a2<<" "<<a3<<" "<<a4<<" "<<a5<<" "<<a6<<" "<<D<<" "<<radius<<std::endl;
  std::cout<<"cos1= "<<(pow(radius,2)-pow(Radc[trk],2)-pow(Xc[trk],2)-pow(Yc[trk],2)-a1+a2+a3+a4-a6*sqrt(D))/2/Radc[trk]/Xc[trk]<<std::endl;
  std::cout<<"cos2= "<<(pow(radius,2)-pow(Radc[trk],2)-pow(Xc[trk],2)-pow(Yc[trk],2)-a1+a2+a3+a4-a6*sqrt(D))/2/Radc[trk]/Xc[trk]<<std::endl;
  std::cout<<"T1,2,3,4 "<<t[0]<<" "<<t[1]<<" "<<t[2]<<" "<<t[3]<<std::endl;
  std::cout<<"return the value "<<t[solution]<<" t0="<<t0<<std::endl;
  **/
  return t[solution];

}


