#include "ana.hpp"
#include "../../utils/h10.h"

#include "TH1.h"
#include "TH2.h"
#include "TLegend.h"
#include "TFractionFitter.h"
#include "TMath.h"
#include "TLimit.h"
#include "TLimitDataSource.h"
#include "TConfidenceLevel.h"
#include "TPaveText.h"
#include "TLatex.h"

#include "TString.h"
#include "TRegexp.h"
#include "TPRegexp.h"
#include "THStack.h"

#include <string>
#include <sstream>
#include <iomanip>
#include <vector>

/***********************************************************
 These functins create the final results plots that have 
 been blessed by the NEMO collaboration.
***********************************************************/

// VERSION: 09.09.16

// GLOBAL VARIABLES TO USE
////////////////////////////////////////////////////////////
string header  = "NEMO-3";
Int_t TITLES   =             0;  // SET HIST TITLE OR NOT
Int_t LEG_HEAD =             1;  // SET LEGEND TITLE OR NOT

Double_t XLS =           0.050;  // X-AXIS LABEL SIZE
Double_t XTS =           0.050;  // X-AXIS TITLE SIZE
Double_t XTO =             1.1;  // X-AXIS TITLE OFFSET

Double_t YLS =           0.048;  // Y-AXIS LABEL SIZE
Double_t YTS =           0.048;  // Y-AXIS TITLE SIZE
Double_t YTO =             1.6;  // Y-AXIS TITLE OFFSET

Int_t DPS   =               20;  // DATA POINT STYLE
Size_t DMS  =      Size_t(1.0);  // DATA MARKER SIZE
Width_t DLW =     Width_t(2.0);  // DATA LINE WIDTH
Width_t LW  =     Width_t(2.4);  // HISTO LINE WIDTH

const Double_t LX2 =      0.98;  // LEGEND X2 COORD.
const Double_t LY2 =      0.98;  // LEGEND Y2 COORD.
const Double_t LX1 =  LX2-0.28;  // LEGEND X1 COORD.
const Double_t LY1 =       LY2;  // LEGEND Y1 COORD.
const Double_t LCS =     0.020;  // LEGEND CHAR SPACING
const Double_t LLS =     0.045;  // LEGEND LINE SPACING
const Double_t LTS =    0.0350;  // LEGEND TEXT SIZE
////////////////////////////////////////////////////////////

string printHalfStatSyst(Float_t halflife, Float_t staterr, 
			 Float_t systerr, Int_t precision){
  ostringstream oss;
  Int_t base = (Int_t)TMath::Log10(halflife);
  Int_t statp=precision-(base-(Int_t)TMath::Log10(staterr));
  Int_t systp=precision-(base-(Int_t)TMath::Log10(systerr));
  
  oss<<"#Tau^{2#beta2#nu}_{1/2}  =  "<<setprecision(precision)
     <<halflife/pow(10.,base)
     <<"  #pm  "<<setprecision(statp)<<staterr/pow(10.,base)
     <<"(stat)  #pm  "<<setprecision(systp)<<systerr/pow(10.,base)
     <<"(syst)  #times10^{"<<base<<"}";
  
  return oss.str();
}

Double_t chi2Calc(TH1* mydat,TH1* total){
  //return mydat->Chi2Test(total,"UUCHI2");  // sometimes this
  return mydat->Chi2Test(total,"UWCHI2");  // default
  //return mydat->Chi2Test(total,"WWCHI2");
  //return mydat->Chi2Test(total,"NORMCHI2");
}

Double_t chi2Man(TH1* mydat,TH1* total){
  
  Double_t mychi2=0;
  Float_t dat=0,tot=0;
  Float_t daterr=0,toterr=0;
  Int_t bins=mydat->GetNbinsX();
  for(Int_t i=1;i<=bins;i++){
    dat=tot=daterr=toterr=0;
    dat=mydat->GetBinContent(i);
    tot=total->GetBinContent(i);
    daterr=mydat->GetBinError(i);
    toterr=total->GetBinError(i);
    if(mydat->GetBinContent(i)>0){
    //if(mydat->GetBinContent(i)!=0){
      mychi2 += (pow(dat-tot,2))/(pow(daterr,2)+pow(toterr,2));
    }
  }
  
  cout<<"\nManual = "<<setprecision(1)<<mychi2
      <<"\n    UU = "<<mydat->Chi2Test(total,"UUCHI2")
      <<"\n    UW = "<<mydat->Chi2Test(total,"UWCHI2")
      <<"\n    WW = "<<mydat->Chi2Test(total,"WWCHI2")
      <<"\n  NORM = "<<mydat->Chi2Test(total,"NORMCHI2")
      <<endl;
  
  return mychi2;
}

/////////////////////////////////////////////////////////////////////
Double_t chi2Range(TH1* mydat,TH1* total,Double_t xmin,Double_t xmax){
  Double_t wd=mydat->GetBinWidth(1);
  Int_t bmin=Int_t(floor(xmin/wd));
  Int_t bmax=Int_t(ceil(xmax/wd));
  
  Double_t mychi2=0;
  Float_t dat=0,tot=0;
  Float_t daterr=0,toterr=0;
  Int_t bins=mydat->GetNbinsX();
  for(Int_t i=bmin+1;i<=bmax;i++){
    dat=tot=daterr=toterr=0;
    dat=mydat->GetBinContent(i);
    tot=total->GetBinContent(i);
    daterr=mydat->GetBinError(i);
    toterr=total->GetBinError(i);
    if(mydat->GetBinContent(i)>0){
    //if(mydat->GetBinContent(i)!=0){
      mychi2 += (pow(dat-tot,2))/(pow(daterr,2)+pow(toterr,2));
    }
  }
  return mychi2;
}
Int_t ndfRange(TH1* mydat,TH1* total,Double_t xmin,Double_t xmax){
  Double_t wd=mydat->GetBinWidth(1);
  Int_t bmin=Int_t(floor(xmin/wd));
  Int_t bmax=Int_t(ceil(xmax/wd));
  Int_t ndf=0;
  for(Int_t i=bmin+1;i<=bmax;i++){
    if(mydat->GetBinContent(i)>0) ndf++;
    //if(mydat->GetBinContent(i)!=0) ndf++;
  }
  
  return ndf;
}
//////////////////////////////////////////////////////////////////////

Int_t ndfCalc(TH1* mydat){
  Int_t ndf=0;
  Int_t bins=mydat->GetNbinsX();
  for(Int_t i=1;i<=bins;i++){
    if(mydat->GetBinContent(i)>0) ndf++;
    //if(mydat->GetBinContent(i)!=0) ndf++;
  }
  return ndf;  
}

Int_t round(Float_t value){
  Int_t intvalue=TMath::FloorNint(value);
  Float_t remainder=value-intvalue;
  if(remainder<0.50) return intvalue;
  if(remainder>=0.50) return intvalue+1;
}

Color_t ncolor(Int_t n=0){
  switch(n){
  case 0:return  kWhite;    case 1:return  kBlack;  case 2:return  kGray+2;
  case 3:return  kBlue+1;   case 4:return  kRed+1;  case 5:return  kGreen+1;
  case 6:return  kMagenta+2;case 7:return  kCyan+2; case 8:return  kOrange+2;
  case 9:return  kPink+2;   case 10:return kAzure+2;case 11:return kSpring+2;
  case 12:return kRed-4;    case 13:return kBlue-4; case 14:return kGreen-4;
  case 15:return kMagenta-4;case 16:return kCyan-4; case 17:return kOrange-4;
  case 18:return kPink-4;   case 19:return kAzure-4;case 20:return kSpring-4;
  default:return 40;
  }
}

Bool_t sfind(string s1, string s2){
  return (s1.find(s2) != string::npos);
}

Bool_t prfind(string s11, string s22){
  TString s1=TString(s11);
  TPRegexp s2=TPRegexp(s22);
  return s1.Contains(s2);
}

Int_t sveto(string s1, vector<string> v1){
  for(Int_t i=0; i<v1.size(); i++){
    if(prfind(s1,v1[i])) return 1;
  }
  return 0;
}

string match_bkg(string one){
  if(prfind(one,"([Aa][Cc]|(228)+){2}"))  return "Ac-228";
  if(prfind(one,"([Bb][Ii]|(212)+){2}"))  return "Bi-212";
  if(prfind(one,"([Tt][Ll]|(208)+){2}"))  return "Tl-208";
  if(prfind(one,"([Bb][Ii]|(214)+){2}"))  return "Bi-214";
  if(prfind(one,"([Pp][Bb]|(214)+){2}"))  return "Pb-214";
  if(prfind(one,"([Pp][Aa]|(234m)+){2}")) return "Pa-234m";
  if(prfind(one,"([Tt][Ll]|(207)+){2}"))  return "Tl-207";
  if(prfind(one,"([Pp][Bb]|(211)+){2}"))  return "Pb-211";
  if(prfind(one,"([Ss][Rr]|(90)+){2}"))   return "Sr-90";
  if(prfind(one,"([Bb][Ii]|(207)+){2}"))  return "Bi-207";
  if(prfind(one,"([Ee][Uu]|(152)+){2}"))  return "Eu-152";
  if(prfind(one,"([Ee][Uu]|(154)+){2}"))  return "Eu-154";
  if(prfind(one,"([Bb][Ii]|(210)+){2}"))  return "Bi-210";
  if(prfind(one,"([Cc][Oo]|(60)+){2}"))   return "Co-60";
  if(prfind(one,"([Yy]|(90)+){2}"))       return "Y-90";
  if(prfind(one,"([Kk]|(40)+){2}"))       return "K-40";
  return "";
}

string match_bb(string one){
  if(prfind(one,"([Zz][Rr]|(96)+){2}"))   return "Zr-96";
  if(prfind(one,"([Cc][Aa]|(48)+){2}"))   return "Ca-48";
  if(prfind(one,"([Nn][Dd]|(150)+){2}"))  return "Nd-150";
  if(prfind(one,"([Cc][Dd]|(116)+){2}"))  return "Cd-116";
  if(prfind(one,"([Tt][Ee]|(130)+){2}"))  return "Te-130";
  if(prfind(one,"([Ss][Ee]|(82)+){2}"))   return "Se-82";
  if(prfind(one,"([Mm][Oo]|(100)+){2}"))  return "Mo-100";
  return "";
}

void setStyle(TH1 *hist){
  hist->GetXaxis()->SetLabelSize(XLS);
  hist->GetXaxis()->SetTitleSize(XTS);
  hist->GetXaxis()->SetTitleOffset(XTO);
  hist->GetYaxis()->SetLabelSize(YLS);
  hist->GetYaxis()->SetTitleSize(YTS);
  hist->GetYaxis()->SetTitleOffset(YTO);
}

Int_t axisType(string s_type){
  if(s_type=="EE"  ||  s_type=="ee" ) return 1;
  if(s_type=="E"   ||  s_type=="e" )  return 2;
  if(s_type=="COS" ||  s_type=="Cos" || s_type=="cos" ) return 3;
  if(s_type=="EG"  ||  s_type=="eg" ) return 4;
  if(s_type=="G"   ||  s_type=="g" )  return 5;
  if(s_type=="CM"  ||  s_type=="cm")  return 6;
  if(s_type=="") return 0;
  return -1;
}

void setData(TH1 *mydat,string s_type){
  if(axisType(s_type)==-1) mydat->SetXTitle(s_type.c_str());
  if(axisType(s_type)==0)  mydat->SetXTitle("");
  if(axisType(s_type)==1)  mydat->SetXTitle("E_{1} + E_{2} (MeV)");
  if(axisType(s_type)==4)  mydat->SetXTitle("E_{e} + E_{#gamma} (MeV)");
  if(axisType(s_type)==2)  mydat->SetXTitle("E_{e} (MeV)");
  if(axisType(s_type)==5)  mydat->SetXTitle("E_{#gamma} (MeV)");
  if(axisType(s_type)==3)  mydat->SetXTitle("Cos(#theta)");
  if(axisType(s_type)==6)  mydat->SetXTitle("cm");
  
  ostringstream oss;
  if(axisType(s_type)==-1) oss<<"counts";
  if(axisType(s_type)==0)  oss<<"counts";
  if(axisType(s_type)==1)  oss<<"counts / "<<mydat->GetBinWidth(0)<<" MeV ";
  if(axisType(s_type)==2)  oss<<"counts / "<<mydat->GetBinWidth(0)<<" MeV ";
  if(axisType(s_type)==3)  oss<<"counts / "<<mydat->GetBinWidth(0)<<"";
  if(axisType(s_type)==4)  oss<<"counts / "<<mydat->GetBinWidth(0)<<" MeV ";
  if(axisType(s_type)==5)  oss<<"counts / "<<mydat->GetBinWidth(0)<<" MeV ";
  if(axisType(s_type)==6)  oss<<"counts / "<<mydat->GetBinWidth(0)<<" cm ";
  mydat->SetYTitle(oss.str().c_str());
  
  mydat->SetStats(0);
  mydat->SetMarkerStyle(DPS);
  mydat->SetMarkerSize(DMS);
  mydat->SetLineColor(kBlack);
  mydat->SetLineWidth(DLW);
}

void setStyle(TLegend *leg){
  leg->SetFillColor(0);
  leg->SetTextSize(LTS);
}

void DrawBreakdown(const char *hname,ana10 *data,bana10 **bgr,Int_t nbgr,
		   bana10 **sig,Int_t nsig,Int_t *hdraw,vector<string> foilbkgs,
		   vector<string> vetobkgs,Double_t xmin,Double_t xmax,
		   Int_t rebin,string s_type,Double_t yscale,string isotope){
  
  //**************************************************************************
  //    THIS IS FOR DRAWING INDIVIDUAL BACKGROUND CONTRIBUTIONS
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // 
  // This can draw just your internal backgrounds, of the sum of externals
  // and sum of neighboor foil contributions, etc.
  // 
  // (options)
  // hdraw     : let's you choose which components you want to draw
  //             hdraw[0] = draw data
  //             hdraw[1] = draw signals
  //             hdraw[2] = draw my foil internals
  //             hdraw[3] = draw other foil bkgs
  //             hdraw[4] = draw external bkgs
  //             hdraw[5] = draw total MC
  // foilbkgs  : the neighboor sources you want to draw
  // vetobkgs  : any internal backgrounds you don't want to be drawn
  // xmin      : X-axis min value
  // xmax      : X-axis max value
  // rebin     : Hist rebinning
  // s_type    : Distribution type, "ee", "e", "cos", "eg", "g"
  // yscale    : Multiply the ymax by this value to adjust the scale
  // isotope   : To name your isotope in Latex for the Legend, ie "^{96}Zr",
  //             "^{48}Ca", "^{150}Nd", "^{116}Cd" ect ...
  // 
  //**************************************************************************
  
  string my_int="zr96";
  const Int_t BI=foilbkgs.size();
    
  Int_t tex=0;
  if(isotope!="") tex=1;
  
  Int_t type=axisType(s_type);
    
  TH1D *mydat=(TH1D*)data->Geth1(hname)->Clone("datas"); 
  TH1D *mybak_ext=(TH1D*)bgr[0]->Geth1(hname)->Clone("ext"); 
  TH1D *mysig;
  if(nsig) mysig=(TH1D*)sig[0]->Geth1(hname)->Clone("sig"); 
  TH1D *total=(TH1D*)bgr[0]->Geth1(hname)->Clone("totals"); 
  
  mybak_ext->Reset();
  if(nsig) mysig->Reset();
  total->Reset();
  
  mydat->Sumw2();
  mybak_ext->Sumw2();
  if(nsig) mysig->Sumw2();
  total->Sumw2();
  
  Bool_t attached=false;
  Float_t events[100];
  Int_t index[100];
  TH1D *myint[100];
  string myname[100];
  Int_t M=0;
  TH1D *bkint[100];
  for(Int_t b=0;b<BI;b++){  
    bkint[b]=(TH1D*)bgr[0]->Geth1(hname)->Clone(foilbkgs[b].c_str());
    bkint[b]->Reset();
    bkint[b]->Sumw2();
  }
  for(Int_t i=0;i<nbgr;i++){
    string bname = bgr[i]->Getname();
    attached=false;
    if(!attached && prfind(bname,my_int) && ! sveto(bname,vetobkgs)){
      if(match_bkg(bname) != ""){
	myname[M]=match_bkg(bname);
      }else myname[M]=bname;
      myint[M]=(TH1D*) bgr[i]->Geth1(hname)->Clone(bname.c_str());
      myint[M]->Sumw2();
      myint[M]->Add(myint[M],bgr[i]->Geth1(hname),0.,bgr[i]->norm);
      events[M]=myint[M]->Integral();
      total->Add(myint[M]);
      M++;
      attached=true;
    }
    for(Int_t b=0;b<BI;b++){
      if(!attached && prfind(bname,foilbkgs[b])){
	bkint[b]->Add(bkint[b],bgr[i]->Geth1(hname),1.,bgr[i]->norm);
	total->Add(bgr[i]->Geth1(hname),bgr[i]->norm);
	attached=true;
      }
    }
    if(!attached){
      mybak_ext->Add(mybak_ext,bgr[i]->Geth1(hname),1.,bgr[i]->norm);  
      total->Add(bgr[i]->Geth1(hname),bgr[i]->norm);
      attached=true;
    }
  }
  
  for(Int_t i=0;i<nsig;i++){
    mysig->Add(sig[i]->Geth1(hname),sig[i]->norm);
    total->Add(sig[i]->Geth1(hname),sig[i]->norm);
  }
  
  const Int_t MI=M;
  TMath::Sort(MI,events,index);
  
  ostringstream ossTitle;
  if(!tex) ossTitle<<"NEMO-3 Preliminary Results";
  if(tex) ossTitle<<"NEMO-3 Preliminary Results for "<<isotope.c_str()<<"";
  mydat->SetTitle(ossTitle.str().c_str());
  if(! TITLES) mydat->SetTitle("");
  
  mydat->Rebin(rebin);
  setData(mydat,s_type);
  
  Double_t ymin=0.5, ymax=0;
  ymax=mydat->GetMaximum()+mydat->GetBinError(mydat->GetMaximumBin());
  if(type!=3) ymax=ymax*1.2;
  if(type==3) ymax=ymax*1.8;
  ymax=ymax*yscale;
  
  if(xmin!=0 || xmax!=0) mydat->SetAxisRange(xmin,xmax,"X");
  mydat->SetAxisRange(ymin,ymax,"Y");
  
  setStyle(mydat);
  if(hdraw[0]) mydat->Draw("pe");
  
  Int_t clr=3;
  if(nsig){
    mysig->Rebin(rebin);
    mysig->SetLineWidth(LW);
    if(xmin!=0 || xmax!=0) mysig->SetAxisRange(xmin,xmax,"X");
    mysig->SetAxisRange(ymin,ymax,"Y");
    if(hdraw[1]){
      mysig->SetLineColor(ncolor(clr));
      clr++;
      mysig->Draw("hist same");
    }  
  }
  for(Int_t m=0;m<MI;m++){
    myint[index[m]]->Rebin(rebin);
    myint[index[m]]->SetLineWidth(LW);
    if(xmin!=0 || xmax!=0) myint[index[m]]->SetAxisRange(xmin,xmax,"X");
    myint[index[m]]->SetAxisRange(ymin,ymax,"Y");
    if(hdraw[2]){
      myint[index[m]]->SetLineColor(ncolor(clr));
      clr++;
      myint[index[m]]->Draw("hist same");
    }
  }  
  for(Int_t b=0;b<BI;b++){
    bkint[b]->Rebin(rebin);
    bkint[b]->SetLineWidth(LW);
    if(xmin!=0 || xmax!=0) bkint[b]->SetAxisRange(xmin,xmax,"X");
    bkint[b]->SetAxisRange(ymin,ymax,"Y");
    if(hdraw[3]){
      bkint[b]->SetLineColor(ncolor(clr));
      clr++;
      bkint[b]->Draw("hist same");
    }
  }
  mybak_ext->Rebin(rebin);
  mybak_ext->SetLineWidth(LW);
  if(xmin!=0 || xmax!=0) mybak_ext->SetAxisRange(xmin,xmax,"X");
  mybak_ext->SetAxisRange(ymin,ymax,"Y");
  if(hdraw[4]){
    mybak_ext->SetLineColor(ncolor(clr));
    clr++;
    mybak_ext->Draw("hist same");
  }
  total->Rebin(rebin);
  total->SetLineColor(12);
  total->SetLineWidth(LW);
  if(xmin!=0 || xmax!=0) total->SetAxisRange(xmin,xmax,"X");
  total->SetAxisRange(ymin,ymax,"Y");
  if(hdraw[5]) total->Draw("hist same");
  if(hdraw[0]) mydat->Draw("pe same");
  
  //Double_t chi2=chi2Calc(mydat,total);
  //Int_t ndf=ndfCalc(mydat);
  Double_t chi2=chi2Range(mydat,total,xmin,xmax);
  Int_t ndf=ndfRange(mydat,total,xmin,xmax);
  
  Double_t base = TMath::FloorNint(TMath::Log10(chi2)) 
    + TMath::FloorNint(TMath::Log10(ndf)) + 2.;
  Double_t LX=LX1;
  LX-=(LCS*base);
  Double_t LY=LY2;
  if(LEG_HEAD) LY-=LLS;
  if(hdraw[0]) LY-=LLS;
  if(hdraw[1]) LY-=LLS;
  if(hdraw[2]) LY-=(LLS*MI);
  if(hdraw[3]) LY-=(LLS*BI);
  if(hdraw[4]) LY-=LLS;
  if(hdraw[5]) LY-=LLS;
  if(LY<0.1) LY=0.1;
  
  cout<<"CHI2 / NDF = "<<setprecision(1)<<chi2<<" / "
      <<setprecision(0)<<ndf<<endl;
  ostringstream chindf;
  chindf<<"chi2/ndf = "<<fixed<<setprecision(1)<<chi2<<" / "
	<<setprecision(0)<<ndf;
  //ostringstream totchi;
  //totchi<<"Total ("<<fixed<<setprecision(1)<<chi2<<" / "
  //     <<setprecision(0)<<ndf<<")";
  
  TLegend *leg;
  leg=new TLegend(LX,LY,LX2,LY2);
  setStyle(leg);
  if(LEG_HEAD) leg->SetHeader(header.c_str());
  if(hdraw[0]) leg->AddEntry(mydat,"data","lpe");
  if(nsig && hdraw[1]) leg->AddEntry(mysig,"signal","l");
  for(Int_t m=0;m<MI;m++){
    if(hdraw[2]) leg->AddEntry(myint[index[m]],(myname[index[m]]).c_str(),"l");
  }
  for(Int_t b=0;b<BI;b++){
    if(hdraw[3]) leg->AddEntry(bkint[b],("int "+match_bb(foilbkgs[b])).c_str(),"l");
  }
  if(hdraw[4]) leg->AddEntry(mybak_ext,"externals","l");
  if(hdraw[5]) leg->AddEntry(total,"Total","l");
  //if(hdraw[5]) leg->AddEntry(total,totchi.str().c_str(),"l");
  leg->AddEntry(total,chindf.str().c_str(),"l");
  leg->Draw("same");
  
}


void DrawBkgConf(const char *hname,ana10 *data,bana10 **bgr,Int_t nbgr,
		 bana10 **sig,Int_t nsig,Double_t xmin,Double_t xmax,
		 Int_t rebin,string s_type,Double_t yscale,Int_t notused,
		 Int_t order,string isotope){
  
  //**************************************************************************
  //    THIS IS FOR DRAWING ADDED BACKGROUND CONTRIBUTIONS
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // 
  // This will draw the backgrounds added up with the data and the signal from 
  // your fit. This format has been BLESSED by the nemo collaboration to use.
  // 
  // (options)
  // xmin      : X-axis min value
  // xmax      : X-axis max value
  // rebin     : Hist rebinning
  // s_type    : Distribution type, "ee", "e", "cos", "eg", "g"
  // yscale    : Multiply the ymax by this value to adjust the scale
  // order     : Changes the order of the hist stacking, either "0" or "1" 
  //             flipps the internal and external order, signal is always 
  //             stacked on the top.
  // isotope   : To name your isotope in Latex for the Legend, ie "^{96}Zr",
  //             "^{48}Ca", "^{150}Nd", "^{116}Cd" ect ...
  // 
  //**************************************************************************
  
  Int_t tex=0;
  if(isotope!="") tex=1;
  
  Int_t type=axisType(s_type);
    
  TH1D *mydat=(TH1D*)data->Geth1(hname)->Clone("datas"); 
  TH1D *mybak_ext=(TH1D*)bgr[0]->Geth1(hname)->Clone("ext"); 
  TH1D *mybak_int=(TH1D*)bgr[0]->Geth1(hname)->Clone("int"); 
  TH1D *mysig;
  if(nsig) mysig=(TH1D*)sig[0]->Geth1(hname)->Clone("mysig"); 
  TH1D *total=(TH1D*)bgr[0]->Geth1(hname)->Clone("totals"); 
  
  mybak_ext->Reset();
  mybak_int->Reset();
  if(nsig) mysig->Reset();
  total->Reset();
  
  mydat->Sumw2();
  mybak_ext->Sumw2();
  mybak_int->Sumw2();
  if(nsig) mysig->Sumw2();
  total->Sumw2();
  
  // deal with the datas
  Double_t dat_events=0,dat_error=0;    
  for(Int_t i=0;i<=mydat->GetNbinsX();i++){
    if(mydat->GetBinContent(i)>0) dat_events+=mydat->GetBinContent(i);
    dat_error+=TMath::Power(mydat->GetBinError(i),2);
  }
  dat_error=TMath::Sqrt(dat_error);
  
  // split up the backgrounds
  Double_t bak_ext_evn=0,bak_int_evn=0,bak_gas_evn=0;
  Double_t bak_ext_err=0,bak_int_err=0,bak_gas_err=0;
  Bool_t attached=false;
  for(Int_t i=0;i<nbgr;i++){
    string bname = bgr[i]->Getname();
    attached=false;
    
    if(!attached && (sfind(bname,"cu")||sfind(bname,"te")||sfind(bname,"zr")
		     ||sfind(bname,"nd")||sfind(bname,"ca")||sfind(bname,"cd")
		     ||sfind(bname,"se")||sfind(bname,"mo"))){
      //cout<<"INT = "<<bname<<endl;
      mybak_int->Add(bgr[i]->Geth1(hname),bgr[i]->norm);
      bak_int_evn+=bgr[i]->GetExpected(hname);
      bak_int_err+=TMath::Power(bgr[i]->GetErrorExpected(hname),2);
      attached=true;
    }
    else if(!attached && (sfind(bname,"ex")||sfind(bname,"pm")
			  ||sfind(bname,"co")||sfind(bname,"wire")
			  ||sfind(bname,"sc")||sfind(bname,"sf")
			  ||sfind(bname,"sw")||sfind(bname,"gas"))){
      //cout<<"EXT = "<<bname<<endl;
      mybak_ext->Add(bgr[i]->Geth1(hname),bgr[i]->norm);
      bak_ext_evn+=bgr[i]->GetExpected(hname);
      bak_ext_err+=TMath::Power(bgr[i]->GetErrorExpected(hname),2);
      attached=true;
    }
    else if(!attached){
      cout<<"\n DrawBkgConf WARNING: component "<<bname
	  <<" was not recognised as being internal, external or gas \n\n"<<endl;
    }
  }
  
  bak_ext_err=TMath::Sqrt(bak_ext_err);
  bak_int_err=TMath::Sqrt(bak_int_err);
  
  // deal with the signals
  Double_t sig_events=0,sig_error=0;
  if(nsig){
    for(Int_t i=0;i<nsig;i++){
      mysig->Add(sig[i]->Geth1(hname),sig[i]->norm);
      sig_events+=sig[i]->GetExpected(hname);
      sig_error+=TMath::Power(sig[i]->GetErrorExpected(hname),2);
    }
    sig_error=TMath::Sqrt(sig_error);
  }
  
  // sum up total events from bkgs and sigs
  Double_t tot_events=0,tot_error=0;
  if(nsig){
    for(Int_t i=0;i<nsig;i++){
      total->Add(sig[i]->Geth1(hname),sig[i]->norm);
      tot_events+=sig[i]->GetExpected(hname);
      tot_error+=TMath::Power(sig[i]->GetErrorExpected(hname),2);
    }
  }
  for(Int_t i=0;i<nbgr;i++){
    total->Add(bgr[i]->Geth1(hname),bgr[i]->norm);
    tot_events+=bgr[i]->GetExpected(hname);
    tot_error+=TMath::Power(bgr[i]->GetErrorExpected(hname),2);
  }
  tot_error=TMath::Sqrt(tot_error);
  
  mydat->Rebin(rebin);
  mybak_ext->Rebin(rebin);
  mybak_int->Rebin(rebin);
  if(nsig) mysig->Rebin(rebin);
  total->Rebin(rebin);
  
  ostringstream ossTitle;
  if(!tex) ossTitle<<"NEMO-3 Preliminary Results";
  if(tex) ossTitle<<"NEMO-3 Preliminary Results for "<<isotope.c_str()<<"";
  mydat->SetTitle(ossTitle.str().c_str());
  if(! TITLES) mydat->SetTitle("");
  
  setData(mydat,s_type);
  
  Double_t ymin=0.5, ymax=0;
  ymax=mydat->GetMaximum()+mydat->GetBinError(mydat->GetMaximumBin());
  if(type!=3) ymax=ymax*1.2;
  if(type==3) ymax=ymax*1.8;
  ymax=ymax*yscale;
  
  if(xmin!=0 || xmax!=0) mydat->SetAxisRange(xmin,xmax,"X");
  mydat->SetAxisRange(ymin,ymax,"Y");
  
  setStyle(mydat);
  mydat->Draw("pe");
  
  // external settings
  mybak_ext->SetLineColor(kGreen-3);
  mybak_ext->SetLineWidth(LW);
  mybak_ext->SetFillStyle(3004);
  mybak_ext->SetFillColor(kGreen-3);
  if(xmin!=0 || xmax!=0) mybak_ext->SetAxisRange(xmin,xmax,"X");
  mybak_ext->SetAxisRange(ymin,ymax,"Y");
  
  // internal settings
  mybak_int->SetLineColor(kRed-3);
  mybak_int->SetLineWidth(LW);
  mybak_int->SetFillStyle(3005);
  mybak_int->SetFillColor(kRed-3);
  if(xmin!=0 || xmax!=0) mybak_int->SetAxisRange(xmin,xmax,"X");
  mybak_int->SetAxisRange(ymin,ymax,"Y");
  
  // signal settings
  if(nsig){
    mysig->SetLineColor(kBlue-3);
    mysig->SetLineWidth(LW);
    mysig->SetFillStyle(3006);
    mysig->SetFillColor(kBlue-3);
    if(xmin!=0 || xmax!=0) mysig->SetAxisRange(xmin,xmax,"X");
    mysig->SetAxisRange(ymin,ymax,"Y");
  }
  
  THStack *st=new THStack("st","st");    
  if(order){
    st->Add(mybak_ext);
    st->Add(mybak_int);
    if(nsig) st->Add(mysig);
  }else{
    st->Add(mybak_int);
    st->Add(mybak_ext);
    if(nsig) st->Add(mysig);
  }
  st->Draw("hist same");
  
  total->SetLineColor(12);
  total->SetLineWidth(LW);
  if(xmin!=0 || xmax!=0) total->SetAxisRange(xmin,xmax,"X");
  total->SetAxisRange(ymin,ymax,"Y");
  
  total->Draw("hist same");
  mydat->Draw("pe same");
  
  ostringstream oss1,oss2,oss3,oss4,oss5,oss6;
  oss1<<"data = "<<round(dat_events)<<"";
  if(!tex && nsig) oss5<<"signal = "<<round(sig_events)
		       <<" #pm "<<round(sig_error)<<"";
  if(tex && nsig) oss5<<isotope.c_str()<<" 2#nu#beta#beta = "
		      <<round(sig_events)<<" #pm "<<round(sig_error)<<"";
  oss4<<"int bkgs = "<<round(bak_int_evn)<<" #pm "<<round(bak_int_err)<<"";
  oss2<<"ext bkgs = "<<round(bak_ext_evn)<<" #pm "<<round(bak_ext_err)<<"";
  oss6<<"Total = "<<round(tot_events)<<" #pm "<<round(tot_error)<<"";
  
  //Double_t chi2=chi2Calc(mydat,total);
  //Int_t ndf=ndfCalc(mydat);
  Double_t chi2=chi2Range(mydat,total,xmin,xmax);
  Int_t ndf=ndfRange(mydat,total,xmin,xmax);
  
  Double_t base = TMath::FloorNint(TMath::Log10(tot_events))
    + TMath::FloorNint(TMath::Log10(tot_error)) + 2;
  Double_t LX=LX1;
  LX-=(LCS*base);
  if(tex) LX-=(LCS*2.);
  
  Double_t LY=LY2;
  if(LEG_HEAD) LY-=LLS;
  LY-=(LLS*4.); // data + ints + exts + chi2
  if(nsig) LY-=(LLS*nsig);
  if(LY<0.0) LY=0.0;
  
  cout<<"CHI2 / NDF = "<<setprecision(1)<<chi2<<" / "
      <<setprecision(0)<<ndf<<endl;
  ostringstream chindf;
  chindf<<"chi2/ndf = "<<fixed<<setprecision(1)<<chi2<<" / "
	<<setprecision(0)<<ndf;
  
  TLegend *leg;
  leg=new TLegend(LX,LY,LX2,LY2);
  setStyle(leg);
  if(LEG_HEAD) leg->SetHeader(header.c_str());
  leg->AddEntry(mydat,oss1.str().c_str(),"lpe");
  if(nsig) leg->AddEntry(mysig,oss5.str().c_str(),"f");
  leg->AddEntry(mybak_int,oss4.str().c_str(),"f");
  leg->AddEntry(mybak_ext,oss2.str().c_str(),"f");
  leg->AddEntry(total,oss6.str().c_str(),"l");
  leg->AddEntry(total,chindf.str().c_str(),"l");
  leg->Draw("same");
  
}


void DrawBkgAdd(const char *hname,ana10 *data,bana10 **bgr,Int_t nbgr,
		bana10 **sig,Int_t nsig,Double_t xmin,Double_t xmax,
		Int_t rebin,string s_type,string isotope){
  
  //**************************************************************************
  // THIS IS FOR DRAWING BACKGROUND ADDED HISTOGRAMS
  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // 
  // (flags)
  // s_type: "ee", "e", "cos", "eg", "g"
  //  
  // isotope: "^{96}Zr", "^{48}Ca", "^{150}Nd",  etc...
  //
  // this will draw the backgrounds added up with the data and the signal from 
  // your fit. this format has been BLESSED by the nemo collaboration to use.
  // 
  //**************************************************************************
  
  if(!nsig){
    std::cout<<"\n\n - NEED A SIGNAL - \n\n"<<std::endl;
    return;  
  }
  
  Int_t type=axisType(s_type);
  Int_t tex=0;
  if(isotope!="") tex=1;
    
  Int_t fillstyle=3003;
  
  TH1D *mydat=(TH1D*)data->Geth1(hname)->Clone("datas"); 
  TH1D *mybak=(TH1D*)bgr[0]->Geth1(hname)->Clone("bkgs added"); 
  TH1D *mysig=(TH1D*)sig[0]->Geth1(hname)->Clone("2b2n signal"); 
  TH1D *total=(TH1D*)sig[0]->Geth1(hname)->Clone("total mc"); 
  
  mybak->Reset();
  mysig->Reset();
  total->Reset();
  
  mydat->Sumw2();
  mybak->Sumw2();
  mysig->Sumw2();
  total->Sumw2();
  
  // deal with the datas
  Double_t dat_events=0,dat_error=0;    
  for(Int_t i=0;i<=mydat->GetNbinsX();i++){
    //if(mydat->GetBinContent(i)>0) dat_events+=mydat->GetBinContent(i);
    dat_events+=mydat->GetBinContent(i);
    dat_error+=TMath::Power(mydat->GetBinError(i),2);
  }
  dat_error=TMath::Sqrt(dat_error);
    
  // deal with the backgrounds
  Double_t bak_events=0,bak_error=0;
  for(Int_t i=0;i<nbgr;i++){
    mybak->Add(bgr[i]->Geth1(hname),bgr[i]->norm);
    bak_events+=bgr[i]->GetExpected(hname);
    bak_error+=TMath::Power(bgr[i]->GetErrorExpected(hname),2);
  }
  bak_error=TMath::Sqrt(bak_error);
  
  // deal with the signals
  Double_t sig_events=0,sig_error=0;
  for(Int_t i=0;i<nsig;i++){
    mysig->Add(sig[i]->Geth1(hname),sig[i]->norm);
    sig_events+=sig[i]->GetExpected(hname);
    sig_error+=TMath::Power(sig[i]->GetErrorExpected(hname),2);
  }
  sig_error=TMath::Sqrt(sig_error);
  
  // sum up total events from bkgs and sigs
  Double_t tot_events=0,tot_error=0;
  for(Int_t i=0;i<nsig;i++){
    total->Add(sig[i]->Geth1(hname),sig[i]->norm);
    tot_events+=sig[i]->GetExpected(hname);
    tot_error+=TMath::Power(sig[i]->GetErrorExpected(hname),2);
  }
  for(Int_t i=0;i<nbgr;i++){
    total->Add(bgr[i]->Geth1(hname),bgr[i]->norm);
    tot_events+=bgr[i]->GetExpected(hname);
    tot_error+=TMath::Power(bgr[i]->GetErrorExpected(hname),2);
  }
  tot_error=TMath::Sqrt(tot_error);
  
  mydat->Rebin(rebin);
  mybak->Rebin(rebin);
  mysig->Rebin(rebin);
  total->Rebin(rebin);
  
  ostringstream ossTitle;
  if(!tex) ossTitle<<"NEMO-3 Preliminary Results";
  if(tex) ossTitle<<"NEMO-3 Preliminary Results for "<<isotope.c_str()<<"";
  mydat->SetTitle(ossTitle.str().c_str());
  if(! TITLES) mydat->SetTitle("");
    
  setData(mydat,s_type);
  
  Double_t ymin=0.5, ymax=0;
  ymax=mydat->GetMaximum()+mydat->GetBinError(mydat->GetMaximumBin());
  if(type!=3) ymax=ymax*1.2;
  if(type==3) ymax=ymax*1.2;
  
  if(xmin!=0 || xmax!=0) mydat->SetAxisRange(xmin,xmax,"X");
  mydat->SetAxisRange(ymin,ymax,"Y");
  
  setStyle(mydat);
  mydat->Draw("pe");
    
  //mybak->SetFillStyle(fillstyle);
  mybak->SetFillStyle(3004);
  mybak->SetLineColor(kRed-3);
  mybak->SetLineWidth(LW);
  mybak->SetFillColor(kRed-3);
  if(xmin!=0 || xmax!=0) mybak->SetAxisRange(xmin,xmax,"X");
  mybak->SetAxisRange(ymin,ymax,"Y");
  mybak->Draw("hist same");
  
  //mysig->SetFillStyle(fillstyle);
  mysig->SetFillStyle(3005);
  mysig->SetLineColor(kBlue-3);
  mysig->SetLineWidth(LW);
  mysig->SetFillColor(kBlue-3);
  if(xmin!=0 || xmax!=0) mysig->SetAxisRange(xmin,xmax,"X");
  mysig->SetAxisRange(ymin,ymax,"Y");
  mysig->Draw("hist same");
  
  total->SetFillStyle(1001);
  total->SetLineColor(12);
  total->SetLineWidth(LW);
  if(xmin!=0 || xmax!=0) total->SetAxisRange(xmin,xmax,"X");
  total->SetAxisRange(ymin,ymax,"Y");
  total->Draw("hist same");
  
  mydat->Draw("pe same");
  
  ostringstream oss1,oss2,oss3,oss4;
  oss1<<"data = "<<round(dat_events)<<"";
  oss2<<"sum of bkgs = "<<round(bak_events)<<" #pm "<<round(bak_error)<<"";
  if(!tex) oss3<<"signal = "<<round(sig_events)<<" #pm "
	       <<round(sig_error)<<"";
  if(tex) oss3<<isotope.c_str()<<" 2#nu#beta#beta = "
	      <<round(sig_events)<<" #pm "<<round(sig_error)<<"";
  oss4<<"Total = "<<round(tot_events)<<" #pm "<<round(tot_error)<<"";
  
  Double_t chi2=chi2Calc(mydat,total);
  Int_t ndf=ndfCalc(mydat);
  //Double_t chi2=chi2Range(mydat,total,xmin,xmax);
  //Int_t ndf=ndfRange(mydat,total,xmin,xmax);
  
  Double_t base = TMath::FloorNint(TMath::Log10(tot_events))
    + TMath::FloorNint(TMath::Log10(tot_error)) + 2;
  Double_t LX=LX1;
  LX-=(LCS*base);
  if(tex) LX-=(LCS*2.);
    
  Double_t LY=LY2;
  if(LEG_HEAD) LY-=LLS;
  LY-=(LLS*3.); // data + bkgs + chi2
  if(nsig) LY-=(LLS*nsig);
  if(LY<0.0) LY=0.0;
  
  cout<<"CHI2 / NDF = "<<setprecision(1)<<chi2<<" / "
      <<setprecision(0)<<ndf<<endl;
  ostringstream chindf;
  chindf<<"chi2/ndf = "<<fixed<<setprecision(1)<<chi2<<" / "
	<<setprecision(0)<<ndf;
    
  TLegend *leg;
  leg=new TLegend(LX,LY,LX2,LY2);
  setStyle(leg);
  if(LEG_HEAD) leg->SetHeader(header.c_str());
  leg->AddEntry(mydat,oss1.str().c_str(),"lpe");
  leg->AddEntry(mysig,oss3.str().c_str(),"f");
  leg->AddEntry(mybak,oss2.str().c_str(),"f");
  leg->AddEntry(total,oss4.str().c_str(),"l");
  leg->AddEntry(total,chindf.str().c_str(),"l");
  leg->Draw("same");
  
}


void DrawBkgSub(const char *hname,ana10 *data,bana10 **bgr,Int_t nbgr,
		bana10 **sig,Int_t nsig,Double_t xmin,Double_t xmax,
		Int_t rebin,string s_type,string isotope){

  //**************************************************************************
  // THIS IS FOR DRAWING BACKGROUND SUBTRACTED HISTOGRAMS
  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // 
  // (flags)
  // s_type: "ee", "e", "cos", "eg", "g"
  //  
  // isotope: "^{96}Zr", "^{48}Ca", "^{150}Nd",  etc...
  //
  // this will draw the backgrounds added up with the data and the signal from 
  // your fit. this format has been BLESSED by the nemo collaboration to use.
  // 
  //**************************************************************************
  
  if(!nsig){
    std::cout<<"\n\n - NEED A SIGNAL - \n\n"<<std::endl;
    return;  
  }
  
  Int_t type=axisType(s_type);
  Int_t tex=0;
  if(isotope!="") tex=1;
    
  Int_t fillstyle=3003;
  
  TH1D *mydat=(TH1D*)data->Geth1(hname)->Clone("datas"); 
  TH1D *mybak=(TH1D*)bgr[0]->Geth1(hname)->Clone("bkgs added"); 
  TH1D *mysig=(TH1D*)sig[0]->Geth1(hname)->Clone("2b2n signal"); 
  //TH1D *total=(TH1D*)sig[0]->Geth1(hname)->Clone("total mc"); 
  
  mybak->Reset();
  mysig->Reset();
  //total->Reset();
  
  mydat->Sumw2();
  mybak->Sumw2();
  mysig->Sumw2();
  //total->Sumw2();
  
  // deal with the datas
  Double_t dat_events=0,dat_error=0;    
  for(Int_t i=1;i<=mydat->GetNbinsX();i++){
    if(mydat->GetBinContent(i)>0) dat_events+=mydat->GetBinContent(i);
    dat_error+=TMath::Power(mydat->GetBinError(i),2);
  }
  dat_error=TMath::Sqrt(dat_error);
  
  // deal with the backgrounds
  Double_t bak_events=0,bak_error=0;
  for(Int_t i=0;i<nbgr;i++){
    mybak->Add(bgr[i]->Geth1(hname),bgr[i]->norm);
    bak_events+=bgr[i]->GetExpected(hname);
    bak_error+=TMath::Power(bgr[i]->GetErrorExpected(hname),2);
  }
  bak_error=TMath::Sqrt(bak_error);
  
  // deal with the signals
  Double_t sig_events=0,sig_error=0;
  for(Int_t i=0;i<nsig;i++){
    mysig->Add(sig[i]->Geth1(hname),sig[i]->norm);
    sig_events+=sig[i]->GetExpected(hname);
    sig_error+=TMath::Power(sig[i]->GetErrorExpected(hname),2);
  }
  sig_error=TMath::Sqrt(sig_error);
  
  // sum up total events from bkgs and sigs
  /*
  Double_t tot_events=0,tot_error=0;
  for(Int_t i=0;i<nsig;i++){
    total->Add(total,sig[i]->Geth1(hname),1.,sig[i]->norm);
    tot_events+=sig[i]->GetExpected(hname);
    tot_error+=TMath::Power(sig[i]->GetErrorExpected(hname),2);
  }
  for(Int_t i=0;i<nbgr;i++){
    total->Add(total,bgr[i]->Geth1(hname),1.,bgr[i]->norm);
    tot_events+=bgr[i]->GetExpected(hname);
    tot_error+=TMath::Power(bgr[i]->GetErrorExpected(hname),2);
  }
  tot_error=TMath::Sqrt(tot_error);
  */
  
  mydat->Rebin(rebin);
  mybak->Rebin(rebin);
  mysig->Rebin(rebin);
  //total->Rebin(rebin);
  
  ostringstream ossTitle;
  if(!tex) ossTitle<<"NEMO-3 Preliminary Results";
  if(tex) ossTitle<<"NEMO-3 Preliminary Results for "<<isotope.c_str()<<"";
  mydat->SetTitle(ossTitle.str().c_str());
  if(! TITLES) mydat->SetTitle("");
  
  setData(mydat,s_type);
  mydat->Add(mybak,-1.);
  
  // Force errors to be sqrt(N) in each bin
  ///////////////////////////////////////////////////////////
  for(Int_t i=1; i<=mysig->GetNbinsX(); i++){
    if(mysig->GetBinContent(i)<0) mysig->SetBinContent(i,0);
    mysig->SetBinError(i,sqrt(mysig->GetBinContent(i)));
  }
  
  Double_t ymin=0.5, ymax=0;
  ymax=mydat->GetMaximum()+mydat->GetBinError(mydat->GetMaximumBin());
  if(type!=3) ymax=ymax*1.2;
  if(type==3) ymax=ymax*1.2;
  
  if(xmin!=0 || xmax!=0) mydat->SetAxisRange(xmin,xmax,"X");
  mydat->SetAxisRange(ymin,ymax,"Y");
  
  setStyle(mydat);
  mydat->Draw("pe");
  
  //mysig->SetFillStyle(fillstyle);
  mysig->SetFillStyle(3005);
  mysig->SetLineColor(kBlue-3);
  mysig->SetLineWidth(LW);
  mysig->SetFillColor(kBlue-3);
  if(xmin!=0 || xmax!=0) mysig->SetAxisRange(xmin,xmax,"X");
  mysig->SetAxisRange(ymin,ymax,"Y");
  mysig->Draw("hist same");
  
  mydat->Draw("pe same");
  
  ostringstream oss1,oss2,oss3,oss4;
  oss1<<"data - bkgs = "<<round(dat_events-bak_events)<<" #pm "
      <<round(TMath::Sqrt(dat_events+(bak_error*bak_error)))<<"";
  if(!tex) oss3<<"signal = "<<round(sig_events)<<" #pm "
	       <<round(sig_error)<<"";
  if(tex) oss3<<isotope.c_str()<<" 2#nu#beta#beta = "
	      <<round(sig_events)<<" #pm "<<round(sig_error)<<"";
  
  //Double_t chi2=chi2Man(mydat,mysig);
  //Double_t chi2=chi2Calc(mydat,mysig);
  Double_t chi2=mydat->Chi2Test(mysig,"UUCHI2");
  //Int_t ndf=ndfCalc(mydat);
  //Double_t chi2=chi2Range(mydat,mysig,xmin,xmax);
  Int_t ndf=ndfRange(mydat,mysig,xmin,xmax);
  
  Double_t base = TMath::FloorNint(TMath::Log10(sig_events))
    + TMath::FloorNint(TMath::Log10(sig_error)) + 2;
  Double_t LX=LX1;
  LX-=(LCS*base);
  if(tex) LX-=(LCS*2.);
    
  Double_t LY=LY2;
  if(LEG_HEAD) LY-=LLS;
  LY-=(LLS*2.); // data + chi2
  if(nsig) LY-=(LLS*nsig);
  if(LY<0.0) LY=0.0;
  
  cout<<"CHI2 / NDF = "<<setprecision(1)<<chi2<<" / "
      <<setprecision(0)<<ndf<<endl;
  ostringstream chindf;
  chindf<<"chi2/ndf = "<<fixed<<setprecision(1)<<chi2<<" / "
	<<setprecision(0)<<ndf;
  
  TLegend *leg;
  leg=new TLegend(LX,LY,LX2,LY2);
  setStyle(leg);
  if(LEG_HEAD) leg->SetHeader(header.c_str());
  leg->AddEntry(mydat,oss1.str().c_str(),"lpe");
  leg->AddEntry(mysig,oss3.str().c_str(),"f");
  leg->AddEntry(mysig,chindf.str().c_str(),"l");
  leg->Draw("same");
  
}


void DrawLim(const char *hname,ana10 *data,bana10 **bgr,Int_t nbgr,
	     bana10 **sig,Int_t nsig,bana10 **lim,Int_t nlim,Double_t xmin,
	     Double_t xmax,Int_t rebin,Double_t limEvents,string s_type,
	     string isotope,string decaymode){

  //**************************************************************************
  // THIS IS FOR DRAWING BACKGROUND ADDED HISTOGRAMS WITH LIMITS
  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // 
  // (flags)
  // 
  // limEvents: the number of signal events at 90% CL
  // 
  // s_type: to set axis lables, ("ee", "e", "cos", "eg", "g")
  //  
  // isotope: legend lable in latex, ("^{96}Zr", "^{48}Ca", "^{150}Nd",
  //                                  "^{116}Cd" ...)
  // 
  // decaymode: limit you are producing, ("m1", "m2", "m6" ...)
  // 
  // this will draw the backgrounds added up with the data and the signal from 
  // your fit and the limit MC. this format has not yet been blessed by the 
  // nemo collaboration to use.
  // 
  //**************************************************************************
  
  if(!nsig>=1){
    std::cout<<"\n\n - NEED A SIGNAL - \n\n"<<std::endl;
    return;  
  }
  
  if(!nlim>=1){
    std::cout<<"\n\n - NEED LIMIT MC - \n\n"<<std::endl;
    return;  
  }
  
  Int_t type=axisType(s_type);
  Int_t tex=0,tex2=0;
  if(isotope!="") tex=1;
  if(decaymode!="") tex2=1;
  
  Int_t fillstyle=3003;
  
  TH1D *mydat=(TH1D*)data->Geth1(hname)->Clone("datas"); 
  TH1D *mybak=(TH1D*)bgr[0]->Geth1(hname)->Clone("bkgs_mc"); 
  TH1D *mysig=(TH1D*)sig[0]->Geth1(hname)->Clone("signal_mc"); 
  TH1D *mylim=(TH1D*)lim[0]->Geth1(hname)->Clone("limit_mc"); 
  TH1D *total=(TH1D*)bgr[0]->Geth1(hname)->Clone("total_mc"); 
  
  mybak->Reset();
  mysig->Reset();
  //mylim->Reset();
  total->Reset();
  
  mydat->Sumw2();
  mybak->Sumw2();
  mysig->Sumw2();
  total->Sumw2();
  
  // deal with the datas
  Double_t dat_events=0,dat_error=0;    
  for(Int_t i=0;i<=mydat->GetNbinsX();i++){
    if(mydat->GetBinContent(i)>0) dat_events+=mydat->GetBinContent(i);
    dat_error+=TMath::Power(mydat->GetBinError(i),2);
  }
  dat_error=TMath::Sqrt(dat_error);
  
  // deal with the backgrounds
  Double_t bak_events=0,bak_error=0;
  for(Int_t i=0;i<nbgr;i++){
    mybak->Add(bgr[i]->Geth1(hname),bgr[i]->norm);
    bak_events+=bgr[i]->GetExpected(hname);
    bak_error+=TMath::Power(bgr[i]->GetErrorExpected(hname),2);
  }
  bak_error=TMath::Sqrt(bak_error);
  
  // deal with the signals
  Double_t sig_events=0,sig_error=0;
  for(Int_t i=0;i<nsig;i++){
    mysig->Add(sig[i]->Geth1(hname),sig[i]->norm);
    sig_events+=sig[i]->GetExpected(hname);
    sig_error+=TMath::Power(sig[i]->GetErrorExpected(hname),2);
  }
  sig_error=TMath::Sqrt(sig_error);
  
  // deal with the limits
  Double_t lim_events=0;
  if(limEvents!=-1) mylim->Scale(limEvents/mylim->Integral());
  lim_events=mylim->Integral();
  
  // sum up total events from bkgs and sigs
  Double_t tot_events=0,tot_error=0;
  for(Int_t i=0;i<nsig;i++){
    total->Add(sig[i]->Geth1(hname),sig[i]->norm);
    tot_events+=sig[i]->GetExpected(hname);
    tot_error+=TMath::Power(sig[i]->GetErrorExpected(hname),2);
  }
  for(Int_t i=0;i<nbgr;i++){
    total->Add(bgr[i]->Geth1(hname),bgr[i]->norm);
    tot_events+=bgr[i]->GetExpected(hname);
    tot_error+=TMath::Power(bgr[i]->GetErrorExpected(hname),2);
  }
  //total->Add(total,mylim,1.,1.);
  //tot_events+=mylim->Integral();
  tot_error=TMath::Sqrt(tot_error);
  
  mydat->Rebin(rebin);
  mybak->Rebin(rebin);
  mysig->Rebin(rebin);
  mylim->Rebin(rebin);
  total->Rebin(rebin);
  
  ostringstream ossTitle;
  if(!tex) ossTitle<<"NEMO-3 Preliminary Results";
  if(tex) ossTitle<<"NEMO-3 Preliminary Results for "<<isotope.c_str()<<"";
  mydat->SetTitle(ossTitle.str().c_str());
  if(! TITLES) mydat->SetTitle("");
  
  setData(mydat,s_type);
  
  const Double_t ymin=0.1;
  Double_t ymax=0;
  ymax=(mydat->GetMaximum()+mydat->GetBinError(mydat->GetMaximumBin()))*11.0;
  
  if(xmin!=0 || xmax!=0) mydat->SetAxisRange(xmin,xmax,"X");
  mydat->SetAxisRange(ymin,ymax,"Y");
  
  setStyle(mydat);
  mydat->Draw("pe");
  
  mybak->SetFillStyle(fillstyle);
  mybak->SetLineColor(kRed-3);
  mybak->SetLineWidth(LW);
  mybak->SetFillColor(kRed-3);
  if(xmin!=0 || xmax!=0) mybak->SetAxisRange(xmin,xmax,"X");
  mybak->SetAxisRange(ymin,ymax,"Y");
  mybak->Draw("hist same");
  
  mysig->SetFillStyle(fillstyle);
  mysig->SetLineColor(kBlue-3);
  mysig->SetLineWidth(LW);
  mysig->SetFillColor(kBlue-3);
  if(xmin!=0 || xmax!=0) mysig->SetAxisRange(xmin,xmax,"X");
  mysig->SetAxisRange(ymin,ymax,"Y");
  mysig->Draw("hist same");
  
  mylim->SetFillStyle(fillstyle);
  mylim->SetLineColor(kGreen-3);
  mylim->SetLineWidth(LW);
  mylim->SetFillColor(kGreen-3);
  if(xmin!=0 || xmax!=0) mylim->SetAxisRange(xmin,xmax,"X");
  mylim->SetAxisRange(ymin,ymax,"Y");
  mylim->Draw("hist same");
  
  total->SetFillStyle(1001);
  total->SetLineColor(12);
  total->SetLineWidth(LW);
  if(xmin!=0 || xmax!=0) total->SetAxisRange(xmin,xmax,"X");
  total->SetAxisRange(ymin,ymax,"Y");
  total->Draw("hist same");
  
  mydat->Draw("pe same");
  
  ostringstream oss1,oss2,oss3,oss4,oss5;
  oss1<<"data = "<<round(dat_events)<<"";
  oss2<<"sum of bkgs = "<<round(bak_events)<<" #pm "<<round(bak_error)<<"";
  if(!tex) oss3<<"signal = "<<round(sig_events)<<" #pm "
	       <<round(sig_error)<<"";
  if(tex) oss3<<isotope.c_str()<<" 2#nu#beta#beta = "
	      <<round(sig_events)<<" #pm "<<round(sig_error)<<"";
  oss4<<"0#nu limit < "<<fixed<<setprecision(2)<<lim_events<<"";
  oss5<<"Total = "<<round(tot_events)<<" #pm "<<round(tot_error)<<"";
  
  Double_t base = TMath::FloorNint(TMath::Log10(tot_events))
    + TMath::FloorNint(TMath::Log10(tot_error)) + 2;
  Double_t LX=LX1;
  LX-=(LCS*base);
  if(tex) LX-=(LCS*2.);
    
  Double_t LY=LY2;
  if(LEG_HEAD) LY-=LLS;
  LY-=(LLS*3.); // data + bkgs +lim
  if(nsig) LY-=(LLS*nsig);
  if(LY<0.0) LY=0.0;
  
  TLegend *leg;
  leg=new TLegend(LX,LY,LX2,LY2);
  setStyle(leg);
  if(LEG_HEAD) leg->SetHeader(header.c_str());
  leg->AddEntry(mydat,oss1.str().c_str(),"lpe");
  leg->AddEntry(mysig,oss3.str().c_str(),"f");
  leg->AddEntry(mybak,oss2.str().c_str(),"f");
  leg->AddEntry(mylim,oss4.str().c_str(),"f");
  leg->Draw("same");
  
}

