#include "ana.hpp"

#include <string>
#include <sstream>
#include <iomanip>
#include "../../utils/h10.h"
#include "TH2.h"
#include "TH1.h"
#include "TLegend.h"
#include "TFractionFitter.h"
#include "TMath.h"
#include "TLimit.h"
#include "TLimitDataSource.h"
#include "TConfidenceLevel.h"
#include "TPaveText.h"
#include "TLatex.h"


// colour scheme
//////////////////////////////////////////////////////////////////////////
Color_t palitra(Int_t n=0){
  switch(n){
  case 0:return  kWhite;     case 1:return  kBlack;   case 2:return  kGray+2;  
  case 3:return  kRed+1;     case 4:return  kBlue+1;  case 5:return  kGreen+1;
  case 6:return  kMagenta+2; case 7:return  kCyan+2;  case 8:return  kOrange+2;
  case 9:return  kPink+2;    case 10:return kAzure+2; case 11:return kSpring+2;
  case 12:return kRed-4;     case 13:return kBlue-4;  case 14:return kGreen-4;
  case 15:return kMagenta-4; case 16:return kCyan-4;  case 17:return kOrange-4;
  case 18:return kPink-4;    case 19:return kAzure-4; case 20:return kSpring-4;
  default:return 40;
  }
}
//////////////////////////////////////////////////////////////////////////


//Drawing utilities
///////////////////////////////////////////////////////////////////////////
void Drawall_d1(const char * hname, ana10 *data,bana10 **bgr, Int_t nbgr){
  // Draws hname from all bgr histos at the same time;
  gROOT->SetStyle("Plain"); // added by matt kauer
  TH1D *sum=(TH1D*)data->Geth1(hname)->Clone("bsum"); 
  TH1D *temp[1000]; 
  sum->Reset();
  //TLegend *leg= new TLegend(0.5,0.45,0.9,0.82);
  TLegend *leg= new TLegend(0.45,0.6,0.995,0.995);
  leg->SetFillColor(0);
  Int_t bcolor=3;
  data->Geth1(hname)->SetMinimum(0);
  data->Geth1(hname)->Draw("e");
  ostringstream oss2;
  oss2<<data->Getname()<<" ("<<data->Geth1(hname)->Integral()<<") ";
  
  leg->AddEntry(data->Geth1(hname),oss2.str().c_str(),"l");
  for(Int_t i=0;i<nbgr;i++){
    temp[i]=(TH1D* ) bgr[i]->Geth1(hname)->Clone(string("temp"+i).c_str());
    temp[i]->Add(temp[i],bgr[i]->Geth1(hname),0,bgr[i]->norm);
    temp[i]->SetTitleSize(0.06,"X");
    temp[i]->SetTitleSize(0.06,"Y");
    temp[i]->SetLineWidth(2);
    temp[i]->SetLineColor(palitra(i+bcolor));
    temp[i]->Draw("same");
    ostringstream oss;
    oss<<bgr[i]->Getname()<<" ("<<temp[i]->Integral()<<") A="<<bgr[i]->activity<<" ";
    if (temp[i]->Integral()>0) leg->AddEntry(temp[i],oss.str().c_str(),"l");
    sum->Add(sum,temp[i],1.,1.);
  }
    sum->SetTitleSize(0.06,"X");
    sum->SetTitleSize(0.06,"Y");
    sum->SetLineWidth(2);
    sum->SetLineColor(palitra(2));
    sum->Draw("hist,same"); 
    ostringstream oss1;
    oss1<<"Total MC ("<<sum->Integral()<<") ";
    leg->AddEntry(sum,oss1.str().c_str(),"l");
    leg->Draw();
}

//*******************************************************************
// This function Draws plot for hname histo
// It draws data, backgournds,sum of all backgrounds and legend
// Input parameters:
// hname -- name of the histogram to draw
// data,bgr,nbgr,sig,nsig -- pointers to data, background and signal components in the memory
// kLegend -- boolean flag to draw the legend or not (default = true)
// kSum -- boolean flag to draw one MC component on top of another. If false than MC components plots are superimposed on on another.(default=false)
// cumulative -- cumulative plot, each bin actually equal to sum of all higher bins in the original histogram. (default=false)
// blow -- low boundary for the plot. If 0 use histograms low boumdary (default=0)
// bup -- upper boundary for the plot. If 0 use histogram upper boundary (default=0)
// rebin -- rebin histogram (defoult = 1, no rebinning)
//*******************************************************************
void Drawall_d1(const char *hname, ana10 *data, bana10 **bgr, Int_t nbgr, bana10 ** sig, Int_t nsig, Int_t kLegend, Bool_t kSum, Bool_t cumulative, Double_t blow, Double_t bup, Int_t rebin, Int_t fillstyle, Int_t printnum, Double_t yscale){
  
  
  gROOT->SetStyle("Plain");
  
  TH1D *sum=(TH1D*)bgr[0]->Geth1(hname)->Clone("bsum"); 
  sum->Reset();
  sum->Sumw2();
  sum->Rebin(rebin);
  
  Float_t mcetot=0;
  Float_t mctot=0;
  TH1D *temp[1000]; 
  Float_t a[1000];
  Int_t index[1000];
  string legends[1000];
  
  Double_t lxmax=0.996;
  Double_t lymax=0.996;
  Double_t lymin=lymax-(printnum*0.04)-0.12;
  //if(printnum>10) lymin=0.59;
  TLegend *leg;
  if(kLegend==1) leg= new TLegend(0.35,lymin,lxmax,lymax);
  if(kLegend==2) leg= new TLegend(0.60,lymin,lxmax,lymax);
  leg->SetFillColor(0);
  leg->SetTextSize(0.035);
  
  Int_t bcolor=3;
  TH1D *datah;
  Float_t di;
  Int_t lx,ux; //lower and upper limit in bins

  if(data){
    if(cumulative){
      datah=(TH1D *)GetCumulativeDown(data->Geth1(hname)->Rebin(rebin,"newdata"),bup);
    }else{
      datah=(TH1D *)data->Geth1(hname)->Clone("newdata");
    }
    datah->Sumw2();
    if(blow!=0){
      Double_t ledges[10000];
      datah->GetLowEdge(ledges);
      lx=0;
      ux=0;
      for(Int_t k=0;k<=datah->GetNbinsX()+1;k++){
	if(blow<ledges[k] && lx==0) lx = k-1;
	if(bup<ledges[k] && ux==0) ux = k-1;
      }
      if(ux==0) ux=datah->GetNbinsX()+1;
      
    }else{
      lx=0;
      ux=datah->GetNbinsX()+1;
    }
    if(cumulative){
      di = datah->GetBinContent(lx);
    }else{
      di = datah->Integral(lx,ux);
      datah->Rebin(rebin);
    }
    if(blow) datah->GetXaxis()->SetRangeUser(blow,bup);
    
    datah->SetMinimum(0.01);
    datah->SetMarkerStyle(20);
    datah->SetMarkerSize(0.8);
    datah->SetLineWidth(1);
    datah->Draw("ep");
    
    ostringstream oss;
    oss<<data->Getname()<<" ("<<di<<") ";
  
    leg->AddEntry(datah,oss.str().c_str(),"lpe");
  }
  
  TH1D *sum2=(TH1D*)datah->Clone();
  sum2->Reset();
  sum2->Sumw2();
  
  for(Int_t i=0;i<nbgr;i++){
    if(cumulative){
      temp[i]=(TH1D *)GetCumulativeDown(bgr[i]->Geth1(hname)->Rebin(rebin,"newtemp"),bup);
    }else{
      temp[i]=(TH1D* ) bgr[i]->Geth1(hname)->Clone(string("temp"+i).c_str());
      temp[i]->Rebin(rebin);
    }
    temp[i]->Sumw2();
    //temp[i]->Add(temp[i],bgr[i]->Geth1(hname),0,bgr[i]->norm);
    temp[i]->Scale(bgr[i]->norm);
    temp[i]->SetTitleSize(0.06,"X");
    temp[i]->SetTitleSize(0.06,"Y");
    temp[i]->SetLineWidth(2);
    Float_t ni;
    ni = bgr[i]->GetExpected(hname,lx,ux);
    Float_t dni=0;
    if(bgr[i]->norm!=0 and ni!=0) dni=bgr[i]->GetErrorExpected(hname,lx,ux);
    mcetot+=dni*dni;
    mctot+=ni;
    //cout<<"Draw "<<bgr[i]->Getname()<<" A="<<bgr[i]->norm<<" ea="<<bgr[i]->eactivity<<" ni = "<<temp[i]->Integral()<<" dni="<<dni<<endl; 
    ostringstream oss1;
    if(kLegend==1) oss1<<bgr[i]->Getname()<<" "<<printenumber(ni,dni,2)<<" "<<"A="<<printenumber(bgr[i]->activity,bgr[i]->eactivity,2);
    if(kLegend==2) oss1<<bgr[i]->Getname();
    a[i]=ni;
    legends[i]=oss1.str();
    sum->Add(sum,temp[i],1.,1.);
  }
  Int_t j=nbgr;
  for(Int_t i=0;i<nsig;i++){
    if(cumulative){
      temp[j]=(TH1D *)GetCumulativeDown(sig[i]->Geth1(hname)->Rebin(rebin,"newtemp"),bup);      
    }else{
      temp[j]=(TH1D* ) sig[i]->Geth1(hname)->Clone(string("temp"+j).c_str());
      temp[j]->Rebin(rebin);
    }
    temp[j]->Sumw2();
    //temp[j]->Add(temp[j],sig[i]->Geth1(hname),0,sig[i]->norm);
    temp[j]->Scale(sig[i]->norm);
    temp[j]->SetTitleSize(0.06,"X");
    temp[j]->SetTitleSize(0.06,"Y");
    temp[j]->SetLineWidth(2);
    Float_t ni;
    ni = sig[i]->GetExpected(hname,lx,ux);
    Float_t dni=0;
    if(ni!=0 && sig[i]->activity!=0) dni=sig[i]->GetErrorExpected(hname,lx,ux);
    mcetot+=dni*dni;
    mctot+=ni;
    ostringstream oss2;
    if(kLegend==1) oss2<<sig[i]->Getname()<<" "<<printenumber(ni,dni,2)<<" \n"<<" A="<<printenumber(sig[i]->activity,sig[i]->eactivity,2);
    if(kLegend==2) oss2<<sig[i]->Getname();
    a[j]=ni;
    legends[j]=oss2.str();
    sum->Add(sum,temp[j],1.,1.);
    j++;
  }
  TMath::Sort(nbgr+nsig,a,index);
  if(!kSum){
    if(printnum>nbgr+nsig) printnum=nbgr+nsig;
    for(Int_t i=0;i<nbgr+nsig;i++){
      temp[index[i]]->SetLineColor(palitra(i+bcolor));
      //temp[index[i]]->Draw("hist same l");
      if(i<printnum && a[index[i]]>0){
	temp[index[i]]->Draw("hist same");
	leg->AddEntry(temp[index[i]],legends[index[i]].c_str(),"l");
      }else{
	if(a[index[i]]>0) sum2->Add(temp[index[i]]);
      }
    }
    sum2->SetLineColor(palitra(printnum+bcolor));
    sum2->SetLineWidth(2);
    //sum2->Draw("hist same");
    //leg->AddEntry(sum2,"sum of other bkgs","l");
  }
  
  if(kSum){
    if(printnum>nbgr+nsig) printnum=nbgr+nsig;
    for(Int_t i=nbgr+nsig-2;i>=0;i--){
      temp[index[i]]->Add(temp[index[i+1]]);
    }
    for(Int_t i=0;i<nbgr+nsig;i++){
      temp[index[i]]->SetLineColor(palitra(i+bcolor));
      temp[index[i]]->SetFillColor(palitra(i+bcolor));
      temp[index[i]]->SetFillStyle(fillstyle);
      //temp[index[i]]->Draw("hist same");
      if(i<printnum && a[index[i]]>0){
	temp[index[i]]->Draw("hist same");
	leg->AddEntry(temp[index[i]],legends[index[i]].c_str(),"f");
      }
    }
  }
  
  sum->SetTitleSize(0.06,"X");
  sum->SetTitleSize(0.06,"Y");
  sum->SetLineWidth(2);
  sum->SetLineColor(palitra(2));
  if(kSum && data){
    datah->Draw("ep same");
  }
  sum->Draw("hist same"); 
  datah->Draw("ep same");
  
  Double_t ymax=1.1*sum->GetMaximum();
  if(datah->GetMaximum() > sum->GetMaximum()) ymax=1.1*datah->GetMaximum();
  datah->SetAxisRange(0.1,ymax*yscale,"Y");
  
  
  /// Calcultae chi2/ndf for data/MC comparison
  //////////////////////////////////////////////////
  Float_t chi2=0;
  Int_t ndf=0;
  /*
  if(blow!=0){
    Double_t ledges[10000];
    datah->GetLowEdge(ledges);
    lx=0;
    ux=0;
    for(Int_t k=0;k<=datah->GetNbinsX()+1;k++){
      if(blow<ledges[k] && lx==0) lx = k-1;
      if(bup<ledges[k] && ux==0) ux = k-1;
    }
    if(ux==0) ux=datah->GetNbinsX()+1;
    
  }else{
    lx=0;
    ux=datah->GetNbinsX()+1;
  }
  for(Int_t i=lx;i<=ux;i++){
    Float_t derr=datah->GetBinError(i);
    Float_t mcerr=sum->GetBinError(i);
    if(datah->GetBinContent(i)>0 && !cumulative){
      chi2+=(datah->GetBinContent(i)-sum->GetBinContent(i))*(datah->GetBinContent(i)-sum->GetBinContent(i))/(derr*derr+mcerr*mcerr);
      ndf++;
    }
  }
  */
  //////////////////////////////////////////////////
  /// here's a quicker way to return the total chi2
  for(Int_t i=1;i<=datah->GetNbinsX();i++){
    if(datah->GetBinContent(i)>0) ndf++;
  }
  chi2=datah->Chi2Test(sum,"UWCHI2");
  //////////////////////////////////////////////////
  
  ostringstream oss3;
  oss3<<"Total MC = "<<fixed<<setprecision(0)<<mctot<<" #pm "<<Int_t(TMath::Sqrt(mcetot));
  leg->AddEntry(sum,oss3.str().c_str(),"l");
  
  ostringstream oss4;
  cout<<hname<<" ==> chi2/ndf = "<<fixed<<setprecision(1)<<chi2<<" / "<<setprecision(0)<<ndf<<endl;
  oss4<<"chi2/ndf = "<<fixed<<setprecision(1)<<chi2<<" / "<<setprecision(0)<<ndf;
  leg->AddEntry(sum,oss4.str().c_str(),"l");
  
  
  if(kLegend) leg->Draw();
  
}

//*******************************************************************
//
//   Wraper for Drawall_d1 function to be used with 
//   multi channels analysis
//
//*******************************************************************

void Drawall_d1_Ch(const char * hname, bana10 *data[],bana10 *bgr[][MAXCOMP], Int_t nbgr[], bana10 *sig[][MAXCOMP], Int_t nsig[],string namechnl[],Int_t nchnls,Bool_t kLegend, Bool_t kSum,Bool_t cumulative, Double_t blow,Double_t bup,Int_t rebin,Int_t fillstyle){

  //The idea: create a new bana10 instance for every MC and data component as a sum of representations in all channels
 
  bana10 * newdata;
  bana10 * newbgr[MAXCOMP];
  Int_t newnbgr = 0;
  bana10 * newsig[MAXCOMP];
  Int_t newnsig=0;


  Draw_Ch_merge(hname,data,bgr,nbgr,sig,nsig,namechnl,nchnls,newdata,newbgr,newnbgr,newsig,newnsig);


  Drawall_d1(hname, newdata,newbgr, newnbgr, newsig, newnsig,kLegend,kSum,cumulative,blow,bup,rebin,fillstyle);

}



//*******************************************************************
// This function Draws plot for hname histo
// It draws data, backgournds,sum of all backgrounds and legend
// Input parameters:
// hname -- name of the histogram to draw
// data,bgr,nbgr,sig,nsig -- pointers to data, background and signal components in the memory
// kLegend -- boolean flag to draw the legend or not (default = true)
// kSum -- boolean flag to draw one MC component on top of another. If false than MC components plots are superimposed on on another.(default=false)
// cumulative -- cumulative plot, each bin actually equal to sum of all higher bins in the original histogram. (default=false)
// blow -- low boundary for the plot. If 0 use histograms low boumdary (default=0)
// bup -- upper boundary for the plot. If 0 use histogram upper boundary (default=0)
// rebin -- rebin histogram (defoult = 1, no rebinning)
//*******************************************************************
void Drawconf_d1(const char *hname,ana10 *data,bana10 **bgr,Int_t nbgr,bana10 ** sig,Int_t nsig,Bool_t kLegend,Bool_t kSum,Bool_t cumulative,Double_t blow,Double_t bup,Int_t rebin,Int_t fillstyle,Bool_t draw_0n,bana10 **lim,Int_t nlim){
  
  if(cumulative){
    cout<<"Drawconf_d1 has no support for cumulative option."<<endl;
    return;
  }

  vector<TH1D *> highlight;
  vector<string> highlightname;

  //total background sum
  gROOT->SetStyle("Plain");
  TH1D * h; //dummy pointer
  TH1D *sumb=(TH1D*)bgr[0]->Geth1(hname)->Clone("bsum"); 
  //total signal sum
  TH1D *sums=(TH1D*)bgr[0]->Geth1(hname)->Clone("ssum"); 

  //group bgr into external, wire chamber and internal
  TH1D *sumex=(TH1D*)bgr[0]->Geth1(hname)->Clone("exsum"); 
  TH1D *sumrn=(TH1D*)bgr[0]->Geth1(hname)->Clone("rnsum"); 
  TH1D *sumin=(TH1D*)bgr[0]->Geth1(hname)->Clone("insum"); 

  //expected number of events
  Float_t mcetot=0;
  Float_t mctot=0;

  Float_t exetot=0;
  Float_t extot=0;
  Float_t rnetot=0;
  Float_t rntot=0;
  Float_t inetot=0;
  Float_t intot=0;

  Float_t sigetot=0;
  Float_t sigtot=0;


  TH1D *temp[1000]; 
  Float_t a[1000];
  Int_t index[1000];
  string legends[1000];

  //initialize all the histogramms
  sumb->Reset();
  sumb->Rebin(rebin);
  sums->Reset();
  sums->Rebin(rebin);
  sumex->Reset();
  sumex->Rebin(rebin);
  sumin->Reset();
  sumin->Rebin(rebin);
  sumrn->Reset();
  sumrn->Rebin(rebin);


  TLegend *leg= new TLegend(0.55,0.55,0.90,0.90);
  leg->SetFillColor(0);
  Int_t bcolor=3;

  TH1D * datah;
  Float_t di;
  Int_t lx,ux; //lower and upper limit in bins

  //deal with data histogram first
  if(data){
    datah=(TH1D *)data->Geth1(hname)->Clone("newdata");
    if(blow!=0){
      Double_t ledges[10000];
      datah->GetLowEdge(ledges);
      lx=0;
      ux=0;
      for(Int_t k=0;k<=datah->GetNbinsX()+1;k++){
	if(blow<ledges[k] && lx==0) lx = k-1;
	if(bup<ledges[k] && ux==0) ux = k-1;
      }
      if(ux==0) ux=datah->GetNbinsX()+1;
      
    }else{
      lx=0;
      ux=datah->GetNbinsX()+1;
    }
    di = datah->Integral(lx,ux);
    datah->Rebin(rebin);
    if(blow) datah->GetXaxis()->SetRangeUser(blow,bup);
  }
  datah->SetMinimum(0.);

  //loop through the backgrounds now
  for(Int_t i=0;i<nbgr;i++){
    temp[i]=(TH1D* ) bgr[i]->Geth1(hname)->Clone(string("temp"+i).c_str());
    temp[i]->Rebin(rebin);

    temp[i]->Scale(bgr[i]->norm);
    temp[i]->SetTitleSize(0.06,"X");
    temp[i]->SetTitleSize(0.06,"Y");
    temp[i]->SetLineWidth(4);

    Float_t ni;
    ni = bgr[i]->GetExpected(hname,lx,ux);
    Float_t dni=0;
    if(bgr[i]->norm!=0 and ni!=0) dni=bgr[i]->GetErrorExpected(hname,lx,ux);

    string bname = bgr[i]->Getname();
    Bool_t attached=false;

    mcetot+=dni*dni;
    mctot+=ni;
    sumb->Add(sumb,temp[i],1.,1.);
    if(!attached & (stest(bname,"ex") || stest(bname,"pm") || stest(bname,"sci") || stest(bname,"ss"))){
      //external background
      sumex->Add(sumex,temp[i],1.,1.);
      exetot+=dni*dni;
      extot+=ni;
      attached = true;
    }
    if(!attached & (stest(bname,"sw") || stest(bname,"sf") || stest(bname,"gas") || stest(bname,"wire"))){
      //tracking chamber background
      sumrn->Add(sumrn,temp[i],1.,1.);
      rnetot+=dni*dni;
      rntot+=ni;
      attached = true;
    }   
     if(!attached &(stest(bname,"foil") || stest(bname,"inbg") || stest(bname,"met") || stest(bname,"com") || stest(bname,"sen") || stest(bname,"seo") || stest(bname,"cu") || stest(bname,"myl") || stest(bname,"te") || stest(bname,"zr96") || stest(bname,"nd150") || stest(bname,"ca48")) ){
      //internal background
      sumin->Add(sumin,temp[i],1.,1.);
      inetot+=dni*dni;
      intot+=ni;
      attached = true;
    }
    if(!attached & (stest(bname,"2b2n"))){
      cout<<"Cloning 2b2n" <<endl;
      h = (TH1D*)bgr[0]->Geth1(hname)->Clone("2b2n"); 
      h->Reset();
      h->Rebin(rebin);
      h->Add(h,temp[i],1.,1.);
      highlight.push_back(h);
      string nn;
      if(stest(bname,"te")){
	nn = "^{130}Te 2#beta2#nu";
      }else if(stest(bname,"se")){
	nn = "^{82}Se 2#beta2#nu";
      }else{
	nn = "2#beta2#nu";
      }
      highlightname.push_back(nn);
      attached = true;
      cout<<"Cloned" <<endl;
    }
    if(!attached){
      cout<<"DRAWCONF WARNING: component "<<bname<<" was not recognised as internal, external or wire chamber"<<endl;
    }

  }
  Int_t j=nbgr;
  for(Int_t i=0;i<nsig;i++){
    temp[j]=(TH1D* ) sig[i]->Geth1(hname)->Clone(string("temp"+j).c_str());
    temp[j]->Rebin(rebin);

    temp[j]->Scale(sig[i]->norm);
    temp[j]->SetTitleSize(0.06,"X");
    temp[j]->SetTitleSize(0.06,"Y");
    temp[j]->SetLineWidth(4);

    Float_t ni;
    ni = sig[i]->GetExpected(hname,lx,ux);
    Float_t dni=0;
    if(ni!=0 && sig[i]->activity!=0) dni=sig[i]->GetErrorExpected(hname,lx,ux);

    sigetot+=dni*dni;
    sigtot+=ni;

    mcetot+=dni*dni;
    mctot+=ni;

    string bname = sig[i]->Getname();
    Bool_t attached=false;

    sums->Add(sums,temp[j],1.,1.);

    if(!attached & (stest(bname,"2b2n"))){
      cout<<"Cloning 2b2n" <<endl;
      h = (TH1D*)bgr[0]->Geth1(hname)->Clone("2b2n"); 
      h->Reset();
      h->Rebin(rebin);
      h->Add(h,temp[j],1.,1.);
      highlight.push_back(h);
      string nn;
      if(stest(bname,"te")){
	nn = "^{130}Te 2#beta2#nu";
      }else if(stest(bname,"se")){
	nn = "^{82}Se 2#beta2#nu";
      }else{
	nn = "2#beta2#nu";
      }
      highlightname.push_back(nn);
      attached = true;
      cout<<"Cloned" <<endl;
    }
    if(!attached & (stest(bname,"bi210"))){
      h = (TH1D*)bgr[0]->Geth1(hname)->Clone("b210"); 
      h->Reset();
      h->Rebin(rebin);
      h->Add(h,temp[j],1.,1.);
      highlight.push_back(h);
      string nn = "surface ^{210}Bi";
      highlightname.push_back(nn);
      attached = true;
    }
    if(!attached & (stest(bname,"234"))){
      cout<<"Cloning Pa234" <<endl;
      h = (TH1D*)bgr[0]->Geth1(hname)->Clone("p234"); 
      h->Reset();
      h->Rebin(rebin);
      h->Add(h,temp[j],1.,1.);
      highlight.push_back(h);
      string nn = "internal ^{234m}Pa";
      highlightname.push_back(nn);
      attached = true;
      cout<<"Cloned" <<endl;
    }
    if(!attached & (stest(bname,"k40"))){
      cout<<"Cloning K40" <<endl;
      h = (TH1D*)bgr[0]->Geth1(hname)->Clone("k40"); 
      h->Reset();
      h->Rebin(rebin);
      h->Add(h,temp[j],1.,1.);
      highlight.push_back(h);
      string nn = "internal ^{40}K";
      highlightname.push_back(nn);
      attached = true;
      cout<<"Cloned" <<endl;
    }
    if(!attached & (stest(bname,"ex") || stest(bname,"pm") || stest(bname,"sci") || stest(bname,"ss"))){
      //external background
      sumex->Add(sumex,temp[j],1.,1.);
      exetot+=dni*dni;
      extot+=ni;
      attached = true;
    }
    if(!attached & (stest(bname,"sw") || stest(bname,"sf") || stest(bname,"gas") || stest(bname,"wire"))){
      //tracking chamber background
      sumrn->Add(sumrn,temp[j],1.,1.);
      rnetot+=dni*dni;
      rntot+=ni;
      attached = true;
    }   
    if(!attached &(stest(bname,"foil") || stest(bname,"inbg") || stest(bname,"met") || stest(bname,"com") || stest(bname,"sen") || stest(bname,"seo") || stest(bname,"cu") || stest(bname,"myl") || stest(bname,"te") || stest(bname,"zr96") || stest(bname,"nd150") || stest(bname,"ca48"))){
      //internal background
      sumin->Add(sumin,temp[j],1.,1.);
      inetot+=dni*dni;
      intot+=ni;
      attached = true;
    }
    if(!attached){
      cout<<"DRAWCONF WARNING: component "<<bname<<" was not recognised as internal, external or wire chamber"<<endl;
    }
    j++;
  }

  sumrn->Add(sumex,1);
  sumin->Add(sumrn,1);

  if(kSum){
    sums->Add(sumb,1);
  }else{
    //datah->Add(sumb,-1.);
  }
  for(size_t i=0;i<highlight.size();i++){
    cout<<"Summing hl "<<i<<endl;
    if(i==0){
      highlight[i]->Add(sumin,1);
    }else{
      highlight[i]->Add(highlight[i-1]);
    }    
  }
  
  cout<<"Start drawing "<<endl;

  if(kSum){
    leg->AddEntry(datah,"raw data","lpe");
  }else{
    //    leg->AddEntry(datah,"data, bgr subtracted");
    leg->AddEntry(datah,"data","lpe");
  }
  leg->SetTextSize(0.04);

  datah->SetLineWidth(3);
  datah->SetMinimum(0.01);
  datah->SetMarkerStyle(20);
  datah->SetMarkerSize(1.);
  datah->SetStats(0);
  datah->SetTitle("");
  datah->GetYaxis()->SetTitleOffset(1.2);
  //  datah->SetXTitle("E_{#beta}, MeV");
  //  datah->SetYTitle("Counts");
  datah->SetMinimum(0.);
  datah->Draw("ep");
  
  for(size_t i=0;i<highlight.size();i++){
    cout<<"setting highlights"<<endl;
    highlight[i]->SetTitleSize(0.06,"X");
    highlight[i]->SetTitleSize(0.06,"Y");
    highlight[i]->SetLineWidth(2);
    highlight[i]->SetLineColor(1);
    highlight[i]->SetFillStyle(fillstyle);
    highlight[i]->SetFillColor(palitra(i+4));
    leg->AddEntry(highlight[i],highlightname[i].c_str(),"f");
    
  }

  /**
  sums->SetTitleSize(0.06,"X");
  sums->SetTitleSize(0.06,"Y");
  sums->SetLineWidth(2);
  sums->SetLineColor(1);
  sums->SetFillStyle(3344);
  sums->SetFillColor(8);
  //  leg->AddEntry(sums,"2#beta2#nu signal ","f");
  **/

  sumex->SetTitleSize(0.06,"X");
  sumex->SetTitleSize(0.06,"Y");
  sumex->SetLineWidth(2);
  sumex->SetFillStyle(fillstyle);
  sumex->SetFillColor(palitra(1));
  sumex->SetLineColor(1);
  sumex->SetLineStyle(2);
  leg->AddEntry(sumex,"external Bkg","f");

  sumrn->SetTitleSize(0.06,"X");
  sumrn->SetTitleSize(0.06,"Y");
  sumrn->SetLineWidth(2);
  sumrn->SetFillStyle(fillstyle);
  sumrn->SetFillColor(palitra(2));
  sumrn->SetLineColor(1);
  sumrn->SetLineStyle(2);
  leg->AddEntry(sumrn,"gg chamber Bkg","f");

  sumin->SetTitleSize(0.06,"X");
  sumin->SetTitleSize(0.06,"Y");
  sumin->SetLineWidth(2);
  sumin->SetFillStyle(fillstyle);
  sumin->SetFillColor(palitra(3));
  sumin->SetLineColor(1);
  sumin->SetLineStyle(2);
  leg->AddEntry(sumin,"internal Bkg","f");


  // sums->Draw("hist,same"); 
  for(Int_t i=highlight.size()-1;i>=0;i--){
    cout<<"Draw hl "<<i<<endl;
    highlight[i]->Draw("hist,same");
  }
 
  sumin->Draw("hist,same"); 
  sumrn->Draw("hist,same"); 
  sumex->Draw("hist,same"); 

  if(draw_0n && nlim > 0){
    TH1D *sum0n=(TH1D*)lim[0]->Geth1(hname)->Clone("0n_sum"); 

    Float_t A=1;
    Int_t isrc[10],nisrc;
    Float_t scale_t;
    /**
    //MO
    isrc[0]=0;
    isrc[1]=1;
    nisrc=2;
    scale_t = 3E+23;
    Float_t T = A2HalfLife(A,isrc,nisrc);
    A = T/scale_t;
    sum0n->Scale(A);

    **/
    sum0n->SetLineColor(14);
    sum0n->SetLineWidth(4);
    sum0n->Rebin(rebin);


    Float_t max_0n=0;
    Float_t max_d = 0;
    //estimate scaling
    for(int i=lx/rebin;i<=ux/rebin;i++){
      if(datah->GetBinContent(i)>max_d) max_d = datah->GetBinContent(i);
      if(sum0n->GetBinContent(i)>max_0n) max_0n= sum0n->GetBinContent(i);
    }
    scale_t = max_d/5/max_0n;
    cout<<"DRAWCONF: max_d max_0n scale "<<max_d<<" "<<max_0n<<" "<<scale_t<<endl;
    sum0n->Scale(scale_t);

    sum0n->Draw("hist,C,same");
    leg->AddEntry(sum0n,"2#beta0#nu spectrum","l");
}

  datah->Draw("epsame");


  leg->SetLineWidth(-1);
  if(kLegend) leg->Draw();
  //  TPaveText * pt = new TPaveText(0.55,0.65,0.98,0.90,"NDCR");

  /**
  pt->SetFillColor(0);
  pt->AddText("^{130}Te 454 g, 534 days");
  pt->AddText("607 events");
  pt->AddText("109 2#beta2#nu events");
  pt->AddText("S/B = 0.25");
  //  pt->SetTextAlign(3);
  pt->SetTextFont(12);
  //  pt->Draw();

  **/
}
//*******************************************************************
//
//   Wraper for Drawconf_d1 function to be used with 
//   multi channels analysis
//
//*******************************************************************

void Drawconf_d1_Ch(const char * hname, bana10 *data[],bana10 *bgr[][MAXCOMP], Int_t nbgr[], bana10 *sig[][MAXCOMP], Int_t nsig[],string namechnl[],Int_t nchnls,Bool_t kLegend, Bool_t kSum,Bool_t cumulative, Double_t blow,Double_t bup,Int_t rebin,Int_t fillstyle,Bool_t draw_0n,bana10 **lim,Int_t nlim){

  //The idea: create a new bana10 instance for every MC and data component as a sum of representations in all channels

  bana10 * newdata;
  bana10 * newbgr[MAXCOMP];
  Int_t newnbgr = 0;
  bana10 * newsig[MAXCOMP];
  Int_t newnsig=0;

  Draw_Ch_merge(hname, data,bgr,nbgr,sig,nsig,namechnl,nchnls,newdata,newbgr,newnbgr,newsig,newnsig);

  Drawconf_d1(hname, newdata,newbgr, newnbgr, newsig, newnsig,kLegend,kSum,cumulative,blow,bup,rebin,fillstyle,draw_0n,lim,nlim);

}


void Drawall_d2x(const char * hname, ana10 *data,bana10 **bgr, Int_t nbgr, bana10 ** sig, Int_t nsig,Bool_t kLegend,Int_t rebin){
  // Draws hname from all bgr histos at the same time;
  gROOT->SetStyle("Plain");
  TH1D *sum=(TH1D*)data->Geth2(hname)->ProjectionX((string("bsumt")+hname).c_str(),0,1000)->Clone((string("bsum")+hname).c_str()); 
  TH1D *temp[1000]; 
  Float_t a[1000];
  Int_t index[1000];
  string legends[1000];
  sum->Reset();
  //TLegend *leg= new TLegend(0.5,0.45,0.9,0.82);
  TLegend *leg= new TLegend(0.45,0.6,0.995,0.995);
  leg->SetFillColor(0);
  Int_t bcolor=3;
  if(data){
    string dtmp=string(hname)+string("x");
    TH1D * dtmpx=data->Geth2(hname)->ProjectionX(dtmp.c_str(),0,1000);
    dtmpx->SetMinimum(0);
    dtmpx->Rebin(rebin)->Draw("e");
    ostringstream oss2;
    oss2<<data->Getname()<<" ("<<dtmpx->Integral()<<") ";
  
    leg->AddEntry(dtmpx,oss2.str().c_str(),"l");
  }
  for(Int_t i=0;i<nbgr;i++){
    temp[i]=(TH1D* ) bgr[i]->Geth2(hname)->ProjectionX((string("temp")+hname+(char)i).c_str(),0,1000);
    temp[i]->Scale(bgr[i]->norm);
    temp[i]->SetTitleSize(0.06,"X");
    temp[i]->SetTitleSize(0.06,"Y");
    temp[i]->SetLineWidth(2);
    a[i]=temp[i]->Integral();
    //    temp[i]->Draw("same");
    ostringstream oss;
    oss<<bgr[i]->Getname()<<" ("<<temp[i]->Integral()<<") A="<<setprecision(4)<<bgr[i]->activity<<" ";
    //    if (temp[i]->Integral()>0) leg->AddEntry(temp[i],oss.str().c_str(),"l");
    legends[i]=oss.str();
    sum->Add(temp[i],1.);
    //cout<<"Draw x projection total sum = "<<sum->Integral()<<endl;
    //cout<<"Draw x temp "<<i<<" = "<<temp[i]->Integral()<<endl;
  }
  Int_t j=nbgr;
  for(Int_t i=0;i<nsig;i++){
    temp[j]=(TH1D*) sig[i]->Geth2(hname)->ProjectionX((string("temp")+hname+(char)j).c_str(),0,1000);
    temp[j]->Scale(sig[i]->norm);
    temp[j]->SetTitleSize(0.06,"X");
    temp[j]->SetTitleSize(0.06,"Y");
    temp[j]->SetLineWidth(2);
    a[j]=temp[j]->Integral();
    //    temp[j]->Draw("same");
    ostringstream oss;
    oss<<sig[i]->Getname()<<" ("<<temp[j]->Integral()<<") A="<<setprecision(4)<<sig[i]->activity<<" ";
    //    if (temp[j]->Integral()>0) leg->AddEntry(temp[j],oss.str().c_str(),"l");
    legends[j]=oss.str();
    sum->Add(temp[j],1.);
    //cout<<"Draw x projection total sum = "<<sum->Integral()<<endl;
    //cout<<"Draw x temp "<<j<<" = "<<temp[j]->Integral()<<endl;
    j++;
  }
  
  TMath::Sort(nbgr+nsig,a,index);
  for(Int_t i=0;i<nbgr+nsig;i++){
    temp[index[i]]->SetLineColor(palitra(i+bcolor));
    temp[index[i]]->Rebin(rebin)->Draw("same");
    if(a[index[i]]>0 && i<10){
      leg->AddEntry(temp[index[i]],legends[index[i]].c_str(),"l");
    }
  }

    sum->SetTitleSize(0.06,"X");
    sum->SetTitleSize(0.06,"Y");
    sum->SetLineWidth(2);
    sum->SetLineColor(palitra(2));
    sum->Rebin(rebin)->Draw("hist,same"); 
    ostringstream oss1;
    oss1<<"Total MC ("<<sum->Integral()<<") ";
    leg->AddEntry(sum,oss1.str().c_str(),"l");
    if(kLegend) leg->Draw();
}

void Drawall_d2y(const char * hname, ana10 *data,bana10 **bgr, Int_t nbgr, bana10 ** sig, Int_t nsig,Bool_t kLegend,Int_t rebin){
  // Draws hname from all bgr histos at the same time;
  gROOT->SetStyle("Plain");
  TH1D *sum=(TH1D*)bgr[0]->Geth2(hname)->ProjectionY((string("bsumty")+hname).c_str(),0,1000)->Clone((string("bsumy")+hname).c_str()); 
  TH1D *temp[1000]; 
  TH1D * dtmpx;
  Float_t a[1000];
  Int_t index[1000];
  string legends[1000];
  sum->Reset();
  //TLegend *leg= new TLegend(0.5,0.45,0.9,0.82);
  TLegend *leg= new TLegend(0.45,0.6,0.995,0.995);
  leg->SetFillColor(0);
  Int_t bcolor=3;
  if(data){
    string dtmp=string(hname)+string("y");
    dtmpx=data->Geth2(hname)->ProjectionY(dtmp.c_str(),0,1000);
    dtmpx->SetMinimum(0);
    dtmpx->Rebin(rebin)->Draw("e");
    ostringstream oss2;
    oss2<<data->Getname()<<" ("<<dtmpx->Integral()<<") ";
  
    leg->AddEntry(dtmpx,oss2.str().c_str(),"l");
  }
  for(Int_t i=0;i<nbgr;i++){
    temp[i]=(TH1D* ) bgr[i]->Geth2(hname)->ProjectionY((string("ytemp")+hname+(char)i).c_str(),0,1000);
    temp[i]->Scale(bgr[i]->norm);
    temp[i]->SetTitleSize(0.06,"X");
    temp[i]->SetTitleSize(0.06,"Y");
    temp[i]->SetLineWidth(2);
    a[i]=temp[i]->Integral();
    //    temp[i]->Draw("same");
    ostringstream oss;
    oss<<bgr[i]->Getname()<<" ("<<temp[i]->Integral()<<") A="<<setprecision(4)<<bgr[i]->activity<<" ";
    //    if (temp[i]->Integral()>0) leg->AddEntry(temp[i],oss.str().c_str(),"l");
    legends[i]=oss.str();
    sum->Add(sum,temp[i],1.,1.);
  }
  Int_t j=nbgr;
  for(Int_t i=0;i<nsig;i++){
    temp[j]=(TH1D*) sig[i]->Geth2(hname)->ProjectionY((string("ytemp")+hname+(char)j).c_str(),0,1000);
    temp[j]->Scale(sig[i]->norm);
    temp[j]->SetTitleSize(0.06,"X");
    temp[j]->SetTitleSize(0.06,"Y");
    temp[j]->SetLineWidth(2);
    a[j]=temp[j]->Integral();
    //    temp[j]->Draw("same");
    ostringstream oss;
    oss<<sig[i]->Getname()<<" ("<<temp[j]->Integral()<<") A="<<setprecision(4)<<sig[i]->activity<<" ";
    //    if (temp[j]->Integral()>0) leg->AddEntry(temp[j],oss.str().c_str(),"l");
    legends[j]=oss.str();
    sum->Add(sum,temp[j],1.,1.);
    j++;
  }
  
  TMath::Sort(nbgr+nsig,a,index);
  for(Int_t i=0;i<nbgr+nsig;i++){
    temp[index[i]]->SetLineColor(palitra(i+bcolor));
    temp[index[i]]->Rebin(rebin)->Draw("same");
    if(a[index[i]]>0 && i<10){
      leg->AddEntry(temp[index[i]],legends[index[i]].c_str(),"l");
    }
  }

    sum->SetTitleSize(0.06,"X");
    sum->SetTitleSize(0.06,"Y");
    sum->SetLineWidth(2);
    sum->SetLineColor(palitra(2));
    sum->Rebin(rebin)->Draw("hist,same"); 
    ostringstream oss1;
    oss1<<"Total MC ("<<sum->Integral()<<") ";
    leg->AddEntry(sum,oss1.str().c_str(),"l");
    if(kLegend) leg->Draw();
}

void Drawdiff_d1(const char * hname, ana10 *data,bana10 **bgr, Int_t nbgr, bana10 ** sig, Int_t nsig,Bool_t kLegend, int rebin){
  // Draws hname from all bgr histos at the same time;
  gROOT->SetStyle("Plain");
  TH1D *sum=(TH1D*)bgr[0]->Geth1(hname)->Clone("bsum"); 
  TH1D *temp[1000]; 
  sum->Reset();
  Int_t bcolor=2;

  for(Int_t i=0;i<nbgr;i++){
    temp[i]=(TH1D* ) bgr[i]->Geth1(hname)->Clone(string("temp"+i).c_str());
    temp[i]->Add(temp[i],bgr[i]->Geth1(hname),0,bgr[i]->norm);
    sum->Add(sum,temp[i],1.,1.);
  }
  Int_t j=nbgr;
  for(Int_t i=0;i<nsig;i++){
    temp[j]=(TH1D* ) sig[i]->Geth1(hname)->Clone(string("temp"+j).c_str());
    temp[j]->Add(temp[j],sig[i]->Geth1(hname),0,sig[i]->norm);
    sum->Add(sum,temp[j],1.,1.);
    j++;
  }
  
    TH1D *diff=(TH1D*)data->Geth1(hname)->Clone("dsum"); 
    diff->Add(sum,-1.);

    diff->Rebin(rebin);
    sum->Rebin(rebin);
    diff->Divide(sum);

    diff->SetTitleSize(0.06,"X");
    diff->SetTitleSize(0.06,"Y");
    diff->SetLineWidth(2);
    diff->SetLineColor(nbgr+nsig+bcolor);
    diff->Draw(); 
}

string printenumber(Float_t number, Float_t error, Int_t precision){
  ostringstream oss;
  //find main base first
  Int_t base = Int_t(TMath::Log10(number));
  if(number==0) base=0;
  if(base<0) base=base-1;
  if(error!=0){
    if(base>=0 && base<2){
      Int_t errp=precision-(base-Int_t(TMath::Log10(error)));
      if(errp<=0){
	precision=precision-errp+1;
	errp=1;
      }
      oss<<fixed<<setprecision(precision)<<number<<" #pm "<<setprecision(errp)<<error;
    }else{
      Int_t errp=precision-(base-Int_t(TMath::Log10(error)));
      if(errp<=0){
	precision=precision-errp+1;
	errp=1;
      }
      oss<<"("<<fixed<<setprecision(precision)<<number/TMath::Power(10.,base)<<" #pm "<<setprecision(errp)<<error/TMath::Power(10.,base)<<") #times10^{"<<base<<"}";
    }
  }else{
    if(base>=0 && base<3){
      oss<<fixed<<setprecision(precision)<<number;
    }else{
      oss<<fixed<<setprecision(precision)<<number/pow(10.,base)<<" #times10^{"<<base<<"}";
    }
  }
  return oss.str();
}


void Drawsn_d1(const char * hname, ana10 *data,bana10 **bgr, Int_t nbgr, bana10 ** sig, Int_t nsig,Bool_t kLegend, Bool_t kSum,Bool_t cumulative, Double_t blow,Double_t bup,Int_t rebin){

  gROOT->SetStyle("Plain");
  TH1D *sum=(TH1D*)bgr[0]->Geth1(hname)->Clone("bsum"); 
  Float_t mcetot=0;
  Float_t mctot=0;
  TH1D *temp[1000]; 
  Float_t a[1000];
  Int_t index[1000];
  string legends[1000];
  sum->Reset();
  sum->Rebin(rebin);

  //TLegend *leg= new TLegend(0.45,0.45,0.98,0.82);
  TLegend *leg= new TLegend(0.45,0.7,0.995,0.995);
  leg->SetFillColor(0);
  Int_t bcolor=3;
  TH1D * datah;
  Float_t di;
  Int_t lx,ux; //lower and upper limit in bins

  if(data){
    if(cumulative){
      datah=(TH1D *)GetCumulativeDown(data->Geth1(hname)->Rebin(rebin,"newdata"),bup);
    }else{
      datah=(TH1D *)data->Geth1(hname)->Clone("newdata");
    }
    if(blow!=0){
      Double_t ledges[10000];
      datah->GetLowEdge(ledges);
      lx=0;
      ux=0;
      for(Int_t k=0;k<datah->GetNbinsX()+1;k++){
	if(blow<ledges[k] && lx==0) lx = k-1;
	if(bup<ledges[k] && ux==0) ux = k-1;
      }
      if(ux==0) ux=datah->GetNbinsX()+1;
      
    }else{
      lx=0;
      ux=datah->GetNbinsX()+1;
    }
    if(cumulative){
      di = datah->GetBinContent(lx);
    }else{
      di = datah->Integral(lx,ux);
      datah->Rebin(rebin);
    }
    if(blow) datah->GetXaxis()->SetRangeUser(blow,bup);
    datah->SetMinimum(0.01);
    datah->SetMarkerStyle(20);
    
    ostringstream oss2;
    oss2<<data->Getname()<<" ("<<di<<") ";
  
  }
  for(Int_t i=0;i<nbgr;i++){
    //      cout<<"Draw "<<bgr[i]->Getname()<<" A="<<bgr[i]->activity<<" ea="<<bgr[i]->eactivity<<endl; 
    if(cumulative){
      temp[i]=(TH1D *)GetCumulativeDown(bgr[i]->Geth1(hname)->Rebin(rebin,"newtemp"),bup);
    }else{
      temp[i]=(TH1D* ) bgr[i]->Geth1(hname)->Clone(string("temp"+i).c_str());
      temp[i]->Rebin(rebin);
    }
    //    temp[i]->Add(temp[i],bgr[i]->Geth1(hname),0,bgr[i]->norm);
    temp[i]->Scale(bgr[i]->norm);
    temp[i]->SetTitleSize(0.045,"X");
    temp[i]->SetTitleSize(0.045,"Y");
    temp[i]->SetLabelSize(0.045,"X");
    temp[i]->SetLabelSize(0.045,"Y");
    temp[i]->SetLineWidth(2);
    ostringstream oss;
    Float_t ni;
    ni = bgr[i]->GetExpected(hname,lx,ux);
    Float_t dni=0;
    if(bgr[i]->norm!=0 and ni!=0) dni=bgr[i]->GetErrorExpected(hname,lx,ux);
    mcetot+=dni*dni;
    mctot+=ni;
    if(stest(bgr[i]->Getname(),"2b2n")){
      oss<<"^{82}Se #beta#beta2#nu";
    }
    if(stest(bgr[i]->Getname(),"tl208")){
      oss<<"^{208}Tl";
    }
    if(stest(bgr[i]->Getname(),"bi214")){
      oss<<"^{214}Bi";
    }

    a[i]=ni;
    legends[i]=oss.str();
    sum->Add(sum,temp[i],1.,1.);
  }
  Int_t j=nbgr;
  for(Int_t i=0;i<nsig;i++){
    if(cumulative){
      temp[j]=(TH1D *)GetCumulativeDown(sig[i]->Geth1(hname)->Rebin(rebin,"newtemp"),bup);      
    }else{
      temp[j]=(TH1D* ) sig[i]->Geth1(hname)->Clone(string("temp"+j).c_str());
      temp[j]->Rebin(rebin);
    }
    //    temp[j]->Add(temp[j],sig[i]->Geth1(hname),0,sig[i]->norm);
    temp[j]->SetLabelSize(0.05,"X");
    temp[j]->SetLabelSize(0.05,"Y");
    temp[j]->Scale(sig[i]->norm);
    temp[j]->SetTitleSize(0.02,"X");
    temp[j]->SetTitleSize(0.02,"Y");
    temp[j]->SetLineWidth(3);
    ostringstream oss;
    Float_t ni;
    ni = sig[i]->GetExpected(hname,lx,ux);
    Float_t dni=0;
    if(ni!=0 && sig[i]->activity!=0) dni=sig[i]->GetErrorExpected(hname,lx,ux);
    mcetot+=dni*dni;
    mctot+=ni;
    if(sig[i]->Getname() =="2b2n"){
      oss<<"^{82}Se #beta#beta2#nu";
    }
    a[j]=ni;
    legends[j]=oss.str();
    sum->Add(sum,temp[j],1.,1.);
    j++;
  }
  TMath::Sort(nbgr+nsig,a,index);
  if(!kSum){
    for(Int_t i=0;i<nbgr+nsig;i++){
      temp[index[i]]->SetLineColor(palitra(i+bcolor));
      temp[index[i]]->Draw("same");
      if(i<10 && a[index[i]]>0){
	leg->AddEntry(temp[index[i]],legends[index[i]].c_str(),"l");
    }
    }
  }
  if(kSum){
    for(Int_t i=nbgr+nsig-2;i>=0;i--){
      temp[index[i]]->Add(temp[index[i+1]]);
    }
    for(Int_t i=0;i<nbgr+nsig;i++){
      temp[index[i]]->SetLineColor(palitra(i+bcolor));
      temp[index[i]]->SetFillColor(palitra(i+bcolor));
      temp[index[i]]->SetFillStyle(1001); // 1001
      if(blow) temp[index[i]]->GetXaxis()->SetRangeUser(blow,bup);
      temp[index[i]]->Draw("hist,same");
      if(i<10 && a[index[i]]>0){
	leg->AddEntry(temp[index[i]],legends[index[i]].c_str(),"f");
      }
    }
  }
  sum->Draw("hist,same"); 
  
  ostringstream oss1;
  oss1 <<"Total MC ";

  leg->AddEntry(sum,oss1.str().c_str(),"l");
  if(kLegend) leg->Draw();
}


//*****************************************************************
//
// Utility function to merge several channels into one to 
// draw them all together by other Draw functions
//
//*****************************************************************

void Draw_Ch_merge(const char * hname, bana10 *data[],bana10 *bgr[][MAXCOMP], Int_t nbgr[], bana10 *sig[][MAXCOMP], Int_t nsig[],string namechnl[],Int_t nchnls,bana10* & newdata,bana10* newbgr[MAXCOMP],Int_t& newnbgr, bana10 * newsig[MAXCOMP], Int_t& newnsig){

  //The idea: create a new bana10 instance for every MC and data component as a sum of representations in all channels


  newnbgr=0;
  newnsig=0;
  map<string,bool> newamap;
  map<string,bool>::iterator itr;
 
  for(size_t i = 0; i<nchnls; i++){
    //deal with the backgrounds first
    string chname = namechnl[i];
    for(size_t j = 0; j<nbgr[i];j++){
      string bname = bgr[i][j]->Getname();
      bname.erase(0,chname.length());
      if ( newamap.find(bname) == newamap.end()){ 
	//	cout << "Coping component " <<bname << endl;
	//create a new bana10 instance
	bana10 * nb = new bana10(bname.c_str(),1,1);
	//copy histogram definition
	TH1D* bh = bgr[i][j]->Geth1(hname);
	nb->Addh1(hname,bh);
	bh = nb->Geth1(hname);
	bh->Reset();
	Int_t nbins = bh->GetNbinsX();
	Double_t h[10000];
	Double_t eh[10000];		
	double totactivity = 0;
	//loop over all channels again
	for(size_t k = 0; k<nchnls;k++){
	  for(size_t l = 0;l<nbgr[k];l++){
	    if(string(bgr[k][l]->Getname()).find(bname) != string::npos){
	      TH1D* htemp = bgr[k][l]->Geth1(hname);
	      if(htemp->GetNbinsX() != nbins){
		cout<<"Drawall_d1Ch ERROR: histogram size is different in different channels, cant draw!"<<endl;
		return;
	      }
	      //	      cout << "Adding component " <<bgr[k][l]->Getname()<<" A="<<bgr[k][l]->activity<<" htemp = "<<htemp->Integral()<< endl;
	      for(size_t n = 0;n<nbins+2;n++){
		if(totactivity==0){
		  h[n]=0;
		  eh[n]=0;
		}
		h[n]+=htemp->GetBinContent(n)*bgr[k][l]->activity;
		eh[n]+=pow(htemp->GetBinError(n)*bgr[k][l]->activity,2)+pow(htemp->GetBinContent(n)*bgr[k][l]->eactivity,2);		
		//		cout<<"h ("<<n<<")="<<h[n]<<endl;
	      }
	      totactivity +=bgr[k][l]->activity;
	      
	    }
	  }
	}//end of looping over all channels
	//Fill the histogram
	if(totactivity == 0) totactivity = 1;
	totactivity /=nchnls;
	for(size_t n = 0;n<nbins+2;n++){
	  bh->SetBinContent(n,h[n]/totactivity);
	  bh->SetBinError(n,sqrt(eh[n])/totactivity);
	}
	nb->setnormalization(totactivity,0);
	//store the background in the list
	newamap[bname] = true;
	newbgr[newnbgr] = nb;
	newnbgr++;
	// 	cout<<"bgr "<<bname<<" bana10 "<<nb<<" nbgr "<<newnbgr<<" bh="<<bh->Integral()<<" totact="<<totactivity<<endl;
     }
      
    }//end of dealing with backgrounds
    //signals now
    for(size_t j = 0; j<nsig[i];j++){
      string bname = sig[i][j]->Getname();
      bname.erase(0,chname.length());
      if ( newamap.find(bname) == newamap.end()){ 
	cout << "Coping component" <<bname << endl;
	double totalactivity = 0;
	//create a new bana10 instance
	bana10 * nb = new bana10(bname.c_str(),1,1);
	//copy histogram definition
	TH1D* bh = sig[i][j]->Geth1(hname);
	nb->Addh1(hname,bh);
	bh = nb->Geth1(hname);
	bh->Reset();
	Int_t nbins = bh->GetNbinsX();
	Double_t h[10000];
	Double_t eh[10000];		
	//loop over all channels again
	for(size_t k = 0; k<nchnls;k++){
	  for(size_t l = 0; l<nsig[k];l++){
	    if(string(sig[k][l]->Getname()).find(bname) != string::npos){
	      TH1D* htemp = sig[k][l]->Geth1(hname);
	      if(htemp->GetNbinsX() != nbins){
		cout<<"Drawall_d1Ch ERROR: "<<hname<<" size is different in channels for component "<<bname<<", cant draw!"<<endl;
		return;
	      }
	      for(size_t n = 0;n<nbins+2;n++){
		if(totalactivity==0){
		  h[n]=0;
		  eh[n]=0;
		}
		h[n]+=htemp->GetBinContent(n)*sig[k][l]->activity;
		eh[n]+=pow(htemp->GetBinError(n)*sig[k][l]->activity,2)+pow(htemp->GetBinContent(n)*sig[k][l]->eactivity,2);
	      }
	      totalactivity+=sig[k][l]->activity;

	    }
	  }
	}//end of looping over all channels
	//Fill the histogram
	if(totalactivity==0) totalactivity = 1;
	totalactivity /= nchnls;
	for(size_t n = 0;n<nbins+2;n++){
	  bh->SetBinContent(n,h[n]/totalactivity);
	  bh->SetBinError(n,sqrt(eh[n])/totalactivity);
	}
	//store the background in the list
	nb->setnormalization(totalactivity,0);
	newamap[bname] = true;
	newsig[newnsig] = nb;
	newnsig++;
      }
      
    }//end of dealing with signals


  }//end of nchnls loop

  //Finally the data
  string bname = data[0]->Getname();
  bname.erase(0,namechnl[0].length());
  if ( newamap.find(bname) == newamap.end()){ 
    cout << "Coping component" <<bname << endl;
    //create a new bana10 instance
    bana10 * nb = new bana10(bname.c_str(),1,1);
    //copy histogram definition
    TH1D* bh = data[0]->Geth1(hname);
    nb->Addh1(hname,bh);
    bh = nb->Geth1(hname);
    bh->Reset();
    Int_t nbins = bh->GetNbinsX();
    cout<<"datat nbins = "<<nbins<<endl;    
    Double_t h[10000];
    Double_t eh[10000];		
    //loop over all channels again
    for(size_t k = 0;k<nchnls;k++){
	if(true){
	  TH1D* htemp = data[k]->Geth1(hname);
	  if(htemp->GetNbinsX() != nbins){
	    cout<<"Drawall_d1Ch ERROR: "<<hname<<" size is different in channels for component "<<bname<<", cant draw!"<<endl;
	    return;
	  }
	  for(size_t n = 0;n<nbins+2;n++){
	    if(k==0){
	      h[n]=0;
	      eh[n]=0;
	    }
	    h[n]+=htemp->GetBinContent(n);
	    eh[n]+=pow(htemp->GetBinError(n),2);
	  }
	  
	}
    }//end of looping over all channels
    //Fill the histogram
    for(size_t n = 0;n<nbins+2;n++){
      bh->SetBinContent(n,h[n]);
      bh->SetBinError(n,sqrt(eh[n]));
    }
    //store the background in the list
    newamap[bname] = true;
    newdata = nb;
  }
      


}

