#include "ana.hpp"

#include <string>
#include <sstream>
#include <iomanip>
#include "../../utils/h10.h"
#include "TH2.h"
#include "TH1.h"
#include "TLegend.h"
#include "TMath.h"



//************************************************************
// Calculate Likelihood ratio (test statistics) for LLH analysis
// X = Psb(data)/Pb(data)
// Psb -- probability to observe data at s+b hypothesis
// Pb -- probability to observe data at b only hypothesis
// both PDF shoud be extebdable!
//************************************************************
Double_t X(RooAbsPdf &bspdf, RooAbsPdf &bpdf,RooDataSet &data){
  Double_t Xstat=0.1;
    RooNLLVar *nllsb=new RooNLLVar("nllsb","nllsb",bspdf,data,RooCmdArg("Extended",kTRUE));
    RooNLLVar *nllb=new RooNLLVar("nllb","nllb",bpdf,data,RooCmdArg("Extended",kTRUE));

    Xstat = -nllsb->getVal(data.get())+nllb->getVal(data.get());//NLL returns -Log(L)!!!
    delete nllsb;
    delete nllb;
    return Xstat;
}

Double_t X1(RooAbsPdf &sbpdf, RooAbsPdf &bpdf,RooDataSet &data){
  Double_t Xstat=0.0;

    RooFormulaVar llratio("llratio","log(sbmodel/bmodel)","log(@0/@1)",RooArgList(sbpdf,bpdf));
    data.addColumn(sbpdf);
    data.addColumn(bpdf);
    data.addColumn(llratio);
  
  for(Int_t i=0;i<data.numEntries();i++){
          data.get(i);
       Xstat+=llratio.getVal();
  }

    Xstat-=sbpdf.extendedTerm((Int_t)data.sumEntries(),data.get());
    Xstat+=bpdf.extendedTerm((Int_t)data.sumEntries(),data.get());
  return Xstat ;
}
//************************************************************
//
//Calculates CL for given s+b and hyp in comparison with observed data 
//Use MC method
//
//************************************************************
Double_t LLHCL(RooAbsPdf *bpdf,RooAbsPdf *sbpdf,RooDataSet *expdata,const Double_t precision){
  
  Bool_t trace = false;  

  if(trace){
    cout<<"LLHCL PDF's. s+b pdf, norm = "<<sbpdf->expectedEvents(0)<<endl;
    sbpdf->Print("v");
    cout<<"LLHCL PDF's. b pdf, norm = "<<bpdf->expectedEvents(0)<<endl;
    bpdf->Print("v");
    expdata->Print("v");
  }


  RooDataSet *_expdata = (RooDataSet *)expdata->Clone("expcopy");
  Double_t Xobs = X(*sbpdf,*bpdf,*_expdata);
  delete _expdata;

  Long_t NCLSB=0;
  Long_t NCLB=0;

  Long_t Ntoy=1.65/precision/precision/0.9;// with 99% probability MC CL value will lay within precision

  if(trace) cout<<"Xobs="<<Xobs<<endl;

  RooArgSet *obsvar = sbpdf->getDependents(expdata);
  if(Ntoy>1000) Ntoy=1000;
  Ntoy=300;
  cout<<"LLHCL "<<Ntoy<<" MC to generate"<<endl;

  RooDataSet *toyb,*toysb;
  for(Long_t i=0;i<Ntoy;i++){
    //generate toy MC
    
    Double_t Xb=0,Xsb=0;
    toyb=bpdf->generate(*obsvar,0,RooCmdArg("Extended",kTRUE));
    toysb=sbpdf->generate(*obsvar,0,RooCmdArg("Extended",kTRUE));

    Xb = X(*sbpdf,*bpdf,*toyb);
    if(Xb<Xobs) NCLB++;
    Xsb = X(*sbpdf,*bpdf,*toysb);
    if(Xsb<Xobs) NCLSB++;
    if(i%10==0){
      cout<<i;
      if(trace){
	cout<<" Xsb="<<Xsb<<" Xb="<<Xb<<endl;
	toysb->Print();
      }else{
	cout<<endl;
      }
    }
    delete toyb;
    delete toysb;
  }
  cout<<endl;
  if(trace){
    cout<<" NCLSB="<<NCLSB<<" NCLB="<<NCLB<<endl;;
  }
  Double_t LL = ((Double_t)NCLSB)*1.0/NCLB;
  delete obsvar;
  return LL;
}

//***********************************************************
//
//Calculate signal limit at given CL with LLH method
//
//***********************************************************
Double_t LLHlimit(TList *bpdfs, TList *bnorms,RooAbsPdf *spdf,RooDataSet *expdata,const Double_t CL, const Double_t eps,const Double_t guesslimit){

  Bool_t trace = true;
  RooRealVar signorm("signorm","signorm",0,0,1000000);
  RooArgList argbpdfs(*bpdfs,"Backgrounds");
  RooArgList  argnbpdfs(*bnorms,"N Backgrounds");
  RooAbsPdf *bpdf = new RooAddPdf("bmodel","bmodel",argbpdfs,argnbpdfs);

  argbpdfs.add(*spdf);
  argnbpdfs.add(signorm);
  RooAbsPdf *sbpdf;
  Double_t a=0.01;
  Double_t b=guesslimit;
  Double_t c,CLa,CLb,CLc;
  Double_t precision=0.1;
  signorm = a;
  sbpdf= new RooAddPdf("sbmodel","sbmodel",argbpdfs,argnbpdfs);  
  CLa=1;
  delete sbpdf;
  if(trace) cout<<" a="<<a<<" CLa="<<CLa;

  signorm = b;
  sbpdf= new RooAddPdf("sbmodel","sbmodel",argbpdfs,argnbpdfs);  
  CLb=0;
  delete sbpdf;
  if(trace) cout<<" b="<<b<<" CLb="<<CLb;
  

  c=a+(b-a)/2.;
  signorm = c;
  sbpdf= new RooAddPdf("sbmodel","sbmodel",argbpdfs,argnbpdfs);  


  CLc=LLHCL(bpdf,sbpdf,expdata,precision);
  delete sbpdf;
  if(trace) cout<<"c="<<c<<" CLc="<<CLc;
  if(CLc==CL) return c;

  while(fabs(a-b)>2*eps){
    if(CLc<CL){
      b=c;
      CLb=CLc;
      precision=0.1;
      if(CLb>CL-0.1) precision=0.1;
      if(CLb>CL-0.05) precision = 0.05;
      if(CLb>CL-0.02) precision = 0.02;
    }else{
      a=c;
      CLa=CLc;
      precision=0.1;
      if(CLa<CL+0.1) precision=0.1;
      if(CLa<CL+0.05) precision = 0.05;
      if(CLa<CL+0.02) precision = 0.02;
    }
    c=a+(b-a)*(CLa-CL)/(CLa-CLb);
    signorm = c;
    sbpdf= new RooAddPdf("sbmodel","sbmodel",argbpdfs,argnbpdfs);  
    CLc=LLHCL(bpdf,sbpdf,expdata,precision);
    delete sbpdf;
    if(trace) cout<<"c="<<c<<" CLc="<<CLc;
 }
     if(CLc<CL){
      b=c;
      CLb=CLc;
    }else{
      a=c;
      CLa=CLc;
    }
    c=(a+b)/2;
 return c;
}

Double_t RooLLHLimit(const char * name, bana10 *data,bana10 **bgr,Int_t nbgr,bana10 ** sig,Int_t nsig,bana10 ** lim, Int_t nlim,Double_t xmin){

  //xmin -- lower boundary for pdf from wich to start the fit
  if(nlim==0){
    cout<<"RooLLHLimit: ERROR, no components to set limit at specified"<<endl;
    return 10.;
  }   

  if(!data->Checkh1(name)){
    cout<<"RooLLHLimit: ERROR, histogram "<<name<<" was not booked in uloop.C"<<endl;
    return 10.;
  }


    RooRealVar mes("Esum","Electrons energy sum, MeV",3.,xmin,4.001); 
    RooRealVar mes2("Ediff","|E1-E2|, MeV",0.,3.8); 
    TList bgrpdfs,bgrnorms;


    Int_t nbins=data->Geth1(name)->GetNbinsX();
    Double_t lbound = data->Geth1(name)->GetBinLowEdge(1);
    Double_t binsize=data->Geth1(name)->GetBinWidth(nbins-1);

    for(Int_t i=0;i<nbgr;i++){
      bgr[i]->ConstructPdf(mes,true);
      bgrpdfs.Add(bgr[i]->d1pdf);
      bgrnorms.Add(new RooRealVar(("N"+bgr[i]->name).c_str(),("N"+bgr[i]->name).c_str(),bgr[i]->GetExpected(name,1+(Int_t)(xmin/binsize-lbound/binsize),200),0,10000000));
    }

    lim[0]->ConstructPdf(mes,true);
  
    data->ConstructSet(mes,true);

    //Calculate LH ratio ststistics
  
    RooDataSet *expdata = (RooDataSet* )data->Esum_set;
    Double_t limitLLH=10.;
    limitLLH= LLHlimit( &bgrpdfs, &bgrnorms,lim[0]->d1pdf,expdata,0.1,0.3,30.);
    return limitLLH;
}
