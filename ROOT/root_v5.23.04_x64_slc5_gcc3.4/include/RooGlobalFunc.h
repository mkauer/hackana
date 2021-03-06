/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 *    File: $Id: RooGlobalFunc.h,v 1.14 2007/07/16 21:04:28 wouter Exp $
 * Authors:                                                                  *
 *   WV, Wouter Verkerke, UC Santa Barbara, verkerke@slac.stanford.edu       *
 *   DK, David Kirkby,    UC Irvine,         dkirkby@uci.edu                 *
 *                                                                           *
 * Copyright (c) 2000-2005, Regents of the University of California          *
 *                          and Stanford University. All rights reserved.    *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/
#ifndef ROO_GLOBAL_FUNC
#define ROO_GLOBAL_FUNC

#include "RooCmdArg.h"
#include "RooDataHist.h"
#include "RooAbsPdf.h"
#include "RooRealConstant.h"
#include "RooMsgService.h"

class RooAbsData ;
class RooArgSet ;
class RooCategory ;
class RooAbsReal ;
class RooAbsBinning ;
class RooAbsPdf ;
class RooConstVar ;
class RooRealVar ;
class RooAbsCategory ;
class TF1 ;
class TF2 ;
class TF3 ;

namespace RooFit {

// RooAbsReal::plotOn arguments
RooCmdArg DrawOption(const char* opt) ;
RooCmdArg Normalization(Double_t scaleFactor) ;
RooCmdArg Slice(const RooArgSet& sliceSet) ;
RooCmdArg Slice(RooCategory& cat, const char* label) ;
RooCmdArg Project(const RooArgSet& projSet) ;
RooCmdArg ProjWData(const RooAbsData& projData, Bool_t binData=kFALSE) ;
RooCmdArg ProjWData(const RooArgSet& projSet, const RooAbsData& projData, Bool_t binData=kFALSE) ;
RooCmdArg Asymmetry(const RooCategory& cat) ;
RooCmdArg Precision(Double_t prec) ;
RooCmdArg ShiftToZero() ;
RooCmdArg Range(const char* rangeName, Bool_t adjustNorm=kTRUE) ;
RooCmdArg Range(Double_t lo, Double_t hi, Bool_t adjustNorm=kTRUE) ;
RooCmdArg NormRange(const char* rangeNameList) ;
RooCmdArg VLines() ;
RooCmdArg LineColor(Color_t color) ;
RooCmdArg LineStyle(Style_t style) ;
RooCmdArg LineWidth(Width_t width) ;
RooCmdArg FillColor(Color_t color) ;
RooCmdArg FillStyle(Style_t style) ;
RooCmdArg ProjectionRange(const char* rangeName) ;
RooCmdArg Name(const char* name) ;
RooCmdArg Invisible() ;
RooCmdArg AddTo(const char* name, double wgtSel=1.0, double wgtOther=1.0) ;
RooCmdArg EvalErrorValue(Double_t value) ;
RooCmdArg MoveToBack()  ;

// RooAbsPdf::plotOn arguments
RooCmdArg Normalization(Double_t scaleFactor, RooAbsPdf::ScaleType scaleType) ;
RooCmdArg Components(const RooArgSet& compSet) ;
RooCmdArg Components(const char* compSpec) ;

// RooAbsData::plotOn arguments
RooCmdArg Cut(const char* cutSpec) ;
RooCmdArg Cut(const RooFormulaVar& cutVar) ;
RooCmdArg Binning(const RooAbsBinning& binning) ;
RooCmdArg Binning(const char* binningName) ;
RooCmdArg Binning(Int_t nBins, Double_t xlo=0., Double_t xhi=0.) ;
RooCmdArg MarkerStyle(Style_t style) ;
RooCmdArg MarkerSize(Size_t size) ;
RooCmdArg MarkerColor(Color_t color) ;
RooCmdArg CutRange(const char* rangeName) ;
RooCmdArg XErrorSize(Double_t width) ;
RooCmdArg RefreshNorm() ;
RooCmdArg Efficiency(const RooCategory& cat) ;
RooCmdArg Rescale(Double_t factor) ;

// RooDataHist::ctor arguments
RooCmdArg Weight(Double_t wgt) ;
RooCmdArg Index(RooCategory& icat) ;
RooCmdArg Import(const char* state, TH1& histo) ;
RooCmdArg Import(TH1& histo) ;

// RooDataSet::ctor arguments
RooCmdArg WeightVar(const char* name) ;
RooCmdArg WeightVar(const RooRealVar& arg) ;
RooCmdArg Import(const char* state, RooDataSet& data) ;
RooCmdArg Import(RooDataSet& data) ;
RooCmdArg Import(TTree& tree) ;
RooCmdArg StoreError(const RooArgSet& aset) ; 
RooCmdArg StoreAsymError(const RooArgSet& aset) ; 

// RooChi2Var::ctor arguments
RooCmdArg Extended(Bool_t flag=kTRUE) ;
RooCmdArg DataError(RooDataHist::ErrorType) ;
RooCmdArg NumCPU(Int_t nCPU, Bool_t interleave=kFALSE) ;

// RooAbsPdf::printLatex arguments
RooCmdArg Columns(Int_t ncol) ;
RooCmdArg OutputFile(const char* fileName) ;
RooCmdArg Format(const char* format, Int_t sigDigit) ;
RooCmdArg Format(const char* what, const RooCmdArg& arg1=RooCmdArg::none(), const RooCmdArg& arg2=RooCmdArg::none(),
                 const RooCmdArg& arg3=RooCmdArg::none(),const RooCmdArg& arg4=RooCmdArg::none(),
                 const RooCmdArg& arg5=RooCmdArg::none(),const RooCmdArg& arg6=RooCmdArg::none(),
                 const RooCmdArg& arg7=RooCmdArg::none(),const RooCmdArg& arg8=RooCmdArg::none()) ;
RooCmdArg Sibling(const RooAbsCollection& sibling) ;

// RooAbsRealLValue::frame arguments
RooCmdArg Title(const char* name) ;
RooCmdArg Bins(Int_t nbin) ;
RooCmdArg AutoSymRange(const RooAbsData& data, Double_t marginFactor=0.1) ;
RooCmdArg AutoRange(const RooAbsData& data, Double_t marginFactor=0.1) ;

// RooAbsData::createHistogram arguments
RooCmdArg AutoSymBinning(Int_t nbins=100, Double_t marginFactor=0.1) ; 
RooCmdArg AutoBinning(Int_t nbins=100, Double_t marginFactor=0.1) ;


// RooAbsData::reduce arguments
RooCmdArg SelectVars(const RooArgSet& vars) ;
RooCmdArg EventRange(Int_t nStart, Int_t nStop) ;

// RooAbsPdf::fitTo arguments
RooCmdArg FitOptions(const char* opts) ;
RooCmdArg Optimize(Bool_t flag=kTRUE) ;
RooCmdArg ProjectedObservables(const RooArgSet& set) ; // obsolete, for backward compatibility
RooCmdArg ConditionalObservables(const RooArgSet& set) ;
RooCmdArg Verbose(Bool_t flag=kTRUE) ;
RooCmdArg Save(Bool_t flag=kTRUE) ;
RooCmdArg Timer(Bool_t flag=kTRUE) ;
RooCmdArg PrintLevel(Int_t code) ;
RooCmdArg Warnings(Bool_t flag=kTRUE) ; 
RooCmdArg Strategy(Int_t code) ;
RooCmdArg InitialHesse(Bool_t flag=kTRUE) ;
RooCmdArg Hesse(Bool_t flag=kTRUE) ;
RooCmdArg Minos(Bool_t flag=kTRUE) ;
RooCmdArg Minos(const RooArgSet& minosArgs) ;
RooCmdArg SplitRange(Bool_t flag=kTRUE) ;
RooCmdArg SumCoefRange(const char* rangeName) ;
RooCmdArg Constrain(const RooArgSet& params) ;
RooCmdArg ExternalConstraints(const RooArgSet& constraintPdfs) ;
RooCmdArg PrintEvalErrors(Int_t numErrors) ;
RooCmdArg EvalErrorWall(Bool_t flag) ;
RooCmdArg SumW2Error(Bool_t flag) ;
RooCmdArg CloneData(Bool_t flag) ;
RooCmdArg Integrate(Bool_t flag) ;

// RooAbsPdf::paramOn arguments
RooCmdArg Label(const char* str) ;
RooCmdArg Layout(Double_t xmin, Double_t xmax=0.99, Double_t ymin=0.95) ;
RooCmdArg Parameters(const RooArgSet& params) ;
RooCmdArg ShowConstants(Bool_t flag=kTRUE) ;

// RooTreeData::statOn arguments
RooCmdArg What(const char* str) ;

// RooProdPdf::ctor arguments
RooCmdArg Conditional(const RooArgSet& pdfSet, const RooArgSet& depSet, Bool_t depsAreCond=kFALSE) ;

// RooAbsPdf::generate arguments
RooCmdArg ProtoData(const RooDataSet& protoData, Bool_t randomizeOrder=kFALSE, Bool_t resample=kFALSE) ;
RooCmdArg NumEvents(Int_t numEvents) ;

// RooAbsRealLValue::createHistogram arguments
RooCmdArg YVar(const RooAbsRealLValue& var, const RooCmdArg& arg=RooCmdArg::none()) ;
RooCmdArg ZVar(const RooAbsRealLValue& var, const RooCmdArg& arg=RooCmdArg::none()) ;
RooCmdArg AxisLabel(const char* name) ;
RooCmdArg Scaling(Bool_t flag) ; 

// RooAbsReal::createIntegral arguments
RooCmdArg NormSet(const RooArgSet& nset) ;
RooCmdArg NumIntConfig(const RooNumIntConfig& cfg) ;

// RooMCStudy::ctor arguments
RooCmdArg Silence(Bool_t flag=kTRUE) ;
RooCmdArg FitModel(RooAbsPdf& pdf) ;
RooCmdArg FitOptions(const RooCmdArg& arg1                ,const RooCmdArg& arg2=RooCmdArg::none(),
                     const RooCmdArg& arg3=RooCmdArg::none(),const RooCmdArg& arg4=RooCmdArg::none(),
                     const RooCmdArg& arg5=RooCmdArg::none(),const RooCmdArg& arg6=RooCmdArg::none()) ;
RooCmdArg Binned(Bool_t flag=kTRUE) ;

// RooMCStudy::plot* arguments
RooCmdArg Frame(const RooCmdArg& arg1                ,const RooCmdArg& arg2=RooCmdArg::none(),
                const RooCmdArg& arg3=RooCmdArg::none(),const RooCmdArg& arg4=RooCmdArg::none(),
                const RooCmdArg& arg5=RooCmdArg::none(),const RooCmdArg& arg6=RooCmdArg::none()) ;
RooCmdArg FrameBins(Int_t nbins) ;
RooCmdArg FrameRange(Double_t xlo, Double_t xhi) ;
RooCmdArg FitGauss(Bool_t flag=kTRUE) ;

// RooRealVar::format arguments
RooCmdArg AutoPrecision(Int_t ndigit=2) ;
RooCmdArg FixedPrecision(Int_t ndigit=2) ;
RooCmdArg TLatexStyle(Bool_t flag=kTRUE) ;
RooCmdArg LatexStyle(Bool_t flag=kTRUE) ;
RooCmdArg LatexTableStyle(Bool_t flag=kTRUE) ;
RooCmdArg VerbatimName(Bool_t flag=kTRUE) ;

// RooMsgService::addReportingStream arguments
RooCmdArg Topic(Int_t topic) ;
RooCmdArg ObjectName(const char* name) ;
RooCmdArg ClassName(const char* name) ;
RooCmdArg BaseClassName(const char* name) ;
RooCmdArg TagName(const char* name) ;
RooCmdArg OutputStream(ostream& os) ;
RooCmdArg Prefix(Bool_t flag) ;
RooCmdArg Color(Color_t color) ;


// RooWorkspace::import() arguments
RooCmdArg RenameConflictNodes(const char* suffix) ;
RooCmdArg RenameAllNodes(const char* suffix) ; 
RooCmdArg RenameAllVariables(const char* suffix) ; 
RooCmdArg RenameVariable(const char* inputName, const char* outputName) ;
RooCmdArg Rename(const char* suffix) ;
RooCmdArg RecycleConflictNodes(Bool_t flag=kTRUE) ;

// RooSimCloneTool::build() arguments
RooCmdArg SplitParam(const char* varname, const char* catname) ;
RooCmdArg SplitParam(const RooRealVar& var, const RooAbsCategory& cat) ;
RooCmdArg SplitParamConstrained(const char* varname, const char* catname, const char* rsname) ;
RooCmdArg SplitParamConstrained(const RooRealVar& var, const RooAbsCategory& cat, const char* rsname) ;
RooCmdArg Restrict(const char* catName, const char* stateNameList) ;

// RooAbsPdf::createCdf() arguments
RooCmdArg SupNormSet(const RooArgSet& nset) ; 
RooCmdArg ScanParameters(Int_t nbins,Int_t intOrder) ;
RooCmdArg ScanNumCdf() ;
RooCmdArg ScanAllCdf() ;
RooCmdArg ScanNoCdf() ;

RooConstVar& RooConst(Double_t val) ; 

}

namespace RooFitShortHand {

RooArgSet S(const RooAbsArg& v1) ;
RooArgSet S(const RooAbsArg& v1, const RooAbsArg& v2) ;
RooArgSet S(const RooAbsArg& v1, const RooAbsArg& v2, const RooAbsArg& v3) ;
RooArgSet S(const RooAbsArg& v1, const RooAbsArg& v2, const RooAbsArg& v3, const RooAbsArg& v4) ;
RooArgSet S(const RooAbsArg& v1, const RooAbsArg& v2, const RooAbsArg& v3, const RooAbsArg& v4, const RooAbsArg& v5) ;
RooArgSet S(const RooAbsArg& v1, const RooAbsArg& v2, const RooAbsArg& v3, const RooAbsArg& v4, const RooAbsArg& v5, 
            const RooAbsArg& v6) ;
RooArgSet S(const RooAbsArg& v1, const RooAbsArg& v2, const RooAbsArg& v3, const RooAbsArg& v4, const RooAbsArg& v5, 
            const RooAbsArg& v6, const RooAbsArg& v7) ;
RooArgSet S(const RooAbsArg& v1, const RooAbsArg& v2, const RooAbsArg& v3, const RooAbsArg& v4, const RooAbsArg& v5, 
            const RooAbsArg& v6, const RooAbsArg& v7, const RooAbsArg& v8) ;
RooArgSet S(const RooAbsArg& v1, const RooAbsArg& v2, const RooAbsArg& v3, const RooAbsArg& v4, const RooAbsArg& v5, 
            const RooAbsArg& v6, const RooAbsArg& v7, const RooAbsArg& v8, const RooAbsArg& v9) ;

RooArgList L(const RooAbsArg& v1) ;
RooArgList L(const RooAbsArg& v1, const RooAbsArg& v2) ;
RooArgList L(const RooAbsArg& v1, const RooAbsArg& v2, const RooAbsArg& v3) ;
RooArgList L(const RooAbsArg& v1, const RooAbsArg& v2, const RooAbsArg& v3, const RooAbsArg& v4) ;
RooArgList L(const RooAbsArg& v1, const RooAbsArg& v2, const RooAbsArg& v3, const RooAbsArg& v4, const RooAbsArg& v5) ;
RooArgList L(const RooAbsArg& v1, const RooAbsArg& v2, const RooAbsArg& v3, const RooAbsArg& v4, const RooAbsArg& v5, 
             const RooAbsArg& v6) ;
RooArgList L(const RooAbsArg& v1, const RooAbsArg& v2, const RooAbsArg& v3, const RooAbsArg& v4, const RooAbsArg& v5, 
             const RooAbsArg& v6, const RooAbsArg& v7) ;
RooArgList L(const RooAbsArg& v1, const RooAbsArg& v2, const RooAbsArg& v3, const RooAbsArg& v4, const RooAbsArg& v5, 
             const RooAbsArg& v6, const RooAbsArg& v7, const RooAbsArg& v8) ;
RooArgList L(const RooAbsArg& v1, const RooAbsArg& v2, const RooAbsArg& v3, const RooAbsArg& v4, const RooAbsArg& v5, 
             const RooAbsArg& v6, const RooAbsArg& v7, const RooAbsArg& v8, const RooAbsArg& v9) ;

RooConstVar& C(Double_t value) ;

} // End namespace ShortHand

#endif
