#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TInterpreter.h>
#include <TString.h>
#include <TObjString.h>
#include <TObjArray.h>
#include <TMath.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH1D.h>
#include <TF1.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TDatabasePDG.h>
#include <TGraph.h>

#include "AliAODRecoDecayHF.h"
#include "AliRDHFCuts.h"
#include "AliRDHFCutsDplustoKpipi.h"
#include "AliRDHFCutsDStartoKpipi.h"
#include "AliRDHFCutsDstoKKpi.h"
#include "AliRDHFCutsD0toKpi.h"
#include "AliHFMassFitter.h"
#include "AliNormalizationCounter.h"
#endif


// MACRO to perform fits to D meson invariant mass spectra
// and store raw yields and cut object into a root output file
// Origin: F. Prino (prino@to.infn.it)
// D0: C. Bianchin (cbianchi@pd.infn.it)



//
enum {kD0toKpi, kDplusKpipi, kDStarD0pi, kDsKKpi, kDStarD0piMult, kDStarD0piToyMC};
enum {kBoth, kParticleOnly, kAntiParticleOnly};
enum {kExpo=0, kLinear, kPol2, kThrExpo=5};
enum {kGaus=0, kDoubleGaus};


// Common variables: to be configured by the user
const Int_t nPtBins=1;
const Int_t nRebin=3;
//Double_t ptlims[nPtBins+1]={0.,9999.};//{0.,0.5,1.,2.,3.,4.,5.,6.,8.,12.,16.,20.,24};
Int_t rebin[nRebin]={1,2,3};//};//{30,10,10,10,10,10,10,10,10,10,10,10};
Int_t firstUsedBin[nRebin]={51,51,51};//};//{-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1}; // -1 uses all bins, >=1 to set the lower bin to be accepted from original histo
Double_t InitSigma[nRebin]={0.00065,0.00065,0.00065};//};//,0.00065,0.00065};//0.00045;
Int_t multbin = 5;
Int_t ptbin = 4;
multbin--;
ptbin--;
TString hmassName = Form("hmass_%d_%d",multbin,ptbin);
//hmassName += multbin;
//hmassName += "_";
//hmassName += ptbin;

TString hmassNameInt = Form("hmass_0_%d",ptbin);
//hmassNameInt += multbin;

TString hsignalName = Form("hSignal_%d",multbin);
//hsignalName += multbin;

TString hsignalNameInt = "hSignal_0";
TString suffix="";


//const Int_t nPtBins=7;//6;
//Double_t ptlims[nPtBins+1]={1.,2.,3.,4.,5.,6.,8.,12.};
//Int_t rebin[nPtBins+1]={8,6,10,10,10,10,10,10}; //for looser cuts
//Int_t rebin[nPtBins+1]={10,10,10,10,10,10,10,10}; //for Chiara's cuts
Int_t typeb[2]={kThrExpo,kExpo};
Int_t types=kGaus;
Int_t optPartAntiPart=kBoth;
Int_t factor4refl=0;
const Int_t nMassMin=3;
const Int_t nMassMax=9;
Double_t minMassForFit[nMassMin]={0.139,0.14.0.141};//};//};//};
Double_t maxMassForFit[nMassMax]={0.197,0.188,0.187,0.186,0.184,0.182,0.178,0.176,0.170};
//multbin_3 ptbin 4{0.197,0.189,0.188,0.187,0.186,0.182,0.179,0.178,0.177,0.173,0.172,0.170};
//multbin_3 ptbin_3{0.199,0.194,0.189,0.188,0.185,0.184,0.183,0.181,0.18,0.179,0.178,0.176,0.175,0.174,0.173,0.172,0.170};
//multbin_2 ptbin_3{0.199,0.194,0.192,0.189,0.188,0.185,0.184,0.183,0.181,0.18,0.179,0.178,0.176,0.175,0.174,0.173,0.172,0.170};
//multbin_2 ptbin_2{0.19,0.189,0.188,0.183,0.182,0.181,0.18,0.175,0.171,0.170};
//multbin_2 ptbin_1{0.170,0.171,0.172,0.174,0.175,0.177,0.179,0.181,0.182,0.184,0.185,0.188,0.189,0.19,0.191,0.192,0.193,0.194,0.198};
//multbin_3 ptbin_1{0.170,0.171,0.172,0.174,0.175,0.176,0.182,0.184,0.185,0.188,0.189,0.192,0.194,0.197};
//multbin_4 ptbin_1{0.191,0.188,0.184,0.182,0.179,0.176,0.175,0.174,0.172,0.170};
//multbin_5 ptbin_1{0.197,0.193,0.191,0.19,0.189,0.179,0.175,0.174,0.172,0.171,0.170};
Double_t minBin = 0.07;
Double_t maxBin = 0.1;
//Int_t nSigmaRangeForCounting=3;
Float_t massRangeForCounting=0.00135; // GeV
TH2F* hPtMass=0x0;
Double_t nEventsForNorm=0.;
Bool_t fixPeakSigma = kFALSE;
//for D0only
Bool_t cutsappliedondistr=kFALSE;
const Int_t nsamples=1;//3;
//Int_t nevents[nsamples]={8.5635859e+07 /*LHC10b+cpart*/,8.9700624e+07/*LHC10dpart*/};
//Int_t nevents[nsamples]={9.0374946e+07 /*LHC10b+c*/,9.7593785e+07/*LHC10d*/};
//Int_t nevents[nsamples]={1.1777812e+08 /*LHC10dnewTPCpid*/,9.7593785e+07/*LHC10d*/};
//Int_t nevents[nsamples]={1.1777812e+08 /*LHC10dnewTPCpid*/,9.7593785e+07/*LHC10d*/,9.0374946e+07 /*LHC10b+c*/};
Int_t nevents[nsamples]={300};//8.1611334400000000e+08 /*LHC16l*/ /*LHC10dnewTPCpid*/}; //for compare fit
//Int_t nevents[nsamples]={1. /*LHC10dnewTPCpid*/,1 /*LHC10dnewTPCpidrescale*/};
// Functions

Bool_t LoadDplusHistos(TObjArray* listFiles, TH1F** hMass);
Bool_t LoadDsHistos(TObjArray* listFiles, TH1F** hMass);
Bool_t LoadD0toKpiHistos(TObjArray* listFiles, TH1F** hMass);
Bool_t LoadDstarD0piHistos(TString* fileName, TH1F** hMass);
Bool_t LoadDstarD0piHistosMult(TString fileName, TH1F** hMass, Float_t& centralvalue);
Bool_t LoadDstarD0piHistosIntMult(TString fileNameInt, TH1F** hMass, Float_t& centralvalue);
Bool_t LoadDstarToyMC(TString fileName, TH1F** hMass);

TH1F* RebinHisto(TH1F* hOrig, Int_t reb, Int_t firstUse=-1);



void YieldExtractionSystematic(Int_t analysisType=kDStarD0piMult,
		    TString fileName="Good_224_RawYield_Mult_Both_pp13TeV_BCin3.0Sigma_ThrExpo.root",
        TString fileNameInt="Good_RawYield_IntMult_Both_pp13TeV_BCin3.0Sigma_ThrExpo.root",
	       ){
  //

  gInterpreter->ExecuteMacro("$ALICE_PHYSICS/PWGHF/vertexingHF/macros/LoadLibraries.C");
  gStyle->SetOptTitle(0);

  /*
  TObjArray* listFiles=new TObjArray();
  if(fileNameb!="") listFiles->AddLast(new TObjString(fileNameb.Data()));
  if(fileNamec!="") listFiles->AddLast(new TObjString(fileNamec.Data()));
  if(fileNamed!="") listFiles->AddLast(new TObjString(fileNamed.Data()));
  if(fileNamee!="") listFiles->AddLast(new TObjString(fileNamee.Data()));
  if(listFiles->GetEntries()==0){
    printf("Missing file names in input\n");
    return;
  }
  */

  TH1F** hmass=new TH1F*[nRebin];
  TH1F** hmassint=new TH1F*[nRebin];
  for(Int_t i=0;i<nRebin;i++) {
    hmass[i]=0x0;
    hmassint[i]=0x0;
  }

  Float_t massD, massD0_fDstar, centralvaluemult, centralvalueint;
  Bool_t retCode;
  if(analysisType==kD0toKpi){
    retCode=LoadD0toKpiHistos(listFiles,hmass);
    massD=TDatabasePDG::Instance()->GetParticle(421)->Mass();
  }
  else if(analysisType==kDplusKpipi){
    retCode=LoadDplusHistos(listFiles,hmass);
    massD=TDatabasePDG::Instance()->GetParticle(411)->Mass();
  }
 else if(analysisType==kDStarD0pi){
    retCode=LoadDstarD0piHistos(listFiles,hmass);
    massD=TDatabasePDG::Instance()->GetParticle(413)->Mass();
    massD0_fDstar=TDatabasePDG::Instance()->GetParticle(421)->Mass();
    massD =massD-massD0_fDstar;
  }
  else if(analysisType==kDStarD0piMult){
    retCode=LoadDstarD0piHistosMult(fileName,hmass,centralvaluemult)*LoadDstarD0piHistosIntMult(fileNameInt,hmassint,centralvalueint);
    massD=TDatabasePDG::Instance()->GetParticle(413)->Mass();
    massD0_fDstar=TDatabasePDG::Instance()->GetParticle(421)->Mass();
    massD =massD-massD0_fDstar;
  }
  else if(analysisType==kDStarD0piToyMC){
    retCode=LoadDstarToyMC(fileName,hmass);
    massD=TDatabasePDG::Instance()->GetParticle(413)->Mass();
    massD0_fDstar=TDatabasePDG::Instance()->GetParticle(421)->Mass();
    massD =massD-massD0_fDstar;
  }
  else if(analysisType==kDsKKpi){
    retCode=LoadDsHistos(listFiles,hmass);
    massD=TDatabasePDG::Instance()->GetParticle(431)->Mass();
  }
  else{
    printf("Wrong analysis type parameter\n");
    return;
  }
  if(!retCode){
    printf("ERROR in reading input files\n");
    return;
  } 
  
  const Int_t totalvar=nRebin*2*nMassMax*nMassMin;

  TH1D* hCntSig1=new TH1D("hCntSig1","hCntSig1",totalvar,0.5,totalvar+0.5);
  TH1D* hCntSig2=new TH1D("hCntSig2","hCntSig2",totalvar,0.5,totalvar+0.5);
  TH1D* hNDiffCntSig1=new TH1D("hNDiffCntSig1","hNDiffCntSig1",totalvar,0.5,totalvar+0.5);
  TH1D* hNDiffCntSig2=new TH1D("hNDiffCntSig2","hNDiffCntSig2",totalvar,0.5,totalvar+0.5);
  TH1D* hSignal=new TH1D("hSignal","hSignal",totalvar,0.5,totalvar+0.5);
  TH1D* hRelErrSig=new TH1D("hRelErrSig","hRelErrSig",totalvar,0.5,totalvar+0.5);
  TH1D* hInvSignif=new TH1D("hInvSignif","hInvSignif",totalvar,0.5,totalvar+0.5);
  TH1D* hBackground=new TH1D("hBackground","hBackground",totalvar,0.5,totalvar+0.5);
  TH1D* hBackgroundNormSigma=new TH1D("hBackgroundNormSigma","hBackgroundNormSigma",totalvar,0.5,totalvar+0.5);
  TH1D* hSignificance=new TH1D("hSignificance","hSignificance",totalvar,0.5,totalvar+0.5);
  TH1D* hMass=new TH1D("hMass","hMass",totalvar,0.5,totalvar+0.5);
  TH1D* hSigma=new TH1D("hSigma","hSigma",totalvar,0.5,totalvar+0.5);
  TH1D* hCntRMS=new TH1D("hCntRMS","hCntRMS",totalvar,0.5,totalvar+0.5);
  TH1I* hCntRMSBin=new TH1I("hCntRMSBin","hCntRMSBin",30,30,60);
  TH1D* hSigDistribution=new TH1D("hSigDistribution","hSigDistribution",100,minBin,maxBin);


  Int_t nMassBins=hmass[0]->GetNbinsX();
  /*
  Double_t hmin=TMath::Max(minMassForFit,hmass[0]->GetBinLowEdge(2));
  Double_t hmax=TMath::Min(maxMassForFit,hmass[0]->GetBinLowEdge(nMassBins-2)+hmass[0]->GetBinWidth(nMassBins-2));
  */
  Float_t minBinSum=hmass[0]->FindBin(massD-massRangeForCounting);
  Float_t maxBinSum=hmass[0]->FindBin(massD+massRangeForCounting);
  Int_t iPad=1;

  TF1* funBckStore1=0x0;
  TF1* funBckStore2=0x0;
  TF1* funBckStore3=0x0;

  AliHFMassFitter** fitter=new AliHFMassFitter*[totalvar];
  AliHFMassFitter** fitterint=new AliHFMassFitter*[totalvar];
  Double_t arrchisquare[totalvar];
  Double_t arrchisquareint[totalvar];
  /*
  TCanvas* c1= new TCanvas("c1","MassSpectra",1000,800);
  Int_t nx=3, ny=2;
  if(totalvar>6){
    nx=4;
    ny=3;
  }
  if(totalvar>12){
    nx=5;
    ny=4;
  }
  c1->Divide(nx,ny);
  */
  Double_t sig,errsig,s,errs,b,errb,sigint,errsigint,sint,errsint,bint,errbint,dvar[totalvar];
  Int_t rebinnew[nRebin];
  Int_t icount=0,xcount=0;
  for (Int_t iMinMass=0; iMinMass<nMassMin; iMinMass++) {
  for (Int_t iMaxMass=0; iMaxMass<nMassMax; iMaxMass++) {
  for (Int_t itypeb=0; itypeb<2; itypeb++) {
  for (Int_t iBin=0; iBin<nRebin; iBin++){
    //c1->cd(iPad++);
    dvar[icount]=icount;
    Int_t origNbins=hmass[0]->GetNbinsX();
    TH1F* hRebinned=RebinHisto(hmass[0],rebin[iBin],firstUsedBin[iBin]);
    TH1F* hRebinnedInt=RebinHisto(hmassint[0],rebin[iBin],firstUsedBin[iBin]);
    Double_t hmin=TMath::Max(minMassForFit[iMinMass],hRebinned->GetBinLowEdge(2));
    Double_t hmax=TMath::Min(maxMassForFit[iMaxMass],hRebinned->GetBinLowEdge(hRebinned->GetNbinsX()));
    fitter[icount]=new AliHFMassFitter(hRebinned,hmin,hmax,1,typeb[itypeb],types);
    fitterint[icount]= new AliHFMassFitter(hRebinnedInt,hmin,hmax,1,typeb[itypeb],types);
    rebinnew[iBin]=origNbins/fitter[icount]->GetBinN();
    fitter[icount]->SetReflectionSigmaFactor(factor4refl);
    fitterint[icount]->SetReflectionSigmaFactor(factor4refl);
    fitter[icount]->SetInitialGaussianMean(massD);
    fitterint[icount]->SetInitialGaussianMean(massD);
    //if(analysisType==kDStarD0pi || ) 
    fitter[icount]->SetInitialGaussianSigma(InitSigma[iBin]);
    fitterint[icount]->SetInitialGaussianSigma(InitSigma[iBin]);
    cout << "before fitting" << endl;
    Bool_t out=fitter[icount]->MassFitter(0);
    Bool_t outint=fitterint[icount]->MassFitter(0);
    /*if(!out) {
      fitter[icount]->GetHistoClone()->Draw();
      continue;
    }*/
    Double_t mass=fitter[icount]->GetMean();
    Double_t massint=fitterint[icount]->GetMean();
    Double_t sigma=fitter[icount]->GetSigma();
    Double_t sigmaint=fitterint[icount]->GetSigma();
    Double_t fitprob=fitter[icount]->GetFitProbability();
    arrchisquare[icount]=fitter[icount]->GetReducedChiSquare();
    arrchisquareint[icount]=fitterint[icount]->GetReducedChiSquare();
    TF1* fB1=fitter[icount]->GetBackgroundFullRangeFunc();
    TF1* fB2=fitter[icount]->GetBackgroundRecalcFunc();
    TF1* fM=fitter[icount]->GetMassFunc();
    TF1* fB1int=fitterint[icount]->GetBackgroundFullRangeFunc();
    TF1* fB2int=fitterint[icount]->GetBackgroundRecalcFunc();
    TF1* fMint=fitterint[icount]->GetMassFunc();
    
    if(icount==0 && fB1) funBckStore1=new TF1(*fB1);
    if(icount==0 && fB2) funBckStore2=new TF1(*fB2);
    if(icount==0 && fM) funBckStore3=new TF1(*fM);

    if(icount==0 && fB1int) funBckStore1int=new TF1(*fB1int);
    if(icount==0 && fB2int) funBckStore2int=new TF1(*fB2int);
    if(icount==0 && fMint) funBckStore3int=new TF1(*fMint);
    //fitter[icount]->DrawHere(gPad);    
    fitter[icount]->Signal(3,s,errs);
    fitter[icount]->Background(3,b,errb);
    fitter[icount]->Significance(3,sig,errsig);

    fitterint[icount]->Signal(3,sint,errsint);
    fitterint[icount]->Background(3,bint,errbint);
    fitterint[icount]->Significance(3,sigint,errsigint);

    Double_t ry=fitter[icount]->GetRawYield();
    Double_t ery=fitter[icount]->GetRawYieldError();

    Double_t ryint=fitterint[icount]->GetRawYield();
    Double_t eryint=fitterint[icount]->GetRawYieldError();

    Float_t cntSig1=0.;
    Float_t cntSig2=0.;
    Float_t cntErr=0.;
    Float_t cntRMS=0.;
    for(Int_t iMB=minBinSum; iMB<=maxBinSum; iMB++){
      Float_t bkg1=fB1 ? fB1->Eval(hmass[0]->GetBinCenter(iMB))/rebinnew[iBin] : 0;
      Float_t bkg2=fB2 ? fB2->Eval(hmass[0]->GetBinCenter(iMB))/rebinnew[iBin] : 0;
      cntSig1+=(hmass[0]->GetBinContent(iMB)-bkg1);
      cntSig2+=(hmass[0]->GetBinContent(iMB)-bkg2);
      cntErr+=(hmass[0]->GetBinContent(iMB));
      cntRMS+=TMath::Power((hmass[0]->GetBinContent(iMB)-bkg1),2.);
    }
    /*
    hCntSig1->SetBinContent(icount+1,cntSig1);
    hCntSig1->SetBinError(icount+1,TMath::Sqrt(cntErr));
    hNDiffCntSig1->SetBinContent(icount+1,(s-cntSig1)/s);
    hNDiffCntSig1->SetBinError(icount+1,TMath::Sqrt(cntErr)/s);
    hCntSig2->SetBinContent(icount+1,cntSig2);
    hNDiffCntSig2->SetBinContent(icount+1,(s-cntSig2)/s);
    hNDiffCntSig2->SetBinError(icount+1,TMath::Sqrt(cntErr)/s);
    hCntSig2->SetBinError(icount+1,TMath::Sqrt(cntErr));
    */
    if (sigma>2*sigmaint/3 && sigma<3*sigmaint/2 && sig>3 && fitprob>0.05) {
      hSignal->SetBinContent(xcount+1,ry/ryint);
      hSignal->SetBinError(xcount+1,0);
      hSigDistribution->Fill(ry/ryint);
      xcount++;
    }
    /*
    hRelErrSig->SetBinContent(icount+1,errs/s);
    hInvSignif->SetBinContent(icount+1,1/sig);
    hInvSignif->SetBinError(icount+1,errsig/(sig*sig));
    hBackground->SetBinContent(icount+1,b); //consider sigma
    hBackground->SetBinError(icount+1,errb);
    hBackgroundNormSigma->SetBinContent(icount+1,b/(3*fitter[iBin]->GetSigma())*(3*0.012)); //consider sigma
    hBackgroundNormSigma->SetBinError(icount+1,errb);
    hSignificance->SetBinContent(icount+1,sig);
    hSignificance->SetBinError(icount+1,errsig);
    hMass->SetBinContent(icount+1,mass);
    hMass->SetBinError(icount+1,0.0001);
    hSigma->SetBinContent(icount+1,sigma);
    hSigma->SetBinError(icount+1,0.0001);
    hCntRMS->SetBinContent(icount+1,TMath::Sqrt(cntRMS/(maxBinSum-minBinSum+1)));//fitter[icount]->GetMean()
    hCntRMS->SetBinError(icount+1,0.0001);
    hCntRMSBin->Fill(TMath::Sqrt(cntRMS/(maxBinSum-minBinSum+1)));*/
    icount++;
  }
  }
  }
  }

  /*
  c1->cd(1); // is some cases the fitting function of 1st bin get lost
  funBckStore1->Draw("same");
  funBckStore2->Draw("same");
  funBckStore3->Draw("same");
  */
/*
  TCanvas *cpar=new TCanvas("cpar","Fit params",1200,600);
  cpar->Divide(2,1);
  cpar->cd(1);
  hMass->SetMarkerStyle(20);
  hMass->Draw("PE");
  hMass->GetXaxis()->SetTitle("Variation");
  hMass->GetXaxis()->SetTitle("Mass (GeV/c^{2})");
  cpar->cd(2);
  hSigma->SetMarkerStyle(20);
  hSigma->Draw("PE");
  hSigma->GetXaxis()->SetTitle("Variation");
  hSigma->GetXaxis()->SetTitle("Sigma (GeV/c^{2})");
*/
  
  cout << centralvaluemult << " and " << centralvalueint << endl;
  centralvaluemult/=centralvalueint;
  Float_t updist = hSignal->GetBinContent(hSignal->GetMaximumBin()) - centralvaluemult;
  Float_t lowdist=1., lowdistemp;
  for (Int_t i=0;i<totalvar;i++) {
    if (hSignal->GetBinContent(i+1)!=0) {
      lowdistemp = lowdist;
      lowdist = hSignal->GetBinContent(i+1);
      lowdist = TMath::Min(lowdist,lowdistemp);
    }
  }
  lowdist = centralvaluemult - lowdist;
  Float_t sysuncplus = 100*updist/centralvaluemult;
  Float_t sysuncminus = 100*lowdist/centralvaluemult;
  cout << "Systematic uncertainty = " << "+" << sysuncplus << "%" << "-" << sysuncminus << "%" << endl;
  gStyle->SetOptStat(0);
  TF1 *fCentralValue = new TF1("fCentralValue",Form("%f",centralvaluemult),0,totalvar);
  TCanvas* csig=new TCanvas("csig","Results",1200,600);
  csig->Divide(2,1);
  csig->cd(1);
  hSignal->SetMarkerStyle(20);
  hSignal->SetMarkerColor(4);
  hSignal->SetLineColor(4);
  hSignal->GetXaxis()->SetTitle("Trials");
  hSignal->GetYaxis()->SetTitle("\\mathcurl{R}");
  hSignal->Draw();
  fCentralValue->Draw("SAME");
  /*
  hCntSig1->SetMarkerStyle(26);
  hCntSig1->SetMarkerColor(2);
  hCntSig1->SetLineColor(2);
  hCntSig1->Draw("PSAME");
  hCntSig2->SetMarkerStyle(29);
  hCntSig2->SetMarkerColor(kGray+1);
  hCntSig2->SetLineColor(kGray+1);
  hCntSig2->Draw("PSAME");
  */

  /*
  TLegend* leg=new TLegend(0.4,0.7,0.89,0.89);
  leg->SetFillColor(0);
  TLegendEntry* ent=leg->AddEntry(hSignal,"From Fit","PL");
  ent->SetTextColor(hSignal->GetMarkerColor());
  
  ent=leg->AddEntry(hCntSig1,"From Counting1","PL");
  ent->SetTextColor(hCntSig1->GetMarkerColor());
  ent=leg->AddEntry(hCntSig2,"From Counting2","PL");
  ent->SetTextColor(hCntSig2->GetMarkerColor());
  
  leg->Draw();
  */
  csig->cd(2);
  
  hSigDistribution->Draw();
  hSigDistribution->GetXaxis()->SetTitle("\\mathcurl{R}");
  hSigDistribution->GetYaxis()->SetTitle("count");
  
  /*
  hBackground->SetMarkerStyle(20);
  hBackground->Draw("P");
  hBackground->GetXaxis()->SetTitle("Trials");
  hBackground->GetYaxis()->SetTitle("Background");
  */
  /*
  csig->cd(3);
  hSignificance->SetMarkerStyle(20);
  hSignificance->Draw("P");
  hSignificance->GetXaxis()->SetTitle("Variation");
  hSignificance->GetYaxis()->SetTitle("Significance");
*/
/*
  TCanvas* cDiffS=new TCanvas("cDiffS","Difference",1200,600);
  cDiffS->Divide(2,1);
  cDiffS->cd(1);
  hRelErrSig->SetMarkerStyle(20); //fullcircle
  hRelErrSig->SetTitleOffset(1.2);  
  hRelErrSig->SetTitle("Rel Error from Fit;Variation;Signal Relative Error");
  hRelErrSig->Draw("P");
  hInvSignif->SetMarkerStyle(21); //fullsquare
  hInvSignif->SetMarkerColor(kMagenta+1);
  hInvSignif->SetLineColor(kMagenta+1);
  hInvSignif->Draw("PSAMES");
  TLegend* leg2=new TLegend(0.4,0.7,0.89,0.89);
  leg2->SetFillColor(0);
  TLegendEntry* ent2=leg2->AddEntry(hRelErrSig,"From Fit","P");
  ent2->SetTextColor(hRelErrSig->GetMarkerColor());
  ent2=leg2->AddEntry(hInvSignif,"1/Significance","PL");
  ent2->SetTextColor(hInvSignif->GetMarkerColor());
  leg2->Draw();

  cDiffS->cd(2);
  hNDiffCntSig1->SetMarkerStyle(26);
  hNDiffCntSig1->SetMarkerColor(2);
  hNDiffCntSig1->SetLineColor(2);
  hNDiffCntSig1->SetTitle("Cmp Fit-Count;Variation;(S_{fit}-S_{count})/S_{fit}");
  hNDiffCntSig1->Draw("P");
  hNDiffCntSig2->SetMarkerStyle(29);
  hNDiffCntSig2->SetMarkerColor(kGray+1);
  hNDiffCntSig2->SetLineColor(kGray+1);
  hNDiffCntSig2->Draw("PSAME");
  TLegend* leg1=new TLegend(0.4,0.7,0.89,0.89);
  leg1->SetFillColor(0);
  TLegendEntry* ent1=leg1->AddEntry(hNDiffCntSig1,"From Counting1","PL");
  ent1->SetTextColor(hNDiffCntSig1->GetMarkerColor());
  ent1=leg1->AddEntry(hNDiffCntSig2,"From Counting2","PL");
  ent1->SetTextColor(hNDiffCntSig2->GetMarkerColor());
  leg1->Draw();

  Int_t d=0;

  TGraph* grReducedChiSquare=new TGraph(totalvar,dvar,arrchisquare);
  grReducedChiSquare->SetName("grReducedChiSquare");
  grReducedChiSquare->SetTitle("Reduced Chi2;Variation;#tilde{#chi}^{2}");
  TCanvas *cChi2=new TCanvas("cChi2","reduced chi square",600,600);
  cChi2->cd();
  grReducedChiSquare->SetMarkerStyle(21);
  grReducedChiSquare->Draw("AP");

  TCanvas* cbkgNormSigma=new TCanvas("cbkgNormSigma","Background normalized to sigma",400,600);
  cbkgNormSigma->cd();
  hBackgroundNormSigma->SetMarkerStyle(20);
  hBackgroundNormSigma->Draw("P");
  hBackgroundNormSigma->GetXaxis()->SetTitle("Variation");
  hBackgroundNormSigma->GetYaxis()->SetTitle("Background #times 3 #times 0.012/ (3 #times #sigma)");
  hBackgroundNormSigma->Draw();

  TCanvas* cCntRMS=new TCanvas("cCntRMS","RMS of the multitrials",1200,600);
  cCntRMS->Divide(2,1);
  cCntRMS->cd(1);
  hCntRMS->SetMarkerStyle(20);
  hCntRMS->Draw("P");
  hCntRMS->GetXaxis()->SetTitle("Variation");
  hCntRMS->GetYaxis()->SetTitle("RMS");

  cCntRMS->cd(2);
  hCntRMSBin->Draw("");
  hCntRMSBin->GetXaxis()->SetTitle("RMS");
  hCntRMSBin->GetYaxis()->SetTitle("Count");
*/

  TString partname="Both";
  if(optPartAntiPart==kParticleOnly) {
    if(analysisType==kD0toKpi) partname="D0";
    if(analysisType==kDplusKpipi) partname="Dplus";
    if(analysisType==kDsKKpi) partname="Dsplus";
  }
  if(optPartAntiPart==kAntiParticleOnly) {
    if(analysisType==kD0toKpi) partname="D0bar";
    if(analysisType==kDplusKpipi) partname="Dminus";
    if(analysisType==kDsKKpi) partname="Dsminus";
  }

  printf("Events for norm = %f\n",nEventsForNorm);
  TH1F* hNEvents=new TH1F("hNEvents","",1,0.,1.);
  hNEvents->SetBinContent(1,nEventsForNorm);

  TFile* outf=new TFile(Form("RawYield%s.root",partname.Data()),"recreate");
  outf->cd();
  hNEvents->Write();
  hMass->Write();
  hSigma->Write();
  hCntSig1->Write();
  hCntSig2->Write();
  hNDiffCntSig1->Write();
  hNDiffCntSig2->Write();
  hSignal->Write();
  hRelErrSig->Write();
  hInvSignif->Write();
  //hBackground->Write();
  hBackgroundNormSigma->Write();
  hSignificance->Write();
  //grReducedChiSquare->Write();
  //hPtMass->Write();
  outf->Close();
}


Bool_t LoadDplusHistos(TObjArray* listFiles, TH1F** hMass){

  Int_t nFiles=listFiles->GetEntries();
  TList **hlist=new TList*[nFiles];
  AliRDHFCutsDplustoKpipi** dcuts=new AliRDHFCutsDplustoKpipi*[nFiles];

  Int_t nReadFiles=0;
  for(Int_t iFile=0; iFile<nFiles; iFile++){
    TString fName=((TObjString*)listFiles->At(iFile))->GetString();    
    TFile *f=TFile::Open(fName.Data());
    if(!f){
      printf("ERROR: file %s does not exist\n",fName.Data());
      continue;
    }
    printf("Open File %s\n",f->GetName()); 
    TDirectory *dir = (TDirectory*)f->Get("PWG3_D2H_InvMassDplus");
    if(!dir){
      printf("ERROR: directory PWG3_D2H_InvMassDplus not found in %s\n",fName.Data());
      continue;
    }
    hlist[nReadFiles]=(TList*)dir->Get(Form("coutputDplus%s",suffix.Data()));
    TList *listcut = (TList*)dir->Get(Form("coutputDplusCuts%s",suffix.Data()));
    dcuts[nReadFiles]=(AliRDHFCutsDplustoKpipi*)listcut->At(0);
    if(nReadFiles>0){
      Bool_t sameCuts=dcuts[nReadFiles]->CompareCuts(dcuts[0]);
      if(!sameCuts){
	printf("ERROR: Cut objects do not match\n");
	return kFALSE;
      }
    }
    AliNormalizationCounter* c=(AliNormalizationCounter*)dir->Get(Form("coutputDplusNorm%s",suffix.Data()));
    printf("Events for normalization = %f\n",c->GetNEventsForNorm());
    nEventsForNorm+=c->GetNEventsForNorm();
    nReadFiles++;
  }
  if(nReadFiles<nFiles){
    printf("WARNING: not all requested files have been found\n");
  }

  Int_t nPtBinsCuts=dcuts[0]->GetNPtBins();
  printf("Number of pt bins for cut object = %d\n",nPtBins);
  Float_t *ptlimsCuts=dcuts[0]->GetPtBinLimits();
  ptlimsCuts[nPtBinsCuts]=ptlimsCuts[nPtBinsCuts-1]+4.;

  Int_t iFinBin=0;
  for(Int_t i=0;i<nPtBinsCuts;i++){
    if(ptlimsCuts[i]>=ptlims[iFinBin+1]) iFinBin+=1; 
    if(iFinBin>nPtBins) break;
    if(ptlimsCuts[i]>=ptlims[iFinBin] && 
       ptlimsCuts[i+1]<=ptlims[iFinBin+1]){
      for(Int_t iFile=0; iFile<nReadFiles; iFile++){
	TString histoName;
	if(optPartAntiPart==kBoth) histoName.Form("hMassPt%dTC",i);
	else if(optPartAntiPart==kParticleOnly) histoName.Form("hMassPt%dTCPlus",i);
	else if(optPartAntiPart==kAntiParticleOnly) histoName.Form("hMassPt%dTCMinus",i);
	TH1F* htemp=(TH1F*)hlist[iFile]->FindObject(histoName.Data());
	if(!htemp){
	  printf("ERROR: Histogram %s not found\n",histoName.Data());
	  return kFALSE;
	}
	if(!hMass[iFinBin]){
	  hMass[iFinBin]=new TH1F(*htemp);
	}else{
	  hMass[iFinBin]->Add(htemp);
	}
      }
    }
  }
  TString partname="Both";
  if(optPartAntiPart==kParticleOnly) partname="Dplus";
  if(optPartAntiPart==kAntiParticleOnly) partname="Dminus";

  for(Int_t iFile=0; iFile<nReadFiles; iFile++){
    TH2F* htemp2=(TH2F*)hlist[iFile]->FindObject("hPtVsMassTC");
    if(iFile==0){
      hPtMass=new TH2F(*htemp2);
    }else{
      hPtMass->Add(htemp2);
    }
  }

  TFile* outf=new TFile(Form("RawYield%s.root",partname.Data()),"recreate");
  outf->cd();
  dcuts[0]->Write();
  outf->Close();

  return kTRUE;

}


Bool_t LoadDsHistos(TObjArray* listFiles, TH1F** hMass){

  Int_t nFiles=listFiles->GetEntries();
  TList **hlist=new TList*[nFiles];
  AliRDHFCutsDstoKKpi** dcuts=new AliRDHFCutsDstoKKpi*[nFiles];

  Int_t nReadFiles=0;
  for(Int_t iFile=0; iFile<nFiles; iFile++){
    TString fName=((TObjString*)listFiles->At(iFile))->GetString();    
    TFile *f=TFile::Open(fName.Data());
    if(!f){
      printf("ERROR: file %s does not exist\n",fName.Data());
      continue;
    }
    printf("Open File %s\n",f->GetName()); 
    TDirectory *dir = (TDirectory*)f->Get("PWG3_D2H_InvMassDs");
    if(!dir){
      printf("ERROR: directory PWG3_D2H_InvMassDs not found in %s\n",fName.Data());
      continue;
    }
    hlist[nReadFiles]=(TList*)dir->Get("coutputDs");
    TList *listcut = (TList*)dir->Get("coutputDsCuts");
    dcuts[nReadFiles]=(AliRDHFCutsDstoKKpi*)listcut->At(0);
    cout<< dcuts[nReadFiles]<<endl;
    if(nReadFiles>0){
      Bool_t sameCuts=dcuts[nReadFiles]->CompareCuts(dcuts[0]);
      if(!sameCuts){
	printf("ERROR: Cut objects do not match\n");
	return kFALSE;
      }
    }
    nReadFiles++;
  }
  if(nReadFiles<nFiles){
    printf("WARNING: not all requested files have been found\n");
  }

  Int_t nPtBinsCuts=dcuts[0]->GetNPtBins();
  printf("Number of pt bins for cut object = %d\n",nPtBins);
  Float_t *ptlimsCuts=dcuts[0]->GetPtBinLimits();
  ptlimsCuts[nPtBinsCuts]=ptlimsCuts[nPtBinsCuts-1]+4.;

  Int_t iFinBin=0;
  for(Int_t i=0;i<nPtBinsCuts;i++){
    if(ptlimsCuts[i]>=ptlims[iFinBin+1]) iFinBin+=1; 
    if(iFinBin>nPtBins) break;
    if(ptlimsCuts[i]>=ptlims[iFinBin] && 
       ptlimsCuts[i+1]<=ptlims[iFinBin+1]){
      for(Int_t iFile=0; iFile<nReadFiles; iFile++){
	TString histoName;
	if(optPartAntiPart==kBoth) histoName.Form("hMassAllPt%dphi",i);
	else if(optPartAntiPart==kParticleOnly){
	  printf("Particle/Antiparticle not yet enabled for Ds");
	  histoName.Form("hMassAllPt%dphi",i);
	}
	else if(optPartAntiPart==kAntiParticleOnly){
	  printf("Particle/Antiparticle not yet enabled for Ds");
	  histoName.Form("hMassAllPt%dphi",i);
	}
	TH1F* htemp=(TH1F*)hlist[iFile]->FindObject(histoName.Data());
	if(!htemp){
	  printf("ERROR: Histogram %s not found\n",histoName.Data());
	  return kFALSE;
	}
	if(!hMass[iFinBin]){
	  hMass[iFinBin]=new TH1F(*htemp);
	}else{
	  hMass[iFinBin]->Add(htemp);
	}
      }
    }
  }
  TString partname="Both";
  if(optPartAntiPart==kParticleOnly) partname="Both";
  if(optPartAntiPart==kAntiParticleOnly) partname="Both";

  for(Int_t iFile=0; iFile<nReadFiles; iFile++){
    TH2F* htemp2=(TH2F*)hlist[iFile]->FindObject("hPtVsMassPhi");
    if(iFile==0){
      hPtMass=new TH2F(*htemp2);
    }else{
      hPtMass->Add(htemp2);
    }
  }

  TFile* outf=new TFile(Form("RawYield%s.root",partname.Data()),"recreate");
  outf->cd();
  dcuts[0]->Write();
  outf->Close();

  return kTRUE;

}


Bool_t LoadD0toKpiHistos(TObjArray* listFiles, TH1F** hMass){
  //
  Int_t nFiles=listFiles->GetEntries();
  TList **hlist=new TList*[nFiles];
  AliRDHFCutsD0toKpi** dcuts=new AliRDHFCutsD0toKpi*[nFiles];

  Int_t nReadFiles=0;
  for(Int_t iFile=0; iFile<nFiles; iFile++){
    TString fName=((TObjString*)listFiles->At(iFile))->GetString();    
    TFile *f=TFile::Open(fName.Data());
    if(!f){
      printf("ERROR: file %s does not exist\n",fName.Data());
      continue;
    }
    printf("Open File %s\n",f->GetName());

    TString dirname="PWG3_D2H_D0InvMass";
    if(optPartAntiPart==kParticleOnly) dirname+="D0";
    else if(optPartAntiPart==kAntiParticleOnly) dirname+="D0bar";
    if(cutsappliedondistr) dirname+="C";
    TDirectory *dir = (TDirectory*)f->Get(dirname);
    if(!dir){
      printf("ERROR: directory %s not found in %s\n",dirname.Data(),fName.Data());
      continue;
    }
    TString listmassname="coutputmassD0Mass";
    if(optPartAntiPart==kParticleOnly) listmassname+="D0";
    else if(optPartAntiPart==kAntiParticleOnly) listmassname+="D0bar";
    if(cutsappliedondistr) listmassname+="C";

    hlist[nReadFiles]=(TList*)dir->Get(listmassname);

    TString cutsobjname="cutsD0";
    if(optPartAntiPart==kParticleOnly) cutsobjname+="D0";
    else if(optPartAntiPart==kAntiParticleOnly) cutsobjname+="D0bar";
    if(cutsappliedondistr) cutsobjname+="C";

    dcuts[nReadFiles]=(AliRDHFCutsD0toKpi*)dir->Get(cutsobjname);
    if(nReadFiles>0){
      Bool_t sameCuts=dcuts[nReadFiles]->CompareCuts(dcuts[0]);
      if(!sameCuts){
	printf("ERROR: Cut objects do not match\n");
	return kFALSE;
      }
    }
    nReadFiles++;
  }
  if(nReadFiles<nFiles){
    printf("WARNING: not all requested files have been found\n");
    if (nReadFiles==0) {
      printf("ERROR: Any file/dir found\n");
      return kFALSE;
    }
  }

  Int_t nPtBinsCuts=dcuts[0]->GetNPtBins();
  printf("Number of pt bins for cut object = %d\n",nPtBins);
  Float_t *ptlimsCuts=dcuts[0]->GetPtBinLimits();
  ptlimsCuts[nPtBinsCuts]=ptlimsCuts[nPtBinsCuts-1]+4.;

  Int_t iFinBin=0;
  for(Int_t i=0;i<nPtBinsCuts;i++){
    if(ptlimsCuts[i]>=ptlims[iFinBin+1]) iFinBin+=1; 
    if(iFinBin>nPtBins) break;
    if(ptlimsCuts[i]>=ptlims[iFinBin] && 
       ptlimsCuts[i+1]<=ptlims[iFinBin+1]){
      for(Int_t iFile=0; iFile<nReadFiles; iFile++){
	TH1F* htemp=(TH1F*)hlist[iFile]->FindObject(Form("histMass_%d",i));
	if(!hMass[iFinBin]){
	  hMass[iFinBin]=new TH1F(*htemp);
	}else{
	  hMass[iFinBin]->Add(htemp);
	}
      }
    }
  }
  TString partname="Both";
  if(optPartAntiPart==kParticleOnly) partname="D0";
  if(optPartAntiPart==kAntiParticleOnly) partname="D0bar";

  TFile* outf=new TFile(Form("RawYield%s.root",partname.Data()),"recreate");
  outf->cd();
  dcuts[0]->Write();
  outf->Close();
  return kTRUE;
}
/////////////////////Multiplicity
Bool_t LoadDstarD0piHistosMult(TString fileName, TH1F** hMass, Float_t& centralvalue){
  //
  TList *hlist=new TList();
  TFile *f=TFile::Open(fileName.Data());
    if(!f){
      printf("ERROR: file %s does not exist\n",fileName.Data());
      continue;
    }
  printf("Open File %s\n",f->GetName());
  //AliRDHFCutsDStartoKpipi** dcuts=new AliRDHFCutsDStartoKpipi*[nFiles];

    //TString fName=((TObjString*)listFiles->At(iFile))->GetString(); 
    
  //TString listmassname="listMass";
  f->cd();
  
  //hlist=(TList*)f->Get(listmassname);
  //Int_t nPtnPtBins;
  //printf("Number of pt bins for cut object = %d\n",nPtBins);
  //Float_t *ptlimsCuts=dcuts[0]->GetPtBinLimits();
  //ptlimsCuts[nPtBinsCuts]=ptlimsCuts[nPtBinsCuts-1]+4.;




  //Int_t iFinBin=0;
  TH1F* htemp=(TH1F*)f->Get(hmassName.Data());
  TH1F* hmasstemp=(TH1F*)f->Get(hsignalName.Data());
  centralvalue=hmasstemp->GetBinContent(ptbin+1);
  hMass[0]=new TH1F(*htemp);

  //TString partname="Both";
  //if(optPartAntiPart==kParticleOnly) partname="DStar";
  //if(optPartAntiPart==kAntiParticleOnly) partname="DStarbar";

  //TFile* outf=new TFile(Form("RawYield%s.root",partname.Data()),"recreate");
  //outf->cd();
  //dcuts[0]->Write();
  //outf->Close();
  return kTRUE;
}

Bool_t LoadDstarD0piHistosIntMult(TString fileNameInt, TH1F** hMass, Float_t& centralvalue){
  //
  TList *hlist=new TList();
  TFile *f=TFile::Open(fileNameInt.Data());
    if(!f){
      printf("ERROR: file %s does not exist\n",fileNameInt.Data());
      continue;
    }
  printf("Open File %s\n",f->GetName());
  //AliRDHFCutsDStartoKpipi** dcuts=new AliRDHFCutsDStartoKpipi*[nFiles];

    //TString fName=((TObjString*)listFiles->At(iFile))->GetString(); 
    
  //TString listmassname="listMass";
  f->cd();
  
  //hlist=(TList*)f->Get(listmassname);
  //Int_t nPtnPtBins;
  //printf("Number of pt bins for cut object = %d\n",nPtBins);
  //Float_t *ptlimsCuts=dcuts[0]->GetPtBinLimits();
  //ptlimsCuts[nPtBinsCuts]=ptlimsCuts[nPtBinsCuts-1]+4.;




  //Int_t iFinBin=0;
  TH1F* htemp=(TH1F*)f->Get(hmassNameInt.Data()));
  TH1F* hmasstemp=(TH1F*)f->Get(hsignalNameInt.Data());
  centralvalue=hmasstemp->GetBinContent(ptbin+1);
  hMass[0]=new TH1F(*htemp);

  //TString partname="Both";
  //if(optPartAntiPart==kParticleOnly) partname="DStar";
  //if(optPartAntiPart==kAntiParticleOnly) partname="DStarbar";

  //TFile* outf=new TFile(Form("RawYield%s.root",partname.Data()),"recreate");
  //outf->cd();
  //dcuts[0]->Write();
  //outf->Close();
  return kTRUE;
}

Bool_t LoadDstarToyMC(TString fileName, TH1F** hMass){
  //
  TList *hlist=new TList();
  TFile *f=TFile::Open(fileName.Data());
    if(!f){
      printf("ERROR: file %s does not exist\n",fName.Data());
      continue;
    }
  printf("Open File %s\n",f->GetName());
  f->cd();
  TH1F* htemp=(TH1F*)f->Get(Form("histos",0));
  hMass[0]=new TH1F(*htemp);

  return kTRUE;
}

Bool_t LoadDstarD0piHistos(TObjArray* listFiles, TH1F** hMass){
  //
  Int_t nFiles=listFiles->GetEntries();
  TList **hlist=new TList*[nFiles];
  AliRDHFCutsDStartoKpipi** dcuts=new AliRDHFCutsDStartoKpipi*[nFiles];
  Int_t nReadFiles=0;
  for(Int_t iFile=0; iFile<nFiles; iFile++){
    TString fName=((TObjString*)listFiles->At(iFile))->GetString();    
    TFile *f=TFile::Open(fName.Data());
    if(!f){
      printf("ERROR: file %s does not exist\n",fName.Data());
      continue;
    }
    printf("Open File %s\n",f->GetName());
    TString dirname="PWG3_D2H_DStarSpectra";
    TDirectory *dir = (TDirectory*)f->Get(dirname);
    if(!dir){
      printf("ERROR: directory %s not found in %s\n",dirname.Data(),fName.Data());
      continue;
    }
    TString listmassname="DStarPID00";

    hlist[nReadFiles]=(TList*)dir->Get(listmassname);
    TString cutsobjname="DStartoKpipiCuts";
    dcuts[nReadFiles]=(AliRDHFCutsDStartoKpipi*)dir->Get(cutsobjname);
    if(nReadFiles>0){
      Bool_t sameCuts=dcuts[nReadFiles]->CompareCuts(dcuts[0]);
      if(!sameCuts){
	printf("ERROR: Cut objects do not match\n");
	return kFALSE;
      }
    }
    nReadFiles++;
  }
  if(nReadFiles<nFiles){
    printf("WARNING: not all requested files have been found\n");
    if (nReadFiles==0) {
      printf("ERROR: Any file/dir found\n");
      return kFALSE;
    }
  }
  Int_t nPtBinsCuts=dcuts[0]->GetNPtBins();
  printf("Number of pt bins for cut object = %d\n",nPtBins);
  Float_t *ptlimsCuts=dcuts[0]->GetPtBinLimits();
  ptlimsCuts[nPtBinsCuts]=ptlimsCuts[nPtBinsCuts-1]+4.;




  Int_t iFinBin=0;
  for(Int_t i=0;i<nPtBinsCuts;i++){
    if(ptlimsCuts[i]>=ptlims[iFinBin+1]) iFinBin+=1; 
    if(iFinBin>nPtBins) break;
    if(ptlimsCuts[i]>=ptlims[iFinBin] && 
       ptlimsCuts[i+1]<=ptlims[iFinBin+1]){
      for(Int_t iFile=0; iFile<nReadFiles; iFile++){
	TH1F* htemp=(TH1F*)hlist[iFile]->FindObject(Form("histDeltaMass_%d",i));
	if(!hMass[iFinBin]){
	  hMass[iFinBin]=new TH1F(*htemp);
	}else{
	  hMass[iFinBin]->Add(htemp);
	}
      }
    }
  }
  TString partname="Both";
  if(optPartAntiPart==kParticleOnly) partname="DStar";
  if(optPartAntiPart==kAntiParticleOnly) partname="DStarbar";

  TFile* outf=new TFile(Form("RawYield%s.root",partname.Data()),"recreate");
  outf->cd();
  dcuts[0]->Write();
  outf->Close();
  return kTRUE;
}


void CompareFitTypes(TString* paths, TString* legtext,Int_t ncmp=3,TString* filenameYield=0x0){
  //read ncmp RawYield.roots and draw them together
  //set the global variable nevents before running
  //arguments:
  // - paths= vector of ncmp dimension with the paths of the file RawYield.root
  // - legtext= vector of ncmp dimension with the label for the legend
  // - ncmp= number of files to compare (default is 3)
  // - filenameYield= optional ncmp-dimensional array with the difference between the names of the files to be compared (useful if the 2 files are in the same directory but have different names)

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetFrameFillColor(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetFrameBorderMode(0);

  if(!filenameYield) filenameYield=new TString[ncmp];

  for(Int_t k=0;k<ncmp;k++){
    if(!filenameYield) filenameYield[k]="RawYield.root";
    filenameYield[k].Prepend(paths[k]);
  }
  
  TCanvas* cSig=new TCanvas("cSig","Results",1200,600);
  cSig->Divide(3,1);
  TCanvas* cBkgN=new TCanvas("cBkgN","Background normalized to sigma",400,600);
  TCanvas* cDiffS=new TCanvas("cDiffS","Difference",1200,600);
  cDiffS->Divide(2,1);
  TCanvas *cChi2=new TCanvas("cChi2","reduced chi square",600,600);

  TLegend* leg1=new TLegend(0.4,0.7,0.89,0.89);
  leg1->SetFillColor(0);
  TLegend* leg2=(TLegend*)leg1->Clone();
  TLegend* leg3=(TLegend*)leg1->Clone();
  TLegend* leg4=new TLegend(0.4,0.6,0.8,0.89);
  leg4->SetFillColor(0);

  for(Int_t iTypes=0;iTypes<ncmp;iTypes++){
    TFile* fin=new TFile(filenameYield[iTypes]);
    if(!fin){
      printf("WARNING: %s not found",filenameYield[iTypes].Data());
      continue;
    }

    TH1F* hSignal=(TH1F*)fin->Get("hSignal");
    TH1F* hBackground=(TH1F*)fin->Get("hBackground");
    TH1F* hBackgroundNormSigma=(TH1F*)fin->Get("hBackgroundNormSigma");
    TH1F* hSignificance=(TH1F*)fin->Get("hSignificance");
    hSignal->SetName(Form("%s%d",hSignal->GetName(),iTypes));
    hBackground->SetName(Form("%s%d",hBackground->GetName(),iTypes));
    hBackgroundNormSigma->SetName(Form("%s%d",hBackgroundNormSigma->GetName(),iTypes));
    hSignificance->SetName(Form("%s%d",hSignificance->GetName(),iTypes));

    hSignal->SetMarkerColor(iTypes+2);
    hSignal->SetLineColor(iTypes+2);
    hBackground->SetMarkerColor(iTypes+2);
    hBackground->SetLineColor(iTypes+2);
    hBackgroundNormSigma->SetMarkerColor(iTypes+2);
    hBackgroundNormSigma->SetLineColor(iTypes+2);
    hSignificance->SetMarkerColor(iTypes+2);
    hSignificance->SetLineColor(iTypes+2);

    TLegendEntry* ent4=leg4->AddEntry(hSignal,Form("%s",legtext[iTypes].Data()),"PL");
    ent4->SetTextColor(hSignal->GetMarkerColor());

    if(ncmp==nsamples){
      printf("Info: Normalizing signal, background and significance to the number of events\n");
      hSignal->Scale(1./nevents[iTypes]);
      hBackground->Scale(1./nevents[iTypes]);
      hBackgroundNormSigma->Scale(1./nevents[iTypes]);
      hSignificance->Scale(1./TMath::Sqrt(nevents[iTypes]));
    }

    if(iTypes==0){
      cSig->cd(1);
      hSignal->DrawClone("P");
      cSig->cd(2);
      hBackground->DrawClone("P");
      cSig->cd(3);
      hSignificance->DrawClone("P");
      cBkgN->cd();
      hBackgroundNormSigma->DrawClone("P");
    } else{
      cSig->cd(1);
      hSignal->DrawClone("Psames");
      cSig->cd(2);
      hBackground->DrawClone("Psames");
      cSig->cd(3);
      hSignificance->DrawClone("Psames");
      cBkgN->cd();
      hBackgroundNormSigma->DrawClone("Psames");
    }

    TH1F* hRelErrSig=(TH1F*)fin->Get("hRelErrSig");
    TH1F* hInvSignif=(TH1F*)fin->Get("hInvSignif");
    hRelErrSig->SetName(Form("%s%d",hRelErrSig->GetName(),iTypes));
    hInvSignif->SetName(Form("%s%d",hInvSignif->GetName(),iTypes));

    hRelErrSig->SetMarkerColor(iTypes+2);
    hRelErrSig->SetLineColor(iTypes+2);
    hInvSignif->SetMarkerColor(iTypes+2);
    hInvSignif->SetLineColor(iTypes+2);

    TLegendEntry* ent1=leg1->AddEntry(hRelErrSig,Form("From Fit (%s)",legtext[iTypes].Data()),"P");
    ent1->SetTextColor(hRelErrSig->GetMarkerColor());
    ent1=leg1->AddEntry(hInvSignif,Form("1/Significance (%s)",legtext[iTypes].Data()),"PL");
    ent1->SetTextColor(hInvSignif->GetMarkerColor());

    cDiffS->cd(1);
    if(iTypes==0){
      hRelErrSig->DrawClone("P");
      hInvSignif->DrawClone();
    } else{
      hRelErrSig->DrawClone("Psames");
      hInvSignif->DrawClone("sames");
    }

    TH1F* hNDiffCntSig1=(TH1F*)fin->Get("hNDiffCntSig1");
    TH1F* hNDiffCntSig2=(TH1F*)fin->Get("hNDiffCntSig2");
    hNDiffCntSig1->SetName(Form("%s%d",hNDiffCntSig1->GetName(),iTypes));
    hNDiffCntSig2->SetName(Form("%s%d",hNDiffCntSig2->GetName(),iTypes));

    hNDiffCntSig1->SetMarkerColor(iTypes+2);
    hNDiffCntSig1->SetLineColor(iTypes+2);
    hNDiffCntSig2->SetMarkerColor(iTypes+2);
    hNDiffCntSig2->SetLineColor(iTypes+2);
    TLegendEntry* ent2=leg2->AddEntry(hNDiffCntSig1,Form("From Counting1 (%s)",legtext[iTypes].Data()),"PL");
    ent2->SetTextColor(hNDiffCntSig1->GetMarkerColor());
    ent2=leg2->AddEntry(hNDiffCntSig2,Form("From Counting2 (%s)",legtext[iTypes].Data()),"PL");
    ent2->SetTextColor(hNDiffCntSig2->GetMarkerColor());

    cDiffS->cd(2);
    if(iTypes==0){
      hNDiffCntSig1->DrawClone();
      hNDiffCntSig2->DrawClone();
    }else{
     hNDiffCntSig1->DrawClone("sames");
     hNDiffCntSig2->DrawClone("sames");
    }

    TGraph* grReducedChiSquare=(TGraph*)fin->Get("grReducedChiSquare");
    grReducedChiSquare->SetName(Form("%s%d",grReducedChiSquare->GetName(),iTypes));

    grReducedChiSquare->SetMarkerColor(iTypes+2);
    grReducedChiSquare->SetLineColor(iTypes+2);
    TLegendEntry* ent3=leg3->AddEntry(grReducedChiSquare,Form("%s",legtext[iTypes].Data()),"PL");
    ent3->SetTextColor(grReducedChiSquare->GetMarkerColor());

    cChi2->cd();
    if(iTypes==0) grReducedChiSquare->DrawClone("AP");
    else grReducedChiSquare->DrawClone("P");
  }

  cSig->cd(1);
  leg4->Draw();

  cDiffS->cd(1);
  leg1->Draw();

  cDiffS->cd(2);
  leg2->Draw();

  cChi2->cd();
  leg3->Draw();

  TFile* fout=new TFile("ComparisonRawYield.root","RECREATE");
  fout->cd();
  cDiffS->Write();
  cChi2->Write();
  fout->Close();
}



TH1F* RebinHisto(TH1F* hOrig, Int_t reb, Int_t firstUse){
  // Rebin histogram, from bin firstUse to lastUse
  // Use all bins if firstUse=-1

  Int_t nBinOrig=hOrig->GetNbinsX();
  Int_t firstBinOrig=1;
  Int_t lastBinOrig=nBinOrig;
  Int_t nBinOrigUsed=nBinOrig;
  Int_t nBinFinal=nBinOrig/reb;
  if(firstUse>=1){ 
    firstBinOrig=firstUse;
    nBinFinal=(nBinOrig-firstUse+1)/reb;
    nBinOrigUsed=nBinFinal*reb;
    lastBinOrig=firstBinOrig+nBinOrigUsed-1;
  }else{
    Int_t exc=nBinOrigUsed%reb;
    if(exc!=0){
      nBinOrigUsed-=exc;
      firstBinOrig+=exc/2;
      lastBinOrig=firstBinOrig+nBinOrigUsed-1;
    }
  }

  printf("Rebin from %d bins to %d bins -- Used bins=%d in range %d-%d\n",nBinOrig,nBinFinal,nBinOrigUsed,firstBinOrig,lastBinOrig);
  Float_t lowLim=hOrig->GetXaxis()->GetBinLowEdge(firstBinOrig);
  Float_t hiLim=hOrig->GetXaxis()->GetBinUpEdge(lastBinOrig);
  TH1F* hRebin=new TH1F(Form("%s-rebin",hOrig->GetName()),hOrig->GetTitle(),nBinFinal,lowLim,hiLim);
  Int_t lastSummed=firstBinOrig-1;
  for(Int_t iBin=1;iBin<=nBinFinal; iBin++){
    Float_t sum=0.;
    for(Int_t iOrigBin=0;iOrigBin<reb;iOrigBin++){
      sum+=hOrig->GetBinContent(lastSummed+1);
      lastSummed++;
    }
    hRebin->SetBinContent(iBin,sum);
  }
  return hRebin;
}
