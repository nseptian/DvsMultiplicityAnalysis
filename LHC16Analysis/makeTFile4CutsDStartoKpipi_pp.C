// gSystem->SetIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_ROOT/ITS -I$ALICE_ROOT/TPC -I$ALICE_ROOT/CONTAINERS -I$ALICE_ROOT/STEER -I$ALICE_ROOT/TRD -I$ALICE_ROOT/macros -I$ALICE_ROOT/ANALYSIS -I$ALICE_ROOT/PWG3 -I$ALICE_ROOT/PWG3/vertexingHF -g");

#include <Riostream.h>
#include <TFile.h>
#include "AliRDHFCutsDStartoKpipi.h"
#include <TClonesArray.h>
#include <TParameter.h>


//Use:
//Set hard coded commentet with //set this!!
// root[] .L makeTFile4CutsDStar.....C++
// root[] makeInputAliAnalysisTaskSED0Mass()
// root[] makeInputAliAnalysisTaskSESignificanceMaximization()
//similar macros for the other D mesons

//Author: Alessandro Grelli, a.grelli@uu.nl


//macro to make a .root file which contains an AliRDHFCutsDStartoKpipi for AliAnalysisTaskSEDStarSpectra task and CF task
void makeInputAliAnalysisTaskSEDStarSpectra(const char *set_cuts="pp13TeV"){

  AliRDHFCutsDStartoKpipi* RDHFDStartoKpipi=new AliRDHFCutsDStartoKpipi();

RDHFDStartoKpipi->SetUseCentrality(kFALSE);
  //Centrality selection
 //RDHFDStartoKpipi->SetTriggerClass("CVHMSH2");
 RDHFDStartoKpipi->SetTriggerClass("CVHMV0M");
 //RDHFDStartoKpipi->SetTriggerClass("");
 //RDHFDStartoKpipi->SetTriggerMask(0);
 //RDHFDStartoKpipi->SetTriggerMask(AliVEvent::kINT7);
 RDHFDStartoKpipi->SetTriggerMask(AliVEvent::kHighMultV0);
 //RDHFDStartoKpipi->SetTriggerClass("CINT7");
 RDHFDStartoKpipi->SetUseOnlyOneTrigger(kFALSE);
 RDHFDStartoKpipi->SetUsePhysicsSelection(kTRUE);
// RDHFDStartoKpipi->SetTriggerMask(AliVEvent::kEMC1 | AliVvent::kEMC7);
// RDHFDStartoKpipi->SetTriggerClass("CEMC");


  RDHFDStartoKpipi->SetName("DStartoKpipiCuts");
  RDHFDStartoKpipi->SetTitle("Cuts for D* analysis");

  AliESDtrackCuts* esdTrackCuts=new AliESDtrackCuts();
  esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
  //default
  esdTrackCuts->SetRequireTPCRefit(kTRUE);
  esdTrackCuts->SetRequireITSRefit(kTRUE);
  esdTrackCuts->SetEtaRange(-0.8,0.8);
  esdTrackCuts->SetMinNClustersTPC(70);
  esdTrackCuts->SetMaxChi2PerClusterTPC(2);
  esdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
  esdTrackCuts->SetMinNClustersITS(2);
  esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny); //kFirst; kAny I selected myself
  esdTrackCuts->SetMinDCAToVertexXY(0.);
  esdTrackCuts->SetPtRange(0.3,1.e10); //later has to be edited, too large minimum Pt?

 // default is kBoth, otherwise kAny
  //esdTrackCuts->SetMinDCAToVertexXY(0.);
 // esdTrackCuts->SetPtRange(0.3,1.e10);

 // soft pion pre-selections
  AliESDtrackCuts* esdSoftPicuts=new AliESDtrackCuts();
  esdSoftPicuts->SetRequireSigmaToVertex(kFALSE);
  //default
  esdSoftPicuts->SetRequireTPCRefit(kFALSE);//////////////////////////////
  esdSoftPicuts->SetRequireITSRefit(kTRUE);
  //esdSoftPicuts->SetMinNClustersTPC(40);
  //esdSoftPicuts->SetMinNClustersITS(3);
  esdSoftPicuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny); //test d0 asimmetry ////////////////////////////
  esdSoftPicuts->SetPtRange(0.0,1.e10);

  // set pre selections
  RDHFDStartoKpipi->AddTrackCuts(esdTrackCuts);
  RDHFDStartoKpipi->AddTrackCutsSoftPi(esdSoftPicuts);

  const Int_t nvars=16;
  const Int_t nptbins=13;

  Float_t* ptbins;
  ptbins=new Float_t[nptbins+1];
  ptbins[0]=0.5;
  ptbins[1]=1.;
  ptbins[2]=2.;
  ptbins[3]=3.;
  ptbins[4]=4.;
  ptbins[5]=5.;
  ptbins[6]=6.;
  ptbins[7]=7.;
  ptbins[8]=8.;
  ptbins[9]=12.;
  ptbins[10]=16.;
  ptbins[11]=24.;
  ptbins[12]=36.;
  ptbins[13]=9999.;
  //ptbins[14]=64.;
  RDHFDStartoKpipi->SetPtBins(nptbins+1,ptbins);

  Float_t** rdcutsvalmine;
  rdcutsvalmine=new Float_t*[nvars];
  for(Int_t iv=0;iv<nvars;iv++){
    rdcutsvalmine[iv]=new Float_t[nptbins];
  }

  if(set_cuts=="pp13TeV"){
          //040612
     //0.5-1
    rdcutsvalmine[0][0]=0.03;
    rdcutsvalmine[1][0]=0.05;//0.1;//
    rdcutsvalmine[2][0]=0.9;
    rdcutsvalmine[3][0]=0.45;
    rdcutsvalmine[4][0]=0.45;
    rdcutsvalmine[5][0]=0.1;
    rdcutsvalmine[6][0]=0.1;
    rdcutsvalmine[7][0]=-0.000075;;//-0.00005; //-0.0001; //was 0.01
    rdcutsvalmine[8][0]=0.8;
    rdcutsvalmine[9][0]=0.3;
    rdcutsvalmine[10][0]=0.015;//0.3;
    rdcutsvalmine[11][0]=0.03;
    rdcutsvalmine[12][0]=0.3;//0.1; //was 100
    rdcutsvalmine[13][0]=1;
    rdcutsvalmine[14][0]=0.;
    rdcutsvalmine[15][0]=1.5;
    
    //1-2

    rdcutsvalmine[0][1]=0.03;
    rdcutsvalmine[1][1]=0.0315;//0.035;  //0.04
    rdcutsvalmine[2][1]=0.8;
    rdcutsvalmine[3][1]=0.5;
    rdcutsvalmine[4][1]=0.5;
    rdcutsvalmine[5][1]=0.1;  //0.04
    rdcutsvalmine[6][1]=0.1;
    rdcutsvalmine[7][1]=-0.00033;//-0.0008;//-0.0005;//-0.0002;//-0.00015; // -0.0003
    rdcutsvalmine[8][1]=0.865;
    rdcutsvalmine[9][1]=0.3;
    rdcutsvalmine[10][1]=0.15;
    rdcutsvalmine[11][1]=0.05;
    rdcutsvalmine[12][1]=0.5;
    rdcutsvalmine[13][1]=1;
    rdcutsvalmine[14][1]=0.; ///-1
    rdcutsvalmine[15][1]=4.;
    /*
    //040612
     //0.5-1
    rdcutsvalmine[0][0]=0.050;
    rdcutsvalmine[1][0]=0.05;
    rdcutsvalmine[2][0]=0.9;
    rdcutsvalmine[3][0]=0.45;
    rdcutsvalmine[4][0]=0.45;
    rdcutsvalmine[5][0]=0.1;
    rdcutsvalmine[6][0]=0.1;
    rdcutsvalmine[7][0]=-0.000075;//-0.00005; //-0.0001; //was 0.01
    rdcutsvalmine[8][0]=0.8;
    rdcutsvalmine[9][0]=0.3;
    rdcutsvalmine[10][0]=0.3;
    rdcutsvalmine[11][0]=0.03;
    rdcutsvalmine[12][0]=0.1; //was 100
    rdcutsvalmine[13][0]=1;
    rdcutsvalmine[14][0]=0.85;
    rdcutsvalmine[15][0]=1.5;
    //1-2
    rdcutsvalmine[0][1]=0.032;
    rdcutsvalmine[1][1]=0.035;  //0.04
    rdcutsvalmine[2][1]=0.9;
    rdcutsvalmine[3][1]=0.6;
    rdcutsvalmine[4][1]=0.6;
    rdcutsvalmine[5][1]=0.1;  //0.04
    rdcutsvalmine[6][1]=0.1;
    rdcutsvalmine[7][1]=-0.00025;//-0.0005;//-0.0002;//-0.00015; // -0.0003
    rdcutsvalmine[8][1]=0.82;
    rdcutsvalmine[9][1]=0.3;
    rdcutsvalmine[10][1]=0.3;
    rdcutsvalmine[11][1]=0.05;
    rdcutsvalmine[12][1]=0.2;
    rdcutsvalmine[13][1]=1;
    rdcutsvalmine[14][1]=0.88; ///-1
    rdcutsvalmine[15][1]=3;
    */
//2-3
    rdcutsvalmine[0][2]=0.032;
    rdcutsvalmine[1][2]=0.027;
    rdcutsvalmine[2][2]=0.8;
    rdcutsvalmine[3][2]=1.0;
    rdcutsvalmine[4][2]=1.0;
    rdcutsvalmine[5][2]=0.1;
    rdcutsvalmine[6][2]=0.1;
    rdcutsvalmine[7][2]=-0.00010;//-0.00019;
    rdcutsvalmine[8][2]=0.9;
    rdcutsvalmine[9][2]=0.3;
    rdcutsvalmine[10][2]=0.15;
    rdcutsvalmine[11][2]=0.05;
    rdcutsvalmine[12][2]=0.5;
    rdcutsvalmine[13][2]=1;
    rdcutsvalmine[14][2]=0.;
    rdcutsvalmine[15][2]=4.;

    //3-4
    rdcutsvalmine[0][3]=0.032;
    rdcutsvalmine[1][3]=0.03375;
    rdcutsvalmine[2][3]=0.8;
    rdcutsvalmine[3][3]=1.0;
    rdcutsvalmine[4][3]=1.0;
    rdcutsvalmine[5][3]=0.1;
    rdcutsvalmine[6][3]=0.1;
    rdcutsvalmine[7][3]=-0.00010;
    rdcutsvalmine[8][3]=0.83;
    rdcutsvalmine[9][3]=0.3;
    rdcutsvalmine[10][3]=0.15;
    rdcutsvalmine[11][3]=0.05;
    rdcutsvalmine[12][3]=0.5;
    rdcutsvalmine[13][3]=1;
    rdcutsvalmine[14][3]=0.;
    rdcutsvalmine[15][3]=0.;

    //4-5
    rdcutsvalmine[0][4]=0.032;
    rdcutsvalmine[1][4]=0.042;
    rdcutsvalmine[2][4]=0.9;
    rdcutsvalmine[3][4]=1.0;
    rdcutsvalmine[4][4]=1.0;
    rdcutsvalmine[5][4]=0.1;
    rdcutsvalmine[6][4]=0.1;
    rdcutsvalmine[7][4]=-0.000028;//0.00005; //-12
    rdcutsvalmine[8][4]=0.81; //0.88
    rdcutsvalmine[9][4]=0.3;
    rdcutsvalmine[10][4]=0.1;
    rdcutsvalmine[11][4]=0.05; //0.05
    rdcutsvalmine[12][4]=100;
    rdcutsvalmine[13][4]=1;
    rdcutsvalmine[14][4]=0.;
    rdcutsvalmine[15][4]=0.;

//5-6
    rdcutsvalmine[0][5]=0.036;
    rdcutsvalmine[1][5]=0.05;
    rdcutsvalmine[2][5]=1.0;
    rdcutsvalmine[3][5]=1.0;
    rdcutsvalmine[4][5]=1.0;
    rdcutsvalmine[5][5]=0.09;
    rdcutsvalmine[6][5]=0.09;
    rdcutsvalmine[7][5]=0.000055; //-1
    rdcutsvalmine[8][5]=0.79;
    rdcutsvalmine[9][5]=0.3;
    rdcutsvalmine[10][5]=0.1;
    rdcutsvalmine[11][5]=0.05;
    rdcutsvalmine[12][5]=100;
    rdcutsvalmine[13][5]=1;
    rdcutsvalmine[14][5]=0.;
    rdcutsvalmine[15][5]=0.;

//6-7
    rdcutsvalmine[0][6]=0.036;
    rdcutsvalmine[1][6]=0.1;
    rdcutsvalmine[2][6]=1.0;
    rdcutsvalmine[3][6]=1;  //0.8
    rdcutsvalmine[4][6]=1;
    rdcutsvalmine[5][6]=0.1;
    rdcutsvalmine[6][6]=0.1;
    rdcutsvalmine[7][6]=0.001;
    rdcutsvalmine[8][6]=0.7;//0.97;
    rdcutsvalmine[9][6]=0.3;
    rdcutsvalmine[10][6]=0.1;
    rdcutsvalmine[11][6]=0.05;
    rdcutsvalmine[12][6]=100;
    rdcutsvalmine[13][6]=1.;
    rdcutsvalmine[14][6]=0.;
    rdcutsvalmine[15][6]=0.;

//7-8
    rdcutsvalmine[0][7]=0.036;
    rdcutsvalmine[1][7]=0.1;
    rdcutsvalmine[2][7]=1.0;
    rdcutsvalmine[3][7]=1.0;//0.6;
    rdcutsvalmine[4][7]=1.0;
    rdcutsvalmine[5][7]=0.1;
    rdcutsvalmine[6][7]=0.1;
    rdcutsvalmine[7][7]=0.00025;
    rdcutsvalmine[8][7]=0.8;//0.97;
    rdcutsvalmine[9][7]=0.3;
    rdcutsvalmine[10][7]=0.1;
    rdcutsvalmine[11][7]=0.05;
    rdcutsvalmine[12][7]=100;
    rdcutsvalmine[13][7]=1.0;
    rdcutsvalmine[14][7]=0.;
    rdcutsvalmine[15][7]=0.;

    //8-12
    rdcutsvalmine[0][8]=0.05;
    rdcutsvalmine[1][8]=0.105;
    rdcutsvalmine[2][8]=1.0;
    rdcutsvalmine[3][8]=1.0;
    rdcutsvalmine[4][8]=1.0;
    rdcutsvalmine[5][8]=0.1;
    rdcutsvalmine[6][8]=0.1;
    rdcutsvalmine[7][8]=0.01;
    rdcutsvalmine[8][8]=0.68;//0.96;
    rdcutsvalmine[9][8]=0.3;
    rdcutsvalmine[10][8]=0.1;
    rdcutsvalmine[11][8]=0.05;
    rdcutsvalmine[12][8]=100;
    rdcutsvalmine[13][8]=1.0;
    rdcutsvalmine[14][8]=0.;
    rdcutsvalmine[15][8]=0.;//4.;

//12-16     230312:

    rdcutsvalmine[0][9]=0.094;
    rdcutsvalmine[1][9]=0.1;
    rdcutsvalmine[2][9]=1.0;
    rdcutsvalmine[3][9]=0.3;
    rdcutsvalmine[4][9]=0.3;
    rdcutsvalmine[5][9]=0.2;
    rdcutsvalmine[6][9]=0.2;
    rdcutsvalmine[7][9]=0.01;
    rdcutsvalmine[8][9]=0.6;//0.96;
    rdcutsvalmine[9][9]=0.3;
    rdcutsvalmine[10][9]=0.15;
    rdcutsvalmine[11][9]=0.45;
    rdcutsvalmine[12][9]=100;
    rdcutsvalmine[13][9]=1.0;
    rdcutsvalmine[14][9]=0.;
    rdcutsvalmine[15][9]=0.;//4.;

//12-16 (ptbin11)
    rdcutsvalmine[0][10]=0.094;
    rdcutsvalmine[1][10]=0.1;//0.035;
    rdcutsvalmine[2][10]=1.0;
    rdcutsvalmine[3][10]=0.3;
    rdcutsvalmine[4][10]=0.3;
    rdcutsvalmine[5][10]=0.2;
    rdcutsvalmine[6][10]=0.2;
    rdcutsvalmine[7][10]=0.01;
    rdcutsvalmine[8][10]=0.6;
    rdcutsvalmine[9][10]=0.3;
    rdcutsvalmine[10][10]=0.15;
    rdcutsvalmine[11][10]=0.45;
    rdcutsvalmine[12][10]=100.;
    rdcutsvalmine[13][10]=1;
    rdcutsvalmine[14][10]=0.;
    rdcutsvalmine[15][10]=0.;//3;

    //16-24 (ptbin12)
    rdcutsvalmine[0][11]=0.7;
    rdcutsvalmine[1][11]=0.15;
    rdcutsvalmine[2][11]=1.0;
    rdcutsvalmine[3][11]=0.3;
    rdcutsvalmine[4][11]=0.3;
    rdcutsvalmine[5][11]=0.5;
    rdcutsvalmine[6][11]=0.5;
    rdcutsvalmine[7][11]=0.01; //was 0.04;
    rdcutsvalmine[8][11]=0.6;
    rdcutsvalmine[9][11]=0.3;
    rdcutsvalmine[10][11]=0.3;
    rdcutsvalmine[11][11]=0.05;
    rdcutsvalmine[12][11]=100.;
    rdcutsvalmine[13][11]=1;
    rdcutsvalmine[14][11]=0.;
    rdcutsvalmine[15][11]=0.;//2.0;

//24-36 (ptbin13)
    rdcutsvalmine[0][12]=0.7;
    rdcutsvalmine[1][12]=0.125;
    rdcutsvalmine[2][12]=1.0;
    rdcutsvalmine[3][12]=0.3;
    rdcutsvalmine[4][12]=0.3;
    rdcutsvalmine[5][12]=0.5;
    rdcutsvalmine[6][12]=0.5;
    rdcutsvalmine[7][12]=0.01;//0.03;
    rdcutsvalmine[8][12]=0.6;//0.93;
    rdcutsvalmine[9][12]=0.3;
    rdcutsvalmine[10][12]=0.1;
    rdcutsvalmine[11][12]=0.05;
    rdcutsvalmine[12][12]=10000.;
    rdcutsvalmine[13][12]=0.5;
    rdcutsvalmine[14][12]=-1.;
    rdcutsvalmine[15][12]=0;

//>36 (ptbin14)
    rdcutsvalmine[0][13]=0.7;
    rdcutsvalmine[1][13]=0.15;
    rdcutsvalmine[2][13]=1.0;
    rdcutsvalmine[3][13]=0.3;
    rdcutsvalmine[4][13]=0.3;
    rdcutsvalmine[5][13]=0.5;
    rdcutsvalmine[6][13]=0.5;
    rdcutsvalmine[7][13]=0.01;//0.35;
    rdcutsvalmine[8][13]=0.6;//0.95;
    rdcutsvalmine[9][13]=0.3;
    rdcutsvalmine[10][13]=0.1;
    rdcutsvalmine[11][13]=0.05;
    rdcutsvalmine[12][13]=100000.;
    rdcutsvalmine[13][13]=0.5;
    rdcutsvalmine[14][13]=-1.;
    rdcutsvalmine[15][13]=0;

 //************************************************************************************

/*
       //0-0.5
    rdcutsvalmine[0][0]=0.032;
    rdcutsvalmine[1][0]=0.04;
    rdcutsvalmine[2][0]=0.8;
    rdcutsvalmine[3][0]=1.0;//1.333;
    rdcutsvalmine[4][0]=1.0;
    rdcutsvalmine[5][0]=0.1;
    rdcutsvalmine[6][0]=0.1;
    rdcutsvalmine[7][0]=-0.0002;
    rdcutsvalmine[8][0]=0.99;
    rdcutsvalmine[9][0]=0.3;
    rdcutsvalmine[10][0]=0.15;
    rdcutsvalmine[11][0]=0.05;
    rdcutsvalmine[12][0]=0.5;
    rdcutsvalmine[13][0]=0.5;
    rdcutsvalmine[14][0]=0.998;
    rdcutsvalmine[15][0]=6.0;

     //0.5-1
    rdcutsvalmine[0][1]=0.032;
    rdcutsvalmine[1][1]=0.04;
    rdcutsvalmine[2][1]=0.8;
    rdcutsvalmine[3][1]=1.0;//1.333;
    rdcutsvalmine[4][1]=1.0;
    rdcutsvalmine[5][1]=0.1;
    rdcutsvalmine[6][1]=0.1;
    rdcutsvalmine[7][1]=-0.0002;
    rdcutsvalmine[8][1]=0.99;
    rdcutsvalmine[9][1]=0.3;
    rdcutsvalmine[10][1]=0.15;
    rdcutsvalmine[11][1]=0.05;
    rdcutsvalmine[12][1]=0.5;
    rdcutsvalmine[13][1]=0.5;
    rdcutsvalmine[14][1]=0.998;
    rdcutsvalmine[15][1]=6.0;

    //1-2
    rdcutsvalmine[0][2]=0.032;
    rdcutsvalmine[1][2]=0.0133;
    rdcutsvalmine[2][2]=0.8;
    rdcutsvalmine[3][2]=1.0;
    rdcutsvalmine[4][2]=1.0;
    rdcutsvalmine[5][2]=0.1;
    rdcutsvalmine[6][2]=0.1;
    rdcutsvalmine[7][2]=-0.00005;
    rdcutsvalmine[8][2]=0.9973;
    rdcutsvalmine[9][2]=0.3;
    rdcutsvalmine[10][2]=0.15;
    rdcutsvalmine[11][2]=0.05;
    rdcutsvalmine[12][2]=0.5;
    rdcutsvalmine[13][2]=0.5;
    rdcutsvalmine[14][2]=0.998;
    rdcutsvalmine[15][2]=8.0;
    //2-3
    rdcutsvalmine[0][3]=0.032;
    rdcutsvalmine[1][3]=0.05;
    rdcutsvalmine[2][3]=0.8;
    rdcutsvalmine[3][3]=1.0;//1.333;
    rdcutsvalmine[4][3]=1.0;
    rdcutsvalmine[5][3]=0.1;
    rdcutsvalmine[6][3]=0.1;
    rdcutsvalmine[7][3]=-0.0001;
    rdcutsvalmine[8][3]=0.988;
    rdcutsvalmine[9][3]=0.3;
    rdcutsvalmine[10][3]=0.15;
    rdcutsvalmine[11][3]=0.05;
    rdcutsvalmine[12][3]=0.5;
    rdcutsvalmine[13][3]=0.5;
    rdcutsvalmine[14][3]=0.998;
    rdcutsvalmine[15][3]=7.0;

    //3-4
    rdcutsvalmine[0][4]=0.7;
    rdcutsvalmine[1][4]=0.0225;
    rdcutsvalmine[2][4]=0.7;
    rdcutsvalmine[3][4]=0.9;
    rdcutsvalmine[4][4]=0.9;
    rdcutsvalmine[5][4]=0.1;
    rdcutsvalmine[6][4]=0.1;
    rdcutsvalmine[7][4]=-0.000433;
    rdcutsvalmine[8][4]=0.9973;
    rdcutsvalmine[9][4]=0.3;
    rdcutsvalmine[10][4]=0.1;
    rdcutsvalmine[11][4]=0.05;
    rdcutsvalmine[12][4]=100.;
    rdcutsvalmine[13][4]=8.17;

    //4-5
    rdcutsvalmine[0][5]=0.032;
    rdcutsvalmine[1][5]=0.025;
    rdcutsvalmine[2][5]=0.8;
    rdcutsvalmine[3][5]=1.;
    rdcutsvalmine[4][5]=1.;
    rdcutsvalmine[5][5]=0.1;
    rdcutsvalmine[6][5]=0.1;
    rdcutsvalmine[7][5]=-0.00025;
    rdcutsvalmine[8][5]=0.994;
    rdcutsvalmine[9][5]=0.20;
    rdcutsvalmine[10][5]=0.15;
    rdcutsvalmine[11][5]=0.05;
    rdcutsvalmine[12][5]=10;
    rdcutsvalmine[13][5]=0.5;
    rdcutsvalmine[14][5]=0.998;
    rdcutsvalmine[15][5]=6.;

    //5-6
    rdcutsvalmine[0][6]=0.04;
    rdcutsvalmine[1][6]=0.0125;
    rdcutsvalmine[2][6]=1.;
    rdcutsvalmine[3][6]=1.;
    rdcutsvalmine[4][6]=1.;
    rdcutsvalmine[5][6]=0.075;
    rdcutsvalmine[6][6]=0.075;
    rdcutsvalmine[7][6]=-0.0002;
    rdcutsvalmine[8][6]=0.989;
    rdcutsvalmine[9][6]=0.15;
    rdcutsvalmine[10][6]=0.15;
    rdcutsvalmine[11][6]=0.2;
    rdcutsvalmine[12][6]=10.;
    rdcutsvalmine[13][6]=0.5;
    rdcutsvalmine[14][6]=0.998;
    rdcutsvalmine[15][6]=8.;

    //6-7
    rdcutsvalmine[0][7]=0.036;
    rdcutsvalmine[1][7]=0.02;
    rdcutsvalmine[2][7]=1.0;
    rdcutsvalmine[3][7]=1.;
    rdcutsvalmine[4][7]=1.;
    rdcutsvalmine[5][7]=0.1;
    rdcutsvalmine[6][7]=0.1;
    rdcutsvalmine[7][7]=-0.00013;
    rdcutsvalmine[8][7]=0.99;
    rdcutsvalmine[9][7]=0.3;
    rdcutsvalmine[10][7]=0.1;
    rdcutsvalmine[11][7]=0.2;
    rdcutsvalmine[12][7]=100.;
    rdcutsvalmine[13][7]=0.5;
    rdcutsvalmine[14][7]=0.998;
    rdcutsvalmine[15][7]=5.75;

    //7-8
    rdcutsvalmine[0][8]=0.036;
    rdcutsvalmine[1][8]=0.03;
    rdcutsvalmine[2][8]=1.0;
    rdcutsvalmine[3][8]=1.;
    rdcutsvalmine[4][8]=1.;
    rdcutsvalmine[5][8]=0.12;
    rdcutsvalmine[6][8]=0.12;
    rdcutsvalmine[7][8]=-0.000163;
    rdcutsvalmine[8][8]=0.989;
    rdcutsvalmine[9][8]=0.3;
    rdcutsvalmine[10][8]=0.15;
    rdcutsvalmine[11][8]=0.2;
    rdcutsvalmine[12][8]=100.;
    rdcutsvalmine[13][8]=0.5;
    rdcutsvalmine[14][8]=0.998;
    rdcutsvalmine[15][8]=1.;

    //8-10
    rdcutsvalmine[0][9]=0.055;
    rdcutsvalmine[1][9]=0.0125;
    rdcutsvalmine[2][9]=1.0;
    rdcutsvalmine[3][9]=0.9;
    rdcutsvalmine[4][9]=0.9;
    rdcutsvalmine[5][9]=0.15;
    rdcutsvalmine[6][9]=0.15;
    rdcutsvalmine[7][9]=-0.000067;
    rdcutsvalmine[8][9]=0.9967;
    rdcutsvalmine[9][9]=0.3;
    rdcutsvalmine[10][9]=0.15;
    rdcutsvalmine[11][9]=0.34;
    rdcutsvalmine[12][9]=100.;
    rdcutsvalmine[13][9]=0.5;
    rdcutsvalmine[14][9]=0.998;
    rdcutsvalmine[15][9]=7.33;

    //10-12     230312:
    rdcutsvalmine[0][10]=0.055;
    rdcutsvalmine[1][10]=0.04;
    rdcutsvalmine[2][10]=1.0;
    rdcutsvalmine[3][10]=0.9;
    rdcutsvalmine[4][10]=0.9;
    rdcutsvalmine[5][10]=0.15;
    rdcutsvalmine[6][10]=0.15;
    rdcutsvalmine[7][10]=0.0001;
    rdcutsvalmine[8][10]=0.975;
    rdcutsvalmine[9][10]=0.3;
    rdcutsvalmine[10][10]=0.15;
    rdcutsvalmine[11][10]=0.35;
    rdcutsvalmine[12][10]=100.;
    rdcutsvalmine[13][10]=0.5;
    rdcutsvalmine[14][10]=0.998;
    rdcutsvalmine[15][10]=1.;


    //12-16 (ptbin11)
    rdcutsvalmine[0][11]=0.074;
    rdcutsvalmine[1][11]=0.03;
    rdcutsvalmine[2][11]=0.8;
    rdcutsvalmine[3][11]=0.3;
    rdcutsvalmine[4][11]=0.3;
    rdcutsvalmine[5][11]=0.2;
    rdcutsvalmine[6][11]=0.2;
    rdcutsvalmine[7][11]=0.0001;
    rdcutsvalmine[8][11]=0.965;
    rdcutsvalmine[9][11]=0.3;
    rdcutsvalmine[10][11]=0.1;
    rdcutsvalmine[11][11]=0.45;
    rdcutsvalmine[12][11]=100.;
    rdcutsvalmine[13][11]=1;
    rdcutsvalmine[14][11]=0.998;
    rdcutsvalmine[15][11]=4.25;



    //16-20 (ptbin12)
    rdcutsvalmine[0][12]=0.074;
    rdcutsvalmine[1][12]=0.02;
    rdcutsvalmine[2][12]=1.0;
    rdcutsvalmine[3][12]=0.30;
    rdcutsvalmine[4][12]=0.30;
    rdcutsvalmine[5][12]=0.50;
    rdcutsvalmine[6][12]=0.50;
    rdcutsvalmine[7][12]=0.004;
    rdcutsvalmine[8][12]=0.9583;
    rdcutsvalmine[9][12]=0.1;
    rdcutsvalmine[10][12]=0.1;
    rdcutsvalmine[11][12]=0.05;
    rdcutsvalmine[12][12]=100.;
    rdcutsvalmine[13][12]=1;
    rdcutsvalmine[14][12]=0.9;
    rdcutsvalmine[15][12]=4.0;

    //20-24 (ptbin13)
    rdcutsvalmine[0][13]=0.074;
    rdcutsvalmine[1][13]=0.0433;
    rdcutsvalmine[2][13]=1.0;
    rdcutsvalmine[3][13]=0.30;
    rdcutsvalmine[4][13]=0.30;
    rdcutsvalmine[5][13]=0.50;
    rdcutsvalmine[6][13]=0.50;
    rdcutsvalmine[7][13]=0.004;
    rdcutsvalmine[8][13]=0.7;
    rdcutsvalmine[9][13]=0.1;
    rdcutsvalmine[10][13]=0.1;
    rdcutsvalmine[11][13]=0.05;
    rdcutsvalmine[12][13]=100.;
    rdcutsvalmine[13][13]=1;
    rdcutsvalmine[14][13]=0.9;
    rdcutsvalmine[15][13]=5.33;

    //>24 (ptbin14)
    rdcutsvalmine[0][14]=0.074;
    rdcutsvalmine[1][14]=0.8;
    rdcutsvalmine[2][14]=1.0;
    rdcutsvalmine[3][14]=0.30;
    rdcutsvalmine[4][14]=0.30;
    rdcutsvalmine[5][14]=0.50;
    rdcutsvalmine[6][14]=0.50;
    rdcutsvalmine[7][14]=0.0033;
    rdcutsvalmine[8][14]=0.9417;
    rdcutsvalmine[9][14]=0.1;
    rdcutsvalmine[10][14]=0.1;
    rdcutsvalmine[11][14]=0.05;
    rdcutsvalmine[12][14]=100.;
    rdcutsvalmine[13][14]=1;
    rdcutsvalmine[14][14]=0.9;
    rdcutsvalmine[15][14]=0.67;
*/

//vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
/*
// Cuts 07C_2C (as 07C_2B, then compared with Alessandro Cuts )
    //0-0.5
    rdcutsvalmine[0][0]=0.032; 		// inv mass
    rdcutsvalmine[1][0]=0.02;		// dca
    rdcutsvalmine[2][0]=0.8;		// cos theta star
    rdcutsvalmine[3][0]=1.333;		// ptK
    rdcutsvalmine[4][0]=1.0;		// pt PI
    rdcutsvalmine[5][0]=0.07;		//  d0 K
    rdcutsvalmine[6][0]=0.07;		// d0 Pi
    rdcutsvalmine[7][0]=-0.000362;	// d0d0
    rdcutsvalmine[8][0]=0.99;		// costheta point
    rdcutsvalmine[9][0]=0.3;            // inv mass half width of D*
    rdcutsvalmine[10][0]=0.15;		// half width of MKpiPi -Md0
    rdcutsvalmine[11][0]=0.05;		// pt min of soft pion
    rdcutsvalmine[12][0]=0.5;		// pt max of soft pion
    rdcutsvalmine[13][0]=0.5;		// theta between soft pion and d0 plane
    rdcutsvalmine[14][0]=0.998;		// cos theta XY
    rdcutsvalmine[15][0]=7.0;		// norm decay length

     //0.5-1
    rdcutsvalmine[0][1]=0.032;
    rdcutsvalmine[1][1]=0.02;
    rdcutsvalmine[2][1]=0.8;
    rdcutsvalmine[3][1]=1.333;
    rdcutsvalmine[4][1]=1.0;
    rdcutsvalmine[5][1]=0.1;		//0.07;
    rdcutsvalmine[6][1]=0.1;		//0.07;
    rdcutsvalmine[7][1]=-0.000362;
    rdcutsvalmine[8][1]=0.99;
    rdcutsvalmine[9][1]=0.3;
    rdcutsvalmine[10][1]=0.15;
    rdcutsvalmine[11][1]=0.05;
    rdcutsvalmine[12][1]=0.5;
    rdcutsvalmine[13][1]=0.5;
    rdcutsvalmine[14][1]=0.998;
    rdcutsvalmine[15][1]=7.0;

    //1-2
    rdcutsvalmine[0][2]=0.032;
    rdcutsvalmine[1][2]=0.025;
    rdcutsvalmine[2][2]=0.8;
    rdcutsvalmine[3][2]=1.0;
    rdcutsvalmine[4][2]=1.0;
    rdcutsvalmine[5][2]=0.1;
    rdcutsvalmine[6][2]=0.1;
    rdcutsvalmine[7][2]=-0.0002;
    rdcutsvalmine[8][2]=0.995;
    rdcutsvalmine[9][2]=0.3;
    rdcutsvalmine[10][2]=0.15;
    rdcutsvalmine[11][2]=0.05;
    rdcutsvalmine[12][2]=0.5;
    rdcutsvalmine[13][2]=0.5;
    rdcutsvalmine[14][2]=0.998;
    rdcutsvalmine[15][2]=9;
    //2-3
    rdcutsvalmine[0][3]=0.032;
    rdcutsvalmine[1][3]=0.025;
    rdcutsvalmine[2][3]=0.8;
    rdcutsvalmine[3][3]=1.0;
    rdcutsvalmine[4][3]=1.0;
    rdcutsvalmine[5][3]=0.1;
    rdcutsvalmine[6][3]=0.1;
    rdcutsvalmine[7][3]=-0.0003;
    rdcutsvalmine[8][3]=0.996;
    rdcutsvalmine[9][3]=0.3;
    rdcutsvalmine[10][3]=0.15;
    rdcutsvalmine[11][3]=0.05;
    rdcutsvalmine[12][3]=0.5;
    rdcutsvalmine[13][3]=0.5;
    rdcutsvalmine[14][3]=0.998;
    rdcutsvalmine[15][3]=8.0;

    //3-4
    rdcutsvalmine[0][4]=0.032;
    rdcutsvalmine[1][4]=0.025;
    rdcutsvalmine[2][4]=0.8;
    rdcutsvalmine[3][4]=1.0;
    rdcutsvalmine[4][4]=1.0;
    rdcutsvalmine[5][4]=0.1;
    rdcutsvalmine[6][4]=0.1;
    rdcutsvalmine[7][4]=-0.000362;
    rdcutsvalmine[8][4]=0.996;
    rdcutsvalmine[9][4]=0.3;
    rdcutsvalmine[10][4]=0.15;
    rdcutsvalmine[11][4]=0.05;
    rdcutsvalmine[12][4]=0.5;
    rdcutsvalmine[13][4]=0.5;
    rdcutsvalmine[14][4]=0.998;
    rdcutsvalmine[15][4]=7.0;

    //4-5
    rdcutsvalmine[0][5]=0.032;
    rdcutsvalmine[1][5]=0.02;
    rdcutsvalmine[2][5]=0.8;
    rdcutsvalmine[3][5]=1.;
    rdcutsvalmine[4][5]=1.;
    rdcutsvalmine[5][5]=0.1;
    rdcutsvalmine[6][5]=0.1;
    rdcutsvalmine[7][5]=-0.00028;
    rdcutsvalmine[8][5]=0.996;
    rdcutsvalmine[9][5]=0.2;
    rdcutsvalmine[10][5]=0.15;
    rdcutsvalmine[11][5]=0.05;
    rdcutsvalmine[12][5]=10;
    rdcutsvalmine[13][5]=0.5;
    rdcutsvalmine[14][5]=0.998;
    rdcutsvalmine[15][5]=6.5;

    //5-6
    rdcutsvalmine[0][6]=0.04;
    rdcutsvalmine[1][6]=0.02;
    rdcutsvalmine[2][6]=1.;
    rdcutsvalmine[3][6]=1.;
    rdcutsvalmine[4][6]=1.;
    rdcutsvalmine[5][6]=0.1;
    rdcutsvalmine[6][6]=0.1;
    rdcutsvalmine[7][6]=-0.00025;
    rdcutsvalmine[8][6]=0.994;
    rdcutsvalmine[9][6]=0.15;
    rdcutsvalmine[10][6]=0.15;
    rdcutsvalmine[11][6]=0.2;
    rdcutsvalmine[12][6]=10.;
    rdcutsvalmine[13][6]=0.5;
    rdcutsvalmine[14][6]=0.998;
    rdcutsvalmine[15][6]=6.0;

    //6-7
    rdcutsvalmine[0][7]=0.036;
    rdcutsvalmine[1][7]=0.02;
    rdcutsvalmine[2][7]=1.0;
    rdcutsvalmine[3][7]=1.;
    rdcutsvalmine[4][7]=1.;
    rdcutsvalmine[5][7]=0.1;
    rdcutsvalmine[6][7]=0.1;
    rdcutsvalmine[7][7]=-0.0002;
    rdcutsvalmine[8][7]=0.994;
    rdcutsvalmine[9][7]=0.3;
    rdcutsvalmine[10][7]=0.1;
    rdcutsvalmine[11][7]=0.2;
    rdcutsvalmine[12][7]=100.;
    rdcutsvalmine[13][7]=0.5;
    rdcutsvalmine[14][7]=0.998;
    rdcutsvalmine[15][7]=5.5;

    //7-8
    rdcutsvalmine[0][8]=0.036;
    rdcutsvalmine[1][8]=0.02;
    rdcutsvalmine[2][8]=1.0;
    rdcutsvalmine[3][8]=1.;
    rdcutsvalmine[4][8]=1.;
    rdcutsvalmine[5][8]=0.12;
    rdcutsvalmine[6][8]=0.12;
    rdcutsvalmine[7][8]=-0.000165;
    rdcutsvalmine[8][8]=0.992;
    rdcutsvalmine[9][8]=0.3;
    rdcutsvalmine[10][8]=0.15;
    rdcutsvalmine[11][8]=0.2;
    rdcutsvalmine[12][8]=100.;
    rdcutsvalmine[13][8]=0.5;
    rdcutsvalmine[14][8]=0.998;
    rdcutsvalmine[15][8]=5.;

    //8-10
    rdcutsvalmine[0][9]=0.055;
    rdcutsvalmine[1][9]=0.02;
    rdcutsvalmine[2][9]=1.0;
    rdcutsvalmine[3][9]=0.9;
    rdcutsvalmine[4][9]=0.9;
    rdcutsvalmine[5][9]=0.15;
    rdcutsvalmine[6][9]=0.15;
    rdcutsvalmine[7][9]=-0.000065;
    rdcutsvalmine[8][9]=0.992;
    rdcutsvalmine[9][9]=0.3;
    rdcutsvalmine[10][9]=0.15;
    rdcutsvalmine[11][9]=0.34;
    rdcutsvalmine[12][9]=100.;
    rdcutsvalmine[13][9]=0.5;
    rdcutsvalmine[14][9]=0.998;
    rdcutsvalmine[15][9]=4.5;

    //10-12
    rdcutsvalmine[0][10]=0.055;
    rdcutsvalmine[1][10]=0.02;
    rdcutsvalmine[2][10]=1.0;
    rdcutsvalmine[3][10]=0.9;
    rdcutsvalmine[4][10]=0.9;
    rdcutsvalmine[5][10]=0.15;
    rdcutsvalmine[6][10]=0.15;
    rdcutsvalmine[7][10]=0.0;
    rdcutsvalmine[8][10]=0.98;
    rdcutsvalmine[9][10]=0.3;
    rdcutsvalmine[10][10]=0.15;
    rdcutsvalmine[11][10]=0.35;
    rdcutsvalmine[12][10]=100.;
    rdcutsvalmine[13][10]=0.5;
    rdcutsvalmine[14][10]=0.998;
    rdcutsvalmine[15][10]=4.;


    //12-16 (ptbin11)
    rdcutsvalmine[0][11]=0.074;
    rdcutsvalmine[1][11]=0.025;
    rdcutsvalmine[2][11]=0.8;
    rdcutsvalmine[3][11]=0.3;
    rdcutsvalmine[4][11]=0.3;
    rdcutsvalmine[5][11]=0.2;
    rdcutsvalmine[6][11]=0.2;
    rdcutsvalmine[7][11]=0.0000;
    rdcutsvalmine[8][11]=0.97;
    rdcutsvalmine[9][11]=0.3;
    rdcutsvalmine[10][11]=0.1;
    rdcutsvalmine[11][11]=0.45;
    rdcutsvalmine[12][11]=100.;
    rdcutsvalmine[13][11]=1;
    rdcutsvalmine[14][11]=0.998;
    rdcutsvalmine[15][11]=4.0;



    //16-20 (ptbin12)
    rdcutsvalmine[0][12]=0.074;
    rdcutsvalmine[1][12]=0.035;
    rdcutsvalmine[2][12]=1.0;
    rdcutsvalmine[3][12]=0.30;
    rdcutsvalmine[4][12]=0.30;
    rdcutsvalmine[5][12]=0.50;
    rdcutsvalmine[6][12]=0.50;
    rdcutsvalmine[7][12]=0.0033;
    rdcutsvalmine[8][12]=0.95;
    rdcutsvalmine[9][12]=0.1;
    rdcutsvalmine[10][12]=0.1;
    rdcutsvalmine[11][12]=0.05;
    rdcutsvalmine[12][12]=100.;
    rdcutsvalmine[13][12]=1;
    rdcutsvalmine[14][12]=0.99;
    rdcutsvalmine[15][12]=3.5;

    //20-24 (ptbin13)
    rdcutsvalmine[0][13]=0.074;
    rdcutsvalmine[1][13]=0.04;
    rdcutsvalmine[2][13]=1.0;
    rdcutsvalmine[3][13]=0.30;
    rdcutsvalmine[4][13]=0.30;
    rdcutsvalmine[5][13]=0.50;
    rdcutsvalmine[6][13]=0.50;
    rdcutsvalmine[7][13]=0.0033;
    rdcutsvalmine[8][13]=0.85;
    rdcutsvalmine[9][13]=0.1;
    rdcutsvalmine[10][13]=0.1;
    rdcutsvalmine[11][13]=0.05;
    rdcutsvalmine[12][13]=100.;
    rdcutsvalmine[13][13]=1;
    rdcutsvalmine[14][13]=0.99;		//0.9;
    rdcutsvalmine[15][13]=3.5;		//0.0

    //>25 (ptbin14)
    rdcutsvalmine[0][14]=0.074;
    rdcutsvalmine[1][14]=0.8;
    rdcutsvalmine[2][14]=1.0;
    rdcutsvalmine[3][14]=0.30;
    rdcutsvalmine[4][14]=0.30;
    rdcutsvalmine[5][14]=0.50;
    rdcutsvalmine[6][14]=0.50;
    rdcutsvalmine[7][14]=0.0033;
    rdcutsvalmine[8][14]=0.7;
    rdcutsvalmine[9][14]=0.1;
    rdcutsvalmine[10][14]=0.1;
    rdcutsvalmine[11][14]=0.05;
    rdcutsvalmine[12][14]=100.;
    rdcutsvalmine[13][14]=1;
    rdcutsvalmine[14][14]=0.9;
    rdcutsvalmine[15][14]=0.5;
*/

/*
// Cuts 07C_2B (corresponds to 07C from train19, with corrections form 160412_CutsOpt.root and cuts 'smoothened')
    //0-0.5
    rdcutsvalmine[0][0]=0.032; 		// inv mass
    rdcutsvalmine[1][0]=0.02;		// dca
    rdcutsvalmine[2][0]=0.8;		// cos theta star
    rdcutsvalmine[3][0]=1.333;		// ptK
    rdcutsvalmine[4][0]=1.0;		// pt PI
    rdcutsvalmine[5][0]=0.07;		//  d0 K
    rdcutsvalmine[6][0]=0.07;		// d0 Pi
    rdcutsvalmine[7][0]=-0.000362;	// d0d0
    rdcutsvalmine[8][0]=0.99;		// costheta point
    rdcutsvalmine[9][0]=0.3;            // inv mass half width of D*
    rdcutsvalmine[10][0]=0.15;		// half width of MKpiPi -Md0
    rdcutsvalmine[11][0]=0.05;		// pt min of soft pion
    rdcutsvalmine[12][0]=0.5;		// pt max of soft pion
    rdcutsvalmine[13][0]=0.5;		// theta between soft pion and d0 plane
    rdcutsvalmine[14][0]=0.998;		// cos theta XY
    rdcutsvalmine[15][0]=7.0;		// norm decay length

     //0.5-1
    rdcutsvalmine[0][1]=0.032;
    rdcutsvalmine[1][1]=0.02;
    rdcutsvalmine[2][1]=0.8;
    rdcutsvalmine[3][1]=1.333;
    rdcutsvalmine[4][1]=1.0;
    rdcutsvalmine[5][1]=0.1;		//0.07;
    rdcutsvalmine[6][1]=0.1;		//0.07;
    rdcutsvalmine[7][1]=-0.000362;
    rdcutsvalmine[8][1]=0.99;
    rdcutsvalmine[9][1]=0.3;
    rdcutsvalmine[10][1]=0.15;
    rdcutsvalmine[11][1]=0.05;
    rdcutsvalmine[12][1]=0.5;
    rdcutsvalmine[13][1]=0.5;
    rdcutsvalmine[14][1]=0.998;
    rdcutsvalmine[15][1]=7.0;

    //1-2
    rdcutsvalmine[0][2]=0.032;
    rdcutsvalmine[1][2]=0.025;
    rdcutsvalmine[2][2]=0.8;
    rdcutsvalmine[3][2]=1.0;
    rdcutsvalmine[4][2]=1.0;
    rdcutsvalmine[5][2]=0.1;
    rdcutsvalmine[6][2]=0.1;
    rdcutsvalmine[7][2]=-0.0002;
    rdcutsvalmine[8][2]=0.995;
    rdcutsvalmine[9][2]=0.3;
    rdcutsvalmine[10][2]=0.15;
    rdcutsvalmine[11][2]=0.05;
    rdcutsvalmine[12][2]=0.5;
    rdcutsvalmine[13][2]=0.5;
    rdcutsvalmine[14][2]=0.998;
    rdcutsvalmine[15][2]=9;
    //2-3
    rdcutsvalmine[0][3]=0.032;
    rdcutsvalmine[1][3]=0.025;
    rdcutsvalmine[2][3]=0.8;
    rdcutsvalmine[3][3]=1.0;
    rdcutsvalmine[4][3]=1.0;
    rdcutsvalmine[5][3]=0.1;
    rdcutsvalmine[6][3]=0.1;
    rdcutsvalmine[7][3]=-0.0003;
    rdcutsvalmine[8][3]=0.996;
    rdcutsvalmine[9][3]=0.3;
    rdcutsvalmine[10][3]=0.15;
    rdcutsvalmine[11][3]=0.05;
    rdcutsvalmine[12][3]=0.5;
    rdcutsvalmine[13][3]=0.5;
    rdcutsvalmine[14][3]=0.998;
    rdcutsvalmine[15][3]=8.0;

    //3-4
    rdcutsvalmine[0][4]=0.032;
    rdcutsvalmine[1][4]=0.025;
    rdcutsvalmine[2][4]=0.8;
    rdcutsvalmine[3][4]=1.0;
    rdcutsvalmine[4][4]=1.0;
    rdcutsvalmine[5][4]=0.1;
    rdcutsvalmine[6][4]=0.1;
    rdcutsvalmine[7][4]=-0.000362;
    rdcutsvalmine[8][4]=0.996;
    rdcutsvalmine[9][4]=0.3;
    rdcutsvalmine[10][4]=0.15;
    rdcutsvalmine[11][4]=0.05;
    rdcutsvalmine[12][4]=0.5;
    rdcutsvalmine[13][4]=0.5;
    rdcutsvalmine[14][4]=0.998;
    rdcutsvalmine[15][4]=7.0;

    //4-5
    rdcutsvalmine[0][5]=0.032;
    rdcutsvalmine[1][5]=0.02;
    rdcutsvalmine[2][5]=0.8;
    rdcutsvalmine[3][5]=1.;
    rdcutsvalmine[4][5]=1.;
    rdcutsvalmine[5][5]=0.1;
    rdcutsvalmine[6][5]=0.1;
    rdcutsvalmine[7][5]=-0.00028;
    rdcutsvalmine[8][5]=0.996;
    rdcutsvalmine[9][5]=0.2;
    rdcutsvalmine[10][5]=0.15;
    rdcutsvalmine[11][5]=0.05;
    rdcutsvalmine[12][5]=10;
    rdcutsvalmine[13][5]=0.5;
    rdcutsvalmine[14][5]=0.998;
    rdcutsvalmine[15][5]=6.5;

    //5-6
    rdcutsvalmine[0][6]=0.04;
    rdcutsvalmine[1][6]=0.02;
    rdcutsvalmine[2][6]=1.;
    rdcutsvalmine[3][6]=1.;
    rdcutsvalmine[4][6]=1.;
    rdcutsvalmine[5][6]=0.1;
    rdcutsvalmine[6][6]=0.1;
    rdcutsvalmine[7][6]=-0.00025;
    rdcutsvalmine[8][6]=0.994;
    rdcutsvalmine[9][6]=0.15;
    rdcutsvalmine[10][6]=0.15;
    rdcutsvalmine[11][6]=0.2;
    rdcutsvalmine[12][6]=10.;
    rdcutsvalmine[13][6]=0.5;
    rdcutsvalmine[14][6]=0.998;
    rdcutsvalmine[15][6]=6.0;

    //6-7
    rdcutsvalmine[0][7]=0.036;
    rdcutsvalmine[1][7]=0.02;
    rdcutsvalmine[2][7]=1.0;
    rdcutsvalmine[3][7]=1.;
    rdcutsvalmine[4][7]=1.;
    rdcutsvalmine[5][7]=0.1;
    rdcutsvalmine[6][7]=0.1;
    rdcutsvalmine[7][7]=-0.0002;
    rdcutsvalmine[8][7]=0.994;
    rdcutsvalmine[9][7]=0.3;
    rdcutsvalmine[10][7]=0.1;
    rdcutsvalmine[11][7]=0.2;
    rdcutsvalmine[12][7]=100.;
    rdcutsvalmine[13][7]=0.5;
    rdcutsvalmine[14][7]=0.998;
    rdcutsvalmine[15][7]=5.5;

    //7-8
    rdcutsvalmine[0][8]=0.036;
    rdcutsvalmine[1][8]=0.02;
    rdcutsvalmine[2][8]=1.0;
    rdcutsvalmine[3][8]=1.;
    rdcutsvalmine[4][8]=1.;
    rdcutsvalmine[5][8]=0.12;
    rdcutsvalmine[6][8]=0.12;
    rdcutsvalmine[7][8]=-0.000165;
    rdcutsvalmine[8][8]=0.992;
    rdcutsvalmine[9][8]=0.3;
    rdcutsvalmine[10][8]=0.15;
    rdcutsvalmine[11][8]=0.2;
    rdcutsvalmine[12][8]=100.;
    rdcutsvalmine[13][8]=0.5;
    rdcutsvalmine[14][8]=0.998;
    rdcutsvalmine[15][8]=5.;

    //8-10
    rdcutsvalmine[0][9]=0.055;
    rdcutsvalmine[1][9]=0.02;
    rdcutsvalmine[2][9]=1.0;
    rdcutsvalmine[3][9]=0.9;
    rdcutsvalmine[4][9]=0.9;
    rdcutsvalmine[5][9]=0.15;
    rdcutsvalmine[6][9]=0.15;
    rdcutsvalmine[7][9]=-0.000065;
    rdcutsvalmine[8][9]=0.992;
    rdcutsvalmine[9][9]=0.3;
    rdcutsvalmine[10][9]=0.15;
    rdcutsvalmine[11][9]=0.34;
    rdcutsvalmine[12][9]=100.;
    rdcutsvalmine[13][9]=0.5;
    rdcutsvalmine[14][9]=0.998;
    rdcutsvalmine[15][9]=4.5;

    //10-12
    rdcutsvalmine[0][10]=0.055;
    rdcutsvalmine[1][10]=0.02;
    rdcutsvalmine[2][10]=1.0;
    rdcutsvalmine[3][10]=0.9;
    rdcutsvalmine[4][10]=0.9;
    rdcutsvalmine[5][10]=0.15;
    rdcutsvalmine[6][10]=0.15;
    rdcutsvalmine[7][10]=0.0;
    rdcutsvalmine[8][10]=0.98;
    rdcutsvalmine[9][10]=0.3;
    rdcutsvalmine[10][10]=0.15;
    rdcutsvalmine[11][10]=0.35;
    rdcutsvalmine[12][10]=100.;
    rdcutsvalmine[13][10]=0.5;
    rdcutsvalmine[14][10]=0.998;
    rdcutsvalmine[15][10]=4.;


    //12-16 (ptbin11)
    rdcutsvalmine[0][11]=0.074;
    rdcutsvalmine[1][11]=0.025;
    rdcutsvalmine[2][11]=0.8;
    rdcutsvalmine[3][11]=0.3;
    rdcutsvalmine[4][11]=0.3;
    rdcutsvalmine[5][11]=0.2;
    rdcutsvalmine[6][11]=0.2;
    rdcutsvalmine[7][11]=0.0000;
    rdcutsvalmine[8][11]=0.97;
    rdcutsvalmine[9][11]=0.3;
    rdcutsvalmine[10][11]=0.1;
    rdcutsvalmine[11][11]=0.45;
    rdcutsvalmine[12][11]=100.;
    rdcutsvalmine[13][11]=1;
    rdcutsvalmine[14][11]=0.998;
    rdcutsvalmine[15][11]=4.0;



    //16-20 (ptbin12)
    rdcutsvalmine[0][12]=0.074;
    rdcutsvalmine[1][12]=0.035;
    rdcutsvalmine[2][12]=1.0;
    rdcutsvalmine[3][12]=0.30;
    rdcutsvalmine[4][12]=0.30;
    rdcutsvalmine[5][12]=0.50;
    rdcutsvalmine[6][12]=0.50;
    rdcutsvalmine[7][12]=0.0033;
    rdcutsvalmine[8][12]=0.95;
    rdcutsvalmine[9][12]=0.1;
    rdcutsvalmine[10][12]=0.1;
    rdcutsvalmine[11][12]=0.05;
    rdcutsvalmine[12][12]=100.;
    rdcutsvalmine[13][12]=1;
    rdcutsvalmine[14][12]=0.99;
    rdcutsvalmine[15][12]=3.5;

    //20-24 (ptbin13)
    rdcutsvalmine[0][13]=0.074;
    rdcutsvalmine[1][13]=0.04;
    rdcutsvalmine[2][13]=1.0;
    rdcutsvalmine[3][13]=0.30;
    rdcutsvalmine[4][13]=0.30;
    rdcutsvalmine[5][13]=0.50;
    rdcutsvalmine[6][13]=0.50;
    rdcutsvalmine[7][13]=0.0033;
    rdcutsvalmine[8][13]=0.85;
    rdcutsvalmine[9][13]=0.1;
    rdcutsvalmine[10][13]=0.1;
    rdcutsvalmine[11][13]=0.05;
    rdcutsvalmine[12][13]=100.;
    rdcutsvalmine[13][13]=1;
    rdcutsvalmine[14][13]=0.99;		//0.9;
    rdcutsvalmine[15][13]=3.5;		//0.0

    //>25 (ptbin14)
    rdcutsvalmine[0][14]=0.074;
    rdcutsvalmine[1][14]=0.8;
    rdcutsvalmine[2][14]=1.0;
    rdcutsvalmine[3][14]=0.30;
    rdcutsvalmine[4][14]=0.30;
    rdcutsvalmine[5][14]=0.50;
    rdcutsvalmine[6][14]=0.50;
    rdcutsvalmine[7][14]=0.0033;
    rdcutsvalmine[8][14]=0.7;
    rdcutsvalmine[9][14]=0.1;
    rdcutsvalmine[10][14]=0.1;
    rdcutsvalmine[11][14]=0.05;
    rdcutsvalmine[12][14]=100.;
    rdcutsvalmine[13][14]=1;
    rdcutsvalmine[14][14]=0.9;
    rdcutsvalmine[15][14]=0.5;
*/

/*

// Cuts 07C_2A (corresponds to 07C from train19, with corrections form 160412_CutsOpt.root)
    //0-0.5
    rdcutsvalmine[0][0]=0.032; 		// inv mass
    rdcutsvalmine[1][0]=0.02;		// dca
    rdcutsvalmine[2][0]=0.8;		// cos theta star
    rdcutsvalmine[3][0]=1.333;		// ptK
    rdcutsvalmine[4][0]=1.0;		// pt PI
    rdcutsvalmine[5][0]=0.07;		//  d0 K
    rdcutsvalmine[6][0]=0.07;		// d0 Pi
    rdcutsvalmine[7][0]=-0.000362;	// d0d0
    rdcutsvalmine[8][0]=0.99;		// costheta point
    rdcutsvalmine[9][0]=0.3;            // inv mass half width of D*
    rdcutsvalmine[10][0]=0.15;		// half width of MKpiPi -Md0
    rdcutsvalmine[11][0]=0.05;		// pt min of soft pion
    rdcutsvalmine[12][0]=0.5;		// pt max of soft pion
    rdcutsvalmine[13][0]=0.5;		// theta between soft pion and d0 plane
    rdcutsvalmine[14][0]=0.998;		// cos theta XY
    rdcutsvalmine[15][0]=7.0;		// norm decay length

     //0.5-1
    rdcutsvalmine[0][1]=0.032;
    rdcutsvalmine[1][1]=0.02;
    rdcutsvalmine[2][1]=0.8;
    rdcutsvalmine[3][1]=1.333;
    rdcutsvalmine[4][1]=1.0;
    rdcutsvalmine[5][1]=0.1;		//0.07;
    rdcutsvalmine[6][1]=0.1;		//0.07;
    rdcutsvalmine[7][1]=-0.000362;
    rdcutsvalmine[8][1]=0.99;
    rdcutsvalmine[9][1]=0.3;
    rdcutsvalmine[10][1]=0.15;
    rdcutsvalmine[11][1]=0.05;
    rdcutsvalmine[12][1]=0.5;
    rdcutsvalmine[13][1]=0.5;
    rdcutsvalmine[14][1]=0.998;
    rdcutsvalmine[15][1]=7.0;

    //1-2
    rdcutsvalmine[0][2]=0.032;
    rdcutsvalmine[1][2]=0.015;
    rdcutsvalmine[2][2]=0.8;
    rdcutsvalmine[3][2]=1.0;
    rdcutsvalmine[4][2]=1.0;
    rdcutsvalmine[5][2]=0.1;
    rdcutsvalmine[6][2]=0.1;
    rdcutsvalmine[7][2]=-0.0002;
    rdcutsvalmine[8][2]=0.995;
    rdcutsvalmine[9][2]=0.3;
    rdcutsvalmine[10][2]=0.15;
    rdcutsvalmine[11][2]=0.05;
    rdcutsvalmine[12][2]=0.5;
    rdcutsvalmine[13][2]=0.5;
    rdcutsvalmine[14][2]=0.998;
    rdcutsvalmine[15][2]=9;
    //2-3
    rdcutsvalmine[0][3]=0.032;
    rdcutsvalmine[1][3]=0.035;
    rdcutsvalmine[2][3]=0.8;
    rdcutsvalmine[3][3]=1.0;
    rdcutsvalmine[4][3]=1.0;
    rdcutsvalmine[5][3]=0.1;
    rdcutsvalmine[6][3]=0.1;
    rdcutsvalmine[7][3]=-0.0003;
    rdcutsvalmine[8][3]=0.991;
    rdcutsvalmine[9][3]=0.3;
    rdcutsvalmine[10][3]=0.15;
    rdcutsvalmine[11][3]=0.05;
    rdcutsvalmine[12][3]=0.5;
    rdcutsvalmine[13][3]=0.5;
    rdcutsvalmine[14][3]=0.998;
    rdcutsvalmine[15][3]=8.0;

    //3-4
    rdcutsvalmine[0][4]=0.032;
    rdcutsvalmine[1][4]=0.025;
    rdcutsvalmine[2][4]=0.8;
    rdcutsvalmine[3][4]=1.0;
    rdcutsvalmine[4][4]=1.0;
    rdcutsvalmine[5][4]=0.1;
    rdcutsvalmine[6][4]=0.1;
    rdcutsvalmine[7][4]=-0.000362;
    rdcutsvalmine[8][4]=0.996;
    rdcutsvalmine[9][4]=0.3;
    rdcutsvalmine[10][4]=0.15;
    rdcutsvalmine[11][4]=0.05;
    rdcutsvalmine[12][4]=0.5;
    rdcutsvalmine[13][4]=0.5;
    rdcutsvalmine[14][4]=0.998;
    rdcutsvalmine[15][4]=7.0;

    //4-5
    rdcutsvalmine[0][5]=0.032;
    rdcutsvalmine[1][5]=0.02;
    rdcutsvalmine[2][5]=0.8;
    rdcutsvalmine[3][5]=1.;
    rdcutsvalmine[4][5]=1.;
    rdcutsvalmine[5][5]=0.1;
    rdcutsvalmine[6][5]=0.1;
    rdcutsvalmine[7][5]=-0.000267;
    rdcutsvalmine[8][5]=0.996;
    rdcutsvalmine[9][5]=0.2;
    rdcutsvalmine[10][5]=0.15;
    rdcutsvalmine[11][5]=0.05;
    rdcutsvalmine[12][5]=10;
    rdcutsvalmine[13][5]=0.5;
    rdcutsvalmine[14][5]=0.998;
    rdcutsvalmine[15][5]=6.5;

    //5-6
    rdcutsvalmine[0][6]=0.04;
    rdcutsvalmine[1][6]=0.02;
    rdcutsvalmine[2][6]=1.;
    rdcutsvalmine[3][6]=1.;
    rdcutsvalmine[4][6]=1.;
    rdcutsvalmine[5][6]=0.1;
    rdcutsvalmine[6][6]=0.1;
    rdcutsvalmine[7][6]=-0.00025;
    rdcutsvalmine[8][6]=0.99;
    rdcutsvalmine[9][6]=0.15;
    rdcutsvalmine[10][6]=0.15;
    rdcutsvalmine[11][6]=0.2;
    rdcutsvalmine[12][6]=10.;
    rdcutsvalmine[13][6]=0.5;
    rdcutsvalmine[14][6]=0.998;
    rdcutsvalmine[15][6]=7.0;

    //6-7
    rdcutsvalmine[0][7]=0.036;
    rdcutsvalmine[1][7]=0.023;
    rdcutsvalmine[2][7]=1.0;
    rdcutsvalmine[3][7]=1.;
    rdcutsvalmine[4][7]=1.;
    rdcutsvalmine[5][7]=0.1;
    rdcutsvalmine[6][7]=0.1;
    rdcutsvalmine[7][7]=-0.0002;
    rdcutsvalmine[8][7]=0.994;
    rdcutsvalmine[9][7]=0.3;
    rdcutsvalmine[10][7]=0.1;
    rdcutsvalmine[11][7]=0.2;
    rdcutsvalmine[12][7]=100.;
    rdcutsvalmine[13][7]=0.5;
    rdcutsvalmine[14][7]=0.998;
    rdcutsvalmine[15][7]=5.5;

    //7-8
    rdcutsvalmine[0][8]=0.036;
    rdcutsvalmine[1][8]=0.02;
    rdcutsvalmine[2][8]=1.0;
    rdcutsvalmine[3][8]=1.;
    rdcutsvalmine[4][8]=1.;
    rdcutsvalmine[5][8]=0.12;
    rdcutsvalmine[6][8]=0.12;
    rdcutsvalmine[7][8]=-0.000165;
    rdcutsvalmine[8][8]=0.989;
    rdcutsvalmine[9][8]=0.3;
    rdcutsvalmine[10][8]=0.15;
    rdcutsvalmine[11][8]=0.2;
    rdcutsvalmine[12][8]=100.;
    rdcutsvalmine[13][8]=0.5;
    rdcutsvalmine[14][8]=0.998;
    rdcutsvalmine[15][8]=3.5;

    //8-10
    rdcutsvalmine[0][9]=0.055;
    rdcutsvalmine[1][9]=0.02;
    rdcutsvalmine[2][9]=1.0;
    rdcutsvalmine[3][9]=1.0;
    rdcutsvalmine[4][9]=1.0;
    rdcutsvalmine[5][9]=0.15;
    rdcutsvalmine[6][9]=0.15;
    rdcutsvalmine[7][9]=-0.000065;
    rdcutsvalmine[8][9]=0.995;
    rdcutsvalmine[9][9]=0.3;
    rdcutsvalmine[10][9]=0.15;
    rdcutsvalmine[11][9]=0.34;
    rdcutsvalmine[12][9]=100.;
    rdcutsvalmine[13][9]=0.5;
    rdcutsvalmine[14][9]=0.998;
    rdcutsvalmine[15][9]=5.5;

    //10-12
    rdcutsvalmine[0][10]=0.055;
    rdcutsvalmine[1][10]=0.02;
    rdcutsvalmine[2][10]=1.0;
    rdcutsvalmine[3][10]=1.0;
    rdcutsvalmine[4][10]=1.0;
    rdcutsvalmine[5][10]=0.15;
    rdcutsvalmine[6][10]=0.15;
    rdcutsvalmine[7][10]=0.0;
    rdcutsvalmine[8][10]=0.98;
    rdcutsvalmine[9][10]=0.3;
    rdcutsvalmine[10][10]=0.15;
    rdcutsvalmine[11][10]=0.35;
    rdcutsvalmine[12][10]=100.;
    rdcutsvalmine[13][10]=0.5;
    rdcutsvalmine[14][10]=0.998;
    rdcutsvalmine[15][10]=2.5;


    //12-16 (ptbin11)
    rdcutsvalmine[0][11]=0.074;
    rdcutsvalmine[1][11]=0.035;
    rdcutsvalmine[2][11]=0.8;
    rdcutsvalmine[3][11]=0.3;//1.0; This has been added as a hunch
    rdcutsvalmine[4][11]=0.3;//1.0;
    rdcutsvalmine[5][11]=0.2;
    rdcutsvalmine[6][11]=0.2;
    rdcutsvalmine[7][11]=0.0000;
    rdcutsvalmine[8][11]=0.97;
    rdcutsvalmine[9][11]=0.3;
    rdcutsvalmine[10][11]=0.1;
    rdcutsvalmine[11][11]=0.45;
    rdcutsvalmine[12][11]=100.;
    rdcutsvalmine[13][11]=1;
    rdcutsvalmine[14][11]=0.998;
    rdcutsvalmine[15][11]=4.0;



    //16-20 (ptbin12)
    rdcutsvalmine[0][12]=0.074;
    rdcutsvalmine[1][12]=0.025;
    rdcutsvalmine[2][12]=1.0;
    rdcutsvalmine[3][12]=0.30;
    rdcutsvalmine[4][12]=0.30;
    rdcutsvalmine[5][12]=0.50;
    rdcutsvalmine[6][12]=0.50;
    rdcutsvalmine[7][12]=0.0033;
    rdcutsvalmine[8][12]=0.9;
    rdcutsvalmine[9][12]=0.1;
    rdcutsvalmine[10][12]=0.1;
    rdcutsvalmine[11][12]=0.05;
    rdcutsvalmine[12][12]=100.;
    rdcutsvalmine[13][12]=1;
    rdcutsvalmine[14][12]=0.99;
    rdcutsvalmine[15][12]=3.5;

    //20-24 (ptbin13)
    rdcutsvalmine[0][13]=0.074;
    rdcutsvalmine[1][13]=0.055;
    rdcutsvalmine[2][13]=1.0;
    rdcutsvalmine[3][13]=0.30;
    rdcutsvalmine[4][13]=0.30;
    rdcutsvalmine[5][13]=0.50;
    rdcutsvalmine[6][13]=0.50;
    rdcutsvalmine[7][13]=0.0033;
    rdcutsvalmine[8][13]=0.85;
    rdcutsvalmine[9][13]=0.1;
    rdcutsvalmine[10][13]=0.1;
    rdcutsvalmine[11][13]=0.05;
    rdcutsvalmine[12][13]=100.;
    rdcutsvalmine[13][13]=1;
    rdcutsvalmine[14][13]=0.99;		//0.9;
    rdcutsvalmine[15][13]=3.5;		//0.0

    //>25 (ptbin14)
    rdcutsvalmine[0][14]=0.074;
    rdcutsvalmine[1][14]=0.8;
    rdcutsvalmine[2][14]=1.0;
    rdcutsvalmine[3][14]=0.30;
    rdcutsvalmine[4][14]=0.30;
    rdcutsvalmine[5][14]=0.50;
    rdcutsvalmine[6][14]=0.50;
    rdcutsvalmine[7][14]=0.0033;
    rdcutsvalmine[8][14]=0.7;
    rdcutsvalmine[9][14]=0.1;
    rdcutsvalmine[10][14]=0.1;
    rdcutsvalmine[11][14]=0.05;
    rdcutsvalmine[12][14]=100.;
    rdcutsvalmine[13][14]=1;
    rdcutsvalmine[14][14]=0.9;
    rdcutsvalmine[15][14]=0.5;
*/
/*
// Cuts 07C, stronger than 07A everywhere
    //0-0.5
    rdcutsvalmine[0][0]=0.032; 		// inv mass
    rdcutsvalmine[1][0]=0.02;		// dca
    rdcutsvalmine[2][0]=0.8;		// cos theta star
    rdcutsvalmine[3][0]=1.333;		// ptK
    rdcutsvalmine[4][0]=1.0;		// pt PI
    rdcutsvalmine[5][0]=0.07;		//  d0 K
    rdcutsvalmine[6][0]=0.07;		// d0 Pi
    rdcutsvalmine[7][0]=-0.000362;	// d0d0
    rdcutsvalmine[8][0]=0.99;		// costheta point
    rdcutsvalmine[9][0]=0.3;            // inv mass half width of D*
    rdcutsvalmine[10][0]=0.15;		// half width of MKpiPi -Md0
    rdcutsvalmine[11][0]=0.05;		// pt min of soft pion
    rdcutsvalmine[12][0]=0.5;		// pt max of soft pion
    rdcutsvalmine[13][0]=0.5;		// theta between soft pion and d0 plane
    rdcutsvalmine[14][0]=0.998;		// cos theta XY
    rdcutsvalmine[15][0]=7.0;		// norm decay length

     //0.5-1
    rdcutsvalmine[0][1]=0.032;
    rdcutsvalmine[1][1]=0.02;
    rdcutsvalmine[2][1]=0.8;
    rdcutsvalmine[3][1]=1.333;
    rdcutsvalmine[4][1]=1.0;
    rdcutsvalmine[5][1]=0.1;		//0.07;
    rdcutsvalmine[6][1]=0.1;		//0.07;
    rdcutsvalmine[7][1]=-0.000362;
    rdcutsvalmine[8][1]=0.99;
    rdcutsvalmine[9][1]=0.3;
    rdcutsvalmine[10][1]=0.15;
    rdcutsvalmine[11][1]=0.05;
    rdcutsvalmine[12][1]=0.5;
    rdcutsvalmine[13][1]=0.5;
    rdcutsvalmine[14][1]=0.998;
    rdcutsvalmine[15][1]=7.0;

    //1-2
    rdcutsvalmine[0][2]=0.032;
    rdcutsvalmine[1][2]=0.02;
    rdcutsvalmine[2][2]=0.8;
    rdcutsvalmine[3][2]=1.0;
    rdcutsvalmine[4][2]=1.0;
    rdcutsvalmine[5][2]=0.1;		//0.07;
    rdcutsvalmine[6][2]=0.1;		//0.07;
    rdcutsvalmine[7][2]=-0.000362;
    rdcutsvalmine[8][2]=0.996;
    rdcutsvalmine[9][2]=0.3;
    rdcutsvalmine[10][2]=0.15;
    rdcutsvalmine[11][2]=0.05;
    rdcutsvalmine[12][2]=0.5;
    rdcutsvalmine[13][2]=0.5;
    rdcutsvalmine[14][2]=0.998;
    rdcutsvalmine[15][2]=8.5;
    //2-3
    rdcutsvalmine[0][3]=0.032;
    rdcutsvalmine[1][3]=0.02;
    rdcutsvalmine[2][3]=0.8;
    rdcutsvalmine[3][3]=1.0;
    rdcutsvalmine[4][3]=1.0;
    rdcutsvalmine[5][3]=0.1;		//0.07;
    rdcutsvalmine[6][3]=0.1;		//0.07;
    rdcutsvalmine[7][3]=-0.000362;
    rdcutsvalmine[8][3]=0.996;
    rdcutsvalmine[9][3]=0.3;
    rdcutsvalmine[10][3]=0.15;
    rdcutsvalmine[11][3]=0.05;
    rdcutsvalmine[12][3]=0.5;
    rdcutsvalmine[13][3]=0.5;
    rdcutsvalmine[14][3]=0.998;
    rdcutsvalmine[15][3]=8.0;

    //3-4
    rdcutsvalmine[0][4]=0.032;
    rdcutsvalmine[1][4]=0.02;
    rdcutsvalmine[2][4]=0.8;
    rdcutsvalmine[3][4]=1.0;		//1.333;
    rdcutsvalmine[4][4]=1.0;
    rdcutsvalmine[5][4]=0.1;		//0.07;
    rdcutsvalmine[6][4]=0.1;		//0.07;
    rdcutsvalmine[7][4]=-0.000362;
    rdcutsvalmine[8][4]=0.996;		//0.99;
    rdcutsvalmine[9][4]=0.3;
    rdcutsvalmine[10][4]=0.15;
    rdcutsvalmine[11][4]=0.05;
    rdcutsvalmine[12][4]=0.5;
    rdcutsvalmine[13][4]=0.5;
    rdcutsvalmine[14][4]=0.998;
    rdcutsvalmine[15][4]=8.0;		//7.0;

    //4-5
    rdcutsvalmine[0][5]=0.032;
    rdcutsvalmine[1][5]=0.025;
    rdcutsvalmine[2][5]=0.8;
    rdcutsvalmine[3][5]=1.;
    rdcutsvalmine[4][5]=1.;
    rdcutsvalmine[5][5]=0.1;		//0.075;
    rdcutsvalmine[6][5]=0.1;		//0.075;
    rdcutsvalmine[7][5]=-0.000285;
    rdcutsvalmine[8][5]=0.995;		//0.99;
    rdcutsvalmine[9][5]=0.15;
    rdcutsvalmine[10][5]=0.15;
    rdcutsvalmine[11][5]=0.05;
    rdcutsvalmine[12][5]=10;
    rdcutsvalmine[13][5]=0.5;
    rdcutsvalmine[14][5]=0.998;
    rdcutsvalmine[15][5]=7.0;		//6.;

    //5-6
    rdcutsvalmine[0][6]=0.04;
    rdcutsvalmine[1][6]=0.03;
    rdcutsvalmine[2][6]=1.;
    rdcutsvalmine[3][6]=1.;
    rdcutsvalmine[4][6]=1.;
    rdcutsvalmine[5][6]=0.1;		//0.075;
    rdcutsvalmine[6][6]=0.1;		//0.075;
    rdcutsvalmine[7][6]=-0.00028;
    rdcutsvalmine[8][6]=0.993;		//0.99;
    rdcutsvalmine[9][6]=0.15;
    rdcutsvalmine[10][6]=0.15;
    rdcutsvalmine[11][6]=0.2;
    rdcutsvalmine[12][6]=10.;
    rdcutsvalmine[13][6]=0.5;
    rdcutsvalmine[14][6]=0.998;
    rdcutsvalmine[15][6]=6.0;		//5.;

    //6-7
    rdcutsvalmine[0][7]=0.036;
    rdcutsvalmine[1][7]=0.03750;
    rdcutsvalmine[2][7]=1.0;
    rdcutsvalmine[3][7]=1.;
    rdcutsvalmine[4][7]=1.;
    rdcutsvalmine[5][7]=0.1;
    rdcutsvalmine[6][7]=0.1;
    rdcutsvalmine[7][7]=-0.00022;
    rdcutsvalmine[8][7]=0.99;		//0.985;
    rdcutsvalmine[9][7]=0.3;
    rdcutsvalmine[10][7]=0.1;
    rdcutsvalmine[11][7]=0.2;
    rdcutsvalmine[12][7]=100.;
    rdcutsvalmine[13][7]=0.5;
    rdcutsvalmine[14][7]=0.998;
    rdcutsvalmine[15][7]=5.0;		//4.0;

    //7-8
    rdcutsvalmine[0][8]=0.036;
    rdcutsvalmine[1][8]=0.03750;
    rdcutsvalmine[2][8]=1.0;
    rdcutsvalmine[3][8]=1.;
    rdcutsvalmine[4][8]=1.;
    rdcutsvalmine[5][8]=0.12;
    rdcutsvalmine[6][8]=0.12;
    rdcutsvalmine[7][8]=-0.000165;
    rdcutsvalmine[8][8]=0.99;		//0.98;
    rdcutsvalmine[9][8]=0.3;
    rdcutsvalmine[10][8]=0.15;
    rdcutsvalmine[11][8]=0.2;
    rdcutsvalmine[12][8]=100.;
    rdcutsvalmine[13][8]=0.5;
    rdcutsvalmine[14][8]=0.998;
    rdcutsvalmine[15][8]=5.0;		//4.;

    //8-10
    rdcutsvalmine[0][9]=0.055;
    rdcutsvalmine[1][9]=0.04;
    rdcutsvalmine[2][9]=1.0;
    rdcutsvalmine[3][9]=1.0;
    rdcutsvalmine[4][9]=1.0;
    rdcutsvalmine[5][9]=0.15;
    rdcutsvalmine[6][9]=0.15;
    rdcutsvalmine[7][9]=-0.00006;
    rdcutsvalmine[8][9]=0.99;		//0.98;
    rdcutsvalmine[9][9]=0.3;
    rdcutsvalmine[10][9]=0.15;
    rdcutsvalmine[11][9]=0.34;
    rdcutsvalmine[12][9]=100.;
    rdcutsvalmine[13][9]=0.5;
    rdcutsvalmine[14][9]=0.998;
    rdcutsvalmine[15][9]=5.0;		//4.;

    //10-12
    rdcutsvalmine[0][10]=0.055;
    rdcutsvalmine[1][10]=0.03;
    rdcutsvalmine[2][10]=1.0;
    rdcutsvalmine[3][10]=1.0;
    rdcutsvalmine[4][10]=1.0;
    rdcutsvalmine[5][10]=0.15;
    rdcutsvalmine[6][10]=0.15;
    rdcutsvalmine[7][10]=-0.00004;
    rdcutsvalmine[8][10]=0.99;
    rdcutsvalmine[9][10]=0.3;
    rdcutsvalmine[10][10]=0.15;
    rdcutsvalmine[11][10]=0.35;
    rdcutsvalmine[12][10]=100.;
    rdcutsvalmine[13][10]=0.5;
    rdcutsvalmine[14][10]=0.998;
    rdcutsvalmine[15][10]=5.0;		//4.;


    //12-16 (ptbin11)
    rdcutsvalmine[0][11]=0.074;
    rdcutsvalmine[1][11]=0.025;
    rdcutsvalmine[2][11]=0.8;
    rdcutsvalmine[3][11]=1.0;
    rdcutsvalmine[4][11]=1.0;
    rdcutsvalmine[5][11]=0.2;
    rdcutsvalmine[6][11]=0.2;
    rdcutsvalmine[7][11]=-0.000023;
    rdcutsvalmine[8][11]=0.97;
    rdcutsvalmine[9][11]=0.3;
    rdcutsvalmine[10][11]=0.1;
    rdcutsvalmine[11][11]=0.45;
    rdcutsvalmine[12][11]=100.;
    rdcutsvalmine[13][11]=1;
    rdcutsvalmine[14][11]=0.998;
    rdcutsvalmine[15][11]=4.0;		//3.;



    //16-20 (ptbin12)
    rdcutsvalmine[0][12]=0.074;
    rdcutsvalmine[1][12]=0.04;		//0.8;
    rdcutsvalmine[2][12]=1.0;
    rdcutsvalmine[3][12]=0.30;
    rdcutsvalmine[4][12]=0.30;
    rdcutsvalmine[5][12]=0.50;
    rdcutsvalmine[6][12]=0.50;
    rdcutsvalmine[7][12]=0.0;		//0.0033;
    rdcutsvalmine[8][12]=0.9;		//0.7;
    rdcutsvalmine[9][12]=0.1;
    rdcutsvalmine[10][12]=0.1;
    rdcutsvalmine[11][12]=0.05;
    rdcutsvalmine[12][12]=100.;
    rdcutsvalmine[13][12]=1;
    rdcutsvalmine[14][12]=0.99;		//0.9;
    rdcutsvalmine[15][12]=1.5;		//0.0;

    //20-24 (ptbin13)
    rdcutsvalmine[0][13]=0.074;
    rdcutsvalmine[1][13]=0.05;		//0.8;
    rdcutsvalmine[2][13]=1.0;
    rdcutsvalmine[3][13]=0.30;
    rdcutsvalmine[4][13]=0.30;
    rdcutsvalmine[5][13]=0.50;
    rdcutsvalmine[6][13]=0.50;
    rdcutsvalmine[7][13]=0.0001;	//0.0033;
    rdcutsvalmine[8][13]=0.8;		//0.7;
    rdcutsvalmine[9][13]=0.1;
    rdcutsvalmine[10][13]=0.1;
    rdcutsvalmine[11][13]=0.05;
    rdcutsvalmine[12][13]=100.;
    rdcutsvalmine[13][13]=1;
    rdcutsvalmine[14][13]=0.99;		//0.9;
    rdcutsvalmine[15][13]=1.0;		//0.0

    //>25 (ptbin14)
    rdcutsvalmine[0][14]=0.074;
    rdcutsvalmine[1][14]=0.8;
    rdcutsvalmine[2][14]=1.0;
    rdcutsvalmine[3][14]=0.30;
    rdcutsvalmine[4][14]=0.30;
    rdcutsvalmine[5][14]=0.50;
    rdcutsvalmine[6][14]=0.50;
    rdcutsvalmine[7][14]=0.0033;
    rdcutsvalmine[8][14]=0.7;
    rdcutsvalmine[9][14]=0.1;
    rdcutsvalmine[10][14]=0.1;
    rdcutsvalmine[11][14]=0.05;
    rdcutsvalmine[12][14]=100.;
    rdcutsvalmine[13][14]=1;
    rdcutsvalmine[14][14]=0.9;
    rdcutsvalmine[15][14]=0.0;
*/

/*

// Cuts 07B, stronger than 07A in bins up to 4 GeV/c and [16,24]
    //0-0.5
    rdcutsvalmine[0][0]=0.032;
    rdcutsvalmine[1][0]=0.02;
    rdcutsvalmine[2][0]=0.8;
    rdcutsvalmine[3][0]=1.333;
    rdcutsvalmine[4][0]=1.0;
    rdcutsvalmine[5][0]=0.07;
    rdcutsvalmine[6][0]=0.07;
    rdcutsvalmine[7][0]=-0.000362;
    rdcutsvalmine[8][0]=0.99;
    rdcutsvalmine[9][0]=0.3;
    rdcutsvalmine[10][0]=0.15;
    rdcutsvalmine[11][0]=0.05;
    rdcutsvalmine[12][0]=0.5;
    rdcutsvalmine[13][0]=0.5;
    rdcutsvalmine[14][0]=0.998;
    rdcutsvalmine[15][0]=7.0;

     //0.5-1
    rdcutsvalmine[0][1]=0.032;
    rdcutsvalmine[1][1]=0.02;
    rdcutsvalmine[2][1]=0.8;
    rdcutsvalmine[3][1]=1.333;
    rdcutsvalmine[4][1]=1.0;
    rdcutsvalmine[5][1]=0.07;
    rdcutsvalmine[6][1]=0.07;
    rdcutsvalmine[7][1]=-0.000362;
    rdcutsvalmine[8][1]=0.99;
    rdcutsvalmine[9][1]=0.3;
    rdcutsvalmine[10][1]=0.15;
    rdcutsvalmine[11][1]=0.05;
    rdcutsvalmine[12][1]=0.5;
    rdcutsvalmine[13][1]=0.5;
    rdcutsvalmine[14][1]=0.998;
    rdcutsvalmine[15][1]=7.0;

    //1-2
    rdcutsvalmine[0][2]=0.032;
    rdcutsvalmine[1][2]=0.02;
    rdcutsvalmine[2][2]=0.8;
    rdcutsvalmine[3][2]=1.0;
    rdcutsvalmine[4][2]=1.0;
    rdcutsvalmine[5][2]=0.07;
    rdcutsvalmine[6][2]=0.07;
    rdcutsvalmine[7][2]=-0.000362;
    rdcutsvalmine[8][2]=0.996;
    rdcutsvalmine[9][2]=0.3;
    rdcutsvalmine[10][2]=0.15;
    rdcutsvalmine[11][2]=0.05;
    rdcutsvalmine[12][2]=0.5;
    rdcutsvalmine[13][2]=0.5;
    rdcutsvalmine[14][2]=0.998;
    rdcutsvalmine[15][2]=8.5;
    //2-3
    rdcutsvalmine[0][3]=0.032;
    rdcutsvalmine[1][3]=0.02;
    rdcutsvalmine[2][3]=0.8;
    rdcutsvalmine[3][3]=1.0;
    rdcutsvalmine[4][3]=1.0;
    rdcutsvalmine[5][3]=0.07;
    rdcutsvalmine[6][3]=0.07;
    rdcutsvalmine[7][3]=-0.000362;
    rdcutsvalmine[8][3]=0.995;
    rdcutsvalmine[9][3]=0.3;
    rdcutsvalmine[10][3]=0.15;
    rdcutsvalmine[11][3]=0.05;
    rdcutsvalmine[12][3]=0.5;
    rdcutsvalmine[13][3]=0.5;
    rdcutsvalmine[14][3]=0.998;
    rdcutsvalmine[15][3]=8.0;

    //3-4
    rdcutsvalmine[0][4]=0.032;
    rdcutsvalmine[1][4]=0.02;
    rdcutsvalmine[2][4]=0.8;
    rdcutsvalmine[3][4]=1.0;		//1.333;
    rdcutsvalmine[4][4]=1.0;
    rdcutsvalmine[5][4]=0.07;
    rdcutsvalmine[6][4]=0.07;
    rdcutsvalmine[7][4]=-0.000362;
    rdcutsvalmine[8][4]=0.994;		//0.99;
    rdcutsvalmine[9][4]=0.3;
    rdcutsvalmine[10][4]=0.15;
    rdcutsvalmine[11][4]=0.05;
    rdcutsvalmine[12][4]=0.5;
    rdcutsvalmine[13][4]=0.5;
    rdcutsvalmine[14][4]=0.998;
    rdcutsvalmine[15][4]=7.0;

    //4-5
    rdcutsvalmine[0][5]=0.032;
    rdcutsvalmine[1][5]=0.025;
    rdcutsvalmine[2][5]=0.8;
    rdcutsvalmine[3][5]=1.;
    rdcutsvalmine[4][5]=1.;
    rdcutsvalmine[5][5]=0.075;
    rdcutsvalmine[6][5]=0.075;
    rdcutsvalmine[7][5]=-0.000285;
    rdcutsvalmine[8][5]=0.99;
    rdcutsvalmine[9][5]=0.15;
    rdcutsvalmine[10][5]=0.15;
    rdcutsvalmine[11][5]=0.05;
    rdcutsvalmine[12][5]=10;
    rdcutsvalmine[13][5]=0.5;
    rdcutsvalmine[14][5]=0.998;
    rdcutsvalmine[15][5]=6.;

    //5-6
    rdcutsvalmine[0][6]=0.04;
    rdcutsvalmine[1][6]=0.03;
    rdcutsvalmine[2][6]=1.;
    rdcutsvalmine[3][6]=1.;
    rdcutsvalmine[4][6]=1.;
    rdcutsvalmine[5][6]=0.075;
    rdcutsvalmine[6][6]=0.075;
    rdcutsvalmine[7][6]=-0.00028;
    rdcutsvalmine[8][6]=0.99;
    rdcutsvalmine[9][6]=0.15;
    rdcutsvalmine[10][6]=0.15;
    rdcutsvalmine[11][6]=0.2;
    rdcutsvalmine[12][6]=10.;
    rdcutsvalmine[13][6]=0.5;
    rdcutsvalmine[14][6]=0.998;
    rdcutsvalmine[15][6]=5.;

    //6-7
    rdcutsvalmine[0][7]=0.036;
    rdcutsvalmine[1][7]=0.03750;
    rdcutsvalmine[2][7]=1.0;
    rdcutsvalmine[3][7]=1.;
    rdcutsvalmine[4][7]=1.;
    rdcutsvalmine[5][7]=0.1;
    rdcutsvalmine[6][7]=0.1;
    rdcutsvalmine[7][7]=-0.00022;
    rdcutsvalmine[8][7]=0.985;
    rdcutsvalmine[9][7]=0.3;
    rdcutsvalmine[10][7]=0.1;
    rdcutsvalmine[11][7]=0.2;
    rdcutsvalmine[12][7]=100.;
    rdcutsvalmine[13][7]=0.5;
    rdcutsvalmine[14][7]=0.998;
    rdcutsvalmine[15][7]=4.0;

    //7-8
    rdcutsvalmine[0][8]=0.036;
    rdcutsvalmine[1][8]=0.03750;
    rdcutsvalmine[2][8]=1.0;
    rdcutsvalmine[3][8]=1.;
    rdcutsvalmine[4][8]=1.;
    rdcutsvalmine[5][8]=0.12;
    rdcutsvalmine[6][8]=0.12;
    rdcutsvalmine[7][8]=-0.000165;
    rdcutsvalmine[8][8]=0.98;
    rdcutsvalmine[9][8]=0.3;
    rdcutsvalmine[10][8]=0.15;
    rdcutsvalmine[11][8]=0.2;
    rdcutsvalmine[12][8]=100.;
    rdcutsvalmine[13][8]=0.5;
    rdcutsvalmine[14][8]=0.998;
    rdcutsvalmine[15][8]=4.;

    //8-10
    rdcutsvalmine[0][9]=0.055;
    rdcutsvalmine[1][9]=0.04;
    rdcutsvalmine[2][9]=1.0;
    rdcutsvalmine[3][9]=1.0;
    rdcutsvalmine[4][9]=1.0;
    rdcutsvalmine[5][9]=0.15;
    rdcutsvalmine[6][9]=0.15;
    rdcutsvalmine[7][9]=-0.00006;
    rdcutsvalmine[8][9]=0.98;
    rdcutsvalmine[9][9]=0.3;
    rdcutsvalmine[10][9]=0.15;
    rdcutsvalmine[11][9]=0.34;
    rdcutsvalmine[12][9]=100.;
    rdcutsvalmine[13][9]=0.5;
    rdcutsvalmine[14][9]=0.998;
    rdcutsvalmine[15][9]=4.;

    //10-12     230312:
    rdcutsvalmine[0][10]=0.055;
    rdcutsvalmine[1][10]=0.03;
    rdcutsvalmine[2][10]=1.0;
    rdcutsvalmine[3][10]=1.0;
    rdcutsvalmine[4][10]=1.0;
    rdcutsvalmine[5][10]=0.15;
    rdcutsvalmine[6][10]=0.15;
    rdcutsvalmine[7][10]=-0.00004;
    rdcutsvalmine[8][10]=0.99;
    rdcutsvalmine[9][10]=0.3;
    rdcutsvalmine[10][10]=0.15;
    rdcutsvalmine[11][10]=0.35;
    rdcutsvalmine[12][10]=100.;
    rdcutsvalmine[13][10]=0.5;
    rdcutsvalmine[14][10]=0.998;
    rdcutsvalmine[15][10]=4.;


    //12-16 (ptbin11)
    rdcutsvalmine[0][11]=0.074;
    rdcutsvalmine[1][11]=0.025;
    rdcutsvalmine[2][11]=0.8;
    rdcutsvalmine[3][11]=1.0;
    rdcutsvalmine[4][11]=1.0;
    rdcutsvalmine[5][11]=0.2;
    rdcutsvalmine[6][11]=0.2;
    rdcutsvalmine[7][11]=-0.000023;
    rdcutsvalmine[8][11]=0.97;
    rdcutsvalmine[9][11]=0.3;
    rdcutsvalmine[10][11]=0.1;
    rdcutsvalmine[11][11]=0.45;
    rdcutsvalmine[12][11]=100.;
    rdcutsvalmine[13][11]=1;
    rdcutsvalmine[14][11]=0.998;
    rdcutsvalmine[15][11]=3.;



    //16-20 (ptbin12)
    rdcutsvalmine[0][12]=0.074;
    rdcutsvalmine[1][12]=0.04;		//0.8;
    rdcutsvalmine[2][12]=1.0;
    rdcutsvalmine[3][12]=0.30;
    rdcutsvalmine[4][12]=0.30;
    rdcutsvalmine[5][12]=0.50;
    rdcutsvalmine[6][12]=0.50;
    rdcutsvalmine[7][12]=0.0;		//0.0033;
    rdcutsvalmine[8][12]=0.8;		//0.7;
    rdcutsvalmine[9][12]=0.1;
    rdcutsvalmine[10][12]=0.1;
    rdcutsvalmine[11][12]=0.05;
    rdcutsvalmine[12][12]=100.;
    rdcutsvalmine[13][12]=1;
    rdcutsvalmine[14][12]=0.98;		//0.9;
    rdcutsvalmine[15][12]=1.;		//0.0;

    //20-24 (ptbin13)
    rdcutsvalmine[0][13]=0.074;
    rdcutsvalmine[1][13]=0.04;		//0.8;
    rdcutsvalmine[2][13]=1.0;
    rdcutsvalmine[3][13]=0.30;
    rdcutsvalmine[4][13]=0.30;
    rdcutsvalmine[5][13]=0.50;
    rdcutsvalmine[6][13]=0.50;
    rdcutsvalmine[7][13]=0.0001;	//0.0033;
    rdcutsvalmine[8][13]=0.75;		//0.7;
    rdcutsvalmine[9][13]=0.1;
    rdcutsvalmine[10][13]=0.1;
    rdcutsvalmine[11][13]=0.05;
    rdcutsvalmine[12][13]=100.;
    rdcutsvalmine[13][13]=1;
    rdcutsvalmine[14][13]=0.95;		//0.9;
    rdcutsvalmine[15][13]=0.5;		//0.0

    //>25 (ptbin14)
    rdcutsvalmine[0][14]=0.074;
    rdcutsvalmine[1][14]=0.8;
    rdcutsvalmine[2][14]=1.0;
    rdcutsvalmine[3][14]=0.30;
    rdcutsvalmine[4][14]=0.30;
    rdcutsvalmine[5][14]=0.50;
    rdcutsvalmine[6][14]=0.50;
    rdcutsvalmine[7][14]=0.0033;
    rdcutsvalmine[8][14]=0.7;
    rdcutsvalmine[9][14]=0.1;
    rdcutsvalmine[10][14]=0.1;
    rdcutsvalmine[11][14]=0.05;
    rdcutsvalmine[12][14]=100.;
    rdcutsvalmine[13][14]=1;
    rdcutsvalmine[14][14]=0.9;
    rdcutsvalmine[15][14]=0.0;
*/
/*


// Cuts 07A (corresponds to 100412_train17.root 07, except to to 3 GeV/c)
    //0-0.5
    rdcutsvalmine[0][0]=0.032;
    rdcutsvalmine[1][0]=0.02;
    rdcutsvalmine[2][0]=0.8;
    rdcutsvalmine[3][0]=1.333;
    rdcutsvalmine[4][0]=1.0;
    rdcutsvalmine[5][0]=0.07;
    rdcutsvalmine[6][0]=0.07;
    rdcutsvalmine[7][0]=-0.000362;
    rdcutsvalmine[8][0]=0.99;
    rdcutsvalmine[9][0]=0.3;
    rdcutsvalmine[10][0]=0.15;
    rdcutsvalmine[11][0]=0.05;
    rdcutsvalmine[12][0]=0.5;
    rdcutsvalmine[13][0]=0.5;
    rdcutsvalmine[14][0]=0.998;
    rdcutsvalmine[15][0]=7.0;

     //0.5-1
    rdcutsvalmine[0][1]=0.032;
    rdcutsvalmine[1][1]=0.02;
    rdcutsvalmine[2][1]=0.8;
    rdcutsvalmine[3][1]=1.333;
    rdcutsvalmine[4][1]=1.0;
    rdcutsvalmine[5][1]=0.07;
    rdcutsvalmine[6][1]=0.07;
    rdcutsvalmine[7][1]=-0.000362;
    rdcutsvalmine[8][1]=0.99;
    rdcutsvalmine[9][1]=0.3;
    rdcutsvalmine[10][1]=0.15;
    rdcutsvalmine[11][1]=0.05;
    rdcutsvalmine[12][1]=0.5;
    rdcutsvalmine[13][1]=0.5;
    rdcutsvalmine[14][1]=0.998;
    rdcutsvalmine[15][1]=7.0;

    //1-2
    rdcutsvalmine[0][2]=0.032;
    rdcutsvalmine[1][2]=0.02;
    rdcutsvalmine[2][2]=0.8;
    rdcutsvalmine[3][2]=1.333;
    rdcutsvalmine[4][2]=1.0;
    rdcutsvalmine[5][2]=0.07;
    rdcutsvalmine[6][2]=0.07;
    rdcutsvalmine[7][2]=-0.000362;
    rdcutsvalmine[8][2]=0.99;
    rdcutsvalmine[9][2]=0.3;
    rdcutsvalmine[10][2]=0.15;
    rdcutsvalmine[11][2]=0.05;
    rdcutsvalmine[12][2]=0.5;
    rdcutsvalmine[13][2]=0.5;
    rdcutsvalmine[14][2]=0.998;
    rdcutsvalmine[15][2]=7.0;
    //2-3
    rdcutsvalmine[0][3]=0.032;
    rdcutsvalmine[1][3]=0.02;
    rdcutsvalmine[2][3]=0.8;
    rdcutsvalmine[3][3]=1.333;
    rdcutsvalmine[4][3]=1.0;
    rdcutsvalmine[5][3]=0.07;
    rdcutsvalmine[6][3]=0.07;
    rdcutsvalmine[7][3]=-0.000362;
    rdcutsvalmine[8][3]=0.99;
    rdcutsvalmine[9][3]=0.3;
    rdcutsvalmine[10][3]=0.15;
    rdcutsvalmine[11][3]=0.05;
    rdcutsvalmine[12][3]=0.5;
    rdcutsvalmine[13][3]=0.5;
    rdcutsvalmine[14][3]=0.998;
    rdcutsvalmine[15][3]=7.0;

    //3-4
    rdcutsvalmine[0][4]=0.032;
    rdcutsvalmine[1][4]=0.02;
    rdcutsvalmine[2][4]=0.8;
    rdcutsvalmine[3][4]=1.333;
    rdcutsvalmine[4][4]=1.0;
    rdcutsvalmine[5][4]=0.07;
    rdcutsvalmine[6][4]=0.07;
    rdcutsvalmine[7][4]=-0.000362;
    rdcutsvalmine[8][4]=0.99;
    rdcutsvalmine[9][4]=0.3;
    rdcutsvalmine[10][4]=0.15;
    rdcutsvalmine[11][4]=0.05;
    rdcutsvalmine[12][4]=0.5;
    rdcutsvalmine[13][4]=0.5;
    rdcutsvalmine[14][4]=0.998;
    rdcutsvalmine[15][4]=7.0;

    //4-5
    rdcutsvalmine[0][5]=0.032;
    rdcutsvalmine[1][5]=0.025;
    rdcutsvalmine[2][5]=0.8;
    rdcutsvalmine[3][5]=1.;
    rdcutsvalmine[4][5]=1.;
    rdcutsvalmine[5][5]=0.075;
    rdcutsvalmine[6][5]=0.075;
    rdcutsvalmine[7][5]=-0.000285;
    rdcutsvalmine[8][5]=0.99;
    rdcutsvalmine[9][5]=0.15;
    rdcutsvalmine[10][5]=0.15;
    rdcutsvalmine[11][5]=0.05;
    rdcutsvalmine[12][5]=10;
    rdcutsvalmine[13][5]=0.5;
    rdcutsvalmine[14][5]=0.998;
    rdcutsvalmine[15][5]=6.;

    //5-6
    rdcutsvalmine[0][6]=0.04;
    rdcutsvalmine[1][6]=0.03;
    rdcutsvalmine[2][6]=1.;
    rdcutsvalmine[3][6]=1.;
    rdcutsvalmine[4][6]=1.;
    rdcutsvalmine[5][6]=0.075;
    rdcutsvalmine[6][6]=0.075;
    rdcutsvalmine[7][6]=-0.00028;
    rdcutsvalmine[8][6]=0.99;
    rdcutsvalmine[9][6]=0.15;
    rdcutsvalmine[10][6]=0.15;
    rdcutsvalmine[11][6]=0.2;
    rdcutsvalmine[12][6]=10.;
    rdcutsvalmine[13][6]=0.5;
    rdcutsvalmine[14][6]=0.998;
    rdcutsvalmine[15][6]=5.;

    //6-7
    rdcutsvalmine[0][7]=0.036;
    rdcutsvalmine[1][7]=0.03750;
    rdcutsvalmine[2][7]=1.0;
    rdcutsvalmine[3][7]=1.;
    rdcutsvalmine[4][7]=1.;
    rdcutsvalmine[5][7]=0.1;
    rdcutsvalmine[6][7]=0.1;
    rdcutsvalmine[7][7]=-0.00022;
    rdcutsvalmine[8][7]=0.985;
    rdcutsvalmine[9][7]=0.3;
    rdcutsvalmine[10][7]=0.1;
    rdcutsvalmine[11][7]=0.2;
    rdcutsvalmine[12][7]=100.;
    rdcutsvalmine[13][7]=0.5;
    rdcutsvalmine[14][7]=0.998;
    rdcutsvalmine[15][7]=4.0;

    //7-8
    rdcutsvalmine[0][8]=0.036;
    rdcutsvalmine[1][8]=0.03750;
    rdcutsvalmine[2][8]=1.0;
    rdcutsvalmine[3][8]=1.;
    rdcutsvalmine[4][8]=1.;
    rdcutsvalmine[5][8]=0.12;
    rdcutsvalmine[6][8]=0.12;
    rdcutsvalmine[7][8]=-0.000165;
    rdcutsvalmine[8][8]=0.98;
    rdcutsvalmine[9][8]=0.3;
    rdcutsvalmine[10][8]=0.15;
    rdcutsvalmine[11][8]=0.2;
    rdcutsvalmine[12][8]=100.;
    rdcutsvalmine[13][8]=0.5;
    rdcutsvalmine[14][8]=0.998;
    rdcutsvalmine[15][8]=4.;

    //8-10		230312: changed to cut of bin10 to increase signal
    rdcutsvalmine[0][9]=0.055;
    rdcutsvalmine[1][9]=0.04;
    rdcutsvalmine[2][9]=1.0;
    rdcutsvalmine[3][9]=1.0;
    rdcutsvalmine[4][9]=1.0;
    rdcutsvalmine[5][9]=0.15;
    rdcutsvalmine[6][9]=0.15;
    rdcutsvalmine[7][9]=-0.00006;
    rdcutsvalmine[8][9]=0.98;
    rdcutsvalmine[9][9]=0.3;
    rdcutsvalmine[10][9]=0.15;
    rdcutsvalmine[11][9]=0.34;
    rdcutsvalmine[12][9]=100.;
    rdcutsvalmine[13][9]=0.5;
    rdcutsvalmine[14][9]=0.998;
    rdcutsvalmine[15][9]=4.;

    //10-12     230312:
    rdcutsvalmine[0][10]=0.055;
    rdcutsvalmine[1][10]=0.03;
    rdcutsvalmine[2][10]=1.0;
    rdcutsvalmine[3][10]=1.0;
    rdcutsvalmine[4][10]=1.0;
    rdcutsvalmine[5][10]=0.15;
    rdcutsvalmine[6][10]=0.15;
    rdcutsvalmine[7][10]=-0.00004;
    rdcutsvalmine[8][10]=0.99;
    rdcutsvalmine[9][10]=0.3;
    rdcutsvalmine[10][10]=0.15;
    rdcutsvalmine[11][10]=0.35;
    rdcutsvalmine[12][10]=100.;
    rdcutsvalmine[13][10]=0.5;
    rdcutsvalmine[14][10]=0.998;
    rdcutsvalmine[15][10]=4.;


    //12-16 (ptbin11)
    rdcutsvalmine[0][11]=0.074;
    rdcutsvalmine[1][11]=0.025;
    rdcutsvalmine[2][11]=0.8;
    rdcutsvalmine[3][11]=1.0;
    rdcutsvalmine[4][11]=1.0;
    rdcutsvalmine[5][11]=0.2;
    rdcutsvalmine[6][11]=0.2;
    rdcutsvalmine[7][11]=-0.000023;
    rdcutsvalmine[8][11]=0.97;
    rdcutsvalmine[9][11]=0.3;
    rdcutsvalmine[10][11]=0.1;
    rdcutsvalmine[11][11]=0.45;
    rdcutsvalmine[12][11]=100.;
    rdcutsvalmine[13][11]=1;
    rdcutsvalmine[14][11]=0.998;
    rdcutsvalmine[15][11]=3.;



    //16-20 (ptbin12)
    rdcutsvalmine[0][12]=0.074;
    rdcutsvalmine[1][12]=0.8;
    rdcutsvalmine[2][12]=1.0;
    rdcutsvalmine[3][12]=0.30;
    rdcutsvalmine[4][12]=0.30;
    rdcutsvalmine[5][12]=0.50;
    rdcutsvalmine[6][12]=0.50;
    rdcutsvalmine[7][12]=0.0033;
    rdcutsvalmine[8][12]=0.7;
    rdcutsvalmine[9][12]=0.1;
    rdcutsvalmine[10][12]=0.1;
    rdcutsvalmine[11][12]=0.05;
    rdcutsvalmine[12][12]=100.;
    rdcutsvalmine[13][12]=1;
    rdcutsvalmine[14][12]=0.9;
    rdcutsvalmine[15][12]=0.0;

    //20-24 (ptbin13)
    rdcutsvalmine[0][13]=0.074;
    rdcutsvalmine[1][13]=0.8;
    rdcutsvalmine[2][13]=1.0;
    rdcutsvalmine[3][13]=0.30;
    rdcutsvalmine[4][13]=0.30;
    rdcutsvalmine[5][13]=0.50;
    rdcutsvalmine[6][13]=0.50;
    rdcutsvalmine[7][13]=0.0033;
    rdcutsvalmine[8][13]=0.7;
    rdcutsvalmine[9][13]=0.1;
    rdcutsvalmine[10][13]=0.1;
    rdcutsvalmine[11][13]=0.05;
    rdcutsvalmine[12][13]=100.;
    rdcutsvalmine[13][13]=1;
    rdcutsvalmine[14][13]=0.9;
    rdcutsvalmine[15][13]=0.0;

    //>25 (ptbin14)
    rdcutsvalmine[0][14]=0.074;
    rdcutsvalmine[1][14]=0.8;
    rdcutsvalmine[2][14]=1.0;
    rdcutsvalmine[3][14]=0.30;
    rdcutsvalmine[4][14]=0.30;
    rdcutsvalmine[5][14]=0.50;
    rdcutsvalmine[6][14]=0.50;
    rdcutsvalmine[7][14]=0.0033;
    rdcutsvalmine[8][14]=0.7;
    rdcutsvalmine[9][14]=0.1;
    rdcutsvalmine[10][14]=0.1;
    rdcutsvalmine[11][14]=0.05;
    rdcutsvalmine[12][14]=100.;
    rdcutsvalmine[13][14]=1;
    rdcutsvalmine[14][14]=0.9;
    rdcutsvalmine[15][14]=0.0;
*/





//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^







/*
// Cuts 07A with everything relaxed (corresponds to 100412_train17.root 07B)
    //0-0.5
    rdcutsvalmine[0][0]=0.026;
    rdcutsvalmine[1][0]=0.022;
    rdcutsvalmine[2][0]=0.7;
    rdcutsvalmine[3][0]=0.21;
    rdcutsvalmine[4][0]=0.21;
    rdcutsvalmine[5][0]=0.05;
    rdcutsvalmine[6][0]=0.05;
    rdcutsvalmine[7][0]=0.0004;
    rdcutsvalmine[8][0]=0.6;
    rdcutsvalmine[9][0]=0.3;
    rdcutsvalmine[10][0]=0.1;
    rdcutsvalmine[11][0]=0.05;
    rdcutsvalmine[12][0]=100.;
    rdcutsvalmine[13][0]=0.5;
    rdcutsvalmine[14][0]=-1.;
    rdcutsvalmine[15][0]=0.;

    //0.5-1
    rdcutsvalmine[0][1]=0.026;
    rdcutsvalmine[1][1]=0.04;
    rdcutsvalmine[2][1]=0.8;
    rdcutsvalmine[3][1]=0.5;
    rdcutsvalmine[4][1]=0.5;
    rdcutsvalmine[5][1]=0.08;
    rdcutsvalmine[6][1]=0.08;
    rdcutsvalmine[7][1]=0.0003;
    rdcutsvalmine[8][1]=0.65;
    rdcutsvalmine[9][1]=0.3;
    rdcutsvalmine[10][1]=0.1;
    rdcutsvalmine[11][1]=0.05;
    rdcutsvalmine[12][1]=100.;
    rdcutsvalmine[13][1]=0.5;
    rdcutsvalmine[14][1]=-1.;
    rdcutsvalmine[15][1]=0.;

    //1-2
    rdcutsvalmine[0][2]=0.022;
    rdcutsvalmine[1][2]=0.04;
    rdcutsvalmine[2][2]=0.8;
    rdcutsvalmine[3][2]=0.7;
    rdcutsvalmine[4][2]=0.7;
    rdcutsvalmine[5][2]=0.08;
    rdcutsvalmine[6][2]=0.08;
    rdcutsvalmine[7][2]=0.00014;
    rdcutsvalmine[8][2]=0.75;
    rdcutsvalmine[9][2]=0.3;
    rdcutsvalmine[10][2]=0.1;
    rdcutsvalmine[11][2]=0.05;
    rdcutsvalmine[12][2]=100.;
    rdcutsvalmine[13][2]=0.5;
    rdcutsvalmine[14][2]=-1.;
    rdcutsvalmine[15][2]=0.;
    //2-3
    rdcutsvalmine[0][3]=0.024;
    rdcutsvalmine[1][3]=0.025;
    rdcutsvalmine[2][3]=0.8;
    rdcutsvalmine[3][3]=0.7;
    rdcutsvalmine[4][3]=0.7;
    rdcutsvalmine[5][3]=0.08;
    rdcutsvalmine[6][3]=0.08;
    rdcutsvalmine[7][3]=0.00014;
    rdcutsvalmine[8][3]=0.85;
    rdcutsvalmine[9][3]=0.3;
    rdcutsvalmine[10][3]=0.1;
    rdcutsvalmine[11][3]=0.05;
    rdcutsvalmine[12][3]=100.;
    rdcutsvalmine[13][3]=0.5;
    rdcutsvalmine[14][3]=-1.;
    rdcutsvalmine[15][3]=0.;

    //3-4
    rdcutsvalmine[0][4]=0.032;
    rdcutsvalmine[1][4]=0.025;
    rdcutsvalmine[2][4]=1.0;
    rdcutsvalmine[3][4]=1.;
    rdcutsvalmine[4][4]=1.;
    rdcutsvalmine[5][4]=0.1;
    rdcutsvalmine[6][4]=0.1;
    rdcutsvalmine[7][4]=-0.000362;
    rdcutsvalmine[8][4]=0.99;
    rdcutsvalmine[9][4]=0.3;
    rdcutsvalmine[10][4]=0.15;
    rdcutsvalmine[11][4]=0.05;
    rdcutsvalmine[12][4]=0.5;
    rdcutsvalmine[13][4]=0.5;
    rdcutsvalmine[14][4]=0.998;
    rdcutsvalmine[15][4]=6.0;

    //4-5
    rdcutsvalmine[0][5]=0.032;
    rdcutsvalmine[1][5]=0.025; //0.1
    rdcutsvalmine[2][5]=1.0;
    rdcutsvalmine[3][5]=1.;
    rdcutsvalmine[4][5]=1.;
    rdcutsvalmine[5][5]=0.1;
    rdcutsvalmine[6][5]=0.1;
    rdcutsvalmine[7][5]=-0.000285;
    rdcutsvalmine[8][5]=0.99;
    rdcutsvalmine[9][5]=0.15;
    rdcutsvalmine[10][5]=0.15;
    rdcutsvalmine[11][5]=0.05;
    rdcutsvalmine[12][5]=10;
    rdcutsvalmine[13][5]=0.5;
    rdcutsvalmine[14][5]=0.998;
    rdcutsvalmine[15][5]=5.;

    //5-6
    rdcutsvalmine[0][6]=0.04;
    rdcutsvalmine[1][6]=0.03;
    rdcutsvalmine[2][6]=1.;
    rdcutsvalmine[3][6]=.8;
    rdcutsvalmine[4][6]=.8;
    rdcutsvalmine[5][6]=0.1;
    rdcutsvalmine[6][6]=0.1;
    rdcutsvalmine[7][6]=-0.00028;
    rdcutsvalmine[8][6]=0.99;
    rdcutsvalmine[9][6]=0.15;
    rdcutsvalmine[10][6]=0.15;
    rdcutsvalmine[11][6]=0.2;
    rdcutsvalmine[12][6]=10.;
    rdcutsvalmine[13][6]=0.5;
    rdcutsvalmine[14][6]=0.998;
    rdcutsvalmine[15][6]=4.;



    //6-7
    rdcutsvalmine[0][7]=0.036;
    rdcutsvalmine[1][7]=0.04;
    rdcutsvalmine[2][7]=1.0;
    rdcutsvalmine[3][7]=0.8;
    rdcutsvalmine[4][7]=0.8;
    rdcutsvalmine[5][7]=0.1;
    rdcutsvalmine[6][7]=0.1;
    rdcutsvalmine[7][7]=-0.00022;
    rdcutsvalmine[8][7]=0.985;
    rdcutsvalmine[9][7]=0.3;
    rdcutsvalmine[10][7]=0.1;
    rdcutsvalmine[11][7]=0.2;
    rdcutsvalmine[12][7]=100.;
    rdcutsvalmine[13][7]=0.5;
    rdcutsvalmine[14][7]=0.998;
    rdcutsvalmine[15][7]=3.0;

    //7-8
    rdcutsvalmine[0][8]=0.036;
    rdcutsvalmine[1][8]=0.04;
    rdcutsvalmine[2][8]=1.0;
    rdcutsvalmine[3][8]=0.8;
    rdcutsvalmine[4][8]=0.8;
    rdcutsvalmine[5][8]=0.12;
    rdcutsvalmine[6][8]=0.12;
    rdcutsvalmine[7][8]=-0.00012;
    rdcutsvalmine[8][8]=0.98;
    rdcutsvalmine[9][8]=0.3;
    rdcutsvalmine[10][8]=0.15;
    rdcutsvalmine[11][8]=0.2;
    rdcutsvalmine[12][8]=100.;
    rdcutsvalmine[13][8]=0.5;
    rdcutsvalmine[14][8]=0.998;
    rdcutsvalmine[15][8]=3.;

    //8-10		230312: changed to cut of bin10 to increase signal
    rdcutsvalmine[0][9]=0.055;
    rdcutsvalmine[1][9]=0.04;
    rdcutsvalmine[2][9]=1.0;
    rdcutsvalmine[3][9]=0.8;			//0.9;
    rdcutsvalmine[4][9]=0.8;			//0.9;
    rdcutsvalmine[5][9]=0.15;
    rdcutsvalmine[6][9]=0.15;
    rdcutsvalmine[7][9]=-0.00004;		//-0.00008;
    rdcutsvalmine[8][9]=0.98;
    rdcutsvalmine[9][9]=0.3;
    rdcutsvalmine[10][9]=0.15;
    rdcutsvalmine[11][9]=0.34;
    rdcutsvalmine[12][9]=100.;
    rdcutsvalmine[13][9]=0.5;
    rdcutsvalmine[14][9]=0.998;
    rdcutsvalmine[15][9]=3.;			//4.8;

    //10-12     230312:
    rdcutsvalmine[0][10]=0.055;
    rdcutsvalmine[1][10]=0.03;
    rdcutsvalmine[2][10]=1.0;
    rdcutsvalmine[3][10]=0.7;
    rdcutsvalmine[4][10]=0.7;
    rdcutsvalmine[5][10]=0.15;
    rdcutsvalmine[6][10]=0.15;
    rdcutsvalmine[7][10]=-0.00002;
    rdcutsvalmine[8][10]=0.99;
    rdcutsvalmine[9][10]=0.3;
    rdcutsvalmine[10][10]=0.15;
    rdcutsvalmine[11][10]=0.35;
    rdcutsvalmine[12][10]=100.;
    rdcutsvalmine[13][10]=0.5;
    rdcutsvalmine[14][10]=0.998;
    rdcutsvalmine[15][10]=3.;


    //12-16 (ptbin11)
    rdcutsvalmine[0][11]=0.074;
    rdcutsvalmine[1][11]=0.025;
    rdcutsvalmine[2][11]=0.8;
    rdcutsvalmine[3][11]=0.7;
    rdcutsvalmine[4][11]=0.7;
    rdcutsvalmine[5][11]=0.2;
    rdcutsvalmine[6][11]=0.2;
    rdcutsvalmine[7][11]=-0.000015;
    rdcutsvalmine[8][11]=0.97;
    rdcutsvalmine[9][11]=0.3;
    rdcutsvalmine[10][11]=0.1;
    rdcutsvalmine[11][11]=0.45;
    rdcutsvalmine[12][11]=100.;
    rdcutsvalmine[13][11]=1;
    rdcutsvalmine[14][11]=0.998;
    rdcutsvalmine[15][11]=2.;//5.;



    //16-20 (ptbin12)
    rdcutsvalmine[0][12]=0.074;
    rdcutsvalmine[1][12]=0.8;
    rdcutsvalmine[2][12]=1.0;
    rdcutsvalmine[3][12]=0.30;
    rdcutsvalmine[4][12]=0.30;
    rdcutsvalmine[5][12]=0.50;
    rdcutsvalmine[6][12]=0.50;
    rdcutsvalmine[7][12]=0.0033;
    rdcutsvalmine[8][12]=0.7;
    rdcutsvalmine[9][12]=0.1;
    rdcutsvalmine[10][12]=0.1;
    rdcutsvalmine[11][12]=0.05;
    rdcutsvalmine[12][12]=100.;
    rdcutsvalmine[13][12]=1;
    rdcutsvalmine[14][12]=0.9;
    rdcutsvalmine[15][12]=0.0;

    //20-24 (ptbin13)
    rdcutsvalmine[0][13]=0.074;
    rdcutsvalmine[1][13]=0.8;
    rdcutsvalmine[2][13]=1.0;
    rdcutsvalmine[3][13]=0.30;
    rdcutsvalmine[4][13]=0.30;
    rdcutsvalmine[5][13]=0.50;
    rdcutsvalmine[6][13]=0.50;
    rdcutsvalmine[7][13]=0.0033;
    rdcutsvalmine[8][13]=0.7;
    rdcutsvalmine[9][13]=0.1;
    rdcutsvalmine[10][13]=0.1;
    rdcutsvalmine[11][13]=0.05;
    rdcutsvalmine[12][13]=100.;
    rdcutsvalmine[13][13]=1;
    rdcutsvalmine[14][13]=0.9;
    rdcutsvalmine[15][13]=0.0;

    //>25 (ptbin14)
    rdcutsvalmine[0][14]=0.074;
    rdcutsvalmine[1][14]=0.8;
    rdcutsvalmine[2][14]=1.0;
    rdcutsvalmine[3][14]=0.30;
    rdcutsvalmine[4][14]=0.30;
    rdcutsvalmine[5][14]=0.50;
    rdcutsvalmine[6][14]=0.50;
    rdcutsvalmine[7][14]=0.0033;
    rdcutsvalmine[8][14]=0.7;
    rdcutsvalmine[9][14]=0.1;
    rdcutsvalmine[10][14]=0.1;
    rdcutsvalmine[11][14]=0.05;
    rdcutsvalmine[12][14]=100.;
    rdcutsvalmine[13][14]=1;
    rdcutsvalmine[14][14]=0.9;
    rdcutsvalmine[15][14]=0.0;
*/

/*
// Cuts 21_03_12 with further relaxiations
    //0-0.5
    rdcutsvalmine[0][0]=0.026;
    rdcutsvalmine[1][0]=0.022;
    rdcutsvalmine[2][0]=0.7;
    rdcutsvalmine[3][0]=0.21;
    rdcutsvalmine[4][0]=0.21;
    rdcutsvalmine[5][0]=0.05;
    rdcutsvalmine[6][0]=0.05;
    rdcutsvalmine[7][0]=0.0004;
    rdcutsvalmine[8][0]=0.6;
    rdcutsvalmine[9][0]=0.3;
    rdcutsvalmine[10][0]=0.1;
    rdcutsvalmine[11][0]=0.05;
    rdcutsvalmine[12][0]=100.;
    rdcutsvalmine[13][0]=0.5;
    rdcutsvalmine[14][0]=-1.;
    rdcutsvalmine[15][0]=0.;

    //0.5-1
    rdcutsvalmine[0][1]=0.026;
    rdcutsvalmine[1][1]=0.04;
    rdcutsvalmine[2][1]=0.8;
    rdcutsvalmine[3][1]=0.5;
    rdcutsvalmine[4][1]=0.5;
    rdcutsvalmine[5][1]=0.08;
    rdcutsvalmine[6][1]=0.08;
    rdcutsvalmine[7][1]=0.0003;
    rdcutsvalmine[8][1]=0.65;
    rdcutsvalmine[9][1]=0.3;
    rdcutsvalmine[10][1]=0.1;
    rdcutsvalmine[11][1]=0.05;
    rdcutsvalmine[12][1]=100.;
    rdcutsvalmine[13][1]=0.5;
    rdcutsvalmine[14][1]=-1.;
    rdcutsvalmine[15][1]=0.;

    //1-2
    rdcutsvalmine[0][2]=0.022;
    rdcutsvalmine[1][2]=0.04;
    rdcutsvalmine[2][2]=0.8;
    rdcutsvalmine[3][2]=0.7;
    rdcutsvalmine[4][2]=0.7;
    rdcutsvalmine[5][2]=0.08;
    rdcutsvalmine[6][2]=0.08;
    rdcutsvalmine[7][2]=0.00014;
    rdcutsvalmine[8][2]=0.75;
    rdcutsvalmine[9][2]=0.3;
    rdcutsvalmine[10][2]=0.1;
    rdcutsvalmine[11][2]=0.05;
    rdcutsvalmine[12][2]=100.;
    rdcutsvalmine[13][2]=0.5;
    rdcutsvalmine[14][2]=-1.;
    rdcutsvalmine[15][2]=0.;
    //2-3
    rdcutsvalmine[0][3]=0.024;
    rdcutsvalmine[1][3]=0.025;
    rdcutsvalmine[2][3]=0.8;
    rdcutsvalmine[3][3]=0.7;
    rdcutsvalmine[4][3]=0.7;
    rdcutsvalmine[5][3]=0.08;
    rdcutsvalmine[6][3]=0.08;
    rdcutsvalmine[7][3]=0.00014;
    rdcutsvalmine[8][3]=0.85;
    rdcutsvalmine[9][3]=0.3;
    rdcutsvalmine[10][3]=0.1;
    rdcutsvalmine[11][3]=0.05;
    rdcutsvalmine[12][3]=100.;
    rdcutsvalmine[13][3]=0.5;
    rdcutsvalmine[14][3]=-1.;
    rdcutsvalmine[15][3]=0.;

    //3-4
  //rdcutsvalmine[0][4]=0.032;
    //rdcutsvalmine[1][4]=0.02;//0.01250;		//0.1;//0.025;
    //rdcutsvalmine[2][4]=1.0;			//1.0;//0.7;
    //rdcutsvalmine[3][4]=1.2;//1.333;		//0.5;//1.1;
    //rdcutsvalmine[4][4]=1.0;			//0.5;//1.1;
    //rdcutsvalmine[5][4]=0.1;			//0.1;//0.08;
    //rdcutsvalmine[6][4]=0.1;//0.05833;		//0.1;//0.08;
    //rdcutsvalmine[7][4]=-0.000285;		//-0.00017;
    //rdcutsvalmine[8][4]=0.99;			//0.98;
    //rdcutsvalmine[9][4]=0.3;
    //rdcutsvalmine[10][4]=0.15;
    //rdcutsvalmine[11][4]=0.05;
    //rdcutsvalmine[12][4]=0.5;
    //rdcutsvalmine[13][4]=0.5;
    //rdcutsvalmine[14][4]=0.998;			//0.99;
    //rdcutsvalmine[15][4]=7.0;			//5.;
  //3-4  			Gives high background
    rdcutsvalmine[0][4]=0.032;
    rdcutsvalmine[1][4]=0.02;						//0.1;//0.025;
    rdcutsvalmine[2][4]=1.0;						//1.0;//0.7;
    rdcutsvalmine[3][4]=1.333;						//0.5;//1.1;
    rdcutsvalmine[4][4]=1.0;						//0.5;//1.1;
    rdcutsvalmine[5][4]=0.1;						//0.1;//0.08;
    rdcutsvalmine[6][4]=0.1;						//0.1;//0.08;
    rdcutsvalmine[7][4]=-0.000285;					//-0.00017;
    rdcutsvalmine[8][4]=0.98;						//0.98;
    rdcutsvalmine[9][4]=0.3;
    rdcutsvalmine[10][4]=0.15;
    rdcutsvalmine[11][4]=0.05;
    rdcutsvalmine[12][4]=0.5;
    rdcutsvalmine[13][4]=0.5;
    rdcutsvalmine[14][4]=0.994;			//0.99;
    rdcutsvalmine[15][4]=5.0;			//5.;


    //4-5
    //rdcutsvalmine[0][5]=0.032;
    //rdcutsvalmine[1][5]=0.025;
    //rdcutsvalmine[2][5]=0.7;
    //rdcutsvalmine[3][5]=1.1;//1.2;
    //rdcutsvalmine[4][5]=1.1;//1.2;
    //rdcutsvalmine[5][5]=0.1;//0.09;
    //rdcutsvalmine[6][5]=0.1;//0.09;
    //rdcutsvalmine[7][5]=-0.000337;
    //rdcutsvalmine[8][5]=0.992;
    //rdcutsvalmine[9][5]=0.15;
    //rdcutsvalmine[10][5]=0.15;
    //rdcutsvalmine[11][5]=0.05;
    //rdcutsvalmine[12][5]=10;
    //rdcutsvalmine[13][5]=0.5;
    //rdcutsvalmine[14][5]=0.998;//0.996;
    //rdcutsvalmine[15][5]=8.0;//8.5;

    //5-6
    //rdcutsvalmine[0][6]=0.04;
    //rdcutsvalmine[1][6]=0.027;
    //rdcutsvalmine[2][6]=0.8;
    //rdcutsvalmine[3][6]=1.1;
    //rdcutsvalmine[4][6]=1.1;
    //rdcutsvalmine[5][6]=0.1;
    //rdcutsvalmine[6][6]=0.1;
    //rdcutsvalmine[7][6]=-0.0003;
    //rdcutsvalmine[8][6]=0.99;
    //rdcutsvalmine[9][6]=0.15;
    //rdcutsvalmine[10][6]=0.15;
    //rdcutsvalmine[11][6]=0.2;
    //rdcutsvalmine[12][6]=10.;
    //rdcutsvalmine[13][6]=0.5;
    //rdcutsvalmine[14][6]=0.998;//0.995;
    //rdcutsvalmine[15][6]=5.0;//6.0;

    //4-5
    rdcutsvalmine[0][5]=0.032;
    rdcutsvalmine[1][5]=0.0375; //0.1
    rdcutsvalmine[2][5]=0.5625;
    rdcutsvalmine[3][5]=1.;
    rdcutsvalmine[4][5]=1.;
    rdcutsvalmine[5][5]=0.075;
    rdcutsvalmine[6][5]=0.0583;
    rdcutsvalmine[7][5]=-0.000362;
    rdcutsvalmine[8][5]=0.98;
    rdcutsvalmine[9][5]=0.15;
    rdcutsvalmine[10][5]=0.15;
    rdcutsvalmine[11][5]=0.05;
    rdcutsvalmine[12][5]=10;
    rdcutsvalmine[13][5]=0.5;
    rdcutsvalmine[14][5]=0.993;
    rdcutsvalmine[15][5]=6.;

    //5-6
    rdcutsvalmine[0][6]=0.04;
    rdcutsvalmine[1][6]=0.0375;
    rdcutsvalmine[2][6]=1.;
    rdcutsvalmine[3][6]=1.33;
    rdcutsvalmine[4][6]=1.167;
    rdcutsvalmine[5][6]=0.075;
    rdcutsvalmine[6][6]=0.075;
    rdcutsvalmine[7][6]=-0.000287;
    rdcutsvalmine[8][6]=0.98;
    rdcutsvalmine[9][6]=0.15;
    rdcutsvalmine[10][6]=0.15;
    rdcutsvalmine[11][6]=0.2;
    rdcutsvalmine[12][6]=10.;
    rdcutsvalmine[13][6]=0.5;
    rdcutsvalmine[14][6]=0.993;
    rdcutsvalmine[15][6]=6;


   //6-7
    //rdcutsvalmine[0][7]=0.036;
    //rdcutsvalmine[1][7]=0.027;//0.03750;
    //rdcutsvalmine[2][7]=1.0;
    //rdcutsvalmine[3][7]=1.2;//1.333;
    //rdcutsvalmine[4][7]=1.2;//1.333;
    //rdcutsvalmine[5][7]=0.075;
    //rdcutsvalmine[6][7]=0.075;
    //rdcutsvalmine[7][7]=-0.000243;
    //rdcutsvalmine[8][7]=0.98;//0.986;
    //rdcutsvalmine[9][7]=0.3;
    //rdcutsvalmine[10][7]=0.1;
    //rdcutsvalmine[11][7]=0.2;
    //rdcutsvalmine[12][7]=100.;
    //rdcutsvalmine[13][7]=0.5;
    //rdcutsvalmine[14][7]=0.998;//0.992;
    //rdcutsvalmine[15][7]=4.0;

    //7-8
    //rdcutsvalmine[0][8]=0.036;
    //rdcutsvalmine[1][8]=0.027;//0.03750;
    //rdcutsvalmine[2][8]=1.0;
    //rdcutsvalmine[3][8]=0.833;
    //rdcutsvalmine[4][8]=0.833;
    //rdcutsvalmine[5][8]=0.2;
    //rdcutsvalmine[6][8]=0.2;
    //rdcutsvalmine[7][8]=-0.000163;
    //rdcutsvalmine[8][8]=0.967;
    //rdcutsvalmine[9][8]=0.3;
    //rdcutsvalmine[10][8]=0.15;
    //rdcutsvalmine[11][8]=0.2;
    //rdcutsvalmine[12][8]=100.;
    //rdcutsvalmine[13][8]=0.5;
    //rdcutsvalmine[14][8]=0.995;
    //rdcutsvalmine[15][8]=4.0;//9.;

    //6-7
    rdcutsvalmine[0][7]=0.036;
    rdcutsvalmine[1][7]=0.03750;
    rdcutsvalmine[2][7]=1.0;
    rdcutsvalmine[3][7]=1.333;
    rdcutsvalmine[4][7]=1.333;
    rdcutsvalmine[5][7]=0.075;
    rdcutsvalmine[6][7]=0.075;
    rdcutsvalmine[7][7]=-0.000243;
    rdcutsvalmine[8][7]=0.98;
    rdcutsvalmine[9][7]=0.3;
    rdcutsvalmine[10][7]=0.1;
    rdcutsvalmine[11][7]=0.2;
    rdcutsvalmine[12][7]=100.;
    rdcutsvalmine[13][7]=0.5;
    rdcutsvalmine[14][7]=0.994;
    rdcutsvalmine[15][7]=4.0;

    //7-8
    rdcutsvalmine[0][8]=0.036;
    rdcutsvalmine[1][8]=0.0375;
    rdcutsvalmine[2][8]=1.;
    rdcutsvalmine[3][8]=0.833;
    rdcutsvalmine[4][8]=1.333;
    rdcutsvalmine[5][8]=0.1;
    rdcutsvalmine[6][8]=0.1;
    rdcutsvalmine[7][8]=-0.000163;
    rdcutsvalmine[8][8]=0.97;
    rdcutsvalmine[9][8]=0.3;
    rdcutsvalmine[10][8]=0.15;
    rdcutsvalmine[11][8]=0.2;
    rdcutsvalmine[12][8]=100.;
    rdcutsvalmine[13][8]=0.5;
    rdcutsvalmine[14][8]=0.994;
    rdcutsvalmine[15][8]=4.;

    //8-10
    //rdcutsvalmine[0][9]=0.055;
    //rdcutsvalmine[1][9]=0.03;
    //rdcutsvalmine[2][9]=1.0;
    //rdcutsvalmine[3][9]=0.9;
    //rdcutsvalmine[4][9]=0.9;
    //rdcutsvalmine[5][9]=0.15;
    //rdcutsvalmine[6][9]=0.15;
    //rdcutsvalmine[7][9]=-0.00008;	//-0.00001;
    //rdcutsvalmine[8][9]=0.94;//0.99;		//0.92;
    //rdcutsvalmine[9][9]=0.3;
    //rdcutsvalmine[10][9]=0.15;
    //rdcutsvalmine[11][9]=0.34;
    //rdcutsvalmine[12][9]=100.;
    //rdcutsvalmine[13][9]=0.5;
    //rdcutsvalmine[14][9]=0.998;
    //rdcutsvalmine[15][9]=4.0;//4.8;

    //10-12
    //rdcutsvalmine[0][10]=0.055;
    //rdcutsvalmine[1][10]=0.03;
    //rdcutsvalmine[2][10]=1.0;
    //rdcutsvalmine[3][10]=1.0;
    //rdcutsvalmine[4][10]=1.0;
    //rdcutsvalmine[5][10]=0.15;
    //rdcutsvalmine[6][10]=0.15;
    //rdcutsvalmine[7][10]=-0.00004;//-0.00001;
    //rdcutsvalmine[8][10]=0.99;//0.92;
    //rdcutsvalmine[9][10]=0.3;
    //rdcutsvalmine[10][10]=0.15;
    //rdcutsvalmine[11][10]=0.35;
    //rdcutsvalmine[12][10]=100.;
    //rdcutsvalmine[13][10]=0.5;
    //rdcutsvalmine[14][10]=0.998;
    //rdcutsvalmine[15][10]=4.;

   //8-10
    rdcutsvalmine[0][9]=0.055;
    rdcutsvalmine[1][9]=0.03;
    rdcutsvalmine[2][9]=1.0;
    rdcutsvalmine[3][9]=0.9;
    rdcutsvalmine[4][9]=0.9;
    rdcutsvalmine[5][9]=0.15;
    rdcutsvalmine[6][9]=0.15;
    rdcutsvalmine[7][9]=-0.00008;//-0.00001;
    rdcutsvalmine[8][9]=0.95;//0.92;
    rdcutsvalmine[9][9]=0.3;
    rdcutsvalmine[10][9]=0.15;
    rdcutsvalmine[11][9]=0.34;
    rdcutsvalmine[12][9]=100.;
    rdcutsvalmine[13][9]=0.5;
    rdcutsvalmine[14][9]=0.992;
    rdcutsvalmine[15][9]=4.;

    //10-12
    rdcutsvalmine[0][10]=0.055;
    rdcutsvalmine[1][10]=0.03;
    rdcutsvalmine[2][10]=1.0;
    rdcutsvalmine[3][10]=1.0;
    rdcutsvalmine[4][10]=1.0;
    rdcutsvalmine[5][10]=0.15;
    rdcutsvalmine[6][10]=0.15;
    rdcutsvalmine[7][10]=-0.00004;//-0.00001;
    rdcutsvalmine[8][10]=0.95;//0.92;
    rdcutsvalmine[9][10]=0.3;
    rdcutsvalmine[10][10]=0.15;
    rdcutsvalmine[11][10]=0.35;
    rdcutsvalmine[12][10]=100.;
    rdcutsvalmine[13][10]=0.5;
    rdcutsvalmine[14][10]=0.992;
    rdcutsvalmine[15][10]=4.;

    //12-16 (ptbin11)
    rdcutsvalmine[0][11]=0.074;
    rdcutsvalmine[1][11]=0.04; ////////////////////////////0.04;//0.03750;
    rdcutsvalmine[2][11]=1.;//////////////////////////1.0;
    rdcutsvalmine[3][11]=1.0;//1.30;
    rdcutsvalmine[4][11]=1.0;//1.30;
    rdcutsvalmine[5][11]=0.2;
    rdcutsvalmine[6][11]=0.2;
    rdcutsvalmine[7][11]=-0.000023;
    rdcutsvalmine[8][11]=0.96;
    rdcutsvalmine[9][11]=0.3;
    rdcutsvalmine[10][11]=0.1;
    rdcutsvalmine[11][11]=0.45;
    rdcutsvalmine[12][11]=100.;
    rdcutsvalmine[13][11]=1;
    rdcutsvalmine[14][11]=0.996;////////////////////////0.996;
    rdcutsvalmine[15][11]=4.;//5.;


    //16-20 (ptbin12)
    //rdcutsvalmine[0][12]=0.074;
    //rdcutsvalmine[1][12]=0.8;
    //rdcutsvalmine[2][12]=1.0;
    //rdcutsvalmine[3][12]=0.30;
    //rdcutsvalmine[4][12]=0.30;
    //rdcutsvalmine[5][12]=0.60;
    //rdcutsvalmine[6][12]=0.60;
    //rdcutsvalmine[7][12]=0.0025;
    //rdcutsvalmine[8][12]=0.8;
    //rdcutsvalmine[9][12]=0.3;
    //rdcutsvalmine[10][12]=0.15;
    //rdcutsvalmine[11][12]=0.5;
    //rdcutsvalmine[12][12]=100.;
    //rdcutsvalmine[13][12]=1;
    //rdcutsvalmine[14][12]=0.995;//0.99;
    //rdcutsvalmine[15][12]=0.;

    //20-24
    //rdcutsvalmine[0][13]=0.074;
    //rdcutsvalmine[1][13]=0.8;
    //rdcutsvalmine[2][13]=1.0;
    //rdcutsvalmine[3][13]=0.30;
    //rdcutsvalmine[4][13]=0.30;
    //rdcutsvalmine[5][13]=0.60;
    //rdcutsvalmine[6][13]=0.60;
    //rdcutsvalmine[7][13]=0.0025;
    //rdcutsvalmine[8][13]=0.8;
    //rdcutsvalmine[9][13]=0.3;
    //rdcutsvalmine[10][13]=0.15;
    //rdcutsvalmine[11][13]=0.5;
    //rdcutsvalmine[12][13]=100.;
    //rdcutsvalmine[13][13]=1;
    //rdcutsvalmine[14][13]=0.99;
    //rdcutsvalmine[15][13]=0.;

    //16-20 (ptbin12)
    rdcutsvalmine[0][12]=0.074;
    rdcutsvalmine[1][12]=0.8;
    rdcutsvalmine[2][12]=1.0;
    rdcutsvalmine[3][12]=0.30;
    rdcutsvalmine[4][12]=0.30;
    rdcutsvalmine[5][12]=0.50;
    rdcutsvalmine[6][12]=0.50;
    rdcutsvalmine[7][12]=0.0033;
    rdcutsvalmine[8][12]=0.7;
    rdcutsvalmine[9][12]=0.1;
    rdcutsvalmine[10][12]=0.1;
    rdcutsvalmine[11][12]=0.05;
    rdcutsvalmine[12][12]=100.;
    rdcutsvalmine[13][12]=1;
    rdcutsvalmine[14][12]=0.9;
    rdcutsvalmine[15][12]=0.0;

    //20-24 (ptbin13)
    rdcutsvalmine[0][13]=0.074;
    rdcutsvalmine[1][13]=0.8;
    rdcutsvalmine[2][13]=1.0;
    rdcutsvalmine[3][13]=0.30;
    rdcutsvalmine[4][13]=0.30;
    rdcutsvalmine[5][13]=0.50;
    rdcutsvalmine[6][13]=0.50;
    rdcutsvalmine[7][13]=0.0033;
    rdcutsvalmine[8][13]=0.7;
    rdcutsvalmine[9][13]=0.1;
    rdcutsvalmine[10][13]=0.1;
    rdcutsvalmine[11][13]=0.05;
    rdcutsvalmine[12][13]=100.;
    rdcutsvalmine[13][13]=1;
    rdcutsvalmine[14][13]=0.9;
    rdcutsvalmine[15][13]=0.0;

    //>25 (ptbin14)
    rdcutsvalmine[0][14]=0.074;
    rdcutsvalmine[1][14]=0.8;
    rdcutsvalmine[2][14]=1.0;
    rdcutsvalmine[3][14]=0.30;
    rdcutsvalmine[4][14]=0.30;
    rdcutsvalmine[5][14]=0.50;
    rdcutsvalmine[6][14]=0.50;
    rdcutsvalmine[7][14]=0.0033;
    rdcutsvalmine[8][14]=0.7;
    rdcutsvalmine[9][14]=0.1;
    rdcutsvalmine[10][14]=0.1;
    rdcutsvalmine[11][14]=0.05;
    rdcutsvalmine[12][14]=100.;
    rdcutsvalmine[13][14]=1;
    rdcutsvalmine[14][14]=0.9;
    rdcutsvalmine[15][14]=0.0;
*/
/*
// Cuts 21_03_12
    //0-0.5
    rdcutsvalmine[0][0]=0.026;
    rdcutsvalmine[1][0]=0.022;
    rdcutsvalmine[2][0]=0.7;
    rdcutsvalmine[3][0]=0.21;
    rdcutsvalmine[4][0]=0.21;
    rdcutsvalmine[5][0]=0.05;
    rdcutsvalmine[6][0]=0.05;
    rdcutsvalmine[7][0]=0.0004;
    rdcutsvalmine[8][0]=0.6;
    rdcutsvalmine[9][0]=0.3;
    rdcutsvalmine[10][0]=0.1;
    rdcutsvalmine[11][0]=0.05;
    rdcutsvalmine[12][0]=100.;
    rdcutsvalmine[13][0]=0.5;
    rdcutsvalmine[14][0]=-1.;
    rdcutsvalmine[15][0]=0.;

    //0.5-1
    rdcutsvalmine[0][1]=0.026;
    rdcutsvalmine[1][1]=0.04;
    rdcutsvalmine[2][1]=0.8;
    rdcutsvalmine[3][1]=0.5;
    rdcutsvalmine[4][1]=0.5;
    rdcutsvalmine[5][1]=0.08;
    rdcutsvalmine[6][1]=0.08;
    rdcutsvalmine[7][1]=0.0003;
    rdcutsvalmine[8][1]=0.65;
    rdcutsvalmine[9][1]=0.3;
    rdcutsvalmine[10][1]=0.1;
    rdcutsvalmine[11][1]=0.05;
    rdcutsvalmine[12][1]=100.;
    rdcutsvalmine[13][1]=0.5;
    rdcutsvalmine[14][1]=-1.;
    rdcutsvalmine[15][1]=0.;

    //1-2
    rdcutsvalmine[0][2]=0.022;
    rdcutsvalmine[1][2]=0.04;
    rdcutsvalmine[2][2]=0.8;
    rdcutsvalmine[3][2]=0.7;
    rdcutsvalmine[4][2]=0.7;
    rdcutsvalmine[5][2]=0.08;
    rdcutsvalmine[6][2]=0.08;
    rdcutsvalmine[7][2]=0.00014;
    rdcutsvalmine[8][2]=0.75;
    rdcutsvalmine[9][2]=0.3;
    rdcutsvalmine[10][2]=0.1;
    rdcutsvalmine[11][2]=0.05;
    rdcutsvalmine[12][2]=100.;
    rdcutsvalmine[13][2]=0.5;
    rdcutsvalmine[14][2]=-1.;
    rdcutsvalmine[15][2]=0.;
    //2-3
    rdcutsvalmine[0][3]=0.024;
    rdcutsvalmine[1][3]=0.025;
    rdcutsvalmine[2][3]=0.8;
    rdcutsvalmine[3][3]=0.7;
    rdcutsvalmine[4][3]=0.7;
    rdcutsvalmine[5][3]=0.08;
    rdcutsvalmine[6][3]=0.08;
    rdcutsvalmine[7][3]=0.00014;
    rdcutsvalmine[8][3]=0.85;
    rdcutsvalmine[9][3]=0.3;
    rdcutsvalmine[10][3]=0.1;
    rdcutsvalmine[11][3]=0.05;
    rdcutsvalmine[12][3]=100.;
    rdcutsvalmine[13][3]=0.5;
    rdcutsvalmine[14][3]=-1.;
    rdcutsvalmine[15][3]=0.;

    //3-4
    rdcutsvalmine[0][4]=0.032;
    rdcutsvalmine[1][4]=0.01250;	//0.1;//0.025;
    rdcutsvalmine[2][4]=1.0;			//1.0;//0.7;
    rdcutsvalmine[3][4]=1.333;		//0.5;//1.1;
    rdcutsvalmine[4][4]=1.0;			//0.5;//1.1;
    rdcutsvalmine[5][4]=0.1;			//0.1;//0.08;
    rdcutsvalmine[6][4]=0.05833;		//0.1;//0.08;
    rdcutsvalmine[7][4]=-0.000285;		//-0.00017;
    rdcutsvalmine[8][4]=0.99;			//0.98;
    rdcutsvalmine[9][4]=0.3;
    rdcutsvalmine[10][4]=0.15;
    rdcutsvalmine[11][4]=0.05;
    rdcutsvalmine[12][4]=0.5;
    rdcutsvalmine[13][4]=0.5;
    rdcutsvalmine[14][4]=0.998;			//0.99;
    rdcutsvalmine[15][4]=7.0;			//5.;

    //4-5
    rdcutsvalmine[0][5]=0.032;
    rdcutsvalmine[1][5]=0.0375; //0.1
    rdcutsvalmine[2][5]=0.5625;
    rdcutsvalmine[3][5]=1.;
    rdcutsvalmine[4][5]=1.;
    rdcutsvalmine[5][5]=0.075;
    rdcutsvalmine[6][5]=0.0583;
    rdcutsvalmine[7][5]=-0.000362;
    rdcutsvalmine[8][5]=0.99;
    rdcutsvalmine[9][5]=0.15;
    rdcutsvalmine[10][5]=0.15;
    rdcutsvalmine[11][5]=0.05;
    rdcutsvalmine[12][5]=10;
    rdcutsvalmine[13][5]=0.5;
    rdcutsvalmine[14][5]=0.99;
    rdcutsvalmine[15][5]=8.;

    //5-6
    rdcutsvalmine[0][6]=0.04;
    rdcutsvalmine[1][6]=0.025;
    rdcutsvalmine[2][6]=1.;
    rdcutsvalmine[3][6]=1.33;
    rdcutsvalmine[4][6]=1.167;
    rdcutsvalmine[5][6]=0.075;
    rdcutsvalmine[6][6]=0.075;
    rdcutsvalmine[7][6]=-0.000287;
    rdcutsvalmine[8][6]=0.99;
    rdcutsvalmine[9][6]=0.15;
    rdcutsvalmine[10][6]=0.15;
    rdcutsvalmine[11][6]=0.2;
    rdcutsvalmine[12][6]=10.;
    rdcutsvalmine[13][6]=0.5;
    rdcutsvalmine[14][6]=0.993;
    rdcutsvalmine[15][6]=6.5;

    //6-7
    rdcutsvalmine[0][7]=0.036;
    rdcutsvalmine[1][7]=0.03750;
    rdcutsvalmine[2][7]=1.0;
    rdcutsvalmine[3][7]=1.333;
    rdcutsvalmine[4][7]=1.333;
    rdcutsvalmine[5][7]=0.075;
    rdcutsvalmine[6][7]=0.075;
    rdcutsvalmine[7][7]=-0.000243;
    rdcutsvalmine[8][7]=0.986;
    rdcutsvalmine[9][7]=0.3;
    rdcutsvalmine[10][7]=0.1;
    rdcutsvalmine[11][7]=0.2;
    rdcutsvalmine[12][7]=100.;
    rdcutsvalmine[13][7]=0.5;
    rdcutsvalmine[14][7]=0.992;
    rdcutsvalmine[15][7]=4.0;

    //7-8
    rdcutsvalmine[0][8]=0.036;
    rdcutsvalmine[1][8]=0.0125;
    rdcutsvalmine[2][8]=0.5625;
    rdcutsvalmine[3][8]=0.833;
    rdcutsvalmine[4][8]=1.333;
    rdcutsvalmine[5][8]=0.075;
    rdcutsvalmine[6][8]=0.2;
    rdcutsvalmine[7][8]=-0.000163;
    rdcutsvalmine[8][8]=0.967;
    rdcutsvalmine[9][8]=0.3;
    rdcutsvalmine[10][8]=0.15;
    rdcutsvalmine[11][8]=0.2;
    rdcutsvalmine[12][8]=100.;
    rdcutsvalmine[13][8]=0.5;
    rdcutsvalmine[14][8]=0.997;
    rdcutsvalmine[15][8]=9.;

    //8-10		230312: changed to cut of bin10 to increase signal
    rdcutsvalmine[0][9]=0.055;
    rdcutsvalmine[1][9]=0.03;
    rdcutsvalmine[2][9]=1.0;
    rdcutsvalmine[3][9]=1.0;			//0.9;
    rdcutsvalmine[4][9]=1.0;			//0.9;
    rdcutsvalmine[5][9]=0.15;
    rdcutsvalmine[6][9]=0.15;
    rdcutsvalmine[7][9]=-0.00004;		//-0.00008;
    rdcutsvalmine[8][9]=0.99;
    rdcutsvalmine[9][9]=0.3;
    rdcutsvalmine[10][9]=0.15;
    rdcutsvalmine[11][9]=0.34;
    rdcutsvalmine[12][9]=100.;
    rdcutsvalmine[13][9]=0.5;
    rdcutsvalmine[14][9]=0.998;
    rdcutsvalmine[15][9]=4.;			//4.8;

    //10-12     230312:
    rdcutsvalmine[0][10]=0.055;
    rdcutsvalmine[1][10]=0.03;
    rdcutsvalmine[2][10]=1.0;
    rdcutsvalmine[3][10]=1.0;
    rdcutsvalmine[4][10]=1.0;
    rdcutsvalmine[5][10]=0.15;
    rdcutsvalmine[6][10]=0.15;
    rdcutsvalmine[7][10]=-0.00004;
    rdcutsvalmine[8][10]=0.99;
    rdcutsvalmine[9][10]=0.3;
    rdcutsvalmine[10][10]=0.15;
    rdcutsvalmine[11][10]=0.35;
    rdcutsvalmine[12][10]=100.;
    rdcutsvalmine[13][10]=0.5;
    rdcutsvalmine[14][10]=0.998;
    rdcutsvalmine[15][10]=4.;

    //12-16 (ptbin11)
    //rdcutsvalmine[0][11]=0.074;
    //rdcutsvalmine[1][11]=0.0125;
    //rdcutsvalmine[2][11]=0.5625;
    //rdcutsvalmine[3][11]=1.30;
    //rdcutsvalmine[4][11]=1.30;
    //rdcutsvalmine[5][11]=0.2;
    //rdcutsvalmine[6][11]=0.2;
    //rdcutsvalmine[7][11]=-0.000023;
    //rdcutsvalmine[8][11]=0.97;
    //rdcutsvalmine[9][11]=0.3;
    //rdcutsvalmine[10][11]=0.1;
    //rdcutsvalmine[11][11]=0.45;
    //rdcutsvalmine[12][11]=100.;
    //rdcutsvalmine[13][11]=1;
    //rdcutsvalmine[14][11]=0.999;
    //rdcutsvalmine[15][11]=5.;

    //12-16 (ptbin11)
    rdcutsvalmine[0][11]=0.074;
    rdcutsvalmine[1][11]=0.025; ////////////////////////////0.04;//0.03750;
    rdcutsvalmine[2][11]=0.8;//////////////////////////1.0;
    rdcutsvalmine[3][11]=1.0;//1.30;
    rdcutsvalmine[4][11]=1.0;//1.30;
    rdcutsvalmine[5][11]=0.2;
    rdcutsvalmine[6][11]=0.2;
    rdcutsvalmine[7][11]=-0.000023;
    rdcutsvalmine[8][11]=0.97;
    rdcutsvalmine[9][11]=0.3;
    rdcutsvalmine[10][11]=0.1;
    rdcutsvalmine[11][11]=0.45;
    rdcutsvalmine[12][11]=100.;
    rdcutsvalmine[13][11]=1;
    rdcutsvalmine[14][11]=0.998;////////////////////////0.996;
    rdcutsvalmine[15][11]=4.;//5.;


    //16-20 (ptbin12)
    //rdcutsvalmine[0][12]=0.074;
    //rdcutsvalmine[1][12]=0.8;
    //rdcutsvalmine[2][12]=0.5625;
    //rdcutsvalmine[3][12]=0.30;
    //rdcutsvalmine[4][12]=0.30;
    //rdcutsvalmine[5][12]=0.60;
    //rdcutsvalmine[6][12]=0.60;
    //rdcutsvalmine[7][12]=0.0002;
    //rdcutsvalmine[8][12]=0.97;
    //rdcutsvalmine[9][12]=0.3;
    //rdcutsvalmine[10][12]=0.15;
    //rdcutsvalmine[11][12]=0.5;
    //rdcutsvalmine[12][12]=100.;
    //rdcutsvalmine[13][12]=1;
    //rdcutsvalmine[14][12]=0.998;
    //rdcutsvalmine[15][12]=5.;


    //20-24
    //rdcutsvalmine[0][13]=0.074;
    //rdcutsvalmine[1][13]=0.8;
    //rdcutsvalmine[2][13]=1.0;
    //rdcutsvalmine[3][13]=0.30;
    //rdcutsvalmine[4][13]=0.30;
    //rdcutsvalmine[5][13]=0.60;
    //rdcutsvalmine[6][13]=0.60;
    //rdcutsvalmine[7][13]=0.0025;
    //rdcutsvalmine[8][13]=0.8;
    //rdcutsvalmine[9][13]=0.3;
    //rdcutsvalmine[10][13]=0.15;
    //rdcutsvalmine[11][13]=0.5;
    //rdcutsvalmine[12][13]=100.;
    //rdcutsvalmine[13][13]=1;
    //rdcutsvalmine[14][13]=0.99;
    //rdcutsvalmine[15][13]=0.;

    //16-20 (ptbin12)
    rdcutsvalmine[0][12]=0.074;
    rdcutsvalmine[1][12]=0.8;
    rdcutsvalmine[2][12]=1.0;
    rdcutsvalmine[3][12]=0.30;
    rdcutsvalmine[4][12]=0.30;
    rdcutsvalmine[5][12]=0.50;
    rdcutsvalmine[6][12]=0.50;
    rdcutsvalmine[7][12]=0.0033;
    rdcutsvalmine[8][12]=0.7;
    rdcutsvalmine[9][12]=0.1;
    rdcutsvalmine[10][12]=0.1;
    rdcutsvalmine[11][12]=0.05;
    rdcutsvalmine[12][12]=100.;
    rdcutsvalmine[13][12]=1;
    rdcutsvalmine[14][12]=0.9;
    rdcutsvalmine[15][12]=0.0;

    //20-24 (ptbin13)
    rdcutsvalmine[0][13]=0.074;
    rdcutsvalmine[1][13]=0.8;
    rdcutsvalmine[2][13]=1.0;
    rdcutsvalmine[3][13]=0.30;
    rdcutsvalmine[4][13]=0.30;
    rdcutsvalmine[5][13]=0.50;
    rdcutsvalmine[6][13]=0.50;
    rdcutsvalmine[7][13]=0.0033;
    rdcutsvalmine[8][13]=0.7;
    rdcutsvalmine[9][13]=0.1;
    rdcutsvalmine[10][13]=0.1;
    rdcutsvalmine[11][13]=0.05;
    rdcutsvalmine[12][13]=100.;
    rdcutsvalmine[13][13]=1;
    rdcutsvalmine[14][13]=0.9;
    rdcutsvalmine[15][13]=0.0;

    //>25 (ptbin14)
    rdcutsvalmine[0][14]=0.074;
    rdcutsvalmine[1][14]=0.8;
    rdcutsvalmine[2][14]=1.0;
    rdcutsvalmine[3][14]=0.30;
    rdcutsvalmine[4][14]=0.30;
    rdcutsvalmine[5][14]=0.50;
    rdcutsvalmine[6][14]=0.50;
    rdcutsvalmine[7][14]=0.0033;
    rdcutsvalmine[8][14]=0.7;
    rdcutsvalmine[9][14]=0.1;
    rdcutsvalmine[10][14]=0.1;
    rdcutsvalmine[11][14]=0.05;
    rdcutsvalmine[12][14]=100.;
    rdcutsvalmine[13][14]=1;
    rdcutsvalmine[14][14]=0.9;
    rdcutsvalmine[15][14]=0.0;
*/
/*
// Cuts 08_03_12 Collection Raoul with Alessandro Modifications
    //0-0.5
    rdcutsvalmine[0][0]=0.026;
    rdcutsvalmine[1][0]=0.022;
    rdcutsvalmine[2][0]=0.7;
    rdcutsvalmine[3][0]=0.21;
    rdcutsvalmine[4][0]=0.21;
    rdcutsvalmine[5][0]=0.05;
    rdcutsvalmine[6][0]=0.05;
    rdcutsvalmine[7][0]=0.0004;
    rdcutsvalmine[8][0]=0.6;
    rdcutsvalmine[9][0]=0.3;
    rdcutsvalmine[10][0]=0.1;
    rdcutsvalmine[11][0]=0.05;
    rdcutsvalmine[12][0]=100.;
    rdcutsvalmine[13][0]=0.5;
    rdcutsvalmine[14][0]=-1.;
    rdcutsvalmine[15][0]=0.;

    //0.5-1
    rdcutsvalmine[0][1]=0.026;
    rdcutsvalmine[1][1]=0.04;
    rdcutsvalmine[2][1]=0.8;
    rdcutsvalmine[3][1]=0.5;
    rdcutsvalmine[4][1]=0.5;
    rdcutsvalmine[5][1]=0.08;
    rdcutsvalmine[6][1]=0.08;
    rdcutsvalmine[7][1]=0.0003;
    rdcutsvalmine[8][1]=0.65;
    rdcutsvalmine[9][1]=0.3;
    rdcutsvalmine[10][1]=0.1;
    rdcutsvalmine[11][1]=0.05;
    rdcutsvalmine[12][1]=100.;
    rdcutsvalmine[13][1]=0.5;
    rdcutsvalmine[14][1]=-1.;
    rdcutsvalmine[15][1]=0.;

    //1-2
    rdcutsvalmine[0][2]=0.022;
    rdcutsvalmine[1][2]=0.04;
    rdcutsvalmine[2][2]=0.8;
    rdcutsvalmine[3][2]=0.7;
    rdcutsvalmine[4][2]=0.7;
    rdcutsvalmine[5][2]=0.08;
    rdcutsvalmine[6][2]=0.08;
    rdcutsvalmine[7][2]=0.00014;
    rdcutsvalmine[8][2]=0.75;
    rdcutsvalmine[9][2]=0.3;
    rdcutsvalmine[10][2]=0.1;
    rdcutsvalmine[11][2]=0.05;
    rdcutsvalmine[12][2]=100.;
    rdcutsvalmine[13][2]=0.5;
    rdcutsvalmine[14][2]=-1.;
    rdcutsvalmine[15][2]=0.;
    //2-3
    rdcutsvalmine[0][3]=0.024;
    rdcutsvalmine[1][3]=0.025;
    rdcutsvalmine[2][3]=0.8;
    rdcutsvalmine[3][3]=0.7;
    rdcutsvalmine[4][3]=0.7;
    rdcutsvalmine[5][3]=0.08;
    rdcutsvalmine[6][3]=0.08;
    rdcutsvalmine[7][3]=0.00014;
    rdcutsvalmine[8][3]=0.85;
    rdcutsvalmine[9][3]=0.3;
    rdcutsvalmine[10][3]=0.1;
    rdcutsvalmine[11][3]=0.05;
    rdcutsvalmine[12][3]=100.;
    rdcutsvalmine[13][3]=0.5;
    rdcutsvalmine[14][3]=-1.;
    rdcutsvalmine[15][3]=0.;

    //3-4
    rdcutsvalmine[0][4]=0.032;
    rdcutsvalmine[1][4]=0.02;//0.01250;		//0.1;//0.025;
    rdcutsvalmine[2][4]=1.0;			//1.0;//0.7;
    rdcutsvalmine[3][4]=1.2;//1.333;		//0.5;//1.1;
    rdcutsvalmine[4][4]=1.0;			//0.5;//1.1;
    rdcutsvalmine[5][4]=0.1;			//0.1;//0.08;
    rdcutsvalmine[6][4]=0.1;//0.05833;		//0.1;//0.08;
    rdcutsvalmine[7][4]=-0.000285;		//-0.00017;
    rdcutsvalmine[8][4]=0.99;			//0.98;
    rdcutsvalmine[9][4]=0.3;
    rdcutsvalmine[10][4]=0.15;
    rdcutsvalmine[11][4]=0.05;
    rdcutsvalmine[12][4]=0.5;
    rdcutsvalmine[13][4]=0.5;
    rdcutsvalmine[14][4]=0.998;			//0.99;
    rdcutsvalmine[15][4]=7.0;			//5.;

    //4-5
    rdcutsvalmine[0][5]=0.032;
    rdcutsvalmine[1][5]=0.025;
    rdcutsvalmine[2][5]=0.7;
    rdcutsvalmine[3][5]=1.1;//1.2;
    rdcutsvalmine[4][5]=1.1;//1.2;
    rdcutsvalmine[5][5]=0.1;//0.09;
    rdcutsvalmine[6][5]=0.1;//0.09;
    rdcutsvalmine[7][5]=-0.000337;
    rdcutsvalmine[8][5]=0.992;
    rdcutsvalmine[9][5]=0.15;
    rdcutsvalmine[10][5]=0.15;
    rdcutsvalmine[11][5]=0.05;
    rdcutsvalmine[12][5]=10;
    rdcutsvalmine[13][5]=0.5;
    rdcutsvalmine[14][5]=0.998;//0.996;
    rdcutsvalmine[15][5]=8.0;//8.5;

    //5-6
    rdcutsvalmine[0][6]=0.04;
    rdcutsvalmine[1][6]=0.027;
    rdcutsvalmine[2][6]=0.8;
    rdcutsvalmine[3][6]=1.1;
    rdcutsvalmine[4][6]=1.1;
    rdcutsvalmine[5][6]=0.1;
    rdcutsvalmine[6][6]=0.1;
    rdcutsvalmine[7][6]=-0.0003;
    rdcutsvalmine[8][6]=0.99;
    rdcutsvalmine[9][6]=0.15;
    rdcutsvalmine[10][6]=0.15;
    rdcutsvalmine[11][6]=0.2;
    rdcutsvalmine[12][6]=10.;
    rdcutsvalmine[13][6]=0.5;
    rdcutsvalmine[14][6]=0.998;//0.995;
    rdcutsvalmine[15][6]=5.0;//6.0;

    //6-7
    rdcutsvalmine[0][7]=0.036;
    rdcutsvalmine[1][7]=0.027;//0.03750;
    rdcutsvalmine[2][7]=1.0;
    rdcutsvalmine[3][7]=1.2;//1.333;
    rdcutsvalmine[4][7]=1.2;//1.333;
    rdcutsvalmine[5][7]=0.075;
    rdcutsvalmine[6][7]=0.075;
    rdcutsvalmine[7][7]=-0.000243;
    rdcutsvalmine[8][7]=0.98;//0.986;
    rdcutsvalmine[9][7]=0.3;
    rdcutsvalmine[10][7]=0.1;
    rdcutsvalmine[11][7]=0.2;
    rdcutsvalmine[12][7]=100.;
    rdcutsvalmine[13][7]=0.5;
    rdcutsvalmine[14][7]=0.998;//0.992;
    rdcutsvalmine[15][7]=4.0;

    //7-8
    rdcutsvalmine[0][8]=0.036;
    rdcutsvalmine[1][8]=0.027;//0.03750;
    rdcutsvalmine[2][8]=1.0;
    rdcutsvalmine[3][8]=0.833;
    rdcutsvalmine[4][8]=0.833;
    rdcutsvalmine[5][8]=0.2;
    rdcutsvalmine[6][8]=0.2;
    rdcutsvalmine[7][8]=-0.000163;
    rdcutsvalmine[8][8]=0.967;
    rdcutsvalmine[9][8]=0.3;
    rdcutsvalmine[10][8]=0.15;
    rdcutsvalmine[11][8]=0.2;
    rdcutsvalmine[12][8]=100.;
    rdcutsvalmine[13][8]=0.5;
    rdcutsvalmine[14][8]=0.995;
    rdcutsvalmine[15][8]=4.0;//9.;

    //8-10
    rdcutsvalmine[0][9]=0.055;
    rdcutsvalmine[1][9]=0.03;
    rdcutsvalmine[2][9]=1.0;
    rdcutsvalmine[3][9]=0.9;
    rdcutsvalmine[4][9]=0.9;
    rdcutsvalmine[5][9]=0.15;
    rdcutsvalmine[6][9]=0.15;
    rdcutsvalmine[7][9]=-0.00008;	//-0.00001;
    rdcutsvalmine[8][9]=0.94;//0.99;		//0.92;
    rdcutsvalmine[9][9]=0.3;
    rdcutsvalmine[10][9]=0.15;
    rdcutsvalmine[11][9]=0.34;
    rdcutsvalmine[12][9]=100.;
    rdcutsvalmine[13][9]=0.5;
     rdcutsvalmine[14][9]=0.998;
    rdcutsvalmine[15][9]=4.0;//4.8;

    //10-12
    rdcutsvalmine[0][10]=0.055;
    rdcutsvalmine[1][10]=0.03;
    rdcutsvalmine[2][10]=1.0;
    rdcutsvalmine[3][10]=1.0;
    rdcutsvalmine[4][10]=1.0;
    rdcutsvalmine[5][10]=0.15;
    rdcutsvalmine[6][10]=0.15;
    rdcutsvalmine[7][10]=-0.00004;//-0.00001;
    rdcutsvalmine[8][10]=0.99;//0.92;
    rdcutsvalmine[9][10]=0.3;
    rdcutsvalmine[10][10]=0.15;
    rdcutsvalmine[11][10]=0.35;
    rdcutsvalmine[12][10]=100.;
    rdcutsvalmine[13][10]=0.5;
     rdcutsvalmine[14][10]=0.998;
    rdcutsvalmine[15][10]=4.;

    //12-16 (ptbin11)
    rdcutsvalmine[0][11]=0.074;
    rdcutsvalmine[1][11]=0.025; ////////////////////////////0.04;//0.03750;
    rdcutsvalmine[2][11]=0.8;//////////////////////////1.0;
    rdcutsvalmine[3][11]=1.0;//1.30;
    rdcutsvalmine[4][11]=1.0;//1.30;
    rdcutsvalmine[5][11]=0.2;
    rdcutsvalmine[6][11]=0.2;
    rdcutsvalmine[7][11]=-0.000023;
    rdcutsvalmine[8][11]=0.97;
    rdcutsvalmine[9][11]=0.3;
    rdcutsvalmine[10][11]=0.1;
    rdcutsvalmine[11][11]=0.45;
    rdcutsvalmine[12][11]=100.;
    rdcutsvalmine[13][11]=1;
    rdcutsvalmine[14][11]=0.998;////////////////////////0.996;
    rdcutsvalmine[15][11]=4.;//5.;

    //16-20 (ptbin12)
    rdcutsvalmine[0][12]=0.074;
    rdcutsvalmine[1][12]=0.8;
    rdcutsvalmine[2][12]=1.0;
    rdcutsvalmine[3][12]=0.30;
    rdcutsvalmine[4][12]=0.30;
    rdcutsvalmine[5][12]=0.60;
    rdcutsvalmine[6][12]=0.60;
    rdcutsvalmine[7][12]=0.0025;
    rdcutsvalmine[8][12]=0.8;
    rdcutsvalmine[9][12]=0.3;
    rdcutsvalmine[10][12]=0.15;
    rdcutsvalmine[11][12]=0.5;
    rdcutsvalmine[12][12]=100.;
    rdcutsvalmine[13][12]=1;
    rdcutsvalmine[14][12]=0.995;//0.99;
    rdcutsvalmine[15][12]=0.;


    //20-24
    rdcutsvalmine[0][13]=0.074;
    rdcutsvalmine[1][13]=0.8;
    rdcutsvalmine[2][13]=1.0;
    rdcutsvalmine[3][13]=0.30;
    rdcutsvalmine[4][13]=0.30;
    rdcutsvalmine[5][13]=0.60;
    rdcutsvalmine[6][13]=0.60;
    rdcutsvalmine[7][13]=0.0025;
    rdcutsvalmine[8][13]=0.8;
    rdcutsvalmine[9][13]=0.3;
    rdcutsvalmine[10][13]=0.15;
    rdcutsvalmine[11][13]=0.5;
    rdcutsvalmine[12][13]=100.;
    rdcutsvalmine[13][13]=1;
    rdcutsvalmine[14][13]=0.99;
    rdcutsvalmine[15][13]=0.;

    //>25 (ptbin14)
    rdcutsvalmine[0][14]=0.074;
    rdcutsvalmine[1][14]=0.8;
    rdcutsvalmine[2][14]=1.0;
    rdcutsvalmine[3][14]=0.30;
    rdcutsvalmine[4][14]=0.30;
    rdcutsvalmine[5][14]=0.50;
    rdcutsvalmine[6][14]=0.50;
    rdcutsvalmine[7][14]=0.0033;
    rdcutsvalmine[8][14]=0.7;
    rdcutsvalmine[9][14]=0.1;
    rdcutsvalmine[10][14]=0.1;
    rdcutsvalmine[11][14]=0.05;
    rdcutsvalmine[12][14]=100.;
    rdcutsvalmine[13][14]=1;
    rdcutsvalmine[14][14]=0.9;
    rdcutsvalmine[15][14]=0.0;
*/
/*
// Cuts 07_03_12 Collection Raoul
    //0-0.5
    rdcutsvalmine[0][0]=0.026;
    rdcutsvalmine[1][0]=0.022;
    rdcutsvalmine[2][0]=0.7;
    rdcutsvalmine[3][0]=0.21;
    rdcutsvalmine[4][0]=0.21;
    rdcutsvalmine[5][0]=0.05;
    rdcutsvalmine[6][0]=0.05;
    rdcutsvalmine[7][0]=0.0004;
    rdcutsvalmine[8][0]=0.6;
    rdcutsvalmine[9][0]=0.3;
    rdcutsvalmine[10][0]=0.1;
    rdcutsvalmine[11][0]=0.05;
    rdcutsvalmine[12][0]=100.;
    rdcutsvalmine[13][0]=0.5;
    rdcutsvalmine[14][0]=-1.;
    rdcutsvalmine[15][0]=0.;

    //0.5-1
    rdcutsvalmine[0][1]=0.026;
    rdcutsvalmine[1][1]=0.04;
    rdcutsvalmine[2][1]=0.8;
    rdcutsvalmine[3][1]=0.5;
    rdcutsvalmine[4][1]=0.5;
    rdcutsvalmine[5][1]=0.08;
    rdcutsvalmine[6][1]=0.08;
    rdcutsvalmine[7][1]=0.0003;
    rdcutsvalmine[8][1]=0.65;
    rdcutsvalmine[9][1]=0.3;
    rdcutsvalmine[10][1]=0.1;
    rdcutsvalmine[11][1]=0.05;
    rdcutsvalmine[12][1]=100.;
    rdcutsvalmine[13][1]=0.5;
    rdcutsvalmine[14][1]=-1.;
    rdcutsvalmine[15][1]=0.;

    //1-2
    rdcutsvalmine[0][2]=0.022;
    rdcutsvalmine[1][2]=0.04;
    rdcutsvalmine[2][2]=0.8;
    rdcutsvalmine[3][2]=0.7;
    rdcutsvalmine[4][2]=0.7;
    rdcutsvalmine[5][2]=0.08;
    rdcutsvalmine[6][2]=0.08;
    rdcutsvalmine[7][2]=0.00014;
    rdcutsvalmine[8][2]=0.75;
    rdcutsvalmine[9][2]=0.3;
    rdcutsvalmine[10][2]=0.1;
    rdcutsvalmine[11][2]=0.05;
    rdcutsvalmine[12][2]=100.;
    rdcutsvalmine[13][2]=0.5;
    rdcutsvalmine[14][2]=-1.;
    rdcutsvalmine[15][2]=0.;
    //2-3
    rdcutsvalmine[0][3]=0.024;
    rdcutsvalmine[1][3]=0.025;
    rdcutsvalmine[2][3]=0.8;
    rdcutsvalmine[3][3]=0.7;
    rdcutsvalmine[4][3]=0.7;
    rdcutsvalmine[5][3]=0.08;
    rdcutsvalmine[6][3]=0.08;
    rdcutsvalmine[7][3]=0.00014;
    rdcutsvalmine[8][3]=0.85;
    rdcutsvalmine[9][3]=0.3;
    rdcutsvalmine[10][3]=0.1;
    rdcutsvalmine[11][3]=0.05;
    rdcutsvalmine[12][3]=100.;
    rdcutsvalmine[13][3]=0.5;
    rdcutsvalmine[14][3]=-1.;
    rdcutsvalmine[15][3]=0.;

    //3-4
    rdcutsvalmine[0][4]=0.032;
    rdcutsvalmine[1][4]=0.01250;	//0.1;//0.025;
    rdcutsvalmine[2][4]=1.0;			//1.0;//0.7;
    rdcutsvalmine[3][4]=1.333;		//0.5;//1.1;
    rdcutsvalmine[4][4]=1.0;			//0.5;//1.1;
    rdcutsvalmine[5][4]=0.1;			//0.1;//0.08;
    rdcutsvalmine[6][4]=0.05833;		//0.1;//0.08;
    rdcutsvalmine[7][4]=-0.000285;		//-0.00017;
    rdcutsvalmine[8][4]=0.99;			//0.98;
    rdcutsvalmine[9][4]=0.3;
    rdcutsvalmine[10][4]=0.15;
    rdcutsvalmine[11][4]=0.05;
    rdcutsvalmine[12][4]=0.5;
    rdcutsvalmine[13][4]=0.5;
    rdcutsvalmine[14][4]=0.998;			//0.99;
    rdcutsvalmine[15][4]=7.0;			//5.;

    //4-5
    rdcutsvalmine[0][5]=0.032;
    rdcutsvalmine[1][5]=0.0375; //0.1
    rdcutsvalmine[2][5]=0.5625;
    rdcutsvalmine[3][5]=1.;
    rdcutsvalmine[4][5]=1.;
    rdcutsvalmine[5][5]=0.075;
    rdcutsvalmine[6][5]=0.0583;
    rdcutsvalmine[7][5]=-0.000362;
    rdcutsvalmine[8][5]=0.99;
    rdcutsvalmine[9][5]=0.15;
    rdcutsvalmine[10][5]=0.15;
    rdcutsvalmine[11][5]=0.05;
    rdcutsvalmine[12][5]=10;
    rdcutsvalmine[13][5]=0.5;
    rdcutsvalmine[14][5]=0.99;
    rdcutsvalmine[15][5]=8.;

    //5-6
    rdcutsvalmine[0][6]=0.04;
    rdcutsvalmine[1][6]=0.025;
    rdcutsvalmine[2][6]=1.;
    rdcutsvalmine[3][6]=1.33;
    rdcutsvalmine[4][6]=1.167;
    rdcutsvalmine[5][6]=0.075;
    rdcutsvalmine[6][6]=0.075;
    rdcutsvalmine[7][6]=-0.000287;
    rdcutsvalmine[8][6]=0.99;
    rdcutsvalmine[9][6]=0.15;
    rdcutsvalmine[10][6]=0.15;
    rdcutsvalmine[11][6]=0.2;
    rdcutsvalmine[12][6]=10.;
    rdcutsvalmine[13][6]=0.5;
    rdcutsvalmine[14][6]=0.993;
    rdcutsvalmine[15][6]=6.5;

    //6-7
    rdcutsvalmine[0][7]=0.036;
    rdcutsvalmine[1][7]=0.03750;
    rdcutsvalmine[2][7]=1.0;
    rdcutsvalmine[3][7]=1.333;
    rdcutsvalmine[4][7]=1.333;
    rdcutsvalmine[5][7]=0.075;
    rdcutsvalmine[6][7]=0.075;
    rdcutsvalmine[7][7]=-0.000243;
    rdcutsvalmine[8][7]=0.986;
    rdcutsvalmine[9][7]=0.3;
    rdcutsvalmine[10][7]=0.1;
    rdcutsvalmine[11][7]=0.2;
    rdcutsvalmine[12][7]=100.;
    rdcutsvalmine[13][7]=0.5;
    rdcutsvalmine[14][7]=0.992;
    rdcutsvalmine[15][7]=4.0;

    //7-8
    rdcutsvalmine[0][8]=0.036;
    rdcutsvalmine[1][8]=0.0125;
    rdcutsvalmine[2][8]=0.5625;
    rdcutsvalmine[3][8]=0.833;
    rdcutsvalmine[4][8]=1.333;
    rdcutsvalmine[5][8]=0.075;
    rdcutsvalmine[6][8]=0.2;
    rdcutsvalmine[7][8]=-0.000163;
    rdcutsvalmine[8][8]=0.967;
    rdcutsvalmine[9][8]=0.3;
    rdcutsvalmine[10][8]=0.15;
    rdcutsvalmine[11][8]=0.2;
    rdcutsvalmine[12][8]=100.;
    rdcutsvalmine[13][8]=0.5;
    rdcutsvalmine[14][8]=0.997;
    rdcutsvalmine[15][8]=9.;

    //8-10
    rdcutsvalmine[0][9]=0.055;
    rdcutsvalmine[1][9]=0.03;
    rdcutsvalmine[2][9]=1.0;
    rdcutsvalmine[3][9]=0.9;
    rdcutsvalmine[4][9]=0.9;
    rdcutsvalmine[5][9]=0.15;
    rdcutsvalmine[6][9]=0.15;
    rdcutsvalmine[7][9]=-0.00008;//-0.00001;
    rdcutsvalmine[8][9]=0.99;//0.92;
    rdcutsvalmine[9][9]=0.3;
    rdcutsvalmine[10][9]=0.15;
    rdcutsvalmine[11][9]=0.34;
    rdcutsvalmine[12][9]=100.;
    rdcutsvalmine[13][9]=0.5;
    rdcutsvalmine[14][9]=0.998;
    rdcutsvalmine[15][9]=4.8;

    //10-12
    rdcutsvalmine[0][10]=0.055;
    rdcutsvalmine[1][10]=0.03;
    rdcutsvalmine[2][10]=1.0;
    rdcutsvalmine[3][10]=1.0;
    rdcutsvalmine[4][10]=1.0;
    rdcutsvalmine[5][10]=0.15;
    rdcutsvalmine[6][10]=0.15;
    rdcutsvalmine[7][10]=-0.00004;//-0.00001;
    rdcutsvalmine[8][10]=0.99;//0.92;
    rdcutsvalmine[9][10]=0.3;
    rdcutsvalmine[10][10]=0.15;
    rdcutsvalmine[11][10]=0.35;
    rdcutsvalmine[12][10]=100.;
    rdcutsvalmine[13][10]=0.5;
    rdcutsvalmine[14][10]=0.998;
    rdcutsvalmine[15][10]=4.;

    //12-16 (ptbin11)
    rdcutsvalmine[0][11]=0.074;
    rdcutsvalmine[1][11]=0.0125;
    rdcutsvalmine[2][11]=0.5625;
    rdcutsvalmine[3][11]=1.30;
    rdcutsvalmine[4][11]=1.30;
    rdcutsvalmine[5][11]=0.2;
    rdcutsvalmine[6][11]=0.2;
    rdcutsvalmine[7][11]=-0.000023;
    rdcutsvalmine[8][11]=0.97;
    rdcutsvalmine[9][11]=0.3;
    rdcutsvalmine[10][11]=0.1;
    rdcutsvalmine[11][11]=0.45;
    rdcutsvalmine[12][11]=100.;
    rdcutsvalmine[13][11]=1;
    rdcutsvalmine[14][11]=0.999;
    rdcutsvalmine[15][11]=5.;

    //16-20 (ptbin11)
    rdcutsvalmine[0][12]=0.074;
    rdcutsvalmine[1][12]=0.8;
    rdcutsvalmine[2][12]=0.5625;
    rdcutsvalmine[3][12]=0.30;
    rdcutsvalmine[4][12]=0.30;
    rdcutsvalmine[5][12]=0.60;
    rdcutsvalmine[6][12]=0.60;
    rdcutsvalmine[7][12]=0.0002;
    rdcutsvalmine[8][12]=0.97;
    rdcutsvalmine[9][12]=0.3;
    rdcutsvalmine[10][12]=0.15;
    rdcutsvalmine[11][12]=0.5;
    rdcutsvalmine[12][12]=100.;
    rdcutsvalmine[13][12]=1;
    rdcutsvalmine[14][12]=0.998;
    rdcutsvalmine[15][12]=5.;


    //20-24
    rdcutsvalmine[0][13]=0.074;
    rdcutsvalmine[1][13]=0.8;
    rdcutsvalmine[2][13]=1.0;
    rdcutsvalmine[3][13]=0.30;
    rdcutsvalmine[4][13]=0.30;
    rdcutsvalmine[5][13]=0.60;
    rdcutsvalmine[6][13]=0.60;
    rdcutsvalmine[7][13]=0.0025;
    rdcutsvalmine[8][13]=0.8;
    rdcutsvalmine[9][13]=0.3;
    rdcutsvalmine[10][13]=0.15;
    rdcutsvalmine[11][13]=0.5;
    rdcutsvalmine[12][13]=100.;
    rdcutsvalmine[13][13]=1;
    rdcutsvalmine[14][13]=0.99;
    rdcutsvalmine[15][13]=0.;

    //>25 (ptbin14)
    rdcutsvalmine[0][14]=0.074;
    rdcutsvalmine[1][14]=0.8;
    rdcutsvalmine[2][14]=1.0;
    rdcutsvalmine[3][14]=0.30;
    rdcutsvalmine[4][14]=0.30;
    rdcutsvalmine[5][14]=0.50;
    rdcutsvalmine[6][14]=0.50;
    rdcutsvalmine[7][14]=0.0033;
    rdcutsvalmine[8][14]=0.7;
    rdcutsvalmine[9][14]=0.1;
    rdcutsvalmine[10][14]=0.1;
    rdcutsvalmine[11][14]=0.05;
    rdcutsvalmine[12][14]=100.;
    rdcutsvalmine[13][14]=1;
    rdcutsvalmine[14][14]=0.9;
    rdcutsvalmine[15][14]=0.0;
*/
/*
// Cuts 04_03_12 (very strong, then loosened  DCA, CTS, pt, d0 and CPAxy)  -> 040312.root
    //0-0.5
    rdcutsvalmine[0][0]=0.026;
    rdcutsvalmine[1][0]=0.022;
    rdcutsvalmine[2][0]=0.7;
    rdcutsvalmine[3][0]=0.21;
    rdcutsvalmine[4][0]=0.21;
    rdcutsvalmine[5][0]=0.05;
    rdcutsvalmine[6][0]=0.05;
    rdcutsvalmine[7][0]=0.0004;
    rdcutsvalmine[8][0]=0.6;
    rdcutsvalmine[9][0]=0.3;
    rdcutsvalmine[10][0]=0.1;
    rdcutsvalmine[11][0]=0.05;
    rdcutsvalmine[12][0]=100.;
    rdcutsvalmine[13][0]=0.5;
    rdcutsvalmine[14][0]=-1.;
    rdcutsvalmine[15][0]=0.;

    //0.5-1
    rdcutsvalmine[0][1]=0.026;
    rdcutsvalmine[1][1]=0.04;
    rdcutsvalmine[2][1]=0.8;
    rdcutsvalmine[3][1]=0.5;
    rdcutsvalmine[4][1]=0.5;
    rdcutsvalmine[5][1]=0.08;
    rdcutsvalmine[6][1]=0.08;
    rdcutsvalmine[7][1]=0.0003;
    rdcutsvalmine[8][1]=0.65;
    rdcutsvalmine[9][1]=0.3;
    rdcutsvalmine[10][1]=0.1;
    rdcutsvalmine[11][1]=0.05;
    rdcutsvalmine[12][1]=100.;
    rdcutsvalmine[13][1]=0.5;
    rdcutsvalmine[14][1]=-1.;
    rdcutsvalmine[15][1]=0.;

    //1-2
    rdcutsvalmine[0][2]=0.022;
    rdcutsvalmine[1][2]=0.04;
    rdcutsvalmine[2][2]=0.8;
    rdcutsvalmine[3][2]=0.7;
    rdcutsvalmine[4][2]=0.7;
    rdcutsvalmine[5][2]=0.08;
    rdcutsvalmine[6][2]=0.08;
    rdcutsvalmine[7][2]=0.00014;
    rdcutsvalmine[8][2]=0.75;
    rdcutsvalmine[9][2]=0.3;
    rdcutsvalmine[10][2]=0.1;
    rdcutsvalmine[11][2]=0.05;
    rdcutsvalmine[12][2]=100.;
    rdcutsvalmine[13][2]=0.5;
    rdcutsvalmine[14][2]=-1.;
    rdcutsvalmine[15][2]=0.;
    //2-3
    rdcutsvalmine[0][3]=0.024;
    rdcutsvalmine[1][3]=0.025;
    rdcutsvalmine[2][3]=0.8;
    rdcutsvalmine[3][3]=0.7;
    rdcutsvalmine[4][3]=0.7;
    rdcutsvalmine[5][3]=0.08;
    rdcutsvalmine[6][3]=0.08;
    rdcutsvalmine[7][3]=0.00014;
    rdcutsvalmine[8][3]=0.85;
    rdcutsvalmine[9][3]=0.3;
    rdcutsvalmine[10][3]=0.1;
    rdcutsvalmine[11][3]=0.05;
    rdcutsvalmine[12][3]=100.;
    rdcutsvalmine[13][3]=0.5;
    rdcutsvalmine[14][3]=-1.;
    rdcutsvalmine[15][3]=0.;
*/
/*
    //3-4
    rdcutsvalmine[0][4]=0.032;
    rdcutsvalmine[1][4]=0.03750;
    rdcutsvalmine[2][4]=1.0;
    rdcutsvalmine[3][4]=1.0;
    rdcutsvalmine[4][4]=1.0;
    rdcutsvalmine[5][4]=0.1;
    rdcutsvalmine[6][4]=0.1; ////////////////////////////////////////////////////
    rdcutsvalmine[7][4]=-0.0002;
    rdcutsvalmine[8][4]=0.99;
    rdcutsvalmine[9][4]=0.3;
    rdcutsvalmine[10][4]=0.15;
    rdcutsvalmine[11][4]=0.05;
    rdcutsvalmine[12][4]=0.5;
    rdcutsvalmine[13][4]=0.5;
    rdcutsvalmine[14][4]=0.995;
    rdcutsvalmine[15][4]=4.0;

    //4-5
    rdcutsvalmine[0][5]=0.032;
    rdcutsvalmine[1][5]=0.03750;
    rdcutsvalmine[2][5]=1.0;
    rdcutsvalmine[3][5]=1.0;
    rdcutsvalmine[4][5]=1.0;
    rdcutsvalmine[5][5]=0.075;
    rdcutsvalmine[6][5]=0.075;
    rdcutsvalmine[7][5]=-0.00025;
    rdcutsvalmine[8][5]=0.99;
    rdcutsvalmine[9][5]=0.15;
    rdcutsvalmine[10][5]=0.15;
    rdcutsvalmine[11][5]=0.05;
    rdcutsvalmine[12][5]=10;
    rdcutsvalmine[13][5]=0.5;
    rdcutsvalmine[14][5]=0.99;
    rdcutsvalmine[15][5]=5.0;

    //5-6 (ptbin6)
    rdcutsvalmine[0][6]=0.04;
    rdcutsvalmine[1][6]=0.03750;
    rdcutsvalmine[2][6]=1.0;
    rdcutsvalmine[3][6]=1.167;
    rdcutsvalmine[4][6]=1.167;
    rdcutsvalmine[5][6]=0.075;
    rdcutsvalmine[6][6]=0.075;
    rdcutsvalmine[7][6]=-0.000287;
    rdcutsvalmine[8][6]=0.99;
    rdcutsvalmine[9][6]=0.15;
    rdcutsvalmine[10][6]=0.15;
    rdcutsvalmine[11][6]=0.2;
    rdcutsvalmine[12][6]=10.;
    rdcutsvalmine[13][6]=0.5;
    rdcutsvalmine[14][6]=0.993;
    rdcutsvalmine[15][6]=6.5;

    //6-7
    rdcutsvalmine[0][7]=0.036;
    rdcutsvalmine[1][7]=0.03750;
    rdcutsvalmine[2][7]=1.0;
    rdcutsvalmine[3][7]=1.333;
    rdcutsvalmine[4][7]=1.333;
    rdcutsvalmine[5][7]=0.075;
    rdcutsvalmine[6][7]=0.075;
    rdcutsvalmine[7][7]=-0.0002;
    rdcutsvalmine[8][7]=0.98;
    rdcutsvalmine[9][7]=0.3;
    rdcutsvalmine[10][7]=0.1;
    rdcutsvalmine[11][7]=0.2;
    rdcutsvalmine[12][7]=100.;
    rdcutsvalmine[13][7]=0.5;
    rdcutsvalmine[14][7]=0.99;
    rdcutsvalmine[15][7]=4.0;

    //7-8
    rdcutsvalmine[0][8]=0.036;
    rdcutsvalmine[1][8]=0.03750;
    rdcutsvalmine[2][8]=1.0;
    rdcutsvalmine[3][8]=0.833;
    rdcutsvalmine[4][8]=0.833;
    rdcutsvalmine[5][8]=0.2;
    rdcutsvalmine[6][8]=0.2;
    rdcutsvalmine[7][8]=-0.000163;
    rdcutsvalmine[8][8]=0.96;
    rdcutsvalmine[9][8]=0.3;
    rdcutsvalmine[10][8]=0.15;
    rdcutsvalmine[11][8]=0.2;
    rdcutsvalmine[12][8]=100.;
    rdcutsvalmine[13][8]=0.5;
    rdcutsvalmine[14][8]=0.995;
    rdcutsvalmine[15][8]=3.;

    //8-10 (ptbin9)
    rdcutsvalmine[0][9]=0.055;
    rdcutsvalmine[1][9]=0.03750;
    rdcutsvalmine[2][9]=1.0;
    rdcutsvalmine[3][9]=1.;
    rdcutsvalmine[4][9]=1.;
    rdcutsvalmine[5][9]=0.2;
    rdcutsvalmine[6][9]=0.2;
    rdcutsvalmine[7][9]=-0.0001;
    rdcutsvalmine[8][9]=0.98;
    rdcutsvalmine[9][9]=0.3;
    rdcutsvalmine[10][9]=0.15;
    rdcutsvalmine[11][9]=0.34;
    rdcutsvalmine[12][9]=100.;
    rdcutsvalmine[13][9]=0.5;
    rdcutsvalmine[14][9]=0.992;
    rdcutsvalmine[15][9]=3.;

    //10-12
    rdcutsvalmine[0][10]=0.055;
    rdcutsvalmine[1][10]=0.03750;
    rdcutsvalmine[2][10]=1.0;
    rdcutsvalmine[3][10]=1.333;
    rdcutsvalmine[4][10]=1.333;
    rdcutsvalmine[5][10]=0.2;
    rdcutsvalmine[6][10]=0.2;
    rdcutsvalmine[7][10]=-0.000050;
    rdcutsvalmine[8][10]=0.97;
    rdcutsvalmine[9][10]=0.3;
    rdcutsvalmine[10][10]=0.15;
    rdcutsvalmine[11][10]=0.35;
    rdcutsvalmine[12][10]=100.;
    rdcutsvalmine[13][10]=0.5;
    rdcutsvalmine[14][10]=0.992;
    rdcutsvalmine[15][10]=5.;

    //12-16 (ptbin11)
    rdcutsvalmine[0][11]=0.074;
    rdcutsvalmine[1][11]=0.03750;
    rdcutsvalmine[2][11]=1.0;
    rdcutsvalmine[3][11]=1.30;
    rdcutsvalmine[4][11]=1.30;
    rdcutsvalmine[5][11]=0.2;
    rdcutsvalmine[6][11]=0.2;
    rdcutsvalmine[7][11]=-0.000023;
    rdcutsvalmine[8][11]=0.97;
    rdcutsvalmine[9][11]=0.3;
    rdcutsvalmine[10][11]=0.1;
    rdcutsvalmine[11][11]=0.45;
    rdcutsvalmine[12][11]=100.;
    rdcutsvalmine[13][11]=1;
    rdcutsvalmine[14][11]=0.996;
    rdcutsvalmine[15][11]=5.;

    //16-20  (ptbin12)
    rdcutsvalmine[0][12]=0.074;
    rdcutsvalmine[1][12]=0.8;
    rdcutsvalmine[2][12]=1.0;
    rdcutsvalmine[3][12]=0.30;
    rdcutsvalmine[4][12]=0.30;
    rdcutsvalmine[5][12]=0.60;
    rdcutsvalmine[6][12]=0.60;
    rdcutsvalmine[7][12]=0.0001;
    rdcutsvalmine[8][12]=0.96;
    rdcutsvalmine[9][12]=0.3;
    rdcutsvalmine[10][12]=0.15;
    rdcutsvalmine[11][12]=0.5;
    rdcutsvalmine[12][12]=100.;
    rdcutsvalmine[13][12]=1;
    rdcutsvalmine[14][12]=0.992;
    rdcutsvalmine[15][12]=4.;

    //20-24
    rdcutsvalmine[0][13]=0.074;
    rdcutsvalmine[1][13]=0.8;
    rdcutsvalmine[2][13]=1.0;
    rdcutsvalmine[3][13]=0.30;
    rdcutsvalmine[4][13]=0.30;
    rdcutsvalmine[5][13]=0.60;
    rdcutsvalmine[6][13]=0.60;
    rdcutsvalmine[7][13]=0.003;
    rdcutsvalmine[8][13]=0.7;
    rdcutsvalmine[9][13]=0.3;
    rdcutsvalmine[10][13]=0.15;
    rdcutsvalmine[11][13]=0.5;
    rdcutsvalmine[12][13]=100.;
    rdcutsvalmine[13][13]=1;
    rdcutsvalmine[14][13]=0.99;
    rdcutsvalmine[15][13]=0.;

    //>25 (ptbin14)
    rdcutsvalmine[0][14]=0.074;
    rdcutsvalmine[1][14]=0.8;
    rdcutsvalmine[2][14]=1.0;
    rdcutsvalmine[3][14]=0.30;
    rdcutsvalmine[4][14]=0.30;
    rdcutsvalmine[5][14]=0.50;
    rdcutsvalmine[6][14]=0.50;
    rdcutsvalmine[7][14]=0.033;
    rdcutsvalmine[8][14]=0.7;
    rdcutsvalmine[9][14]=0.1;
    rdcutsvalmine[10][14]=0.1;
    rdcutsvalmine[11][14]=0.05;
    rdcutsvalmine[12][14]=100.;
    rdcutsvalmine[13][14]=1;
    rdcutsvalmine[14][14]=0.9;
    rdcutsvalmine[15][14]=0.0;
*/
/*
// Cuts 01_03_12 (very strong, then loosened  DCA, CTS, pt, d0 and CPAxy) -> 010312.root
    //0-0.5
    rdcutsvalmine[0][0]=0.026;
    rdcutsvalmine[1][0]=0.022;
    rdcutsvalmine[2][0]=0.7;
    rdcutsvalmine[3][0]=0.21;
    rdcutsvalmine[4][0]=0.21;
    rdcutsvalmine[5][0]=0.05;
    rdcutsvalmine[6][0]=0.05;
    rdcutsvalmine[7][0]=0.0004;
    rdcutsvalmine[8][0]=0.6;
    rdcutsvalmine[9][0]=0.3;
    rdcutsvalmine[10][0]=0.1;
    rdcutsvalmine[11][0]=0.05;
    rdcutsvalmine[12][0]=100.;
    rdcutsvalmine[13][0]=0.5;
    rdcutsvalmine[14][0]=-1.;
    rdcutsvalmine[15][0]=0.;

    //0.5-1
    rdcutsvalmine[0][1]=0.026;
    rdcutsvalmine[1][1]=0.04;
    rdcutsvalmine[2][1]=0.8;
    rdcutsvalmine[3][1]=0.5;
    rdcutsvalmine[4][1]=0.5;
    rdcutsvalmine[5][1]=0.08;
    rdcutsvalmine[6][1]=0.08;
    rdcutsvalmine[7][1]=0.0003;
    rdcutsvalmine[8][1]=0.65;
    rdcutsvalmine[9][1]=0.3;
    rdcutsvalmine[10][1]=0.1;
    rdcutsvalmine[11][1]=0.05;
    rdcutsvalmine[12][1]=100.;
    rdcutsvalmine[13][1]=0.5;
    rdcutsvalmine[14][1]=-1.;
    rdcutsvalmine[15][1]=0.;

    //1-2
    rdcutsvalmine[0][2]=0.022;
    rdcutsvalmine[1][2]=0.04;
    rdcutsvalmine[2][2]=0.8;
    rdcutsvalmine[3][2]=0.7;
    rdcutsvalmine[4][2]=0.7;
    rdcutsvalmine[5][2]=0.08;
    rdcutsvalmine[6][2]=0.08;
    rdcutsvalmine[7][2]=0.00014;
    rdcutsvalmine[8][2]=0.75;
    rdcutsvalmine[9][2]=0.3;
    rdcutsvalmine[10][2]=0.1;
    rdcutsvalmine[11][2]=0.05;
    rdcutsvalmine[12][2]=100.;
    rdcutsvalmine[13][2]=0.5;
    rdcutsvalmine[14][2]=-1.;
    rdcutsvalmine[15][2]=0.;
    //2-3
    rdcutsvalmine[0][3]=0.024;
    rdcutsvalmine[1][3]=0.025;
    rdcutsvalmine[2][3]=0.8;
    rdcutsvalmine[3][3]=0.7;
    rdcutsvalmine[4][3]=0.7;
    rdcutsvalmine[5][3]=0.08;
    rdcutsvalmine[6][3]=0.08;
    rdcutsvalmine[7][3]=0.00014;
    rdcutsvalmine[8][3]=0.85;
    rdcutsvalmine[9][3]=0.3;
    rdcutsvalmine[10][3]=0.1;
    rdcutsvalmine[11][3]=0.05;
    rdcutsvalmine[12][3]=100.;
    rdcutsvalmine[13][3]=0.5;
    rdcutsvalmine[14][3]=-1.;
    rdcutsvalmine[15][3]=0.;


    //3-4
    rdcutsvalmine[0][4]=0.032;
    rdcutsvalmine[1][4]=0.03750;
    rdcutsvalmine[2][4]=1.0;
    rdcutsvalmine[3][4]=1.0;
    rdcutsvalmine[4][4]=1.0;
    rdcutsvalmine[5][4]=0.1;
    rdcutsvalmine[6][4]=0.1; ////////////////////////////////////////////////////
    rdcutsvalmine[7][4]=-0.000285;
    rdcutsvalmine[8][4]=0.99;
    rdcutsvalmine[9][4]=0.3;
    rdcutsvalmine[10][4]=0.15;
    rdcutsvalmine[11][4]=0.05;
    rdcutsvalmine[12][4]=0.5;
    rdcutsvalmine[13][4]=0.5;
    rdcutsvalmine[14][4]=0.995;
    rdcutsvalmine[15][4]=7.0;

    //4-5
    rdcutsvalmine[0][5]=0.032;
    rdcutsvalmine[1][5]=0.03750;
    rdcutsvalmine[2][5]=1.0;
    rdcutsvalmine[3][5]=1.0;
    rdcutsvalmine[4][5]=1.0;
    rdcutsvalmine[5][5]=0.075;
    rdcutsvalmine[6][5]=0.075;
    rdcutsvalmine[7][5]=-0.000362;
    rdcutsvalmine[8][5]=0.99;
    rdcutsvalmine[9][5]=0.15;
    rdcutsvalmine[10][5]=0.15;
    rdcutsvalmine[11][5]=0.05;
    rdcutsvalmine[12][5]=10;
    rdcutsvalmine[13][5]=0.5;
    rdcutsvalmine[14][5]=0.99;
    rdcutsvalmine[15][5]=8.0;

    //5-6 (ptbin6)
    rdcutsvalmine[0][6]=0.04;
    rdcutsvalmine[1][6]=0.03750;
    rdcutsvalmine[2][6]=1.0;
    rdcutsvalmine[3][6]=1.167;
    rdcutsvalmine[4][6]=1.167;
    rdcutsvalmine[5][6]=0.075;
    rdcutsvalmine[6][6]=0.075;
    rdcutsvalmine[7][6]=-0.000287;
    rdcutsvalmine[8][6]=0.99;
    rdcutsvalmine[9][6]=0.15;
    rdcutsvalmine[10][6]=0.15;
    rdcutsvalmine[11][6]=0.2;
    rdcutsvalmine[12][6]=10.;
    rdcutsvalmine[13][6]=0.5;
    rdcutsvalmine[14][6]=0.993;
    rdcutsvalmine[15][6]=6.5;
*/
/*    //6-7
    rdcutsvalmine[0][7]=0.036;
    rdcutsvalmine[1][7]=0.03750;
    rdcutsvalmine[2][7]=1.0;
    rdcutsvalmine[3][7]=1.333;
    rdcutsvalmine[4][7]=1.333;
    rdcutsvalmine[5][7]=0.075;
    rdcutsvalmine[6][7]=0.075;
    rdcutsvalmine[7][7]=-0.000243;
    rdcutsvalmine[8][7]=0.986;
    rdcutsvalmine[9][7]=0.3;
    rdcutsvalmine[10][7]=0.1;
    rdcutsvalmine[11][7]=0.2;
    rdcutsvalmine[12][7]=100.;
    rdcutsvalmine[13][7]=0.5;
    rdcutsvalmine[14][7]=0.992;
    rdcutsvalmine[15][7]=4.0;

    //7-8
    rdcutsvalmine[0][8]=0.036;
    rdcutsvalmine[1][8]=0.03750;
    rdcutsvalmine[2][8]=1.0;
    rdcutsvalmine[3][8]=0.833;
    rdcutsvalmine[4][8]=0.833;
    rdcutsvalmine[5][8]=0.2;
    rdcutsvalmine[6][8]=0.2;
    rdcutsvalmine[7][8]=-0.000163;
    rdcutsvalmine[8][8]=0.967;
    rdcutsvalmine[9][8]=0.3;
    rdcutsvalmine[10][8]=0.15;
    rdcutsvalmine[11][8]=0.2;
    rdcutsvalmine[12][8]=100.;
    rdcutsvalmine[13][8]=0.5;
    rdcutsvalmine[14][8]=0.995;
    rdcutsvalmine[15][8]=9.;
*/
/*
    //8-10 (ptbin9)
    rdcutsvalmine[0][9]=0.055;
    rdcutsvalmine[1][9]=0.03750;
    rdcutsvalmine[2][9]=1.0;
    rdcutsvalmine[3][9]=1.;
    rdcutsvalmine[4][9]=1.;
    rdcutsvalmine[5][9]=0.2;
    rdcutsvalmine[6][9]=0.2;
    rdcutsvalmine[7][9]=-0.000163;
    rdcutsvalmine[8][9]=0.987;
    rdcutsvalmine[9][9]=0.3;
    rdcutsvalmine[10][9]=0.15;
    rdcutsvalmine[11][9]=0.34;
    rdcutsvalmine[12][9]=100.;
    rdcutsvalmine[13][9]=0.5;
    rdcutsvalmine[14][9]=0.996;
    rdcutsvalmine[15][9]=5;

    //10-12
    rdcutsvalmine[0][10]=0.055;
    rdcutsvalmine[1][10]=0.03750;
    rdcutsvalmine[2][10]=1.0;
    rdcutsvalmine[3][10]=1.333;
    rdcutsvalmine[4][10]=1.333;
    rdcutsvalmine[5][10]=0.2;
    rdcutsvalmine[6][10]=0.2;
    rdcutsvalmine[7][10]=-0.000100;
    rdcutsvalmine[8][10]=0.97;
    rdcutsvalmine[9][10]=0.3;
    rdcutsvalmine[10][10]=0.15;
    rdcutsvalmine[11][10]=0.35;
    rdcutsvalmine[12][10]=100.;
    rdcutsvalmine[13][10]=0.5;
    rdcutsvalmine[14][10]=0.996;
    rdcutsvalmine[15][10]=7.;
*/
/*
    //12-16 (ptbin11)
    rdcutsvalmine[0][11]=0.074;
    rdcutsvalmine[1][11]=0.03750;
    rdcutsvalmine[2][11]=1.0;
    rdcutsvalmine[3][11]=1.30;
    rdcutsvalmine[4][11]=1.30;
    rdcutsvalmine[5][11]=0.2;
    rdcutsvalmine[6][11]=0.2;
    rdcutsvalmine[7][11]=-0.000023;
    rdcutsvalmine[8][11]=0.97;
    rdcutsvalmine[9][11]=0.3;
    rdcutsvalmine[10][11]=0.1;
    rdcutsvalmine[11][11]=0.45;
    rdcutsvalmine[12][11]=100.;
    rdcutsvalmine[13][11]=1;
    rdcutsvalmine[14][11]=0.996;
    rdcutsvalmine[15][11]=5.;
*/

/*
    //16-20  (ptbin12)
    //rdcutsvalmine[0][12]=0.074;
    //rdcutsvalmine[1][12]=0.8;
    //rdcutsvalmine[2][12]=1.0;
    //rdcutsvalmine[3][12]=0.30;
    //rdcutsvalmine[4][12]=0.30;
    //rdcutsvalmine[5][12]=0.60;
    //rdcutsvalmine[6][12]=0.60;
    //rdcutsvalmine[7][12]=0.0002;
    //rdcutsvalmine[8][12]=0.97;
    //rdcutsvalmine[9][12]=0.3;
    //rdcutsvalmine[10][12]=0.15;
    //rdcutsvalmine[11][12]=0.5;
    //rdcutsvalmine[12][12]=100.;
    //rdcutsvalmine[13][12]=1;
    //rdcutsvalmine[14][12]=0.996;
    //rdcutsvalmine[15][12]=5.;
    rdcutsvalmine[0][12]=0.074;
    rdcutsvalmine[1][12]=0.8;
    rdcutsvalmine[2][12]=1.0;
    rdcutsvalmine[3][12]=0.30;
    rdcutsvalmine[4][12]=0.30;
    rdcutsvalmine[5][12]=0.60;
    rdcutsvalmine[6][12]=0.60;
    rdcutsvalmine[7][12]=0.0025;
    rdcutsvalmine[8][12]=0.8;
    rdcutsvalmine[9][12]=0.3;
    rdcutsvalmine[10][12]=0.15;
    rdcutsvalmine[11][12]=0.5;
    rdcutsvalmine[12][12]=100.;
    rdcutsvalmine[13][12]=1;
    rdcutsvalmine[14][12]=0.99;
    rdcutsvalmine[15][12]=0.;


    //20-24
    rdcutsvalmine[0][13]=0.074;
    rdcutsvalmine[1][13]=0.8;
    rdcutsvalmine[2][13]=1.0;
    rdcutsvalmine[3][13]=0.30;
    rdcutsvalmine[4][13]=0.30;
    rdcutsvalmine[5][13]=0.60;
    rdcutsvalmine[6][13]=0.60;
    rdcutsvalmine[7][13]=0.0025;
    rdcutsvalmine[8][13]=0.8;
    rdcutsvalmine[9][13]=0.3;
    rdcutsvalmine[10][13]=0.15;
    rdcutsvalmine[11][13]=0.5;
    rdcutsvalmine[12][13]=100.;
    rdcutsvalmine[13][13]=1;
    rdcutsvalmine[14][13]=0.99;
    rdcutsvalmine[15][13]=0.;

    //>25 (ptbin14)
    rdcutsvalmine[0][14]=0.074;
    rdcutsvalmine[1][14]=0.8;
    rdcutsvalmine[2][14]=1.0;
    rdcutsvalmine[3][14]=0.30;
    rdcutsvalmine[4][14]=0.30;
    rdcutsvalmine[5][14]=0.50;
    rdcutsvalmine[6][14]=0.50;
    rdcutsvalmine[7][14]=0.0033;
    rdcutsvalmine[8][14]=0.7;
    rdcutsvalmine[9][14]=0.1;
    rdcutsvalmine[10][14]=0.1;
    rdcutsvalmine[11][14]=0.05;
    rdcutsvalmine[12][14]=100.;
    rdcutsvalmine[13][14]=1;
    rdcutsvalmine[14][14]=0.9;
    rdcutsvalmine[15][14]=0.0;
*/
/*
// Cuts 27_02_12 (very strong) 	-> 290212.root
    //0-0.5
    rdcutsvalmine[0][0]=0.026;
    rdcutsvalmine[1][0]=0.022;
    rdcutsvalmine[2][0]=0.7;
    rdcutsvalmine[3][0]=0.21;
    rdcutsvalmine[4][0]=0.21;
    rdcutsvalmine[5][0]=0.05;
    rdcutsvalmine[6][0]=0.05;
    rdcutsvalmine[7][0]=0.0004;
    rdcutsvalmine[8][0]=0.6;
    rdcutsvalmine[9][0]=0.3;
    rdcutsvalmine[10][0]=0.1;
    rdcutsvalmine[11][0]=0.05;
    rdcutsvalmine[12][0]=100.;
    rdcutsvalmine[13][0]=0.5;
    rdcutsvalmine[14][0]=-1.;
    rdcutsvalmine[15][0]=0.;

    //0.5-1
    rdcutsvalmine[0][1]=0.026;
    rdcutsvalmine[1][1]=0.04;
    rdcutsvalmine[2][1]=0.8;
    rdcutsvalmine[3][1]=0.5;
    rdcutsvalmine[4][1]=0.5;
    rdcutsvalmine[5][1]=0.08;
    rdcutsvalmine[6][1]=0.08;
    rdcutsvalmine[7][1]=0.0003;
    rdcutsvalmine[8][1]=0.65;
    rdcutsvalmine[9][1]=0.3;
    rdcutsvalmine[10][1]=0.1;
    rdcutsvalmine[11][1]=0.05;
    rdcutsvalmine[12][1]=100.;
    rdcutsvalmine[13][1]=0.5;
    rdcutsvalmine[14][1]=-1.;
    rdcutsvalmine[15][1]=0.;

    //1-2
    rdcutsvalmine[0][2]=0.022;
    rdcutsvalmine[1][2]=0.04;
    rdcutsvalmine[2][2]=0.8;
    rdcutsvalmine[3][2]=0.7;
    rdcutsvalmine[4][2]=0.7;
    rdcutsvalmine[5][2]=0.08;
    rdcutsvalmine[6][2]=0.08;
    rdcutsvalmine[7][2]=0.00014;
    rdcutsvalmine[8][2]=0.75;
    rdcutsvalmine[9][2]=0.3;
    rdcutsvalmine[10][2]=0.1;
    rdcutsvalmine[11][2]=0.05;
    rdcutsvalmine[12][2]=100.;
    rdcutsvalmine[13][2]=0.5;
    rdcutsvalmine[14][2]=-1.;
    rdcutsvalmine[15][2]=0.;
    //2-3
    rdcutsvalmine[0][3]=0.024;
    rdcutsvalmine[1][3]=0.025;
    rdcutsvalmine[2][3]=0.8;
    rdcutsvalmine[3][3]=0.7;
    rdcutsvalmine[4][3]=0.7;
    rdcutsvalmine[5][3]=0.08;
    rdcutsvalmine[6][3]=0.08;
    rdcutsvalmine[7][3]=0.00014;
    rdcutsvalmine[8][3]=0.85;
    rdcutsvalmine[9][3]=0.3;
    rdcutsvalmine[10][3]=0.1;
    rdcutsvalmine[11][3]=0.05;
    rdcutsvalmine[12][3]=100.;
    rdcutsvalmine[13][3]=0.5;
    rdcutsvalmine[14][3]=-1.;
    rdcutsvalmine[15][3]=0.;

*/
/*
    //3-4
    rdcutsvalmine[0][4]=0.032;
    rdcutsvalmine[1][4]=0.01250;	//0.1;//0.025;
    rdcutsvalmine[2][4]=1.0;			//1.0;//0.7;
    rdcutsvalmine[3][4]=1.333;		//0.5;//1.1;
    rdcutsvalmine[4][4]=1.0;			//0.5;//1.1;
    rdcutsvalmine[5][4]=0.1;			//0.1;//0.08;
    rdcutsvalmine[6][4]=0.05833;		//0.1;//0.08;
    rdcutsvalmine[7][4]=-0.000285;		//-0.00017;
    rdcutsvalmine[8][4]=0.99;			//0.98;
    rdcutsvalmine[9][4]=0.3;
    rdcutsvalmine[10][4]=0.15;
    rdcutsvalmine[11][4]=0.05;
    rdcutsvalmine[12][4]=0.5;
    rdcutsvalmine[13][4]=0.5;
    rdcutsvalmine[14][4]=0.998;			//0.99;
    rdcutsvalmine[15][4]=7.0;			//5.;
*/
/*
    //4-5
    rdcutsvalmine[0][5]=0.032;
    rdcutsvalmine[1][5]=0.03750;//0.1;//0.025;
    rdcutsvalmine[2][5]=0.5625;//1.0;//0.7;
    rdcutsvalmine[3][5]=1.0;//0.5;//1.2;
    rdcutsvalmine[4][5]=1.0;//0.5;//1.2;
    rdcutsvalmine[5][5]=0.075;//0.1;//0.09;
    rdcutsvalmine[6][5]=0.0583;//0.1;//0.09;
    rdcutsvalmine[7][5]=-0.000362;//-0.00017;
    rdcutsvalmine[8][5]=0.99;//0.98;
    rdcutsvalmine[9][5]=0.15;
    rdcutsvalmine[10][5]=0.15;
    rdcutsvalmine[11][5]=0.05;
    rdcutsvalmine[12][5]=10;
    rdcutsvalmine[13][5]=0.5;
    rdcutsvalmine[14][5]=0.99;//0.99;
    rdcutsvalmine[15][5]=8.0;//0.;

    //5-6 (ptbin6)
    rdcutsvalmine[0][6]=0.04;
    rdcutsvalmine[1][6]=0.025;//0.1;//0.027;
    rdcutsvalmine[2][6]=1.0;//1.0;//0.8;
    rdcutsvalmine[3][6]=1.33;//0.5;//1.1;
    rdcutsvalmine[4][6]=1.167;//0.5;//1.1;
    rdcutsvalmine[5][6]=0.075;//0.2;//0.1;
    rdcutsvalmine[6][6]=0.075;//0.2;//0.1;
    rdcutsvalmine[7][6]=-0.000287;//-0.00016;
    rdcutsvalmine[8][6]=0.99;//0.98;
    rdcutsvalmine[9][6]=0.15;
    rdcutsvalmine[10][6]=0.15;
    rdcutsvalmine[11][6]=0.2;
    rdcutsvalmine[12][6]=10.;
    rdcutsvalmine[13][6]=0.5;
    rdcutsvalmine[14][6]=0.993;//0.99;
    rdcutsvalmine[15][6]=6.5;//0;

    //6-7
    rdcutsvalmine[0][7]=0.036;
    rdcutsvalmine[1][7]=0.03750;//0.1;//0.027;
    rdcutsvalmine[2][7]=1.0;//1.0;//0.8;
    rdcutsvalmine[3][7]=1.333;//0.5;//0.9;
    rdcutsvalmine[4][7]=1.333;//0.5;//0.9;
    rdcutsvalmine[5][7]=0.075;//0.2;//0.1;
    rdcutsvalmine[6][7]=0.075;//0.2;//0.1;
    rdcutsvalmine[7][7]=-0.000243;//-0.00019;
    rdcutsvalmine[8][7]=0.986;//0.98;
    rdcutsvalmine[9][7]=0.3;
    rdcutsvalmine[10][7]=0.1;
    rdcutsvalmine[11][7]=0.2;
    rdcutsvalmine[12][7]=100.;
    rdcutsvalmine[13][7]=0.5;
    rdcutsvalmine[14][7]=0.992;//0.99;
    rdcutsvalmine[15][7]=4.0;//0.;

    //7-8
    rdcutsvalmine[0][8]=0.036;
    rdcutsvalmine[1][8]=0.01250;//0.1;//0.027;
    rdcutsvalmine[2][8]=0.5625;//1.0;
    rdcutsvalmine[3][8]=0.833;//0.5;//0.9;
    rdcutsvalmine[4][8]=1.333;//0.5;//0.9;
    rdcutsvalmine[5][8]=0.075;//0.2;//0.1;
    rdcutsvalmine[6][8]=0.2;//0.2;//0.1;
    rdcutsvalmine[7][8]=-0.000163;//-0.00012;
    rdcutsvalmine[8][8]=0.967;//0.95;
    rdcutsvalmine[9][8]=0.3;
    rdcutsvalmine[10][8]=0.15;
    rdcutsvalmine[11][8]=0.2;
    rdcutsvalmine[12][8]=100.;
    rdcutsvalmine[13][8]=0.5;
    rdcutsvalmine[14][8]=0.997;//0.995;
    rdcutsvalmine[15][8]=9.;//5.;

    //8-10 (ptbin9)
    rdcutsvalmine[0][9]=0.055;
    rdcutsvalmine[1][9]=0.03750;//0.1;//0.03;
    rdcutsvalmine[2][9]=1.0;//1.0;
    rdcutsvalmine[3][9]=1.33;//0.5;//0.9;
    rdcutsvalmine[4][9]=1.1667;//0.5;//0.9;
    rdcutsvalmine[5][9]=0.2;//0.2;//0.15;
    rdcutsvalmine[6][9]=0.2;//0.2;//0.15;
    rdcutsvalmine[7][9]=-0.000163;//-0.00008;
    rdcutsvalmine[8][9]=0.987;//0.96;
    rdcutsvalmine[9][9]=0.3;
    rdcutsvalmine[10][9]=0.15;
    rdcutsvalmine[11][9]=0.34;
    rdcutsvalmine[12][9]=100.;
    rdcutsvalmine[13][9]=0.5;
    rdcutsvalmine[14][9]=0.998;//0.995;
    rdcutsvalmine[15][9]=5;//0;

    //10-12
    rdcutsvalmine[0][10]=0.055;
    rdcutsvalmine[1][10]=0.01250;//0.1;//0.03;
    rdcutsvalmine[2][10]=0.625;//1.0;
    rdcutsvalmine[3][10]=1.333;//0.5;//1.0;
    rdcutsvalmine[4][10]=1.333;//0.5;//1.0;
    rdcutsvalmine[5][10]=0.2;//0.2;//0.15;
    rdcutsvalmine[6][10]=0.2;//0.2;//0.15;
    rdcutsvalmine[7][10]=-0.000100;//-0.00000;//-0.00001;
    rdcutsvalmine[8][10]=0.97;//0.97;//0.92;
    rdcutsvalmine[9][10]=0.3;
    rdcutsvalmine[10][10]=0.15;
    rdcutsvalmine[11][10]=0.35;
    rdcutsvalmine[12][10]=100.;
    rdcutsvalmine[13][10]=0.5;
    rdcutsvalmine[14][10]=0.998;//0.992;
    rdcutsvalmine[15][10]=7.;//4.;

    //12-16 (ptbin11)
    rdcutsvalmine[0][11]=0.074;
    rdcutsvalmine[1][11]=0.01250;//0.1;//0.07;
    rdcutsvalmine[2][11]=0.5625;//1.0;
    rdcutsvalmine[3][11]=1.30;//0.3;
    rdcutsvalmine[4][11]=1.30;//0.3;
    rdcutsvalmine[5][11]=0.2;//0.2;//0.15;
    rdcutsvalmine[6][11]=0.2;//0.2;//0.15;
    rdcutsvalmine[7][11]=-0.000023;//-0.0000001;//0.0001;
    rdcutsvalmine[8][11]=0.97;//0.94;//0.87;
    rdcutsvalmine[9][11]=0.3;
    rdcutsvalmine[10][11]=0.1;
    rdcutsvalmine[11][11]=0.45;
    rdcutsvalmine[12][11]=100.;
    rdcutsvalmine[13][11]=1;
    rdcutsvalmine[14][11]=0.999;//0.995;
    rdcutsvalmine[15][11]=5.;//2;

    //16-20  (ptbin12)
    rdcutsvalmine[0][12]=0.074;
    rdcutsvalmine[1][12]=0.8;//0.8;
    rdcutsvalmine[2][12]=0.56250;//1.0;
    rdcutsvalmine[3][12]=0.30;//0.3;
    rdcutsvalmine[4][12]=0.30;//0.3;
    rdcutsvalmine[5][12]=0.60;//0.6;
    rdcutsvalmine[6][12]=0.60;//0.6;
    rdcutsvalmine[7][12]=0.0002;//0.0002;//0.06;
    rdcutsvalmine[8][12]=0.97;//0.93;//0.1;
    rdcutsvalmine[9][12]=0.3;
    rdcutsvalmine[10][12]=0.15;
    rdcutsvalmine[11][12]=0.5;
    rdcutsvalmine[12][12]=100.;
    rdcutsvalmine[13][12]=1;
    rdcutsvalmine[14][12]=0.998;//0.992;
    rdcutsvalmine[15][12]=5.;//0.;

    //20-24
    rdcutsvalmine[0][13]=0.074;
    rdcutsvalmine[1][13]=0.8;//0.8;
    rdcutsvalmine[2][13]=1.0;//1.0;
    rdcutsvalmine[3][13]=0.30;//0.3;
    rdcutsvalmine[4][13]=0.30;//0.3;
    rdcutsvalmine[5][13]=0.60;//0.6;
    rdcutsvalmine[6][13]=0.60;//0.6;
    rdcutsvalmine[7][13]=0.0025;//0.0025;//0.01;
    rdcutsvalmine[8][13]=0.8;//0.8;//0.1;
    rdcutsvalmine[9][13]=0.3;
    rdcutsvalmine[10][13]=0.15;
    rdcutsvalmine[11][13]=0.5;
    rdcutsvalmine[12][13]=100.;
    rdcutsvalmine[13][13]=1;
    rdcutsvalmine[14][13]=0.99;//0.99;
    rdcutsvalmine[15][13]=0.;//0.;

    //>25 (ptbin14)
    rdcutsvalmine[0][14]=0.074;
    rdcutsvalmine[1][14]=0.6;//0.6;
    rdcutsvalmine[2][14]=0.6875;//1.0;
    rdcutsvalmine[3][14]=0.30;//.4;
    rdcutsvalmine[4][14]=0.30;//.4;
    rdcutsvalmine[5][14]=0.50;//0.5;
    rdcutsvalmine[6][14]=0.50;//0.5;
    rdcutsvalmine[7][14]=0.0033;//0.033;//0.1;
    rdcutsvalmine[8][14]=0.7;//0.7;//0.7;
    rdcutsvalmine[9][14]=0.1;
    rdcutsvalmine[10][14]=0.1;
    rdcutsvalmine[11][14]=0.05;
    rdcutsvalmine[12][14]=100.;
    rdcutsvalmine[13][14]=1;
    rdcutsvalmine[14][14]=0.9;//0.9;
    rdcutsvalmine[15][14]=0.0;//0.;
*/

/*
//  Cuts 20-02-12		-> 220212.root
    //0-0.5
    rdcutsvalmine[0][0]=0.026;
    rdcutsvalmine[1][0]=0.022;
    rdcutsvalmine[2][0]=0.7;
    rdcutsvalmine[3][0]=0.21;
    rdcutsvalmine[4][0]=0.21;
    rdcutsvalmine[5][0]=0.05;
    rdcutsvalmine[6][0]=0.05;
    rdcutsvalmine[7][0]=0.0004;
    rdcutsvalmine[8][0]=0.6;
    rdcutsvalmine[9][0]=0.3;
    rdcutsvalmine[10][0]=0.1;
    rdcutsvalmine[11][0]=0.05;
    rdcutsvalmine[12][0]=100.;
    rdcutsvalmine[13][0]=0.5;
    rdcutsvalmine[14][0]=-1.;
    rdcutsvalmine[15][0]=0.;

    //0.5-1
    rdcutsvalmine[0][1]=0.026;
    rdcutsvalmine[1][1]=0.04;
    rdcutsvalmine[2][1]=0.8;
    rdcutsvalmine[3][1]=0.5;
    rdcutsvalmine[4][1]=0.5;
    rdcutsvalmine[5][1]=0.08;
    rdcutsvalmine[6][1]=0.08;
    rdcutsvalmine[7][1]=0.0003;
    rdcutsvalmine[8][1]=0.65;
    rdcutsvalmine[9][1]=0.3;
    rdcutsvalmine[10][1]=0.1;
    rdcutsvalmine[11][1]=0.05;
    rdcutsvalmine[12][1]=100.;
    rdcutsvalmine[13][1]=0.5;
    rdcutsvalmine[14][1]=-1.;
    rdcutsvalmine[15][1]=0.;

    //1-2
    rdcutsvalmine[0][2]=0.022;
    rdcutsvalmine[1][2]=0.04;
    rdcutsvalmine[2][2]=0.8;
    rdcutsvalmine[3][2]=0.7;
    rdcutsvalmine[4][2]=0.7;
    rdcutsvalmine[5][2]=0.08;
    rdcutsvalmine[6][2]=0.08;
    rdcutsvalmine[7][2]=0.00014;
    rdcutsvalmine[8][2]=0.75;
    rdcutsvalmine[9][2]=0.3;
    rdcutsvalmine[10][2]=0.1;
    rdcutsvalmine[11][2]=0.05;
    rdcutsvalmine[12][2]=100.;
    rdcutsvalmine[13][2]=0.5;
    rdcutsvalmine[14][2]=-1.;
    rdcutsvalmine[15][2]=0.;
    //2-3
    rdcutsvalmine[0][3]=0.024;
    rdcutsvalmine[1][3]=0.025;
    rdcutsvalmine[2][3]=0.8;
    rdcutsvalmine[3][3]=0.7;
    rdcutsvalmine[4][3]=0.7;
    rdcutsvalmine[5][3]=0.08;
    rdcutsvalmine[6][3]=0.08;
    rdcutsvalmine[7][3]=0.00014;
    rdcutsvalmine[8][3]=0.85;
    rdcutsvalmine[9][3]=0.3;
    rdcutsvalmine[10][3]=0.1;
    rdcutsvalmine[11][3]=0.05;
    rdcutsvalmine[12][3]=100.;
    rdcutsvalmine[13][3]=0.5;
    rdcutsvalmine[14][3]=-1.;
    rdcutsvalmine[15][3]=0.;


    //3-4
    rdcutsvalmine[0][4]=0.032;
    rdcutsvalmine[1][4]=0.025; // 0.8
    rdcutsvalmine[2][4]=0.7;
    rdcutsvalmine[3][4]=1.1;
    rdcutsvalmine[4][4]=1.1;
    rdcutsvalmine[5][4]=0.08;
    rdcutsvalmine[6][4]=0.08;
    rdcutsvalmine[7][4]=-0.00027;
    rdcutsvalmine[8][4]=0.992;
    rdcutsvalmine[9][4]=0.3;
    rdcutsvalmine[10][4]=0.15;
    rdcutsvalmine[11][4]=0.05;
    rdcutsvalmine[12][4]=0.5;
    rdcutsvalmine[13][4]=0.5;
    rdcutsvalmine[14][4]=0.998;
    rdcutsvalmine[15][4]=9.75;
*/
/*
    //4-5
    rdcutsvalmine[0][5]=0.032;
    rdcutsvalmine[1][5]=0.025; //0.1
    rdcutsvalmine[2][5]=0.7;
    rdcutsvalmine[3][5]=1.2;
    rdcutsvalmine[4][5]=1.2;
    rdcutsvalmine[5][5]=0.09;
    rdcutsvalmine[6][5]=0.09;
    rdcutsvalmine[7][5]=-0.000337;
    rdcutsvalmine[8][5]=0.992;
    rdcutsvalmine[9][5]=0.15;
    rdcutsvalmine[10][5]=0.15;
    rdcutsvalmine[11][5]=0.05;
    rdcutsvalmine[12][5]=10;
    rdcutsvalmine[13][5]=0.5;
    rdcutsvalmine[14][5]=0.996;
    rdcutsvalmine[15][5]=8.5;

    //5-6
    rdcutsvalmine[0][6]=0.04;
    rdcutsvalmine[1][6]=0.027;
    rdcutsvalmine[2][6]=0.8;
    rdcutsvalmine[3][6]=1.1;
    rdcutsvalmine[4][6]=1.1;
    rdcutsvalmine[5][6]=0.1;
    rdcutsvalmine[6][6]=0.1;
    rdcutsvalmine[7][6]=-0.0003;
    rdcutsvalmine[8][6]=0.99;
    rdcutsvalmine[9][6]=0.15;
    rdcutsvalmine[10][6]=0.15;
    rdcutsvalmine[11][6]=0.2;
    rdcutsvalmine[12][6]=10.;
    rdcutsvalmine[13][6]=0.5;
    rdcutsvalmine[14][6]=0.995;
    rdcutsvalmine[15][6]=6.0;
*/
/*
    //6-7
    rdcutsvalmine[0][7]=0.036;
    rdcutsvalmine[1][7]=0.027;
    rdcutsvalmine[2][7]=0.8;
    rdcutsvalmine[3][7]=0.9;
    rdcutsvalmine[4][7]=0.9;
    rdcutsvalmine[5][7]=0.1;
    rdcutsvalmine[6][7]=0.1;
    rdcutsvalmine[7][7]=-0.000282;
    rdcutsvalmine[8][7]=0.994;
    rdcutsvalmine[9][7]=0.3;
    rdcutsvalmine[10][7]=0.1;
    rdcutsvalmine[11][7]=0.2;
    rdcutsvalmine[12][7]=100.;
    rdcutsvalmine[13][7]=0.5;
    rdcutsvalmine[14][7]=0.996;
    rdcutsvalmine[15][7]=5.8;

    //7-8
    rdcutsvalmine[0][8]=0.036;
    rdcutsvalmine[1][8]=0.027;
    rdcutsvalmine[2][8]=1;
    rdcutsvalmine[3][8]=0.9;
    rdcutsvalmine[4][8]=0.9;
    rdcutsvalmine[5][8]=0.1;
    rdcutsvalmine[6][8]=0.1;
    rdcutsvalmine[7][8]=-0.00019;
    rdcutsvalmine[8][8]=0.97;
    rdcutsvalmine[9][8]=0.3;
    rdcutsvalmine[10][8]=0.15;
    rdcutsvalmine[11][8]=0.2;
    rdcutsvalmine[12][8]=100.;
    rdcutsvalmine[13][8]=0.5;
    rdcutsvalmine[14][8]=0.998;
    rdcutsvalmine[15][8]=7.5;

    //8-10
    rdcutsvalmine[0][9]=0.055;
    rdcutsvalmine[1][9]=0.03;
    rdcutsvalmine[2][9]=1.0;
    rdcutsvalmine[3][9]=0.9;
    rdcutsvalmine[4][9]=0.9;
    rdcutsvalmine[5][9]=0.15;
    rdcutsvalmine[6][9]=0.15;
    rdcutsvalmine[7][9]=-0.000138;//-0.00001;
    rdcutsvalmine[8][9]=0.985;//0.92;
    rdcutsvalmine[9][9]=0.3;
    rdcutsvalmine[10][9]=0.15;
    rdcutsvalmine[11][9]=0.34;
    rdcutsvalmine[12][9]=100.;
    rdcutsvalmine[13][9]=0.5;
     rdcutsvalmine[14][9]=0.998;
    rdcutsvalmine[15][9]=6.5;

    //10-12
    rdcutsvalmine[0][10]=0.055;
    rdcutsvalmine[1][10]=0.03;
    rdcutsvalmine[2][10]=1.0;
    rdcutsvalmine[3][10]=1.0;
    rdcutsvalmine[4][10]=1.0;
    rdcutsvalmine[5][10]=0.15;
    rdcutsvalmine[6][10]=0.15;
    rdcutsvalmine[7][10]=-0.00004;//-0.00001;
    rdcutsvalmine[8][10]=0.99;//0.92;
    rdcutsvalmine[9][10]=0.3;
    rdcutsvalmine[10][10]=0.15;
    rdcutsvalmine[11][10]=0.35;
    rdcutsvalmine[12][10]=100.;
    rdcutsvalmine[13][10]=0.5;
     rdcutsvalmine[14][10]=0.996;
    rdcutsvalmine[15][10]=5.5;

    //12-16
    rdcutsvalmine[0][11]=0.074;
    rdcutsvalmine[1][11]=0.07;
    rdcutsvalmine[2][11]=1.0;
    rdcutsvalmine[3][11]=0.3;
    rdcutsvalmine[4][11]=0.3;
    rdcutsvalmine[5][11]=0.15;
    rdcutsvalmine[6][11]=0.15;
    rdcutsvalmine[7][11]=-0.000067;//0.0001;
    rdcutsvalmine[8][11]=0.97;//0.87;
    rdcutsvalmine[9][11]=0.3;
    rdcutsvalmine[10][11]=0.1;
    rdcutsvalmine[11][11]=0.45;
    rdcutsvalmine[12][11]=100.;
    rdcutsvalmine[13][11]=1;
     rdcutsvalmine[14][11]=0.998;
    rdcutsvalmine[15][11]=6.5;

    //16-20
    rdcutsvalmine[0][12]=0.074;
    rdcutsvalmine[1][12]=0.8;
    rdcutsvalmine[2][12]=1.0;
    rdcutsvalmine[3][12]=.3;
    rdcutsvalmine[4][12]=.3;
    rdcutsvalmine[5][12]=0.6;
    rdcutsvalmine[6][12]=0.6;
    rdcutsvalmine[7][12]=0.00015;//0.06;
    rdcutsvalmine[8][12]=0.97;//0.1;
    rdcutsvalmine[9][12]=0.3;
    rdcutsvalmine[10][12]=0.15;
    rdcutsvalmine[11][12]=0.5;
    rdcutsvalmine[12][12]=100.;
    rdcutsvalmine[13][12]=1;
     rdcutsvalmine[14][12]=0.998;
    rdcutsvalmine[15][12]=1.;

    //20-24
    rdcutsvalmine[0][13]=0.074;
    rdcutsvalmine[1][13]=0.8;
    rdcutsvalmine[2][13]=1.0;
    rdcutsvalmine[3][13]=.3;
    rdcutsvalmine[4][13]=.3;
    rdcutsvalmine[5][13]=0.6;
    rdcutsvalmine[6][13]=0.6;
    rdcutsvalmine[7][13]=0.00233;//0.01;
    rdcutsvalmine[8][13]=0.89;//0.1;
    rdcutsvalmine[9][13]=0.3;
    rdcutsvalmine[10][13]=0.15;
    rdcutsvalmine[11][13]=0.5;
    rdcutsvalmine[12][13]=100.;
    rdcutsvalmine[13][13]=1;
    rdcutsvalmine[14][13]=0.995;
    rdcutsvalmine[15][13]=0.;

    //>25
    rdcutsvalmine[0][14]=0.074;
    rdcutsvalmine[1][14]=0.6;
    rdcutsvalmine[2][14]=1.0;
    rdcutsvalmine[3][14]=.4;
    rdcutsvalmine[4][14]=.4;
    rdcutsvalmine[5][14]=0.5;
    rdcutsvalmine[6][14]=0.5;
    rdcutsvalmine[7][14]=0.033;//0.1;
    rdcutsvalmine[8][14]=0.8;//0.7;
    rdcutsvalmine[9][14]=0.1;
    rdcutsvalmine[10][14]=0.1;
    rdcutsvalmine[11][14]=0.05;
    rdcutsvalmine[12][14]=100.;
    rdcutsvalmine[13][14]=1;
    rdcutsvalmine[14][14]=0.9;
    rdcutsvalmine[15][14]=0.;
*/
/*
//  Cuts 15-02-12 (cuts for 0-20%)  -> 200212.root
    //0-0.5
    rdcutsvalmine[0][0]=0.026;
    rdcutsvalmine[1][0]=0.022;
    rdcutsvalmine[2][0]=0.7;
    rdcutsvalmine[3][0]=0.21;
    rdcutsvalmine[4][0]=0.21;
    rdcutsvalmine[5][0]=0.05;
    rdcutsvalmine[6][0]=0.05;
    rdcutsvalmine[7][0]=0.0004;
    rdcutsvalmine[8][0]=0.6;
    rdcutsvalmine[9][0]=0.3;
    rdcutsvalmine[10][0]=0.1;
    rdcutsvalmine[11][0]=0.05;
    rdcutsvalmine[12][0]=100.;
    rdcutsvalmine[13][0]=0.5;
    rdcutsvalmine[14][0]=-1.;
    rdcutsvalmine[15][0]=0.;

    //0.5-1
    rdcutsvalmine[0][1]=0.026;
    rdcutsvalmine[1][1]=0.04;
    rdcutsvalmine[2][1]=0.8;
    rdcutsvalmine[3][1]=0.5;
    rdcutsvalmine[4][1]=0.5;
    rdcutsvalmine[5][1]=0.08;
    rdcutsvalmine[6][1]=0.08;
    rdcutsvalmine[7][1]=0.0003;
    rdcutsvalmine[8][1]=0.65;
    rdcutsvalmine[9][1]=0.3;
    rdcutsvalmine[10][1]=0.1;
    rdcutsvalmine[11][1]=0.05;
    rdcutsvalmine[12][1]=100.;
    rdcutsvalmine[13][1]=0.5;
    rdcutsvalmine[14][1]=-1.;
    rdcutsvalmine[15][1]=0.;

    //1-2
    rdcutsvalmine[0][2]=0.022;
    rdcutsvalmine[1][2]=0.04;
    rdcutsvalmine[2][2]=0.8;
    rdcutsvalmine[3][2]=0.7;
    rdcutsvalmine[4][2]=0.7;
    rdcutsvalmine[5][2]=0.08;
    rdcutsvalmine[6][2]=0.08;
    rdcutsvalmine[7][2]=0.00014;
    rdcutsvalmine[8][2]=0.75;
    rdcutsvalmine[9][2]=0.3;
    rdcutsvalmine[10][2]=0.1;
    rdcutsvalmine[11][2]=0.05;
    rdcutsvalmine[12][2]=100.;
    rdcutsvalmine[13][2]=0.5;
    rdcutsvalmine[14][2]=-1.;
    rdcutsvalmine[15][2]=0.;
    //2-3
    rdcutsvalmine[0][3]=0.024;
    rdcutsvalmine[1][3]=0.025;
    rdcutsvalmine[2][3]=0.8;
    rdcutsvalmine[3][3]=0.7;
    rdcutsvalmine[4][3]=0.7;
    rdcutsvalmine[5][3]=0.08;
    rdcutsvalmine[6][3]=0.08;
    rdcutsvalmine[7][3]=0.00014;
    rdcutsvalmine[8][3]=0.85;
    rdcutsvalmine[9][3]=0.3;
    rdcutsvalmine[10][3]=0.1;
    rdcutsvalmine[11][3]=0.05;
    rdcutsvalmine[12][3]=100.;
    rdcutsvalmine[13][3]=0.5;
    rdcutsvalmine[14][3]=-1.;
    rdcutsvalmine[15][3]=0.;


    //3-4
    rdcutsvalmine[0][4]=0.032;
    rdcutsvalmine[1][4]=0.025; // 0.8
    rdcutsvalmine[2][4]=0.7;
    rdcutsvalmine[3][4]=1.1;
    rdcutsvalmine[4][4]=1.1;
    rdcutsvalmine[5][4]=0.08;
    rdcutsvalmine[6][4]=0.08;
    rdcutsvalmine[7][4]=-0.00027;
    rdcutsvalmine[8][4]=0.995;
    rdcutsvalmine[9][4]=0.3;
    rdcutsvalmine[10][4]=0.15;
    rdcutsvalmine[11][4]=0.05;
    rdcutsvalmine[12][4]=0.5;
    rdcutsvalmine[13][4]=0.5;
    rdcutsvalmine[14][4]=0.998;
    rdcutsvalmine[15][4]=9.5;

    //4-5
    rdcutsvalmine[0][5]=0.032;
    rdcutsvalmine[1][5]=0.025; //0.1
    rdcutsvalmine[2][5]=0.7;
    rdcutsvalmine[3][5]=1.2;
    rdcutsvalmine[4][5]=1.2;
    rdcutsvalmine[5][5]=0.09;
    rdcutsvalmine[6][5]=0.09;
    rdcutsvalmine[7][5]=-0.00027;
    rdcutsvalmine[8][5]=0.998;
    rdcutsvalmine[9][5]=0.15;
    rdcutsvalmine[10][5]=0.15;
    rdcutsvalmine[11][5]=0.05;
    rdcutsvalmine[12][5]=10;
    rdcutsvalmine[13][5]=0.5;
    rdcutsvalmine[14][5]=0.998;
    rdcutsvalmine[15][5]=6.;

    //5-6
    rdcutsvalmine[0][6]=0.04;
    rdcutsvalmine[1][6]=0.027;
    rdcutsvalmine[2][6]=0.8;
    rdcutsvalmine[3][6]=1.1;
    rdcutsvalmine[4][6]=1.1;
    rdcutsvalmine[5][6]=0.1;
    rdcutsvalmine[6][6]=0.1;
    rdcutsvalmine[7][6]=-0.0002;
    rdcutsvalmine[8][6]=0.996;
    rdcutsvalmine[9][6]=0.15;
    rdcutsvalmine[10][6]=0.15;
    rdcutsvalmine[11][6]=0.2;
    rdcutsvalmine[12][6]=10.;
    rdcutsvalmine[13][6]=0.5;
    rdcutsvalmine[14][6]=0.998;
    rdcutsvalmine[15][6]=5.2;

    //6-7
    rdcutsvalmine[0][7]=0.036;
    rdcutsvalmine[1][7]=0.027;
    rdcutsvalmine[2][7]=0.8;
    rdcutsvalmine[3][7]=0.9;
    rdcutsvalmine[4][7]=0.9;
    rdcutsvalmine[5][7]=0.1;
    rdcutsvalmine[6][7]=0.1;
    rdcutsvalmine[7][7]=-0.00019;
    rdcutsvalmine[8][7]=0.99;
    rdcutsvalmine[9][7]=0.3;
    rdcutsvalmine[10][7]=0.1;
    rdcutsvalmine[11][7]=0.2;
    rdcutsvalmine[12][7]=100.;
    rdcutsvalmine[13][7]=0.5;
    rdcutsvalmine[14][7]=0.998;
    rdcutsvalmine[15][7]=5.;

    //7-8
    rdcutsvalmine[0][8]=0.036;
    rdcutsvalmine[1][8]=0.027;
    rdcutsvalmine[2][8]=1;
    rdcutsvalmine[3][8]=0.9;
    rdcutsvalmine[4][8]=0.9;
    rdcutsvalmine[5][8]=0.1;
    rdcutsvalmine[6][8]=0.1;
    rdcutsvalmine[7][8]=-0.00016;
    rdcutsvalmine[8][8]=0.99;
    rdcutsvalmine[9][8]=0.3;
    rdcutsvalmine[10][8]=0.15;
    rdcutsvalmine[11][8]=0.2;
    rdcutsvalmine[12][8]=100.;
    rdcutsvalmine[13][8]=0.5;
    rdcutsvalmine[14][8]=0.998;
    rdcutsvalmine[15][8]=5.;

    //8-10
    rdcutsvalmine[0][9]=0.055;
    rdcutsvalmine[1][9]=0.03;
    rdcutsvalmine[2][9]=1.0;
    rdcutsvalmine[3][9]=0.9;
    rdcutsvalmine[4][9]=0.9;
    rdcutsvalmine[5][9]=0.15;
    rdcutsvalmine[6][9]=0.15;
    rdcutsvalmine[7][9]=-0.00008;//-0.00001;
    rdcutsvalmine[8][9]=0.99;//0.92;
    rdcutsvalmine[9][9]=0.3;
    rdcutsvalmine[10][9]=0.15;
    rdcutsvalmine[11][9]=0.34;
    rdcutsvalmine[12][9]=100.;
    rdcutsvalmine[13][9]=0.5;
     rdcutsvalmine[14][9]=0.998;
    rdcutsvalmine[15][9]=4.8;

    //10-12
    rdcutsvalmine[0][10]=0.055;
    rdcutsvalmine[1][10]=0.03;
    rdcutsvalmine[2][10]=1.0;
    rdcutsvalmine[3][10]=1.0;
    rdcutsvalmine[4][10]=1.0;
    rdcutsvalmine[5][10]=0.15;
    rdcutsvalmine[6][10]=0.15;
    rdcutsvalmine[7][10]=-0.00004;//-0.00001;
    rdcutsvalmine[8][10]=0.99;//0.92;
    rdcutsvalmine[9][10]=0.3;
    rdcutsvalmine[10][10]=0.15;
    rdcutsvalmine[11][10]=0.35;
    rdcutsvalmine[12][10]=100.;
    rdcutsvalmine[13][10]=0.5;
     rdcutsvalmine[14][10]=0.998;
    rdcutsvalmine[15][10]=4.;

    //12-16
    rdcutsvalmine[0][11]=0.074;
    rdcutsvalmine[1][11]=0.07;
    rdcutsvalmine[2][11]=1.0;
    rdcutsvalmine[3][11]=0.3;
    rdcutsvalmine[4][11]=0.3;
    rdcutsvalmine[5][11]=0.15;
    rdcutsvalmine[6][11]=0.15;
    rdcutsvalmine[7][11]=-0.0000001;//0.0001;
    rdcutsvalmine[8][11]=0.94;//0.87;
    rdcutsvalmine[9][11]=0.3;
    rdcutsvalmine[10][11]=0.1;
    rdcutsvalmine[11][11]=0.45;
    rdcutsvalmine[12][11]=100.;
    rdcutsvalmine[13][11]=1;
     rdcutsvalmine[14][11]=0.996;
    rdcutsvalmine[15][11]=4;

    //16-20
    rdcutsvalmine[0][12]=0.074;
    rdcutsvalmine[1][12]=0.8;
    rdcutsvalmine[2][12]=1.0;
    rdcutsvalmine[3][12]=.3;
    rdcutsvalmine[4][12]=.3;
    rdcutsvalmine[5][12]=0.6;
    rdcutsvalmine[6][12]=0.6;
    rdcutsvalmine[7][12]=0.00015;//0.06;
    rdcutsvalmine[8][12]=0.93;//0.1;
    rdcutsvalmine[9][12]=0.3;
    rdcutsvalmine[10][12]=0.15;
    rdcutsvalmine[11][12]=0.5;
    rdcutsvalmine[12][12]=100.;
    rdcutsvalmine[13][12]=1;
     rdcutsvalmine[14][12]=0.996;
    rdcutsvalmine[15][12]=0.;

    //20-24
    rdcutsvalmine[0][13]=0.074;
    rdcutsvalmine[1][13]=0.8;
    rdcutsvalmine[2][13]=1.0;
    rdcutsvalmine[3][13]=.3;
    rdcutsvalmine[4][13]=.3;
    rdcutsvalmine[5][13]=0.6;
    rdcutsvalmine[6][13]=0.6;
    rdcutsvalmine[7][13]=0.00233;//0.01;
    rdcutsvalmine[8][13]=0.7;//0.1;
    rdcutsvalmine[9][13]=0.3;
    rdcutsvalmine[10][13]=0.15;
    rdcutsvalmine[11][13]=0.5;
    rdcutsvalmine[12][13]=100.;
    rdcutsvalmine[13][13]=1;
    rdcutsvalmine[14][13]=0.99;
    rdcutsvalmine[15][13]=0.;

    //>25
    rdcutsvalmine[0][14]=0.074;
    rdcutsvalmine[1][14]=0.6;
    rdcutsvalmine[2][14]=1.0;
    rdcutsvalmine[3][14]=.4;
    rdcutsvalmine[4][14]=.4;
    rdcutsvalmine[5][14]=0.5;
    rdcutsvalmine[6][14]=0.5;
    rdcutsvalmine[7][14]=0.033;//0.1;
    rdcutsvalmine[8][14]=0.7;//0.7;
    rdcutsvalmine[9][14]=0.1;
    rdcutsvalmine[10][14]=0.1;
    rdcutsvalmine[11][14]=0.05;
    rdcutsvalmine[12][14]=100.;
    rdcutsvalmine[13][14]=1;
    rdcutsvalmine[14][14]=0.9;
    rdcutsvalmine[15][14]=0.;
*/

  }
  if(set_cuts=="heidelberg"){

    //0-0.5
    rdcutsvalmine[0][0]=0.7;
    rdcutsvalmine[1][0]=0.03;
    rdcutsvalmine[2][0]=0.7;
    rdcutsvalmine[3][0]=0.8;
    rdcutsvalmine[4][0]=0.8;
    rdcutsvalmine[5][0]=0.1;
    rdcutsvalmine[6][0]=0.1;
    rdcutsvalmine[7][0]=-0.00002;
    rdcutsvalmine[8][0]=0.9;
    rdcutsvalmine[9][0]=0.3;
    rdcutsvalmine[10][0]=0.1;
    rdcutsvalmine[11][0]=0.05;
    rdcutsvalmine[12][0]=100.;
    rdcutsvalmine[13][0]=0.5;
    //0.5-1
    rdcutsvalmine[0][1]=0.7;
    rdcutsvalmine[1][1]=0.03;
    rdcutsvalmine[2][1]=0.7;
    rdcutsvalmine[3][1]=0.8;
    rdcutsvalmine[4][1]=0.8;
    rdcutsvalmine[5][1]=0.1;
    rdcutsvalmine[6][1]=0.1;
    rdcutsvalmine[7][1]=-0.00002;
    rdcutsvalmine[8][1]=0.9;
    rdcutsvalmine[9][1]=0.3;
    rdcutsvalmine[10][1]=0.1;
    rdcutsvalmine[11][1]=0.05;
    rdcutsvalmine[12][1]=100.;
    rdcutsvalmine[13][1]=0.5;

    //1-2
    rdcutsvalmine[0][2]=0.7;
    rdcutsvalmine[1][2]=0.03;
    rdcutsvalmine[2][2]=0.7;
    rdcutsvalmine[3][2]=0.8;
    rdcutsvalmine[4][2]=0.8;
    rdcutsvalmine[5][2]=0.1;
    rdcutsvalmine[6][2]=0.1;
    rdcutsvalmine[7][2]=-0.00002;
    rdcutsvalmine[8][2]=0.9;
    rdcutsvalmine[9][2]=0.3;
    rdcutsvalmine[10][2]=0.1;
    rdcutsvalmine[11][2]=0.05;
    rdcutsvalmine[12][2]=100.;
    rdcutsvalmine[13][2]=0.5;
    //2-3
    rdcutsvalmine[0][3]=0.7;
    rdcutsvalmine[1][3]=0.03;
    rdcutsvalmine[2][3]=0.7;
    rdcutsvalmine[3][3]=0.8;
    rdcutsvalmine[4][3]=0.8;
    rdcutsvalmine[5][3]=0.1;
    rdcutsvalmine[6][3]=0.1;
    rdcutsvalmine[7][3]=-0.00002;
    rdcutsvalmine[8][3]=0.9;
    rdcutsvalmine[9][3]=0.3;
    rdcutsvalmine[10][3]=0.1;
    rdcutsvalmine[11][3]=0.05;
    rdcutsvalmine[12][3]=100.;
    rdcutsvalmine[13][3]=0.5;
    //3-4
    rdcutsvalmine[0][4]=0.7;
    rdcutsvalmine[1][4]=0.03;
    rdcutsvalmine[2][4]=0.7;
    rdcutsvalmine[3][4]=0.9;
    rdcutsvalmine[4][4]=0.9;
    rdcutsvalmine[5][4]=0.1;
    rdcutsvalmine[6][4]=0.1;
    rdcutsvalmine[7][4]=0.000002;
    rdcutsvalmine[8][4]=0.8;
    rdcutsvalmine[9][4]=0.3;
    rdcutsvalmine[10][4]=0.1;
    rdcutsvalmine[11][4]=0.05;
    rdcutsvalmine[12][4]=100.;
    rdcutsvalmine[13][4]=0.5;
    //4-5
    rdcutsvalmine[0][5]=0.7;
    rdcutsvalmine[1][5]=0.03;
    rdcutsvalmine[2][5]=0.7;
    rdcutsvalmine[3][5]=0.9;
    rdcutsvalmine[4][5]=0.9;
    rdcutsvalmine[5][5]=0.1;
    rdcutsvalmine[6][5]=0.1;
    rdcutsvalmine[7][5]=0.000002;
    rdcutsvalmine[8][5]=0.8;
    rdcutsvalmine[9][5]=0.3;
    rdcutsvalmine[10][5]=0.1;
    rdcutsvalmine[11][5]=0.05;
    rdcutsvalmine[12][5]=100.;
    rdcutsvalmine[13][5]=0.5;
    //5-6
    rdcutsvalmine[0][6]=0.7;
    rdcutsvalmine[1][6]=0.03;
    rdcutsvalmine[2][6]=0.7;
    rdcutsvalmine[3][6]=1.0;
    rdcutsvalmine[4][6]=1.0;
    rdcutsvalmine[5][6]=0.1;
    rdcutsvalmine[6][6]=0.1;
    rdcutsvalmine[7][6]=0.000002;
    rdcutsvalmine[8][6]=0.8;
    rdcutsvalmine[9][6]=0.3;
    rdcutsvalmine[10][6]=0.1;
    rdcutsvalmine[11][6]=0.05;
    rdcutsvalmine[12][6]=100.;
    rdcutsvalmine[13][6]=0.5;
    //6-8
    rdcutsvalmine[0][7]=0.7;
    rdcutsvalmine[1][7]=0.03;
    rdcutsvalmine[2][7]=0.7;
    rdcutsvalmine[3][7]=1.0;
    rdcutsvalmine[4][7]=1.0;
    rdcutsvalmine[5][7]=0.1;
    rdcutsvalmine[6][7]=0.1;
    rdcutsvalmine[7][7]=0.000002;
    rdcutsvalmine[8][7]=0.8;
    rdcutsvalmine[9][7]=0.3;
    rdcutsvalmine[10][7]=0.1;
    rdcutsvalmine[11][7]=0.05;
    rdcutsvalmine[12][7]=100.;
    rdcutsvalmine[13][7]=0.5;
    //8-12
    rdcutsvalmine[0][8]=0.7;
    rdcutsvalmine[1][8]=0.03;
    rdcutsvalmine[2][8]=0.7;
    rdcutsvalmine[3][8]=1.0;
    rdcutsvalmine[4][8]=1.0;
    rdcutsvalmine[5][8]=0.1;
    rdcutsvalmine[6][8]=0.1;
    rdcutsvalmine[7][8]=0.000002;
    rdcutsvalmine[8][8]=0.8;
    rdcutsvalmine[9][8]=0.3;
    rdcutsvalmine[10][8]=0.1;
    rdcutsvalmine[11][8]=0.05;
    rdcutsvalmine[12][8]=100.;
    rdcutsvalmine[13][8]=0.5;
    //12-16
    rdcutsvalmine[0][9]=0.7;
    rdcutsvalmine[1][9]=0.03;
    rdcutsvalmine[2][9]=0.7;
    rdcutsvalmine[3][9]=1.0;
    rdcutsvalmine[4][9]=1.0;
    rdcutsvalmine[5][9]=0.1;
    rdcutsvalmine[6][9]=0.1;
    rdcutsvalmine[7][9]=0.000002;
    rdcutsvalmine[8][9]=0.8;
    rdcutsvalmine[9][9]=0.3;
    rdcutsvalmine[10][9]=0.1;
    rdcutsvalmine[11][9]=0.05;
    rdcutsvalmine[12][9]=100.;
    rdcutsvalmine[13][9]=0.5;
    //16-20
    rdcutsvalmine[0][10]=0.7;
    rdcutsvalmine[1][10]=0.03;
    rdcutsvalmine[2][10]=0.7;
    rdcutsvalmine[3][10]=1.0;
    rdcutsvalmine[4][10]=1.0;
    rdcutsvalmine[5][10]=0.1;
    rdcutsvalmine[6][10]=0.1;
    rdcutsvalmine[7][10]=0.000002;
    rdcutsvalmine[8][10]=0.8;
    rdcutsvalmine[9][10]=0.3;
    rdcutsvalmine[10][10]=0.1;
    rdcutsvalmine[11][10]=0.05;
    rdcutsvalmine[12][10]=100.;
    rdcutsvalmine[13][10]=0.5;
    //20-24
    rdcutsvalmine[0][11]=0.7;
    rdcutsvalmine[1][11]=0.03;
    rdcutsvalmine[2][11]=0.7;
    rdcutsvalmine[3][11]=1.0;
    rdcutsvalmine[4][11]=1.0;
    rdcutsvalmine[5][11]=0.1;
    rdcutsvalmine[6][11]=0.1;
    rdcutsvalmine[7][11]=0.000002;
    rdcutsvalmine[8][11]=0.8;
    rdcutsvalmine[9][11]=0.3;
    rdcutsvalmine[10][11]=0.1;
    rdcutsvalmine[11][11]=0.05;
    rdcutsvalmine[12][11]=100.;
    rdcutsvalmine[13][11]=0.5;
    //>24
    rdcutsvalmine[0][12]=0.7;
    rdcutsvalmine[1][12]=0.03;
    rdcutsvalmine[2][12]=0.7;
    rdcutsvalmine[3][12]=1.0;
    rdcutsvalmine[4][12]=1.0;
    rdcutsvalmine[5][12]=0.1;
    rdcutsvalmine[6][12]=0.1;
    rdcutsvalmine[7][12]=0.000002;
    rdcutsvalmine[8][12]=0.8;
    rdcutsvalmine[9][12]=0.3;
    rdcutsvalmine[10][12]=0.1;
    rdcutsvalmine[11][12]=0.05;
    rdcutsvalmine[12][12]=100.;
    rdcutsvalmine[13][12]=0.5;
  }

  RDHFDStartoKpipi->SetCuts(nvars,nptbins,rdcutsvalmine);

   Bool_t pidflag=kTRUE;
  RDHFDStartoKpipi->SetUsePID(pidflag);
  if(pidflag) cout<<"PID is used"<<endl;
  else cout<<"PID is not used"<<endl;

  AliAODPidHF* pidObj=new AliAODPidHF();

  //kaon
  Double_t sigmasK[5]={3.,3.,3.,3.,3.};
  pidObj->SetSigma(sigmasK);
  pidObj->SetAsym(kTRUE);
  pidObj->SetMatch(1);
  pidObj->SetTPC(kTRUE);
  pidObj->SetTOF(kTRUE);
  Double_t plimK[2]={0.5,0.8};
  pidObj->SetPLimit(plimK,2);
  //pidObjK->SetTOFdecide(kTRUE);

  RDHFDStartoKpipi->SetPidHF(pidObj);

  // PID SETTINGS
  // AliAODPidHF* pidObj=new AliAODPidHF();
  // pidObj->SetName("pid4DSatr");
  // // Int_t mode=1;
  // Double_t priors[5]={0.01,0.001,0.3,0.3,0.3};
  // pidObj->SetPriors(priors);
  // pidObj->SetMatch(mode);
  // pidObj->SetSigma(0,3); // TPC ----> previously it was 3sig, but this reduces background
  // pidObj->SetSigma(3,3); // TOF
  // pidObj->SetTPC(kTRUE);
  // pidObj->SetTOF(kTRUE);
  // pidObj->SetMaxTrackMomForCombinedPID(3.);
  // RDHFDStartoKpipi->SetPidHF(pidObj);

  //activate pileup rejection
  RDHFDStartoKpipi->SetOptPileup(AliRDHFCuts::kRejectPileupEvent);

  //TFile *fFlat=TFile::Open("./CentrDistrBins005.root","READ");
  // TCanvas *c=fFlat->Get("cintegral");
  // TH1F *hfl=(TH1F*)c->FindObject("hint");
  // RDHFDStartoKpipi->SetHistoForCentralityFlattening(hfl,40.,50.,0.,0);
  /*
  RDHFDStartoKpipi->PrintAll();

  TFile* fout=new TFile("DStartoKpipiCuts.root","recreate");   //set this!!
  fout->cd();
  RDHFDStartoKpipi->Write();
  fout->Close();
  */
  // TFile *fFlat=TFile::Open("./CentrDistrBins005.root","READ");
  // TCanvas *c=fFlat->Get("cintegral");
  // TH1F *hfl=(TH1F*)c->FindObject("hint");
  // RDHFDStartoKpipi->SetHistoForCentralityFlattening(hfl,0.,10.,0.,0);

  RDHFDStartoKpipi->PrintAll();

  TFile* fout=new TFile("DStartoKpipiCuts_pp13TeV_kHighMultV0.root","recreate");   //set this!!
  fout->cd();
  RDHFDStartoKpipi->Write();
  fout->Close();

}

//macro to make a .root file (for significance maximization) which contains an AliRDHFCutsDStartoKpipi with loose set of cuts  and TParameter with the tighest value of these cuts
// copy form D0 ... NOT TESTED YIET ... to be tested!!

void makeInputAliAnalysisTaskSEDstarSignificanceMaximization(){

  AliRDHFCutsDStartoKpipi* RDHFDStartoKpipi=new AliRDHFCutsDStartoKpipi();

//Trigger when running over data
  RDHFDStartoKpipi->SetTriggerClass("");
  RDHFDStartoKpipi->ResetMaskAndEnableMBTrigger();
  RDHFDStartoKpipi->EnableCentralTrigger();
  RDHFDStartoKpipi->EnableSemiCentralTrigger();

  //Centrality selection
  RDHFDStartoKpipi->SetUseCentrality(kTRUE);
  RDHFDStartoKpipi->SetMinCentrality(40.);
  RDHFDStartoKpipi->SetMaxCentrality(50.);
  //RDHFDStartoKpipi->SetMinCentrality(0.);
 // RDHFDStartoKpipi->SetMaxCentrality(20.);

  RDHFDStartoKpipi->SetName("DStartoKpipiCuts");
  RDHFDStartoKpipi->SetTitle("Cuts for D* analysis");

  AliESDtrackCuts* esdTrackCuts=new AliESDtrackCuts();
  esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
  //default
  esdTrackCuts->SetRequireTPCRefit(kTRUE);
  esdTrackCuts->SetRequireITSRefit(kTRUE);
  esdTrackCuts->SetMinNClustersTPC(70);
  esdTrackCuts->SetMinNClustersITS(2);
  esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny); //kFirst; kAny I selected myself
  esdTrackCuts->SetMinDCAToVertexXY(0.);
  esdTrackCuts->SetPtRange(0.3,1.e10);

 // default is kBoth, otherwise kAny
  esdTrackCuts->SetMinDCAToVertexXY(0.);
 // esdTrackCuts->SetPtRange(0.3,1.e10);

 // soft pion pre-selections
  AliESDtrackCuts* esdSoftPicuts=new AliESDtrackCuts();
  esdSoftPicuts->SetRequireSigmaToVertex(kFALSE);
  //default
  esdSoftPicuts->SetRequireTPCRefit(kTRUE);//////////////////////////////
  esdSoftPicuts->SetRequireITSRefit(kTRUE);
  //esdSoftPicuts->SetMinNClustersTPC(40);
  //esdSoftPicuts->SetMinNClustersITS(3);
  esdSoftPicuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny); //test d0 asimmetry ////////////////////////////
  esdSoftPicuts->SetPtRange(0.0,1.e10);

  // set pre selections
  RDHFDStartoKpipi->AddTrackCuts(esdTrackCuts);
  RDHFDStartoKpipi->AddTrackCutsSoftPi(esdSoftPicuts);

  const Int_t nvars=16;
  const Int_t nptbins=12;

  Float_t* ptbins;
  ptbins=new Float_t[nptbins+1];
  ptbins[0]=1.;
  ptbins[1]=2.;
  ptbins[2]=3.;
  ptbins[3]=4.;
  ptbins[4]=5.;
  ptbins[5]=6.;
  ptbins[6]=7.;
  ptbins[7]=8.;
  ptbins[8]=10.;
  ptbins[9]=12.;
  ptbins[10]=16.;
  ptbins[11]=20.;
  ptbins[12]=24.;
  RDHFDStartoKpipi->SetPtBins(nptbins+1,ptbins);

  Float_t** rdcutsvalmine;
  rdcutsvalmine=new Float_t*[nvars];
  for(Int_t iv=0;iv<nvars;iv++){
    rdcutsvalmine[iv]=new Float_t[nptbins];
  }
      //
    //1-2
    rdcutsvalmine[0][0]=0.032;
    rdcutsvalmine[1][0]=0.0225;
    rdcutsvalmine[2][0]=0.8;
    rdcutsvalmine[3][0]=1.0;
    rdcutsvalmine[4][0]=1.0;
    rdcutsvalmine[5][0]=0.1;
    rdcutsvalmine[6][0]=0.1;
    rdcutsvalmine[7][0]=-0.0002;
    rdcutsvalmine[8][0]=0.97;
    rdcutsvalmine[9][0]=0.3;
    rdcutsvalmine[10][0]=0.15;
    rdcutsvalmine[11][0]=0.05;
    rdcutsvalmine[12][0]=0.5;
    rdcutsvalmine[13][0]=0.5;
    rdcutsvalmine[14][0]=0.998;
    rdcutsvalmine[15][0]=6.0;

//2-3
    rdcutsvalmine[0][1]=0.032;
    rdcutsvalmine[1][1]=0.045;
    rdcutsvalmine[2][1]=0.8;
    rdcutsvalmine[3][1]=1.0;
    rdcutsvalmine[4][1]=1.0;
    rdcutsvalmine[5][1]=0.1;
    rdcutsvalmine[6][1]=0.1;
    rdcutsvalmine[7][1]=-0.00015;
    rdcutsvalmine[8][1]=0.96;
    rdcutsvalmine[9][1]=0.3;
    rdcutsvalmine[10][1]=0.15;
    rdcutsvalmine[11][1]=0.05;
    rdcutsvalmine[12][1]=0.5;
    rdcutsvalmine[13][1]=0.5;
    rdcutsvalmine[14][1]=0.998;
    rdcutsvalmine[15][1]=3.5;

    //3-4
    rdcutsvalmine[0][2]=0.032;
    rdcutsvalmine[1][2]=0.045;
    rdcutsvalmine[2][2]=0.8;
    rdcutsvalmine[3][2]=1.0;
    rdcutsvalmine[4][2]=1.0;
    rdcutsvalmine[5][2]=0.1;
    rdcutsvalmine[6][2]=0.1;
    rdcutsvalmine[7][2]=-0.0001;
    rdcutsvalmine[8][2]=0.95;
    rdcutsvalmine[9][2]=0.3;
    rdcutsvalmine[10][2]=0.15;
    rdcutsvalmine[11][2]=0.05;
    rdcutsvalmine[12][2]=0.5;
    rdcutsvalmine[13][2]=1;
    rdcutsvalmine[14][2]=0.998;
    rdcutsvalmine[15][2]=3.;

    //4-5
    rdcutsvalmine[0][3]=0.032;
    rdcutsvalmine[1][3]=0.035;
    rdcutsvalmine[2][3]=0.8;
    rdcutsvalmine[3][3]=1.;
    rdcutsvalmine[4][3]=1.;
    rdcutsvalmine[5][3]=0.1;
    rdcutsvalmine[6][3]=0.1;
    rdcutsvalmine[7][3]=-0.00011;
    rdcutsvalmine[8][3]=0.94;
    rdcutsvalmine[9][3]=0.20;
    rdcutsvalmine[10][3]=0.15;
    rdcutsvalmine[11][3]=0.05;
    rdcutsvalmine[12][3]=10;
    rdcutsvalmine[13][3]=1;
    rdcutsvalmine[14][3]=0.998;
    rdcutsvalmine[15][3]=3.;

//5-6
    rdcutsvalmine[0][4]=0.036;
    rdcutsvalmine[1][4]=0.035;
    rdcutsvalmine[2][4]=1.;
    rdcutsvalmine[3][4]=0.8;
    rdcutsvalmine[4][4]=0.8;
    rdcutsvalmine[5][4]=0.1;
    rdcutsvalmine[6][4]=0.1;
    rdcutsvalmine[7][4]=-0.0001;
    rdcutsvalmine[8][4]=0.95;
    rdcutsvalmine[9][4]=0.3;
    rdcutsvalmine[10][4]=0.15;
    rdcutsvalmine[11][4]=0.2;
    rdcutsvalmine[12][4]=10.;
    rdcutsvalmine[13][4]=1;
    rdcutsvalmine[14][4]=0.998;
    rdcutsvalmine[15][4]=3.;

//6-7
    rdcutsvalmine[0][5]=0.036;
    rdcutsvalmine[1][5]=0.035;
    rdcutsvalmine[2][5]=1.0;
    rdcutsvalmine[3][5]=0.8;
    rdcutsvalmine[4][5]=0.8;
    rdcutsvalmine[5][5]=0.1;
    rdcutsvalmine[6][5]=0.1;
    rdcutsvalmine[7][5]=-0.00004;
    rdcutsvalmine[8][5]=0.88;//0.97;
    rdcutsvalmine[9][5]=0.3;
    rdcutsvalmine[10][5]=0.15;
    rdcutsvalmine[11][5]=0.2;
    rdcutsvalmine[12][5]=100.;
    rdcutsvalmine[13][5]=1;
    rdcutsvalmine[14][5]=0.996;
    rdcutsvalmine[15][5]=1.;

//7-8
    rdcutsvalmine[0][6]=0.04;
    rdcutsvalmine[1][6]=0.035;
    rdcutsvalmine[2][6]=1.0;
    rdcutsvalmine[3][6]=0.7;
    rdcutsvalmine[4][6]=0.7;
    rdcutsvalmine[5][6]=0.1;
    rdcutsvalmine[6][6]=0.1;
    rdcutsvalmine[7][6]=-0.00001;
    rdcutsvalmine[8][6]=0.88;//0.97;
    rdcutsvalmine[9][6]=0.3;
    rdcutsvalmine[10][6]=0.15;
    rdcutsvalmine[11][6]=0.2;
    rdcutsvalmine[12][6]=100.;
    rdcutsvalmine[13][6]=1.0;
    rdcutsvalmine[14][6]=0.996;
    rdcutsvalmine[15][6]=0;

    //8-10
    rdcutsvalmine[0][7]=0.055;
    rdcutsvalmine[1][7]=0.04;
    rdcutsvalmine[2][7]=1.0;
    rdcutsvalmine[3][7]=0.7;
    rdcutsvalmine[4][7]=0.7;
    rdcutsvalmine[5][7]=0.1;
    rdcutsvalmine[6][7]=0.1;
    rdcutsvalmine[7][7]=-0.0000033;
    rdcutsvalmine[8][7]=0.88;//0.96;
    rdcutsvalmine[9][7]=0.3;
    rdcutsvalmine[10][7]=0.15;
    rdcutsvalmine[11][7]=0.3;
    rdcutsvalmine[12][7]=100.;
    rdcutsvalmine[13][7]=1.0;
    rdcutsvalmine[14][7]=0.995;
    rdcutsvalmine[15][7]=0;//4.;

//10-12     230312:
    rdcutsvalmine[0][8]=0.060;
    rdcutsvalmine[1][8]=0.042;
    rdcutsvalmine[2][8]=1.0;
    rdcutsvalmine[3][8]=0.65;
    rdcutsvalmine[4][8]=0.65;
    rdcutsvalmine[5][8]=0.1;
    rdcutsvalmine[6][8]=0.1;
    rdcutsvalmine[7][8]=-0.000001;
    rdcutsvalmine[8][8]=0.8;//0.96;
    rdcutsvalmine[9][8]=0.3;
    rdcutsvalmine[10][8]=0.15;
    rdcutsvalmine[11][8]=0.35;
    rdcutsvalmine[12][8]=100.;
    rdcutsvalmine[13][8]=1.0;
    rdcutsvalmine[14][8]=0.995;
    rdcutsvalmine[15][8]=0.;//4.;

//12-16 (ptbin11)
    rdcutsvalmine[0][9]=0.074;
    rdcutsvalmine[1][9]=0.08;//0.035;
    rdcutsvalmine[2][9]=0.8;
    rdcutsvalmine[3][9]=0.3;
    rdcutsvalmine[4][9]=0.3;
    rdcutsvalmine[5][9]=0.2;
    rdcutsvalmine[6][9]=0.2;
    rdcutsvalmine[7][9]=0.0005;
    rdcutsvalmine[8][9]=0.8;
    rdcutsvalmine[9][9]=0.3;
    rdcutsvalmine[10][9]=0.15;
    rdcutsvalmine[11][9]=0.45;
    rdcutsvalmine[12][9]=100.;
    rdcutsvalmine[13][9]=1;
    rdcutsvalmine[14][9]=0.99;
    rdcutsvalmine[15][9]=0.0;//3;

    //16-20 (ptbin12)
    rdcutsvalmine[0][10]=0.074;
    rdcutsvalmine[1][10]=0.03;
    rdcutsvalmine[2][10]=1.0;
    rdcutsvalmine[3][10]=0.30;
    rdcutsvalmine[4][10]=0.30;
    rdcutsvalmine[5][10]=0.50;
    rdcutsvalmine[6][10]=0.50;
    rdcutsvalmine[7][10]=0.004;
    rdcutsvalmine[8][10]=0.85;
    rdcutsvalmine[9][10]=0.3;
    rdcutsvalmine[10][10]=0.15;
    rdcutsvalmine[11][10]=0.05;
    rdcutsvalmine[12][10]=100.;
    rdcutsvalmine[13][10]=1;
    rdcutsvalmine[14][10]=0.9;
    rdcutsvalmine[15][10]=0.0;//2.0;

//20-24 (ptbin13)
    rdcutsvalmine[0][11]=0.074;
    rdcutsvalmine[1][11]=0.09;
    rdcutsvalmine[2][11]=1.0;
    rdcutsvalmine[3][11]=0.30;
    rdcutsvalmine[4][11]=0.30;
    rdcutsvalmine[5][11]=0.50;
    rdcutsvalmine[6][11]=0.50;
    rdcutsvalmine[7][11]=0.003;
    rdcutsvalmine[8][11]=0.85;//0.93;
    rdcutsvalmine[9][11]=0.3;
    rdcutsvalmine[10][11]=0.15;
    rdcutsvalmine[11][11]=0.05;
    rdcutsvalmine[12][11]=100.;
    rdcutsvalmine[13][11]=1;
    rdcutsvalmine[14][11]=0.88;
    rdcutsvalmine[15][11]=0;

//>24 (ptbin14)
    rdcutsvalmine[0][12]=0.074;
    rdcutsvalmine[1][12]=0.1;
    rdcutsvalmine[2][12]=1.0;
    rdcutsvalmine[3][12]=0.30;
    rdcutsvalmine[4][12]=0.30;
    rdcutsvalmine[5][12]=0.50;
    rdcutsvalmine[6][12]=0.50;
    rdcutsvalmine[7][12]=0.0035;
    rdcutsvalmine[8][12]=0.7;//0.95;
    rdcutsvalmine[9][12]=0.3;
    rdcutsvalmine[10][12]=0.15;
    rdcutsvalmine[11][12]=0.05;
    rdcutsvalmine[12][12]=100.;
    rdcutsvalmine[13][12]=1;
    rdcutsvalmine[14][12]=0.9;
    rdcutsvalmine[15][12]=0;


  RDHFDStartoKpipi->SetCuts(nvars,nptbins,rdcutsvalmine);


  Int_t nvarsforopt=RDHFDStartoKpipi->GetNVarsForOpt();
  const Int_t dim=4; //set this!!
  const Bool_t *boolforopt;
  boolforopt=new Bool_t[nvars];
  if(dim>nvarsforopt){
    cout<<"Number of variables for optimization has probably changed, check and edit accordingly"<<endl;
    return;
  } else {
    if(dim==nvarsforopt){
      boolforopt=RDHFDStartoKpipi->GetVarsForOpt();
    }else{
      TString *names;
      names=new TString[nvars];
      TString answer="";
      Int_t checktrue=0;
      names=RDHFDStartoKpipi->GetVarNames();
      for(Int_t i=0;i<nvars;i++){
	cout<<names[i]<<" for opt? (y/n)"<<endl;
	cin>>answer;
	if(answer=="y") {
	  boolforopt[i]=kTRUE;
	  checktrue++;
	}
	else boolforopt[i]=kFALSE;
      }
      if (checktrue!=dim) {
	cout<<"Error! You set "<<checktrue<<" kTRUE instead of "<<dim<<endl;
	return;
      }
      RDHFDStartoKpipi->SetVarsForOpt(dim,boolforopt);
    }
  }


  Float_t tighterval[dim][nptbins];

  //ptK
  //ptPi
  //MD0
  //theta

// [1, 2]
  tighterval[0][0]=0.02;
  tighterval[1][0]=-0.00035;
  tighterval[2][0]=0.995;
  tighterval[3][0]=7;

// [2,3]
  tighterval[0][1]=0.02;
  tighterval[1][1]=-0.00022;
  tighterval[2][1]=0.994;
  tighterval[3][1]=5;

// [3, 4]
  tighterval[0][2]=0.01;
  tighterval[1][2]=-0.00018;
  tighterval[2][2]=0.993;
  tighterval[3][2]=5;

// [4, 5]
  tighterval[0][3]=0.02;
  tighterval[1][3]=-0.00017;
  tighterval[2][3]=0.99;
  tighterval[3][3]=5;



// [5, 6]
  tighterval[0][4]=0.02;
  tighterval[1][4]=-0.00017;
  tighterval[2][4]=0.99;
  tighterval[3][4]=5;

// [6, 7]
  tighterval[0][5]=0.02;
  tighterval[1][5]=-0.00016;
  tighterval[2][5]=2.0;
  tighterval[3][5]=5;

// [7, 8]
  tighterval[0][6]=0.02;
  tighterval[1][6]=-0.00016;
  tighterval[2][6]=0.99;
  tighterval[3][6]=4.5;

 // [8, 10]
  tighterval[0][7]=0.02;
  tighterval[1][7]=-0.0001;
  tighterval[2][7]=0.99;
  tighterval[3][7]=4.5;

// [10, 12]
  tighterval[0][8]=0.03;
  tighterval[1][8]=-0.0001;
  tighterval[2][8]=0.99;
  tighterval[3][8]=4;

// [12, 16]
  tighterval[0][9]=0.03;
  tighterval[1][9]=-0.00006;
  tighterval[2][9]=0.99;
  tighterval[3][9]=3.5;


// [16, 20]
  tighterval[0][10]=0.03;
  tighterval[1][10]=-0.0001;
  tighterval[2][10]=0.99;
  tighterval[3][10]=2;


// [20, 24]
  tighterval[0][11]=0.03;
  tighterval[1][11]=-0.0001;
  tighterval[2][11]=0.99;
  tighterval[3][11]=2;


// >24
  tighterval[0][12]=0.03;
  tighterval[1][12]=-0.0001;
  tighterval[2][12]=0.99;
  tighterval[3][12]=2;


  TString name="";
  Int_t arrdim=dim*nptbins;
  cout<<"Will save "<<arrdim<<" TParameter<float>"<<endl;
  TClonesArray max("TParameter<float>",arrdim);
  for(Int_t ival=0;ival<dim;ival++){
    for(Int_t jpt=0;jpt<nptbins;jpt++){
      name=Form("par%dptbin%d",ival,jpt);
      cout<<"Setting "<<name.Data()<<" to "<<tighterval[ival][jpt]<<endl;
      new(max[jpt*dim+ival])TParameter<float>(name.Data(),tighterval[ival][jpt]);
    }
  }

  Bool_t pidflag=kTRUE;
  RDHFDStartoKpipi->SetUsePID(pidflag);
  if(pidflag) cout<<"PID is used"<<endl;
  else cout<<"PID is not used"<<endl;

  // PID SETTINGS
  AliAODPidHF* pidObj=new AliAODPidHF();
  // pidObj->SetName("pid4DSatr");
  Int_t mode=1;
  Double_t priors[5]={0.01,0.001,0.3,0.3,0.3};
  pidObj->SetPriors(priors);
  pidObj->SetMatch(mode);
  pidObj->SetSigma(0,3); // TPC ----> previously it was 3sig, but this reduces background
  pidObj->SetSigma(3,3); // TOF
  pidObj->SetTPC(kTRUE);
  pidObj->SetTOF(kTRUE);
  RDHFDStartoKpipi->SetPidHF(pidObj);


  TFile *fFlat=TFile::Open("./CentrDistrBins005.root","READ");
  TCanvas *c=fFlat->Get("cintegral");
  TH1F *hfl=(TH1F*)c->FindObject("hint");
  RDHFDStartoKpipi->SetHistoForCentralityFlattening(hfl,40.,50.,0.,0);

  //activate pileup rejection
  //RDHFDStartoKpipi->SetOptPileup(AliRDHFCuts::kRejectPileupEvent);

  TFile* fout=new TFile("cuts4SignifMaximDStar_07.root","recreate");   //set this!!
  fout->cd();
  RDHFDStartoKpipi->Write();
  max.Write();
  fout->Close();

}
