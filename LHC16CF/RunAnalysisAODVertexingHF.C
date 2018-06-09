class AliAnalysisGrid;
class AliAnalysisAlien;

void RunAnalysisAODVertexingHF()
{
  //
  // Test macro for AliAnalysisTaskSE's for heavy-flavour candidates
  // It has the structure of a Analysis Train:
  // - in this macro, change things related to running mode
  //   and input preparation 
  // - add your task using a AddTaskXXX macro 
  //
  // A.Dainese, andrea.dainese@lnl.infn.it
  // "grid" mode added by R.Bala, bala@to.infn.it
  //


  gSystem->SetIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_PHYSICS -I$ALICE_PHYSICS/include -I$ALICE_ROOT/ITS -I$ALICE_ROOT/TPC -I$ALICE_ROOT/CONTAINERS -I$ALICE_ROOT/STEER/STEER -I$ALICE_ROOT/STEER/STEERBase -I$ALICE_ROOT/STEER/ESD -I$ALICE_ROOT/STEER/AOD -I$ALICE_ROOT/TRD -I$ALICE_ROOT/macros -I$ALICE_ROOT/ANALYSIS  -I$ALICE_PHYSICS/OADB -I$ALICE_PHYSICS/PWGHF -I$ALICE_PHYSICS/PWGHF/base -I$ALICE_PHYSICS/PWGHF/vertexingHF -I$ALICE_PHYSICS/PWG/FLOW/Base -I$ALICE_PHYSICS/PWG/FLOW/Tasks -I$ALICE_PHYSICS/PWG -g"); 
  //
  TString trainName = "D2H";
  TString analysisMode = "grid"; // "local", "grid", or "proof"
  TString inputMode    = "list"; // "list", "xml", or "dataset"
  Long64_t nentries=99999999,firstentry=0;
  Bool_t useParFiles=kFALSE;
  Bool_t useAlienPlugin=kTRUE;
  TString pluginmode="terminate";
  TString testfileslistWithPlugin="";
  Bool_t saveProofToAlien=kFALSE;
  TString proofOutdir = "";
  TString loadMacroPath="$ALICE/../../../AliPhysics/PWGHF/vertexingHF/macros/";
  //TString localDataPath="data/";
  //TString loadMacroPath="./"; // this is normally needed for CAF
  //

  if(analysisMode=="grid") {
    // Connect to AliEn
    TGrid::Connect("alien://");
  } else if(analysisMode=="proof") {
    // Connect to the PROOF cluster
    if(inputMode!="dataset") {printf("Input mode must be dataset, for proof analysis\n"); return;}
    gEnv->SetValue("XSec.GSI.DelegProxy","2");
    TProof::Open("alicecaf");
    //TProof::Reset("alicecaf");
    if(saveProofToAlien) {
      TGrid::Connect("alien://");
      if(gGrid) {
	TString homedir = gGrid->GetHomeDirectory();
	TString workdir = homedir + trainName;
	if(!gGrid->Cd(workdir)) {
	  gGrid->Cd(homedir);
	  if(gGrid->Mkdir(workdir)) {
	    gGrid->Cd(trainName);
	    ::Info("VertexingTrain::Connect()", "Directory %s created", gGrid->Pwd());
	  }
	}	   
	gGrid->Mkdir("proof_output");
	gGrid->Cd("proof_output");
	proofOutdir = Form("alien://%s", gGrid->Pwd());
      } 
    }
  }


  // AliRoot libraries
  if(analysisMode=="local" || analysisMode=="grid") {
    TString loadLibraries="LoadLibraries.C"; loadLibraries.Prepend(loadMacroPath.Data());
    gROOT->LoadMacro(loadLibraries.Data());
    LoadLibraries(useParFiles);
    gSystem->Load("libGui");
    gSystem->Load("libRAWDatabase");
    gSystem->Load("libCDB");
    gSystem->Load("libSTEER");
    gSystem->Load("libTRDbase");
    gSystem->Load("libPWGTRD");
  } else if (analysisMode=="proof") {
    gSystem->Load("libTree");
    gSystem->Load("libGeom");
    gSystem->Load("libPhysics");
    gSystem->Load("libVMC");
    gSystem->Load("libMinuit");
    gSystem->Load("libGui");
    gSystem->Load("libRAWDatabase");
    gSystem->Load("libCDB");
    gSystem->Load("libSTEER");
    gSystem->Load("libTRDbase");
    gSystem->Load("libPWGTRD");
    // Enable the needed packages
    //gProof->ClearPackages();
    TString parDir="/afs/cern.ch/user/d/dainesea/code/";
    TString parFile;
    if(!useParFiles) {
      gProof->UploadPackage("AF-v4-17");
      gProof->EnablePackage("AF-v4-17");
      // --- Enable the PWGHFvertexingHF Package
      parFile="PWGHFvertexingHF.par"; parFile.Prepend(parDir.Data());
      gProof->UploadPackage(parFile.Data());
      gProof->EnablePackage("PWGHFvertexingHF");
    } else {
      // --- Enable the STEERBase Package
      parFile="STEERBase.par"; parFile.Prepend(parDir.Data());
      gProof->UploadPackage(parFile.Data());
      gProof->EnablePackage("STEERBase");
      // --- Enable the ESD Package
      parFile="ESD.par"; parFile.Prepend(parDir.Data());
      gProof->UploadPackage(parFile.Data());
      gProof->EnablePackage("ESD");
      // --- Enable the AOD Package
      parFile="AOD.par"; parFile.Prepend(parDir.Data());
      gProof->UploadPackage(parFile.Data());
      gProof->EnablePackage("AOD");
      // --- Enable the ANALYSIS Package
      parFile="ANALYSIS.par"; parFile.Prepend(parDir.Data());
      gProof->UploadPackage(parFile.Data());
      gProof->EnablePackage("ANALYSIS");
      // --- Enable the ANALYSISalice Package
      parFile="ANALYSISalice.par"; parFile.Prepend(parDir.Data());
      gProof->UploadPackage(parFile.Data());
      gProof->EnablePackage("ANALYSISalice");
      // --- Enable the CORRFW Package
      parFile="CORRFW.par"; parFile.Prepend(parDir.Data());
      gProof->UploadPackage(parFile.Data());
      gProof->EnablePackage("CORRFW");
      // --- Enable the PWGHFbase Package
      parFile="PWGHFbase.par"; parFile.Prepend(parDir.Data());
      gProof->UploadPackage(parFile.Data());
      gProof->EnablePackage("PWGHFbase");
      // --- Enable the PWGHFvertexingHF Package
      parFile="PWGHFvertexingHF.par"; parFile.Prepend(parDir.Data());
      gProof->UploadPackage(parFile.Data());
      gProof->EnablePackage("PWGHFvertexingHF");
    }
    gProof->ShowEnabledPackages(); // show a list of enabled packages
  }


  // Create Alien plugin, if requested
  if(useAlienPlugin) {  
    //    if(analysisMode!="grid") {printf("Analysis mode must be grid, to use alien plugin\n"); return;}
    AliAnalysisGrid *alienHandler = CreateAlienHandler(pluginmode,useParFiles,testfileslistWithPlugin);  
    if(!alienHandler) return;
  }


  //-------------------------------------------------------------------
  // Prepare input
  TChain *chainAOD = 0;
  
  TString dataset; // for proof

  if(!useAlienPlugin) {
    //TString makeAODInputChain="../MakeAODInputChain.C"; makeAODInputChain.Prepend(loadMacroPath.Data());
    if(inputMode=="list") {
      // Local files
      TChain *chainAOD = new TChain("aodTree");
      TChain *chainAODfriend = new TChain("aodTree");
      chainAOD->Add("data/AliAOD.root");
      chainAODfriend->Add("data/AliAOD.VertexingHF.root");
      chainAOD->AddFriend(chainAODfriend);
      //gROOT->LoadMacro(makeAODInputChain.Data());
      //chainAOD = MakeAODInputChain();// with this it reads ./AliAOD.root and ./AliAOD.VertexingHF.root
      //chainAOD = MakeAODInputChain("alien://alice/cern.ch/user/n/nseptian/newtrain/out_lhc08x/180100/",1,1);
      printf("ENTRIES %d\n",chainAOD->GetEntries());
    } else if(inputMode=="xml") {
      // xml
      gROOT->LoadMacro(makeAODInputChain.Data());
      chainAOD = MakeAODInputChain("collection_aod.xml","collection_aodHF.xml");
    } else if(inputMode=="dataset") {
      // CAF dataset
      //gProof->ShowDataSets();
      dataset="/ITS/dainesea/AODVertexingHF_LHC08x_180100";
    }
  }

  // Create the analysis manager
  AliAnalysisManager *mgr  = new AliAnalysisManager("My Manager","My Manager");
  mgr->SetDebugLevel(10);
  // Connect plug-in to the analysis manager
  if(useAlienPlugin) mgr->SetGridHandler(alienHandler);

  // Input
  AliAODInputHandler *inputHandler = new AliAODInputHandler("handler","handler for D2H");
  if(analysisMode=="proof" ) {
    inputHandler->AddFriend("./AliAOD.VertexingHF.root");
    //inputHandler->AddFriend("deltas/AliAOD.VertexingHF.root");
    if(saveProofToAlien) mgr->SetSpecialOutputLocation(proofOutdir);
  }
  mgr->SetInputEventHandler(inputHandler);
  //-------------------------------------------------------------------

  
  //-------------------------------------------------------------------
  // Analysis tasks (wagons of the train)   
  //
  // First add the task for the PID response setting
 gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
 AliAnalysisTaskSE *setupTask = AddTaskPIDResponse(kTRUE,kTRUE,kTRUE,"1");
/*
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDqa.C");
  AliAnalysisTaskPIDqa *pidQA = AddTaskPIDqa();
*/

  TString taskName;
  
  ////// ADD THE FULL D2H TRAIN
  /*taskName="../AddD2HTrain.C"; taskName.Prepend(loadMacroPath.Data());
  gROOT->LoadMacro(taskName.Data());
  Bool_t readMC=kFALSE;
  AddD2HTrain(readMC);//,1,0,0,0,0,0,0,0,0,0,0);*/
  
  ////// OR ADD INDIVIDUAL TASKS
  /*
  taskName="AddTaskHFQA.C"; taskName.Prepend(loadMacroPath.Data());
  gROOT->LoadMacro(taskName.Data());
  AliAnalysisTaskSEHFQA *qaTask = AddTaskHFQA(2,"DStartoKpipiCuts.root",kFALSE, kFALSE, 0 , "test",kTRUE,kTRUE,kFALSE, kTRUE, kFALSE,kTRUE);
*/
  /*taskName="AddTaskHFQA.C"; taskName.Prepend(loadMacroPath.Data());
  gROOT->LoadMacro(taskName.Data());
  AliAnalysisTaskSEHFQA *qaTask = AddTaskHFQA(2,"DStartoKpipiCuts.root",kFALSE, kFALSE, 0 , "QApp",kTRUE,kTRUE,kTRUE, kTRUE, kFALSE,kTRUE);
*/
//AddTaskHFQA(AliAnalysisTaskSEHFQA::DecChannel ch,TString filecutsname="DStartoKpipiCuts.root",Bool_t readMC=kFALSE, Bool_t simplemode=kFALSE, Int_t system=0 /*0=pp, 1=PbPb*/, TString finDirname="",Bool_t trackon=kTRUE,Bool_t pidon=kTRUE,Bool_t centralityon=kTRUE, Bool_t eventselon=kTRUE, Bool_t flowobson=kFALSE, Bool_t filldistribforeffcheckson=kFALSE)

/*
   taskName="AddTaskVertexingHF.C"; taskName.Prepend(loadMacroPath.Data());
   gROOT->LoadMacro(taskName.Data());
   AliAnalysisTaskSEVertexingHF *hfTask = AddTaskVertexingHF();
*/
  /*  taskName="AddTaskCompareHF.C"; taskName.Prepend(loadMacroPath.Data());
    gROOT->LoadMacro(taskName.Data());
    AliAnalysisTaskSECompareHF *cmpTask = AddTaskCompareHF();
   
    taskName="AddTaskD0Mass.C"; taskName.Prepend(loadMacroPath.Data());
    gROOT->LoadMacro(taskName.Data());
    AliAnalysisTaskSED0Mass *d0massTask = AddTaskD0Mass();
    AliAnalysisTaskSED0Mass *d0massLikeSignTask = AddTaskD0Mass(1); 

    taskName="AddTaskDStarSpectra.C"; taskName.Prepend(loadMacroPath.Data());
    gROOT->LoadMacro(taskName.Data());
    AliAnalysisTaskSEDStarSpectra *dStarSpectraTask = AddTaskDStarSpectra(0,0,100,"DStartoKpipiCuts.root","nseptian",kFALSE,kFALSE);
*/

/*
    gROOT->LoadMacro("$PWD/AliAnalysisMyTaskSEDvsMultiplicity.cxx++g");
    taskName="AddTaskDvsMultiplicity.C"; taskName.Prepend(loadMacroPath.Data());
    gROOT->LoadMacro(taskName.Data());
    AliAnalysisMyTaskSEDvsMultiplicity *dMultTask = AddTaskDvsMultiplicity(0,kTRUE,0,413,"pp13TeV_kINT7","DStartoKpipiCuts.root","DStartoKpipiCuts","estimators.root",11.98,kFALSE,0,AliAnalysisMyTaskSEDvsMultiplicity::kNtrk10,AliAnalysisMyTaskSEDvsMultiplicity::kEta10,kFALSE,16);
*/

    gROOT->LoadMacro("MyAliCFTaskVertexingHF.cxx++g");
    taskName="MyAddTaskCFVertexingHFCascade.C";// taskName.Prepend(loadMacroPath.Data());
    gROOT->LoadMacro(taskName.Data());
    MyAliCFTaskVertexingHF *CFtask = AddTaskCFVertexingHFCascade("DStartoKpipiCuts.root","DStartoKpipiCuts","pp13TeV_kINT7",MyAliCFTaskVertexingHF::kCheetah,kFALSE,kFALSE,413,2,kTRUE,kFALSE,kFALSE,kFALSE,kFALSE,"estimators.root",MyAliCFTaskVertexingHF::kNtrk10,kTRUE,kFALSE,11.98,kFALSE);
    
    /*
    taskName="AddTaskCFDStar.C"; taskName.Prepend(loadMacroPath.Data()); //test addtaskcfdstar
    gROOT->LoadMacro(taskName.Data());
    AliCFTaskForDStarAnalysis *dstarCFTask = AddTaskCFDStar();
    */
    /*
    taskName="AddTaskDStarCharmFraction.C"; taskName.Prepend(loadMacroPath.Data());
    gROOT->LoadMacro(taskName.Data());
    AliAnalysisTaskSEDStarCharmFraction *dstarCharmFractionTask = AddTaskDStarCharmFraction("DStartoKpipiCuts.root",kFALSE,"",3.,6.,9.,kFALSE,0.,"AnalysisResults.root");
    */
/*
    taskName="AddTaskD0Mass.C"; taskName.Prepend(loadMacroPath.Data());
    gROOT->LoadMacro(taskName.Data());
    AliAnalysisTaskSED0Mass *d0massTask = AddTaskD0Mass(0,kFALSE,kFALSE,kFALSE,0,0,0,0,"Loose","D0toKpiCutsppNoNoRecVtx(null)PileupRej.root","D0toKpiCuts",kFALSE,false,false,false,false,false,true,1);

*/  	
/*
    taskName="AddTaskDplus.C"; taskName.Prepend(loadMacroPath.Data());
    gROOT->LoadMacro(taskName.Data());
    AliAnalysisTaskSEDplus *dplusTask = AddTaskDplus();

    taskName="AddTaskDs.C"; taskName.Prepend(loadMacroPath.Data());
    gROOT->LoadMacro(taskName.Data());
    AliAnalysisTaskSEDs *dsTask = AddTaskDs();
    
    //taskName="AddTaskSelectHF.C"; taskName.Prepend(loadMacroPath.Data());
    //gROOT->LoadMacro(taskName.Data());
    //AliAnalysisTaskSESelectHF *seleTask = AddTaskSelectHF();
    
    taskName="AddTaskBkgLikeSignD0.C"; taskName.Prepend(loadMacroPath.Data());
    gROOT->LoadMacro(taskName.Data());
    AliAnalysisTaskSEBkgLikeSignD0 *lsD0Task = AddTaskBkgLikeSignD0();
    
    taskName="AddTaskCFMultiVarMultiStep.C"; taskName.Prepend(loadMacroPath.Data());
    gROOT->LoadMacro(taskName.Data());
    AliCFHeavyFlavourTaskMultiVarMultiStep *cfmvmsTask = AddTaskCFMultiVarMultiStep();
    

    taskName="AddTaskSECharmFraction.C"; 
    taskName.Prepend(loadMacroPath.Data());
    gROOT->LoadMacro(taskName.Data());
    Int_t switchMC[5]={0,0,0,0,0};
    Int_t ppPbPb=1;// 0 for pp, 1 for PbPb, used to siwtch on/off the removal of daughters from the primary vertex
    AliAnalysisTaskSECharmFraction *cFractTask = AddTaskSECharmFraction("standard",switchMC,readMC,kTRUE,kFALSE,"D0toKpiCharmFractCuts.root","c",ppPbPb);
    // arguments: filename,switchMC,readmc,usepid,likesign,cutfilename,containerprefix

    
    // attach a private task (not committed)
    // (the files MyTask.h MyTask.cxx AddMyTask.C have to be declared in plugin
    // configuration, see below)
    
    if(analysisMode.Data()=="proof") {
    gProof->LoadMacro("MyTask.cxx++g");
    } else {
    gROOT->LoadMacro("MyTask.cxx++g");
    }
    gROOT->LoadMacro("AddMyTask.C");
    MyTask *myTask = AddMyTask();
    
    
    if(analysisMode.Data()=="proof") {
    gProof->LoadMacro("AliDStarJets.cxx++g");
    } else {
    gROOT->LoadMacro("AliDStarJets.cxx++g");
    }
    gROOT->LoadMacro("AddTaskDStarJets.C");
    AliDStarJets *myTask = AddTaskDStarJets();
  */
  //-------------------------------------------------------------------
  
  //
  // Run the analysis
  //    
  if(chainAOD) printf("CHAIN HAS %d ENTRIES\n",(Int_t)chainAOD->GetEntries());
  
  if(!mgr->InitAnalysis()) return;
  mgr->PrintStatus();
  if(analysisMode=="grid" && !useAlienPlugin) analysisMode="local";
  if(analysisMode!="proof") {
    mgr->StartAnalysis(analysisMode.Data(),chainAOD,nentries,firstentry);
  } else {
    // proof
    mgr->StartAnalysis(analysisMode.Data(),dataset.Data(),nentries,firstentry);
  }
  
  return;
}
//_____________________________________________________________________________
//
AliAnalysisGrid* CreateAlienHandler(TString pluginmode,Bool_t useParFiles=kFALSE, TString testfileslistWithPlugin="")
{
  // Check if user has a valid token, otherwise make one. This has limitations.
  // One can always follow the standard procedure of calling alien-token-init then
  //   source /tmp/gclient_env_$UID in the current shell.
   AliAnalysisAlien *plugin = new AliAnalysisAlien();
   // Set the run mode (can be "full", "test", "offline", "submit" or "terminate")
   plugin->SetRunMode(pluginmode.Data());
   plugin->SetUser("nseptian");
   // Set versions of used packages
   plugin->SetAPIVersion("V1.1x");
   //plugin->SetROOTVersion("v5-34-30-alice6-2");
   //plugin->SetAliROOTVersion("v5-08-13m-cookdedx-1");
   plugin->SetNtestFiles(1);
   plugin->SetAliPhysicsVersion("vAN-20180310-1");//("v5-09-02-01-1");
   gROOT->LoadMacro("../AddGoodRuns.C");

   // Declare input data to be processed.
   //************************************************
   // Set data file list to test on local mode
   //************************************************  
   //plugin->SetFileForTestMode(testfileslistWithPlugin.Data());

   //************************************************
   // Set data search pattern for DATA
   //************************************************  
   //Method 1: To create automatically xml through plugin /*/ 
   //plugin->SetGridDataDir("/alice/data/2016/LHC16l"); // specify LHC period
   // OR plugin->SetFriendChainName("deltas/AliAOD.VertexingHF.root");
   // Adds only the good runs from the Monalisa Run Condition Table
   // More than one period can be added but the period name has to be removed from GridDataDir (to be tested)
   Int_t totruns=0;
   //totruns += AddGoodRuns(plugin,"LHC10b"); // specify LHC period
   //totruns += AddGoodRuns(plugin,"LHC10c"); // specify LHC period
   //totruns += AddGoodRuns(plugin,"LHC16k");
   //totruns += AddGoodRuns(plugin,"LHC16l");

   plugin->SetGridDataDir("/alice/sim/2017/LHC17d20a1");
   plugin->SetDataPattern("/AOD/*/*AliAOD.root"); // specify AOD set
   plugin->SetFriendChainName("./AliAOD.VertexingHF.root");
   totruns += AddGoodRuns(plugin,"LHC16k","LHC17d20a1");
   plugin->SetNrunsPerMaster(totruns);
   // Method 2: Declare existing data files (e.g xml collections)

   //plugin->AddDataFile("/alice/cern.ch/user/r/rbala/000168068_000170593.xml");
   //  plugin->SetDataPattern("*AliAOD.root");
   //  plugin->SetFriendChainName("./AliAOD.VertexingHF.root"); 

   //************************************************
   // Set data search pattern for MONTECARLO
   //************************************************
    
   //plugin->SetGridDataDir("/alice/sim/2016/LHC16j2a1"); // specify MC sample
   //plugin->SetDataPattern("/AOD/*/*AliAOD.root"); // specify AOD set
   //plugin->SetFriendChainName("./AliAOD.VertexingHF.root");
   // OR plugin->SetFriendChainName("deltas/AliAOD.VertexingHF.root");
   // Adds only the good runs from the Monalisa Run Condition Table 
   // More than one period can be added!
   //Int_t totruns=0;
   //totruns += AddGoodRuns(plugin,"LHC16l","LHC16j2a2"); // specify LHC period for anchor runs; and the name of the MC production
   //totruns += AddGoodRuns(plugin,"LHC16k","LHC16j2a1"); // specify LHC period for anchor runs;  and the name of the MC production
   //totruns += AddGoodRuns(plugin,"LHC10d","LHC10f7"); // specify LHC period for anchor runs;  and the name of the MC production
   //plugin->SetNrunsPerMaster(totruns);
   
   //
   // Define alien work directory where all files will be copied. Relative to alien $HOME.
   plugin->SetGridWorkingDir("DStarvsMultiplicityAnalysisMC");
   // Name of executable
   plugin->SetExecutable("DStarvsMultiplicityAnalysisLHC16kMC.sh");
   // Declare alien output directory. Relative to working directory.
   plugin->SetGridOutputDir("LHC17d20a1_full_CorrectionMyDStarvsMultMC_pp13TeV_promptD"); // In this case will be $HOME/work/output
   // Declare the analysis source files names separated by blancs. To be compiled runtime
   // using ACLiC on the worker nodes.
   plugin->SetAnalysisSource("MyAliCFTaskVertexingHF.cxx");
   // Declare all libraries (other than the default ones for the framework. These will be
   // loaded by the generated analysis macro. Add all extra files (task .cxx/.h) here.
   plugin->SetAdditionalLibs("libPWGflowBase.so libPWGflowTasks.so libPWGHFbase.so libPWGHFvertexingHF.so libGui.so libRAWDatabase.so libCDB.so libSTEER.so libTRDbase.so libPWGTRD.so MyAliCFTaskVertexingHF.cxx MyAliCFTaskVertexingHF.h");
  
   // use par files
   if(useParFiles) {
     plugin->EnablePackage("STEERBase.par");
     plugin->EnablePackage("ESD.par");
     plugin->EnablePackage("AOD.par");
     plugin->EnablePackage("ANALYSIS.par");
     plugin->EnablePackage("OADB.par");
     plugin->EnablePackage("ANALYSISalice.par");
     plugin->EnablePackage("CORRFW.par");
     plugin->EnablePackage("PWGHFbase.par");
     plugin->EnablePackage("PWGHFvertexingHF.par");
   }
   plugin->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_PHYSICS -I$ALICE_PHYSICS/include -I$ALICE_ROOT/ITS -I$ALICE_ROOT/TPC -I$ALICE_ROOT/CONTAINERS -I$ALICE_ROOT/STEER/STEER -I$ALICE_ROOT/STEER/STEERBase -I$ALICE_ROOT/STEER/ESD -I$ALICE_ROOT/STEER/AOD -I$ALICE_ROOT/TRD -I$ALICE_ROOT/macros -I$ALICE_ROOT/ANALYSIS  -I$ALICE_PHYSICS/OADB -I$ALICE_PHYSICS/PWGHF -I$ALICE_PHYSICS/PWGHF/base -I$ALICE_PHYSICS/PWGHF/vertexingHF -I$ALICE_PHYSICS/PWG/FLOW/Base -I$ALICE_PHYSICS/PWG/FLOW/Tasks -I$ALICE_PHYSICS/PWG -g"); 

   plugin->SetDefaultOutputs(kTRUE);
   // merging via jdl
   plugin->SetMergeViaJDL(kFALSE);
   plugin->SetOneStageMerging(kTRUE);
   plugin->SetMaxMergeStages(2);


   // Optionally set a name for the generated analysis macro (default MyAnalysis.C)
   plugin->SetAnalysisMacro("AnalysisDStarvsMultiplicityMC.C");
   plugin->SetTTL(86000);
   //plugin->SetRunPrefix("000");
   //plugin->AddRunNumber(256504);
   //plugin->SetNMasterJobs(5);
   plugin->SetSplitMaxInputFileNumber(50);
   //plugin->SetRunRange(256504,258537);
   // Optionally set maximum number of input files/subjob (default 100, put 0 to ignore)
   // Optionally modify the name of the generated JDL (default analysis.jdl)
   plugin->SetJDLName("TaskHF.jdl");

   return plugin;
}