
void runGridTrain( Bool_t kUseMC = kFALSE,
		   Int_t runNumber = 117112,
		   TString dataSet = "LHC10b",
		   TString pass = "pass3",   
		   Bool_t mergeJDL = kFALSE, //you should use the kTRUE when runnign full mode 
		   TString gridMode = "test", 
		   TString workDir = "", //set your working directory
		   TString rootVer = "v5-34-08-6", 
		   TString alirootVer = "vAN-20141012")
{  
  // Load common libraries
  gSystem->Load("libCore.so");  
  gSystem->Load("libTree.so");
  gSystem->Load("libGeom.so");
  gSystem->Load("libVMC.so");
  gSystem->Load("libPhysics.so");

  gSystem->Load("libSTEERBase");
  gSystem->Load("libESD");
  gSystem->Load("libAOD");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice"); 
  gSystem->Load("libCORRFW");  

  //Use AliRoot includes to compile our task

  gROOT->ProcessLine(".include $ALICE_ROOT/include");
             
  gSystem->Load("libPWGmuon");
  gSystem->Load("libPWGTools");  

  //You will also need to load the jet-analysis framework
  //consult Twiki or people from PWG-JE

  // Create and configure the alien handler plugin
  gROOT->LoadMacro("CreateAlienHandler.C");
  AliAnalysisGrid *alienHandler = CreateAlienHandler(runNumber, dataSet, pass, gridMode, kUseMC, rootVer, alirootVer, mergeJDL, workDir);  
  if (!alienHandler) return;

  // Create the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("testAnalysis");

  // Connect plug-in to the analysis manager
  mgr->SetGridHandler(alienHandler);
  
  //In case of Train Analysis
  AliESDInputHandler *esdHandler = new AliESDInputHandler();
  //  esdHandler->SetReadFriends(kFALSE);
  mgr->SetInputEventHandler(esdHandler);

  if(kUseMC)
    {
      AliMCEventHandler* mcHandler = new AliMCEventHandler();
      mgr->SetMCtruthEventHandler(mcHandler);
      //mcHandler->SetReadTR(kTRUE);
    }
  
  //Needed for pID
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
  AliAnalysisTaskPIDResponse *taskpidresponse = AddTaskPIDResponse(kUseMC);

  // Here load your analyses classe

  // Enable debug printouts
  mgr->SetDebugLevel(2);

  if (!mgr->InitAnalysis())
    return;

  mgr->PrintStatus();
  // Start analysis in grid.
  mgr->StartAnalysis("grid");
};
