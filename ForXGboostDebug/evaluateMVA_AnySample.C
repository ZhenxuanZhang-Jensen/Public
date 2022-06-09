#include "TROOT.h"
#include "TKey.h"
#include "TFile.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TObjArray.h"
#include "THStack.h"
#include "TLegend.h"
#include "TEfficiency.h"
#include "TGraphAsymmErrors.h"
#include "TF1.h"
#include "TMath.h"
#include "TCut.h"
#include "TPaletteAxis.h"

#include "TMVA/Reader.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <string>
#include <algorithm>

void evaluateMVA_AnySample(){
  
  //Any Test sample: data, Hgg signal, bkg
  string treeFileName_sig = "/eos/user/z/zhenxuan/Hgg_mass/MiniTree/output_sig125.root";
  string treeFileName_pp = "/eos/user/z/zhenxuan/Hgg_mass/MiniTree/New_MCpp_DataDriven_QCD_SFs_sEoEWgt_2DpTWgt.root";
  string treeFileName_pfff = "/eos/user/z/zhenxuan/Hgg_mass/MiniTree/New_DataDriven_QCD_SFs_sEoEWgt_2DpTWgt.root";
  // string outputFileName = "/eos/user/z/zhenxuan/Hgg_mass/MiniTree/TMVA_diphoMVA.root"; // For TMVA
  string outputFileName = "/eos/user/z/zhenxuan/Hgg_mass/MiniTree/XGboost_diphoMVA_t8.root"; // For XGboost

 //Weight file
  // string weightfile = "/afs/cern.ch/user/z/zhenxuan/CMSSW_10_6_20/src/Hgg_mass/BDT/XGboost/convert_pkl2xml/convert_v3/DiphotonXGboost_pro2.weights.xml";
  string weightfile = "/afs/cern.ch/user/z/zhenxuan/CMSSW_10_6_20/src/Hgg_mass/BDT/XGboost/weights.xml";
  // string weightfile = "/afs/cern.ch/user/z/zhenxuan/CMSSW_10_6_20/src/Hgg_mass/BDT/XGboost/convert_pkl2xml/convert_v3/DiphotonXGboost_v7.weights.xml";
  // string weightfile = "/afs/cern.ch/user/z/zhenxuan/CMSSW_10_6_20/src/Hgg_mass/BDT/XGboost/altDiphoModel_IncUL2017.xml";
  // string weightfile = "/afs/cern.ch/user/z/zhenxuan/CMSSW_10_6_20/src/Hgg_mass/BDT/XGboost/XGboost_diphoMVA_testv1.weights.xml";

  ///
  TFile *treeFile_sig = new TFile(treeFileName_sig.c_str());
  TFile *treeFile_pp = new TFile(treeFileName_pp.c_str());
  TFile *treeFile_pfff = new TFile(treeFileName_pfff.c_str());

  TTree *t_sig = (TTree*)treeFile_sig->Get("/tagsDumper/trees/ggh_125_13TeV_UntaggedTag");
  TTree *t_pp = (TTree*)treeFile_pp->Get("pp");
  TTree *t_pfff = (TTree*)treeFile_pfff->Get("DataDriven_QCD");

  Long64_t nEntries_sig = t_sig->GetEntries();
  Long64_t nEntries_pp = t_pp->GetEntries();
  Long64_t nEntries_pfff = t_pfff->GetEntries();

  //output file create
  TFile *outputFile = new TFile (outputFileName.c_str(),"RECREATE");
  TTree *outputTree_sig = t_sig->CloneTree(0); //new TTree("sigPhoTree","sigPhoTree");
  TTree *outputTree_pp = t_pp->CloneTree(0); //new TTree("bkgPhoTree","bkgPhoTree");
  TTree *outputTree_pfff = t_pfff->CloneTree(0); //new TTree("bkgPhoTree","bkgPhoTree");

  //input vars: or equivalent variables in the test sample
  Float_t leadmva;
  Float_t subleadmva;
  Float_t leadptom;
  Float_t subleadptom;
  Float_t leadeta;
  Float_t subleadeta;
  Float_t CosPhi;
  Float_t vtxprob;
  Float_t sigmarv;
  Float_t sigmawv;





  //output mva
  Float_t DiphotonMVA_self;
  std::vector<float> MVA_vec;

  outputTree_sig->Branch("DiphotonMVA_self",&DiphotonMVA_self);
  outputTree_pp->Branch("DiphotonMVA_self",&DiphotonMVA_self);
  outputTree_pfff->Branch("DiphotonMVA_self",&DiphotonMVA_self);
  
  //====MVA==
  TMVA::Reader *diphotonMva = new TMVA::Reader("!V");
  //Load the var in the same order as training or in weight file, with the first var name used in TMVA training and the second from the test samples, can be different name; then book the MVA 
  diphotonMva->AddVariable("leadmva",&leadmva);
  diphotonMva->AddVariable("subleadmva",&subleadmva);
  diphotonMva->AddVariable("leadptom",&leadptom);
  diphotonMva->AddVariable("subleadptom",&subleadptom);
  diphotonMva->AddVariable("leadeta",&leadeta);
  diphotonMva->AddVariable("subleadeta",&subleadeta);
  diphotonMva->AddVariable("CosPhi",&CosPhi);
  diphotonMva->AddVariable("vtxprob",&vtxprob);
  diphotonMva->AddVariable("sigmarv",&sigmarv);
  diphotonMva->AddVariable("sigmawv",&sigmawv);
  diphotonMva->BookMVA("BDT",weightfile.c_str());
 
  //varibles in input sig trees
  t_sig->SetBranchAddress("leadmva",&leadmva);
  t_sig->SetBranchAddress("subleadmva",&subleadmva);
  t_sig->SetBranchAddress("leadptom",&leadptom);
  t_sig->SetBranchAddress("subleadptom",&subleadptom);
  t_sig->SetBranchAddress("leadeta",&leadeta);
  t_sig->SetBranchAddress("subleadeta",&subleadeta);
  t_sig->SetBranchAddress("CosPhi",&CosPhi);
  t_sig->SetBranchAddress("vtxprob",&vtxprob);
  t_sig->SetBranchAddress("sigmarv",&sigmarv);
  t_sig->SetBranchAddress("sigmawv",&sigmawv);
 

  //varibles in input bkg trees
  t_pp->SetBranchAddress("leadmva",&leadmva);
  t_pp->SetBranchAddress("subleadmva",&subleadmva);
  t_pp->SetBranchAddress("leadptom",&leadptom);
  t_pp->SetBranchAddress("subleadptom",&subleadptom);
  t_pp->SetBranchAddress("leadeta",&leadeta);
  t_pp->SetBranchAddress("subleadeta",&subleadeta);
  t_pp->SetBranchAddress("CosPhi",&CosPhi);
  t_pp->SetBranchAddress("vtxprob",&vtxprob);
  t_pp->SetBranchAddress("sigmarv",&sigmarv);
  t_pp->SetBranchAddress("sigmawv",&sigmawv);

  t_pfff->SetBranchAddress("leadmva",&leadmva);
  t_pfff->SetBranchAddress("subleadmva",&subleadmva);
  t_pfff->SetBranchAddress("leadptom",&leadptom);
  t_pfff->SetBranchAddress("subleadptom",&subleadptom);
  t_pfff->SetBranchAddress("leadeta",&leadeta);
  t_pfff->SetBranchAddress("subleadeta",&subleadeta);
  t_pfff->SetBranchAddress("CosPhi",&CosPhi);
  t_pfff->SetBranchAddress("vtxprob",&vtxprob);
  t_pfff->SetBranchAddress("sigmarv",&sigmarv);
  t_pfff->SetBranchAddress("sigmawv",&sigmawv);
 
  for(int i = 0; i < nEntries_sig; i++){
  // for(int i = 0; i < 100000; i++){
    t_sig->GetEntry(i);
    DiphotonMVA_self = -999.;
    // MVA_vec = diphotonMva->EvaluateMulticlass( "BDT" ); 
    DiphotonMVA_self = diphotonMva->EvaluateMVA( "BDT" ); 
    // DiphotonMVA_self = MVA_vec[0];
    if (i < 10){
      std::cout << "debug value sig mva scores:" << DiphotonMVA_self << std::endl;
    }
    outputTree_sig->Fill();
  }

  cout <<"signal evaluated!" << endl;

  for(int i = 0; i < nEntries_pp; i++){
  // for(int i = 0; i < 100000; i++){
    t_pp->GetEntry(i);
    DiphotonMVA_self = -999.;
    // MVA_vec = diphotonMva->EvaluateMulticlass( "BDT" ); 
    DiphotonMVA_self = diphotonMva->EvaluateMVA( "BDT" ); 
    // DiphotonMVA_self = MVA_vec[1];
    if (i < 10){
      std::cout << "debug value bkg pp mva scores:" << DiphotonMVA_self << std::endl;

    }
    outputTree_pp->Fill();    
  }

  cout <<"background pp evaluated!" << endl;

  for(int i = 0; i < nEntries_pfff; i++){
  // for(int i = 0; i < 100000; i++){
    t_pfff->GetEntry(i);
    DiphotonMVA_self = -999.;
    // MVA_vec = diphotonMva->EvaluateMulticlass( "BDT" ); 
    DiphotonMVA_self = diphotonMva->EvaluateMVA( "BDT" ); 
    // DiphotonMVA_self = MVA_vec[0];
    outputTree_pfff->Fill();    
  }

  cout <<"background pfff evaluated!" << endl;

  outputFile->Write();

}

