#include <iostream>
#include <sstream>
#include <istream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cmath>
#include <functional>
#include <vector>
#include <cassert>
#include "TMath.h"
#include "TRandom.h"

#include "parsePileUpJSON2.h"
#include "SMPJ/AnalysisFW/plugins/AnalysisTemplate.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"

//#include "../Unfolding/RooUnfold-1.1.1/src/RooUnfold.h"
//#include "../Unfolding/RooUnfold-1.1.1/src/RooUnfoldResponse.h"

using namespace std;

//---------------------------- Constructor Of The Class TriggerTurnOn -------------------------- //
Analysis_Template_MC::Analysis_Template_MC(edm::ParameterSet const& cfg)
{
  mFileName       = cfg.getParameter<std::vector<std::string>>            ("filename");
  mTreeName       = cfg.getParameter<std::string>               ("treename");
  mDirName        = cfg.getParameter<std::string>               ("dirname");
  
  mMinPt     = cfg.getParameter<double> ("minPt");
  mYMax      = cfg.getParameter<double> ("ymax");
  mJetID     = cfg.getParameter<int>  ("JetID");
  
  mprintOk   = cfg.getParameter<int>  ("printOk");

 //    mIsMCarlo       = cfg.getUntrackedParameter<bool>             ("isMCarlo");
 //    mJECUncSrcNames = cfg.getParameter<std::vector<std::string> > ("jecUncSrcNames");
  //jecs = new JECs(mIsMCarlo, mGlobalTag, mjettype);
  mGlobalTag        = cfg.getParameter<std::string>               ("pseudoglobaltag");
  mjettype        = cfg.getParameter<std::string>               ("jettype");
  
  mIsMCarlo       = cfg.getUntrackedParameter<bool>             ("isMCarlo");
  mJECUncSrcNames = cfg.getParameter<std::vector<std::string> > ("jecUncSrcNames");
  mJECUncSrc = cfg.getParameter<std::string> ("jecUncSrc");
  
  mMCSlice= cfg.getParameter<int>("MCSlice");
  mPUReweighting= cfg.getUntrackedParameter<bool>("PUReweighting");
  mLowPileUp= cfg.getUntrackedParameter<bool>("LowPileUp");

}

//------------------------------ Declaration Of The Function beginjob() ------------------------//
void Analysis_Template_MC::beginJob()
 {

   //for(unsigned i=0;i<mFileName.size();i++){
   //mInf = TFile::Open(mFileName[i].c_str());
   //mDir = (TDirectoryFile*)mInf->Get(mDirName.c_str());  
   //TTree* mTreeFirst=(TTree*)mDir->Get(mTreeName.c_str());
   //mTree.push_back(mTreeFirst)
   //Event = new QCDEvent();
   //TBranch *branch = mTree[i]->GetBranch("events");
   //branch->SetAddress(&Event);  
   //}
   
   //jecs = new JECs(mIsMCarlo, mGlobalTag, mjettype,mJECUncSrcNames);
   
   //Pt binning
   double Ptbinning[73] = {18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84,97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468,507, 548, 592, 638, 686, 737, 790, 846, 905, 967,1032, 1101, 1172, 1248, 1327, 1410, 1497, 1588, 1684, 1784, 1890, 2000,2116, 2238, 2366, 2500, 2640, 2787, 2941, 3103, 3273, 3450, 3637, 3832,4037, 4252, 4477, 4713, 4961, 5220, 5492, 5777, 6076, 6389, 6717, 7000};

   int Ptbins=72;
   
   //------------------ Histogram Booking --------------------------- //
   num_of_Vtx     = fs->make<TH1F>("num_of_Vtx","num_of_Vtx",100,0.,100.);
   num_of_VtxGood = fs->make<TH1F>("num_of_VtxGood","num_of_VtxGood",100,0.,100.);
   
   mc_pthat = fs->make<TH1F>("mc_pthat","mc_pthat",200,0.,2000.);
   mc_pthat_weighted = fs->make<TH1F>("mc_pthat_weighted","mc_pthat_weighted",100,0.,2000.);
   mc_pthat_weighted->Sumw2();
   
   pt0_GENJet  = fs->make<TH1F>("pt0_GENJet","pt0_GENJet",Ptbins,Ptbinning); pt0_GENJet->Sumw2();
   pt1_GENJet  = fs->make<TH1F>("pt1_GENJet","pt1_GENJet",Ptbins,Ptbinning); pt1_GENJet->Sumw2();
   y0_GENJet = fs->make<TH1F>("y0_GENJet","y0_GENJet",60,-5.,5.); y0_GENJet->Sumw2();
   y1_GENJet = fs->make<TH1F>("y1_GENJet","y1_GENJet",60,-5.,5.); y1_GENJet->Sumw2();
   phi0_GENJet = fs->make<TH1F>("phi0_GENJet","phi0_GENJet",60, -TMath::Pi(),TMath::Pi()); phi0_GENJet->Sumw2();
   phi1_GENJet = fs->make<TH1F>("phi1_GENJet","phi1_GENJet",60, -TMath::Pi(),TMath::Pi()); phi1_GENJet->Sumw2();
   
   pt0_DETJet  = fs->make<TH1F>("pt0_DETJet","pt0_DETJet",Ptbins,Ptbinning); pt0_DETJet->Sumw2();
   pt1_DETJet  = fs->make<TH1F>("pt1_DETJet","pt1_DETJet",Ptbins,Ptbinning); pt1_DETJet->Sumw2();
   pt0_DETJetUncor  = fs->make<TH1F>("pt0_DETJetUncor","pt0_DETJetUncor",Ptbins,Ptbinning); pt0_DETJetUncor->Sumw2();
   pt1_DETJetUncor  = fs->make<TH1F>("pt1_DETJetUncor","pt1_DETJetUncor",Ptbins,Ptbinning); pt1_DETJetUncor->Sumw2();
   y0_DETJet = fs->make<TH1F>("y0_DETJet","y0_DETJet",60,-5.,5.); y0_DETJet->Sumw2();
   y1_DETJet = fs->make<TH1F>("y1_DETJet","y1_DETJet",60,-5.,5.); y1_DETJet->Sumw2();
   phi0_DETJet = fs->make<TH1F>("phi0_DETJet","phi0_DETJet",60, -TMath::Pi(),TMath::Pi()); phi0_DETJet->Sumw2();
   phi1_DETJet = fs->make<TH1F>("phi1_DETJet","phi1_DETJet",60, -TMath::Pi(),TMath::Pi()); phi1_DETJet->Sumw2();

   pt0_DETBTagJet  = fs->make<TH1F>("pt0_DETBTagJet","pt0_DETBTagJet",Ptbins,Ptbinning); pt0_DETBTagJet->Sumw2();
   pt1_DETBTagJet  = fs->make<TH1F>("pt1_DETBTagJet","pt1_DETBTagJet",Ptbins,Ptbinning); pt1_DETBTagJet->Sumw2();
   y0_DETBTagJet = fs->make<TH1F>("y0_DETBTagJet","y0_DETBTagJet",60,-5.,5.); y0_DETBTagJet->Sumw2();
   y1_DETBTagJet = fs->make<TH1F>("y1_DETBTagJet","y1_DETBTagJet",60,-5.,5.); y1_DETBTagJet->Sumw2();

   ptPair_DETBTagJet  = fs->make<TH1F>("ptPair_DETBTagJet","ptPair_DETBTagJet",Ptbins,Ptbinning); ptPair_DETBTagJet->Sumw2();

   pt0_DETBTagJetInclusive  = fs->make<TH1F>("pt0_DETJetInclusive","pt0_DETJetInclusive",Ptbins,Ptbinning); pt0_DETBTagJetInclusive->Sumw2();
   
   Multiplicity_DETJet = fs->make<TH1F>("Multiplicity_DETJet","Multiplicity_DETJet",21, -0.5,20.5); Multiplicity_DETJet->Sumw2();
   Multiplicity_GENJet = fs->make<TH1F>("Multiplicity_GENJet","Multiplicity_GENJet",21, -0.5,20.5); Multiplicity_GENJet->Sumw2();

   Multiplicity_DETJetBTags = fs->make<TH1F>("Multiplicity_DETJetBTags","Multiplicity_DETJetBTags",21, -0.5,20.5); Multiplicity_DETJetBTags->Sumw2();
   Multiplicity_DETJetBTagsInclusive = fs->make<TH1F>("Multiplicity_DETJetBTagsInclusive","Multiplicity_DETJetBTagsInclusive",21, -0.5,20.5); Multiplicity_DETJetBTagsInclusive->Sumw2();
   
   pt0_DETInclJet  = fs->make<TH1F>("pt0_DETInclJet","pt0_DETInclJet",Ptbins,Ptbinning); pt0_DETInclJet->Sumw2();
   pt0_DETInclJetUncor  = fs->make<TH1F>("pt0_DETInclJetUncor","pt0_DETInclJetUncor",Ptbins,Ptbinning); pt0_DETInclJetUncor->Sumw2();
   y0_DETInclJet = fs->make<TH1F>("y0_DETInclJet","y0_DETInclJet",60,-5.,5.); y0_DETInclJet->Sumw2();
   phi0_DETInclJet = fs->make<TH1F>("phi0_DETInclJet","phi0_DETInclJet",60, -TMath::Pi(),TMath::Pi()); phi0_DETInclJet->Sumw2();
   
   PF_MatchedInclusiveJets  = fs->make<TH1F>("PF_MatchedInclusiveJets","PF_MatchedInclusiveJets",Ptbins,Ptbinning); PF_MatchedInclusiveJets->Sumw2(); 
   Gen_MatchedInclusiveJets  = fs->make<TH1F>("Gen_MatchedInclusiveJets","Gen_MatchedInclusiveJets",Ptbins,Ptbinning); Gen_MatchedInclusiveJets->Sumw2(); 
   TwoD_MatchedInclusiveJets   = fs->make<TH2F>("TwoD_MatchedInclusiveJets","TwoD_MatchedInclusiveJets",Ptbins,Ptbinning,Ptbins,Ptbinning); TwoD_MatchedInclusiveJets->Sumw2(); 
   PF_FakeInclusiveJets   = fs->make<TH1F>("PF_FakeInclusiveJets","PF_FakeInclusiveJets",Ptbins,Ptbinning); PF_FakeInclusiveJets->Sumw2(); 
   Gen_MissInclusiveJets   = fs->make<TH1F>("Gen_MissInclusiveJets","Gen_MissInclusiveJets",Ptbins,Ptbinning); Gen_MissInclusiveJets->Sumw2(); 
   
   PileUpVSVertex  = fs->make<TH2F>("PileUpVSVertex","PileUpVSVertex",31,-0.5,30.5,31,-0.5,30.5); PileUpVSVertex->Sumw2(); 

   DeltaR_Jets = fs->make<TH1F>("DeltaR_Jets","DeltaR_Jets",20, 0,10); DeltaR_Jets->Sumw2();
   
   AcceptancePtJets= fs->make<TH1F>("AcceptancePtJets","AcceptancePtJets",Ptbins,Ptbinning); AcceptancePtJets->Sumw2();
   PurityPtJets= fs->make<TH1F>("PurityPtJets","PurityPtJets",Ptbins,Ptbinning); PurityPtJets->Sumw2();
   BackgroundPtJets= fs->make<TH1F>("BackgroundPtJets","BackgroundPtJets",Ptbins,Ptbinning); BackgroundPtJets->Sumw2();
   StabilityPtJets= fs->make<TH1F>("StabilityPtJets","StabilityPtJets",Ptbins,Ptbinning); StabilityPtJets->Sumw2();
   
   AcceptancePtJets_1bin= fs->make<TH1F>("AcceptancePtJets_1bin","AcceptancePtJets_1bin",Ptbins,Ptbinning); 
   PurityPtJets_1bin= fs->make<TH1F>("PurityPtJets_1bin","PurityPtJets_1bin",Ptbins,Ptbinning); 
   AcceptancePtJets_2bin= fs->make<TH1F>("AcceptancePtJets_2bin","AcceptancePtJets_2bin",Ptbins,Ptbinning); 
   PurityPtJets_2bin= fs->make<TH1F>("PurityPtJets_2bin","PurityPtJets_2bin",Ptbins,Ptbinning); 
   AcceptancePtJets_3bin= fs->make<TH1F>("AcceptancePtJets_3bin","AcceptancePtJets_3bin",Ptbins,Ptbinning); 
   PurityPtJets_3bin= fs->make<TH1F>("PurityPtJets_3bin","PurityPtJets_3bin",Ptbins,Ptbinning); 
   AcceptancePtJets_4bin= fs->make<TH1F>("AcceptancePtJets_4bin","AcceptancePtJets_4bin",Ptbins,Ptbinning); 
   PurityPtJets_4bin= fs->make<TH1F>("PurityPtJets_4bin","PurityPtJets_4bin",Ptbins,Ptbinning); 
   AcceptancePtJets_5bin= fs->make<TH1F>("AcceptancePtJets_5bin","AcceptancePtJets_5bin",Ptbins,Ptbinning); 
   PurityPtJets_5bin= fs->make<TH1F>("PurityPtJets_5bin","PurityPtJets_5bin",Ptbins,Ptbinning); 
   AcceptancePtJets_6bin= fs->make<TH1F>("AcceptancePtJets_6bin","AcceptancePtJets_6bin",Ptbins,Ptbinning); 
   PurityPtJets_6bin= fs->make<TH1F>("PurityPtJets_6bin","PurityPtJets_6bin",Ptbins,Ptbinning); 
   AcceptancePtJets_7bin= fs->make<TH1F>("AcceptancePtJets_7bin","AcceptancePtJets_7bin",Ptbins,Ptbinning); 
   PurityPtJets_7bin= fs->make<TH1F>("PurityPtJets_7bin","PurityPtJets_7bin",Ptbins,Ptbinning); 
   
   BackgroundPtJets_1bin= fs->make<TH1F>("BackgroundPtJets_1bin","BackgroundPtJets_1bin",Ptbins,Ptbinning); 
   StabilityPtJets_1bin= fs->make<TH1F>("StabilityPtJets_1bin","StabilityPtJets_1bin",Ptbins,Ptbinning); 
   BackgroundPtJets_2bin= fs->make<TH1F>("BackgroundPtJets_2bin","BackgroundPtJets_2bin",Ptbins,Ptbinning); 
   StabilityPtJets_2bin= fs->make<TH1F>("StabilityPtJets_2bin","StabilityPtJets_2bin",Ptbins,Ptbinning); 
   BackgroundPtJets_3bin= fs->make<TH1F>("BackgroundPtJets_3bin","BackgroundPtJets_3bin",Ptbins,Ptbinning); 
   StabilityPtJets_3bin= fs->make<TH1F>("StabilityPtJets_3bin","StabilityPtJets_3bin",Ptbins,Ptbinning); 
   BackgroundPtJets_4bin= fs->make<TH1F>("BackgroundPtJets_4bin","BackgroundPtJets_4bin",Ptbins,Ptbinning); 
   StabilityPtJets_4bin= fs->make<TH1F>("StabilityPtJets_4bin","StabilityPtJets_4bin",Ptbins,Ptbinning); 
   BackgroundPtJets_5bin= fs->make<TH1F>("BackgroundPtJets_5bin","BackgroundPtJets_5bin",Ptbins,Ptbinning); 
   StabilityPtJets_5bin= fs->make<TH1F>("StabilityPtJets_5bin","StabilityPtJets_5bin",Ptbins,Ptbinning); 
   BackgroundPtJets_6bin= fs->make<TH1F>("BackgroundPtJets_6bin","BackgroundPtJets_6bin",Ptbins,Ptbinning); 
   StabilityPtJets_6bin= fs->make<TH1F>("StabilityPtJets_6bin","StabilityPtJets_6bin",Ptbins,Ptbinning); 
   BackgroundPtJets_7bin= fs->make<TH1F>("BackgroundPtJets_7bin","BackgroundPtJets_7bin",Ptbins,Ptbinning); 
   StabilityPtJets_7bin= fs->make<TH1F>("StabilityPtJets_7bin","StabilityPtJets_7bin",Ptbins,Ptbinning); 
   
   pt_DETInclJet_1bin  = fs->make<TH1F>("pt_DETInclJet_1bin","pt_DETInclJet_1bin",Ptbins,Ptbinning); pt_DETInclJet_1bin->Sumw2();
   pt_DETInclJet_2bin  = fs->make<TH1F>("pt_DETInclJet_2bin","pt_DETInclJet_2bin",Ptbins,Ptbinning); pt_DETInclJet_2bin->Sumw2();
   pt_DETInclJet_3bin  = fs->make<TH1F>("pt_DETInclJet_3bin","pt_DETInclJet_3bin",Ptbins,Ptbinning); pt_DETInclJet_3bin->Sumw2();
   pt_DETInclJet_4bin  = fs->make<TH1F>("pt_DETInclJet_4bin","pt_DETInclJet_4bin",Ptbins,Ptbinning); pt_DETInclJet_4bin->Sumw2();
   pt_DETInclJet_5bin  = fs->make<TH1F>("pt_DETInclJet_5bin","pt_DETInclJet_5bin",Ptbins,Ptbinning); pt_DETInclJet_5bin->Sumw2();
   pt_DETInclJet_6bin  = fs->make<TH1F>("pt_DETInclJet_6bin","pt_DETInclJet_6bin",Ptbins,Ptbinning); pt_DETInclJet_6bin->Sumw2();
   pt_DETInclJet_7bin  = fs->make<TH1F>("pt_DETInclJet_7bin","pt_DETInclJet_7bin",Ptbins,Ptbinning); pt_DETInclJet_7bin->Sumw2();

   pt_DETBTagJet_1bin  = fs->make<TH1F>("pt_DETBTagJet_1bin","pt_DETBTagJet_1bin",Ptbins,Ptbinning); pt_DETBTagJet_1bin->Sumw2();
   pt_DETBTagJet_2bin  = fs->make<TH1F>("pt_DETBTagJet_2bin","pt_DETBTagJet_2bin",Ptbins,Ptbinning); pt_DETBTagJet_2bin->Sumw2();
   pt_DETBTagJet_3bin  = fs->make<TH1F>("pt_DETBTagJet_3bin","pt_DETBTagJet_3bin",Ptbins,Ptbinning); pt_DETBTagJet_3bin->Sumw2();
   pt_DETBTagJet_4bin  = fs->make<TH1F>("pt_DETBTagJet_4bin","pt_DETBTagJet_4bin",Ptbins,Ptbinning); pt_DETBTagJet_4bin->Sumw2();
   pt_DETBTagJet_5bin  = fs->make<TH1F>("pt_DETBTagJet_5bin","pt_DETBTagJet_5bin",Ptbins,Ptbinning); pt_DETBTagJet_5bin->Sumw2();
   pt_DETBTagJet_6bin  = fs->make<TH1F>("pt_DETBTagJet_6bin","pt_DETBTagJet_6bin",Ptbins,Ptbinning); pt_DETBTagJet_6bin->Sumw2();
   pt_DETBTagJet_7bin  = fs->make<TH1F>("pt_DETBTagJet_7bin","pt_DETBTagJet_7bin",Ptbins,Ptbinning); pt_DETBTagJet_7bin->Sumw2();

   pt_DETInclJetCrossSectNorm_1bin  = fs->make<TH1F>("pt_DETInclJetCrossSectNorm_1bin","pt_DETInclJetCrossSectNorm_1bin",Ptbins,Ptbinning); pt_DETInclJetCrossSectNorm_1bin->Sumw2();
   pt_DETInclJetCrossSectNorm_2bin  = fs->make<TH1F>("pt_DETInclJetCrossSectNorm_2bin","pt_DETInclJetCrossSectNorm_2bin",Ptbins,Ptbinning); pt_DETInclJetCrossSectNorm_2bin->Sumw2();
   pt_DETInclJetCrossSectNorm_3bin  = fs->make<TH1F>("pt_DETInclJetCrossSectNorm_3bin","pt_DETInclJetCrossSectNorm_3bin",Ptbins,Ptbinning); pt_DETInclJetCrossSectNorm_3bin->Sumw2();
   pt_DETInclJetCrossSectNorm_4bin  = fs->make<TH1F>("pt_DETInclJetCrossSectNorm_4bin","pt_DETInclJetCrossSectNorm_4bin",Ptbins,Ptbinning); pt_DETInclJetCrossSectNorm_4bin->Sumw2();
   pt_DETInclJetCrossSectNorm_5bin  = fs->make<TH1F>("pt_DETInclJetCrossSectNorm_5bin","pt_DETInclJetCrossSectNorm_5bin",Ptbins,Ptbinning); pt_DETInclJetCrossSectNorm_5bin->Sumw2();
   pt_DETInclJetCrossSectNorm_6bin  = fs->make<TH1F>("pt_DETInclJetCrossSectNorm_6bin","pt_DETInclJetCrossSectNorm_6bin",Ptbins,Ptbinning); pt_DETInclJetCrossSectNorm_6bin->Sumw2();
   pt_DETInclJetCrossSectNorm_7bin  = fs->make<TH1F>("pt_DETInclJetCrossSectNorm_7bin","pt_DETInclJetCrossSectNorm_7bin",Ptbins,Ptbinning); pt_DETInclJetCrossSectNorm_7bin->Sumw2();

   pt_DETInclJet40_1bin  = fs->make<TH1F>("pt_DETInclJet40_1bin","pt_DETInclJet40_1bin",Ptbins,Ptbinning); pt_DETInclJet40_1bin->Sumw2();
   pt_DETInclJet40_2bin  = fs->make<TH1F>("pt_DETInclJet40_2bin","pt_DETInclJet40_2bin",Ptbins,Ptbinning); pt_DETInclJet40_2bin->Sumw2();
   pt_DETInclJet40_3bin  = fs->make<TH1F>("pt_DETInclJet40_3bin","pt_DETInclJet40_3bin",Ptbins,Ptbinning); pt_DETInclJet40_3bin->Sumw2();
   pt_DETInclJet40_4bin  = fs->make<TH1F>("pt_DETInclJet40_4bin","pt_DETInclJet40_4bin",Ptbins,Ptbinning); pt_DETInclJet40_4bin->Sumw2();
   pt_DETInclJet40_5bin  = fs->make<TH1F>("pt_DETInclJet40_5bin","pt_DETInclJet40_5bin",Ptbins,Ptbinning); pt_DETInclJet40_5bin->Sumw2();
   pt_DETInclJet40_6bin  = fs->make<TH1F>("pt_DETInclJet40_6bin","pt_DETInclJet40_6bin",Ptbins,Ptbinning); pt_DETInclJet40_6bin->Sumw2();
   pt_DETInclJet40_7bin  = fs->make<TH1F>("pt_DETInclJet40_7bin","pt_DETInclJet40_7bin",Ptbins,Ptbinning); pt_DETInclJet40_7bin->Sumw2();


   pt_DETInclJet60_1bin  = fs->make<TH1F>("pt_DETInclJet60_1bin","pt_DETInclJet60_1bin",Ptbins,Ptbinning); pt_DETInclJet60_1bin->Sumw2();
   pt_DETInclJet60_2bin  = fs->make<TH1F>("pt_DETInclJet60_2bin","pt_DETInclJet60_2bin",Ptbins,Ptbinning); pt_DETInclJet60_2bin->Sumw2();
   pt_DETInclJet60_3bin  = fs->make<TH1F>("pt_DETInclJet60_3bin","pt_DETInclJet60_3bin",Ptbins,Ptbinning); pt_DETInclJet60_3bin->Sumw2();
   pt_DETInclJet60_4bin  = fs->make<TH1F>("pt_DETInclJet60_4bin","pt_DETInclJet60_4bin",Ptbins,Ptbinning); pt_DETInclJet60_4bin->Sumw2();
   pt_DETInclJet60_5bin  = fs->make<TH1F>("pt_DETInclJet60_5bin","pt_DETInclJet60_5bin",Ptbins,Ptbinning); pt_DETInclJet60_5bin->Sumw2();
   pt_DETInclJet60_6bin  = fs->make<TH1F>("pt_DETInclJet60_6bin","pt_DETInclJet60_6bin",Ptbins,Ptbinning); pt_DETInclJet60_6bin->Sumw2();
   pt_DETInclJet60_7bin  = fs->make<TH1F>("pt_DETInclJet60_7bin","pt_DETInclJet60_7bin",Ptbins,Ptbinning); pt_DETInclJet60_7bin->Sumw2();

   pt_DETInclJet80_1bin  = fs->make<TH1F>("pt_DETInclJet80_1bin","pt_DETInclJet80_1bin",Ptbins,Ptbinning); pt_DETInclJet80_1bin->Sumw2();
   pt_DETInclJet80_2bin  = fs->make<TH1F>("pt_DETInclJet80_2bin","pt_DETInclJet80_2bin",Ptbins,Ptbinning); pt_DETInclJet80_2bin->Sumw2();
   pt_DETInclJet80_3bin  = fs->make<TH1F>("pt_DETInclJet80_3bin","pt_DETInclJet80_3bin",Ptbins,Ptbinning); pt_DETInclJet80_3bin->Sumw2();
   pt_DETInclJet80_4bin  = fs->make<TH1F>("pt_DETInclJet80_4bin","pt_DETInclJet80_4bin",Ptbins,Ptbinning); pt_DETInclJet80_4bin->Sumw2();
   pt_DETInclJet80_5bin  = fs->make<TH1F>("pt_DETInclJet80_5bin","pt_DETInclJet80_5bin",Ptbins,Ptbinning); pt_DETInclJet80_5bin->Sumw2();
   pt_DETInclJet80_6bin  = fs->make<TH1F>("pt_DETInclJet80_6bin","pt_DETInclJet80_6bin",Ptbins,Ptbinning); pt_DETInclJet80_6bin->Sumw2();
   pt_DETInclJet80_7bin  = fs->make<TH1F>("pt_DETInclJet80_7bin","pt_DETInclJet80_7bin",Ptbins,Ptbinning); pt_DETInclJet80_7bin->Sumw2();

   pt_DETInclJet140_1bin  = fs->make<TH1F>("pt_DETInclJet140_1bin","pt_DETInclJet140_1bin",Ptbins,Ptbinning); pt_DETInclJet140_1bin->Sumw2();
   pt_DETInclJet140_2bin  = fs->make<TH1F>("pt_DETInclJet140_2bin","pt_DETInclJet140_2bin",Ptbins,Ptbinning); pt_DETInclJet140_2bin->Sumw2();
   pt_DETInclJet140_3bin  = fs->make<TH1F>("pt_DETInclJet140_3bin","pt_DETInclJet140_3bin",Ptbins,Ptbinning); pt_DETInclJet140_3bin->Sumw2();
   pt_DETInclJet140_4bin  = fs->make<TH1F>("pt_DETInclJet140_4bin","pt_DETInclJet140_4bin",Ptbins,Ptbinning); pt_DETInclJet140_4bin->Sumw2();
   pt_DETInclJet140_5bin  = fs->make<TH1F>("pt_DETInclJet140_5bin","pt_DETInclJet140_5bin",Ptbins,Ptbinning); pt_DETInclJet140_5bin->Sumw2();
   pt_DETInclJet140_6bin  = fs->make<TH1F>("pt_DETInclJet140_6bin","pt_DETInclJet140_6bin",Ptbins,Ptbinning); pt_DETInclJet140_6bin->Sumw2();
   pt_DETInclJet140_7bin  = fs->make<TH1F>("pt_DETInclJet140_7bin","pt_DETInclJet140_7bin",Ptbins,Ptbinning); pt_DETInclJet140_7bin->Sumw2();

   pt_DETInclJet200_1bin  = fs->make<TH1F>("pt_DETInclJet200_1bin","pt_DETInclJet200_1bin",Ptbins,Ptbinning); pt_DETInclJet200_1bin->Sumw2();
   pt_DETInclJet200_2bin  = fs->make<TH1F>("pt_DETInclJet200_2bin","pt_DETInclJet200_2bin",Ptbins,Ptbinning); pt_DETInclJet200_2bin->Sumw2();
   pt_DETInclJet200_3bin  = fs->make<TH1F>("pt_DETInclJet200_3bin","pt_DETInclJet200_3bin",Ptbins,Ptbinning); pt_DETInclJet200_3bin->Sumw2();
   pt_DETInclJet200_4bin  = fs->make<TH1F>("pt_DETInclJet200_4bin","pt_DETInclJet200_4bin",Ptbins,Ptbinning); pt_DETInclJet200_4bin->Sumw2();
   pt_DETInclJet200_5bin  = fs->make<TH1F>("pt_DETInclJet200_5bin","pt_DETInclJet200_5bin",Ptbins,Ptbinning); pt_DETInclJet200_5bin->Sumw2();
   pt_DETInclJet200_6bin  = fs->make<TH1F>("pt_DETInclJet200_6bin","pt_DETInclJet200_6bin",Ptbins,Ptbinning); pt_DETInclJet200_6bin->Sumw2();
   pt_DETInclJet200_7bin  = fs->make<TH1F>("pt_DETInclJet200_7bin","pt_DETInclJet200_7bin",Ptbins,Ptbinning); pt_DETInclJet200_7bin->Sumw2();

   pt_DETInclJet260_1bin  = fs->make<TH1F>("pt_DETInclJet260_1bin","pt_DETInclJet260_1bin",Ptbins,Ptbinning); pt_DETInclJet260_1bin->Sumw2();
   pt_DETInclJet260_2bin  = fs->make<TH1F>("pt_DETInclJet260_2bin","pt_DETInclJet260_2bin",Ptbins,Ptbinning); pt_DETInclJet260_2bin->Sumw2();
   pt_DETInclJet260_3bin  = fs->make<TH1F>("pt_DETInclJet260_3bin","pt_DETInclJet260_3bin",Ptbins,Ptbinning); pt_DETInclJet260_3bin->Sumw2();
   pt_DETInclJet260_4bin  = fs->make<TH1F>("pt_DETInclJet260_4bin","pt_DETInclJet260_4bin",Ptbins,Ptbinning); pt_DETInclJet260_4bin->Sumw2();
   pt_DETInclJet260_5bin  = fs->make<TH1F>("pt_DETInclJet260_5bin","pt_DETInclJet260_5bin",Ptbins,Ptbinning); pt_DETInclJet260_5bin->Sumw2();
   pt_DETInclJet260_6bin  = fs->make<TH1F>("pt_DETInclJet260_6bin","pt_DETInclJet260_6bin",Ptbins,Ptbinning); pt_DETInclJet260_6bin->Sumw2();
   pt_DETInclJet260_7bin  = fs->make<TH1F>("pt_DETInclJet260_7bin","pt_DETInclJet260_7bin",Ptbins,Ptbinning); pt_DETInclJet260_7bin->Sumw2();

   pt_DETInclJet320_1bin  = fs->make<TH1F>("pt_DETInclJet320_1bin","pt_DETInclJet320_1bin",Ptbins,Ptbinning); pt_DETInclJet320_1bin->Sumw2();
   pt_DETInclJet320_2bin  = fs->make<TH1F>("pt_DETInclJet320_2bin","pt_DETInclJet320_2bin",Ptbins,Ptbinning); pt_DETInclJet320_2bin->Sumw2();
   pt_DETInclJet320_3bin  = fs->make<TH1F>("pt_DETInclJet320_3bin","pt_DETInclJet320_3bin",Ptbins,Ptbinning); pt_DETInclJet320_3bin->Sumw2();
   pt_DETInclJet320_4bin  = fs->make<TH1F>("pt_DETInclJet320_4bin","pt_DETInclJet320_4bin",Ptbins,Ptbinning); pt_DETInclJet320_4bin->Sumw2();
   pt_DETInclJet320_5bin  = fs->make<TH1F>("pt_DETInclJet320_5bin","pt_DETInclJet320_5bin",Ptbins,Ptbinning); pt_DETInclJet320_5bin->Sumw2();
   pt_DETInclJet320_6bin  = fs->make<TH1F>("pt_DETInclJet320_6bin","pt_DETInclJet320_6bin",Ptbins,Ptbinning); pt_DETInclJet320_6bin->Sumw2();
   pt_DETInclJet320_7bin  = fs->make<TH1F>("pt_DETInclJet320_7bin","pt_DETInclJet320_7bin",Ptbins,Ptbinning); pt_DETInclJet320_7bin->Sumw2();

   pt_DETInclJet400_1bin  = fs->make<TH1F>("pt_DETInclJet400_1bin","pt_DETInclJet400_1bin",Ptbins,Ptbinning); pt_DETInclJet400_1bin->Sumw2();
   pt_DETInclJet400_2bin  = fs->make<TH1F>("pt_DETInclJet400_2bin","pt_DETInclJet400_2bin",Ptbins,Ptbinning); pt_DETInclJet400_2bin->Sumw2();
   pt_DETInclJet400_3bin  = fs->make<TH1F>("pt_DETInclJet400_3bin","pt_DETInclJet400_3bin",Ptbins,Ptbinning); pt_DETInclJet400_3bin->Sumw2();
   pt_DETInclJet400_4bin  = fs->make<TH1F>("pt_DETInclJet400_4bin","pt_DETInclJet400_4bin",Ptbins,Ptbinning); pt_DETInclJet400_4bin->Sumw2();
   pt_DETInclJet400_5bin  = fs->make<TH1F>("pt_DETInclJet400_5bin","pt_DETInclJet400_5bin",Ptbins,Ptbinning); pt_DETInclJet400_5bin->Sumw2();
   pt_DETInclJet400_6bin  = fs->make<TH1F>("pt_DETInclJet400_6bin","pt_DETInclJet400_6bin",Ptbins,Ptbinning); pt_DETInclJet400_6bin->Sumw2();
   pt_DETInclJet400_7bin  = fs->make<TH1F>("pt_DETInclJet400_7bin","pt_DETInclJet400_7bin",Ptbins,Ptbinning); pt_DETInclJet400_7bin->Sumw2();

   pt_DETInclJet450_1bin  = fs->make<TH1F>("pt_DETInclJet450_1bin","pt_DETInclJet450_1bin",Ptbins,Ptbinning); pt_DETInclJet450_1bin->Sumw2();
   pt_DETInclJet450_2bin  = fs->make<TH1F>("pt_DETInclJet450_2bin","pt_DETInclJet450_2bin",Ptbins,Ptbinning); pt_DETInclJet450_2bin->Sumw2();
   pt_DETInclJet450_3bin  = fs->make<TH1F>("pt_DETInclJet450_3bin","pt_DETInclJet450_3bin",Ptbins,Ptbinning); pt_DETInclJet450_3bin->Sumw2();
   pt_DETInclJet450_4bin  = fs->make<TH1F>("pt_DETInclJet450_4bin","pt_DETInclJet450_4bin",Ptbins,Ptbinning); pt_DETInclJet450_4bin->Sumw2();
   pt_DETInclJet450_5bin  = fs->make<TH1F>("pt_DETInclJet450_5bin","pt_DETInclJet450_5bin",Ptbins,Ptbinning); pt_DETInclJet450_5bin->Sumw2();
   pt_DETInclJet450_6bin  = fs->make<TH1F>("pt_DETInclJet450_6bin","pt_DETInclJet450_6bin",Ptbins,Ptbinning); pt_DETInclJet450_6bin->Sumw2();
   pt_DETInclJet450_7bin  = fs->make<TH1F>("pt_DETInclJet450_7bin","pt_DETInclJet450_7bin",Ptbins,Ptbinning); pt_DETInclJet450_7bin->Sumw2();
  
   pt_GENInclJet_1bin  = fs->make<TH1F>("pt_GENInclJet_1bin","pt_GENInclJet_1bin",Ptbins,Ptbinning); pt_GENInclJet_1bin->Sumw2();
   pt_GENInclJet_2bin  = fs->make<TH1F>("pt_GENInclJet_2bin","pt_GENInclJet_2bin",Ptbins,Ptbinning); pt_GENInclJet_2bin->Sumw2();
   pt_GENInclJet_3bin  = fs->make<TH1F>("pt_GENInclJet_3bin","pt_GENInclJet_3bin",Ptbins,Ptbinning); pt_GENInclJet_3bin->Sumw2();
   pt_GENInclJet_4bin  = fs->make<TH1F>("pt_GENInclJet_4bin","pt_GENInclJet_4bin",Ptbins,Ptbinning); pt_GENInclJet_4bin->Sumw2();
   pt_GENInclJet_5bin  = fs->make<TH1F>("pt_GENInclJet_5bin","pt_GENInclJet_5bin",Ptbins,Ptbinning); pt_GENInclJet_5bin->Sumw2();
   pt_GENInclJet_6bin  = fs->make<TH1F>("pt_GENInclJet_6bin","pt_GENInclJet_6bin",Ptbins,Ptbinning); pt_GENInclJet_6bin->Sumw2();
   pt_GENInclJet_7bin  = fs->make<TH1F>("pt_GENInclJet_7bin","pt_GENInclJet_7bin",Ptbins,Ptbinning); pt_GENInclJet_7bin->Sumw2();

   pt_GENBTagJet_1bin  = fs->make<TH1F>("pt_GENBTagJet_1bin","pt_GENBTagJet_1bin",Ptbins,Ptbinning); pt_GENBTagJet_1bin->Sumw2();
   pt_GENBTagJet_2bin  = fs->make<TH1F>("pt_GENBTagJet_2bin","pt_GENBTagJet_2bin",Ptbins,Ptbinning); pt_GENBTagJet_2bin->Sumw2();
   pt_GENBTagJet_3bin  = fs->make<TH1F>("pt_GENBTagJet_3bin","pt_GENBTagJet_3bin",Ptbins,Ptbinning); pt_GENBTagJet_3bin->Sumw2();
   pt_GENBTagJet_4bin  = fs->make<TH1F>("pt_GENBTagJet_4bin","pt_GENBTagJet_4bin",Ptbins,Ptbinning); pt_GENBTagJet_4bin->Sumw2();
   pt_GENBTagJet_5bin  = fs->make<TH1F>("pt_GENBTagJet_5bin","pt_GENBTagJet_5bin",Ptbins,Ptbinning); pt_GENBTagJet_5bin->Sumw2();
   pt_GENBTagJet_6bin  = fs->make<TH1F>("pt_GENBTagJet_6bin","pt_GENBTagJet_6bin",Ptbins,Ptbinning); pt_GENBTagJet_6bin->Sumw2();
   pt_GENBTagJet_7bin  = fs->make<TH1F>("pt_GENBTagJet_7bin","pt_GENBTagJet_7bin",Ptbins,Ptbinning); pt_GENBTagJet_7bin->Sumw2();

   pt_GENInclJetCrossSectNorm_1bin  = fs->make<TH1F>("pt_GENInclJetCrossSectNorm_1bin","pt_GENInclJetCrossSectNorm_1bin",Ptbins,Ptbinning); pt_GENInclJetCrossSectNorm_1bin->Sumw2();
   pt_GENInclJetCrossSectNorm_2bin  = fs->make<TH1F>("pt_GENInclJetCrossSectNorm_2bin","pt_GENInclJetCrossSectNorm_2bin",Ptbins,Ptbinning); pt_GENInclJetCrossSectNorm_2bin->Sumw2();
   pt_GENInclJetCrossSectNorm_3bin  = fs->make<TH1F>("pt_GENInclJetCrossSectNorm_3bin","pt_GENInclJetCrossSectNorm_3bin",Ptbins,Ptbinning); pt_GENInclJetCrossSectNorm_3bin->Sumw2();
   pt_GENInclJetCrossSectNorm_4bin  = fs->make<TH1F>("pt_GENInclJetCrossSectNorm_4bin","pt_GENInclJetCrossSectNorm_4bin",Ptbins,Ptbinning); pt_GENInclJetCrossSectNorm_4bin->Sumw2();
   pt_GENInclJetCrossSectNorm_5bin  = fs->make<TH1F>("pt_GENInclJetCrossSectNorm_5bin","pt_GENInclJetCrossSectNorm_5bin",Ptbins,Ptbinning); pt_GENInclJetCrossSectNorm_5bin->Sumw2();
   pt_GENInclJetCrossSectNorm_6bin  = fs->make<TH1F>("pt_GENInclJetCrossSectNorm_6bin","pt_GENInclJetCrossSectNorm_6bin",Ptbins,Ptbinning); pt_GENInclJetCrossSectNorm_6bin->Sumw2();
   pt_GENInclJetCrossSectNorm_7bin  = fs->make<TH1F>("pt_GENInclJetCrossSectNorm_7bin","pt_GENInclJetCrossSectNorm_7bin",Ptbins,Ptbinning); pt_GENInclJetCrossSectNorm_7bin->Sumw2();
   
   PF_FakeInclusiveJets_1bin  = fs->make<TH1F>("PF_FakeInclusiveJets_1bin","PF_FakeInclusiveJets_1bin",Ptbins,Ptbinning); PF_FakeInclusiveJets_1bin->Sumw2();
   PF_FakeInclusiveJets_2bin  = fs->make<TH1F>("PF_FakeInclusiveJets_2bin","PF_FakeInclusiveJets_2bin",Ptbins,Ptbinning); PF_FakeInclusiveJets_2bin->Sumw2();
   PF_FakeInclusiveJets_3bin  = fs->make<TH1F>("PF_FakeInclusiveJets_3bin","PF_FakeInclusiveJets_3bin",Ptbins,Ptbinning); PF_FakeInclusiveJets_3bin->Sumw2();
   PF_FakeInclusiveJets_4bin  = fs->make<TH1F>("PF_FakeInclusiveJets_4bin","PF_FakeInclusiveJets_4bin",Ptbins,Ptbinning); PF_FakeInclusiveJets_4bin->Sumw2();
   PF_FakeInclusiveJets_5bin  = fs->make<TH1F>("PF_FakeInclusiveJets_5bin","PF_FakeInclusiveJets_5bin",Ptbins,Ptbinning); PF_FakeInclusiveJets_5bin->Sumw2();
   PF_FakeInclusiveJets_6bin  = fs->make<TH1F>("PF_FakeInclusiveJets_6bin","PF_FakeInclusiveJets_6bin",Ptbins,Ptbinning); PF_FakeInclusiveJets_6bin->Sumw2();
   PF_FakeInclusiveJets_7bin  = fs->make<TH1F>("PF_FakeInclusiveJets_7bin","PF_FakeInclusiveJets_7bin",Ptbins,Ptbinning); PF_FakeInclusiveJets_7bin->Sumw2();
   
   Gen_MissInclusiveJets_1bin  = fs->make<TH1F>("Gen_MissInclusiveJets_1bin","Gen_MissInclusiveJets_1bin",Ptbins,Ptbinning); Gen_MissInclusiveJets_1bin->Sumw2();
   Gen_MissInclusiveJets_2bin  = fs->make<TH1F>("Gen_MissInclusiveJets_2bin","Gen_MissInclusiveJets_2bin",Ptbins,Ptbinning); Gen_MissInclusiveJets_2bin->Sumw2();
   Gen_MissInclusiveJets_3bin  = fs->make<TH1F>("Gen_MissInclusiveJets_3bin","Gen_MissInclusiveJets_3bin",Ptbins,Ptbinning); Gen_MissInclusiveJets_3bin->Sumw2();
   Gen_MissInclusiveJets_4bin  = fs->make<TH1F>("Gen_MissInclusiveJets_4bin","Gen_MissInclusiveJets_4bin",Ptbins,Ptbinning); Gen_MissInclusiveJets_4bin->Sumw2();
   Gen_MissInclusiveJets_5bin  = fs->make<TH1F>("Gen_MissInclusiveJets_5bin","Gen_MissInclusiveJets_5bin",Ptbins,Ptbinning); Gen_MissInclusiveJets_5bin->Sumw2();
   Gen_MissInclusiveJets_6bin  = fs->make<TH1F>("Gen_MissInclusiveJets_6bin","Gen_MissInclusiveJets_6bin",Ptbins,Ptbinning); Gen_MissInclusiveJets_6bin->Sumw2();
   Gen_MissInclusiveJets_7bin  = fs->make<TH1F>("Gen_MissInclusiveJets_7bin","Gen_MissInclusiveJets_7bin",Ptbins,Ptbinning); Gen_MissInclusiveJets_7bin->Sumw2();
   
   PF_MatchedInclusiveJets_1bin  = fs->make<TH1F>("PF_MatchedInclusiveJets_1bin","PF_MatchedInclusiveJets_1bin",Ptbins,Ptbinning); PF_MatchedInclusiveJets_1bin->Sumw2();
   PF_MatchedInclusiveJets_2bin  = fs->make<TH1F>("PF_MatchedInclusiveJets_2bin","PF_MatchedInclusiveJets_2bin",Ptbins,Ptbinning); PF_MatchedInclusiveJets_2bin->Sumw2();
   PF_MatchedInclusiveJets_3bin  = fs->make<TH1F>("PF_MatchedInclusiveJets_3bin","PF_MatchedInclusiveJets_3bin",Ptbins,Ptbinning); PF_MatchedInclusiveJets_3bin->Sumw2();
   PF_MatchedInclusiveJets_4bin  = fs->make<TH1F>("PF_MatchedInclusiveJets_4bin","PF_MatchedInclusiveJets_4bin",Ptbins,Ptbinning); PF_MatchedInclusiveJets_4bin->Sumw2();
   PF_MatchedInclusiveJets_5bin  = fs->make<TH1F>("PF_MatchedInclusiveJets_5bin","PF_MatchedInclusiveJets_5bin",Ptbins,Ptbinning); PF_MatchedInclusiveJets_5bin->Sumw2();
   PF_MatchedInclusiveJets_6bin  = fs->make<TH1F>("PF_MatchedInclusiveJets_6bin","PF_MatchedInclusiveJets_6bin",Ptbins,Ptbinning); PF_MatchedInclusiveJets_6bin->Sumw2();
   PF_MatchedInclusiveJets_7bin  = fs->make<TH1F>("PF_MatchedInclusiveJets_7bin","PF_MatchedInclusiveJets_7bin",Ptbins,Ptbinning); PF_MatchedInclusiveJets_7bin->Sumw2();
   
   Gen_MatchedInclusiveJets_1bin  = fs->make<TH1F>("Gen_MatchedInclusiveJets_1bin","Gen_MatchedInclusiveJets_1bin",Ptbins,Ptbinning); Gen_MatchedInclusiveJets_1bin->Sumw2();
   Gen_MatchedInclusiveJets_2bin  = fs->make<TH1F>("Gen_MatchedInclusiveJets_2bin","Gen_MatchedInclusiveJets_2bin",Ptbins,Ptbinning); Gen_MatchedInclusiveJets_2bin->Sumw2();
   Gen_MatchedInclusiveJets_3bin  = fs->make<TH1F>("Gen_MatchedInclusiveJets_3bin","Gen_MatchedInclusiveJets_3bin",Ptbins,Ptbinning); Gen_MatchedInclusiveJets_3bin->Sumw2();
   Gen_MatchedInclusiveJets_4bin  = fs->make<TH1F>("Gen_MatchedInclusiveJets_4bin","Gen_MatchedInclusiveJets_4bin",Ptbins,Ptbinning); Gen_MatchedInclusiveJets_4bin->Sumw2();
   Gen_MatchedInclusiveJets_5bin  = fs->make<TH1F>("Gen_MatchedInclusiveJets_5bin","Gen_MatchedInclusiveJets_5bin",Ptbins,Ptbinning); Gen_MatchedInclusiveJets_5bin->Sumw2();
   Gen_MatchedInclusiveJets_6bin  = fs->make<TH1F>("Gen_MatchedInclusiveJets_6bin","Gen_MatchedInclusiveJets_6bin",Ptbins,Ptbinning); Gen_MatchedInclusiveJets_6bin->Sumw2();
   Gen_MatchedInclusiveJets_7bin  = fs->make<TH1F>("Gen_MatchedInclusiveJets_7bin","Gen_MatchedInclusiveJets_7bin",Ptbins,Ptbinning); Gen_MatchedInclusiveJets_7bin->Sumw2();
   
   TwoD_MatchedInclusiveJets_1bin   = fs->make<TH2F>("TwoD_MatchedInclusiveJets_1bin","TwoD_MatchedInclusiveJets_1bin",Ptbins,Ptbinning,Ptbins,Ptbinning); TwoD_MatchedInclusiveJets_1bin->Sumw2(); 
   TwoD_MatchedInclusiveJets_2bin   = fs->make<TH2F>("TwoD_MatchedInclusiveJets_2bin","TwoD_MatchedInclusiveJets_2bin",Ptbins,Ptbinning,Ptbins,Ptbinning); TwoD_MatchedInclusiveJets_2bin->Sumw2(); 
   TwoD_MatchedInclusiveJets_3bin   = fs->make<TH2F>("TwoD_MatchedInclusiveJets_3bin","TwoD_MatchedInclusiveJets_3bin",Ptbins,Ptbinning,Ptbins,Ptbinning); TwoD_MatchedInclusiveJets_3bin->Sumw2(); 
   TwoD_MatchedInclusiveJets_4bin   = fs->make<TH2F>("TwoD_MatchedInclusiveJets_4bin","TwoD_MatchedInclusiveJets_4bin",Ptbins,Ptbinning,Ptbins,Ptbinning); TwoD_MatchedInclusiveJets_4bin->Sumw2(); 
   TwoD_MatchedInclusiveJets_5bin   = fs->make<TH2F>("TwoD_MatchedInclusiveJets_5bin","TwoD_MatchedInclusiveJets_5bin",Ptbins,Ptbinning,Ptbins,Ptbinning); TwoD_MatchedInclusiveJets_5bin->Sumw2(); 
   TwoD_MatchedInclusiveJets_6bin   = fs->make<TH2F>("TwoD_MatchedInclusiveJets_6bin","TwoD_MatchedInclusiveJets_6bin",Ptbins,Ptbinning,Ptbins,Ptbinning); TwoD_MatchedInclusiveJets_6bin->Sumw2(); 
   TwoD_MatchedInclusiveJets_7bin   = fs->make<TH2F>("TwoD_MatchedInclusiveJets_7bin","TwoD_MatchedInclusiveJets_7bin",Ptbins,Ptbinning,Ptbins,Ptbinning); TwoD_MatchedInclusiveJets_7bin->Sumw2(); 
   
   Resolution1D = fs->make<TH1F>("Resolution1D","Resolution1D",100,-5,5);
   ResolutionForward1D = fs->make<TH1F>("ResolutionForward1D","ResolutionForward1D",100,-5,5);
   
   ResolutionInclusiveJets = fs->make<TProfile>("ResolutionInclusiveJets","ResolutionInclusiveJets",Ptbins,Ptbinning,0,5);
   ResolutionInclusiveJets_1bin = fs->make<TProfile>("ResolutionInclusiveJets_1bin","ResolutionInclusiveJets_1bin",Ptbins,Ptbinning,0,5);
   ResolutionInclusiveJets_2bin = fs->make<TProfile>("ResolutionInclusiveJets_2bin","ResolutionInclusiveJets_2bin",Ptbins,Ptbinning,0,5);
   ResolutionInclusiveJets_3bin = fs->make<TProfile>("ResolutionInclusiveJets_3bin","ResolutionInclusiveJets_3bin",Ptbins,Ptbinning,0,5);
   ResolutionInclusiveJets_4bin = fs->make<TProfile>("ResolutionInclusiveJets_4bin","ResolutionInclusiveJets_4bin",Ptbins,Ptbinning,0,5);
   ResolutionInclusiveJets_5bin = fs->make<TProfile>("ResolutionInclusiveJets_5bin","ResolutionInclusiveJets_5bin",Ptbins,Ptbinning,0,5);
   ResolutionInclusiveJets_6bin = fs->make<TProfile>("ResolutionInclusiveJets_6bin","ResolutionInclusiveJets_6bin",Ptbins,Ptbinning,0,5);
   ResolutionInclusiveJets_7bin = fs->make<TProfile>("ResolutionInclusiveJets_7bin","ResolutionInclusiveJets_7bin",Ptbins,Ptbinning,0,5);

   ResolutionTagAndProbe = fs->make<TProfile>("ResolutionTagAndProbe","ResolutionTagAndProbe",Ptbins,Ptbinning,0,5);
   
   //trigger efficiency measurement
   TagAndProbeNum= fs->make<TH1F>("TagAndProbeNum","TagAndProbeNum",Ptbins,Ptbinning); 
   TagAndProbeDen= fs->make<TH1F>("TagAndProbeDen","TagAndProbeDen",Ptbins,Ptbinning); 
   TagAndProbeEff= fs->make<TH1F>("TagAndProbeEff","TagAndProbeEff",Ptbins,Ptbinning); 

   hist_leading_pt_emulated_Jet60= fs->make<TH1F>("hist_leading_pt_emulated_Jet60","hist_leading_pt_emulated_Jet60",Ptbins,Ptbinning); 
   hist_leading_pt_all_Jet60= fs->make<TH1F>("hist_leading_pt_all_Jet60","hist_leading_pt_all_Jet60",Ptbins,Ptbinning); 
   hist_leading_eta_emulated_Jet60= fs->make<TH1F>("hist_leading_eta_emulated_Jet60","hist_leading_eta_emulated_Jet60",24,-5.2,5.2); 
   hist_leading_eta_all_Jet60= fs->make<TH1F>("hist_leading_eta_all_Jet60","hist_leading_eta_all_Jet60",24,-5.2,5.2); 
   
   hist_leading_pt_HLT_Jet60U_eff= fs->make<TH1F>("hist_leading_pt_HLT_Jet60U_eff","hist_leading_pt_HLT_Jet60U_eff",Ptbins,Ptbinning); 
   hist_leading_eta_HLT_Jet60U_eff= fs->make<TH1F>("hist_leading_eta_HLT_Jet60U_eff","hist_leading_eta_HLT_Jet60U_eff",24,-5.2,5.2); 
   
   hist_leading_pt_emulated_Jet80= fs->make<TH1F>("hist_leading_pt_emulated_Jet80","hist_leading_pt_emulated_Jet80",Ptbins,Ptbinning);
   hist_leading_pt_all_Jet80= fs->make<TH1F>("hist_leading_pt_all_Jet80","hist_leading_pt_all_Jet80",Ptbins,Ptbinning);
   hist_leading_eta_emulated_Jet80= fs->make<TH1F>("hist_leading_eta_emulated_Jet80","hist_leading_eta_emulated_Jet80",24,-5.2,5.2);
   hist_leading_eta_all_Jet80= fs->make<TH1F>("hist_leading_eta_all_Jet80","hist_leading_eta_all_Jet80",24,-5.2,5.2);
   
   hist_leading_pt_HLT_Jet80U_eff= fs->make<TH1F>("hist_leading_pt_HLT_Jet80U_eff","hist_leading_pt_HLT_Jet80U_eff",Ptbins,Ptbinning);
   hist_leading_eta_HLT_Jet80U_eff= fs->make<TH1F>("hist_leading_eta_HLT_Jet80U_eff","hist_leading_eta_HLT_Jet80U_eff",24,-5.2,5.2);
   
   hist_leading_pt_emulated_Jet140= fs->make<TH1F>("hist_leading_pt_emulated_Jet140","hist_leading_pt_emulated_Jet140",Ptbins,Ptbinning);
   hist_leading_pt_all_Jet140= fs->make<TH1F>("hist_leading_pt_all_Jet140","hist_leading_pt_all_Jet140",Ptbins,Ptbinning);
   hist_leading_eta_emulated_Jet140= fs->make<TH1F>("hist_leading_eta_emulated_Jet140","hist_leading_eta_emulated_Jet140",24,-5.2,5.2);
   hist_leading_eta_all_Jet140= fs->make<TH1F>("hist_leading_eta_all_Jet140","hist_leading_eta_all_Jet140",24,-5.2,5.2);
   
   hist_leading_pt_HLT_Jet140U_eff= fs->make<TH1F>("hist_leading_pt_HLT_Jet140U_eff","hist_leading_pt_HLT_Jet140U_eff",Ptbins,Ptbinning);
   hist_leading_eta_HLT_Jet140U_eff= fs->make<TH1F>("hist_leading_eta_HLT_Jet140U_eff","hist_leading_eta_HLT_Jet140U_eff",24,-5.2,5.2);
   
   hist_leading_pt_emulated_Jet200= fs->make<TH1F>("hist_leading_pt_emulated_Jet200","hist_leading_pt_emulated_Jet200",Ptbins,Ptbinning);
   hist_leading_pt_all_Jet200= fs->make<TH1F>("hist_leading_pt_all_Jet200","hist_leading_pt_all_Jet200",Ptbins,Ptbinning);
   hist_leading_eta_emulated_Jet200= fs->make<TH1F>("hist_leading_eta_emulated_Jet200","hist_leading_eta_emulated_Jet200",24,-5.2,5.2);
   hist_leading_eta_all_Jet200= fs->make<TH1F>("hist_leading_eta_all_Jet200","hist_leading_eta_all_Jet200",24,-5.2,5.2);
   
   hist_leading_pt_HLT_Jet200U_eff= fs->make<TH1F>("hist_leading_pt_HLT_Jet200U_eff","hist_leading_pt_HLT_Jet200U_eff",Ptbins,Ptbinning);
   hist_leading_eta_HLT_Jet200U_eff= fs->make<TH1F>("hist_leading_eta_HLT_Jet200U_eff","hist_leading_eta_HLT_Jet200U_eff",24,-5.2,5.2);

   hist_leading_pt_emulated_Jet260= fs->make<TH1F>("hist_leading_pt_emulated_Jet260","hist_leading_pt_emulated_Jet260",Ptbins,Ptbinning);
   hist_leading_pt_all_Jet260= fs->make<TH1F>("hist_leading_pt_all_Jet260","hist_leading_pt_all_Jet260",Ptbins,Ptbinning);
   hist_leading_eta_emulated_Jet260= fs->make<TH1F>("hist_leading_eta_emulated_Jet260","hist_leading_eta_emulated_Jet260",24,-5.2,5.2);
   hist_leading_eta_all_Jet260= fs->make<TH1F>("hist_leading_eta_all_Jet260","hist_leading_eta_all_Jet260",24,-5.2,5.2);
   
   hist_leading_pt_HLT_Jet260U_eff= fs->make<TH1F>("hist_leading_pt_HLT_Jet260U_eff","hist_leading_pt_HLT_Jet260U_eff",Ptbins,Ptbinning);
   hist_leading_eta_HLT_Jet260U_eff= fs->make<TH1F>("hist_leading_eta_HLT_Jet260U_eff","hist_leading_eta_HLT_Jet260U_eff",24,-5.2,5.2);
   
   hist_leading_pt_emulated_Jet320= fs->make<TH1F>("hist_leading_pt_emulated_Jet320","hist_leading_pt_emulated_Jet320",Ptbins,Ptbinning);
   hist_leading_pt_all_Jet320= fs->make<TH1F>("hist_leading_pt_all_Jet320","hist_leading_pt_all_Jet320",Ptbins,Ptbinning);
   hist_leading_eta_emulated_Jet320= fs->make<TH1F>("hist_leading_eta_emulated_Jet320","hist_leading_eta_emulated_Jet320",24,-5.2,5.2);
   hist_leading_eta_all_Jet320= fs->make<TH1F>("hist_leading_eta_all_Jet320","hist_leading_eta_all_Jet320",24,-5.2,5.2);
   
   hist_leading_pt_HLT_Jet320U_eff= fs->make<TH1F>("hist_leading_pt_HLT_Jet320U_eff","hist_leading_pt_HLT_Jet320U_eff",Ptbins,Ptbinning);
   hist_leading_eta_HLT_Jet320U_eff= fs->make<TH1F>("hist_leading_eta_HLT_Jet320U_eff","hist_leading_eta_HLT_Jet320U_eff",24,-5.2,5.2);
   
   hist_leading_pt_emulated_Jet400= fs->make<TH1F>("hist_leading_pt_emulated_Jet400","hist_leading_pt_emulated_Jet400",Ptbins,Ptbinning);
   hist_leading_pt_all_Jet400= fs->make<TH1F>("hist_leading_pt_all_Jet400","hist_leading_pt_all_Jet400",Ptbins,Ptbinning);
   hist_leading_eta_emulated_Jet400= fs->make<TH1F>("hist_leading_eta_emulated_Jet400","hist_leading_eta_emulated_Jet400",24,-5.2,5.2);
   hist_leading_eta_all_Jet400= fs->make<TH1F>("hist_leading_eta_all_Jet400","hist_leading_eta_all_Jet400",24,-5.2,5.2);
   
   hist_leading_pt_HLT_Jet400U_eff= fs->make<TH1F>("hist_leading_pt_HLT_Jet400U_eff","hist_leading_pt_HLT_Jet400U_eff",Ptbins,Ptbinning);
   hist_leading_eta_HLT_Jet400U_eff= fs->make<TH1F>("hist_leading_eta_HLT_Jet400U_eff","hist_leading_eta_HLT_Jet400U_eff",24,-5.2,5.2);
   
   hist_leading_pt_emulated_Jet450= fs->make<TH1F>("hist_leading_pt_emulated_Jet450","hist_leading_pt_emulated_Jet450",Ptbins,Ptbinning);
   hist_leading_pt_all_Jet450= fs->make<TH1F>("hist_leading_pt_all_Jet450","hist_leading_pt_all_Jet450",Ptbins,Ptbinning);
   hist_leading_eta_emulated_Jet450= fs->make<TH1F>("hist_leading_eta_emulated_Jet450","hist_leading_eta_emulated_Jet450",24,-5.2,5.2);
   hist_leading_eta_all_Jet450= fs->make<TH1F>("hist_leading_eta_all_Jet450","hist_leading_eta_all_Jet450",24,-5.2,5.2);
   
   hist_leading_pt_HLT_Jet450U_eff= fs->make<TH1F>("hist_leading_pt_HLT_Jet450U_eff","hist_leading_pt_HLT_Jet450U_eff",Ptbins,Ptbinning);
   hist_leading_eta_HLT_Jet450U_eff= fs->make<TH1F>("hist_leading_eta_HLT_Jet450U_eff","hist_leading_eta_HLT_Jet450U_eff",24,-5.2,5.2);
   
   hist_leading_pt_emulated_Jet500= fs->make<TH1F>("hist_leading_pt_emulated_Jet500","hist_leading_pt_emulated_Jet500",Ptbins,Ptbinning);
   hist_leading_pt_all_Jet500= fs->make<TH1F>("hist_leading_pt_all_Jet500","hist_leading_pt_all_Jet500",Ptbins,Ptbinning);
   hist_leading_eta_emulated_Jet500= fs->make<TH1F>("hist_leading_eta_emulated_Jet500","hist_leading_eta_emulated_Jet500",24,-5.2,5.2);
   hist_leading_eta_all_Jet500= fs->make<TH1F>("hist_leading_eta_all_Jet500","hist_leading_eta_all_Jet500",24,-5.2,5.2);
   
   hist_leading_pt_HLT_Jet500U_eff= fs->make<TH1F>("hist_leading_pt_HLT_Jet500U_eff","hist_leading_pt_HLT_Jet500U_eff",Ptbins,Ptbinning);
   hist_leading_eta_HLT_Jet500U_eff= fs->make<TH1F>("hist_leading_eta_HLT_Jet500U_eff","hist_leading_eta_HLT_Jet500U_eff",24,-5.2,5.2);
   
   //Jet energy systematic uncertainty
   pt_DETInclJetUP_1bin  = fs->make<TH1F>("pt_DETInclJetUP_1bin","pt_DETInclJetUP_1bin",Ptbins,Ptbinning); pt_DETInclJetUP_1bin->Sumw2();
   pt_DETInclJetUP_2bin  = fs->make<TH1F>("pt_DETInclJetUP_2bin","pt_DETInclJetUP_2bin",Ptbins,Ptbinning); pt_DETInclJetUP_2bin->Sumw2();
   pt_DETInclJetUP_3bin  = fs->make<TH1F>("pt_DETInclJetUP_3bin","pt_DETInclJetUP_3bin",Ptbins,Ptbinning); pt_DETInclJetUP_3bin->Sumw2();
   pt_DETInclJetUP_4bin  = fs->make<TH1F>("pt_DETInclJetUP_4bin","pt_DETInclJetUP_4bin",Ptbins,Ptbinning); pt_DETInclJetUP_4bin->Sumw2();
   pt_DETInclJetUP_5bin  = fs->make<TH1F>("pt_DETInclJetUP_5bin","pt_DETInclJetUP_5bin",Ptbins,Ptbinning); pt_DETInclJetUP_5bin->Sumw2();
   pt_DETInclJetUP_6bin  = fs->make<TH1F>("pt_DETInclJetUP_6bin","pt_DETInclJetUP_6bin",Ptbins,Ptbinning); pt_DETInclJetUP_6bin->Sumw2();
   pt_DETInclJetUP_7bin  = fs->make<TH1F>("pt_DETInclJetUP_7bin","pt_DETInclJetUP_7bin",Ptbins,Ptbinning); pt_DETInclJetUP_7bin->Sumw2();
   
   pt_DETInclJetDOWN_1bin  = fs->make<TH1F>("pt_DETInclJetDOWN_1bin","pt_DETInclJetDOWN_1bin",Ptbins,Ptbinning); pt_DETInclJetDOWN_1bin->Sumw2();
   pt_DETInclJetDOWN_2bin  = fs->make<TH1F>("pt_DETInclJetDOWN_2bin","pt_DETInclJetDOWN_2bin",Ptbins,Ptbinning); pt_DETInclJetDOWN_2bin->Sumw2();
   pt_DETInclJetDOWN_3bin  = fs->make<TH1F>("pt_DETInclJetDOWN_3bin","pt_DETInclJetDOWN_3bin",Ptbins,Ptbinning); pt_DETInclJetDOWN_3bin->Sumw2();
   pt_DETInclJetDOWN_4bin  = fs->make<TH1F>("pt_DETInclJetDOWN_4bin","pt_DETInclJetDOWN_4bin",Ptbins,Ptbinning); pt_DETInclJetDOWN_4bin->Sumw2();
   pt_DETInclJetDOWN_5bin  = fs->make<TH1F>("pt_DETInclJetDOWN_5bin","pt_DETInclJetDOWN_5bin",Ptbins,Ptbinning); pt_DETInclJetDOWN_5bin->Sumw2();
   pt_DETInclJetDOWN_6bin  = fs->make<TH1F>("pt_DETInclJetDOWN_6bin","pt_DETInclJetDOWN_6bin",Ptbins,Ptbinning); pt_DETInclJetDOWN_6bin->Sumw2();
   pt_DETInclJetDOWN_7bin  = fs->make<TH1F>("pt_DETInclJetDOWN_7bin","pt_DETInclJetDOWN_7bin",Ptbins,Ptbinning); pt_DETInclJetDOWN_7bin->Sumw2();

   Chargedhf0_DETJet = fs->make<TH1F>("Chargedhf0_DETJet","Chargedhf0_DETJet",50,0,1); Chargedhf0_DETJet->Sumw2();
   Chargedef0_DETJet = fs->make<TH1F>("Chargedef0_DETJet","Chargedef0_DETJet",50,0,1); Chargedef0_DETJet->Sumw2();
   Neutralhf0_DETJet = fs->make<TH1F>("Neutralhf0_DETJet","Neutralhf0_DETJet",50,0,1); Neutralhf0_DETJet->Sumw2();
   Photonef0_DETJet = fs->make<TH1F>("Photonef0_DETJet","Photonef0_DETJet",50,0,1); Photonef0_DETJet->Sumw2();
   Hadronef0_DETJet = fs->make<TH1F>("Hadronef0_DETJet","Hadronef0_DETJet",50,0,1); Hadronef0_DETJet->Sumw2();
   Muonef0_DETJet = fs->make<TH1F>("Muonef0_DETJet","Muonef0_DETJet",50,0,1); Muonef0_DETJet->Sumw2();
   Electromagneticef0_DETJet = fs->make<TH1F>("Electromagneticef0_DETJet","Electromagneticef0_DETJet",50,0,1); Electromagneticef0_DETJet->Sumw2();
   
   ChargedhMultiplicity0_DETJet = fs->make<TH1F>("ChargedhMultiplicity0_DETJet","ChargedhMultiplicity0_DETJet",100,0,100); ChargedhMultiplicity0_DETJet->Sumw2();
   ChargedeMultiplicity0_DETJet = fs->make<TH1F>("ChargedeMultiplicity0_DETJet","ChargedeMultiplicity0_DETJet",100,0,100); ChargedeMultiplicity0_DETJet->Sumw2();
   NeutralhMultiplicity0_DETJet = fs->make<TH1F>("NeutralhMultiplicity0_DETJet","NeutralhMultiplicity0_DETJet",100,0,100); NeutralhMultiplicity0_DETJet->Sumw2();
   PhotoneMultiplicity0_DETJet = fs->make<TH1F>("PhotoneMultiplicity0_DETJet","PhotoneMultiplicity0_DETJet",100,0,100); PhotoneMultiplicity0_DETJet->Sumw2();
   HadroneMultiplicity0_DETJet = fs->make<TH1F>("HadroneMultiplicity0_DETJet","HadroneMultiplicity0_DETJet",100,0,100); HadroneMultiplicity0_DETJet->Sumw2();
   MuoneMultiplicity0_DETJet = fs->make<TH1F>("MuoneMultiplicity0_DETJet","MuoneMultiplicity0_DETJet",100,0,100); MuoneMultiplicity0_DETJet->Sumw2();
   ElectromagneticeMultiplicity0_DETJet = fs->make<TH1F>("ElectromagneticeMultiplicity0_DETJet","ElectromagneticeMultiplicity0_DETJet",100,0,100); ElectromagneticeMultiplicity0_DETJet->Sumw2();

   TruePileUpMC = fs->make<TH1F>("TruePileUpMC","TruePileUpMC",100,0,100); TruePileUpMC->Sumw2();
   TruePileUpMCInteger = fs->make<TH1F>("TruePileUpMCInteger","TruePileUpMCInteger",100,0,100); TruePileUpMCInteger->Sumw2();

   TruePileUpDataInteger = fs->make<TH1F>("TruePileUpDataInteger","TruePileUpDataInteger",100,0,100); TruePileUpDataInteger->Sumw2();

   MET_DET  = fs->make<TH1F>("MET_DET","MET_DET",500,0,2000); MET_DET->Sumw2();
   METPhi_DET  = fs->make<TH1F>("METPhi_DET","METPhi_DET",35,-7,7); METPhi_DET->Sumw2();
   FractionMET_DET  = fs->make<TH1F>("FractionMET_DET","FractionMET_DET",50,0,1); FractionMET_DET->Sumw2();
   FractionMET_DET_1bin  = fs->make<TH1F>("FractionMET_DET_1bin","FractionMET_DET_1bin",50,0,1); FractionMET_DET_1bin->Sumw2();
   FractionMET_DET_2bin  = fs->make<TH1F>("FractionMET_DET_2bin","FractionMET_DET_2bin",50,0,1); FractionMET_DET_2bin->Sumw2();
   FractionMET_DET_3bin  = fs->make<TH1F>("FractionMET_DET_3bin","FractionMET_DET_3bin",50,0,1); FractionMET_DET_3bin->Sumw2();
   FractionMET_DET_4bin  = fs->make<TH1F>("FractionMET_DET_4bin","FractionMET_DET_4bin",50,0,1); FractionMET_DET_4bin->Sumw2();
   FractionMET_DET_5bin  = fs->make<TH1F>("FractionMET_DET_5bin","FractionMET_DET_5bin",50,0,1); FractionMET_DET_5bin->Sumw2();
   FractionMET_DET_6bin  = fs->make<TH1F>("FractionMET_DET_6bin","FractionMET_DET_6bin",50,0,1); FractionMET_DET_6bin->Sumw2();
   FractionMET_DET_7bin  = fs->make<TH1F>("FractionMET_DET_7bin","FractionMET_DET_7bin",50,0,1); FractionMET_DET_7bin->Sumw2();

   //B-tag Discriminants

   CSV_DET  = fs->make<TH1F>("CSV_DET","CSV_DET",50,0,1); CSV_DET->Sumw2();
   TCHP_DET  = fs->make<TH1F>("TCHP_DET","TCHP_DET",100,0,10); TCHP_DET->Sumw2();
   TCHE_DET  = fs->make<TH1F>("TCHE_DET","TCHE_DET",100,0,10); TCHE_DET->Sumw2();
   SimpleSecVertexHE_DET  = fs->make<TH1F>("SimpleSecVertexHE_DET","SimpleSecVertexHE_DET",50,0,10); SimpleSecVertexHE_DET->Sumw2();
   SimpleSecVertexHP_DET  = fs->make<TH1F>("SimpleSecVertexHP_DET","SimpleSecVertexHP_DET",50,0,10); SimpleSecVertexHP_DET->Sumw2();

   jecs = new JECs(mIsMCarlo, mGlobalTag, mjettype,mJECUncSrc,mJECUncSrcNames);

   //Unfolding
   /*resp_jetpt1etabin = fs->make<RooUnfoldResponse>(pt0_DETJet,pt0_DETJet, "resp_jetpt1etabin", "jetpt1etabin");
     resp_jetpt2etabin = fs->make<RooUnfoldResponse>(pt0_DETJet,pt0_DETJet, "resp_jetpt2etabin", "jetpt2etabin");
     resp_jetpt3etabin = fs->make<RooUnfoldResponse>(pt0_DETJet,pt0_DETJet, "resp_jetpt3etabin", "jetpt3etabin");
     resp_jetpt4etabin = fs->make<RooUnfoldResponse>(pt0_DETJet,pt0_DETJet, "resp_jetpt4etabin", "jetpt4etabin");
     resp_jetpt5etabin = fs->make<RooUnfoldResponse>(pt0_DETJet,pt0_DETJet, "resp_jetpt5etabin", "jetpt5etabin");
     resp_jetpt6etabin = fs->make<RooUnfoldResponse>(pt0_DETJet,pt0_DETJet, "resp_jetpt6etabin", "jetpt6etabin");
     resp_jetpt7etabin = fs->make<RooUnfoldResponse>(pt0_DETJet,pt0_DETJet, "resp_jetpt7etabin", "jetpt7etabin");*/

 } // end of function beginJob()


 //------------------------ endjob() function declaration ---------------------- //
 void Analysis_Template_MC::endJob()
 {

   mInf->Close();
   
 } // closing endJob()





 //--------------------------- analyze() fuction declaration ------------------ //
void Analysis_Template_MC::analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup)
 {

   cout<<" Size "<<mFileName.size()<<endl;
   
   for(unsigned ifile=0;ifile<mFileName.size();ifile++){
     
     mInf = TFile::Open(mFileName[ifile].c_str());
     mDir = (TDirectoryFile*)mInf->Get(mDirName.c_str());  
     mTree=(TTree*)mDir->Get(mTreeName.c_str());

     cout << mFileName[ifile].c_str() << endl;
     
     Event = new QCDEvent();
     TBranch *branch = mTree->GetBranch("events");
     branch->SetAddress(&Event);  

     unsigned NEntries = mTree->GetEntries();
     cout<<"Reading TREE: "<<NEntries<<" events"<<endl;  
     
     int decade = 0 ;
     
     float hweight=1.;  ///Initial value to one
     
     int nJetTrigList = 20;
     int nJetTrig = 20;
  
     TString HLTJet[nJetTrig];
     int ihltj[nJetTrig];
     float prescalej[nJetTrig];
     
     int HLTJetPtN[10] = {40,60,80,140,200,260,320,400,450,500};
     int HLTJetPtNThres[10] = {60,80,140,200,260,320,400,450,500,600};
     int ATrig[10] = {10,10,10,10,10,10,10,10,10,10};
     
     bool l1cut[nJetTrig], hltcut[nJetTrig];
     bool hltPassj[nJetTrig];
     char trigtitle[200];
     
     for (int i=0; i<nJetTrigList/2; i++){
       sprintf(trigtitle, "HLT_PFJet%i_v3",HLTJetPtN[i]);
       HLTJet[i] = trigtitle;
       sprintf(trigtitle, "HLT_PFJet%i_v4",HLTJetPtN[i]);
       HLTJet[i+10] = trigtitle;
     }

     for (int i=0; i<nJetTrig; i++){
       ihltj[i] = -1 ;
     }
     
     ///Trigger infos
     TH1F *hTrigNames = (TH1F*)mInf->Get("ak4/TriggerNames");
     cout<<"Finding trigger mapping: "<<endl;
     //----------- loop over the X-axis labels -----------------
     for(int ibin=0;ibin<hTrigNames->GetNbinsX();ibin++) {
       TString ss(hTrigNames->GetXaxis()->GetBinLabel(ibin+1));
       
       for (int ii=0; ii<nJetTrig; ii++){
	 if (ss == HLTJet[ii]) {
	   ihltj[ii] = ibin;
	   continue;
	 }
       } // for (int ii=0; ii<nJetTrig; ii++)                                                                                                                                 
     }
     
     for (int ij=0; ij<nJetTrig; ij++){
       if (ihltj[ij] == -1) {
	 cout<<"The requested trigger ("<<HLTJet[ij]<<") is not found "<<endl;
	 //         break;                                                                                                                                                            
       }
       else {
	 cout<<HLTJet[ij]<<" --> "<<ihltj[ij]<<endl;
       }
     }
     //
     int random=0;

     for(unsigned  l=0; l<NEntries; l++) {
     //for(unsigned  l=0; l<1000; l++) {
       
       double WeightMC=1;
       unsigned NEntriesNorm=NEntries;

       mTree->GetEntry(l);

       //LOW PILE-UP
       
       double TruePileUpInteger=0;
       if(mIsMCarlo) TruePileUpInteger=Event->evtHdr().pu()/Event->evtHdr().nbx();
              
       if(mIsMCarlo){
	 if(mMCSlice==0){WeightMC=2022100000*Event->evtHdr().weight()/NEntriesNorm;}//flat MC
	 if(mMCSlice==1){
	   
  /*       
         // MadGraph Monte Carlo
          if(ifile < 11){WeightMC=27540000*Event->evtHdr().weight()/NEntriesNorm;}//100-200 Slice
          if(ifile==11){WeightMC=1717000*Event->evtHdr().weight()/NEntriesNorm;}//200-300 Slice
          if(ifile==12){WeightMC=351300*Event->evtHdr().weight()/NEntriesNorm;}//300-500
          if(ifile==13){WeightMC=31630*Event->evtHdr().weight()/NEntriesNorm;}//500-700
          if(ifile==14){WeightMC=6802*Event->evtHdr().weight()/NEntriesNorm;}//700-1000
          if(ifile==15){WeightMC=1206*Event->evtHdr().weight()/NEntriesNorm;}//1000-1500
          if(ifile==16){WeightMC=120.4*Event->evtHdr().weight()/NEntriesNorm;}//1500-2000
          if(ifile==17){WeightMC=25.25*Event->evtHdr().weight()/NEntriesNorm;}//2000-Inf
            
    */   
      
        //      P8 Monte Carlo
       if(ifile==0){WeightMC=140932000*Event->evtHdr().weight()/NEntriesNorm;}//30-50 Slice
	   if(ifile==1){WeightMC=19204300*Event->evtHdr().weight()/NEntriesNorm;}//50-80 Slice
	   if(ifile==2){WeightMC=2762530*Event->evtHdr().weight()/NEntriesNorm;}//80-120 Slice
	   if(ifile==3){WeightMC=471100*Event->evtHdr().weight()/NEntriesNorm;}//120-170
	   if(ifile==4){WeightMC=117276*Event->evtHdr().weight()/NEntriesNorm;}//170-300
	   if(ifile==5){WeightMC=7823*Event->evtHdr().weight()/NEntriesNorm;}//300-470
	   if(ifile==6){WeightMC=648.2*Event->evtHdr().weight()/NEntriesNorm;}//470-600
	   if(ifile==7){WeightMC=186.9*Event->evtHdr().weight()/NEntriesNorm;}//600-800
	   if(ifile==8){WeightMC=32.293*Event->evtHdr().weight()/NEntriesNorm;}//800-1000
	   if(ifile==9){WeightMC=9.4183*Event->evtHdr().weight()/NEntriesNorm;}//1000-1400
	   if(ifile==10){WeightMC=0.84265*Event->evtHdr().weight()/NEntriesNorm;}//1400-1800
	   if(ifile==11){WeightMC=0.114943*Event->evtHdr().weight()/NEntriesNorm;}//1800-2400
	   if(ifile==12){WeightMC=0.00682981*Event->evtHdr().weight()/NEntriesNorm;}//2400-3200
	   if(ifile==13){WeightMC=0.000165445*Event->evtHdr().weight()/NEntriesNorm;}//3200-Inf
	    
    }


	 if(mPUReweighting) {
	   double PUValue=Event->evtHdr().pu()/Event->evtHdr().nbx();

	   //RUNC
	   if(PUValue==0) WeightMC=WeightMC*467.57407;
           if(PUValue==1) WeightMC=WeightMC*30.31823;
           if(PUValue==2) WeightMC=WeightMC*26.30839;
           if(PUValue==3) WeightMC=WeightMC*6.08921;
           if(PUValue==4) WeightMC=WeightMC*0;
           if(PUValue==5) WeightMC=WeightMC*0;
           if(PUValue==6) WeightMC=WeightMC*0.97311;
           if(PUValue==7) WeightMC=WeightMC*0.12387;
           if(PUValue==8) WeightMC=WeightMC*0.00199;
           if(PUValue==9) WeightMC=WeightMC*0;
           if(PUValue==10) WeightMC=WeightMC*0.03264;
           if(PUValue==11) WeightMC=WeightMC*0.05545;
           if(PUValue==12) WeightMC=WeightMC*0.01525;
           if(PUValue==13) WeightMC=WeightMC*0.00746;
           if(PUValue==14) WeightMC=WeightMC*0.00383;
           if(PUValue==15) WeightMC=WeightMC*0.43554;
           if(PUValue==16) WeightMC=WeightMC*1.45626;
           if(PUValue==17) WeightMC=WeightMC*1.07330;
           if(PUValue==18) WeightMC=WeightMC*1.20143;
           if(PUValue==19) WeightMC=WeightMC*1;
           if(PUValue==20) WeightMC=WeightMC*0.98406;
           if(PUValue==21) WeightMC=WeightMC*1.51376;
           if(PUValue==22) WeightMC=WeightMC*1.24476;
           if(PUValue==23) WeightMC=WeightMC*1.25060;
           if(PUValue==24) WeightMC=WeightMC*0.97207;
	   if(PUValue==25) WeightMC=WeightMC*0.47589;
	   if(PUValue>=26) WeightMC=WeightMC*0;
	 }
	 
	 hweight=WeightMC;
	 
       }
       if(hweight<0){ continue; }

       //cout<<"************************* NEW EVENT ******************************************"<<endl;
    
       //----------- progress report -------------
       double progress = 10.0*l/(1.0*NEntriesNorm);
       int k = TMath::FloorNint(progress);
       if (k > decade)
	 cout<<10*k<<" %"<<endl;
       decade = k;
       //----------- read the event --------------
       
       if(!mIsMCarlo){
	 
	 for (int j=0; j<nJetTrig; j++){
	   hltPassj[j] = false ;
	   l1cut[j] = false;
	   hltcut[j] = false;
	   prescalej[j] = 1;
	 }
      
	 for (int j=0; j<nJetTrig; j++){
	
	   if (Event->fired(ihltj[j]) > 0) {
	     hltPassj[j] = true;
	     if(Event->preL1(ihltj[j])>=0){ prescalej[j] = Event->preL1(ihltj[j]) * Event->preHLT(ihltj[j]);}
	     else {prescalej[j] = Event->preHLT(ihltj[j]);}
	     if(mprintOk==4){
	       if(j==3 || j==27) cout<<j<<" "<<HLTJet[j]<<" "<<prescalej[j]<<" PASS **********************************"<<Event->preL1(ihltj[j])<<" "<< Event->preHLT(ihltj[j])<<endl;
	     }

	     //cout<<Event->evtHdr().runNo()<<" "<<Event->evtHdr().event()<<" "<<Event->evtHdr().lumi()<<" "<<HLTJet[j]<<" "<<prescaleCorrect<<endl; 

	   }
	 }
       }
    
       ///pthat & mc_weight
       //int NumberEvent=Event->evtHdr().lumi();
    
       double pthat = Event->evtHdr().pthat();
       double mc_weight = Event->evtHdr().weight();
       if(mprintOk==10) printf("\npthat=%f  mc_weight=%e\n",pthat,mc_weight);
       
       mc_pthat->Fill(pthat);
       mc_pthat_weighted->Fill(pthat,hweight);
     
       unsigned n_genJets = Event->nGenJets();
       int DetJets=0;
       int DetJetsBTags=0;
       int DETjet_ok[100]; for(int ii=0;ii<100;++ii){DETjet_ok[ii]=0;}
       
       jecs->JEC_CHScorrections(Event, Event->nPFJetsCHS(), mIsMCarlo,mJECUncSrcNames);
       //jecs->JEC_corrections(Event, Event->nPFJetsCHS(), mIsMCarlo,mJECUncSrcNames);
       
       if(mIsMCarlo){
      
	 /////////////////////////////////////////////////////////////////////////////////////////////
	 ///Examine GenJets
	 
	 /// Dump all GEN jets
	 if(mprintOk==2){
	   printf("Number of GENJets=%d\n",n_genJets);
	   for(unsigned j=0; j<n_genJets; ++j){
	     printf("j=%2d  pt=%8.3f  y=%6.3f  phi=%6.3f\n",j,Event->genjet(j).pt(),Event->genjet(j).Rapidity(),Event->genjet(j).phi());
	   }
	 }
	 
	 ///Apply Jet cuts.  Very General to all existing GEN Jets
	 int GENjet_ok[100]; for(int ii=0;ii<100;++ii){GENjet_ok[ii]=0;}
	 
	 int GenJets=0;
	 
	 if(mprintOk==1) cout<<"Vertex info: numVtx="<<Event->evtHdr().nVtx()<<"  numVtxGood="<<Event->evtHdr().nVtxGood()<<"  isPVgood()="<<Event->evtHdr().isPVgood()<<endl;
	 
	 num_of_Vtx->Fill(Event->evtHdr().nVtx(),hweight);
	 /// Keep events with PVgood
	 if (Event->evtHdr().isPVgood() != 1) continue;
	 
	 for(unsigned j=0; j< n_genJets; ++j){
	   if(Event->genjet(j).pt()<mMinPt) continue;
	   if(fabs(Event->genjet(j).Rapidity())>mYMax) continue;
	   GENjet_ok[j]=1;
	   
	   //if(Event->genjet(j).pt()>=50){
	     GenJets++;
	     
	     if(fabs(Event->genjet(j).Rapidity())<=0.5) pt_GENInclJet_1bin->Fill(Event->genjet(j).pt(),hweight);
	     if(fabs(Event->genjet(j).Rapidity())>0.5 && fabs(Event->genjet(j).Rapidity())<=1.0) pt_GENInclJet_2bin->Fill(Event->genjet(j).pt(),hweight);
	     if(fabs(Event->genjet(j).Rapidity())>1.0 && fabs(Event->genjet(j).Rapidity())<=1.5) pt_GENInclJet_3bin->Fill(Event->genjet(j).pt(),hweight);
	     if(fabs(Event->genjet(j).Rapidity())>1.5 && fabs(Event->genjet(j).Rapidity())<=2.0) pt_GENInclJet_4bin->Fill(Event->genjet(j).pt(),hweight);
	     if(fabs(Event->genjet(j).Rapidity())>2.0 && fabs(Event->genjet(j).Rapidity())<=2.5) pt_GENInclJet_5bin->Fill(Event->genjet(j).pt(),hweight);
	     if(fabs(Event->genjet(j).Rapidity())>2.5 && fabs(Event->genjet(j).Rapidity())<=3.0) pt_GENInclJet_6bin->Fill(Event->genjet(j).pt(),hweight);
	     if(fabs(Event->genjet(j).Rapidity())>3.2 && fabs(Event->genjet(j).Rapidity())<=4.7) pt_GENInclJet_7bin->Fill(Event->genjet(j).pt(),hweight);
	     
	     if(Event->GenJetFlavour(j)>=4.99 && Event->GenJetFlavour(j)<=5.01){
	       if(fabs(Event->genjet(j).Rapidity())<=0.5) pt_GENBTagJet_1bin->Fill(Event->genjet(j).pt(),hweight);
	       if(fabs(Event->genjet(j).Rapidity())>0.5 && fabs(Event->genjet(j).Rapidity())<=1.0) pt_GENBTagJet_2bin->Fill(Event->genjet(j).pt(),hweight);
	       if(fabs(Event->genjet(j).Rapidity())>1.0 && fabs(Event->genjet(j).Rapidity())<=1.5) pt_GENBTagJet_3bin->Fill(Event->genjet(j).pt(),hweight);
	       if(fabs(Event->genjet(j).Rapidity())>1.5 && fabs(Event->genjet(j).Rapidity())<=2.0) pt_GENBTagJet_4bin->Fill(Event->genjet(j).pt(),hweight);
	       if(fabs(Event->genjet(j).Rapidity())>2.0 && fabs(Event->genjet(j).Rapidity())<=2.5) pt_GENBTagJet_5bin->Fill(Event->genjet(j).pt(),hweight);
	       if(fabs(Event->genjet(j).Rapidity())>2.5 && fabs(Event->genjet(j).Rapidity())<=3.0) pt_GENBTagJet_6bin->Fill(Event->genjet(j).pt(),hweight);
	       if(fabs(Event->genjet(j).Rapidity())>3.2 && fabs(Event->genjet(j).Rapidity())<=4.7) pt_GENBTagJet_7bin->Fill(Event->genjet(j).pt(),hweight);
	     }
	   //}
	 }
	 
	 Multiplicity_GENJet->Fill(GenJets,hweight);
       
	 /// Keep events where leading Jet[0] and Jet[1] survived cuts
	 if((GENjet_ok[0]==1)&&(GENjet_ok[1]==1)) {
	   
	   ///////////////////////////////////// Measurement with Gen Jets ///////////////////////////////////////////////////
	   
	   pt0_GENJet->Fill(Event->genjet(0).pt(),hweight);
	   pt1_GENJet->Fill(Event->genjet(1).pt(),hweight);
	   y0_GENJet->Fill(Event->genjet(0).Rapidity(),hweight);
	   y1_GENJet->Fill(Event->genjet(1).Rapidity(),hweight);
	   phi0_GENJet->Fill(Event->genjet(0).phi(),hweight);
	   phi1_GENJet->Fill(Event->genjet(1).phi(),hweight);
	   
	 } //end of GEN Jets
       
       
	 /////////////////////////////////////////////////////////////////////////////////////////////
	 /// PFJets
	 /////////////////////////////////////////////Vertex!!!!////////////////////////////////////////////////////////
	 
	 /// Vertex selection
	 TruePileUpMC->Fill(Event->evtHdr().trpu(),hweight);
	 PileUpVSVertex->Fill(TruePileUpInteger,Event->evtHdr().nVtx(),hweight);

	 //cout<<Event->evtHdr().pu()<<" "<<Event->evtHdr().ootpuEarly()<<" "<<Event->evtHdr().ootpuLate()<<" "<<Event->evtHdr().intpu()<<" "<<Event->evtHdr().nbx()<<endl;

	 if(mprintOk==1){
	   printf("Number of PFJets=%d\n",Event->nPFJetsCHS());
	   for(unsigned j=0; j<Event->nPFJetsCHS(); ++j){
	     printf("j=%2d ptCor=%8.3f  y=%6.3f  phi=%6.3f\n",j,Event->pfjetchs(j).ptCor(),Event->pfjetchs(j).y(),Event->pfjetchs(j).phi());}
	 }
       
	 ///Apply Jet cuts.  Very Deteral to all existing DET Jets
	  
	 //if (n_PFJets != 4 ) continue; // 4jet event for P8vsMadGraph 
	 //To require a two-jet event
	 //if (Event->pfjetchs(2).ptCor() > 114 ) continue;      // Suggested by Panos
	 
	 //To require a three-jet event
	 //if (Event->pfjetchs(3).ptCor() > 114 ) continue;      // Suggested by Panos
             
	 //To require a four-jet event
	 //if (Event->pfjetchs(4).ptCor() > 114 ) continue;      // Suggested by Panos
	 
	 
	 if (Event->pfjetchs(0).ptCor()/pthat > 3.5 ) continue;      // Suggested by Panos

	 for(unsigned j=0; j< Event->nPFJetsCHS(); ++j){
	   if(Event->pfjetchs(j).ptCor()<mMinPt) continue;
	   if(fabs(Event->pfjetchs(j).y())>mYMax) continue;

	   //if(mLowPileUp) {if (TruePileUpInteger<=5) continue;}
	   
	   //if(Event->pfjetchs(j).ptCor()>=50 && TruePileUpInteger<=10){
	   //if(Event->pfjetchs(j).ptCor() >= 20){

	     DetJets++;
	     DETjet_ok[j]=1;
	     
	     if(Event->pfjetchs(j).csv()>0.679)
		    DetJetsBTags++;

	     pt0_DETInclJet->Fill(Event->pfjetchs(j).ptCor(),hweight);
	     pt0_DETInclJetUncor->Fill(Event->pfjetchs(j).pt(),hweight);
	     y0_DETInclJet->Fill(Event->pfjetchs(j).y(),hweight);
	     phi0_DETInclJet->Fill(Event->pfjetchs(j).phi(),hweight);

	     //new histograms
	     Chargedhf0_DETJet->Fill(Event->pfjetchs(j).chf(),hweight);//charged hadron energy fraction
	     Neutralhf0_DETJet->Fill(Event->pfjetchs(j).nhf(),hweight);//neutral hadron energy fraction
	     Chargedef0_DETJet->Fill(Event->pfjetchs(j).cemf(),hweight);//charged em energy fraction
	     Photonef0_DETJet->Fill(Event->pfjetchs(j).nemf(),hweight);//neutral em energy fraction
	     Hadronef0_DETJet->Fill(Event->pfjetchs(j).hf_hf(),hweight);//hadron fraction energy fraction
	     Electromagneticef0_DETJet->Fill(Event->pfjetchs(j).hf_phf(),hweight);//electromagnetic fraction energy fraction
	     Muonef0_DETJet->Fill(Event->pfjetchs(j).muf(),hweight); //muon fraction energy fraction
	 
	     ChargedhMultiplicity0_DETJet->Fill(Event->pfjetchs(j).chm(),hweight);//charged hadron multiplicity
	     NeutralhMultiplicity0_DETJet->Fill(Event->pfjetchs(j).nhm(),hweight);//neutral hadron multiplicity
	     ChargedeMultiplicity0_DETJet->Fill(Event->pfjetchs(j).elm(),hweight);//electron em multiplicity
	     PhotoneMultiplicity0_DETJet->Fill(Event->pfjetchs(j).phm(),hweight);//neutral em multiplicity
	     HadroneMultiplicity0_DETJet->Fill(Event->pfjetchs(j).hf_hm(),hweight);//hadron fraction multiplicity
	     ElectromagneticeMultiplicity0_DETJet->Fill(Event->pfjetchs(j).hf_phm(),hweight);//electromagnetic fraction multiplicity
	     MuoneMultiplicity0_DETJet->Fill(Event->pfjetchs(j).mum(),hweight); //muon fraction multiplicity
	    
	     if(fabs(Event->pfjetchs(j).y())<=0.5 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) {pt_DETInclJet_1bin->Fill(Event->pfjetchs(j).ptCor(),hweight); pt_DETInclJetUP_1bin->Fill(Event->pfjetchs(j).ptCor()+Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);pt_DETInclJetDOWN_1bin->Fill(Event->pfjetchs(j).ptCor()-Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);}
	     if(fabs(Event->pfjetchs(j).y())>0.5 && fabs(Event->pfjetchs(j).y())<=1.0 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) {pt_DETInclJet_2bin->Fill(Event->pfjetchs(j).ptCor(),hweight);pt_DETInclJetUP_2bin->Fill(Event->pfjetchs(j).ptCor()+Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);pt_DETInclJetDOWN_2bin->Fill(Event->pfjetchs(j).ptCor()-Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);}
	     if(fabs(Event->pfjetchs(j).y())>1.0 && fabs(Event->pfjetchs(j).y())<=1.5 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) {pt_DETInclJet_3bin->Fill(Event->pfjetchs(j).ptCor(),hweight);pt_DETInclJetUP_3bin->Fill(Event->pfjetchs(j).ptCor()+Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);pt_DETInclJetDOWN_3bin->Fill(Event->pfjetchs(j).ptCor()-Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);}
	     if(fabs(Event->pfjetchs(j).y())>1.5 && fabs(Event->pfjetchs(j).y())<=2.0 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) {pt_DETInclJet_4bin->Fill(Event->pfjetchs(j).ptCor(),hweight);pt_DETInclJetUP_4bin->Fill(Event->pfjetchs(j).ptCor()+Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);pt_DETInclJetDOWN_4bin->Fill(Event->pfjetchs(j).ptCor()-Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);}
	     if(fabs(Event->pfjetchs(j).y())>2.0 && fabs(Event->pfjetchs(j).y())<=2.5 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) {pt_DETInclJet_5bin->Fill(Event->pfjetchs(j).ptCor(),hweight);pt_DETInclJetUP_5bin->Fill(Event->pfjetchs(j).ptCor()+Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);pt_DETInclJetDOWN_5bin->Fill(Event->pfjetchs(j).ptCor()-Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);}
	     if(fabs(Event->pfjetchs(j).y())>2.5 && fabs(Event->pfjetchs(j).y())<=3.0 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) {pt_DETInclJet_6bin->Fill(Event->pfjetchs(j).ptCor(),hweight);pt_DETInclJetUP_6bin->Fill(Event->pfjetchs(j).ptCor()+Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);pt_DETInclJetDOWN_6bin->Fill(Event->pfjetchs(j).ptCor()-Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);}
	     if(fabs(Event->pfjetchs(j).y())>3.2 && fabs(Event->pfjetchs(j).y())<=4.7) {
	       if(Event->pfjetchs(j).nemf()<0.90 && Event->pfjetchs(j).ncand()>10 && Event->pfmet().met_o_sumet() < 0.3 )//tight JETID forward region
		 {pt_DETInclJet_7bin->Fill(Event->pfjetchs(j).ptCor(),hweight);pt_DETInclJetUP_7bin->Fill(Event->pfjetchs(j).ptCor()+Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);pt_DETInclJetDOWN_7bin->Fill(Event->pfjetchs(j).ptCor()-Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);}
	     }

	     //UncertaintySourcesEvaluation
	     for(int isrc=0;isrc < (int) mJECUncSrcNames.size();isrc++) {
	       //cout<<isrc<<" "<<Event->pfjetchs(j).uncSrc(isrc)<<endl;
	     }

	     ///B-tagged jets

	     if(fabs(Event->pfjetchs(j).y()) < 2.5){
	       CSV_DET->Fill(Event->pfjetchs(j).csv(),hweight);
	       TCHP_DET->Fill(Event->pfjetchs(j).tchp(),hweight);
	       TCHE_DET->Fill(Event->pfjetchs(j).tche(),hweight);
	       SimpleSecVertexHE_DET->Fill(Event->pfjetchs(j).simplesecvertexhe(),hweight);
	       SimpleSecVertexHP_DET->Fill(Event->pfjetchs(j).simplesecvertexhp(),hweight);
	     }

	     if(Event->pfjetchs(j).csv() > 0.679){
	       if(fabs(Event->pfjetchs(j).y())<=0.5 && Event->pfjetchs(j).tightID()) {pt_DETBTagJet_1bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
	       if(fabs(Event->pfjetchs(j).y())>0.5 && fabs(Event->pfjetchs(j).y())<=1.0 && Event->pfjetchs(j).tightID()) {pt_DETBTagJet_2bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
	       if(fabs(Event->pfjetchs(j).y())>1.0 && fabs(Event->pfjetchs(j).y())<=1.5 && Event->pfjetchs(j).tightID()) {pt_DETBTagJet_3bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
	       if(fabs(Event->pfjetchs(j).y())>1.5 && fabs(Event->pfjetchs(j).y())<=2.0 && Event->pfjetchs(j).tightID()) {pt_DETBTagJet_4bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
	       if(fabs(Event->pfjetchs(j).y())>2.0 && fabs(Event->pfjetchs(j).y())<=2.5 && Event->pfjetchs(j).tightID()) {pt_DETBTagJet_5bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
	       if(fabs(Event->pfjetchs(j).y())>2.5 && fabs(Event->pfjetchs(j).y())<=3.0 && Event->pfjetchs(j).tightID()) {pt_DETBTagJet_6bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
	       if(fabs(Event->pfjetchs(j).y())>3.2 && fabs(Event->pfjetchs(j).y())<=4.7) {
		 if(Event->pfjetchs(j).nemf()<0.90 && Event->pfjetchs(j).ncand()>10)//tight JETID forward region
		   {pt_DETBTagJet_7bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
	       }
	     }
	   }
	 
	 if(DETjet_ok[0]==1){
	   TruePileUpMCInteger->Fill(TruePileUpInteger,hweight);
	   num_of_VtxGood->Fill(Event->evtHdr().nVtxGood(),hweight);
	   
	   Multiplicity_DETJetBTags->Fill(DetJetsBTags,hweight);
	   if(DetJetsBTags>=2){
	     Multiplicity_DETJetBTagsInclusive->Fill(DetJets,hweight);
	   }
	   
	   if(Event->pfjetchs(0).csv()>0.679){
	     pt0_DETBTagJetInclusive->Fill(Event->pfjetchs(0).ptCor(),hweight);
	   }
	   
	   if(fabs(Event->pfjetchs(0).y())<=3.0 && Event->pfjetchs(0).tightID()){
	     MET_DET->Fill(Event->pfmet().met(),hweight); 
	     METPhi_DET->Fill(fabs(Event->pfmet().phi()),hweight); 
	     FractionMET_DET->Fill(Event->pfmet().met_o_sumet(),hweight); 
	   }

       // MET fraction for each rap. bin

	   if(fabs(Event->pfjetchs(0).y()) <= 0.5 && Event->pfjetchs(0).tightID())
	     FractionMET_DET_1bin->Fill(Event->pfmet().met_o_sumet(),hweight); 
       if(fabs(Event->pfjetchs(0).y()) > 0.5 && fabs(Event->pfjetchs(0).y()) <= 1.0  && Event->pfjetchs(0).tightID())
	     FractionMET_DET_2bin->Fill(Event->pfmet().met_o_sumet(),hweight); 
       if(fabs(Event->pfjetchs(0).y()) > 1.0 && fabs(Event->pfjetchs(0).y()) <= 1.5  && Event->pfjetchs(0).tightID())
	     FractionMET_DET_3bin->Fill(Event->pfmet().met_o_sumet(),hweight); 
       if(fabs(Event->pfjetchs(0).y()) > 1.5 && fabs(Event->pfjetchs(0).y()) <= 2.0  && Event->pfjetchs(0).tightID())
	     FractionMET_DET_4bin->Fill(Event->pfmet().met_o_sumet(),hweight); 
       if(fabs(Event->pfjetchs(0).y()) > 2.0 && fabs(Event->pfjetchs(0).y()) <= 2.5  && Event->pfjetchs(0).tightID())
	     FractionMET_DET_5bin->Fill(Event->pfmet().met_o_sumet(),hweight); 
       if(fabs(Event->pfjetchs(0).y()) > 2.5 && fabs(Event->pfjetchs(0).y()) <= 3.0  && Event->pfjetchs(0).tightID())
	     FractionMET_DET_6bin->Fill(Event->pfmet().met_o_sumet(),hweight); 

       

	   if(fabs(Event->pfjetchs(0).y())>3.2){
	     if(Event->pfjetchs(0).nemf()<0.90 && Event->pfjetchs(0).ncand()>10){//tight JETID forward region
	       MET_DET->Fill(Event->pfmet().met(),hweight); 
	       METPhi_DET->Fill(fabs(Event->pfmet().phi()),hweight); 
	       FractionMET_DET->Fill(Event->pfmet().met_o_sumet(),hweight); 
	       FractionMET_DET_7bin->Fill(Event->pfmet().met_o_sumet(),hweight); 
	     }
	   }
	 }
              
	 /// Keep events where leading Jet[0] and Jet[1] survived cuts
	 if((DETjet_ok[0]==1)&&(DETjet_ok[1]==1)) {
	   
	   ///////////////////////////////////// Measurement with Det Jets ///////////////////////////////////////////////////
	   
	   pt0_DETJet->Fill(Event->pfjetchs(0).ptCor(),hweight);
	   pt0_DETJetUncor->Fill(Event->pfjetchs(0).pt(),hweight);
	   pt1_DETJet->Fill(Event->pfjetchs(1).ptCor(),hweight);
	   pt1_DETJetUncor->Fill(Event->pfjetchs(1).pt(),hweight);
	   y0_DETJet->Fill(Event->pfjetchs(0).y(),hweight);
	   y1_DETJet->Fill(Event->pfjetchs(1).y(),hweight);
	   phi0_DETJet->Fill(Event->pfjetchs(0).phi(),hweight);
	   phi1_DETJet->Fill(Event->pfjetchs(1).phi(),hweight);
	 
	   if(Event->pfjetchs(0).csv()>0.679 && Event->pfjetchs(1).csv()>0.679){  
	     pt0_DETBTagJet->Fill(Event->pfjetchs(0).ptCor(),hweight);
	     pt1_DETBTagJet->Fill(Event->pfjetchs(1).ptCor(),hweight);
	     y0_DETBTagJet->Fill(Event->pfjetchs(0).y(),hweight);
	     y1_DETBTagJet->Fill(Event->pfjetchs(1).y(),hweight);
	     reco::Particle::LorentzVector Lead1=Event->pfjetchs(0).p4();
	     reco::Particle::LorentzVector Lead2=Event->pfjetchs(1).p4();
	     reco::Particle::LorentzVector Pair=Lead1+Lead2;
	     //cout<<Lead1.pt()<<" "<<Lead2.pt()<<" "<<Pair.pt()<<endl;
	     ptPair_DETBTagJet->Fill(Pair.pt(),hweight);
	   }

	 }	  
        
	 /// Vertex selection
	 if(mprintOk==1) cout<<"Vertex info: numVtx="<<Event->evtHdr().nVtx()<<"  numVtxGood="<<Event->evtHdr().nVtxGood()<<"  isPVgood()="<<Event->evtHdr().isPVgood()<<"   pfRho="<<Event->evtHdr().pfRho()<<endl;
       
	 /// Dump all jets No cuts only isPVgood()
	 if(mprintOk==1){
	   printf("Number of PFJets=%d\n",Event->nPFJetsCHS());
	   for(unsigned j=0; j<Event->nPFJetsCHS(); ++j){
	     printf("j=%2d  pt=%8.3f ptCor=%8.3f y=%6.3f  phi=%6.3f   cor=%6.3f   tightID=%d\n",j,Event->pfjetchs(j).pt(),Event->pfjetchs(j).ptCor(),Event->pfjetchs(j).y(),Event->pfjetchs(j).phi(),Event->pfjetchs(j).cor(),Event->pfjetchs(j).tightID());
	   }
	 }
       
	 // PURITY AND STABILITY
	 for(unsigned j=0; j< Event->nPFJetsCHS(); ++j){
	   if(Event->pfjetchs(j).ptCor()<mMinPt) continue;
	   if(fabs(Event->pfjetchs(j).y())>mYMax) continue;
	 
	   int histoEntry=0;
	   double matchJet=1000;
	   double matchJetCounter=1000;
	   
	   for(unsigned i=0; i< n_genJets; ++i){
	     if(Event->genjet(i).pt()<mMinPt) continue;
	     if(fabs(Event->genjet(i).Rapidity())>mYMax) continue;
	     
	     double deltaPhiMatch=Event->genjet(i).phi()-Event->pfjetchs(j).phi();
	     if(deltaPhiMatch<-TMath::Pi()) deltaPhiMatch=deltaPhiMatch+2*TMath::Pi();
	     if(deltaPhiMatch>TMath::Pi()) deltaPhiMatch=deltaPhiMatch-2*TMath::Pi();
	     deltaPhiMatch=fabs(deltaPhiMatch);
	     
	     matchJetCounter=sqrt(pow(Event->genjet(i).Rapidity()-Event->pfjetchs(j).y(),2)+pow(deltaPhiMatch,2));
	     
	     if(matchJetCounter<=matchJet) matchJet=matchJetCounter;
	     
	     if (matchJetCounter<0.2 && histoEntry==0){
	       histoEntry=1;
	       PF_MatchedInclusiveJets->Fill(Event->pfjetchs(j).ptCor(),hweight);
	       Gen_MatchedInclusiveJets->Fill(Event->genjet(i).pt(),hweight);
	       TwoD_MatchedInclusiveJets->Fill(Event->pfjetchs(j).ptCor(),Event->genjet(i).pt(),hweight);
	       
	       //Separation in eta bins
	       if(fabs(Event->pfjetchs(j).y())<=0.5 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) { PF_MatchedInclusiveJets_1bin->Fill(Event->pfjetchs(j).ptCor(),hweight); Gen_MatchedInclusiveJets_1bin->Fill(Event->genjet(i).pt(),hweight); TwoD_MatchedInclusiveJets_1bin->Fill(Event->pfjetchs(j).ptCor(),Event->genjet(i).pt(),hweight); } 
	       if(fabs(Event->pfjetchs(j).y())>0.5 && fabs(Event->pfjetchs(j).y())<=1.0 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3)  { PF_MatchedInclusiveJets_2bin->Fill(Event->pfjetchs(j).ptCor(),hweight); Gen_MatchedInclusiveJets_2bin->Fill(Event->genjet(i).pt(),hweight); TwoD_MatchedInclusiveJets_2bin->Fill(Event->pfjetchs(j).ptCor(),Event->genjet(i).pt(),hweight); } 
	       if(fabs(Event->pfjetchs(j).y())>1.0 && fabs(Event->pfjetchs(j).y())<=1.5 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3)  { PF_MatchedInclusiveJets_3bin->Fill(Event->pfjetchs(j).ptCor(),hweight); Gen_MatchedInclusiveJets_3bin->Fill(Event->genjet(i).pt(),hweight); TwoD_MatchedInclusiveJets_3bin->Fill(Event->pfjetchs(j).ptCor(),Event->genjet(i).pt(),hweight); } 
	       if(fabs(Event->pfjetchs(j).y())>1.5 && fabs(Event->pfjetchs(j).y())<=2.0 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3)  { PF_MatchedInclusiveJets_4bin->Fill(Event->pfjetchs(j).ptCor(),hweight); Gen_MatchedInclusiveJets_4bin->Fill(Event->genjet(i).pt(),hweight); TwoD_MatchedInclusiveJets_4bin->Fill(Event->pfjetchs(j).ptCor(),Event->genjet(i).pt(),hweight); } 
	       if(fabs(Event->pfjetchs(j).y())>2.0 && fabs(Event->pfjetchs(j).y())<=2.5 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3)  { PF_MatchedInclusiveJets_5bin->Fill(Event->pfjetchs(j).ptCor(),hweight); Gen_MatchedInclusiveJets_5bin->Fill(Event->genjet(i).pt(),hweight); TwoD_MatchedInclusiveJets_5bin->Fill(Event->pfjetchs(j).ptCor(),Event->genjet(i).pt(),hweight); } 
	       if(fabs(Event->pfjetchs(j).y())>2.5 && fabs(Event->pfjetchs(j).y())<=3.0 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3)  { PF_MatchedInclusiveJets_6bin->Fill(Event->pfjetchs(j).ptCor(),hweight); Gen_MatchedInclusiveJets_6bin->Fill(Event->genjet(i).pt(),hweight); TwoD_MatchedInclusiveJets_6bin->Fill(Event->pfjetchs(j).ptCor(),Event->genjet(i).pt(),hweight); } 
	       if(fabs(Event->pfjetchs(j).y())>3.2 && fabs(Event->pfjetchs(j).y())<=4.7 && Event->pfmet().met_o_sumet() < 0.3)  {
		 if(Event->pfjetchs(j).nemf()<0.90 && Event->pfjetchs(j).ncand()>10){//tight JETID forward region
		   PF_MatchedInclusiveJets_7bin->Fill(Event->pfjetchs(j).ptCor(),hweight); Gen_MatchedInclusiveJets_7bin->Fill(Event->genjet(i).pt(),hweight); TwoD_MatchedInclusiveJets_7bin->Fill(Event->pfjetchs(j).ptCor(),Event->genjet(i).pt(),hweight); } 
	       }
	       
	       //Evaluation of resolution
	       double resolution=fabs(Event->pfjetchs(j).ptCor()-Event->genjet(i).pt())/Event->genjet(i).pt();
	       double resolutionNoAbs=(Event->pfjetchs(j).ptCor()-Event->genjet(i).pt())/Event->genjet(i).pt();
	       Resolution1D->Fill(resolutionNoAbs,hweight);
	       ResolutionInclusiveJets->Fill(Event->pfjetchs(j).ptCor(),resolution,hweight);
	       //Fill that for each eta bins
	       //ResolutionHistoPt1bin[0]->Fill(resolution,hweight);

	       if(fabs(Event->pfjetchs(j).y())<=0.5 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3)  ResolutionInclusiveJets_1bin->Fill(Event->pfjetchs(j).ptCor(),resolution,hweight);
	       if(fabs(Event->pfjetchs(j).y())>0.5 && fabs(Event->pfjetchs(j).y())<=1.0 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) ResolutionInclusiveJets_2bin->Fill(Event->pfjetchs(j).ptCor(),resolution,hweight); 
	       if(fabs(Event->pfjetchs(j).y())>1.0 && fabs(Event->pfjetchs(j).y())<=1.5 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) ResolutionInclusiveJets_3bin->Fill(Event->pfjetchs(j).ptCor(),resolution,hweight);
	       if(fabs(Event->pfjetchs(j).y())>1.5 && fabs(Event->pfjetchs(j).y())<=2.0 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) ResolutionInclusiveJets_4bin->Fill(Event->pfjetchs(j).ptCor(),resolution,hweight);
	       if(fabs(Event->pfjetchs(j).y())>2.0 && fabs(Event->pfjetchs(j).y())<=2.5 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) ResolutionInclusiveJets_5bin->Fill(Event->pfjetchs(j).ptCor(),resolution,hweight); 
	       if(fabs(Event->pfjetchs(j).y())>2.5 && fabs(Event->pfjetchs(j).y())<=3.0 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) ResolutionInclusiveJets_6bin->Fill(Event->pfjetchs(j).ptCor(),resolution,hweight); 
	       if(fabs(Event->pfjetchs(j).y())>3.2 && fabs(Event->pfjetchs(j).y())<=4.7 && Event->pfmet().met_o_sumet() < 0.3){ 
		 if(Event->pfjetchs(j).nemf()<0.90 && Event->pfjetchs(j).ncand()>10){//tight JETID forward region
		   {ResolutionInclusiveJets_7bin->Fill(Event->pfjetchs(j).ptCor(),resolution,hweight); ResolutionForward1D->Fill(resolutionNoAbs,hweight);}
		 }
	       }
	     }
	   }
	   
	   DeltaR_Jets->Fill(matchJet,hweight);
	   
	   if (matchJet>0.2){
	     PF_FakeInclusiveJets->Fill(Event->pfjetchs(j).ptCor(),hweight);
	     
	     if(fabs(Event->pfjetchs(j).y())<=0.5 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) PF_FakeInclusiveJets_1bin->Fill(Event->pfjetchs(j).ptCor(),hweight);
	     if(fabs(Event->pfjetchs(j).y())>0.5 && fabs(Event->pfjetchs(j).y())<=1.0 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) PF_FakeInclusiveJets_2bin->Fill(Event->pfjetchs(j).ptCor(),hweight);
	     if(fabs(Event->pfjetchs(j).y())>1.0 && fabs(Event->pfjetchs(j).y())<=1.5 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) PF_FakeInclusiveJets_3bin->Fill(Event->pfjetchs(j).ptCor(),hweight);
	     if(fabs(Event->pfjetchs(j).y())>1.5 && fabs(Event->pfjetchs(j).y())<=2.0 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) PF_FakeInclusiveJets_4bin->Fill(Event->pfjetchs(j).ptCor(),hweight);
	     if(fabs(Event->pfjetchs(j).y())>2.0 && fabs(Event->pfjetchs(j).y())<=2.5 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) PF_FakeInclusiveJets_5bin->Fill(Event->pfjetchs(j).ptCor(),hweight);
	     if(fabs(Event->pfjetchs(j).y())>2.5 && fabs(Event->pfjetchs(j).y())<=3.0 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) PF_FakeInclusiveJets_6bin->Fill(Event->pfjetchs(j).ptCor(),hweight);
	     if(fabs(Event->pfjetchs(j).y())>3.2 && fabs(Event->pfjetchs(j).y())<=4.7 && Event->pfmet().met_o_sumet() < 0.3 ) PF_FakeInclusiveJets_7bin->Fill(Event->pfjetchs(j).ptCor(),hweight);
	     
	   }
	 }
	 
	 // Missed Events at the Gen Level
	 for(unsigned i=0; i< n_genJets; ++i){
	   if(Event->genjet(i).pt()<mMinPt) continue;
	   if(fabs(Event->genjet(i).Rapidity())>mYMax) continue;
	   double matchJet=1000;
	   double matchJetCounter=1000;
	   
	   for(unsigned j=0; j< Event->nPFJetsCHS(); ++j){
	     if(Event->pfjetchs(j).ptCor()<mMinPt) continue;
	     if(fabs(Event->pfjetchs(j).y())>mYMax) continue;
	   
	     if(fabs(Event->pfjetchs(j).y())<=3.0 && !Event->pfjetchs(j).tightID()) continue;

	     double deltaPhiMatch=Event->genjet(i).phi()-Event->pfjetchs(j).phi();
	     if(deltaPhiMatch<-TMath::Pi()) deltaPhiMatch=deltaPhiMatch+2*TMath::Pi();
	     if(deltaPhiMatch>TMath::Pi()) deltaPhiMatch=deltaPhiMatch-2*TMath::Pi();
	     deltaPhiMatch=fabs(deltaPhiMatch);
	     
	     matchJetCounter=sqrt(pow(Event->genjet(i).Rapidity()-Event->pfjetchs(j).y(),2)+pow(deltaPhiMatch,2));
	     if(matchJetCounter<=matchJet) matchJet=matchJetCounter;
	   }
	   
	   if (matchJet>0.2){
	     Gen_MissInclusiveJets->Fill(Event->genjet(i).pt(),hweight);
	     
	     if(fabs(Event->genjet(i).Rapidity())<=0.5 && Event->pfmet().met_o_sumet() < 0.3) Gen_MissInclusiveJets_1bin->Fill(Event->genjet(i).pt(),hweight);
	     if(fabs(Event->genjet(i).Rapidity())>0.5 && fabs(Event->genjet(i).Rapidity())<=1.0 && Event->pfmet().met_o_sumet() < 0.3) Gen_MissInclusiveJets_2bin->Fill(Event->genjet(i).pt(),hweight);
	     if(fabs(Event->genjet(i).Rapidity())>1.0 && fabs(Event->genjet(i).Rapidity())<=1.5 && Event->pfmet().met_o_sumet() < 0.3) Gen_MissInclusiveJets_3bin->Fill(Event->genjet(i).pt(),hweight);
	     if(fabs(Event->genjet(i).Rapidity())>1.5 && fabs(Event->genjet(i).Rapidity())<=2.0 && Event->pfmet().met_o_sumet() < 0.3) Gen_MissInclusiveJets_4bin->Fill(Event->genjet(i).pt(),hweight);
	     if(fabs(Event->genjet(i).Rapidity())>2.0 && fabs(Event->genjet(i).Rapidity())<=2.5 && Event->pfmet().met_o_sumet() < 0.3) Gen_MissInclusiveJets_5bin->Fill(Event->genjet(i).pt(),hweight);
	     if(fabs(Event->genjet(i).Rapidity())>2.5 && fabs(Event->genjet(i).Rapidity())<=3.0 && Event->pfmet().met_o_sumet() < 0.3) Gen_MissInclusiveJets_6bin->Fill(Event->genjet(i).pt(),hweight);
	     if(fabs(Event->genjet(i).Rapidity())>3.2 && fabs(Event->genjet(i).Rapidity())<=4.7 && Event->pfmet().met_o_sumet() < 0.3) Gen_MissInclusiveJets_7bin->Fill(Event->genjet(i).pt(),hweight);
	     
	     /*if(fabs(Event->genjet(i).y())<=0.5) { resp_jetpt1etabin->Miss(Event->genjet(i).pt(),hweight);}
	       if(fabs(Event->genjet(i).y())>0.5 && fabs(Event->genjet(i).y())<=1.0) { resp_jetpt2etabin->Miss(Event->genjet(i).pt(),hweight);}
	       if(fabs(Event->genjet(i).y())>1.0 && fabs(Event->genjet(i).y())<=1.5) { resp_jetpt3etabin->Miss(Event->genjet(i).pt(),hweight);}
	       if(fabs(Event->genjet(i).y())>1.5 && fabs(Event->genjet(i).y())<=2.0) { resp_jetpt4etabin->Miss(Event->genjet(i).pt(),hweight);}
	       if(fabs(Event->genjet(i).y())>2.0 && fabs(Event->genjet(i).y())<=2.5) { resp_jetpt5etabin->Miss(Event->genjet(i).pt(),hweight);}
	       if(fabs(Event->genjet(i).y())>2.5 && fabs(Event->genjet(i).y())<=3.0) { resp_jetpt6etabin->Miss(Event->genjet(i).pt(),hweight);}
	       if(fabs(Event->genjet(i).y())>3.2 && fabs(Event->genjet(i).y())<=4.7) { resp_jetpt7etabin->Miss(Event->genjet(i).pt(),hweight);}*/
	   }
	 }
       }
     
       // Trigger Ranges code
       if(!mIsMCarlo){
	 
	 double leading_pt = -10.0;
	 double leading_eta = -10.0;

	 num_of_Vtx->Fill(Event->evtHdr().nVtx(),hweight);
	 /// Keep events with PVgood
       
	 if(mprintOk==1){
           for(unsigned j=0; j<Event->nPFJets(); ++j){
	     if(Event->evtHdr().runNo()==251252){
	       if(Event->evtHdr().lumi()==158 && Event->evtHdr().event()==95957128) {
		 cout<<Event->evtHdr().event()<<endl;
		 printf("j=%2d  ptCor=%8.3f  y=%6.3f  phi=%6.3f\n",j,Event->pfjet(j).ptCor(),Event->pfjet(j).y(),Event->pfjet(j).phi());}
	       if(Event->evtHdr().lumi()==158 && Event->evtHdr().event()==95734424) {
		 cout<<Event->evtHdr().event()<<endl;
                 printf("j=%2d  ptCor=%8.3f  y=%6.3f  phi=%6.3f\n",j,Event->pfjet(j).ptCor(),Event->pfjet(j).y(),Event->pfjet(j).phi());}
	       if(Event->evtHdr().lumi()==158 && Event->evtHdr().event()==95800780) {
		 cout<<Event->evtHdr().event()<<endl;
                 printf("j=%2d  ptCor=%8.3f  y=%6.3f  phi=%6.3f\n",j,Event->pfjet(j).ptCor(),Event->pfjet(j).y(),Event->pfjet(j).phi());}
	     }
           }
	 }
	 
	 if (Event->evtHdr().isPVgood() != 1) continue;
	 
	 for (int j=0; j<nJetTrig; j++){
	   
	   hltcut[j]=false;
	   l1cut[j]=false;
	   
	   if (hltPassj[j])
	     {
	       //----------------- L1 Theshold Checking --------------------------//
	       for ( unsigned l1iobj=0; l1iobj<Event->nL1Obj(ihltj[j]); l1iobj++ )
		 {
		   //cout << "L1 object id: "<<l1iobj << " -> " << Event->l1obj(ihltj[j],l1iobj).pt() << " > " << ATrig[j] << endl;
		   if (Event->l1obj(ihltj[j],l1iobj).pt() > ATrig[j]) //threshold for L1 trigger
		     {
		       l1cut[j] = true;
		     } 	       
		 } 
	       
	       //----------------- HLT Theshold Checking --------------------------//
	       for ( unsigned hltiobj=0; hltiobj<Event->nHLTObj(ihltj[j]); hltiobj++ )
		 {
		   //cout << "HLT object id: " << hltiobj << " -> " << Event->hltobj(ihltj[j],hltiobj).pt() << " > " << HLTJetPtNThres[j] << endl;
		   if (Event->hltobj(ihltj[j],hltiobj).pt() > HLTJetPtNThres[j]) //threshold for HLT trigger
		     {
		       hltcut[j] = true;
		     } 	       
		 }
	     }
	 }
     
	 //Trigger tag-and-probe
	 
	 if(Event->nPFJetsCHS()>=2) {
	    int selectEventsLead=0;
	    int selectEventsSubLead=0;
	   
	    if(fabs(Event->pfjetchs(0).eta())<1.3){
	     
	      for ( unsigned hltiobj=0; hltiobj<Event->nHLTObj(ihltj[1]); hltiobj++ )
		{
		  
		  double deltaPhiMatch=Event->pfjetchs(0).phi()-Event->hltobj(ihltj[1],hltiobj).phi();
		  if(deltaPhiMatch<-TMath::Pi()) deltaPhiMatch=deltaPhiMatch+2*TMath::Pi();
		  if(deltaPhiMatch>TMath::Pi()) deltaPhiMatch=deltaPhiMatch-2*TMath::Pi();
		  deltaPhiMatch=fabs(deltaPhiMatch);
		  
		  double matchJetCounter=sqrt(pow(Event->pfjetchs(0).eta()-Event->hltobj(ihltj[1],hltiobj).eta(),2)+pow(deltaPhiMatch,2));

		  if (matchJetCounter<0.5 && Event->hltobj(ihltj[1],hltiobj).pt()>10){	
		    double ResHLT=(Event->pfjetchs(0).ptCor()-Event->hltobj(ihltj[1],hltiobj).pt())/Event->pfjetchs(0).ptCor();
		    ResolutionTagAndProbe->Fill(Event->pfjetchs(0).ptCor(),ResHLT,hweight);		       
		  }
		  
		  if (matchJetCounter<0.5 && Event->hltobj(ihltj[1],hltiobj).pt()>60){	
		    selectEventsLead=1;
		  }
		    
		}
	    }
	    
	    if(fabs(Event->pfjetchs(1).eta())<1.3){
	      for ( unsigned hltiobj=0; hltiobj<Event->nHLTObj(ihltj[1]); hltiobj++ )
		{
		  
		  double deltaPhiMatch=Event->pfjetchs(1).phi()-Event->hltobj(ihltj[1],hltiobj).phi();
		  if(deltaPhiMatch<-TMath::Pi()) deltaPhiMatch=deltaPhiMatch+2*TMath::Pi();
		  if(deltaPhiMatch>TMath::Pi()) deltaPhiMatch=deltaPhiMatch-2*TMath::Pi();
		  deltaPhiMatch=fabs(deltaPhiMatch);
		  
		  double matchJetCounter=sqrt(pow(Event->pfjetchs(1).eta()-Event->hltobj(ihltj[1],hltiobj).eta(),2)+pow(deltaPhiMatch,2));
		  
		  if (matchJetCounter<0.5 && Event->hltobj(ihltj[1],hltiobj).pt()>10){
                    double ResHLT=(Event->pfjetchs(1).ptCor()-Event->hltobj(ihltj[1],hltiobj).pt())/Event->pfjetchs(1).ptCor();
                    ResolutionTagAndProbe->Fill(Event->pfjetchs(1).ptCor(),ResHLT,hweight);
                  }

		  if (matchJetCounter<0.5 && Event->hltobj(ihltj[1],hltiobj).pt()>60) //threshold for HLT trigger                                                             
		    {		     
		      selectEventsSubLead=1;
		    }
		}
	    }
	 
	    int TriggerDenLead=0;
	    int TriggerDenLeadMatch=0;
	    int TriggerNumLead=0;
	    
	    int TriggerDenSubLead=0;
	    int TriggerDenSubLeadMatch=0;
	    int TriggerNumSubLead=0;
	 
	    if(selectEventsLead){
	      
	      if(fabs(Event->pfjetchs(1).eta())<1.5){
		
		double deltaPhiMatch=Event->pfjetchs(0).phi()-Event->pfjetchs(1).phi();
		if(deltaPhiMatch<-TMath::Pi()) deltaPhiMatch=deltaPhiMatch+2*TMath::Pi();
		if(deltaPhiMatch>TMath::Pi()) deltaPhiMatch=deltaPhiMatch-2*TMath::Pi();
		deltaPhiMatch=fabs(deltaPhiMatch);
		
		if(Event->nPFJetsCHS()==2 && deltaPhiMatch>2.7){TriggerDenLead=1;}
		double ptBalance=0.3*(Event->pfjetchs(0).ptCor()+Event->pfjetchs(1).ptCor())/2;
		if(Event->nPFJetsCHS()>2 && deltaPhiMatch>2.7 && Event->pfjetchs(2).ptCor()<ptBalance){TriggerDenLead=1;}
		
		if (TriggerDenLead==1){
		  
		  for ( unsigned hltiobj=0; hltiobj<Event->nHLTObj(ihltj[1]); hltiobj++ )
		    {
		      
		      double deltaPhiMatch2=Event->pfjetchs(1).phi()-Event->hltobj(ihltj[1],hltiobj).phi();
		      if(deltaPhiMatch2<-TMath::Pi()) deltaPhiMatch2=deltaPhiMatch2+2*TMath::Pi();
		      if(deltaPhiMatch2>TMath::Pi()) deltaPhiMatch2=deltaPhiMatch2-2*TMath::Pi();
		      deltaPhiMatch2=fabs(deltaPhiMatch2);
		      
		      double matchJetCounter2=sqrt(pow(Event->pfjetchs(1).eta()-Event->hltobj(ihltj[1],hltiobj).eta(),2)+pow(deltaPhiMatch2,2));
		      
		      if (matchJetCounter2<0.2){
			TriggerDenLeadMatch=1;		 
			if(Event->hltobj(ihltj[1],hltiobj).pt()>60){
			  //cout<<Event->pfjetchs(1).ptCor()<<" "<<Event->hltobj(ihltj[1],hltiobj).pt()<<endl;
			  TriggerNumLead=1;
			}
		      }
		    }
		}
	      }
	    }

	    if(selectEventsSubLead){
	      
	      if(fabs(Event->pfjetchs(0).eta())<1.5){
	       
		double deltaPhiMatch=Event->pfjetchs(0).phi()-Event->pfjetchs(1).phi();
		if(deltaPhiMatch<-TMath::Pi()) deltaPhiMatch=deltaPhiMatch+2*TMath::Pi();
		if(deltaPhiMatch>TMath::Pi()) deltaPhiMatch=deltaPhiMatch-2*TMath::Pi();
		deltaPhiMatch=fabs(deltaPhiMatch);
		
		if(Event->nPFJetsCHS()==2 && deltaPhiMatch>2.7){TriggerDenSubLead=1;}
		double ptBalance=0.3*(Event->pfjetchs(0).ptCor()+Event->pfjetchs(1).ptCor())/2;
		if(Event->nPFJetsCHS()>2 && deltaPhiMatch>2.7 && Event->pfjetchs(2).ptCor()<ptBalance){TriggerDenSubLead=1;}
		
		if (TriggerDenSubLead==1){
		  
		  for ( unsigned hltiobj=0; hltiobj<Event->nHLTObj(ihltj[1]); hltiobj++ )
		    {
		      
		      double deltaPhiMatch2=Event->pfjetchs(0).phi()-Event->hltobj(ihltj[1],hltiobj).phi();
		      if(deltaPhiMatch2<-TMath::Pi()) deltaPhiMatch2=deltaPhiMatch2+2*TMath::Pi();
		      if(deltaPhiMatch2>TMath::Pi()) deltaPhiMatch2=deltaPhiMatch2-2*TMath::Pi();
		      deltaPhiMatch2=fabs(deltaPhiMatch2);
		      
		      double matchJetCounter2=sqrt(pow(Event->pfjetchs(0).eta()-Event->hltobj(ihltj[1],hltiobj).eta(),2)+pow(deltaPhiMatch2,2));
		      
		      if(matchJetCounter2<0.2){
			TriggerDenSubLeadMatch=1;
			
			if(Event->hltobj(ihltj[1],hltiobj).pt()>60){
			  TriggerNumSubLead=1;
			}
		      }
		    }
		}
	      }
	    }

	    if(random==1){
	      random=0;
	      if(TriggerDenLeadMatch){TagAndProbeDen->Fill(Event->pfjetchs(1).ptCor(),hweight);}
	      if(TriggerNumLead){TagAndProbeNum->Fill(Event->pfjetchs(1).ptCor(),hweight);}
	    }

	    if(random==0){
	      random=1;
	      if(TriggerDenSubLeadMatch){TagAndProbeDen->Fill(Event->pfjetchs(0).ptCor(),hweight);}
	      if(TriggerNumSubLead){TagAndProbeNum->Fill(Event->pfjetchs(0).ptCor(),hweight);}
	    }
	 }
	 
	 //end of trigger tag-and-probe 
	 
	 for(unsigned int j=0; j<Event->nPFJetsCHS(); j++) {
	   double pt = Event->pfjetchs(j).ptCor();
	   double eta = Event->pfjetchs(j).eta();
	   
	   if (pt >= 20){
	     if((Event->pfjetchs(j).tightID() && fabs(Event->pfjetchs(j).eta())<=3.0) || (fabs(Event->pfjetchs(j).eta())>3.0 && fabs(Event->pfjetchs(j).eta())<=4.7 && Event->pfjetchs(j).nemf()<0.90 && Event->pfjetchs(j).ncand()>10)) {
	       if (leading_pt < pt) { leading_pt = pt; leading_eta = eta; }
	     }
	   }
	 }
	 
	 hweight=1.;
	 
	 if(fabs(leading_eta)<1.5){
	   if(hltPassj[0]){
	     //if(prescalej[0]>0) hweight=hweight*prescalej[0];
	     //if(prescalej[0]<0) hweight=hweight*(-prescalej[0]);
	     hist_leading_pt_all_Jet60->Fill(leading_pt,hweight);
	     hist_leading_eta_all_Jet60->Fill(leading_eta,hweight);
	   }
	   if (hltPassj[0] && hltcut[0]){
	     //if(prescalej[0]>0) hweight=hweight*prescalej[0];
	     //if(prescalej[0]<0) hweight=hweight*(-prescalej[0]);
	     hist_leading_pt_emulated_Jet60->Fill(leading_pt,hweight);
	     hist_leading_eta_emulated_Jet60->Fill(leading_eta,hweight);
	   }
	   
	   if(hltPassj[1]){
	     //if(prescalej[1]>0) hweight=hweight*prescalej[1];
	     //if(prescalej[1]<0) hweight=hweight*(-prescalej[1]);
	     hist_leading_pt_all_Jet80->Fill(leading_pt,hweight);
	     hist_leading_eta_all_Jet80->Fill(leading_eta,hweight);
	   }
	   if (hltPassj[1] && l1cut[1] && hltcut[1]){
	     //if(prescalej[1]>0) hweight=hweight*prescalej[1];
	     //if(prescalej[1]<0) hweight=hweight*(-prescalej[1]);
	     hist_leading_pt_emulated_Jet80->Fill(leading_pt,hweight);
	     hist_leading_eta_emulated_Jet80->Fill(leading_eta,hweight);
	   }
	   
	   if(hltPassj[2]){
	     //if(prescalej[2]>0) hweight=hweight*prescalej[2];
	     //if(prescalej[2]<0) hweight=hweight*(-prescalej[2]);
	     hist_leading_pt_all_Jet140->Fill(leading_pt,hweight);
	     hist_leading_eta_all_Jet140->Fill(leading_eta,hweight);
	   }
	   if (hltPassj[2] && l1cut[2] && hltcut[2]){
	     //if(prescalej[2]>0) hweight=hweight*prescalej[2];
	     //if(prescalej[2]<0) hweight=hweight*(-prescalej[2]);
	     hist_leading_pt_emulated_Jet140->Fill(leading_pt,hweight);
	     hist_leading_eta_emulated_Jet140->Fill(leading_eta,hweight);
	   }

	   if(hltPassj[3]){
	     //if(prescalej[3]>0) hweight=hweight*prescalej[3];
	     //if(prescalej[3]<0) hweight=hweight*(-prescalej[3]);
	     hist_leading_pt_all_Jet200->Fill(leading_pt,hweight);
	     hist_leading_eta_all_Jet200->Fill(leading_eta,hweight);
	   }
	   if (hltPassj[3] && l1cut[3] && hltcut[3]){
	     //if(prescalej[3]>0) hweight=hweight*prescalej[3];
	     //if(prescalej[3]<0) hweight=hweight*(-prescalej[3]);
	     hist_leading_pt_emulated_Jet200->Fill(leading_pt,hweight);
	     hist_leading_eta_emulated_Jet200->Fill(leading_eta,hweight);
	   }
	   
	   if(hltPassj[4]){
	     //if(prescalej[3]>0) hweight=hweight*prescalej[3];
	     //if(prescalej[3]<0) hweight=hweight*(-prescalej[3]);
	     hist_leading_pt_all_Jet260->Fill(leading_pt,hweight);
	     hist_leading_eta_all_Jet260->Fill(leading_eta,hweight);
	   }
	   if (hltPassj[4] && l1cut[4] && hltcut[4]){
	     //if(prescalej[3]>0) hweight=hweight*prescalej[3];
	     //if(prescalej[3]<0) hweight=hweight*(-prescalej[3]);
	     hist_leading_pt_emulated_Jet260->Fill(leading_pt,hweight);
	     hist_leading_eta_emulated_Jet260->Fill(leading_eta,hweight);
	   }
	   
	   if(hltPassj[5]){
	     //if(prescalej[4]>0) hweight=hweight*prescalej[4];
	     //if(prescalej[4]<0) hweight=hweight*(-prescalej[4]);
	     hist_leading_pt_all_Jet320->Fill(leading_pt,hweight);
	     hist_leading_eta_all_Jet320->Fill(leading_eta,hweight);
	   }
	   if (hltPassj[5] && l1cut[5] && hltcut[5]){
	     //if(prescalej[4]>0) hweight=hweight*prescalej[4];
	     //if(prescalej[4]<0) hweight=hweight*(-prescalej[4]);
	     hist_leading_pt_emulated_Jet320->Fill(leading_pt,hweight);
	     hist_leading_eta_emulated_Jet320->Fill(leading_eta,hweight);
	   }
	   
	   if(hltPassj[6]){
	     //if(prescalej[5]>0) hweight=hweight*prescalej[5];
	     //if(prescalej[5]<0) hweight=hweight*(-prescalej[5]);
	     hist_leading_pt_all_Jet400->Fill(leading_pt,hweight);
	     hist_leading_eta_all_Jet400->Fill(leading_eta,hweight);
	   }
	   if (hltPassj[6] && l1cut[6] && hltcut[6]){
	     //if(prescalej[5]>0) hweight=hweight*prescalej[5];
	     //if(prescalej[5]<0) hweight=hweight*(-prescalej[5]);
	     hist_leading_pt_emulated_Jet400->Fill(leading_pt,hweight);
	     hist_leading_eta_emulated_Jet400->Fill(leading_eta,hweight);
	   }
	   
	   if(hltPassj[7]){
	     //if(prescalej[6]>0) hweight=hweight*prescalej[6];
	     //if(prescalej[6]<0) hweight=hweight*(-prescalej[6]);
	     hist_leading_pt_all_Jet450->Fill(leading_pt,hweight);
	     hist_leading_eta_all_Jet450->Fill(leading_eta,hweight);
	   }
	   if (hltPassj[7] && l1cut[7] && hltcut[7]){
	     //if(prescalej[6]>0) hweight=hweight*prescalej[6];
	     //if(prescalej[6]<0) hweight=hweight*(-prescalej[6]);
	     hist_leading_pt_emulated_Jet450->Fill(leading_pt,hweight);
	     hist_leading_eta_emulated_Jet450->Fill(leading_eta,hweight);
	   }
	   
	   if(hltPassj[8]){
	     //if(prescalej[7]>0) hweight=hweight*prescalej[7];
	     //if(prescalej[7]<0) hweight=hweight*(-prescalej[7]);
	     hist_leading_pt_all_Jet500->Fill(leading_pt,hweight);
	     hist_leading_eta_all_Jet500->Fill(leading_eta,hweight);
	   }
	   if (hltPassj[8] && l1cut[8] && hltcut[8]){
	     //if(prescalej[7]>0) hweight=hweight*prescalej[7];
	     //if(prescalej[7]<0) hweight=hweight*(-prescalej[7]);
	     hist_leading_pt_emulated_Jet500->Fill(leading_pt,hweight);
	     hist_leading_eta_emulated_Jet500->Fill(leading_eta,hweight);
	   } 
	 }
	 
	 hweight=1.;
       
	 for(unsigned int j=0; j<Event->nPFJetsCHS(); j++) {
	   if(Event->pfjetchs(j).ptCor()<mMinPt) continue;
	   if(fabs(Event->pfjetchs(j).y())>mYMax) continue;
	   
	   hweight=1.;

	   //cout<<Event->pfjetchs(j).tche()<<endl;

	   //parsePileUpJSON2();
	   double PileUpData=-10;
	   //PileUpData=getAvgPU(int(Event->evtHdr().runNo()),int(Event->evtHdr().lumi()));
	 
	   //cout<<int(Event->evtHdr().runNo())<<" "<<int(Event->evtHdr().lumi())<<endl;
	   //cout<<getAvgPU(int(Event->evtHdr().runNo()),int(Event->evtHdr().lumi()))<<" "<<int(Event->evtHdr().runNo())<<" "<<int(Event->evtHdr().lumi())<<endl;
    
       if(hltPassj[0] || hltPassj[10]){
	   //if(hltPassj[11]){
	     if(hltPassj[0]){ if(prescalej[0]>0) hweight=hweight*prescalej[0]; if(prescalej[0]<0) hweight=hweight*(-prescalej[0]);}
	     if(hltPassj[10]){ if(prescalej[10]>0) hweight=hweight*prescalej[10]; if(prescalej[10]<0) hweight=hweight*(-prescalej[10]);}

         if(fabs(Event->pfjetchs(j).y())<=0.5 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) {pt_DETInclJet40_1bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
	     if(fabs(Event->pfjetchs(j).y())>=0.5 && fabs(Event->pfjetchs(j).y())<=1.0 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) {pt_DETInclJet40_2bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
	     if(fabs(Event->pfjetchs(j).y())>=1.0 && fabs(Event->pfjetchs(j).y())<=1.5 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) {pt_DETInclJet40_3bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
	     if(fabs(Event->pfjetchs(j).y())>=1.5 && fabs(Event->pfjetchs(j).y())<=2.0 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) {pt_DETInclJet40_4bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
	     if(fabs(Event->pfjetchs(j).y())>=2.0 && fabs(Event->pfjetchs(j).y())<=2.5 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) {pt_DETInclJet40_5bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
	     if(fabs(Event->pfjetchs(j).y())>=2.5 && fabs(Event->pfjetchs(j).y())<=3.0 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) {pt_DETInclJet40_6bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
	     if(fabs(Event->pfjetchs(j).y())>=3.2 && fabs(Event->pfjetchs(j).y())<=4.7){ 
             if(Event->pfjetchs(j).nemf()<0.90 && Event->pfjetchs(j).ncand()>10 && Event->pfmet().met_o_sumet() < 0.3)
                 pt_DETInclJet40_7bin->Fill(Event->pfjetchs(j).ptCor(),hweight);
	        }   
        }    
	   hweight=1.;

	   if(hltPassj[1] || hltPassj[11]){
	   //if(hltPassj[11]){
	     if(hltPassj[1]){ if(prescalej[1]>0) hweight=hweight*prescalej[1]; if(prescalej[1]<0) hweight=hweight*(-prescalej[1]);}
	     if(hltPassj[11]){ if(prescalej[11]>0) hweight=hweight*prescalej[11]; if(prescalej[11]<0) hweight=hweight*(-prescalej[11]);}

         if(fabs(Event->pfjetchs(j).y())<=0.5 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) {pt_DETInclJet60_1bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
	     if(fabs(Event->pfjetchs(j).y())>=0.5 && fabs(Event->pfjetchs(j).y())<=1.0 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) {pt_DETInclJet60_2bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
	     if(fabs(Event->pfjetchs(j).y())>=1.0 && fabs(Event->pfjetchs(j).y())<=1.5 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) {pt_DETInclJet60_3bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
	     if(fabs(Event->pfjetchs(j).y())>=1.5 && fabs(Event->pfjetchs(j).y())<=2.0 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) {pt_DETInclJet60_4bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
	     if(fabs(Event->pfjetchs(j).y())>=2.0 && fabs(Event->pfjetchs(j).y())<=2.5 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) {pt_DETInclJet60_5bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
	     if(fabs(Event->pfjetchs(j).y())>=2.5 && fabs(Event->pfjetchs(j).y())<=3.0 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) {pt_DETInclJet60_6bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
	     if(fabs(Event->pfjetchs(j).y())>=3.2 && fabs(Event->pfjetchs(j).y())<=4.7){ 
             if(Event->pfjetchs(j).nemf()<0.90 && Event->pfjetchs(j).ncand()>10 && Event->pfmet().met_o_sumet() < 0.3)
                 pt_DETInclJet60_7bin->Fill(Event->pfjetchs(j).ptCor(),hweight);
	        }   
        }    
	   hweight=1.;

	   if(hltPassj[2] || hltPassj[12]){
	   //if(hltPassj[12]){
        if(hltPassj[2]){ if(prescalej[2]>0) hweight=hweight*prescalej[2]; if(prescalej[2]<0) hweight=hweight*(-prescalej[2]);}
	    if(hltPassj[12]){ if(prescalej[12]>0) hweight=hweight*prescalej[12]; if(prescalej[12]<0) hweight=hweight*(-prescalej[12]);}

             if(fabs(Event->pfjetchs(j).y())<=0.5 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) {pt_DETInclJet80_1bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
	     if(fabs(Event->pfjetchs(j).y())>=0.5 && fabs(Event->pfjetchs(j).y())<=1.0 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) {pt_DETInclJet80_2bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
	     if(fabs(Event->pfjetchs(j).y())>=1.0 && fabs(Event->pfjetchs(j).y())<=1.5 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) {pt_DETInclJet80_3bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
	     if(fabs(Event->pfjetchs(j).y())>=1.5 && fabs(Event->pfjetchs(j).y())<=2.0 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) {pt_DETInclJet80_4bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
	     if(fabs(Event->pfjetchs(j).y())>=2.0 && fabs(Event->pfjetchs(j).y())<=2.5 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) {pt_DETInclJet80_5bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
	     if(fabs(Event->pfjetchs(j).y())>=2.5 && fabs(Event->pfjetchs(j).y())<=3.0 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) {pt_DETInclJet80_6bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
         if(fabs(Event->pfjetchs(j).y())>=3.2 && fabs(Event->pfjetchs(j).y())<=4.7){ 
             if(Event->pfjetchs(j).nemf()<0.90 && Event->pfjetchs(j).ncand()>10 && Event->pfmet().met_o_sumet() < 0.3)
                 pt_DETInclJet80_7bin->Fill(Event->pfjetchs(j).ptCor(),hweight);
	        }   
       }

	   hweight=1.;

	   if(hltPassj[3] || hltPassj[13]){
	   //if(hltPassj[13]){
        if(hltPassj[3]){ if(prescalej[3]>0) hweight=hweight*prescalej[3]; if(prescalej[3]<0) hweight=hweight*(-prescalej[3]);}
	    if(hltPassj[13]){ if(prescalej[13]>0) hweight=hweight*prescalej[13]; if(prescalej[13]<0) hweight=hweight*(-prescalej[13]);}

             if(fabs(Event->pfjetchs(j).y())<=0.5 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3 ) {pt_DETInclJet140_1bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
	     if(fabs(Event->pfjetchs(j).y())>=0.5 && fabs(Event->pfjetchs(j).y())<=1.0 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) {pt_DETInclJet140_2bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
	     if(fabs(Event->pfjetchs(j).y())>=1.0 && fabs(Event->pfjetchs(j).y())<=1.5 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) {pt_DETInclJet140_3bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
	     if(fabs(Event->pfjetchs(j).y())>=1.5 && fabs(Event->pfjetchs(j).y())<=2.0 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) {pt_DETInclJet140_4bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
	     if(fabs(Event->pfjetchs(j).y())>=2.0 && fabs(Event->pfjetchs(j).y())<=2.5 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) {pt_DETInclJet140_5bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
	     if(fabs(Event->pfjetchs(j).y())>=2.5 && fabs(Event->pfjetchs(j).y())<=3.0 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) {pt_DETInclJet140_6bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
         if(fabs(Event->pfjetchs(j).y())>=3.2 && fabs(Event->pfjetchs(j).y())<=4.7){ 
             if(Event->pfjetchs(j).nemf()<0.90 && Event->pfjetchs(j).ncand()>10 && Event->pfmet().met_o_sumet() < 0.3)
                 pt_DETInclJet140_7bin->Fill(Event->pfjetchs(j).ptCor(),hweight);
	        }         
       }

	   hweight=1.;

	   if(hltPassj[4] || hltPassj[14]){
	   //if(hltPassj[14]){
        if(hltPassj[4]){ if(prescalej[4]>0) hweight=hweight*prescalej[4]; if(prescalej[4]<0) hweight=hweight*(-prescalej[4]);}
	    if(hltPassj[14]){ if(prescalej[14]>0) hweight=hweight*prescalej[14]; if(prescalej[14]<0) hweight=hweight*(-prescalej[14]);}

             if(fabs(Event->pfjetchs(j).y())<=0.5 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) {pt_DETInclJet200_1bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
	     if(fabs(Event->pfjetchs(j).y())>=0.5 && fabs(Event->pfjetchs(j).y())<=1.0 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) {pt_DETInclJet200_2bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
	     if(fabs(Event->pfjetchs(j).y())>=1.0 && fabs(Event->pfjetchs(j).y())<=1.5 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) {pt_DETInclJet200_3bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
	     if(fabs(Event->pfjetchs(j).y())>=1.5 && fabs(Event->pfjetchs(j).y())<=2.0 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) {pt_DETInclJet200_4bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
	     if(fabs(Event->pfjetchs(j).y())>=2.0 && fabs(Event->pfjetchs(j).y())<=2.5 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) {pt_DETInclJet200_5bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
	     if(fabs(Event->pfjetchs(j).y())>=2.5 && fabs(Event->pfjetchs(j).y())<=3.0 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) {pt_DETInclJet200_6bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
         if(fabs(Event->pfjetchs(j).y())>=3.2 && fabs(Event->pfjetchs(j).y())<=4.7){ 
             if(Event->pfjetchs(j).nemf()<0.90 && Event->pfjetchs(j).ncand()>10 && Event->pfmet().met_o_sumet() < 0.3)
                 pt_DETInclJet200_7bin->Fill(Event->pfjetchs(j).ptCor(),hweight);
	        }         

        }

	   hweight=1.;

	   if(hltPassj[5] || hltPassj[15]){
	   //if(hltPassj[15]){
        if(hltPassj[5]){ if(prescalej[5]>0) hweight=hweight*prescalej[5]; if(prescalej[5]<0) hweight=hweight*(-prescalej[5]);}
	    if(hltPassj[15]){ if(prescalej[15]>0) hweight=hweight*prescalej[15]; if(prescalej[15]<0) hweight=hweight*(-prescalej[15]);}

             if(fabs(Event->pfjetchs(j).y())<=0.5 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) {pt_DETInclJet260_1bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
	     if(fabs(Event->pfjetchs(j).y())>=0.5 && fabs(Event->pfjetchs(j).y())<=1.0 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) {pt_DETInclJet260_2bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
	     if(fabs(Event->pfjetchs(j).y())>=1.0 && fabs(Event->pfjetchs(j).y())<=1.5 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) {pt_DETInclJet260_3bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
	     if(fabs(Event->pfjetchs(j).y())>=1.5 && fabs(Event->pfjetchs(j).y())<=2.0 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) {pt_DETInclJet260_4bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
	     if(fabs(Event->pfjetchs(j).y())>=2.0 && fabs(Event->pfjetchs(j).y())<=2.5 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) {pt_DETInclJet260_5bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
	     if(fabs(Event->pfjetchs(j).y())>=2.5 && fabs(Event->pfjetchs(j).y())<=3.0 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) {pt_DETInclJet260_6bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
	     if(fabs(Event->pfjetchs(j).y())>=3.2 && fabs(Event->pfjetchs(j).y())<=4.7 && Event->pfjetchs(j).tightID()) {pt_DETInclJet260_7bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
         if(fabs(Event->pfjetchs(j).y())>=3.2 && fabs(Event->pfjetchs(j).y())<=4.7){ 
             if(Event->pfjetchs(j).nemf()<0.90 && Event->pfjetchs(j).ncand()>10 && Event->pfmet().met_o_sumet() < 0.3)
                 pt_DETInclJet260_7bin->Fill(Event->pfjetchs(j).ptCor(),hweight);
	        }         
 
       }

	   hweight=1.;

	   if(hltPassj[6] || hltPassj[16]){
	   //if(hltPassj[16]){
        if(hltPassj[6]){ if(prescalej[6]>0) hweight=hweight*prescalej[6]; if(prescalej[6]<0) hweight=hweight*(-prescalej[6]);}
	    if(hltPassj[16]){ if(prescalej[16]>0) hweight=hweight*prescalej[16]; if(prescalej[16]<0) hweight=hweight*(-prescalej[16]);}


             if(fabs(Event->pfjetchs(j).y())<=0.5 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) {pt_DETInclJet320_1bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
	     if(fabs(Event->pfjetchs(j).y())>=0.5 && fabs(Event->pfjetchs(j).y())<=1.0 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) {pt_DETInclJet320_2bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
	     if(fabs(Event->pfjetchs(j).y())>=1.0 && fabs(Event->pfjetchs(j).y())<=1.5 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) {pt_DETInclJet320_3bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
	     if(fabs(Event->pfjetchs(j).y())>=1.5 && fabs(Event->pfjetchs(j).y())<=2.0 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) {pt_DETInclJet320_4bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
	     if(fabs(Event->pfjetchs(j).y())>=2.0 && fabs(Event->pfjetchs(j).y())<=2.5 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) {pt_DETInclJet320_5bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
	     if(fabs(Event->pfjetchs(j).y())>=2.5 && fabs(Event->pfjetchs(j).y())<=3.0 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) {pt_DETInclJet320_6bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
         if(fabs(Event->pfjetchs(j).y())>=3.2 && fabs(Event->pfjetchs(j).y())<=4.7){ 
             if(Event->pfjetchs(j).nemf()<0.90 && Event->pfjetchs(j).ncand()>10 && Event->pfmet().met_o_sumet() < 0.3)
                 pt_DETInclJet320_7bin->Fill(Event->pfjetchs(j).ptCor(),hweight);
	        }         

       }

	   hweight=1.;

	   if(hltPassj[7] || hltPassj[17]){
	   //if(hltPassj[17]){
        if(hltPassj[7]){ if(prescalej[7]>0) hweight=hweight*prescalej[7]; if(prescalej[7]<0) hweight=hweight*(-prescalej[7]);}
	    if(hltPassj[17]){ if(prescalej[17]>0) hweight=hweight*prescalej[17]; if(prescalej[17]<0) hweight=hweight*(-prescalej[17]);}

             if(fabs(Event->pfjetchs(j).y())<=0.5 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) {pt_DETInclJet400_1bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
	     if(fabs(Event->pfjetchs(j).y())>=0.5 && fabs(Event->pfjetchs(j).y())<=1.0 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) {pt_DETInclJet400_2bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
	     if(fabs(Event->pfjetchs(j).y())>=1.0 && fabs(Event->pfjetchs(j).y())<=1.5 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) {pt_DETInclJet400_3bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
	     if(fabs(Event->pfjetchs(j).y())>=1.5 && fabs(Event->pfjetchs(j).y())<=2.0 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) {pt_DETInclJet400_4bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
	     if(fabs(Event->pfjetchs(j).y())>=2.0 && fabs(Event->pfjetchs(j).y())<=2.5 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) {pt_DETInclJet400_5bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
	     if(fabs(Event->pfjetchs(j).y())>=2.5 && fabs(Event->pfjetchs(j).y())<=3.0 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) {pt_DETInclJet400_6bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
         if(fabs(Event->pfjetchs(j).y())>=3.2 && fabs(Event->pfjetchs(j).y())<=4.7){ 
             if(Event->pfjetchs(j).nemf()<0.90 && Event->pfjetchs(j).ncand()>10 && Event->pfmet().met_o_sumet() < 0.3)
                 pt_DETInclJet400_7bin->Fill(Event->pfjetchs(j).ptCor(),hweight);
	        }         

        }

	   hweight=1.;

	   if(hltPassj[8] || hltPassj[18]){
	   //if(hltPassj[18]){
        if(hltPassj[8]){ if(prescalej[8]>0) hweight=hweight*prescalej[8]; if(prescalej[8]<0) hweight=hweight*(-prescalej[8]);}
	    if(hltPassj[18]){ if(prescalej[18]>0) hweight=hweight*prescalej[18]; if(prescalej[18]<0) hweight=hweight*(-prescalej[18]);}

             if(fabs(Event->pfjetchs(j).y())<=0.5 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) {pt_DETInclJet450_1bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
	     if(fabs(Event->pfjetchs(j).y())>=0.5 && fabs(Event->pfjetchs(j).y())<=1.0 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) {pt_DETInclJet450_2bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
	     if(fabs(Event->pfjetchs(j).y())>=1.0 && fabs(Event->pfjetchs(j).y())<=1.5 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) {pt_DETInclJet450_3bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
	     if(fabs(Event->pfjetchs(j).y())>=1.5 && fabs(Event->pfjetchs(j).y())<=2.0 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) {pt_DETInclJet450_4bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
	     if(fabs(Event->pfjetchs(j).y())>=2.0 && fabs(Event->pfjetchs(j).y())<=2.5 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) {pt_DETInclJet450_5bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
	     if(fabs(Event->pfjetchs(j).y())>=2.5 && fabs(Event->pfjetchs(j).y())<=3.0 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) {pt_DETInclJet450_6bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
         if(fabs(Event->pfjetchs(j).y())>=3.2 && fabs(Event->pfjetchs(j).y())<=4.7){ 
             if(Event->pfjetchs(j).nemf()<0.90 && Event->pfjetchs(j).ncand()>10 && Event->pfmet().met_o_sumet() < 0.3)
                 pt_DETInclJet450_7bin->Fill(Event->pfjetchs(j).ptCor(),hweight);
	        }         
 
       
       }
	 
	   hweight=1.;

	   if(Event->pfjetchs(j).ptCor()>=mMinPt){
	   
	     if(Event->pfjetchs(0).ptCor()>=mMinPt && Event->pfjetchs(0).ptCor()<74 && (hltPassj[0] || hltPassj[10])){
	     	       
	       if(hltPassj[0]){ if(prescalej[0]>0) hweight=hweight*prescalej[0]; if(prescalej[0]<0) hweight=hweight*(-prescalej[0]);}
	       if(hltPassj[10]){ if(prescalej[10]>0) hweight=hweight*prescalej[10]; if(prescalej[10]<0) hweight=hweight*(-prescalej[10]);}

	       if((fabs(Event->pfjetchs(j).y())<=3.0 && Event->pfjetchs(j).tightID())||(fabs(Event->pfjetchs(j).y())>3.2 && Event->pfjetchs(j).nemf()<0.90 && Event->pfjetchs(j).ncand()>10)){
		 //DetJets++;
		 //DETjet_ok[j]=1;       

		 pt0_DETInclJet->Fill(Event->pfjetchs(j).ptCor(),hweight);
		 pt0_DETInclJetUncor->Fill(Event->pfjetchs(j).pt(),hweight);
		 y0_DETInclJet->Fill(Event->pfjetchs(j).y(),hweight);
		 phi0_DETInclJet->Fill(Event->pfjetchs(j).phi(),hweight);
		 
		 //new histograms
		 Chargedhf0_DETJet->Fill(Event->pfjetchs(j).chf(),hweight);//charged hadron energy fraction
		 Neutralhf0_DETJet->Fill(Event->pfjetchs(j).nhf(),hweight);//neutral hadron energy fraction
		 Chargedef0_DETJet->Fill(Event->pfjetchs(j).cemf(),hweight);//charged em energy fraction
		 Photonef0_DETJet->Fill(Event->pfjetchs(j).nemf(),hweight);//neutral em energy fraction
		 Hadronef0_DETJet->Fill(Event->pfjetchs(j).hf_hf(),hweight);//hadron fraction energy fraction
		 Electromagneticef0_DETJet->Fill(Event->pfjetchs(j).hf_phf(),hweight);//electromagnetic fraction energy fraction
		 Muonef0_DETJet->Fill(Event->pfjetchs(j).muf(),hweight); //muon fraction energy fraction
		 
		 ChargedhMultiplicity0_DETJet->Fill(Event->pfjetchs(j).chm(),hweight);//charged hadron multiplicity
		 NeutralhMultiplicity0_DETJet->Fill(Event->pfjetchs(j).nhm(),hweight);//neutral hadron multiplicity
		 ChargedeMultiplicity0_DETJet->Fill(Event->pfjetchs(j).elm(),hweight);//electron em multiplicity
		 PhotoneMultiplicity0_DETJet->Fill(Event->pfjetchs(j).phm(),hweight);//neutral em multiplicity
		 HadroneMultiplicity0_DETJet->Fill(Event->pfjetchs(j).hf_hm(),hweight);//hadron fraction multiplicity
		 ElectromagneticeMultiplicity0_DETJet->Fill(Event->pfjetchs(j).hf_phm(),hweight);//electromagnetic fraction multiplicity
		 MuoneMultiplicity0_DETJet->Fill(Event->pfjetchs(j).mum(),hweight); //muon fraction multiplicity
	       }

	       if(Event->pfjetchs(j).ptCor()>=10){
		 if(fabs(Event->pfjetchs(j).y())<=0.5 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3 ) {pt_DETInclJet_1bin->Fill(Event->pfjetchs(j).ptCor(),hweight); pt_DETInclJetUP_1bin->Fill(Event->pfjetchs(j).ptCor()+Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);pt_DETInclJetDOWN_1bin->Fill(Event->pfjetchs(j).ptCor()-Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);}
		 if(fabs(Event->pfjetchs(j).y())>0.5 && fabs(Event->pfjetchs(j).y())<=1.0 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) {pt_DETInclJet_2bin->Fill(Event->pfjetchs(j).ptCor(),hweight);pt_DETInclJetUP_2bin->Fill(Event->pfjetchs(j).ptCor()+Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);pt_DETInclJetDOWN_2bin->Fill(Event->pfjetchs(j).ptCor()-Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);}
		 if(fabs(Event->pfjetchs(j).y())>1.0 && fabs(Event->pfjetchs(j).y())<=1.5 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) {pt_DETInclJet_3bin->Fill(Event->pfjetchs(j).ptCor(),hweight);pt_DETInclJetUP_3bin->Fill(Event->pfjetchs(j).ptCor()+Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);pt_DETInclJetDOWN_3bin->Fill(Event->pfjetchs(j).ptCor()-Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);}
		 if(fabs(Event->pfjetchs(j).y())>1.5 && fabs(Event->pfjetchs(j).y())<=2.0 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) {pt_DETInclJet_4bin->Fill(Event->pfjetchs(j).ptCor(),hweight);pt_DETInclJetUP_4bin->Fill(Event->pfjetchs(j).ptCor()+Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);pt_DETInclJetDOWN_4bin->Fill(Event->pfjetchs(j).ptCor()-Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);}
		 if(fabs(Event->pfjetchs(j).y())>2.0 && fabs(Event->pfjetchs(j).y())<=2.5 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) {pt_DETInclJet_5bin->Fill(Event->pfjetchs(j).ptCor(),hweight);pt_DETInclJetUP_5bin->Fill(Event->pfjetchs(j).ptCor()+Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);pt_DETInclJetDOWN_5bin->Fill(Event->pfjetchs(j).ptCor()-Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);}
		 if(fabs(Event->pfjetchs(j).y())>2.5 && fabs(Event->pfjetchs(j).y())<=3.0 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) {pt_DETInclJet_6bin->Fill(Event->pfjetchs(j).ptCor(),hweight);pt_DETInclJetUP_6bin->Fill(Event->pfjetchs(j).ptCor()+Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);pt_DETInclJetDOWN_6bin->Fill(Event->pfjetchs(j).ptCor()-Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);}
		 if(fabs(Event->pfjetchs(j).y())>3.2 && fabs(Event->pfjetchs(j).y())<=4.7){
		   if(Event->pfjetchs(j).nemf()<0.90 && Event->pfjetchs(j).ncand()>10 && Event->pfmet().met_o_sumet() < 0.3 )//tight JETID forward region
		     {pt_DETInclJet_7bin->Fill(Event->pfjetchs(j).ptCor(),hweight);pt_DETInclJetUP_7bin->Fill(Event->pfjetchs(j).ptCor()+Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);pt_DETInclJetDOWN_7bin->Fill(Event->pfjetchs(j).ptCor()-Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);}
		 }
	       
		 if(fabs(Event->pfjetchs(j).y()) < 2.5){
		   CSV_DET->Fill(Event->pfjetchs(j).csv(),hweight);
		   TCHP_DET->Fill(Event->pfjetchs(j).tchp(),hweight);
		   TCHE_DET->Fill(Event->pfjetchs(j).tche(),hweight);
		   SimpleSecVertexHE_DET->Fill(Event->pfjetchs(j).simplesecvertexhe(),hweight);
		   SimpleSecVertexHP_DET->Fill(Event->pfjetchs(j).simplesecvertexhp(),hweight);
		 }

		 if(Event->pfjetchs(j).csv() > 0.679){
		   if(fabs(Event->pfjetchs(j).y())<=0.5 && Event->pfjetchs(j).tightID()) {pt_DETBTagJet_1bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
		   if(fabs(Event->pfjetchs(j).y())>0.5 && fabs(Event->pfjetchs(j).y())<=1.0 && Event->pfjetchs(j).tightID()) {pt_DETBTagJet_2bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
		   if(fabs(Event->pfjetchs(j).y())>1.0 && fabs(Event->pfjetchs(j).y())<=1.5 && Event->pfjetchs(j).tightID()) {pt_DETBTagJet_3bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
		   if(fabs(Event->pfjetchs(j).y())>1.5 && fabs(Event->pfjetchs(j).y())<=2.0 && Event->pfjetchs(j).tightID()) {pt_DETBTagJet_4bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
		   if(fabs(Event->pfjetchs(j).y())>2.0 && fabs(Event->pfjetchs(j).y())<=2.5 && Event->pfjetchs(j).tightID()) {pt_DETBTagJet_5bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
		   if(fabs(Event->pfjetchs(j).y())>2.5 && fabs(Event->pfjetchs(j).y())<=3.0 && Event->pfjetchs(j).tightID()) {pt_DETBTagJet_6bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
		   if(fabs(Event->pfjetchs(j).y())>3.2 && fabs(Event->pfjetchs(j).y())<=4.7) {
		     if(Event->pfjetchs(j).nemf()<0.90 && Event->pfjetchs(j).ncand()>10)//tight JETID forward region
		       {pt_DETBTagJet_7bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
		   }
		 }
	       }
	     }
	     
	     if(Event->pfjetchs(0).ptCor()>=74 && Event->pfjetchs(0).ptCor()<114 && (hltPassj[1] || hltPassj[11])){
	          
	       //if(Event->pfjetchs(0).ptCor()>=137 && Event->pfjetchs(0).ptCor()<200 && hltPassj[1]){
	       if(hltPassj[1]) { if(prescalej[1]>0) hweight=hweight*prescalej[1]; if(prescalej[1]<0) hweight=hweight*(-prescalej[1]);}
	       if(hltPassj[11]) { if(prescalej[11]>0) hweight=hweight*prescalej[11]; if(prescalej[11]<0) hweight=hweight*(-prescalej[11]);}

	       bool switchFillPU=1;
	       if(switchFillPU){
		 TruePileUpDataInteger->Fill(PileUpData,hweight);
		 switchFillPU=0;
	       }

	       if((fabs(Event->pfjetchs(j).y())<=3.0 && Event->pfjetchs(j).tightID())||(fabs(Event->pfjetchs(j).y())>3.2 && Event->pfjetchs(j).nemf()<0.90 && Event->pfjetchs(j).ncand()>10)){
		 DetJets++;
		 DETjet_ok[j]=1;       
		 
		 if(Event->pfjetchs(j).csv()>0.679)
		   DetJetsBTags++;

		 pt0_DETInclJet->Fill(Event->pfjetchs(j).ptCor(),hweight);
		 pt0_DETInclJetUncor->Fill(Event->pfjetchs(j).pt(),hweight);
		 y0_DETInclJet->Fill(Event->pfjetchs(j).y(),hweight);
		 phi0_DETInclJet->Fill(Event->pfjetchs(j).phi(),hweight);
		 
		 //new histograms
		 Chargedhf0_DETJet->Fill(Event->pfjetchs(j).chf(),hweight);//charged hadron energy fraction
		 Neutralhf0_DETJet->Fill(Event->pfjetchs(j).nhf(),hweight);//neutral hadron energy fraction
		 Chargedef0_DETJet->Fill(Event->pfjetchs(j).cemf(),hweight);//charged em energy fraction
		 Photonef0_DETJet->Fill(Event->pfjetchs(j).nemf(),hweight);//neutral em energy fraction
		 Hadronef0_DETJet->Fill(Event->pfjetchs(j).hf_hf(),hweight);//hadron fraction energy fraction
		 Electromagneticef0_DETJet->Fill(Event->pfjetchs(j).hf_phf(),hweight);//electromagnetic fraction energy fraction
		 Muonef0_DETJet->Fill(Event->pfjetchs(j).muf(),hweight); //muon fraction energy fraction
		 
		 ChargedhMultiplicity0_DETJet->Fill(Event->pfjetchs(j).chm(),hweight);//charged hadron multiplicity
		 NeutralhMultiplicity0_DETJet->Fill(Event->pfjetchs(j).nhm(),hweight);//neutral hadron multiplicity
		 ChargedeMultiplicity0_DETJet->Fill(Event->pfjetchs(j).elm(),hweight);//electron em multiplicity
		 PhotoneMultiplicity0_DETJet->Fill(Event->pfjetchs(j).phm(),hweight);//neutral em multiplicity
		 HadroneMultiplicity0_DETJet->Fill(Event->pfjetchs(j).hf_hm(),hweight);//hadron fraction multiplicity
		 ElectromagneticeMultiplicity0_DETJet->Fill(Event->pfjetchs(j).hf_phm(),hweight);//electromagnetic fraction multiplicity
		 MuoneMultiplicity0_DETJet->Fill(Event->pfjetchs(j).mum(),hweight); //muon fraction multiplicity
	       }
	       
	       if(Event->pfjetchs(j).ptCor()>=10){
		 if(fabs(Event->pfjetchs(j).y())<=0.5 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) {pt_DETInclJet_1bin->Fill(Event->pfjetchs(j).ptCor(),hweight); pt_DETInclJetUP_1bin->Fill(Event->pfjetchs(j).ptCor()+Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);pt_DETInclJetDOWN_1bin->Fill(Event->pfjetchs(j).ptCor()-Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);}
		 if(fabs(Event->pfjetchs(j).y())>0.5 && fabs(Event->pfjetchs(j).y())<=1.0 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) {pt_DETInclJet_2bin->Fill(Event->pfjetchs(j).ptCor(),hweight);pt_DETInclJetUP_2bin->Fill(Event->pfjetchs(j).ptCor()+Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);pt_DETInclJetDOWN_2bin->Fill(Event->pfjetchs(j).ptCor()-Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);}
		 if(fabs(Event->pfjetchs(j).y())>1.0 && fabs(Event->pfjetchs(j).y())<=1.5 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) {pt_DETInclJet_3bin->Fill(Event->pfjetchs(j).ptCor(),hweight);pt_DETInclJetUP_3bin->Fill(Event->pfjetchs(j).ptCor()+Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);pt_DETInclJetDOWN_3bin->Fill(Event->pfjetchs(j).ptCor()-Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);}
		 if(fabs(Event->pfjetchs(j).y())>1.5 && fabs(Event->pfjetchs(j).y())<=2.0 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) {pt_DETInclJet_4bin->Fill(Event->pfjetchs(j).ptCor(),hweight);pt_DETInclJetUP_4bin->Fill(Event->pfjetchs(j).ptCor()+Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);pt_DETInclJetDOWN_4bin->Fill(Event->pfjetchs(j).ptCor()-Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);}
		 if(fabs(Event->pfjetchs(j).y())>2.0 && fabs(Event->pfjetchs(j).y())<=2.5 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) {pt_DETInclJet_5bin->Fill(Event->pfjetchs(j).ptCor(),hweight);pt_DETInclJetUP_5bin->Fill(Event->pfjetchs(j).ptCor()+Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);pt_DETInclJetDOWN_5bin->Fill(Event->pfjetchs(j).ptCor()-Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);}
		 if(fabs(Event->pfjetchs(j).y())>2.5 && fabs(Event->pfjetchs(j).y())<=3.0 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) {pt_DETInclJet_6bin->Fill(Event->pfjetchs(j).ptCor(),hweight);pt_DETInclJetUP_6bin->Fill(Event->pfjetchs(j).ptCor()+Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);pt_DETInclJetDOWN_6bin->Fill(Event->pfjetchs(j).ptCor()-Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);}
		 if(fabs(Event->pfjetchs(j).y())>3.2 && fabs(Event->pfjetchs(j).y())<=4.7) {
		   if(Event->pfjetchs(j).nemf()<0.90 && Event->pfjetchs(j).ncand()>10 && Event->pfmet().met_o_sumet() < 0.3)//tight JETID forward region
		     {pt_DETInclJet_7bin->Fill(Event->pfjetchs(j).ptCor(),hweight);pt_DETInclJetUP_7bin->Fill(Event->pfjetchs(j).ptCor()+Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);pt_DETInclJetDOWN_7bin->Fill(Event->pfjetchs(j).ptCor()-Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);}
		 }

		 if(fabs(Event->pfjetchs(j).y()) < 2.5){
		   CSV_DET->Fill(Event->pfjetchs(j).csv(),hweight);
		   TCHP_DET->Fill(Event->pfjetchs(j).tchp(),hweight);
		   TCHE_DET->Fill(Event->pfjetchs(j).tche(),hweight);
		   SimpleSecVertexHE_DET->Fill(Event->pfjetchs(j).simplesecvertexhe(),hweight);
		   SimpleSecVertexHP_DET->Fill(Event->pfjetchs(j).simplesecvertexhp(),hweight);
		 }

		  if(Event->pfjetchs(j).csv() > 0.679){
		   if(fabs(Event->pfjetchs(j).y())<=0.5 && Event->pfjetchs(j).tightID()) {pt_DETBTagJet_1bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
		   if(fabs(Event->pfjetchs(j).y())>0.5 && fabs(Event->pfjetchs(j).y())<=1.0 && Event->pfjetchs(j).tightID()) {pt_DETBTagJet_2bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
		   if(fabs(Event->pfjetchs(j).y())>1.0 && fabs(Event->pfjetchs(j).y())<=1.5 && Event->pfjetchs(j).tightID()) {pt_DETBTagJet_3bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
		   if(fabs(Event->pfjetchs(j).y())>1.5 && fabs(Event->pfjetchs(j).y())<=2.0 && Event->pfjetchs(j).tightID()) {pt_DETBTagJet_4bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
		   if(fabs(Event->pfjetchs(j).y())>2.0 && fabs(Event->pfjetchs(j).y())<=2.5 && Event->pfjetchs(j).tightID()) {pt_DETBTagJet_5bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
		   if(fabs(Event->pfjetchs(j).y())>2.5 && fabs(Event->pfjetchs(j).y())<=3.0 && Event->pfjetchs(j).tightID()) {pt_DETBTagJet_6bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
		   if(fabs(Event->pfjetchs(j).y())>3.2 && fabs(Event->pfjetchs(j).y())<=4.7) {
		     if(Event->pfjetchs(j).nemf()<0.90 && Event->pfjetchs(j).ncand()>10)//tight JETID forward region
		       {pt_DETBTagJet_7bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
		   }
		 }
	       }
	     }
	   
	     if(Event->pfjetchs(0).ptCor()>=114 && Event->pfjetchs(0).ptCor() < 174 && (hltPassj[2] || hltPassj[12])){
	       if(hltPassj[2]) { if(prescalej[2]>0) hweight=hweight*prescalej[2]; if(prescalej[2]<0) hweight=hweight*(-prescalej[2]);}
	       if(hltPassj[12]) { if(prescalej[12]>0) hweight=hweight*prescalej[12]; if(prescalej[12]<0) hweight=hweight*(-prescalej[12]);}

	       bool switchFillPU=1;
               if(switchFillPU){
                 TruePileUpDataInteger->Fill(PileUpData,hweight);
                 switchFillPU=0;
               }


	       if((fabs(Event->pfjetchs(j).y())<=3.0 && Event->pfjetchs(j).tightID())||(fabs(Event->pfjetchs(j).y())>3.2 && Event->pfjetchs(j).nemf()<0.90 && Event->pfjetchs(j).ncand()>10)){
		 DetJets++;
		 DETjet_ok[j]=1;       

		 if(Event->pfjetchs(j).csv()>0.679)
		   DetJetsBTags++;

		 pt0_DETInclJet->Fill(Event->pfjetchs(j).ptCor(),hweight);
		 pt0_DETInclJetUncor->Fill(Event->pfjetchs(j).pt(),hweight);
		 y0_DETInclJet->Fill(Event->pfjetchs(j).y(),hweight);
		 phi0_DETInclJet->Fill(Event->pfjetchs(j).phi(),hweight);
		 
		 //new histograms
		 Chargedhf0_DETJet->Fill(Event->pfjetchs(j).chf(),hweight);//charged hadron energy fraction
		 Neutralhf0_DETJet->Fill(Event->pfjetchs(j).nhf(),hweight);//neutral hadron energy fraction
		 Chargedef0_DETJet->Fill(Event->pfjetchs(j).cemf(),hweight);//charged em energy fraction
		 Photonef0_DETJet->Fill(Event->pfjetchs(j).nemf(),hweight);//neutral em energy fraction
		 Hadronef0_DETJet->Fill(Event->pfjetchs(j).hf_hf(),hweight);//hadron fraction energy fraction
		 Electromagneticef0_DETJet->Fill(Event->pfjetchs(j).hf_phf(),hweight);//electromagnetic fraction energy fraction
		 Muonef0_DETJet->Fill(Event->pfjetchs(j).muf(),hweight); //muon fraction energy fraction
		 
		 ChargedhMultiplicity0_DETJet->Fill(Event->pfjetchs(j).chm(),hweight);//charged hadron multiplicity
		 NeutralhMultiplicity0_DETJet->Fill(Event->pfjetchs(j).nhm(),hweight);//neutral hadron multiplicity
		 ChargedeMultiplicity0_DETJet->Fill(Event->pfjetchs(j).elm(),hweight);//electron em multiplicity
		 PhotoneMultiplicity0_DETJet->Fill(Event->pfjetchs(j).phm(),hweight);//neutral em multiplicity
		 HadroneMultiplicity0_DETJet->Fill(Event->pfjetchs(j).hf_hm(),hweight);//hadron fraction multiplicity
		 ElectromagneticeMultiplicity0_DETJet->Fill(Event->pfjetchs(j).hf_phm(),hweight);//electromagnetic fraction multiplicity
		 MuoneMultiplicity0_DETJet->Fill(Event->pfjetchs(j).mum(),hweight); //muon fraction multiplicity
	       }

	       if(Event->pfjetchs(j).ptCor()>=10){
		 if(fabs(Event->pfjetchs(j).y())<=0.5 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3 ) {pt_DETInclJet_1bin->Fill(Event->pfjetchs(j).ptCor(),hweight); pt_DETInclJetUP_1bin->Fill(Event->pfjetchs(j).ptCor()+Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);pt_DETInclJetDOWN_1bin->Fill(Event->pfjetchs(j).ptCor()-Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);}
		 if(fabs(Event->pfjetchs(j).y())>0.5 && fabs(Event->pfjetchs(j).y())<=1.0 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) {pt_DETInclJet_2bin->Fill(Event->pfjetchs(j).ptCor(),hweight);pt_DETInclJetUP_2bin->Fill(Event->pfjetchs(j).ptCor()+Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);pt_DETInclJetDOWN_2bin->Fill(Event->pfjetchs(j).ptCor()-Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);}
		 if(fabs(Event->pfjetchs(j).y())>1.0 && fabs(Event->pfjetchs(j).y())<=1.5 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) {pt_DETInclJet_3bin->Fill(Event->pfjetchs(j).ptCor(),hweight);pt_DETInclJetUP_3bin->Fill(Event->pfjetchs(j).ptCor()+Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);pt_DETInclJetDOWN_3bin->Fill(Event->pfjetchs(j).ptCor()-Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);}
		 if(fabs(Event->pfjetchs(j).y())>1.5 && fabs(Event->pfjetchs(j).y())<=2.0 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) {pt_DETInclJet_4bin->Fill(Event->pfjetchs(j).ptCor(),hweight);pt_DETInclJetUP_4bin->Fill(Event->pfjetchs(j).ptCor()+Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);pt_DETInclJetDOWN_4bin->Fill(Event->pfjetchs(j).ptCor()-Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);}
		 if(fabs(Event->pfjetchs(j).y())>2.0 && fabs(Event->pfjetchs(j).y())<=2.5 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) {pt_DETInclJet_5bin->Fill(Event->pfjetchs(j).ptCor(),hweight);pt_DETInclJetUP_5bin->Fill(Event->pfjetchs(j).ptCor()+Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);pt_DETInclJetDOWN_5bin->Fill(Event->pfjetchs(j).ptCor()-Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);}
		 if(fabs(Event->pfjetchs(j).y())>2.5 && fabs(Event->pfjetchs(j).y())<=3.0 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) {pt_DETInclJet_6bin->Fill(Event->pfjetchs(j).ptCor(),hweight);pt_DETInclJetUP_6bin->Fill(Event->pfjetchs(j).ptCor()+Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);pt_DETInclJetDOWN_6bin->Fill(Event->pfjetchs(j).ptCor()-Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);}
		 if(fabs(Event->pfjetchs(j).y())>3.2 && fabs(Event->pfjetchs(j).y())<=4.7) {
		   if(Event->pfjetchs(j).nemf()<0.90 && Event->pfjetchs(j).ncand()>10 && Event->pfmet().met_o_sumet() < 0.3)//tight JETID forward region
		     {pt_DETInclJet_7bin->Fill(Event->pfjetchs(j).ptCor(),hweight);pt_DETInclJetUP_7bin->Fill(Event->pfjetchs(j).ptCor()+Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);pt_DETInclJetDOWN_7bin->Fill(Event->pfjetchs(j).ptCor()-Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);}
		 }
		 
		 if(fabs(Event->pfjetchs(j).y()) < 2.5){
		   CSV_DET->Fill(Event->pfjetchs(j).csv(),hweight);
		   TCHP_DET->Fill(Event->pfjetchs(j).tchp(),hweight);
		   TCHE_DET->Fill(Event->pfjetchs(j).tche(),hweight);
		   SimpleSecVertexHE_DET->Fill(Event->pfjetchs(j).simplesecvertexhe(),hweight);
		   SimpleSecVertexHP_DET->Fill(Event->pfjetchs(j).simplesecvertexhp(),hweight);
		 }

		 if(Event->pfjetchs(j).csv() > 0.679){
		   if(fabs(Event->pfjetchs(j).y())<=0.5 && Event->pfjetchs(j).tightID()) {pt_DETBTagJet_1bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
		   if(fabs(Event->pfjetchs(j).y())>0.5 && fabs(Event->pfjetchs(j).y())<=1.0 && Event->pfjetchs(j).tightID()) {pt_DETBTagJet_2bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
		   if(fabs(Event->pfjetchs(j).y())>1.0 && fabs(Event->pfjetchs(j).y())<=1.5 && Event->pfjetchs(j).tightID()) {pt_DETBTagJet_3bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
		   if(fabs(Event->pfjetchs(j).y())>1.5 && fabs(Event->pfjetchs(j).y())<=2.0 && Event->pfjetchs(j).tightID()) {pt_DETBTagJet_4bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
		   if(fabs(Event->pfjetchs(j).y())>2.0 && fabs(Event->pfjetchs(j).y())<=2.5 && Event->pfjetchs(j).tightID()) {pt_DETBTagJet_5bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
		   if(fabs(Event->pfjetchs(j).y())>2.5 && fabs(Event->pfjetchs(j).y())<=3.0 && Event->pfjetchs(j).tightID()) {pt_DETBTagJet_6bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
		   if(fabs(Event->pfjetchs(j).y())>3.2 && fabs(Event->pfjetchs(j).y())<=4.7) {
		     if(Event->pfjetchs(j).nemf()<0.90 && Event->pfjetchs(j).ncand()>10)//tight JETID forward region
		       {pt_DETBTagJet_7bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
		   }
		 }
	       }    
	     }
	     
	     if(Event->pfjetchs(0).ptCor()>=174 && Event->pfjetchs(0).ptCor()<245 && (hltPassj[3] || hltPassj[13])){
	       if(hltPassj[3]) { if(prescalej[3]>0) hweight=hweight*prescalej[3]; if(prescalej[3]<0) hweight=hweight*(-prescalej[3]);}
	       if(hltPassj[13]) { if(prescalej[13]>0) hweight=hweight*prescalej[13]; if(prescalej[13]<0) hweight=hweight*(-prescalej[13]);}

	       bool switchFillPU=1;
               if(switchFillPU){
                 TruePileUpDataInteger->Fill(PileUpData,hweight);
                 switchFillPU=0;
               }

	       if((fabs(Event->pfjetchs(j).y())<=3.0 && Event->pfjetchs(j).tightID())||(fabs(Event->pfjetchs(j).y())>3.2 && Event->pfjetchs(j).nemf()<0.90 && Event->pfjetchs(j).ncand()>10)){
		 DetJets++;
		 DETjet_ok[j]=1;       

		 if(Event->pfjetchs(j).csv()>0.679)
		   DetJetsBTags++;

		 pt0_DETInclJet->Fill(Event->pfjetchs(j).ptCor(),hweight);
		 pt0_DETInclJetUncor->Fill(Event->pfjetchs(j).pt(),hweight);
		 y0_DETInclJet->Fill(Event->pfjetchs(j).y(),hweight);
		 phi0_DETInclJet->Fill(Event->pfjetchs(j).phi(),hweight);
		 
		 //new histograms
		 Chargedhf0_DETJet->Fill(Event->pfjetchs(j).chf(),hweight);//charged hadron energy fraction
		 Neutralhf0_DETJet->Fill(Event->pfjetchs(j).nhf(),hweight);//neutral hadron energy fraction
		 Chargedef0_DETJet->Fill(Event->pfjetchs(j).cemf(),hweight);//charged em energy fraction
		 Photonef0_DETJet->Fill(Event->pfjetchs(j).nemf(),hweight);//neutral em energy fraction
		 Hadronef0_DETJet->Fill(Event->pfjetchs(j).hf_hf(),hweight);//hadron fraction energy fraction
		 Electromagneticef0_DETJet->Fill(Event->pfjetchs(j).hf_phf(),hweight);//electromagnetic fraction energy fraction
		 Muonef0_DETJet->Fill(Event->pfjetchs(j).muf(),hweight); //muon fraction energy fraction
		 
		 ChargedhMultiplicity0_DETJet->Fill(Event->pfjetchs(j).chm(),hweight);//charged hadron multiplicity
		 NeutralhMultiplicity0_DETJet->Fill(Event->pfjetchs(j).nhm(),hweight);//neutral hadron multiplicity
		 ChargedeMultiplicity0_DETJet->Fill(Event->pfjetchs(j).elm(),hweight);//electron em multiplicity
		 PhotoneMultiplicity0_DETJet->Fill(Event->pfjetchs(j).phm(),hweight);//neutral em multiplicity
		 HadroneMultiplicity0_DETJet->Fill(Event->pfjetchs(j).hf_hm(),hweight);//hadron fraction multiplicity
		 ElectromagneticeMultiplicity0_DETJet->Fill(Event->pfjetchs(j).hf_phm(),hweight);//electromagnetic fraction multiplicity
		 MuoneMultiplicity0_DETJet->Fill(Event->pfjetchs(j).mum(),hweight); //muon fraction multiplicity
	       }

	       if(Event->pfjetchs(j).ptCor()>=10){
		 if(fabs(Event->pfjetchs(j).y())<=0.5 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) {pt_DETInclJet_1bin->Fill(Event->pfjetchs(j).ptCor(),hweight); pt_DETInclJetUP_1bin->Fill(Event->pfjetchs(j).ptCor()+Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);pt_DETInclJetDOWN_1bin->Fill(Event->pfjetchs(j).ptCor()-Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);}
		 if(fabs(Event->pfjetchs(j).y())>0.5 && fabs(Event->pfjetchs(j).y())<=1.0 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) {pt_DETInclJet_2bin->Fill(Event->pfjetchs(j).ptCor(),hweight);pt_DETInclJetUP_2bin->Fill(Event->pfjetchs(j).ptCor()+Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);pt_DETInclJetDOWN_2bin->Fill(Event->pfjetchs(j).ptCor()-Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);}
		 if(fabs(Event->pfjetchs(j).y())>1.0 && fabs(Event->pfjetchs(j).y())<=1.5 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) {pt_DETInclJet_3bin->Fill(Event->pfjetchs(j).ptCor(),hweight);pt_DETInclJetUP_3bin->Fill(Event->pfjetchs(j).ptCor()+Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);pt_DETInclJetDOWN_3bin->Fill(Event->pfjetchs(j).ptCor()-Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);}
		 if(fabs(Event->pfjetchs(j).y())>1.5 && fabs(Event->pfjetchs(j).y())<=2.0 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) {pt_DETInclJet_4bin->Fill(Event->pfjetchs(j).ptCor(),hweight);pt_DETInclJetUP_4bin->Fill(Event->pfjetchs(j).ptCor()+Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);pt_DETInclJetDOWN_4bin->Fill(Event->pfjetchs(j).ptCor()-Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);}
		 if(fabs(Event->pfjetchs(j).y())>2.0 && fabs(Event->pfjetchs(j).y())<=2.5 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) {pt_DETInclJet_5bin->Fill(Event->pfjetchs(j).ptCor(),hweight);pt_DETInclJetUP_5bin->Fill(Event->pfjetchs(j).ptCor()+Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);pt_DETInclJetDOWN_5bin->Fill(Event->pfjetchs(j).ptCor()-Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);}
		 if(fabs(Event->pfjetchs(j).y())>2.5 && fabs(Event->pfjetchs(j).y())<=3.0 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) {pt_DETInclJet_6bin->Fill(Event->pfjetchs(j).ptCor(),hweight);pt_DETInclJetUP_6bin->Fill(Event->pfjetchs(j).ptCor()+Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);pt_DETInclJetDOWN_6bin->Fill(Event->pfjetchs(j).ptCor()-Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);}
		 if(fabs(Event->pfjetchs(j).y())>3.2 && fabs(Event->pfjetchs(j).y())<=4.7) {
		   if(Event->pfjetchs(j).nemf()<0.90 && Event->pfjetchs(j).ncand()>10 && Event->pfmet().met_o_sumet() < 0.3)//tight JETID forward region
		     {pt_DETInclJet_7bin->Fill(Event->pfjetchs(j).ptCor(),hweight);pt_DETInclJetUP_7bin->Fill(Event->pfjetchs(j).ptCor()+Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);pt_DETInclJetDOWN_7bin->Fill(Event->pfjetchs(j).ptCor()-Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);}
		 }
		 
		 if(fabs(Event->pfjetchs(j).y()) < 2.5){
		   CSV_DET->Fill(Event->pfjetchs(j).csv(),hweight);
		   TCHP_DET->Fill(Event->pfjetchs(j).tchp(),hweight);
		   TCHE_DET->Fill(Event->pfjetchs(j).tche(),hweight);
		   SimpleSecVertexHE_DET->Fill(Event->pfjetchs(j).simplesecvertexhe(),hweight);
		   SimpleSecVertexHP_DET->Fill(Event->pfjetchs(j).simplesecvertexhp(),hweight);
		 }

		 if(Event->pfjetchs(j).csv() > 0.679){
		   if(fabs(Event->pfjetchs(j).y())<=0.5 && Event->pfjetchs(j).tightID()) {pt_DETBTagJet_1bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
		   if(fabs(Event->pfjetchs(j).y())>0.5 && fabs(Event->pfjetchs(j).y())<=1.0 && Event->pfjetchs(j).tightID()) {pt_DETBTagJet_2bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
		   if(fabs(Event->pfjetchs(j).y())>1.0 && fabs(Event->pfjetchs(j).y())<=1.5 && Event->pfjetchs(j).tightID()) {pt_DETBTagJet_3bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
		   if(fabs(Event->pfjetchs(j).y())>1.5 && fabs(Event->pfjetchs(j).y())<=2.0 && Event->pfjetchs(j).tightID()) {pt_DETBTagJet_4bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
		   if(fabs(Event->pfjetchs(j).y())>2.0 && fabs(Event->pfjetchs(j).y())<=2.5 && Event->pfjetchs(j).tightID()) {pt_DETBTagJet_5bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
		   if(fabs(Event->pfjetchs(j).y())>2.5 && fabs(Event->pfjetchs(j).y())<=3.0 && Event->pfjetchs(j).tightID()) {pt_DETBTagJet_6bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
		   if(fabs(Event->pfjetchs(j).y())>3.2 && fabs(Event->pfjetchs(j).y())<=4.7) {
		     if(Event->pfjetchs(j).nemf()<0.90 && Event->pfjetchs(j).ncand()>10)//tight JETID forward region
		       {pt_DETBTagJet_7bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
		   }
		 }
	       }
	     }
	     
	     if(Event->pfjetchs(0).ptCor()>=245 && Event->pfjetchs(0).ptCor()<300 && (hltPassj[4] || hltPassj[14])){
	       if(hltPassj[4]) { if(prescalej[4]>0) hweight=hweight*prescalej[4]; if(prescalej[4]<0) hweight=hweight*(-prescalej[4]);}
	       if(hltPassj[14]) { if(prescalej[14]>0) hweight=hweight*prescalej[14]; if(prescalej[14]<0) hweight=hweight*(-prescalej[14]);}

	       bool switchFillPU=1;
               if(switchFillPU){
                 TruePileUpDataInteger->Fill(PileUpData,hweight);
                 switchFillPU=0;
               }

	     
	       if((fabs(Event->pfjetchs(j).y())<=3.0 && Event->pfjetchs(j).tightID())||(fabs(Event->pfjetchs(j).y())>3.2 && Event->pfjetchs(j).nemf()<0.90 && Event->pfjetchs(j).ncand()>10)){
		 DetJets++;
		 DETjet_ok[j]=1;       

		 if(Event->pfjetchs(j).csv()>0.679)
		   DetJetsBTags++;

		 pt0_DETInclJet->Fill(Event->pfjetchs(j).ptCor(),hweight);
		 pt0_DETInclJetUncor->Fill(Event->pfjetchs(j).pt(),hweight);
		 y0_DETInclJet->Fill(Event->pfjetchs(j).y(),hweight);
		 phi0_DETInclJet->Fill(Event->pfjetchs(j).phi(),hweight);
		 
		 //new histograms
		 Chargedhf0_DETJet->Fill(Event->pfjetchs(j).chf(),hweight);//charged hadron energy fraction
		 Neutralhf0_DETJet->Fill(Event->pfjetchs(j).nhf(),hweight);//neutral hadron energy fraction
		 Chargedef0_DETJet->Fill(Event->pfjetchs(j).cemf(),hweight);//charged em energy fraction
		 Photonef0_DETJet->Fill(Event->pfjetchs(j).nemf(),hweight);//neutral em energy fraction
		 Hadronef0_DETJet->Fill(Event->pfjetchs(j).hf_hf(),hweight);//hadron fraction energy fraction
		 Electromagneticef0_DETJet->Fill(Event->pfjetchs(j).hf_phf(),hweight);//electromagnetic fraction energy fraction
		 Muonef0_DETJet->Fill(Event->pfjetchs(j).muf(),hweight); //muon fraction energy fraction
		 
		 ChargedhMultiplicity0_DETJet->Fill(Event->pfjetchs(j).chm(),hweight);//charged hadron multiplicity
		 NeutralhMultiplicity0_DETJet->Fill(Event->pfjetchs(j).nhm(),hweight);//neutral hadron multiplicity
		 ChargedeMultiplicity0_DETJet->Fill(Event->pfjetchs(j).elm(),hweight);//electron em multiplicity
		 PhotoneMultiplicity0_DETJet->Fill(Event->pfjetchs(j).phm(),hweight);//neutral em multiplicity
		 HadroneMultiplicity0_DETJet->Fill(Event->pfjetchs(j).hf_hm(),hweight);//hadron fraction multiplicity
		 ElectromagneticeMultiplicity0_DETJet->Fill(Event->pfjetchs(j).hf_phm(),hweight);//electromagnetic fraction multiplicity
		 MuoneMultiplicity0_DETJet->Fill(Event->pfjetchs(j).mum(),hweight); //muon fraction multiplicity
	       }
	       
	       if(Event->pfjetchs(j).ptCor()>=10){
		 if(fabs(Event->pfjetchs(j).y())<=0.5 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) {pt_DETInclJet_1bin->Fill(Event->pfjetchs(j).ptCor(),hweight); pt_DETInclJetUP_1bin->Fill(Event->pfjetchs(j).ptCor()+Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);pt_DETInclJetDOWN_1bin->Fill(Event->pfjetchs(j).ptCor()-Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);}
		 if(fabs(Event->pfjetchs(j).y())>0.5 && fabs(Event->pfjetchs(j).y())<=1.0 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) {pt_DETInclJet_2bin->Fill(Event->pfjetchs(j).ptCor(),hweight);pt_DETInclJetUP_2bin->Fill(Event->pfjetchs(j).ptCor()+Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);pt_DETInclJetDOWN_2bin->Fill(Event->pfjetchs(j).ptCor()-Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);}
		 if(fabs(Event->pfjetchs(j).y())>1.0 && fabs(Event->pfjetchs(j).y())<=1.5 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) {pt_DETInclJet_3bin->Fill(Event->pfjetchs(j).ptCor(),hweight);pt_DETInclJetUP_3bin->Fill(Event->pfjetchs(j).ptCor()+Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);pt_DETInclJetDOWN_3bin->Fill(Event->pfjetchs(j).ptCor()-Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);}
		 if(fabs(Event->pfjetchs(j).y())>1.5 && fabs(Event->pfjetchs(j).y())<=2.0 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) {pt_DETInclJet_4bin->Fill(Event->pfjetchs(j).ptCor(),hweight);pt_DETInclJetUP_4bin->Fill(Event->pfjetchs(j).ptCor()+Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);pt_DETInclJetDOWN_4bin->Fill(Event->pfjetchs(j).ptCor()-Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);}
		 if(fabs(Event->pfjetchs(j).y())>2.0 && fabs(Event->pfjetchs(j).y())<=2.5 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) {pt_DETInclJet_5bin->Fill(Event->pfjetchs(j).ptCor(),hweight);pt_DETInclJetUP_5bin->Fill(Event->pfjetchs(j).ptCor()+Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);pt_DETInclJetDOWN_5bin->Fill(Event->pfjetchs(j).ptCor()-Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);}
		 if(fabs(Event->pfjetchs(j).y())>2.5 && fabs(Event->pfjetchs(j).y())<=3.0 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) {pt_DETInclJet_6bin->Fill(Event->pfjetchs(j).ptCor(),hweight);pt_DETInclJetUP_6bin->Fill(Event->pfjetchs(j).ptCor()+Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);pt_DETInclJetDOWN_6bin->Fill(Event->pfjetchs(j).ptCor()-Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);}
		 if(fabs(Event->pfjetchs(j).y())>3.2 && fabs(Event->pfjetchs(j).y())<=4.7) {
		 if(Event->pfjetchs(j).nemf()<0.90 && Event->pfjetchs(j).ncand()>10 && Event->pfmet().met_o_sumet() < 0.3)//tight JETID forward region
		 {pt_DETInclJet_7bin->Fill(Event->pfjetchs(j).ptCor(),hweight);pt_DETInclJetUP_7bin->Fill(Event->pfjetchs(j).ptCor()+Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);pt_DETInclJetDOWN_7bin->Fill(Event->pfjetchs(j).ptCor()-Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);}
		 }

		 if(fabs(Event->pfjetchs(j).y()) < 2.5){
		   CSV_DET->Fill(Event->pfjetchs(j).csv(),hweight);
		   TCHP_DET->Fill(Event->pfjetchs(j).tchp(),hweight);
		   TCHE_DET->Fill(Event->pfjetchs(j).tche(),hweight);
		   SimpleSecVertexHE_DET->Fill(Event->pfjetchs(j).simplesecvertexhe(),hweight);
		   SimpleSecVertexHP_DET->Fill(Event->pfjetchs(j).simplesecvertexhp(),hweight);
		 }

		 if(Event->pfjetchs(j).csv() > 0.679){
		   if(fabs(Event->pfjetchs(j).y())<=0.5 && Event->pfjetchs(j).tightID()) {pt_DETBTagJet_1bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
		   if(fabs(Event->pfjetchs(j).y())>0.5 && fabs(Event->pfjetchs(j).y())<=1.0 && Event->pfjetchs(j).tightID()) {pt_DETBTagJet_2bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
		   if(fabs(Event->pfjetchs(j).y())>1.0 && fabs(Event->pfjetchs(j).y())<=1.5 && Event->pfjetchs(j).tightID()) {pt_DETBTagJet_3bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
		   if(fabs(Event->pfjetchs(j).y())>1.5 && fabs(Event->pfjetchs(j).y())<=2.0 && Event->pfjetchs(j).tightID()) {pt_DETBTagJet_4bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
		   if(fabs(Event->pfjetchs(j).y())>2.0 && fabs(Event->pfjetchs(j).y())<=2.5 && Event->pfjetchs(j).tightID()) {pt_DETBTagJet_5bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
		   if(fabs(Event->pfjetchs(j).y())>2.5 && fabs(Event->pfjetchs(j).y())<=3.0 && Event->pfjetchs(j).tightID()) {pt_DETBTagJet_6bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
		   if(fabs(Event->pfjetchs(j).y())>3.2 && fabs(Event->pfjetchs(j).y())<=4.7) {
		     if(Event->pfjetchs(j).nemf()<0.90 && Event->pfjetchs(j).ncand()>10)//tight JETID forward region
		       {pt_DETBTagJet_7bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
		   }
		 }
	       }
	     }
	     
	     if(Event->pfjetchs(0).ptCor()>=300 && Event->pfjetchs(0).ptCor()<395 && (hltPassj[5] || hltPassj[15])){
	       //if(Event->pfjetchs(0).ptCor()>=500 && Event->pfjetchs(0).ptCor()<600 && hltPassj[5]){
	       if(hltPassj[5]) { if(prescalej[5]>0) hweight=hweight*prescalej[5]; if(prescalej[5]<0) hweight=hweight*(-prescalej[5]);}
	       if(hltPassj[15]) { if(prescalej[15]>0) hweight=hweight*prescalej[15]; if(prescalej[15]<0) hweight=hweight*(-prescalej[15]);}

	       bool switchFillPU=1;
               if(switchFillPU){
                 TruePileUpDataInteger->Fill(PileUpData,hweight);
                 switchFillPU=0;
               }


	        if((fabs(Event->pfjetchs(j).y())<=3.0 && Event->pfjetchs(j).tightID())||(fabs(Event->pfjetchs(j).y())>3.2 && Event->pfjetchs(j).nemf()<0.90 && Event->pfjetchs(j).ncand()>10)){
		 DetJets++;
		 DETjet_ok[j]=1;       

		 if(Event->pfjetchs(j).csv()>0.679)
		   DetJetsBTags++;

		 pt0_DETInclJet->Fill(Event->pfjetchs(j).ptCor(),hweight);
		 pt0_DETInclJetUncor->Fill(Event->pfjetchs(j).pt(),hweight);
		 y0_DETInclJet->Fill(Event->pfjetchs(j).y(),hweight);
		 phi0_DETInclJet->Fill(Event->pfjetchs(j).phi(),hweight);
		 
		 //new histograms
		 Chargedhf0_DETJet->Fill(Event->pfjetchs(j).chf(),hweight);//charged hadron energy fraction
		 Neutralhf0_DETJet->Fill(Event->pfjetchs(j).nhf(),hweight);//neutral hadron energy fraction
		 Chargedef0_DETJet->Fill(Event->pfjetchs(j).cemf(),hweight);//charged em energy fraction
		 Photonef0_DETJet->Fill(Event->pfjetchs(j).nemf(),hweight);//neutral em energy fraction
		 Hadronef0_DETJet->Fill(Event->pfjetchs(j).hf_hf(),hweight);//hadron fraction energy fraction
		 Electromagneticef0_DETJet->Fill(Event->pfjetchs(j).hf_phf(),hweight);//electromagnetic fraction energy fraction
		 Muonef0_DETJet->Fill(Event->pfjetchs(j).muf(),hweight); //muon fraction energy fraction
		 
		 ChargedhMultiplicity0_DETJet->Fill(Event->pfjetchs(j).chm(),hweight);//charged hadron multiplicity
		 NeutralhMultiplicity0_DETJet->Fill(Event->pfjetchs(j).nhm(),hweight);//neutral hadron multiplicity
		 ChargedeMultiplicity0_DETJet->Fill(Event->pfjetchs(j).elm(),hweight);//electron em multiplicity
		 PhotoneMultiplicity0_DETJet->Fill(Event->pfjetchs(j).phm(),hweight);//neutral em multiplicity
		 HadroneMultiplicity0_DETJet->Fill(Event->pfjetchs(j).hf_hm(),hweight);//hadron fraction multiplicity
		 ElectromagneticeMultiplicity0_DETJet->Fill(Event->pfjetchs(j).hf_phm(),hweight);//electromagnetic fraction multiplicity
		 MuoneMultiplicity0_DETJet->Fill(Event->pfjetchs(j).mum(),hweight); //muon fraction multiplicity
		}

	       if(Event->pfjetchs(j).ptCor()>=10){
		 if(fabs(Event->pfjetchs(j).y())<=0.5 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) {pt_DETInclJet_1bin->Fill(Event->pfjetchs(j).ptCor(),hweight); pt_DETInclJetUP_1bin->Fill(Event->pfjetchs(j).ptCor()+Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);pt_DETInclJetDOWN_1bin->Fill(Event->pfjetchs(j).ptCor()-Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);}
		 if(fabs(Event->pfjetchs(j).y())>0.5 && fabs(Event->pfjetchs(j).y())<=1.0 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) {pt_DETInclJet_2bin->Fill(Event->pfjetchs(j).ptCor(),hweight);pt_DETInclJetUP_2bin->Fill(Event->pfjetchs(j).ptCor()+Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);pt_DETInclJetDOWN_2bin->Fill(Event->pfjetchs(j).ptCor()-Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);}
		 if(fabs(Event->pfjetchs(j).y())>1.0 && fabs(Event->pfjetchs(j).y())<=1.5 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) {pt_DETInclJet_3bin->Fill(Event->pfjetchs(j).ptCor(),hweight);pt_DETInclJetUP_3bin->Fill(Event->pfjetchs(j).ptCor()+Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);pt_DETInclJetDOWN_3bin->Fill(Event->pfjetchs(j).ptCor()-Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);}
		 if(fabs(Event->pfjetchs(j).y())>1.5 && fabs(Event->pfjetchs(j).y())<=2.0 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) {pt_DETInclJet_4bin->Fill(Event->pfjetchs(j).ptCor(),hweight);pt_DETInclJetUP_4bin->Fill(Event->pfjetchs(j).ptCor()+Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);pt_DETInclJetDOWN_4bin->Fill(Event->pfjetchs(j).ptCor()-Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);}
		 if(fabs(Event->pfjetchs(j).y())>2.0 && fabs(Event->pfjetchs(j).y())<=2.5 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) {pt_DETInclJet_5bin->Fill(Event->pfjetchs(j).ptCor(),hweight);pt_DETInclJetUP_5bin->Fill(Event->pfjetchs(j).ptCor()+Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);pt_DETInclJetDOWN_5bin->Fill(Event->pfjetchs(j).ptCor()-Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);}
		 if(fabs(Event->pfjetchs(j).y())>2.5 && fabs(Event->pfjetchs(j).y())<=3.0 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) {pt_DETInclJet_6bin->Fill(Event->pfjetchs(j).ptCor(),hweight);pt_DETInclJetUP_6bin->Fill(Event->pfjetchs(j).ptCor()+Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);pt_DETInclJetDOWN_6bin->Fill(Event->pfjetchs(j).ptCor()-Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);}
		 if(fabs(Event->pfjetchs(j).y())>3.2 && fabs(Event->pfjetchs(j).y())<=4.7) {
		 if(Event->pfjetchs(j).nemf()<0.90 && Event->pfjetchs(j).ncand()>10 && Event->pfmet().met_o_sumet() < 0.3)//tight JETID forward region
		 {pt_DETInclJet_7bin->Fill(Event->pfjetchs(j).ptCor(),hweight);pt_DETInclJetUP_7bin->Fill(Event->pfjetchs(j).ptCor()+Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);pt_DETInclJetDOWN_7bin->Fill(Event->pfjetchs(j).ptCor()-Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);}
		 }

		 if(fabs(Event->pfjetchs(j).y()) < 2.5){
		   CSV_DET->Fill(Event->pfjetchs(j).csv(),hweight);
		   TCHP_DET->Fill(Event->pfjetchs(j).tchp(),hweight);
		   TCHE_DET->Fill(Event->pfjetchs(j).tche(),hweight);
		   SimpleSecVertexHE_DET->Fill(Event->pfjetchs(j).simplesecvertexhe(),hweight);
		   SimpleSecVertexHP_DET->Fill(Event->pfjetchs(j).simplesecvertexhp(),hweight);
		 }

		 if(Event->pfjetchs(j).csv() > 0.679){
		   if(fabs(Event->pfjetchs(j).y())<=0.5 && Event->pfjetchs(j).tightID()) {pt_DETBTagJet_1bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
		   if(fabs(Event->pfjetchs(j).y())>0.5 && fabs(Event->pfjetchs(j).y())<=1.0 && Event->pfjetchs(j).tightID()) {pt_DETBTagJet_2bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
		   if(fabs(Event->pfjetchs(j).y())>1.0 && fabs(Event->pfjetchs(j).y())<=1.5 && Event->pfjetchs(j).tightID()) {pt_DETBTagJet_3bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
		   if(fabs(Event->pfjetchs(j).y())>1.5 && fabs(Event->pfjetchs(j).y())<=2.0 && Event->pfjetchs(j).tightID()) {pt_DETBTagJet_4bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
		   if(fabs(Event->pfjetchs(j).y())>2.0 && fabs(Event->pfjetchs(j).y())<=2.5 && Event->pfjetchs(j).tightID()) {pt_DETBTagJet_5bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
		   if(fabs(Event->pfjetchs(j).y())>2.5 && fabs(Event->pfjetchs(j).y())<=3.0 && Event->pfjetchs(j).tightID()) {pt_DETBTagJet_6bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
		   if(fabs(Event->pfjetchs(j).y())>3.2 && fabs(Event->pfjetchs(j).y())<=4.7) {
		     if(Event->pfjetchs(j).nemf()<0.90 && Event->pfjetchs(j).ncand()>10)//tight JETID forward region
		       {pt_DETBTagJet_7bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
		   }
		 }
	       }
	     }
	     
	     if(Event->pfjetchs(0).ptCor()>=395 && Event->pfjetchs(0).ptCor()<468 && (hltPassj[6] || hltPassj[16])){

	       if(hltPassj[6]) { if(prescalej[6]>0) hweight=hweight*prescalej[6]; if(prescalej[6]<0) hweight=hweight*(-prescalej[6]);}
	       if(hltPassj[16]) { if(prescalej[16]>0) hweight=hweight*prescalej[16]; if(prescalej[16]<0) hweight=hweight*(-prescalej[16]);}

	       bool switchFillPU=1;
               if(switchFillPU){
                 TruePileUpDataInteger->Fill(PileUpData,hweight);
                 switchFillPU=0;
               }


	        if((fabs(Event->pfjetchs(j).y())<=3.0 && Event->pfjetchs(j).tightID())||(fabs(Event->pfjetchs(j).y())>3.2 && Event->pfjetchs(j).nemf()<0.90 && Event->pfjetchs(j).ncand()>10)){
		 DetJets++;
		 DETjet_ok[j]=1;       

		 if(Event->pfjetchs(j).csv()>0.679)
		   DetJetsBTags++;

		 pt0_DETInclJet->Fill(Event->pfjetchs(j).ptCor(),hweight);
		 pt0_DETInclJetUncor->Fill(Event->pfjetchs(j).pt(),hweight);
		 y0_DETInclJet->Fill(Event->pfjetchs(j).y(),hweight);
		 phi0_DETInclJet->Fill(Event->pfjetchs(j).phi(),hweight);
		 
		 //new histograms
		 Chargedhf0_DETJet->Fill(Event->pfjetchs(j).chf(),hweight);//charged hadron energy fraction
		 Neutralhf0_DETJet->Fill(Event->pfjetchs(j).nhf(),hweight);//neutral hadron energy fraction
		 Chargedef0_DETJet->Fill(Event->pfjetchs(j).cemf(),hweight);//charged em energy fraction
		 Photonef0_DETJet->Fill(Event->pfjetchs(j).nemf(),hweight);//neutral em energy fraction
		 Hadronef0_DETJet->Fill(Event->pfjetchs(j).hf_hf(),hweight);//hadron fraction energy fraction
		 Electromagneticef0_DETJet->Fill(Event->pfjetchs(j).hf_phf(),hweight);//electromagnetic fraction energy fraction
		 Muonef0_DETJet->Fill(Event->pfjetchs(j).muf(),hweight); //muon fraction energy fraction
		 
		 ChargedhMultiplicity0_DETJet->Fill(Event->pfjetchs(j).chm(),hweight);//charged hadron multiplicity
		 NeutralhMultiplicity0_DETJet->Fill(Event->pfjetchs(j).nhm(),hweight);//neutral hadron multiplicity
		 ChargedeMultiplicity0_DETJet->Fill(Event->pfjetchs(j).elm(),hweight);//electron em multiplicity
		 PhotoneMultiplicity0_DETJet->Fill(Event->pfjetchs(j).phm(),hweight);//neutral em multiplicity
		 HadroneMultiplicity0_DETJet->Fill(Event->pfjetchs(j).hf_hm(),hweight);//hadron fraction multiplicity
		 ElectromagneticeMultiplicity0_DETJet->Fill(Event->pfjetchs(j).hf_phm(),hweight);//electromagnetic fraction multiplicity
		 MuoneMultiplicity0_DETJet->Fill(Event->pfjetchs(j).mum(),hweight); //muon fraction multiplicity
	       }
	       
	       if(Event->pfjetchs(j).ptCor()>=10){
		 if(fabs(Event->pfjetchs(j).y())<=0.5 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) {pt_DETInclJet_1bin->Fill(Event->pfjetchs(j).ptCor(),hweight); pt_DETInclJetUP_1bin->Fill(Event->pfjetchs(j).ptCor()+Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);pt_DETInclJetDOWN_1bin->Fill(Event->pfjetchs(j).ptCor()-Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);}
		 if(fabs(Event->pfjetchs(j).y())>0.5 && fabs(Event->pfjetchs(j).y())<=1.0 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) {pt_DETInclJet_2bin->Fill(Event->pfjetchs(j).ptCor(),hweight);pt_DETInclJetUP_2bin->Fill(Event->pfjetchs(j).ptCor()+Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);pt_DETInclJetDOWN_2bin->Fill(Event->pfjetchs(j).ptCor()-Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);}
		 if(fabs(Event->pfjetchs(j).y())>1.0 && fabs(Event->pfjetchs(j).y())<=1.5 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) {pt_DETInclJet_3bin->Fill(Event->pfjetchs(j).ptCor(),hweight);pt_DETInclJetUP_3bin->Fill(Event->pfjetchs(j).ptCor()+Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);pt_DETInclJetDOWN_3bin->Fill(Event->pfjetchs(j).ptCor()-Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);}
		 if(fabs(Event->pfjetchs(j).y())>1.5 && fabs(Event->pfjetchs(j).y())<=2.0 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) {pt_DETInclJet_4bin->Fill(Event->pfjetchs(j).ptCor(),hweight);pt_DETInclJetUP_4bin->Fill(Event->pfjetchs(j).ptCor()+Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);pt_DETInclJetDOWN_4bin->Fill(Event->pfjetchs(j).ptCor()-Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);}
		 if(fabs(Event->pfjetchs(j).y())>2.0 && fabs(Event->pfjetchs(j).y())<=2.5 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) {pt_DETInclJet_5bin->Fill(Event->pfjetchs(j).ptCor(),hweight);pt_DETInclJetUP_5bin->Fill(Event->pfjetchs(j).ptCor()+Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);pt_DETInclJetDOWN_5bin->Fill(Event->pfjetchs(j).ptCor()-Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);}
		 if(fabs(Event->pfjetchs(j).y())>2.5 && fabs(Event->pfjetchs(j).y())<=3.0 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) {pt_DETInclJet_6bin->Fill(Event->pfjetchs(j).ptCor(),hweight);pt_DETInclJetUP_6bin->Fill(Event->pfjetchs(j).ptCor()+Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);pt_DETInclJetDOWN_6bin->Fill(Event->pfjetchs(j).ptCor()-Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);}
		 if(fabs(Event->pfjetchs(j).y())>3.2 && fabs(Event->pfjetchs(j).y())<=4.7) {
		 if(Event->pfjetchs(j).nemf()<0.90 && Event->pfjetchs(j).ncand()>10 && Event->pfmet().met_o_sumet() < 0.3)//tight JETID forward region
		 {pt_DETInclJet_7bin->Fill(Event->pfjetchs(j).ptCor(),hweight);pt_DETInclJetUP_7bin->Fill(Event->pfjetchs(j).ptCor()+Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);pt_DETInclJetDOWN_7bin->Fill(Event->pfjetchs(j).ptCor()-Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);}
		 }

		 if(fabs(Event->pfjetchs(j).y()) < 2.5){
		   CSV_DET->Fill(Event->pfjetchs(j).csv(),hweight);
		   TCHP_DET->Fill(Event->pfjetchs(j).tchp(),hweight);
		   TCHE_DET->Fill(Event->pfjetchs(j).tche(),hweight);
		   SimpleSecVertexHE_DET->Fill(Event->pfjetchs(j).simplesecvertexhe(),hweight);
		   SimpleSecVertexHP_DET->Fill(Event->pfjetchs(j).simplesecvertexhp(),hweight);
		 }

		 if(Event->pfjetchs(j).csv() > 0.679){
		   if(fabs(Event->pfjetchs(j).y())<=0.5 && Event->pfjetchs(j).tightID()) {pt_DETBTagJet_1bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
		   if(fabs(Event->pfjetchs(j).y())>0.5 && fabs(Event->pfjetchs(j).y())<=1.0 && Event->pfjetchs(j).tightID()) {pt_DETBTagJet_2bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
		   if(fabs(Event->pfjetchs(j).y())>1.0 && fabs(Event->pfjetchs(j).y())<=1.5 && Event->pfjetchs(j).tightID()) {pt_DETBTagJet_3bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
		   if(fabs(Event->pfjetchs(j).y())>1.5 && fabs(Event->pfjetchs(j).y())<=2.0 && Event->pfjetchs(j).tightID()) {pt_DETBTagJet_4bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
		   if(fabs(Event->pfjetchs(j).y())>2.0 && fabs(Event->pfjetchs(j).y())<=2.5 && Event->pfjetchs(j).tightID()) {pt_DETBTagJet_5bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
		   if(fabs(Event->pfjetchs(j).y())>2.5 && fabs(Event->pfjetchs(j).y())<=3.0 && Event->pfjetchs(j).tightID()) {pt_DETBTagJet_6bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
		   if(fabs(Event->pfjetchs(j).y())>3.2 && fabs(Event->pfjetchs(j).y())<=4.7) {
		     if(Event->pfjetchs(j).nemf()<0.90 && Event->pfjetchs(j).ncand()>10)//tight JETID forward region
		       {pt_DETBTagJet_7bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
		   }
		 }
	       }
	     }
	     
	     if(Event->pfjetchs(0).ptCor()>=468 && Event->pfjetchs(0).ptCor()<507 && (hltPassj[7] || hltPassj[17])){
	       //if(Event->pfjetchs(0).ptCor()>=700 && Event->pfjetchs(0).ptCor()<800 && hltPassj[7]){
	       if(hltPassj[7]) { if(prescalej[7]>0) hweight=hweight*prescalej[7]; if(prescalej[7]<0) hweight=hweight*(-prescalej[7]);}
	       if(hltPassj[17]) { if(prescalej[17]>0) hweight=hweight*prescalej[17]; if(prescalej[17]<0) hweight=hweight*(-prescalej[17]);}

	       bool switchFillPU=1;
               if(switchFillPU){
                 TruePileUpDataInteger->Fill(PileUpData,hweight);
                 switchFillPU=0;
               }


	        if((fabs(Event->pfjetchs(j).y())<=3.0 && Event->pfjetchs(j).tightID())||(fabs(Event->pfjetchs(j).y())>3.2 && Event->pfjetchs(j).nemf()<0.90 && Event->pfjetchs(j).ncand()>10)){
		 DetJets++;
		 DETjet_ok[j]=1;       

		 if(Event->pfjetchs(j).csv()>0.679)
		   DetJetsBTags++;

		 pt0_DETInclJet->Fill(Event->pfjetchs(j).ptCor(),hweight);
		 pt0_DETInclJetUncor->Fill(Event->pfjetchs(j).pt(),hweight);
		 y0_DETInclJet->Fill(Event->pfjetchs(j).y(),hweight);
		 phi0_DETInclJet->Fill(Event->pfjetchs(j).phi(),hweight);
		 
		 //new histograms
		 Chargedhf0_DETJet->Fill(Event->pfjetchs(j).chf(),hweight);//charged hadron energy fraction
		 Neutralhf0_DETJet->Fill(Event->pfjetchs(j).nhf(),hweight);//neutral hadron energy fraction
		 Chargedef0_DETJet->Fill(Event->pfjetchs(j).cemf(),hweight);//charged em energy fraction
		 Photonef0_DETJet->Fill(Event->pfjetchs(j).nemf(),hweight);//neutral em energy fraction
		 Hadronef0_DETJet->Fill(Event->pfjetchs(j).hf_hf(),hweight);//hadron fraction energy fraction
		 Electromagneticef0_DETJet->Fill(Event->pfjetchs(j).hf_phf(),hweight);//electromagnetic fraction energy fraction
		 Muonef0_DETJet->Fill(Event->pfjetchs(j).muf(),hweight); //muon fraction energy fraction
		 
		 ChargedhMultiplicity0_DETJet->Fill(Event->pfjetchs(j).chm(),hweight);//charged hadron multiplicity
		 NeutralhMultiplicity0_DETJet->Fill(Event->pfjetchs(j).nhm(),hweight);//neutral hadron multiplicity
		 ChargedeMultiplicity0_DETJet->Fill(Event->pfjetchs(j).elm(),hweight);//electron em multiplicity
		 PhotoneMultiplicity0_DETJet->Fill(Event->pfjetchs(j).phm(),hweight);//neutral em multiplicity
		 HadroneMultiplicity0_DETJet->Fill(Event->pfjetchs(j).hf_hm(),hweight);//hadron fraction multiplicity
		 ElectromagneticeMultiplicity0_DETJet->Fill(Event->pfjetchs(j).hf_phm(),hweight);//electromagnetic fraction multiplicity
		 MuoneMultiplicity0_DETJet->Fill(Event->pfjetchs(j).mum(),hweight); //muon fraction multiplicity
	       }

	       if(Event->pfjetchs(j).ptCor()>=10){
		 if(fabs(Event->pfjetchs(j).y())<=0.5 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) {pt_DETInclJet_1bin->Fill(Event->pfjetchs(j).ptCor(),hweight); pt_DETInclJetUP_1bin->Fill(Event->pfjetchs(j).ptCor()+Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);pt_DETInclJetDOWN_1bin->Fill(Event->pfjetchs(j).ptCor()-Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);}
		 if(fabs(Event->pfjetchs(j).y())>0.5 && fabs(Event->pfjetchs(j).y())<=1.0 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) {pt_DETInclJet_2bin->Fill(Event->pfjetchs(j).ptCor(),hweight);pt_DETInclJetUP_2bin->Fill(Event->pfjetchs(j).ptCor()+Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);pt_DETInclJetDOWN_2bin->Fill(Event->pfjetchs(j).ptCor()-Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);}
		 if(fabs(Event->pfjetchs(j).y())>1.0 && fabs(Event->pfjetchs(j).y())<=1.5 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) {pt_DETInclJet_3bin->Fill(Event->pfjetchs(j).ptCor(),hweight);pt_DETInclJetUP_3bin->Fill(Event->pfjetchs(j).ptCor()+Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);pt_DETInclJetDOWN_3bin->Fill(Event->pfjetchs(j).ptCor()-Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);}
		 if(fabs(Event->pfjetchs(j).y())>1.5 && fabs(Event->pfjetchs(j).y())<=2.0 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) {pt_DETInclJet_4bin->Fill(Event->pfjetchs(j).ptCor(),hweight);pt_DETInclJetUP_4bin->Fill(Event->pfjetchs(j).ptCor()+Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);pt_DETInclJetDOWN_4bin->Fill(Event->pfjetchs(j).ptCor()-Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);}
		 if(fabs(Event->pfjetchs(j).y())>2.0 && fabs(Event->pfjetchs(j).y())<=2.5 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) {pt_DETInclJet_5bin->Fill(Event->pfjetchs(j).ptCor(),hweight);pt_DETInclJetUP_5bin->Fill(Event->pfjetchs(j).ptCor()+Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);pt_DETInclJetDOWN_5bin->Fill(Event->pfjetchs(j).ptCor()-Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);}
		 if(fabs(Event->pfjetchs(j).y())>2.5 && fabs(Event->pfjetchs(j).y())<=3.0 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) {pt_DETInclJet_6bin->Fill(Event->pfjetchs(j).ptCor(),hweight);pt_DETInclJetUP_6bin->Fill(Event->pfjetchs(j).ptCor()+Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);pt_DETInclJetDOWN_6bin->Fill(Event->pfjetchs(j).ptCor()-Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);}
		 if(fabs(Event->pfjetchs(j).y())>3.2 && fabs(Event->pfjetchs(j).y())<=4.7) {
		 if(Event->pfjetchs(j).nemf()<0.90 && Event->pfjetchs(j).ncand()>10 && Event->pfmet().met_o_sumet() < 0.3)//tight JETID forward region
		 {pt_DETInclJet_7bin->Fill(Event->pfjetchs(j).ptCor(),hweight);pt_DETInclJetUP_7bin->Fill(Event->pfjetchs(j).ptCor()+Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);pt_DETInclJetDOWN_7bin->Fill(Event->pfjetchs(j).ptCor()-Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);}
		 }

		 if(fabs(Event->pfjetchs(j).y()) < 2.5){
		   CSV_DET->Fill(Event->pfjetchs(j).csv(),hweight);
		   TCHP_DET->Fill(Event->pfjetchs(j).tchp(),hweight);
		   TCHE_DET->Fill(Event->pfjetchs(j).tche(),hweight);
		   SimpleSecVertexHE_DET->Fill(Event->pfjetchs(j).simplesecvertexhe(),hweight);
		   SimpleSecVertexHP_DET->Fill(Event->pfjetchs(j).simplesecvertexhp(),hweight);
		 }

		 if(Event->pfjetchs(j).csv() > 0.679){
		   if(fabs(Event->pfjetchs(j).y())<=0.5 && Event->pfjetchs(j).tightID()) {pt_DETBTagJet_1bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
		   if(fabs(Event->pfjetchs(j).y())>0.5 && fabs(Event->pfjetchs(j).y())<=1.0 && Event->pfjetchs(j).tightID()) {pt_DETBTagJet_2bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
		   if(fabs(Event->pfjetchs(j).y())>1.0 && fabs(Event->pfjetchs(j).y())<=1.5 && Event->pfjetchs(j).tightID()) {pt_DETBTagJet_3bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
		   if(fabs(Event->pfjetchs(j).y())>1.5 && fabs(Event->pfjetchs(j).y())<=2.0 && Event->pfjetchs(j).tightID()) {pt_DETBTagJet_4bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
		   if(fabs(Event->pfjetchs(j).y())>2.0 && fabs(Event->pfjetchs(j).y())<=2.5 && Event->pfjetchs(j).tightID()) {pt_DETBTagJet_5bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
		   if(fabs(Event->pfjetchs(j).y())>2.5 && fabs(Event->pfjetchs(j).y())<=3.0 && Event->pfjetchs(j).tightID()) {pt_DETBTagJet_6bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
		   if(fabs(Event->pfjetchs(j).y())>3.2 && fabs(Event->pfjetchs(j).y())<=4.7) {
		     if(Event->pfjetchs(j).nemf()<0.90 && Event->pfjetchs(j).ncand()>10)//tight JETID forward region
		       {pt_DETBTagJet_7bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
		   }
		 }
	       }
	     }

	     if(Event->pfjetchs(0).ptCor()>=507 && (hltPassj[8] || hltPassj[18])){
	       //if(Event->pfjetchs(0).ptCor()>=800 && hltPassj[8]){
	       if(hltPassj[8]){ if(prescalej[8]>0) hweight=hweight*prescalej[8]; if(prescalej[8]<0) hweight=hweight*(-prescalej[8]);}
	       if(hltPassj[18]){ if(prescalej[18]>0) hweight=hweight*prescalej[18]; if(prescalej[18]<0) hweight=hweight*(-prescalej[18]); }
	       
	       bool switchFillPU=1;
               if(switchFillPU){
                 TruePileUpDataInteger->Fill(PileUpData,hweight);
                 switchFillPU=0;
               }


	        if((fabs(Event->pfjetchs(j).y())<=3.0 && Event->pfjetchs(j).tightID())||(fabs(Event->pfjetchs(j).y())>3.2 && Event->pfjetchs(j).nemf()<0.90 && Event->pfjetchs(j).ncand()>10)){
		 DetJets++;
		 DETjet_ok[j]=1;       

		 if(Event->pfjetchs(j).csv()>0.679)
		   DetJetsBTags++;

		 pt0_DETInclJet->Fill(Event->pfjetchs(j).ptCor(),hweight);
		 pt0_DETInclJetUncor->Fill(Event->pfjetchs(j).pt(),hweight);
		 y0_DETInclJet->Fill(Event->pfjetchs(j).y(),hweight);
		 phi0_DETInclJet->Fill(Event->pfjetchs(j).phi(),hweight);
		 
		 //new histograms
		 Chargedhf0_DETJet->Fill(Event->pfjetchs(j).chf(),hweight);//charged hadron energy fraction
		 Neutralhf0_DETJet->Fill(Event->pfjetchs(j).nhf(),hweight);//neutral hadron energy fraction
		 Chargedef0_DETJet->Fill(Event->pfjetchs(j).cemf(),hweight);//charged em energy fraction
		 Photonef0_DETJet->Fill(Event->pfjetchs(j).nemf(),hweight);//neutral em energy fraction
		 Hadronef0_DETJet->Fill(Event->pfjetchs(j).hf_hf(),hweight);//hadron fraction energy fraction
		 Electromagneticef0_DETJet->Fill(Event->pfjetchs(j).hf_phf(),hweight);//electromagnetic fraction energy fraction
		 Muonef0_DETJet->Fill(Event->pfjetchs(j).muf(),hweight); //muon fraction energy fraction
		 
		 ChargedhMultiplicity0_DETJet->Fill(Event->pfjetchs(j).chm(),hweight);//charged hadron multiplicity
		 NeutralhMultiplicity0_DETJet->Fill(Event->pfjetchs(j).nhm(),hweight);//neutral hadron multiplicity
		 ChargedeMultiplicity0_DETJet->Fill(Event->pfjetchs(j).elm(),hweight);//electron em multiplicity
		 PhotoneMultiplicity0_DETJet->Fill(Event->pfjetchs(j).phm(),hweight);//neutral em multiplicity
		 HadroneMultiplicity0_DETJet->Fill(Event->pfjetchs(j).hf_hm(),hweight);//hadron fraction multiplicity
		 ElectromagneticeMultiplicity0_DETJet->Fill(Event->pfjetchs(j).hf_phm(),hweight);//electromagnetic fraction multiplicity
		 MuoneMultiplicity0_DETJet->Fill(Event->pfjetchs(j).mum(),hweight); //muon fraction multiplicity
	       }
	       
	       if(Event->pfjetchs(j).ptCor()>=10){
		 if(fabs(Event->pfjetchs(j).y())<=0.5 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) {pt_DETInclJet_1bin->Fill(Event->pfjetchs(j).ptCor(),hweight); pt_DETInclJetUP_1bin->Fill(Event->pfjetchs(j).ptCor()+Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);pt_DETInclJetDOWN_1bin->Fill(Event->pfjetchs(j).ptCor()-Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);}
		 if(fabs(Event->pfjetchs(j).y())>0.5 && fabs(Event->pfjetchs(j).y())<=1.0 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) {pt_DETInclJet_2bin->Fill(Event->pfjetchs(j).ptCor(),hweight);pt_DETInclJetUP_2bin->Fill(Event->pfjetchs(j).ptCor()+Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);pt_DETInclJetDOWN_2bin->Fill(Event->pfjetchs(j).ptCor()-Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);}
		 if(fabs(Event->pfjetchs(j).y())>1.0 && fabs(Event->pfjetchs(j).y())<=1.5 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) {pt_DETInclJet_3bin->Fill(Event->pfjetchs(j).ptCor(),hweight);pt_DETInclJetUP_3bin->Fill(Event->pfjetchs(j).ptCor()+Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);pt_DETInclJetDOWN_3bin->Fill(Event->pfjetchs(j).ptCor()-Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);}
		 if(fabs(Event->pfjetchs(j).y())>1.5 && fabs(Event->pfjetchs(j).y())<=2.0 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) {pt_DETInclJet_4bin->Fill(Event->pfjetchs(j).ptCor(),hweight);pt_DETInclJetUP_4bin->Fill(Event->pfjetchs(j).ptCor()+Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);pt_DETInclJetDOWN_4bin->Fill(Event->pfjetchs(j).ptCor()-Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);}
		 if(fabs(Event->pfjetchs(j).y())>2.0 && fabs(Event->pfjetchs(j).y())<=2.5 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) {pt_DETInclJet_5bin->Fill(Event->pfjetchs(j).ptCor(),hweight);pt_DETInclJetUP_5bin->Fill(Event->pfjetchs(j).ptCor()+Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);pt_DETInclJetDOWN_5bin->Fill(Event->pfjetchs(j).ptCor()-Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);}
		 if(fabs(Event->pfjetchs(j).y())>2.5 && fabs(Event->pfjetchs(j).y())<=3.0 && Event->pfjetchs(j).tightID() && Event->pfmet().met_o_sumet() < 0.3) {pt_DETInclJet_6bin->Fill(Event->pfjetchs(j).ptCor(),hweight);pt_DETInclJetUP_6bin->Fill(Event->pfjetchs(j).ptCor()+Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);pt_DETInclJetDOWN_6bin->Fill(Event->pfjetchs(j).ptCor()-Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);}
		 if(fabs(Event->pfjetchs(j).y())>3.2 && fabs(Event->pfjetchs(j).y())<=4.7) {
		 if(Event->pfjetchs(j).nemf()<0.90 && Event->pfjetchs(j).ncand()>10 && Event->pfmet().met_o_sumet() < 0.3 )//tight JETID forward region
		 {pt_DETInclJet_7bin->Fill(Event->pfjetchs(j).ptCor(),hweight);pt_DETInclJetUP_7bin->Fill(Event->pfjetchs(j).ptCor()+Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);pt_DETInclJetDOWN_7bin->Fill(Event->pfjetchs(j).ptCor()-Event->pfjetchs(j).ptCor()*Event->pfjetchs(j).uncSrc(0),hweight);}
		 }
		 
		 if(fabs(Event->pfjetchs(j).y()) < 2.5){
		   CSV_DET->Fill(Event->pfjetchs(j).csv(),hweight);
		   TCHP_DET->Fill(Event->pfjetchs(j).tchp(),hweight);
		   TCHE_DET->Fill(Event->pfjetchs(j).tche(),hweight);
		   SimpleSecVertexHE_DET->Fill(Event->pfjetchs(j).simplesecvertexhe(),hweight);
		   SimpleSecVertexHP_DET->Fill(Event->pfjetchs(j).simplesecvertexhp(),hweight);
		 }

		 if(Event->pfjetchs(j).csv() > 0.679){
		   if(fabs(Event->pfjetchs(j).y())<=0.5 && Event->pfjetchs(j).tightID()) {pt_DETBTagJet_1bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
		   if(fabs(Event->pfjetchs(j).y())>0.5 && fabs(Event->pfjetchs(j).y())<=1.0 && Event->pfjetchs(j).tightID()) {pt_DETBTagJet_2bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
		   if(fabs(Event->pfjetchs(j).y())>1.0 && fabs(Event->pfjetchs(j).y())<=1.5 && Event->pfjetchs(j).tightID()) {pt_DETBTagJet_3bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
		   if(fabs(Event->pfjetchs(j).y())>1.5 && fabs(Event->pfjetchs(j).y())<=2.0 && Event->pfjetchs(j).tightID()) {pt_DETBTagJet_4bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
		   if(fabs(Event->pfjetchs(j).y())>2.0 && fabs(Event->pfjetchs(j).y())<=2.5 && Event->pfjetchs(j).tightID()) {pt_DETBTagJet_5bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
		   if(fabs(Event->pfjetchs(j).y())>2.5 && fabs(Event->pfjetchs(j).y())<=3.0 && Event->pfjetchs(j).tightID()) {pt_DETBTagJet_6bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
		   if(fabs(Event->pfjetchs(j).y())>3.2 && fabs(Event->pfjetchs(j).y())<=4.7) {
		     if(Event->pfjetchs(j).nemf()<0.90 && Event->pfjetchs(j).ncand()>10)//tight JETID forward region
		       {pt_DETBTagJet_7bin->Fill(Event->pfjetchs(j).ptCor(),hweight);}
		   }
		 }
	       }
	     }       
	   }	       
	 }
       }
       
       if(DETjet_ok[0]==1 && DETjet_ok[1]==1){

	 pt0_DETJet->Fill(Event->pfjetchs(0).ptCor(),hweight);
	 pt0_DETJetUncor->Fill(Event->pfjetchs(0).pt(),hweight);
	 pt1_DETJet->Fill(Event->pfjetchs(1).ptCor(),hweight);
	 pt1_DETJetUncor->Fill(Event->pfjetchs(1).pt(),hweight);
	 y0_DETJet->Fill(Event->pfjetchs(0).y(),hweight);
	 y1_DETJet->Fill(Event->pfjetchs(1).y(),hweight);
	 phi0_DETJet->Fill(Event->pfjetchs(0).phi(),hweight);
	 phi1_DETJet->Fill(Event->pfjetchs(1).phi(),hweight);
     
	 if(Event->pfjetchs(0).csv()>0.679 && Event->pfjetchs(1).csv()>0.679){  
	   pt0_DETBTagJet->Fill(Event->pfjetchs(0).ptCor(),hweight);
	   pt1_DETBTagJet->Fill(Event->pfjetchs(1).ptCor(),hweight);
	   y0_DETBTagJet->Fill(Event->pfjetchs(0).y(),hweight);
	   y1_DETBTagJet->Fill(Event->pfjetchs(1).y(),hweight);
	   reco::Particle::LorentzVector Lead1=Event->pfjetchs(0).p4();
	   reco::Particle::LorentzVector Lead2=Event->pfjetchs(1).p4();
	   reco::Particle::LorentzVector Pair=Lead1+Lead2;
	   //cout<<Lead1.pt()<<" "<<Lead2.pt()<<" "<<Pair.pt()<<endl;
	   ptPair_DETBTagJet->Fill(Pair.pt(),hweight);
	 }
       }
       
       if(DETjet_ok[0]==1){

	 Multiplicity_DETJet->Fill(DetJets,hweight);
	 num_of_VtxGood->Fill(Event->evtHdr().nVtxGood(),hweight);

	 Multiplicity_DETJetBTags->Fill(DetJetsBTags,hweight);
	 if(DetJetsBTags>=2){
	   Multiplicity_DETJetBTagsInclusive->Fill(DetJets,hweight);
	 }

	 if(Event->pfjetchs(0).csv()>0.679){
	   pt0_DETBTagJetInclusive->Fill(Event->pfjetchs(0).ptCor(),hweight);
	 }

	 if(fabs(Event->pfjetchs(0).y())<=3.0 && Event->pfjetchs(0).tightID()){
	   MET_DET->Fill(Event->pfmet().met(),hweight); 
	   METPhi_DET->Fill(fabs(Event->pfmet().phi()),hweight); 
	   FractionMET_DET->Fill(Event->pfmet().met_o_sumet(),hweight); 
	 }

     // MET fraction for each rap. bin

     if(fabs(Event->pfjetchs(0).y()) <= 0.5 && Event->pfjetchs(0).tightID())
       FractionMET_DET_1bin->Fill(Event->pfmet().met_o_sumet(),hweight); 
     if(fabs(Event->pfjetchs(0).y()) > 0.5 && fabs(Event->pfjetchs(0).y()) <= 1.0  && Event->pfjetchs(0).tightID())
       FractionMET_DET_2bin->Fill(Event->pfmet().met_o_sumet(),hweight); 
     if(fabs(Event->pfjetchs(0).y()) > 1.0 && fabs(Event->pfjetchs(0).y()) <= 1.5  && Event->pfjetchs(0).tightID())
       FractionMET_DET_3bin->Fill(Event->pfmet().met_o_sumet(),hweight); 
     if(fabs(Event->pfjetchs(0).y()) > 1.5 && fabs(Event->pfjetchs(0).y()) <= 2.0  && Event->pfjetchs(0).tightID())
       FractionMET_DET_4bin->Fill(Event->pfmet().met_o_sumet(),hweight); 
     if(fabs(Event->pfjetchs(0).y()) > 2.0 && fabs(Event->pfjetchs(0).y()) <= 2.5  && Event->pfjetchs(0).tightID())
       FractionMET_DET_5bin->Fill(Event->pfmet().met_o_sumet(),hweight); 
     if(fabs(Event->pfjetchs(0).y()) > 2.5 && fabs(Event->pfjetchs(0).y()) <= 3.0  && Event->pfjetchs(0).tightID())
       FractionMET_DET_6bin->Fill(Event->pfmet().met_o_sumet(),hweight); 





	 if(fabs(Event->pfjetchs(0).y())>3.2 && Event->pfjetchs(0).nemf()<0.90 && Event->pfjetchs(0).ncand()>10){
	   MET_DET->Fill(Event->pfmet().met(),hweight); 
	   METPhi_DET->Fill(fabs(Event->pfmet().phi()),hweight); 
	   FractionMET_DET->Fill(Event->pfmet().met_o_sumet(),hweight); 
       FractionMET_DET_7bin->Fill(Event->pfmet().met_o_sumet(),hweight); 
	 }
       }
       
       //if(DetJets==2) {cout<<DetJets<<" Jets: "<<"Lumi "<<Event->evtHdr().lumi()<<" Event "<<Event->evtHdr().event()<<" Run Number "<<Event->evtHdr().runNo()<<endl;}
       //if(DetJets==3) {cout<<DetJets<<" Jets: "<<"Lumi "<<Event->evtHdr().lumi()<<" Event "<<Event->evtHdr().event()<<" Run Number "<<Event->evtHdr().runNo()<<endl;}
       //if(DetJets==4) {cout<<DetJets<<" Jets: "<<"Lumi "<<Event->evtHdr().lumi()<<" Event "<<Event->evtHdr().event()<<" Run Number "<<Event->evtHdr().runNo()<<endl;}
      
       //if(Event->evtHdr().runNo()==251252){
	 //if(Event->evtHdr().lumi()==158) { cout<<Event->evtHdr().event()<<endl;}
	 // if(DetJets==3 && Event->evtHdr().lumi()==158) { cout<<Event->evtHdr().event()<<" "<<DetJets<<" "<<Event->pfjetchs(0).ptCor()<<" "<<Event->pfjetchs(1).ptCor()<<" "<<Event->pfjetchs(2).ptCor()<<" "<<Event->pfjetchs(0).y()<<" "<<Event->pfjetchs(1).y()<<" "<<Event->pfjetchs(2).y()<<endl;}
	 //if(DetJets==2 && Event->evtHdr().lumi()==158) { cout<<Event->evtHdr().event()<<" "<<DetJets<<" "<<Event->pfjetchs(0).ptCor()<<" "<<Event->pfjetchs(1).ptCor()<<" "<<Event->pfjetchs(0).y()<<" "<<Event->pfjetchs(1).y()<<endl;}
	 //if(DetJets==4 && Event->evtHdr().lumi()==158) { cout<<Event->evtHdr().event()<<" "<<DetJets<<" "<<Event->pfjetchs(0).ptCor()<<" "<<Event->pfjetchs(1).ptCor()<<" "<<Event->pfjetchs(2).ptCor()<<" "<<Event->pfjetchs(0).y()<<" "<<Event->pfjetchs(1).y()<<" "<<Event->pfjetchs(2).y()<<endl;}
	 //}
     } // end of event loop
       
     mInf->Close();
   
   }
   
   float EntriesNumAcceptance[100][8];
   float EntriesDenAcceptance[100][8];
   float EntriesNumPurity[100][8];
   float EntriesDenPurity[100][8];
   float EntriesNumBackground[100][8];
   float EntriesDenBackground[100][8];
   float EntriesNumStability[100][8];
   float EntriesDenStability[100][8];

   float Entries[8];
   float Errors[8];

   for(int i=0;i<100;i++){
     for(int j=0;j<8;j++){
       EntriesNumAcceptance[i][j]=0;
       EntriesDenAcceptance[i][j]=0;
       EntriesNumStability[i][j]=0;
       EntriesDenStability[i][j]=0;
       EntriesNumBackground[i][j]=0;
       EntriesDenBackground[i][j]=0;
       EntriesNumPurity[i][j]=0;
       EntriesDenPurity[i][j]=0;
       Entries[j]=0;
       Errors[j]=0;
     }
   }

   //Acceptance

   for(int i=1;i<TwoD_MatchedInclusiveJets->GetXaxis()->GetNbins()+1;i++){
     EntriesDenAcceptance[i][0]=Gen_MatchedInclusiveJets->GetBinContent(i);
     EntriesDenAcceptance[i][1]=pt_GENInclJet_1bin->GetBinContent(i);
     EntriesDenAcceptance[i][2]=pt_GENInclJet_2bin->GetBinContent(i);
     EntriesDenAcceptance[i][3]=pt_GENInclJet_3bin->GetBinContent(i);
     EntriesDenAcceptance[i][4]=pt_GENInclJet_4bin->GetBinContent(i);
     EntriesDenAcceptance[i][5]=pt_GENInclJet_5bin->GetBinContent(i);
     EntriesDenAcceptance[i][6]=pt_GENInclJet_6bin->GetBinContent(i);
     EntriesDenAcceptance[i][7]=pt_GENInclJet_7bin->GetBinContent(i);
     
     for(int j=1;j<TwoD_MatchedInclusiveJets->GetYaxis()->GetNbins()+1;j++){
       EntriesNumAcceptance[i][0]=EntriesNumAcceptance[i][0]+TwoD_MatchedInclusiveJets->GetBinContent(j,i);
       EntriesNumAcceptance[i][1]=EntriesNumAcceptance[i][1]+TwoD_MatchedInclusiveJets_1bin->GetBinContent(j,i);
       EntriesNumAcceptance[i][2]=EntriesNumAcceptance[i][2]+TwoD_MatchedInclusiveJets_2bin->GetBinContent(j,i);
       EntriesNumAcceptance[i][3]=EntriesNumAcceptance[i][3]+TwoD_MatchedInclusiveJets_3bin->GetBinContent(j,i);
       EntriesNumAcceptance[i][4]=EntriesNumAcceptance[i][4]+TwoD_MatchedInclusiveJets_4bin->GetBinContent(j,i);
       EntriesNumAcceptance[i][5]=EntriesNumAcceptance[i][5]+TwoD_MatchedInclusiveJets_5bin->GetBinContent(j,i);
       EntriesNumAcceptance[i][6]=EntriesNumAcceptance[i][6]+TwoD_MatchedInclusiveJets_6bin->GetBinContent(j,i);
       EntriesNumAcceptance[i][7]=EntriesNumAcceptance[i][7]+TwoD_MatchedInclusiveJets_7bin->GetBinContent(j,i);
     }
     
     for(int j=0;j<8;j++){
       if(EntriesDenAcceptance[i][j]!=0){
	 Entries[j]=EntriesNumAcceptance[i][j]/EntriesDenAcceptance[i][j];
       }
     }

     AcceptancePtJets->SetBinContent(i,Entries[0]);
     AcceptancePtJets->SetBinError(i,Errors[0]);
     AcceptancePtJets_1bin->SetBinContent(i,Entries[1]);
     AcceptancePtJets_2bin->SetBinContent(i,Entries[2]);
     AcceptancePtJets_3bin->SetBinContent(i,Entries[3]);
     AcceptancePtJets_4bin->SetBinContent(i,Entries[4]);
     AcceptancePtJets_5bin->SetBinContent(i,Entries[5]);
     AcceptancePtJets_6bin->SetBinContent(i,Entries[6]);
     AcceptancePtJets_7bin->SetBinContent(i,Entries[7]);
     AcceptancePtJets_1bin->SetBinError(i,Errors[0]);
     AcceptancePtJets_2bin->SetBinError(i,Errors[0]);
     AcceptancePtJets_3bin->SetBinError(i,Errors[0]);
     AcceptancePtJets_4bin->SetBinError(i,Errors[0]);
     AcceptancePtJets_5bin->SetBinError(i,Errors[0]);
     AcceptancePtJets_6bin->SetBinError(i,Errors[0]);
     AcceptancePtJets_7bin->SetBinError(i,Errors[0]);
   }
   
   //Purity: matched detector and hadron level/ matched detector level
      
   for(int i=1;i<TwoD_MatchedInclusiveJets->GetXaxis()->GetNbins()+1;i++){
     EntriesNumPurity[i][0]=TwoD_MatchedInclusiveJets->GetBinContent(i,i);

     EntriesNumPurity[i][1]=TwoD_MatchedInclusiveJets_1bin->GetBinContent(i,i);
     EntriesNumPurity[i][2]=TwoD_MatchedInclusiveJets_2bin->GetBinContent(i,i);
     EntriesNumPurity[i][3]=TwoD_MatchedInclusiveJets_3bin->GetBinContent(i,i);
     EntriesNumPurity[i][4]=TwoD_MatchedInclusiveJets_4bin->GetBinContent(i,i);
     EntriesNumPurity[i][5]=TwoD_MatchedInclusiveJets_5bin->GetBinContent(i,i);
     EntriesNumPurity[i][6]=TwoD_MatchedInclusiveJets_6bin->GetBinContent(i,i);
     EntriesNumPurity[i][7]=TwoD_MatchedInclusiveJets_7bin->GetBinContent(i,i);

     for(int j=1;j<TwoD_MatchedInclusiveJets->GetYaxis()->GetNbins()+1;j++){
       EntriesDenPurity[i][0]=EntriesDenPurity[i][0]+TwoD_MatchedInclusiveJets->GetBinContent(j,i);
       EntriesDenPurity[i][1]=EntriesDenPurity[i][1]+TwoD_MatchedInclusiveJets_1bin->GetBinContent(j,i);
       EntriesDenPurity[i][2]=EntriesDenPurity[i][2]+TwoD_MatchedInclusiveJets_2bin->GetBinContent(j,i);
       EntriesDenPurity[i][3]=EntriesDenPurity[i][3]+TwoD_MatchedInclusiveJets_3bin->GetBinContent(j,i);
       EntriesDenPurity[i][4]=EntriesDenPurity[i][4]+TwoD_MatchedInclusiveJets_4bin->GetBinContent(j,i);
       EntriesDenPurity[i][5]=EntriesDenPurity[i][5]+TwoD_MatchedInclusiveJets_5bin->GetBinContent(j,i);
       EntriesDenPurity[i][6]=EntriesDenPurity[i][6]+TwoD_MatchedInclusiveJets_6bin->GetBinContent(j,i);
       EntriesDenPurity[i][7]=EntriesDenPurity[i][7]+TwoD_MatchedInclusiveJets_7bin->GetBinContent(j,i);
     }
     
     for(int j=0;j<8;j++){
       if(EntriesDenPurity[i][j]!=0){
	 Entries[j]=EntriesNumPurity[i][j]/EntriesDenPurity[i][j];
       }
     }

     PurityPtJets->SetBinContent(i,Entries[0]);
     PurityPtJets_1bin->SetBinContent(i,Entries[1]);
     PurityPtJets_2bin->SetBinContent(i,Entries[2]);
     PurityPtJets_3bin->SetBinContent(i,Entries[3]);
     PurityPtJets_4bin->SetBinContent(i,Entries[4]);
     PurityPtJets_5bin->SetBinContent(i,Entries[5]);
     PurityPtJets_6bin->SetBinContent(i,Entries[6]);
     PurityPtJets_7bin->SetBinContent(i,Entries[7]);
     
   }
   
   //Stability:matched detector and hadron level / matched gen level

   for(int i=1;i<TwoD_MatchedInclusiveJets->GetXaxis()->GetNbins()+1;i++){
     EntriesNumStability[i][0]=TwoD_MatchedInclusiveJets->GetBinContent(i,i);

     EntriesNumStability[i][1]=TwoD_MatchedInclusiveJets_1bin->GetBinContent(i,i);
     EntriesNumStability[i][2]=TwoD_MatchedInclusiveJets_2bin->GetBinContent(i,i);
     EntriesNumStability[i][3]=TwoD_MatchedInclusiveJets_3bin->GetBinContent(i,i);
     EntriesNumStability[i][4]=TwoD_MatchedInclusiveJets_4bin->GetBinContent(i,i);
     EntriesNumStability[i][5]=TwoD_MatchedInclusiveJets_5bin->GetBinContent(i,i);
     EntriesNumStability[i][6]=TwoD_MatchedInclusiveJets_6bin->GetBinContent(i,i);
     EntriesNumStability[i][7]=TwoD_MatchedInclusiveJets_7bin->GetBinContent(i,i);

     for(int j=1;j<TwoD_MatchedInclusiveJets->GetYaxis()->GetNbins()+1;j++){
       EntriesDenStability[i][0]=EntriesDenStability[i][0]+TwoD_MatchedInclusiveJets->GetBinContent(i,j);
       EntriesDenStability[i][1]=EntriesDenStability[i][1]+TwoD_MatchedInclusiveJets_1bin->GetBinContent(i,j);
       EntriesDenStability[i][2]=EntriesDenStability[i][2]+TwoD_MatchedInclusiveJets_2bin->GetBinContent(i,j);
       EntriesDenStability[i][3]=EntriesDenStability[i][3]+TwoD_MatchedInclusiveJets_3bin->GetBinContent(i,j);
       EntriesDenStability[i][4]=EntriesDenStability[i][4]+TwoD_MatchedInclusiveJets_4bin->GetBinContent(i,j);
       EntriesDenStability[i][5]=EntriesDenStability[i][5]+TwoD_MatchedInclusiveJets_5bin->GetBinContent(i,j);
       EntriesDenStability[i][6]=EntriesDenStability[i][6]+TwoD_MatchedInclusiveJets_6bin->GetBinContent(i,j);
       EntriesDenStability[i][7]=EntriesDenStability[i][7]+TwoD_MatchedInclusiveJets_7bin->GetBinContent(i,j);
     }

     for(int j=0;j<8;j++){
       if(EntriesDenStability[i][j]!=0){
	 Entries[j]=EntriesNumStability[i][j]/EntriesDenStability[i][j];
       }
     }

     StabilityPtJets->SetBinContent(i,Entries[0]);
     StabilityPtJets_1bin->SetBinContent(i,Entries[1]);
     StabilityPtJets_2bin->SetBinContent(i,Entries[2]);
     StabilityPtJets_3bin->SetBinContent(i,Entries[3]);
     StabilityPtJets_4bin->SetBinContent(i,Entries[4]);
     StabilityPtJets_5bin->SetBinContent(i,Entries[5]);
     StabilityPtJets_6bin->SetBinContent(i,Entries[6]);
     StabilityPtJets_7bin->SetBinContent(i,Entries[7]);
     
   }
  
   //Background

   for(int i=1;i<TwoD_MatchedInclusiveJets->GetXaxis()->GetNbins()+1;i++){
     EntriesDenBackground[i][0]=PF_MatchedInclusiveJets->GetBinContent(i);

     EntriesDenBackground[i][1]=pt_DETInclJet_1bin->GetBinContent(i);
     EntriesDenBackground[i][2]=pt_DETInclJet_2bin->GetBinContent(i);
     EntriesDenBackground[i][3]=pt_DETInclJet_3bin->GetBinContent(i);
     EntriesDenBackground[i][4]=pt_DETInclJet_4bin->GetBinContent(i);
     EntriesDenBackground[i][5]=pt_DETInclJet_5bin->GetBinContent(i);
     EntriesDenBackground[i][6]=pt_DETInclJet_6bin->GetBinContent(i);
     EntriesDenBackground[i][7]=pt_DETInclJet_7bin->GetBinContent(i);

     for(int j=1;j<TwoD_MatchedInclusiveJets->GetYaxis()->GetNbins()+1;j++){
       EntriesNumBackground[i][0]=EntriesNumBackground[i][0]+TwoD_MatchedInclusiveJets->GetBinContent(i,j);

       EntriesNumBackground[i][1]=EntriesNumBackground[i][1]+TwoD_MatchedInclusiveJets_1bin->GetBinContent(i,j);
       EntriesNumBackground[i][2]=EntriesNumBackground[i][2]+TwoD_MatchedInclusiveJets_2bin->GetBinContent(i,j);
       EntriesNumBackground[i][3]=EntriesNumBackground[i][3]+TwoD_MatchedInclusiveJets_3bin->GetBinContent(i,j);
       EntriesNumBackground[i][4]=EntriesNumBackground[i][4]+TwoD_MatchedInclusiveJets_4bin->GetBinContent(i,j);
       EntriesNumBackground[i][5]=EntriesNumBackground[i][5]+TwoD_MatchedInclusiveJets_5bin->GetBinContent(i,j);
       EntriesNumBackground[i][6]=EntriesNumBackground[i][6]+TwoD_MatchedInclusiveJets_6bin->GetBinContent(i,j);
       EntriesNumBackground[i][7]=EntriesNumBackground[i][7]+TwoD_MatchedInclusiveJets_7bin->GetBinContent(i,j);
     }

     for(int j=0;j<8;j++){
       if(EntriesDenBackground[i][j]!=0){
	 Entries[j]=1-EntriesNumBackground[i][j]/EntriesDenBackground[i][j];
       }
     }

     BackgroundPtJets->SetBinContent(i,Entries[0]);
     BackgroundPtJets_1bin->SetBinContent(i,Entries[1]);
     BackgroundPtJets_2bin->SetBinContent(i,Entries[2]);
     BackgroundPtJets_3bin->SetBinContent(i,Entries[3]);
     BackgroundPtJets_4bin->SetBinContent(i,Entries[4]);
     BackgroundPtJets_5bin->SetBinContent(i,Entries[5]);
     BackgroundPtJets_6bin->SetBinContent(i,Entries[6]);
     BackgroundPtJets_7bin->SetBinContent(i,Entries[7]);
     
   }

   int Ptbinwidth[73] = {1,3,3,4,4,5,6,6,7,8,10,10,13,17,19,20,21,22,24,25,27,28,30,32,33,35,38,39,41,44,46,48,51,53,56,59,62,65,69,71,76,79,83,87,91,96,100,106,110,116,122,128,134,140,147,154,162,170,177,187,195,205,215,225,236,248,259,272,285,299,313,328,283};   // first element is by default 1.0

   if(!mIsMCarlo){

     double Entries=0;
     for(int j=1;j<pt_DETInclJet_1bin->GetXaxis()->GetNbins()+1;j++){
       Entries=pt_DETInclJet_1bin->GetBinContent(j)/(2*1119.48*Ptbinwidth[j]);
       pt_DETInclJetCrossSectNorm_1bin->SetBinContent(j,Entries);
       Entries=pt_DETInclJet_2bin->GetBinContent(j)/(2*1119.48*Ptbinwidth[j]);
       pt_DETInclJetCrossSectNorm_2bin->SetBinContent(j,Entries);
       Entries=pt_DETInclJet_3bin->GetBinContent(j)/(2*1119.48*Ptbinwidth[j]);
       pt_DETInclJetCrossSectNorm_3bin->SetBinContent(j,Entries);
       Entries=pt_DETInclJet_4bin->GetBinContent(j)/(2*1119.48*Ptbinwidth[j]);
       pt_DETInclJetCrossSectNorm_4bin->SetBinContent(j,Entries);
       Entries=pt_DETInclJet_5bin->GetBinContent(j)/(2*1119.48*Ptbinwidth[j]);
       pt_DETInclJetCrossSectNorm_5bin->SetBinContent(j,Entries);
       Entries=pt_DETInclJet_6bin->GetBinContent(j)/(2*1119.48*Ptbinwidth[j]);
       pt_DETInclJetCrossSectNorm_6bin->SetBinContent(j,Entries);
       Entries=pt_DETInclJet_7bin->GetBinContent(j)/(2*3358.46*Ptbinwidth[j]);
       pt_DETInclJetCrossSectNorm_7bin->SetBinContent(j,Entries); 
     }
   }

   if(mIsMCarlo){

     double Entries=0;
     for(int j=1;j<pt_DETInclJet_1bin->GetXaxis()->GetNbins()+1;j++){
       Entries=pt_DETInclJet_1bin->GetBinContent(j)/(Ptbinwidth[j]);
       pt_DETInclJetCrossSectNorm_1bin->SetBinContent(j,Entries);
       Entries=pt_DETInclJet_2bin->GetBinContent(j)/(Ptbinwidth[j]);
       pt_DETInclJetCrossSectNorm_2bin->SetBinContent(j,Entries);
       Entries=pt_DETInclJet_3bin->GetBinContent(j)/(Ptbinwidth[j]);
       pt_DETInclJetCrossSectNorm_3bin->SetBinContent(j,Entries);
       Entries=pt_DETInclJet_4bin->GetBinContent(j)/(Ptbinwidth[j]);
       pt_DETInclJetCrossSectNorm_4bin->SetBinContent(j,Entries);
       Entries=pt_DETInclJet_5bin->GetBinContent(j)/(Ptbinwidth[j]);
       pt_DETInclJetCrossSectNorm_5bin->SetBinContent(j,Entries);
       Entries=pt_DETInclJet_6bin->GetBinContent(j)/(Ptbinwidth[j]);
       pt_DETInclJetCrossSectNorm_6bin->SetBinContent(j,Entries);
       Entries=pt_DETInclJet_7bin->GetBinContent(j)/(Ptbinwidth[j]);
       pt_DETInclJetCrossSectNorm_7bin->SetBinContent(j,Entries);

       Entries=pt_GENInclJet_1bin->GetBinContent(j)/(Ptbinwidth[j]);
       pt_GENInclJetCrossSectNorm_1bin->SetBinContent(j,Entries);
       Entries=pt_GENInclJet_2bin->GetBinContent(j)/(Ptbinwidth[j]);
       pt_GENInclJetCrossSectNorm_2bin->SetBinContent(j,Entries);
       Entries=pt_GENInclJet_3bin->GetBinContent(j)/(Ptbinwidth[j]);
       pt_GENInclJetCrossSectNorm_3bin->SetBinContent(j,Entries);
       Entries=pt_GENInclJet_4bin->GetBinContent(j)/(Ptbinwidth[j]);
       pt_GENInclJetCrossSectNorm_4bin->SetBinContent(j,Entries);
       Entries=pt_GENInclJet_5bin->GetBinContent(j)/(Ptbinwidth[j]);
       pt_GENInclJetCrossSectNorm_5bin->SetBinContent(j,Entries);
       Entries=pt_GENInclJet_6bin->GetBinContent(j)/(Ptbinwidth[j]);
       pt_GENInclJetCrossSectNorm_6bin->SetBinContent(j,Entries);
       Entries=pt_GENInclJet_7bin->GetBinContent(j)/(Ptbinwidth[j]);
       pt_GENInclJetCrossSectNorm_7bin->SetBinContent(j,Entries);
     }
   }

   TagAndProbeEff->Divide(TagAndProbeNum,TagAndProbeDen,1.,1.,"B");

   hist_leading_pt_HLT_Jet60U_eff->Divide(hist_leading_pt_emulated_Jet60,hist_leading_pt_all_Jet60,1.,1.,"B");
   hist_leading_eta_HLT_Jet60U_eff->Divide(hist_leading_eta_emulated_Jet60,hist_leading_eta_all_Jet60,1.,1.,"B");
   
   hist_leading_pt_HLT_Jet80U_eff->Divide(hist_leading_pt_emulated_Jet80,hist_leading_pt_all_Jet80,1.,1.,"B");
   hist_leading_eta_HLT_Jet80U_eff->Divide(hist_leading_eta_emulated_Jet80,hist_leading_eta_all_Jet80,1.,1.,"B");

   hist_leading_pt_HLT_Jet140U_eff->Divide(hist_leading_pt_emulated_Jet140,hist_leading_pt_all_Jet140,1.,1.,"B");
   hist_leading_eta_HLT_Jet140U_eff->Divide(hist_leading_eta_emulated_Jet140,hist_leading_eta_all_Jet140,1.,1.,"B");

   hist_leading_pt_HLT_Jet200U_eff->Divide(hist_leading_pt_emulated_Jet200,hist_leading_pt_all_Jet200,1.,1.,"B");
   hist_leading_eta_HLT_Jet200U_eff->Divide(hist_leading_eta_emulated_Jet200,hist_leading_eta_all_Jet200,1.,1.,"B");

   hist_leading_pt_HLT_Jet260U_eff->Divide(hist_leading_pt_emulated_Jet260,hist_leading_pt_all_Jet260,1.,1.,"B");
   hist_leading_eta_HLT_Jet260U_eff->Divide(hist_leading_eta_emulated_Jet260,hist_leading_eta_all_Jet260,1.,1.,"B");

   hist_leading_pt_HLT_Jet320U_eff->Divide(hist_leading_pt_emulated_Jet320,hist_leading_pt_all_Jet320,1.,1.,"B");
   hist_leading_eta_HLT_Jet320U_eff->Divide(hist_leading_eta_emulated_Jet320,hist_leading_eta_all_Jet320,1.,1.,"B");

   hist_leading_pt_HLT_Jet400U_eff->Divide(hist_leading_pt_emulated_Jet400,hist_leading_pt_all_Jet400,1.,1.,"B");
   hist_leading_eta_HLT_Jet400U_eff->Divide(hist_leading_eta_emulated_Jet400,hist_leading_eta_all_Jet400,1.,1.,"B");

   hist_leading_pt_HLT_Jet450U_eff->Divide(hist_leading_pt_emulated_Jet450,hist_leading_pt_all_Jet450,1.,1.,"B");
   hist_leading_eta_HLT_Jet450U_eff->Divide(hist_leading_eta_emulated_Jet450,hist_leading_eta_all_Jet450,1.,1.,"B");

   hist_leading_pt_HLT_Jet500U_eff->Divide(hist_leading_pt_emulated_Jet500,hist_leading_pt_all_Jet500,1.,1.,"B");
   hist_leading_eta_HLT_Jet500U_eff->Divide(hist_leading_eta_emulated_Jet500,hist_leading_eta_all_Jet500,1.,1.,"B");
   
 } // closing analyze() function



Analysis_Template_MC::~Analysis_Template_MC()
{
}


DEFINE_FWK_MODULE(Analysis_Template_MC);

