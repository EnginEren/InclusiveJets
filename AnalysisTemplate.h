#ifndef My_azimuthal_MC_h
#define My_azimuthal_MC_h

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "SMPJ/AnalysisFW/interface/QCDJet.h"
#include "SMPJ/AnalysisFW/interface/QCDEvent.h"
#include "SMPJ/AnalysisFW/interface/QCDEventHdr.h"
#include "SMPJ/AnalysisFW/interface/QCDCaloJet.h"
#include "SMPJ/AnalysisFW/interface/QCDPFJet.h"
#include "SMPJ/AnalysisFW/interface/QCDMET.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "SMPJ/AnalysisFW/plugins/JECs.h"


#include "TTree.h"
#include "TH1F.h"
#include "TFile.h"
#include "TProfile.h"
#include <TMath.h>
#include <boost/shared_ptr.hpp>
using namespace edm;
using namespace std;



class Analysis_Template_MC : public edm::EDAnalyzer
 {

  //typedef reco::Particle::LorentzVector LorentzVector;

  public:
    explicit Analysis_Template_MC(edm::ParameterSet const& cfg);
    virtual void beginJob();
    virtual void analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup);
    virtual void endJob();
    virtual ~Analysis_Template_MC();


  private:
    //---- configurable parameters --------
    std::string mTreeName,mDirName, mGlobalTag, mjettype;
    std::vector<std::string> mFileName;
    double mMinPt, mYMax;
    int    mJetID;        // looseID==1 tightID==2
    int    mprintOk;       // noPrint=0  Print=1
    bool mIsMCarlo;
    std::vector<std::string> mJECUncSrcNames;

    edm::Service<TFileService> fs;
    TTree *mTree;
    TFile *mInf, *mPuf;
    TDirectoryFile *mDir;

    //---- TREE variable --------
    QCDEvent *Event;

    //---- JECs variable --------
    JECs *jecs;

    //--------- Histogram Declaration --------------------//
    // Vertices
    TH1F *num_of_Vtx;
    TH1F *num_of_VtxGood;

    TH1F *mc_pthat;
    TH1F *mc_pthat_weighted;

    // PF Jets
    TH1F *ak7chsPt;
    TH1F *pt0_DETJetHLTj60;
    TH1F *pt1_DETJetHLTj60;
    TH1F *pt0_DETJetHLTHt800;
    TH1F *pt1_DETJetHLTHt800;
    TH1F *pt0_DETJet;
    TH1F *pt1_DETJet;
   
    
    TH1F *Xsect;

    ///Measurement Gen jets
    TH1F *pt0_GENJet;
    TH1F *pt1_GENJet;
    TH1F *y0_GENJet;
    TH1F *y1_GENJet;
    TH1F *phi0_GENJet;
    TH1F *phi1_GENJet;
    TH1F *ak7GenPt;

 };

#endif

