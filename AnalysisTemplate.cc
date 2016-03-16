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
#include "TFile.h"

#include "SMPJ/AnalysisFW/plugins/AnalysisTemplate.h"

#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"

using namespace std;

//---------------------------- Constructor Of The Class TriggerTurnOn -------------------------- //
Analysis_Template_MC::Analysis_Template_MC(edm::ParameterSet const& cfg)
{
     mFileName       = cfg.getParameter<std::vector<std::string>>               ("filename");
     mTreeName       = cfg.getParameter<std::string>               ("treename");
     mDirName        = cfg.getParameter<std::string>               ("dirname");

     mGlobalTag      = cfg.getParameter<std::string>               ("pseudoglobaltag");
     mjettype        = cfg.getParameter<std::string>               ("jettype");

     mMinPt     = cfg.getParameter<double> ("minPt");
     mYMax      = cfg.getParameter<double> ("ymax");
     mJetID     = cfg.getParameter<int>  ("JetID");

     mprintOk   = cfg.getParameter<int>  ("printOk");

     mIsMCarlo       = cfg.getUntrackedParameter<bool>             ("isMCarlo");
     //mJECUncSrcNames = cfg.getParameter<std::vector<std::string> > ("jecUncSrcNames");
}

//------------------------------ Declaration Of The Function beginjob() ------------------------//
void Analysis_Template_MC::beginJob()
 {

     
          
     
     //------------------ Init jec constructor --------------------------- //
     //jecs = new JECs(mIsMCarlo, mGlobalTag, mjettype);

     //------------------ Histogram Booking --------------------------- //

      double Ptbinning[81] = {0, 1, 5, 6, 8, 10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84,97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468,507, 548, 592, 638, 686, 737, 790, 846, 905, 967,1032, 1101, 1172, 1248, 1327, 1410, 1497, 1588, 1684, 1784, 1890, 2000,2116, 2238, 2366, 2500, 2640, 2787, 2941, 3103, 3273, 3450, 3637, 3832,4037, 4252, 4477, 4713, 4961, 5220, 5492, 5777, 6076, 6389, 6717, 7000};
     num_of_Vtx     = fs->make<TH1F>("num_of_Vtx","num_of_Vtx",100,0.,100.);
     num_of_VtxGood = fs->make<TH1F>("num_of_VtxGood","num_of_VtxGood",100,0.,100.);
     ak7chsPt = fs->make<TH1F>("ak7chsPt","ak7chsPt",80,Ptbinning);
     ak7chsPt->Sumw2();
     Xsect = fs->make<TH1F>("Xsect","Cross-Section",80,Ptbinning);
     pt0_DETJetHLTj60  = fs->make<TH1F>("pt0_DETJetHLTj60","pt0_DETJetHLTj60",80,Ptbinning); pt0_DETJetHLTj60->Sumw2();
     pt1_DETJetHLTj60  = fs->make<TH1F>("pt1_DETJetHLTj60","pt1_DETJetHLTj60",80,Ptbinning); pt1_DETJetHLTj60->Sumw2();

     pt0_DETJetHLTHt800  = fs->make<TH1F>("pt0_DETJetHLTHt800","pt0_DETJetHLTHt800",80,Ptbinning); pt1_DETJetHLTHt800->Sumw2();
     pt1_DETJetHLTHt800  = fs->make<TH1F>("pt1_DETJetHLTHt800","pt1_DETJetHLTHt800",80,Ptbinning); pt1_DETJetHLTHt800->Sumw2();

     pt0_DETJet  = fs->make<TH1F>("pt0_DETJet","pt0_DETJet",80,Ptbinning); pt0_DETJet->Sumw2();
     pt1_DETJet  = fs->make<TH1F>("pt1_DETJet","pt1_DETJet",80,Ptbinning); pt1_DETJet->Sumw2();


     mc_pthat = fs->make<TH1F>("mc_pthat","mc_pthat",200,0.,2000.);
     mc_pthat_weighted = fs->make<TH1F>("mc_pthat_weighted","mc_pthat_weighted",200,0.,2000.);
     mc_pthat_weighted->Sumw2();


     pt0_GENJet  = fs->make<TH1F>("pt0_GENJet","pt0_GENJet",80,Ptbinning); pt0_GENJet->Sumw2();
     pt1_GENJet  = fs->make<TH1F>("pt1_GENJet","pt1_GENJet",80,Ptbinning); pt1_GENJet->Sumw2();
     y0_GENJet = fs->make<TH1F>("y0_GENJet","y0_GENJet",60,-3.,3.); y0_GENJet->Sumw2();
     y1_GENJet = fs->make<TH1F>("y1_GENJet","y1_GENJet",60,-3.,3.); y1_GENJet->Sumw2();
     phi0_GENJet = fs->make<TH1F>("phi0_GENJet","phi0_GENJet",60, -TMath::Pi(),TMath::Pi()); phi0_GENJet->Sumw2();
     phi1_GENJet = fs->make<TH1F>("phi1_GENJet","phi1_GENJet",60, -TMath::Pi(),TMath::Pi()); phi1_GENJet->Sumw2();

     ak7GenPt = fs->make<TH1F>("ak7GenPt","ak7GenPt",80,Ptbinning);
     ak7GenPt->Sumw2();

 } // end of function beginJob()





 //------------------------ endjob() function declaration ---------------------- //
 void Analysis_Template_MC::endJob()
 {
   mInf->Close();




 } // closing endJob()





 //--------------------------- analyze() fuction declaration ------------------ //
void Analysis_Template_MC::analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup)
 {
    
     cout << "Size " << mFileName.size() << endl;
     for (unsigned ifile = 0; ifile < mFileName.size(); ifile++) {
            
            cout << "Opening file: " << mFileName[ifile] << endl;
            mInf = TFile::Open(mFileName[ifile].c_str());
            mDir = (TDirectoryFile*)mInf->Get(mDirName.c_str());
            mTree = (TTree*)mDir->Get(mTreeName.c_str());
            Event = new QCDEvent();
            TBranch *branch = mTree->GetBranch("events");
            branch->SetAddress(&Event);

            unsigned NEntries = mTree->GetEntries();
            cout<<"Reading TREE: "<<NEntries<<" events"<<endl;

            int decade = 0 ;

            float hweight=1.;  ///Initial value to one
            //----- TRIGGER Info--------

            int nJetTrig=10, nHTTrig=9;
            int ihltj[nJetTrig];
            int ihltht[nHTTrig];
            TString HLTJet[nJetTrig], HLTHT[nJetTrig];
            char trigtitle[200];
            int HLTJetPtN[10] = {40,60,80,140,200,260,320,400,450,500};
            int HLTHTN[9] = {200,250,300,350,400,475,600,650,800};

            for (int i=0; i<nJetTrig; i++){
                    sprintf(trigtitle, "HLT_PFJet%i_v4",HLTJetPtN[i]);
                    HLTJet[i] = trigtitle;
            }


            for (int i=0; i<nHTTrig; i++){
                    sprintf(trigtitle, "HLT_PFHT%i_v2",HLTHTN[i]);
                    HLTHT[i] = trigtitle;
            }

            for (int i=0; i < nHTTrig; i++)
                         ihltht[i] = -1;

            for (int i=0; i<nJetTrig; i++)
                         ihltj[i] = -1 ;

            TH1F *hTrigNames = (TH1F*)mInf->Get("ak7/TriggerNames");
            cout<<"Finding trigger mapping: "<<endl;
             

            //----------- loop over the X-axis labels -----------------
            for(int ibin=0;ibin<hTrigNames->GetNbinsX();ibin++) {
                    TString ss(hTrigNames->GetXaxis()->GetBinLabel(ibin+1));
                    for (int ii=0; ii<nJetTrig; ii++){
                        if (ss == HLTJet[ii]) {
                            ihltj[ii] = ibin;
                            continue;
                        }
                    }

                    for (int j=0; j < nHTTrig; j++){
                        if (ss == HLTHT[j]){
                            ihltht[j] = ibin;
                            continue;
                        }
                    
                    }
            }


         for (int ij=0; ij<nJetTrig; ij++){
                   if (ihltj[ij] == -1)
                        cout<<"The requested trigger ("<<HLTJet[ij]<<") is not found "<<endl;
                   else 
                        cout<<HLTJet[ij]<<" --> "<<ihltj[ij]<<endl;
         }

        for (int ij=0; ij<nHTTrig; ij++){
                   if (ihltht[ij] == -1)
                        cout<<"The requested trigger ("<<HLTHT[ij]<<") is not found "<<endl;
                   else 
                        cout<<HLTHT[ij]<<" --> "<<ihltht[ij]<<endl;
         }


        for(unsigned  l=0; l<NEntries; l++) {
        //for(unsigned  l=0; l<50; l++) {


        //----------- progress report -------------
        double progress = 10.0*l/(1.0*NEntries);
        int k = TMath::FloorNint(progress);
        if (k > decade)
          cout<<10*k<<" %"<<endl;
        decade = k;
        //----------- read the event --------------
        mTree->GetEntry(l);

        ///pthat & mc_weight
        float pthat = Event->evtHdr().pthat();
        float mc_weight = 1.0; 
        //if(mprintOk==1) printf("\npthat=%f  mc_weight=%e\n",pthat,mc_weight);
        
        
        if (mIsMCarlo){
            if (ifile == 0){mc_weight = 27540000 * Event->evtHdr().weight() / NEntries;}  // slice 100-200 
            if (ifile == 1){mc_weight = 1717000  * Event->evtHdr().weight() / NEntries;}  // slice 200-300
            if (ifile == 2){mc_weight = 351300 * Event->evtHdr().weight() / NEntries;}  // slice 300-500 
            if (ifile == 3){mc_weight = 31630 * Event->evtHdr().weight() / NEntries;}  // slice 500-700
            if (ifile == 4){mc_weight = 6802 * Event->evtHdr().weight() / NEntries;}  // slice 700-1000
            if (ifile == 5){mc_weight = 1206 * Event->evtHdr().weight() / NEntries;}  // slice 1000-1500
            if (ifile == 6){mc_weight = 120.4 * Event->evtHdr().weight() / NEntries;}  // slice 1500-2000
            if (ifile == 7){mc_weight = 25.25 * Event->evtHdr().weight() / NEntries;}  // slice 2000-Inf
        
        }
        
        
        hweight=mc_weight;
        mc_pthat->Fill(pthat);
        mc_pthat_weighted->Fill(pthat,hweight);
        /////////////////////////////////////////////////////////////////////////////////////////////
        ///Examine GenJets
        unsigned n_genJets = Event->nGenJets();


        ///Apply Jet cuts.  Very General to all existing GEN Jets
        int GENjet_ok[100]; for(int ii=0;ii<100;++ii){GENjet_ok[ii]=0;}

        for(unsigned j=0; j< n_genJets; ++j){
            if(Event->genjet(j).pt()<mMinPt) continue;
            //if(fabs(Event->genjet(j).Rapidity())>mYMax) continue;
            GENjet_ok[j]=1;
            ak7GenPt->Fill(Event->genjet(j).pt(),hweight);
        
        }



        /// Keep events where leading Jet[0] and Jet[1] survived cuts
        if((GENjet_ok[0]==1)&&(GENjet_ok[1]==1)) {

           ///////////////////////////////////// Measurement with Gen Jets ///////////////////////////////////////////////////

              //float ptmax_gen=Event->genjet(0).pt();

                  
              pt0_GENJet->Fill(Event->genjet(0).pt(),hweight);
              pt1_GENJet->Fill(Event->genjet(1).pt(),hweight);
              y0_GENJet->Fill(Event->genjet(0).Rapidity(),hweight);
              y1_GENJet->Fill(Event->genjet(1).Rapidity(),hweight);
              phi0_GENJet->Fill(Event->genjet(0).phi(),hweight);
              phi1_GENJet->Fill(Event->genjet(1).phi(),hweight);



        } //end of GEN Jets


        /////////////////////////////////////////////////////////////////////////////////////////////
        /// PFJets
        /////////////////////////////////////////////Vertex!!!!/////////////////////////////////////
        unsigned n_PFJets = Event->nPFJetsCHS();
        int DETjet_ok[100]; for(int ii=0;ii<100;++ii){DETjet_ok[ii]=0;}
        int nJetTrig = 10;
        float prescalej[nJetTrig];
        float prescaleHT[nHTTrig];
        TString HLTJet[nJetTrig];
        bool hltPassj[nJetTrig];
        bool hltPassHT[nHTTrig];
        //jecs->JEC_corrections(Event, n_PFJets, mIsMCarlo);

        //int DETjet_ok[100]; for(int ii=0;ii<100;++ii){DETjet_ok[ii]=0;}

        /// Vertex selection
        //if(mprintOk==1) cout<<"Vertex info: numVtx="<<Event->evtHdr().nVtx()<<"  numVtxGood="<<Event->evtHdr().nVtxGood()<<"  isPVgood()="<<Event->evtHdr().isPVgood()<<"   pfRho="<<Event->evtHdr().pfRho()<<endl;

        num_of_Vtx->Fill(Event->evtHdr().nVtx());
        /// Keep events with PVgood
        if (Event->evtHdr().isPVgood() != 1) continue;

        num_of_VtxGood->Fill(Event->evtHdr().nVtxGood());
        /// Dump all jets No cuts only isPVgood()
        if(mprintOk != 1) continue;
            //printf("Number of PFJets=%d\n",n_PFJets);
        
        
        
        if (!mIsMCarlo){
            
            for (int j=0; j<nJetTrig; j++){
                hltPassj[j] = false ;
                prescalej[j] = 1;
            }
            
            for (int j=0; j<nHTTrig; j++){
                hltPassHT[j] = false ;
                prescaleHT[j] = 1;
            }

            for(int j=0; j < nJetTrig; ++j){
                if (Event->fired(ihltj[j]) > 0) {
                   hltPassj[j] = true;
                   if(Event->preL1(ihltj[j])>=0){ prescalej[j] = Event->preL1(ihltj[j]) * Event->preHLT(ihltj[j]);}
                   else {prescalej[j] = Event->preHLT(ihltj[j]);}
                   //cout<<j<<" "<<HLTJet[j]<<" "<<prescalej[j]<<" PASS **********************************"<<Event->preL1(ihltj[j])<<" "<< Event->preHLT(ihltj[j])<<endl;
                    
            
                }
            
            }
       
            for(int j=0; j < nHTTrig; ++j){
                if (Event->fired(ihltht[j]) > 0) {
                   hltPassHT[j] = true;
                   if(Event->preL1(ihltht[j])>=0){ prescaleHT[j] = Event->preL1(ihltht[j]) * Event->preHLT(ihltht[j]);}
                   else {prescaleHT[j] = Event->preHLT(ihltht[j]);}
                   //cout<<j<<" "<<HLTJet[j]<<" "<<prescalej[j]<<" PASS **********************************"<<Event->preL1(ihltj[j])<<" "<< Event->preHLT(ihltj[j])<<endl;
                    
            
                }
            
            }

        }
       
    
        if ( hltPassHT[8] ){
            hweight = 1.0;  //for PF Jets
            
            if(prescaleHT[8] > 0 ) hweight = hweight * prescaleHT[8];
            if(prescaleHT[8] < 0 ) hweight = hweight * (-prescaleHT[8]);
            
            for(unsigned j=0; j<n_PFJets; ++j){
                if(Event->pfjetchs(j).ptCor() < mMinPt) continue;  //MinPt Cut
                //if(fabs(Event->pfjetchs(j).y()) < 0.5 ){    // for the first rapidity bin. 
                DETjet_ok[j] = 1;
                //ak7chsPt->Fill(Event->pfjetchs(j).ptCor(), hweight);
                //Xsect->Fill(Event->pfjetchs(j).ptCor());
                //}
       
            }
            
            
            if ( (DETjet_ok[0] == 1) && (DETjet_ok[1] == 1)){
              pt0_DETJetHLTHt800->Fill(Event->pfjetchs(0).ptCor(),hweight);
              pt1_DETJetHLTHt800->Fill(Event->pfjetchs(1).ptCor(),hweight);
 
            
            }
            
        
        }
            
        if (  hltPassj[1]  ){
            hweight = 1.0;  //for PF Jets
            if(prescalej[1] > 0 ) hweight = hweight * prescalej[1];
            if(prescalej[1] < 0 ) hweight = hweight * (-prescalej[1]);
            
            //if(prescaleHT[8] > 0 ) hweight = hweight * prescaleHT[8];
            //if(prescaleHT[8] < 0 ) hweight = hweight * (-prescaleHT[8]);
            
            for(unsigned j=0; j<n_PFJets; ++j){
                if(Event->pfjetchs(j).ptCor() < mMinPt) continue;  //MinPt Cut
                //if(fabs(Event->pfjetchs(j).y()) < 0.5 ){    // for the first rapidity bin. 
                DETjet_ok[j] = 1;
                ak7chsPt->Fill(Event->pfjetchs(j).ptCor(), hweight);
                Xsect->Fill(Event->pfjetchs(j).ptCor());
                //}
       
            }
            
            if ( (DETjet_ok[0] == 1) && (DETjet_ok[1] == 1)){
              pt0_DETJetHLTj60->Fill(Event->pfjetchs(0).ptCor(),hweight);
              pt1_DETJetHLTj60->Fill(Event->pfjetchs(1).ptCor(),hweight);
              if ( (Event->pfjetchs(0).ptCor() >= 114.00) && (Event->pfjetchs(0).ptCor() < 133.00) )
                    pt0_DETJet->Fill(Event->pfjetchs(0).ptCor());
                    //pt1_DETJet->Fill(Event->pfjetchs(1).ptCor());
                    
            }

        }      

        if (  hltPassj[2]  ){
            hweight = 1.0;  //for PF Jets
            if(prescalej[2] > 0 ) hweight = hweight * prescalej[2];
            if(prescalej[2] < 0 ) hweight = hweight * (-prescalej[2]);
            
            for(unsigned j=0; j<n_PFJets; ++j){
                if(Event->pfjetchs(j).ptCor() < mMinPt) continue;  //MinPt Cut
                //if(fabs(Event->pfjetchs(j).y()) < 0.5 ){    // for the first rapidity bin. 
                DETjet_ok[j] = 1;
                //}
       
            }
            
            if ( (DETjet_ok[0] == 1) && (DETjet_ok[1] == 1)){
              if ( (Event->pfjetchs(0).ptCor() >= 133.00) && (Event->pfjetchs(0).ptCor() < 220.00) )
                    pt0_DETJet->Fill(Event->pfjetchs(0).ptCor());
                    //pt1_DETJet->Fill(Event->pfjetchs(1).ptCor());
              // Other Stuff      
            }

        }      

         if (  hltPassj[3]  ){
            hweight = 1.0;  //for PF Jets
            if(prescalej[3] > 0 ) hweight = hweight * prescalej[3];
            if(prescalej[3] < 0 ) hweight = hweight * (-prescalej[3]);
            
            for(unsigned j=0; j<n_PFJets; ++j){
                if(Event->pfjetchs(j).ptCor() < mMinPt) continue;  //MinPt Cut
                //if(fabs(Event->pfjetchs(j).y()) < 0.5 ){    // for the first rapidity bin. 
                DETjet_ok[j] = 1;
                //}
       
            }
            
            if ( (DETjet_ok[0] == 1) && (DETjet_ok[1] == 1)){
              if ( (Event->pfjetchs(0).ptCor() >= 220.00) && (Event->pfjetchs(0).ptCor() < 300.00) )
                    pt0_DETJet->Fill(Event->pfjetchs(0).ptCor());
                    //pt1_DETJet->Fill(Event->pfjetchs(1).ptCor());
              // Other Stuff      
            }

        }         
          
        if (  hltPassj[4]  ){
            hweight = 1.0;  //for PF Jets
            if(prescalej[4] > 0 ) hweight = hweight * prescalej[4];
            if(prescalej[4] < 0 ) hweight = hweight * (-prescalej[4]);
            
            for(unsigned j=0; j<n_PFJets; ++j){
                if(Event->pfjetchs(j).ptCor() < mMinPt) continue;  //MinPt Cut
                //if(fabs(Event->pfjetchs(j).y()) < 0.5 ){    // for the first rapidity bin. 
                DETjet_ok[j] = 1;
                //}
       
            }
            
            if ( (DETjet_ok[0] == 1) && (DETjet_ok[1] == 1)){
              if ( (Event->pfjetchs(0).ptCor() >= 300.00) && (Event->pfjetchs(0).ptCor() < 430.00) )
                    pt0_DETJet->Fill(Event->pfjetchs(0).ptCor());
                    //pt1_DETJet->Fill(Event->pfjetchs(1).ptCor());
              // Other Stuff      
            }

        }         
      
        if (  hltPassj[5]  ){
            hweight = 1.0;  //for PF Jets
            if(prescalej[5] > 0 ) hweight = hweight * prescalej[5];
            if(prescalej[5] < 0 ) hweight = hweight * (-prescalej[5]);
            
            for(unsigned j=0; j<n_PFJets; ++j){
                if(Event->pfjetchs(j).ptCor() < mMinPt) continue;  //MinPt Cut
                //if(fabs(Event->pfjetchs(j).y()) < 0.5 ){    // for the first rapidity bin. 
                DETjet_ok[j] = 1;
                //}
       
            }
            
            if ( (DETjet_ok[0] == 1) && (DETjet_ok[1] == 1)){
              if ( (Event->pfjetchs(0).ptCor() >= 430.00) && (Event->pfjetchs(0).ptCor() < 507.00) )
                    pt0_DETJet->Fill(Event->pfjetchs(0).ptCor());
                    //pt1_DETJet->Fill(Event->pfjetchs(1).ptCor());
              // Other Stuff      
            }

        }
        
        if (  hltPassj[6]  ){
            hweight = 1.0;  //for PF Jets
            if(prescalej[6] > 0 ) hweight = hweight * prescalej[6];
            if(prescalej[6] < 0 ) hweight = hweight * (-prescalej[6]);
            
            for(unsigned j=0; j<n_PFJets; ++j){
                if(Event->pfjetchs(j).ptCor() < mMinPt) continue;  //MinPt Cut
                //if(fabs(Event->pfjetchs(j).y()) < 0.5 ){    // for the first rapidity bin. 
                DETjet_ok[j] = 1;
                //}
       
            }
            
            if ( (DETjet_ok[0] == 1) && (DETjet_ok[1] == 1)){
              if ( (Event->pfjetchs(0).ptCor() >= 570.00) && (Event->pfjetchs(0).ptCor() < 638.00) )
                    pt0_DETJet->Fill(Event->pfjetchs(0).ptCor());
                    //pt1_DETJet->Fill(Event->pfjetchs(1).ptCor());
              // Other Stuff      
            }

        }
        
        if (  hltPassj[7]  ){
            hweight = 1.0;  //for PF Jets
            if(prescalej[7] > 0 ) hweight = hweight * prescalej[7];
            if(prescalej[7] < 0 ) hweight = hweight * (-prescalej[7]);
            
            for(unsigned j=0; j<n_PFJets; ++j){
                if(Event->pfjetchs(j).ptCor() < mMinPt) continue;  //MinPt Cut
                //if(fabs(Event->pfjetchs(j).y()) < 0.5 ){    // for the first rapidity bin. 
                DETjet_ok[j] = 1;
                //}
       
            }
            
            if ( (DETjet_ok[0] == 1) && (DETjet_ok[1] == 1)){
              if ( (Event->pfjetchs(0).ptCor() >= 638.00) && (Event->pfjetchs(0).ptCor() < 737.00) )
                    pt0_DETJet->Fill(Event->pfjetchs(0).ptCor());
                    //pt1_DETJet->Fill(Event->pfjetchs(1).ptCor());
              // Other Stuff      
            }

        }
 
        if (  hltPassj[8]  ){
            hweight = 1.0;  //for PF Jets
            if(prescalej[8] > 0 ) hweight = hweight * prescalej[8];
            if(prescalej[8] < 0 ) hweight = hweight * (-prescalej[8]);
            
            for(unsigned j=0; j<n_PFJets; ++j){
                if(Event->pfjetchs(j).ptCor() < mMinPt) continue;  //MinPt Cut
                //if(fabs(Event->pfjetchs(j).y()) < 0.5 ){    // for the first rapidity bin. 
                DETjet_ok[j] = 1;
                //}
       
            }
            
            if ( (DETjet_ok[0] == 1) && (DETjet_ok[1] == 1)){
              if ( Event->pfjetchs(0).ptCor() > 737.00 )
                    pt0_DETJet->Fill(Event->pfjetchs(0).ptCor());
                    //pt1_DETJet->Fill(Event->pfjetchs(1).ptCor());
              // Other Stuff      
            }

        }



        } // end of event loop

    } // end of ifile loop

} // closing analyze() function



Analysis_Template_MC::~Analysis_Template_MC()
{
}


DEFINE_FWK_MODULE(Analysis_Template_MC);

