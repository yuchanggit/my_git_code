#include <map>
#include <vector>
#include <string>
#include <iostream>
#include <algorithm>
#include <TH1D.h>
#include <TH1F.h>
#include <TMath.h>
#include <TFile.h>
#include <TList.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TBranch.h>
#include <TRandom.h>
#include <TProfile.h>
#include <TLorentzVector.h>
#include <TSystemDirectory.h>
#include "untuplizer.h"

using namespace std;

int EventSelectionMuonChannel(TreeReader &data,
  vector<Int_t> *Tight_Muon_index, vector<Int_t> *Loose_Electron_index ,
  vector<Int_t> *Loose_CA8Jet_index, vector<Int_t> *Muon_pass_index,
  vector<Int_t> *leading_and_second_CA8Jet_index,
  vector<Int_t> *Btagged_CA8Jet_index ){

    Int_t    nMu          = data.GetInt("nMu");
    Float_t* muPt         = data.GetPtrFloat("muPt");
    Float_t* muEta        = data.GetPtrFloat("muEta");
    Float_t* muCorrPfIso  = data.GetPtrFloat("muCorrPfIso");

    Int_t    nEle         = data.GetInt("nEle");
    Float_t* elePt        = data.GetPtrFloat("elePt");
    Float_t* eleScEta     = data.GetPtrFloat("eleScEta");

    Int_t    CA8nJet      = data.GetInt("CA8nJet");
    Float_t* CA8jetPt     = data.GetPtrFloat("CA8jetPt");
    Float_t* CA8jetCSV    = data.GetPtrFloat("CA8jetCSV");

    Float_t  pfMetCorrPt  = data.GetFloat("pfMetCorrPt");


  //-----------------------------------------------------------------------------------//      in event selection the return -1 is fail, return 2 is pass 

  //-----------------------------------------------------------------------------------//
  // Muon pT>45, |Eta|<2.1, iso-Muon and the # of iso-Muon >= 1


    for (int i=0;i<Tight_Muon_index->size();i++){

        if ( muPt[ Tight_Muon_index->at(i) ] < 45  ){continue;}
        if ( muEta[ Tight_Muon_index->at(i) ] > 2.1  ){continue;}
        if ( ( muCorrPfIso[ Tight_Muon_index->at(i) ]/muPt[ Tight_Muon_index->at(i) ] ) > 0.05  ){continue;}
        Muon_pass_index->push_back( Tight_Muon_index->at(i) );
    }

    if ( Muon_pass_index->size() < 1 ){return -1;}

  //-----------------------------------------------------------------------------------//
  // reject addtional Muon with pT>10, |Eta|<2.1

    int local_counter1 =0;
 
    for (int i=0; i< Tight_Muon_index->size();i++){
        for (int j=0; j< Muon_pass_index->size();j++){
            if ( Tight_Muon_index->at(i) != Muon_pass_index->at(j) ){
                if( (muPt[ Tight_Muon_index->at(i) ] > 10) && (fabs( muEta[ Tight_Muon_index->at(i)]) < 2.1) )
                {local_counter1 = local_counter1 + 1;}                    
            }
        }    
    }
    
    if ( local_counter1 > 0 ){return -1;}

  //-----------------------------------------------------------------------------------//
  // reject addtional Electron with pT>12, |Eta|<2.1

    int local_counter2 =0;

    for (int i=0; i< Loose_Electron_index->size();i++){
                if( (elePt[ Loose_Electron_index->at(i) ] > 12) && (fabs( eleScEta[ Loose_Electron_index->at(i)]) < 2.1) )
                {local_counter2 = local_counter2 + 1;}
            }

    if ( local_counter2 > 0 ){return -1;}

  //-----------------------------------------------------------------------------------//
  // at least 2 CA8 jets

    if ( Loose_CA8Jet_index->size() < 2 ){return -1;}   

  //-----------------------------------------------------------------------------------//
  // leading jet pT>200, second jet pT>30

    float leading_CA8Jet_pT = -999 ;
    int leading_CA8Jet_index = -1 ;

    for (int i=0; i < Loose_CA8Jet_index->size();i++){
                if( CA8jetPt[ Loose_CA8Jet_index->at(i) ] > leading_CA8Jet_pT )
                        { leading_CA8Jet_pT = CA8jetPt[ Loose_CA8Jet_index->at(i) ];
                          leading_CA8Jet_index = Loose_CA8Jet_index->at(i);
                        }
    }

    float second_CA8Jet_pT = -999 ;
    int second_CA8Jet_index = -1 ;
    
    for (int i=0; i < Loose_CA8Jet_index->size();i++){
            if( Loose_CA8Jet_index->at(i) == leading_CA8Jet_index ){continue;}
                if( CA8jetPt[ Loose_CA8Jet_index->at(i) ] > second_CA8Jet_pT )
                        { second_CA8Jet_pT = CA8jetPt[ Loose_CA8Jet_index->at(i) ];
                          second_CA8Jet_index = Loose_CA8Jet_index->at(i);
                        }
    }
//cout<<"leading jet index:"<<leading_CA8Jet_index <<"second jet index:"<<second_CA8Jet_index<<endl;
//cout<<"hello leading jet"<<endl;

   if ( !( (leading_CA8Jet_pT>200) && (second_CA8Jet_pT>30) ) )
      { return -1;}

//cout<<"leading jet index:"<<leading_CA8Jet_index <<"second jet index:"<<second_CA8Jet_index<<endl;

   leading_and_second_CA8Jet_index->push_back(leading_CA8Jet_index);
   leading_and_second_CA8Jet_index->push_back(second_CA8Jet_index);

  //-----------------------------------------------------------------------------------//
  // at least 1 b-tagged jet

    int local_counter3 =0;

    for (int i=0; i < Loose_CA8Jet_index->size();i++){
            
                if( CA8jetCSV[ Loose_CA8Jet_index->at(i) ] > 0.244 )
                        { local_counter3 = local_counter3 + 1 ;
                          Btagged_CA8Jet_index->
                                       push_back( Loose_CA8Jet_index->at(i) );
			}
    }

    if ( local_counter3  < 1 ){return -1;}

  //-----------------------------------------------------------------------------------//
  // Missing Et>20

    if ( pfMetCorrPt  < 20 ){return -1;}

  //-----------------------------------------------------------------------------------//
  //

  // if an event pass all cuts and go here, then this event will send back 2 and be accepted 
  return 2;

}
