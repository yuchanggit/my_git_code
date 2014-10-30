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

void MuonTightID(TreeReader &data,vector<Int_t> *Tight_Muon_index 
		      ){

    Int_t    nMu          = data.GetInt("nMu");
    Float_t* muPt         = data.GetPtrFloat("muPt");
    Int_t*   isGlobalMuon = data.GetPtrInt("isGlobalMuon");
    Int_t*   isTrackerMuon= data.GetPtrInt("isTrackerMuon");
    Int_t*   muTrkLayers  = data.GetPtrInt("muTrkLayers");
    Int_t*   muPixelHits  = data.GetPtrInt("muPixelHits");
    Int_t*   muHits       = data.GetPtrInt("muHits");
    Int_t*   muMatches    = data.GetPtrInt("muMatches");
    Float_t* mudxy        = data.GetPtrFloat("mudxy");
    Float_t* mudz         = data.GetPtrFloat("mudz");


  for(Int_t i = 0; i < nMu; i++){

       if ( ( isGlobalMuon[i] != 1 ) && ( isTrackerMuon[i] != 1 ) ){continue;}
       if ( muTrkLayers[i] <= 5  ){continue;}
       if ( muPixelHits[i] <= 0  ){continue;}
       if ( muHits[i] <= 0  ){continue;}
       if ( muMatches[i] <= 1  ){continue;}
       if ( mudxy[i] >= 0.2  ){continue;}
       if ( mudz[i] >= 0.5  ){continue;}

       Tight_Muon_index->push_back(i);



    
  }    

  //-----------------------------------------------------------------------------------//



}
