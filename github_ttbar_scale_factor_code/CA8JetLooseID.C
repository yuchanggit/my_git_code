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

void CA8JetLooseID(TreeReader &data,vector<Int_t> *Loose_CA8Jet_index 
		      ){

     Int_t    CA8nJet      = data.GetInt("CA8nJet");
     Float_t* CA8jetPt     = data.GetPtrFloat("CA8jetPt");
     Float_t* CA8jetEta    = data.GetPtrFloat("CA8jetEta");
     Float_t* CA8jetCEmEF  = data.GetPtrFloat("CA8jetCEmEF");
     Float_t* CA8jetCHadEF = data.GetPtrFloat("CA8jetCHadEF");
     Float_t* CA8jetNEmEF  = data.GetPtrFloat("CA8jetNEmEF");
     Float_t* CA8jetNHadEF = data.GetPtrFloat("CA8jetNHadEF");
     Float_t* CA8jetCMulti = data.GetPtrFloat("CA8jetCMulti");
     Int_t*   CA8jetPassID = data.GetPtrInt("CA8jetPassID");

  for(Int_t i = 0; i < CA8nJet; i++){

//      if ( CA8jetPassID[i]!=1  ){continue;}
     if ( CA8jetNHadEF[i]>=0.99  ){continue;}
     if ( CA8jetNEmEF[i]>=0.99 ){continue;}
     if ( fabs( CA8jetEta[i] )<2.4 ){
        if ( CA8jetCHadEF[i]<=0  ){continue;}
        if ( CA8jetCMulti[i]<=0  ){continue;}
        if ( CA8jetCEmEF[i]>=0.990  ){continue;}
     }

       Loose_CA8Jet_index->push_back(i);
    
  }    

  //-----------------------------------------------------------------------------------//



}
