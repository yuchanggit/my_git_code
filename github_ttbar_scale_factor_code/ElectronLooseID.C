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

void ElectronLooseID(TreeReader &data, vector<Int_t> *Loose_Electron_index 
		      ){

  Int_t    nEle           = data.GetInt("nEle"); 
  Int_t*   elePassID      = data.GetPtrInt("elePassID");
  Float_t  eleRho         = data.GetFloat("eleRho");
  Float_t* eleEt          = data.GetPtrFloat("eleEt");
  Float_t* eleScEta       = data.GetPtrFloat("eleScEta");
  Float_t* eleUserTrkIso  = data.GetPtrFloat("eleUserTrkIso");
  Float_t* eleUserCalIso  = data.GetPtrFloat("eleUserCalIso");
  Float_t* elePt          = data.GetPtrFloat("elePt");
  Float_t* eleDelEtaIn    = data.GetPtrFloat("eleDelEtaIn");
  Float_t* eleDelPhiIn    = data.GetPtrFloat("eleDelPhiIn");
  Float_t* eleSigIhIh     = data.GetPtrFloat("eleSigIhIh");
  Float_t* eleHoE         = data.GetPtrFloat("eleHoE");
  Float_t* eleDxy         = data.GetPtrFloat("eleDxy");
  Float_t* eleDz          = data.GetPtrFloat("eleDz");
  Float_t* eleEoverP  	  = data.GetPtrFloat("eleEoverP");
  Float_t* eleCorrPfIso   = data.GetPtrFloat("eleCorrPfIso");
  Int_t*   elePassConv 	  = data.GetPtrInt("elePassConv");
  Float_t* eleMissingHits = data.GetPtrFloat("eleMissingHits");


  for(Int_t i = 0; i < nEle; i++){


//cout<<"nEle: "<<nEle<<endl;

    // barrel selection
    if( fabs(eleScEta[i]) > 0 && fabs(eleScEta[i]) < 1.4442 ){

       if ( fabs( eleDelEtaIn[i] ) > 0.007 ){continue;}
       if ( fabs( eleDelPhiIn[i] ) > 0.15 ){continue;}
       if (  eleSigIhIh[i]  > 0.01 ){continue;}
       if (  eleHoE[i]  > 0.12 ){continue;}
       if ( fabs( eleDxy[i] ) > 0.02 ){continue;}
       if ( fabs( eleDz[i] ) > 0.2 ){continue;}
       if ( fabs( eleEoverP[i] ) > 0.05 ){continue;}
       if ( ( eleCorrPfIso[i]/elePt[i] )  > 0.15 ){continue;}
       if (  elePassConv[i] !=  1 ){continue;}
       if (  eleMissingHits[i]   > 1 ){continue;}

       Loose_Electron_index->push_back(i);

//cout<<"barrel,i:"<< i<<" eta:"<<eleScEta[i]<<endl;

    }

    // endcap selection
    if( fabs(eleScEta[i]) > 1.566 && fabs(eleScEta[i]) < 2.5 ){

       if ( fabs( eleDelEtaIn[i] ) > 0.009 ){continue;}
       if ( fabs( eleDelPhiIn[i] ) > 0.1 ){continue;}
       if (  eleSigIhIh[i]  > 0.03 ){continue;}
       if (  eleHoE[i]  > 0.10 ){continue;}
       if ( fabs( eleDxy[i] ) > 0.02 ){continue;}
       if ( fabs( eleDz[i] ) > 0.2 ){continue;}
       if ( fabs( eleEoverP[i] ) > 0.05 ){continue;}

       if (elePt[i]<20){ 
	     if ( ( eleCorrPfIso[i]/elePt[i] )  > 0.1 ){continue;}}
       else if(elePt[i]>20){ 
	     if ( ( eleCorrPfIso[i]/elePt[i] )  > 0.15){continue;}} 

       if (  elePassConv[i] !=  1 ){continue;}
       if (  eleMissingHits[i]   > 1 ){continue;}

       Loose_Electron_index->push_back(i);

//cout<<"endcap,i:"<< i<<" eta:"<<eleScEta[i]<<endl;
    
    }
  }    

  //-----------------------------------------------------------------------------------//



}
