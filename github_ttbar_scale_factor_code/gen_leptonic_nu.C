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

void gen_leptonic_nu(TreeReader &data,vector<Int_t> *gen_leptonic_nu_index 
		      ){

    Int_t    nGenPar     = data.GetInt("nGenPar");
    Int_t*   genParId    = data.GetPtrInt("genParId");
    Int_t*   genMomParId = data.GetPtrInt("genMomParId");
    Int_t*   genMo1      = data.GetPtrInt("genMo1");
    Int_t*   genParSt    = data.GetPtrInt("genParSt");    
    Float_t* genParE     = data.GetPtrFloat("genParE");
    Float_t* genParPt    = data.GetPtrFloat("genParPt");
    Float_t* genParEta   = data.GetPtrFloat("genParEta");
    Float_t* genParPhi   = data.GetPtrFloat("genParPhi");

  //-----------------------------------------------------------------------------------//

  // finding muon neutrino <- W boson <- top quark 
  // new version to find status=1 neutrino

  for(Int_t i = 0; i < nGenPar; i++){

     if (genParId[i]!=14 && genParId[i]!= -14){continue;}
     if (genParSt[i]!=1){continue;}

//     cout<<"gen paricle index i: "<< i<<endl;                             
//     cout<<"genParId[i]: "<<genParId[i] <<endl;
//     cout<<"genParSt[i]: "<<genParSt[i] <<endl;

     int mom_id = genMomParId[i]; 
     int mom_index = genMo1[i];
     int grandmom_index;
     int W_index=0;
     //set zero to avoid cannot find W and make error later to require W decay from top

     for(int j=0; j > -100; j++){// find W boson (mother of nu) index

//        cout<<"mon index i: "<< mom_index<<endl;
//        cout<<"mom Id: "<< mom_id<<endl;                    
//        cout<<"mom St: "<< genParSt[mom_index]<<endl;

        if(mom_id== 24 || mom_id== -24){

//           cout<<"find W boson when trace back the nu mother"<<endl;
           W_index = mom_index;
           j= -150 ;//make it go out of loop

        }
        else{

//           cout<<"in loop trace nu's mother"<<endl;
           grandmom_index = genMo1[ mom_index];
           mom_index = grandmom_index;
           mom_id = genParId[mom_index];
           
        }
        if(j>20){j=-200;}
        //to end the infinite loop if cannot find W after try 20 times

     } 

//     cout<<"check W find after out of loop?"<<endl;
//     cout<<"W_index: "<<W_index<<endl;
//     cout<<"gen id of W_index: "<<genParId[ W_index ]<<endl;
     
     if ( genParId[ genMo1[ W_index ] ] != 6 && genParId[ genMo1[ W_index ] ] != -6 ){continue;}// require W <- top 

//     cout<<"W's mother is top?"<<endl;
//     cout<<"genParId[ genMo1[ W_index ] ]: "<<genParId[ genMo1[ W_index ] ]<<endl;

//        cout<<"it is leptonic neutrino" <<endl;
//     cout<<"gen paricle index i: "<< i<<endl;                             
//     cout<<"genParId[i]: "<<genParId[i] <<endl;
//     cout<<"genParSt[i]: "<<genParSt[i] <<endl;

        // if it can go to this step, after require all requirement, then save in gen_leptonic_nu_index
        gen_leptonic_nu_index->push_back(i); 
     
  }

  //-----------------------------------------------------------------------------------//
  // finding muon neutrino <- W boson <- top quark   

  // old version to find status=3 neutrino

/*
  for(Int_t i = 0; i < nGenPar; i++){
     if (genParId[i]!=14 && genParId[i]!= -14){continue;}

//     if (genParSt[i]!=3){continue;}
     // status=3 nu's mother is W, grandmother is top 

//     if (genMomParId[i]!=24 && genMomParId[i]!= -24){continue;}
     
//     if ( genParId[ genMo1[ genMo1[i] ] ] != 6 && genParId[ genMo1[ genMo1[i]] ] != -6 ){continue;} 

        cout<<"gen paricle index i: "<< i<<endl;
        cout<<"genParId[i]: "<<genParId[i] <<endl;
        cout<<"genParSt[i]: "<<genParSt[i] <<endl;
        cout<<"genMomParId[i]: "<<genMomParId[i] <<endl;
        cout<<"genMo1[i]: "<<genMo1[i] <<endl;
        cout<<"genMo1[ genMo1[i]]: "<<genMo1[ genMo1[i] ] <<endl;
        cout<<"genParId[ genMo1[ genMo1[i]] ]: "<<genParId[ genMo1[ genMo1[i]] ] <<endl;
//        cout<<"it is leptonic neutrino" <<endl;

        gen_leptonic_nu_index->push_back(i); 

     
  }
*/

  //-----------------------------------------------------------------------------------//



}
