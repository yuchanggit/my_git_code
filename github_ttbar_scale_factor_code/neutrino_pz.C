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
#include <TVector3.h>

using namespace std;

int neutrino_pz(TreeReader &data,vector<Int_t> *Muon_pass_index, 
                 vector<Float_t> *nu_pz , vector<Float_t> *nu_eta){

    Float_t* muPt         = data.GetPtrFloat("muPt");
    Float_t* muEta        = data.GetPtrFloat("muEta");
    Float_t* muPhi        = data.GetPtrFloat("muPhi");
    Float_t* muM          = data.GetPtrFloat("muM");

    Float_t  pfMetCorrPt  = data.GetFloat("pfMetCorrPt");
    Float_t  pfMetCorrPhi = data.GetFloat("pfMetCorrPhi");

  //-----------------------------------------------------------------------------------//
  // calculate neutrino pz  

  float W_mass = 80.3;// W mass = 80.3 GeV

  float inner_product_of_Pt_l_and_MET = 
      muPt[Muon_pass_index->at(0)]*pfMetCorrPt*cos(muPhi[Muon_pass_index->at(0)]-pfMetCorrPhi);
//cout<<"inner_product_of_Pt_l_and_MET: "<<inner_product_of_Pt_l_and_MET<<endl;

  float a = pow(W_mass,2) + 2*inner_product_of_Pt_l_and_MET;
//cout<<"a: "<<a<<endl;

  TLorentzVector mu;
  mu.SetPtEtaPhiM(muPt[Muon_pass_index->at(0)],muEta[Muon_pass_index->at(0)],
                  muPhi[Muon_pass_index->at(0)],muM[Muon_pass_index->at(0)]);

  float mu_pz = mu.Pz();// muon pz
  float mu_en = mu.E();//  muon energy  

  float MET = pfMetCorrPt;// neutrino Pt 

  // use solution of 1 unknow, 2nd order equation to calculate the nu pz
  float argument_inside_sqrt = 16*pow(a,2)*pow(mu_pz,2)-
4*( 4*pow(mu_en,2)-4*pow(mu_pz,2) )*(4*pow(mu_en*MET,2) - pow(a,2) );

  // positive and negative neutrino pz solutions
  float nu_pz_positive;
  float nu_pz_negative;

  // positive and negative neutrino eta solutions
  float nu_eta_positive;
  float nu_eta_negative;  

  float P_tot;

//  cout<<"argument_inside_sqrt: "<<argument_inside_sqrt<<endl;

  if (argument_inside_sqrt <0){return 0;}
  // nu pz calculate fail bcz it cannot be complex number
  else{
   
  nu_pz_positive = ( 4*a*mu_pz + sqrt(argument_inside_sqrt) )/
                  ( 2*( 4*pow(mu_en,2)-4*pow(mu_pz,2) ) ) ;
  nu_pz_negative = ( 4*a*mu_pz - sqrt(argument_inside_sqrt) )/
                   ( 2*( 4*pow(mu_en,2)-4*pow(mu_pz,2) ) ) ;

  nu_pz->push_back(nu_pz_positive);
  nu_pz->push_back(nu_pz_negative);

//cout<<"MET: "<< MET<<endl;

//  nu_eta_positive = - log( tan( atan(MET/nu_pz_positive) / 2 ) );
//  nu_eta_negative = - log( tan( atan(MET/nu_pz_negative) / 2 ) );

  P_tot =  sqrt( pow(MET,2) + pow(nu_pz_positive,2) ) ;
  nu_eta_positive = 0.5 * log( (P_tot+nu_pz_positive)/(P_tot-nu_pz_positive) );

  P_tot =  sqrt( pow(MET,2) + pow(nu_pz_negative,2) ) ;
  nu_eta_negative = 0.5 * log( (P_tot+nu_pz_negative)/(P_tot-nu_pz_negative) );

  nu_eta->push_back( nu_eta_positive );
  nu_eta->push_back( nu_eta_negative );

//cout<<"nu_pz_positive: "<<nu_pz_positive<<" nu_pz_negative: "<<nu_pz_negative<<endl;
//cout<<"nu_eta_positive: "<<nu_eta_positive<<" nu_eta_negative: "<<nu_eta_negative<<endl;

  }
  //-----------------------------------------------------------------------------------//
  // close check for nu pz

/*
  float W_mass_check_positive=
     sqrt( pow( ( mu_en + sqrt(pow(MET,2)+pow(nu_pz_positive,2)) ),2) - 
   ( pow( (mu.Px() + MET*cos(pfMetCorrPhi)),2)+pow( (mu.Py() + MET*sin(pfMetCorrPhi)),2)+
     pow( (mu.Pz() + nu_pz_positive) ,2) ) );

  float W_mass_check_negative=
     sqrt( pow( ( mu_en + sqrt(pow(MET,2)+pow(nu_pz_negative,2)) ),2) -
   ( pow( (mu.Px() + MET*cos(pfMetCorrPhi)),2)+pow( (mu.Py() + MET*sin(pfMetCorrPhi)),2)+ 
     pow( (mu.Pz() + nu_pz_negative) ,2) ) );

cout<<"W_mass_check_positive: "<<W_mass_check_positive <<endl
    <<"W_mass_check_negative: "<<W_mass_check_negative <<endl;
*/  

  //-----------------------------------------------------------------------------------//
  // close check for nu  eta

/*
  TLorentzVector Wboson;

  TLorentzVector nu;
  nu.SetPtEtaPhiM(MET,nu_eta_positive,pfMetCorrPhi,0);

  Wboson = mu + nu;
  cout<<"W mass for positive eat: "<<Wboson.M()<<endl;

//  TLorentzVector nu;
  nu.SetPtEtaPhiM(MET,nu_eta_negative,pfMetCorrPhi,0);

  Wboson = mu + nu;
  cout<<"W mass for negative eat: "<<Wboson.M()<<endl;
*/


}
