#include <map>
#include <vector>
#include <string>
#include <iostream>
#include <algorithm>
#include <TF1.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TBox.h>
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
#include "ElectronLooseID.C"
#include "MuonTightID.C"
#include "CA8JetLooseID.C"
#include "EventSelectionMuonChannel.C"
#include "neutrino_pz.C"
#include "gen_leptonic_nu.C"

using namespace std;

//Double_t fitFunc(Double_t*, Double_t*);
//void reconstructZPrime(TreeReader&, Double_t*, Double_t*);
void ElectronLooseID(TreeReader &, vector<Int_t> *  );
void MuonTightID(TreeReader &, vector<Int_t> *  );
void CA8JetLooseID(TreeReader &, vector<Int_t> *  );
int EventSelectionMuonChannel(TreeReader &, vector<Int_t> *,
 vector<Int_t> *, vector<Int_t> *, vector<Int_t> *,  vector<Int_t> *,
 vector<Int_t> * );
int neutrino_pz(TreeReader &, vector<Int_t> *, vector<Float_t> *, vector<Float_t> *);
void gen_leptonic_nu(TreeReader &,vector<Int_t> *);


//void ZPrimeMass(std::string inputFile){
void TTbarselect(){

//  TreeReader data(inputFile.data());
//  TreeReader data("flattuple_12_1_pNo.root");
  TreeReader data("combine_130_files.root");

  // Declare the histogram of reconstructed mass

  TH1D* h_genParId = new TH1D("h_genParId", "gen particle ID",60 , -30,30 );
  TH1D* h_elePt = new TH1D("h_elePt", "reco Ele Pt no Cuts",20 , 0,200 );
  TH1D* h_muPt = new TH1D("h_muPt", "reco Mu Pt no Cuts",20 , 0,200 );
  TH1D* h_CA8jetPt = new TH1D("h_CA8jetPt", "reco Jet Pt no Cuts",50 , 0,500 );
  TH1D* h_pfMetCorrPt = new TH1D("h_pfMetCorrPt", "reco MET Pt no Cuts",40 , 0,400 );
  TH1D* h_LooseCA8jetPt = new TH1D("h_LooseCA8jetPt", "reco Loose CA8 Jet Pt ",50 , 0,500 );
  TH1D* h_nCA8jet_pass_Evt_cut = new TH1D("h_nCA8jet_pass_Evt_cut", "# of CA8jet pass event cuts",10 , 0,10 );

  TH1D* h_leptonic_top_mass = new TH1D("h_leptonic_top_mass", "leptonic top mass",100 , 0,1000 );
  TH1D* h_hadronic_top_mass = new TH1D("h_hadronic_top_mass", "hadronic top mass",100 , 0,1000 );

//  TH1D* h_deltaR_genNu_recoNu1 = new TH1D("h_deltaR_genNu_recoNu1", "deltaR btw genNu and recoNu of positive pz sol",50 ,0 ,5 );
//  TH1D* h_deltaR_genNu_recoNu2 = new TH1D("h_deltaR_genNu_recoNu2", "deltaR btw genNu and recoNu of negative pz sol",50 ,0 ,5 );
  TH1D* h_deltaR_genNu_recoNu_smaller = new TH1D("h_deltaR_genNu_recoNu_smaller", "choose smaller deltaR of two nu pz sols with gen nu",50 ,0 ,5 );

  TH1D* h_deltaR_genNu_recoLep = new TH1D("h_deltaR_genNu_recoLep", " deltaR BTW genNu and recoLep",50 ,0 ,5 );



int event_pass_counter=0;
double counter1=0;
double counter2=0;
double counter0=0;
double counter3=0;

//--------------------------------------------------------------------//

  // begin of event loop
  
//  for (Long64_t ev = 0; ev < data.GetEntriesFast(); ev++){
  for (Long64_t ev = 0; ev < 100000; ev++){

    if ( ev % 5000 == 0 )
      fprintf(stderr, "Processing event %lli of %lli\n", ev + 1, data.GetEntriesFast());

//        cout<<"ev event: "<<ev<<endl;

    data.GetEntry(ev);
    

    Int_t    nGenPar   = data.GetInt("nGenPar"); 
    Int_t*   genParId  = data.GetPtrInt("genParId");
    Int_t*   genParSt  = data.GetPtrInt("genParSt");
    Float_t* genParPt  = data.GetPtrFloat("genParPt");
    Float_t* genParEta = data.GetPtrFloat("genParEta");
    Float_t* genParPhi = data.GetPtrFloat("genParPhi");
    Float_t* genParM   = data.GetPtrFloat("genParM");
    Float_t* genParE     = data.GetPtrFloat("genParE");
    Int_t*   genMomParId = data.GetPtrInt("genMomParId");
    Int_t*   genMo1      = data.GetPtrInt("genMo1");


    Int_t    nEle      = data.GetInt("nEle");
    Float_t* elePt     = data.GetPtrFloat("elePt");
    Float_t* eleScEta  = data.GetPtrFloat("eleScEta");
    Float_t* eleSigIhIh= data.GetPtrFloat("eleSigIhIh");

    Int_t    nMu          = data.GetInt("nMu");
    Float_t* muPt         = data.GetPtrFloat("muPt");
    Float_t* muEta        = data.GetPtrFloat("muEta");
    Float_t* muPhi        = data.GetPtrFloat("muPhi");
    Float_t* muM          = data.GetPtrFloat("muM");
    Int_t*   isGlobalMuon = data.GetPtrInt("isGlobalMuon");
    Int_t*   isTrackerMuon= data.GetPtrInt("isTrackerMuon");
    Int_t*   muTrkLayers  = data.GetPtrInt("muTrkLayers");
    Int_t*   muPixelHits  = data.GetPtrInt("muPixelHits");
    Int_t*   muHits       = data.GetPtrInt("muHits");
    Int_t*   muMatches    = data.GetPtrInt("muMatches");
    Float_t* mudxy        = data.GetPtrFloat("mudxy");
    Float_t* mudz         = data.GetPtrFloat("mudz");

    Int_t    CA8nJet      = data.GetInt("CA8nJet");
    Float_t* CA8jetPt     = data.GetPtrFloat("CA8jetPt");
    Float_t* CA8jetEta    = data.GetPtrFloat("CA8jetEta");
    Float_t* CA8jetPhi    = data.GetPtrFloat("CA8jetPhi");
    Float_t* CA8jetEn     = data.GetPtrFloat("CA8jetEn");
    Float_t* CA8jetMass   = data.GetPtrFloat("CA8jetMass");
    Float_t* CA8jetCEmEF  = data.GetPtrFloat("CA8jetCEmEF");
    Float_t* CA8jetCHadEF = data.GetPtrFloat("CA8jetCHadEF");
    Float_t* CA8jetNEmEF  = data.GetPtrFloat("CA8jetNEmEF");
    Float_t* CA8jetNHadEF = data.GetPtrFloat("CA8jetNHadEF");
    Float_t* CA8jetCMulti = data.GetPtrFloat("CA8jetCMulti");

    Float_t  pfMetCorrPt  = data.GetFloat("pfMetCorrPt");
    Float_t  pfMetCorrPhi = data.GetFloat("pfMetCorrPhi");

    //-----------------------------------------------------------------------------------// gen particle ID, 

        for(int i=0 ;i<nGenPar ;i++){

          h_genParId->Fill(genParId[i]);}

    //-----------------------------------------------------------------------------------// reco elePt, muPt, CA8jetPt, Met

        for(int i=0 ;i<nEle ;i++){

          h_elePt->Fill(elePt[i]);}

           for(int i=0 ;i<nMu ;i++){

          h_muPt->Fill(muPt[i]);}

        for(int i=0 ;i<CA8nJet ;i++){

          h_CA8jetPt->Fill(CA8jetPt[i]);}

          h_pfMetCorrPt->Fill(pfMetCorrPt);

    //-----------------------------------------------------------------------------------//
    vector<Int_t> Loose_Electron_index;
    ElectronLooseID(data, &Loose_Electron_index );

    vector<Int_t> Tight_Muon_index;
    MuonTightID(data, &Tight_Muon_index );

    vector<Int_t> Loose_CA8Jet_index;
    CA8JetLooseID(data, &Loose_CA8Jet_index );

      for( int i=0 ;i<Loose_CA8Jet_index.size() ;i++){
h_LooseCA8jetPt->Fill( CA8jetPt[ Loose_CA8Jet_index[i] ] );
     }

    //-----------------------------------------------------------------------------------//
    // event selection in muon channel

    vector<Int_t> Muon_pass_index;
//    int *leading_CA8Jet_index;
//    int *second_CA8Jet_index;
    vector<Int_t> leading_and_second_CA8Jet_index;
    vector<Int_t> Btagged_CA8Jet_index;

    int flag1 = -2;// muon channel event selection flag

    flag1 = EventSelectionMuonChannel(data, &Tight_Muon_index, &Loose_Electron_index ,
&Loose_CA8Jet_index, &Muon_pass_index, &leading_and_second_CA8Jet_index,
&Btagged_CA8Jet_index  );
    // flag1 = 2  for pass, and -2 and -1 for fail 


/*    if (flag1==2){
        cout<<"flag1:"<<flag1<<endl;
        event_pass_counter = event_pass_counter + 1 ; 
        h_nCA8jet_pass_Evt_cut->Fill( Loose_CA8Jet_index.size() );
        cout<<"Muon_pass_index size: "<<Muon_pass_index.size()<<endl;
        cout<<"leading CA8jet index: "<<leading_and_second_CA8Jet_index[0]
            <<" second CA8jet index: "<<leading_and_second_CA8Jet_index[1]<<endl;
        cout<<"Btagged_CA8Jet_index size: "<<Btagged_CA8Jet_index.size()<<endl;
        }*/

   //-----------------------------------------------------------------------------------//
   // calculate neutrino pz
   
   vector<Float_t> nu_pz;// neutrino pz 2 solutions   
   vector<Float_t> nu_eta;//neutrino eta of 2 solutions

   if (flag1==2){
      if (Muon_pass_index.size()==1){// present way only calculate nu pz in case of 1 muon 
        neutrino_pz(data, &Muon_pass_index, &nu_pz, &nu_eta);
//        cout<<"nu_pz size: "<<nu_pz.size()<<endl;
          if (nu_pz.size()>0){
//            cout<<"nu_pz_positive: "<<nu_pz[0]<<" nu_eta_positive: "<<nu_eta[0]<<endl;
//            cout<<"nu_pz_negative: "<<nu_pz[1]<<" nu_eta_negative: "<<nu_eta[1]<<endl;
          }
      }
   }

   //-----------------------------------------------------------------------------------// 
   // select gen leptonic neutrino

   // note, we are search in muon channel
    vector<Int_t> gen_leptonic_nu_index;

    gen_leptonic_nu(data, &gen_leptonic_nu_index);

/*    if (gen_leptonic_nu_index.size()>0){
         cout<<"gen_leptonic_nu_index.size(): "<< gen_leptonic_nu_index.size()<<endl;
       for (int i=0; i< gen_leptonic_nu_index.size();i++){

        cout<<"gen paricle index: "<< gen_leptonic_nu_index[i]<<endl;
        cout<<"genParId: "<<genParId[ gen_leptonic_nu_index[i] ] <<endl;
        cout<<"genMomParId: "<<genMomParId[gen_leptonic_nu_index[i]] <<endl;
        cout<<"genParId[ genMo1[ genMo1[i]] ]: "<<genParId[ genMo1[ genMo1[ gen_leptonic_nu_index[i]]] ] <<endl;

        cout<<"genParE: "<<genParE[ gen_leptonic_nu_index[i] ] <<endl;
        cout<<"genParPt: "<<genParPt[ gen_leptonic_nu_index[i] ] <<endl;
        cout<<"genParEta: "<<genParEta[ gen_leptonic_nu_index[i] ] <<endl;
        cout<<"genParPhi: "<<genParPhi[ gen_leptonic_nu_index[i] ] <<endl;

       }
    }
*/
   //-----------------------------------------------------------------------------------//
   //choose one of 2 nu pz has smaller delta R with gen nu 

   TLorentzVector gen_nu;
   TLorentzVector reco_nu0;//neutrino
   TLorentzVector reco_nu1;//neutrino
   float deltaR0;
   float deltaR1;
   int correct_nu_pz_index = -1;
   int flag_match =-1;

   TLorentzVector mu;
   int nu_pz_smaller_deltR_reco_lep_index = -1;

   // require gen_leptonic_nu_index.size()==1 for only one gen-nu, i.e ttbar semi-leptonic decay
   if( ( flag1==2 ) && ( nu_pz.size()!=0 ) && ( gen_leptonic_nu_index.size()==1 ) )
        { gen_nu.SetPtEtaPhiE(genParPt[ gen_leptonic_nu_index[0] ],
                              genParEta[ gen_leptonic_nu_index[0] ],
                              genParPhi[ gen_leptonic_nu_index[0] ],
                              genParE[ gen_leptonic_nu_index[0] ] );

          reco_nu0.SetPtEtaPhiM(pfMetCorrPt,nu_eta[0],pfMetCorrPhi,0);
          reco_nu1.SetPtEtaPhiM(pfMetCorrPt,nu_eta[1],pfMetCorrPhi,0);

          deltaR0 = gen_nu.DeltaR( reco_nu0 ); 
          deltaR1 = gen_nu.DeltaR( reco_nu1 );

          // choose the nu pz has smaller deltaR with gen nu 
          if( deltaR0 > deltaR1 ){ 
              correct_nu_pz_index = 1;
              if( deltaR1<0.2 ){flag_match=1;}
              h_deltaR_genNu_recoNu_smaller->Fill(deltaR1);
          }
          else if( deltaR0 < deltaR1 ){
              correct_nu_pz_index = 0;
              if( deltaR0<0.2 ){flag_match=1;}
              h_deltaR_genNu_recoNu_smaller->Fill(deltaR0);
          }

          // choose one of 2 nu pz has smaller delta R with reco lepton 
          // and see the correct rate of choose it 
          if(Muon_pass_index.size()==1){

              	  mu.SetPtEtaPhiM(   	muPt[Muon_pass_index[0]],
                	                muEta[Muon_pass_index[0]],
                        	        muPhi[Muon_pass_index[0]],
                                	muM[Muon_pass_index[0]]);


          	  deltaR0 = mu.DeltaR( reco_nu0 );
         	  deltaR1 = mu.DeltaR( reco_nu1 );

             	 // choose the nu pz has smaller deltaR with mu 
               if(flag_match==1){  
           	 if( deltaR0 > deltaR1 ){ 
               		 nu_pz_smaller_deltR_reco_lep_index  = 1;
              		  counter1=counter1+1;
                	//h_deltaR_genNu_recoNu_smaller->Fill(deltaR1);
           	 }
           	 else if( deltaR0 < deltaR1 ){
                	nu_pz_smaller_deltR_reco_lep_index = 0;
               		 counter1=counter1+1;
              		  //h_deltaR_genNu_recoNu_smaller->Fill(deltaR0);
           	 }
            	 else if(deltaR0 == deltaR1) {counter3=counter3+1;}
 
       	         if(correct_nu_pz_index == nu_pz_smaller_deltR_reco_lep_index)
           	 {counter2=counter2+1;}
                }//flag_match

            	//deltaR BTW gen nu and reco mu
            	gen_nu.SetPtEtaPhiE(    genParPt[ gen_leptonic_nu_index[0] ],
                	                genParEta[ gen_leptonic_nu_index[0] ],
                        	        genParPhi[ gen_leptonic_nu_index[0] ],
                                	genParE[ gen_leptonic_nu_index[0] ] );

           	 h_deltaR_genNu_recoLep->Fill( mu.DeltaR( gen_nu ) );

            	counter0=counter0+1;

          }//require muon size==0 

        }
   //-----------------------------------------------------------------------------------//
   // choose one of 2 nu pz has smaller delta R with reco lepton 
   // and see the correct rate of choose it 

//   TLorentzVector mu;
//   int nu_pz_smaller_deltR_reco_lep_index = -1;

   // require gen_leptonic_nu_index.size()==1 for only one gen-nu, i.e ttbar semi-leptonic decay
   //the correct_nu_pz_index is defined by require gen_leptonic_nu_index.size()==1
   //require Muon_pass_index.size()==1 aviod error
/*   if( ( flag1==2 ) && ( nu_pz.size()!=0 ) && ( gen_leptonic_nu_index.size()==1 ) && (Muon_pass_index.size()==1) )
        {
             
             mu.SetPtEtaPhiM(muPt[Muon_pass_index[0]],
                             muEta[Muon_pass_index[0]],
                             muPhi[Muon_pass_index[0]],
                             muM[Muon_pass_index[0]]);          


          deltaR0 = mu.DeltaR( reco_nu0 );
          deltaR1 = mu.DeltaR( reco_nu1 );

          // choose the nu pz has smaller deltaR with mu 
          if( deltaR0 > deltaR1 ){ 
              nu_pz_smaller_deltR_reco_lep_index  = 1;
              counter1=counter1+1;
              //h_deltaR_genNu_recoNu_smaller->Fill(deltaR1);
          }
          else if( deltaR0 < deltaR1 ){
              nu_pz_smaller_deltR_reco_lep_index = 0;
              counter1=counter1+1;
              //h_deltaR_genNu_recoNu_smaller->Fill(deltaR0);
          }
          else if(deltaR0 == deltaR1) {counter3=counter3+1;}
 
          if(correct_nu_pz_index == nu_pz_smaller_deltR_reco_lep_index)
          {counter2=counter2+1;}

          //deltaR BTW gen nu and reco mu
          gen_nu.SetPtEtaPhiE(genParPt[ gen_leptonic_nu_index[0] ],
                              genParEta[ gen_leptonic_nu_index[0] ],
                              genParPhi[ gen_leptonic_nu_index[0] ],
                              genParE[ gen_leptonic_nu_index[0] ] );

          h_deltaR_genNu_recoLep->Fill( mu.DeltaR( gen_nu ) );

          counter0=counter0+1;

        }
*/
   //-----------------------------------------------------------------------------------//
   // reconstruct leptonic top

//   TLorentzVector mu;// defined above
   TLorentzVector nu;//neutrino
   TLorentzVector jet;//CA8jet
   TLorentzVector lep_top;

   // pass event selection
   if (flag1==2 && nu_pz.size()!=0 ){// I don't deal with the complex number nu pz

//        cout<<"ev event: "<<ev<<endl;  
//        cout<<"leptonic top"<<endl;

        for (int i=0;i<Muon_pass_index.size();i++){// Muon
            for (int j=0;j<Loose_CA8Jet_index.size();j++){// Jet
                if(Loose_CA8Jet_index[j]==leading_and_second_CA8Jet_index[0]){continue;}
                // remove the leading jet to leptonic top reco

                   for (int k=0; k<Btagged_CA8Jet_index.size();k++ ){
                      if(Loose_CA8Jet_index[j]==Btagged_CA8Jet_index[k]){
                      // require the b-tagged jet

                         for (int l=0;l<nu_pz.size();l++ ){// eta from nu pz 2 solutions
                       
                      mu.SetPtEtaPhiM(muPt[Muon_pass_index[i]],muEta[Muon_pass_index[i]],
                                      muPhi[Muon_pass_index[i]],muM[Muon_pass_index[i]]);

//                      jet.SetPtEtaPhiE(CA8jetPt[Loose_CA8Jet_index[j]],CA8jetEta[Loose_CA8Jet_index[j]],
//                                       CA8jetPhi[Loose_CA8Jet_index[j]],CA8jetEn[Loose_CA8Jet_index[j]]);

                      jet.SetPtEtaPhiM(CA8jetPt[Loose_CA8Jet_index[j]],CA8jetEta[Loose_CA8Jet_index[j]],
                                       CA8jetPhi[Loose_CA8Jet_index[j]],CA8jetMass[Loose_CA8Jet_index[j]]);


                      nu.SetPtEtaPhiM(pfMetCorrPt,nu_eta[l],pfMetCorrPhi,0);

                      lep_top = mu + jet +nu ;
//                      cout<<"i muon: "<<i<<" j jet: "<<j<<" k b-tagged jet: "<<k<<" l nu: "<<l<<endl;
//                      cout<<"top mass:"<<lep_top.M()<<endl;
                      h_leptonic_top_mass->Fill( lep_top.M() );
                         }//nu l loop
                      }//b-tagging if 
                   }//b-taging k loop
            }//jet j loop
        }//muon i loop 
      }//flag & nu_pz size if 

   //-----------------------------------------------------------------------------------//
   // reconstruct hadronic top

   TLorentzVector W_jet;
   TLorentzVector b_jet;
   TLorentzVector had_top;

  if (flag1==2  ){//pass event selection

//        cout<<"ev event: "<<ev<<endl;
//        cout<<"hadronic top"<<endl;        

        W_jet.SetPtEtaPhiM(CA8jetPt[leading_and_second_CA8Jet_index[0]],CA8jetEta[leading_and_second_CA8Jet_index[0]],
                           CA8jetPhi[leading_and_second_CA8Jet_index[0]],CA8jetMass[leading_and_second_CA8Jet_index[0]]);

      

        for (int k=0; k<Btagged_CA8Jet_index.size();k++ ){//require the b-tagged jet
            if( Btagged_CA8Jet_index[k] == leading_and_second_CA8Jet_index[0] ){continue;}
            // remove the leading jet to be b-jet
 
            b_jet.SetPtEtaPhiM(CA8jetPt[Btagged_CA8Jet_index[k]],CA8jetEta[Btagged_CA8Jet_index[k]],
                               CA8jetPhi[Btagged_CA8Jet_index[k]],CA8jetMass[Btagged_CA8Jet_index[k]]);


                      had_top = W_jet +  b_jet ;
//                      cout<<" k b-tagged jet: "<<k<<endl;
//                      cout<<"top mass:"<<had_top.M()<<endl;
                      h_hadronic_top_mass->Fill(had_top.M());
        }//end b-tagging k loop
  }//end if

   //-----------------------------------------------------------------------------------//

  }   // end of event loop

cout<<"counter1: "<<counter1<<endl;
cout<<"counter2: "<<counter2<<endl;
cout<<"counter3(same deltaR): "<<counter3<<endl;
cout<<"correct rate of choosing nu pz has smaller deltaR with reco lepton"<<endl;
cout<<"rate counter2/counter1: "<<counter2/counter1<<endl;

cout<<"counter0: "<<counter0<<endl;

//cout<<"event_pass_counter: "<<event_pass_counter<<endl;

/*    TCanvas * c1 =new TCanvas("c1","",600,600);
    h_genParId->Draw();
    TCanvas * c2 =new TCanvas("c2","",600,600);
    h_elePt->Draw();
    TCanvas * c3 =new TCanvas("c3","",600,600);
    h_muPt->Draw();
    TCanvas * c4 =new TCanvas("c4","",600,600);
    h_CA8jetPt->Draw();
    TCanvas * c5 =new TCanvas("c5","",600,600);
    h_pfMetCorrPt->Draw();
*/
    TCanvas * c6 =new TCanvas("c6","",600,600);
h_deltaR_genNu_recoNu_smaller->Draw();
//h_deltaR_genNu_recoLep->Draw();
//h_LooseCA8jetPt->Draw();
//h_nCA8jet_pass_Evt_cut->Draw();
//h_leptonic_top_mass->Draw();
//h_deltaR_genNu_recoNu1->Draw();
//    TCanvas * c7 =new TCanvas("c7","",600,600);
//h_hadronic_top_mass->Draw();
//h_deltaR_genNu_recoNu2->SetLineColor(2);
//h_deltaR_genNu_recoNu2->SetLineColor(kOrange+10);
//h_deltaR_genNu_recoNu2->Draw("same");
//    c1->Print("genParId.png");
//    c2->Print("elePt.png");
//    c3->Print("muPt.png");
//    c4->Print("CA8jetPt.png");
//    c5->Print("pfMetCorrPt.png");

  // Fitting the mass ratio
  
  
}

