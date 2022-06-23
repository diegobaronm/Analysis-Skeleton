// Header guard to ensure file is imported properly
#ifndef Analysis
#define Analysis

// Include the file that lets the program know about the data
#include "backend/CLoop.h"
#include <iostream>
#include <vector>
#include <algorithm>
//#include <bits/stdc++.h>
#include <utility>

const int run2015Begin = 276262;
const int run2015End   = 284484;

const int run2016Begin = 297730;
const int run2016End   = 311481;

const int run2017Begin = 323427;
const int run2017End   = 341649;

const int run2018Begin = 341649;
const int run2018End   = 364292;


double del_phi(double phi_1, double phi_2){
    double pi=TMath::Pi();
    double phi_1_norm, phi_2_norm;
    if (phi_1<0.0){
        phi_1_norm=phi_1+2*pi;
    }else {
        phi_1_norm=phi_1;
    }

    if (phi_2<0.0){
        phi_2_norm=phi_2+2*pi;
    }else {
        phi_2_norm=phi_2;
    }
    double delta=std::abs(phi_1_norm-phi_2_norm);
    if (delta>pi){
        delta=2*pi-delta;
        delta=std::abs(delta);
    }

    return delta;
}

int is_inside_jets(TLorentzVector * test_jet,TLorentzVector * j1, TLorentzVector * j2){
  double delta_y_j1j2=abs(j1->Rapidity()-j2->Rapidity());
  double delta_y_j1test=abs(j1->Rapidity()-test_jet->Rapidity());
  double delta_y_j2test=abs(j2->Rapidity()-test_jet->Rapidity());
  if(delta_y_j1test>delta_y_j1j2 || delta_y_j2test>delta_y_j1j2){return 0;}
  else{return 1;}
}


void CLoop::Book(double lumFactor) {
    double pi=TMath::Pi();

    h_lep_pt_basic = new TH1F("lep_pt_basic","pT of the light lepton",200,0,200);
    h_lep_pt_basic_dphi = new TH1F("lep_pt_basic_dphi","pT of the light lepton",200,0,200);
    h_lep_pt_basic_dphi_btag = new TH1F("lep_pt_basic_dphi_btag","pT of the light lepton",200,0,200);
    h_lep_pt_basic_dphi_btag_iso = new TH1F("lep_pt_basic_dphi_btag_iso","pT of the light lepton",200,0,200);
    h_lep_pt_basic_dphi_btag_iso_rnn = new TH1F("lep_pt_basic_dphi_btag_iso_rnn","pT of the light lepton",200,0,200);
    h_lep_pt_basic_dphi_btag_iso_rnn_ptl = new TH1F("lep_pt_basic_dphi_btag_iso_rnn_ptl","pT of the light lepton",200,0,200);
    h_lep_pt_basic_dphi_btag_iso_rnn_ptl_nji = new TH1F("lep_pt_basic_dphi_btag_iso_rnn_ptl_nji","pT of the light lepton",200,0,200);
    h_lep_pt_basic_dphi_btag_iso_rnn_ptl_nji_mreco = new TH1F("lep_pt_basic_dphi_btag_iso_rnn_ptl_nji_mreco","pT of the light lepton",200,0,200);
    h_lep_pt_basic_dphi_btag_iso_rnn_ptl_nji_mreco_tpt = new TH1F("lep_pt_basic_dphi_btag_iso_rnn_ptl_nji_mreco_tpt","pT of the light lepton",200,0,200);

    h_delta_phi = new TH1F("delta_phi","Delta phi between tau and lep",32,0,3.2);
    h_n_bjets = new TH1F("n_bjets","Number of b_jets",5,0,5);
    h_lepiso = new TH1F("lepiso","Lep Isolation",2,0,2);
    h_rnn_score_1p = new TH1F("rnn_score_1p","RNN Score 1 prong taus",100,0,1);
    h_rnn_score_3p = new TH1F("rnn_score_3p","RNN Score 3 prong taus",100,0,1);
    h_lep_pt = new TH1F("lep_pt","Lep pT",500,0,500);
    h_tau_pt = new TH1F("tau_pt","Tau pT",500,0,500);
    h_n_jets_interval = new TH1F("n_jets_interval","N jets between rapidity interval",10,0,10);
    h_reco_mass = new TH1F("reco_mass","Reconstructed mass all events",240,0,240);

    // MC Exclusive histogram
    if (lumFactor!=1){
      h_Z_pt_truth_basic = new TH1F("Z_pt_truth_i_basic","Truth ZpT in between events",1000,0,1000);
      h_Z_pt_truth_basic_cuts = new TH1F("Z_pt_truth_i_basic_cuts","Truth ZpT in between events",1000,0,1000);
      h_Z_pt_truth_basic_cuts_tpt = new TH1F("Z_pt_truth_i_basic_cuts_tpt","Truth ZpT in between events",1000,0,1000);

    }

}

void CLoop::Fill(double weight, int z_sample) {
  double pi=TMath::Pi();
  //Charges
  float ql=muon_0_q;
  float qtau=tau_0_q;

  if (ql!=qtau && n_muons==1 && n_taus_rnn_loose>=1){
    
    //angles
    double angle_l_MET=del_phi(muon_0_p4->Phi(),met_reco_p4->Phi());
    double angle_tau_MET=del_phi(tau_0_p4->Phi(),met_reco_p4->Phi());
    double angle=del_phi(tau_0_p4->Phi(),muon_0_p4->Phi());
    //trigger decision
    bool trigger_decision= false;
    bool trigger_match= false;
    if (run_number>= 276262 && run_number<=284484) {
      trigger_decision= bool(HLT_mu20_iloose_L1MU15 | HLT_mu50);
      trigger_match=bool(muTrigMatch_0_HLT_mu20_iloose_L1MU15 | muTrigMatch_0_HLT_mu50);
    } else {
      trigger_decision= bool(HLT_mu26_ivarmedium | HLT_mu50);
      trigger_match=bool(muTrigMatch_0_HLT_mu26_ivarmedium | muTrigMatch_0_HLT_mu50);
    }

    if (trigger_decision  && trigger_match) {

      // RECO mass AND neutrino momentum
      double reco_mass{};

      reco_mass=sqrt(2*(tau_0_p4->Dot(*muon_0_p4)));
  


      // ZpT calculations
      double Z_pt_x=0;
      double Z_pt_y=0;
      double Z_pt=0;
      double truth_z_pt=0.0;

      // truth ZpT definition
      if (z_sample==1 || z_sample==2)
      {
        truth_z_pt=truth_Z_p4->Pt()/1000;
      }

      Z_pt_x=tau_0_p4->Pt()*cos(tau_0_p4->Phi())+muon_0_p4->Pt()*cos(muon_0_p4->Phi());
      Z_pt_y=tau_0_p4->Pt()*sin(tau_0_p4->Phi())+muon_0_p4->Pt()*sin(muon_0_p4->Phi());
      Z_pt=sqrt(Z_pt_x*Z_pt_x+Z_pt_y*Z_pt_y);
      if (z_sample==0){
        truth_z_pt=Z_pt;
      }
  
      


      // NUMBER OF JETS INTERVAL
      int n_jets_interval{};
      n_jets_interval=n_jets_interval+is_inside_jets(ljet_2_p4,ljet_0_p4,ljet_1_p4);
      
      
      // Cuts vector
      vector<int> cuts={0,0,0,0,0,0,0,0};
      // CUTS
      if (angle<=3.2){cuts[0]=1;}
      if(n_bjets_MV2c10_FixedCutBEff_85==0){cuts[1]=1;}
      if(muon_0_iso_TightTrackOnly_FixedRad==1){cuts[2]=1;}
      if(tau_0_n_charged_tracks==1 && tau_0_jet_rnn_score_trans>=0.25){cuts[3]=1;}
      if(tau_0_n_charged_tracks==3 && tau_0_jet_rnn_score_trans>=0.40){cuts[3]=1;}
      if(muon_0_p4->Pt()>=27){cuts[4]=1;}
      if(n_jets_interval==0){cuts[5]=1;}
      if (reco_mass<116 && reco_mass>66){cuts[6]=1;}
      if (tau_0_p4->Pt()>=25){cuts[7]=1;}

      // SUM OF THE VECTOR STORING IF CUTS PASS OR NOT
      int sum{};
      for(auto &j : cuts){sum=sum+j;}

      // FILLING CUTS HISTOGRAMS
      if ((sum-cuts[0])==7) {
        h_delta_phi->Fill(angle,weight);
      }
      if ((sum-cuts[1])==7) {
        h_n_bjets->Fill(n_bjets_MV2c10_FixedCutBEff_85,weight);
      }
      if ((sum-cuts[2])==7) {
        h_lepiso->Fill(muon_0_iso_TightTrackOnly_FixedRad,weight);
      }
      if ((sum-cuts[3])==7) {
        if (tau_0_n_charged_tracks==1){
          h_rnn_score_1p->Fill(tau_0_jet_rnn_score_trans,weight);
        }
        if (tau_0_n_charged_tracks==3){
          h_rnn_score_3p->Fill(tau_0_jet_rnn_score_trans,weight);
        }
      }
      if ((sum-cuts[4])==7) {
        h_lep_pt->Fill(muon_0_p4->Pt(),weight);
      }
      if ((sum-cuts[5])==7) {
        h_n_jets_interval->Fill(n_jets_interval,weight);
      }
      if ((sum-cuts[6])==7) {
        h_reco_mass->Fill(reco_mass,weight);
      }
      if ((sum-cuts[7])==7) {
        h_tau_pt->Fill(tau_0_p4->Pt(),weight);
      }


      // HISTOGRAM FILLING STARTING IN BASIC SELECTION
      if (weight!=1){
        h_Z_pt_truth_basic->Fill(truth_z_pt,weight);
      }




      // ANGLE CUT
      if (cuts[0]==1){

        h_lep_pt_basic_dphi->Fill(tau_0_p4->Pt(),weight);

        /// B TAG CUT
        if (cuts[1]==1){

          h_lep_pt_basic_dphi_btag->Fill(tau_0_p4->Pt(),weight);

          // ISOLATION CUT
          if (cuts[2]==1){

            h_lep_pt_basic_dphi_btag_iso->Fill(tau_0_p4->Pt(),weight);

            // JET RNN SCORE CUT
            if (cuts[3]==1){

              h_lep_pt_basic_dphi_btag_iso_rnn->Fill(tau_0_p4->Pt(),weight);

                // LEPTON PT CUT
              if (cuts[4]==1){

                h_lep_pt_basic_dphi_btag_iso_rnn_ptl->Fill(tau_0_p4->Pt(),weight);

                  // N JETS INTERVAL CUT
                if (cuts[5]==1){

                  h_lep_pt_basic_dphi_btag_iso_rnn_ptl_nji->Fill(tau_0_p4->Pt(),weight);

                  // RECO MASS CUT
                  if(cuts[6]==1){

                    h_lep_pt_basic_dphi_btag_iso_rnn_ptl_nji_mreco->Fill(tau_0_p4->Pt(),weight);
                    if (weight!=1){
                      h_Z_pt_truth_basic_cuts->Fill(truth_z_pt,weight);
                    }

                    //TAU PT CUT
                    if(cuts[7]==1){

                      h_lep_pt_basic_dphi_btag_iso_rnn_ptl_nji_mreco_tpt->Fill(tau_0_p4->Pt(),weight);
                      if (weight!=1){
                        h_Z_pt_truth_basic_cuts_tpt->Fill(truth_z_pt,weight);
                      }

                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}

 // Write your histograms to the output file
void CLoop::Style(double lumFactor) {
 

  h_lep_pt_basic->Write();
  h_lep_pt_basic_dphi->Write();
  h_lep_pt_basic_dphi_btag->Write();
  h_lep_pt_basic_dphi_btag_iso->Write();
  h_lep_pt_basic_dphi_btag_iso_rnn->Write();
  h_lep_pt_basic_dphi_btag_iso_rnn_ptl->Write();
  h_lep_pt_basic_dphi_btag_iso_rnn_ptl_nji->Write();
  h_lep_pt_basic_dphi_btag_iso_rnn_ptl_nji_mreco->Write();
  h_lep_pt_basic_dphi_btag_iso_rnn_ptl_nji_mreco_tpt->Write();


  h_delta_phi->Write();
  h_n_bjets->Write();
  h_lepiso->Write();
  h_rnn_score_1p->Write();
  h_rnn_score_3p->Write();
  h_lep_pt->Write();
  h_n_jets_interval->Write();
  h_reco_mass->Write();
  h_tau_pt->Write();

  if (lumFactor!=1){
    h_Z_pt_truth_basic->Write();
    h_Z_pt_truth_basic_cuts->Write();
    h_Z_pt_truth_basic_cuts_tpt->Write();
  }

}

#endif // End header guard