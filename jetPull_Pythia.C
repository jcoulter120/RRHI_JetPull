// Jennifer Coulter
// February 10th 2016
// Rutgers University, jennifer.coulter@cern.ch
//
// Test macro for plotting the jet pull vector. 
//

//CURRENTLY : This macro attempts to compute the jet pull vector of the leading jet, then calculates the pull angle of that jet with all subjets. 

#include <iostream>
#include <stdio.h>
#include <fstream>
#include <TH1F.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TFile.h>
#include <TTree.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TFile.h>
#include <TStyle.h>
#include <TStopwatch.h>
#include <TRandom3.h>
#include <TChain.h>
#include <TProfile.h>
#include <TStopwatch.h>
#include <TCut.h>
#include <cstdlib>
#include <cmath>
#include "TLegend.h"
#include "TLatex.h"
#include "TMath.h"
#include "TLine.h"
#include "TVector3.h"
#include <TROOT.h>

using namespace std;


// DIVIDE BY BIN WIDTH ==================================================
//creates histograms in terms of variable vs. dN/dT
TH1F* DivideByBinWidth(TH1F * hist, const char * name){

  TH1F* h_return = new TH1F(name, "", hist->GetNbinsX(), 0,1);
  hist->Sumw2(); 
  //loops through all the bins
  for (int i=1;i<=hist->GetNbinsX();++i){
    Float_t bin = hist->GetBinWidth(i);
    Float_t val = hist->GetBinContent(i);
    Float_t valErr = hist->GetBinError(i);
    val = val/bin;
    valErr= valErr/bin;
    h_return->SetBinError(i,valErr);
    h_return->SetBinContent(i, val); 
  }//end bin loop
  return h_return;
}//end rebin function ==================================================

//Function to normalize a vector =======================================
TVector3 Norm(TVector3 v){
  if ( (v(0) == 0) && (v(1) == 0) && (v(2) == 0)) return v; 
  Double_t mag = TMath::Sqrt(v(0)*v(0) + v(1)*v(1) + v(2)*v(2)); 
  v(0) = v(0)/mag;    v(1) = v(1)/mag;   v(2) = v(2)/mag;
  return v; 
}//end normalize =======================================================

//BEGIN MAIN MACRO =====================================================
void jetPull_Pythia(Int_t startfile = 0,
		   Int_t endfile = 44,
		   Int_t jobNumber = 3,
		   int radius = 3,
		   float ptCutLead = 30.0,
		   float ptCut = 30.0,
		   float etaCut = 2.0){

  TH1::SetDefaultSumw2();

  TStopwatch timer;
  bool debug = false;
  
  //DEFINE TREES AND FILES +++++++++++++++++++++++++++++++++++
  TFile * file;
  TFile * weight_file;
  TFile * save_File = new TFile(Form("pythia_jetPull_%d.root", jobNumber),"RECREATE");

  TTree * t;
  TTree * hiEvt;
  TTree * skim;
  TTree * weight;

  //DEFINE HISTOGRAMS +++++++++++++++++++++++++++++++++++
  TH1F * h_jetPull = new TH1F("jetPull", "", 50,0,1);  
  TH1F * h_leadingJetpT = new TH1F("leading_pT", "", 100, 0, 500);
  TH1F * h_pTcut = new TH1F("pTcut", "", 100, 0, 500);
  TH1F * h_nref = new TH1F("nref", "", 12, 0, 12);
  TH1F * h_jetCount = new TH1F("jetCount", "", 12, 0, 12);
  TH1F * h_leadingJeteta = new TH1F("leadingJeteta", "", 60, -2, 2);
  TH1F * h_leadingJetphi = new TH1F("leadingJetphi", "", 60, -3.15, 3.15);
  TH1F * h_weight = new TH1F("weighted", "", 1200, 0, 1200);
  TH1F * h_unweight = new TH1F("unweighted", "", 1200, 0, 1200);
  TH1F * h_etaParticle = new TH1F("etaParticle", "", 60, -3, 3);
  TH1F * h_yParticle = new TH1F("yParticle", "", 60, -3, 3);
  TH2F * h_jetPull_2d = new TH2F("jetPull_2d","Jet pT vs Eta for pT>20",30,-3,+3,60,-3.5,3.5);
  TH2D * h_jetPull_3d = new TH2D("jetPull_3d","phi vs y vs pT" ,30,-3,+3,60,-3.5,3.5);
  TH1F * h_30 = new TH1F("pTcutoff30", "Jet Pull Angle -- with pThat Cutoffs", 60, 0.01, 500);
  TH1F * h_60 = new TH1F("pTcutoff60", "Jet Pull Angle -- with pThat Cutoffs", 60, 0.01, 500);
  TH1F * h_90 = new TH1F("pTcutoff90", "Jet Pull Angle -- with pThat Cutoffs", 60, 0.01, 500);
  TH1F * h_120 = new TH1F("pTcutoff120", "Jet Pull Angle -- with pThat Cutoffs", 60, 0.01, 500);
  TH1F * h_pullAngle = new TH1F("pullAngle", "Jet Pull Angle vs pT", 60, -3.5, 3.5); 

  //DEFINE TREE VARIABLES +++++++++++++++++++++++++++++++++++
  Float_t pt[1000];  
  Float_t eta[1000];  
  Float_t phi[1000];  
  Float_t y[1000];
  Float_t genpt[1000];
  Float_t genphi[1000];
  Float_t geny[1000];
  Float_t geneta[1000];
  
  Int_t nref;
  Int_t ngen;
  Float_t vz;
  Int_t run; 
  Int_t lumi;
  Int_t evt = 0; 
  Int_t isGoodEvt = 0; 
  Int_t halo;
  Int_t noise;
  
  Double_t pThat_weight;
  Float_t pThat;

  //DEFINE VARS FOR CALCULATION +++++++++++++++++++++++++++++
  TVector2 r[1000];
  TVector2 jetPull;
  Int_t leadJet = 0;
  Int_t eventCount = 0;
  Float_t jetPullMag  = 0;
  Float_t pullAngle = 0;
  Float_t dipolarity = 0;

  //LOAD IN FILE NAMES +++++++++++++++++++++++++++++++++++
  string input_file = "Text_Files/pp_MC_HiForest.txt";
  ifstream count(input_file.c_str(), ifstream::in);
  Int_t fileCount = 0;
  string * filename = new string[135];
  
  string line;
  while(getline(count, line)){
    filename[fileCount] = line;
    if (debug) cout << filename[fileCount] << endl; 
    fileCount++;
  }
  count.close();

  //BEGIN FILE LOOP ======================================================
  for(int ifile = startfile; ifile < endfile; ifile++){
    
    string str = filename[ifile];
    file = TFile::Open(str.c_str());
    string w = Form("weights/weights_pp_%d.root", ifile+1);
    cout << w << endl; 
    weight_file = TFile::Open(w.c_str());

    if (debug) cout << "\n **** =========================== New File ================================= **** \n ";
    cout << "File Name: " << filename[ifile] << endl;
    cout << "File Number: " << ifile << "/" << startfile << " to " << endfile << endl;

    //define trees and file +++++++++++++++++
    t = (TTree*)file->Get(Form("ak%dPFJetAnalyzer/t", radius));
    return;
    hiEvt = (TTree*)file->Get("hiEvtAnalyzer/HiTree");
    skim = (TTree*)file->Get("skimanalysis/HltTree");
    weight = (TTree*)weight_file->Get("weights");

    //Set branches of the tree ++++++++++++++ 
    t->SetBranchAddress("jtpt", &pt);
    t->SetBranchAddress("jteta", &eta);
    t->SetBranchAddress("jtphi", &phi);
    t->SetBranchAddress("jty", &y);
    t->SetBranchAddress("nref", &nref);

    t->SetBranchAddress("ngen", &ngen);
    t->SetBranchAddress("genpt", &genpt);
    t->SetBranchAddress("geny", &geny);
    t->SetBranchAddress("genphi", &genphi);
    t->SetBranchAddress("geneta", &geneta);

    hiEvt->SetBranchAddress("vz", &vz);
    hiEvt->SetBranchAddress("evt", &evt);    
    hiEvt->SetBranchAddress("run", &run);    
    hiEvt->SetBranchAddress("lumi", &lumi);    

    skim->SetBranchAddress("pHBHENoiseFilter", &noise);
    skim->SetBranchAddress("pPAcollisionEventSelectionPA",&halo);

    t->SetBranchAddress("pthat", &pThat);
    weight->SetBranchAddress("pthatweight", &pThat_weight);

    t->AddFriend(hiEvt);
    t->AddFriend(skim);
    t->AddFriend(weight);
    
    Long64_t nentries = t->GetEntries();
    nentries = 1000;
    cout << "Events in File: " << nentries << endl;
    eventCount = 0;
    
    //EVENT LOOP  ================================================
    for(Long64_t nentry = 0; nentry<nentries; ++nentry){

      if(nentry%10000 == 0) cout << nentry << endl;
      
      t->GetEntry(nentry);
      skim->GetEntry(nentry);
      hiEvt->GetEntry(nentry);
      weight->GetEntry(nentry);
      
      if(TMath::Abs(vz) > 15 || halo == 0) {
	isGoodEvt = 0;
      	continue;
      }
      
      // get the pt vector for each event which passes you jet selections based on pT and eta
      vector <float> pt_v;
      vector <float> eta_v;
      vector <float> phi_v;
      Double_t ptLead = 0;
      
      for(int ij = 0; ij<nref; ++ij){
	if(pt[ij] > ptCutLead && fabs(eta[ij]) < etaCut) {
	  isGoodEvt = 1;
	  if(pt[ij] > ptLead) {
	    leadJet = ij;
	    ptLead = pt[ij];
	  }
	  if(debug) cout << "LEAD JET: "<<  pt[ij] << " *******************************" << endl;

	}
	if(pt[ij] > ptCut && fabs(eta[ij]) < etaCut) { 
	  pt_v.push_back(pt[ij]);	
	  eta_v.push_back(eta[ij]);	
	  phi_v.push_back(phi[ij]);
	}
      }

      if(debug) cout <<"Total Number of Jets    : "<<nref<<endl;
      if(debug) cout <<"Number of Selected Jets : "<<pt_v.size()<<endl;

      int NJets_Sel = pt_v.size();

      if(NJets_Sel < 2) {
      	if (debug) cout<<"This event had only 1 Jet"<<endl;
	continue;	
      }
      if(!isGoodEvt){  continue;  }

      if(debug) cout<< " \n ******* New Event ******** " << endl;
      if(debug) cout<< " ******* " << nentry << " ******** " << endl;
      
      //reset maximum values
      eventCount++;
      evt++;
      //EVENT LOOP  ====================================================
      //runs through all of the jets, and calculates pull angle with just the leading jet
      //for(Long64_t njet = 0; njet < NJets_Sel; ++njet){	

      //we will only need jet pull for the lead jet
      //if(njet != leadJet) {
      //  continue;
      //}
      
      r[leadJet].Set(0.0,0.0);
      jetPull.Set(0.0,0.0);
      if(debug) cout << "NEW JET: "<<  ngen << " *******************************" << endl;
      
      //PARTICLE LOOP STARTS HERE ===================================
      //calculate jet pull for the individual jet
      for(Long64_t npart = 0; npart < ngen; ++npart){
	
	//these files seem to have a lot of gen particles with pt = 0
	//if ( genpt[npart] == 0) {
	//  continue; 
	//}

	//if( genpt[npart] > pt_v[leadJet] ) {
	  //  continue; 
	  //	}
	
	//Delta R^2 cut
	if((Double_t)(TMath::Power(TMath::Power(phi_v[leadJet] - genphi[npart],2) + TMath::Power(eta_v[leadJet] - geneta[npart],2),2)) > 0.3){
	  continue;
	}
	
	//calculate jet pull
	TVector2 partAxis = TVector2(geneta[npart], genphi[npart]);
	TVector2 jetAxis = TVector2(eta_v[leadJet], phi_v[leadJet]);
	TVector2 newJetPull; 
	r[leadJet].Set((partAxis.X() - jetAxis.X()),(partAxis.Y() - jetAxis.Y()));
	//cout << r[leadJet].X() << " , " << r[leadJet].Y() << endl; 
	
	if(debug) cout << "NEW PART: Lead Jet pT " << pt_v[leadJet] << " Particle pT " << genpt[npart] << " ======> Lead Jet-->Subjet ri" <<  TMath::Sqrt(r[leadJet].X()*r[leadJet].X() + r[leadJet].Y()*r[leadJet].Y()) << endl;
	
	jetPullMag = TMath::Sqrt(r[leadJet].X()*r[leadJet].X() + r[leadJet].Y()*r[leadJet].Y()) * genpt[npart]/pt_v[leadJet];
	if(debug) cout << "JET PULL MAG: " << jetPullMag << endl;
	newJetPull.Set(jetPullMag*r[leadJet].X(), jetPullMag*r[leadJet].Y());
	jetPull.Set(jetPull.X() + newJetPull.X(), jetPull.Y() + newJetPull.Y());
 	if(debug) cout << "JET PULL MAG w/ NEW PARTICLE: " << jetPullMag << endl;
	
	h_yParticle->Fill(geneta[npart], pThat_weight);
	h_etaParticle->Fill(geneta[npart], pThat_weight);
	
      }//PARTICLE LOOP ENDS HERE ===========
      
      jetPullMag = TMath::Sqrt(jetPull.X()*jetPull.X() + jetPull.Y()*jetPull.Y());
      
      //only plot the leading jets
      if(jetPullMag != 0){
	h_jetPull->Fill(jetPullMag, pThat_weight);
	h_jetPull_2d->Fill(jetPull.X(), jetPull.Y());
	h_jetPull_3d->Fill(jetPull.X(), jetPull.Y(), pt_v[leadJet]);
	h_leadingJetpT->Fill(pt_v[leadJet], pThat_weight);
	h_leadingJetphi->Fill(phi_v[leadJet], pThat_weight);
      }

      //}//JET PULL LOOP ENDS HERE ===========================
      
      if(debug) cout << "Number of Jets Selected: " << NJets_Sel << endl;
      
      //calculate the pull angle of each subjet with the leading jet
      //PULL ANGLE LOOP ====================================
      for(Long64_t npull = 0; npull < NJets_Sel; ++npull){
	
	if(debug) cout << "lead jet : " << leadJet << " Current Jet " << npull <<" NSelected : " << NJets_Sel << " Jet Pull " << jetPullMag << endl;

	//prevents double counting of lead jet
	if(npull == leadJet) {
	  continue;
	}
	//prevents us from using jets will jet pull of zero --> why do we get these anyway?
	if(jetPullMag == 0) {
	  continue;
	}
	
	// the vector defining the difference between the lead jet and this particular jet
	TVector2 rpull; 
	rpull.Set(eta_v[leadJet]-eta_v[npull], phi_v[leadJet]-phi_v[npull]);
	
	pullAngle = TMath::ACos((jetPull.X()*rpull.X() + jetPull.Y()*rpull.Y())/(TMath::Sqrt(jetPull.X()*jetPull.X() + jetPull.Y()*jetPull.Y()) * TMath::Sqrt(rpull.X()*rpull.X() + rpull.Y()*rpull.Y())));	
	if(debug) cout << "pt_leadJet " << pt_v[leadJet] << " pt_mainJet " << pt_v[npull] << " Pull Angle: " << pullAngle << "\n arccos inside " << (jetPull.X()*rpull.X() + jetPull.Y()*rpull.Y()) << " / " << (TMath::Sqrt(jetPull.X()*jetPull.X() + jetPull.Y()*jetPull.Y()) * TMath::Sqrt(rpull.X()*rpull.X() + rpull.Y()*rpull.Y())) <<  endl;
	
	//fill pull angle and eta difference (for checking)

	h_pullAngle->Fill(pullAngle, pThat_weight);
	h_leadingJeteta->Fill(eta_v[leadJet] - eta_v[npull]);      
      } // END OF PULL ANGLE LOOP ===============================
      
      //fill values
      h_nref->Fill(nref);
      h_jetCount->Fill(NJets_Sel);
      h_weight->Fill(pThat,pThat_weight);
      h_unweight->Fill(pThat);     
      
      pt_v.clear();
      eta_v.clear();
      phi_v.clear();
      
    }//end of event loop ==========================
    
    gROOT->GetListOfFiles()->Remove(file);
    gROOT->GetListOfFiles()->Remove(weight);
    
    cout << "Events Selected: " << eventCount << endl;
    cout << "File Finished" << endl; 
    
  }//end file loop ================================
  
  //SCALE HISTOGRAMS +++++++++++++++++++++++++++++++++++  
  //CREATE PLOTS VS. dN/dT +++++++++++++++++++++++ (NOT CURRENTLY USED)
  TH1F * h_jetPullScaled = DivideByBinWidth(h_jetPull, "jetPull_scaled");

  save_File->cd();
  
  h_jetPullScaled->Write();
  h_jetPull->Write();
  h_leadingJetpT->Write();
  h_pTcut->Write();
  h_nref->Write();
  h_jetCount->Write();
  h_leadingJeteta->Write();
  h_leadingJetphi->Write();
  h_weight->Write();
  h_unweight->Write();
  h_yParticle->Write();
  h_etaParticle->Write();
  h_jetPull_2d->Write();
  h_jetPull_3d->Write();
  h_pullAngle->Write();
  h_30->Write();
  h_60->Write();
  h_90->Write();
  h_120->Write();
  
  save_File->Write();
  save_File->Close();
  gROOT->GetListOfFiles()->Remove(save_File);

  timer.Stop();
  cout<<"Macro finished: "<<endl;
  
  cout<<"CPU time (min)  = "<<(Float_t)timer.CpuTime()/60<<endl;
  cout<<"Real time (min) = "<<(Float_t)timer.RealTime()/60<<endl;
    
}//end of plot jetPull

 
