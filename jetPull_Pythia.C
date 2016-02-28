// Jennifer Coulter
// February 10th 2016
// Rutgers University, jennifer.coulter@cern.ch
//
// Test macro for plotting the jet pull vector. 
//

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

//PLANE CLASS =======================================================
class Plane{
public:
  TVector3 v1, v2, proj, u1, u2;
  Double_t scalar1, scalar2, mag1, mag2; 
  Plane(TVector3);

  //returns a projection onto the 2D plane 
  TVector3 Projection(TVector3 jaxis){
    //Find the projection of a jet onto this subspace
    if(v1.Mag() == 0) { scalar1 = 0; }   else { scalar1 = jaxis.Dot(v1)/(v1.Dot(v1)); } 
    if(v2.Mag() == 0) { scalar2 = 0; }   else { scalar2 = jaxis.Dot(v2)/(v2.Dot(v2)); } 
    v1 = scalar1*v1;
    v2 = scalar2*v2;
    proj(0) = v1(0) + v2(0);
    proj(1) = v1(1) + v2(1);
    proj(2) = v1(2) + v2(2); 
    
    return proj;
  }//end of projection
};
//plane class constructor
Plane::Plane(TVector3 nT){
  
  //Use TVector3 to find an orthogonal vector and a second vector orthogonal to the first and nT
  v1 = nT.Orthogonal();  v2 = nT.Cross(v1);

  //Normalize, checking for 0 length axes
  if ((v1(0) == 0) && (v1(1) == 0) && (v1(2) == 0)){  v1(0) = 0;    v1(1) = 0;    v1(2) = 0; }
  else { mag1 = v1.Mag();   v1(0) = v1(0)/mag1;    v1(1) = v1(1)/mag1;    v1(2) = v1(2)/mag1; } 
  if ((v2(0) == 0) && (v2(1) == 0) && (v2(2) == 0)){  v2(0) = 0;    v2(1) = 0;    v2(2) = 0; } 
  else { mag2 = v2.Mag();   v2(0) = v2(0)/mag2;    v2(1) = v2(1)/mag2;    v2(2) = v2(2)/mag2; }
}//end plane constructor
//END PLANE CLASS =======================================================

// DIVIDE BY BIN WIDTH ==================================================
//creates histograms in terms of Thrust vs. dN/dT
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
		   Int_t endfile = 22,
		   Int_t jobNumber = 100,
		   int radius = 3,
		   float ptCutLead = 60.0,
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
  TH1F * h_thrust = new TH1F("thrust_unscaled", "", 50,0,1);
  TH1F * h_min = new TH1F("thrust_min", "", 50,0,1);
  TH1F * h_maj = new TH1F("thrust_maj", "", 50,0,1);
  
  TH1F * h_jetPull = new TH1F("jetPull", "", 50,0,1);
  
  TH1F * h_pT = new TH1F("pT", "", 100, 0, 500);
  TH1F * h_pTcut = new TH1F("pTcut", "", 100, 0, 500);
  TH1F * h_nref = new TH1F("nref", "", 12, 0, 12);
  TH1F * h_jetCount = new TH1F("jetCount", "", 12, 0, 12);
  TH1F * h_eta = new TH1F("eta", "", 60, -2, 2);
  TH1F * h_phi = new TH1F("phi", "", 60, -3.15, 3.15);
  TH1F * h_weight = new TH1F("weighted", "", 1200, 0, 1200);
  TH1F * h_unweight = new TH1F("unweighted", "", 1200, 0, 1200);
  TH1F * h_etaParticle = new TH1F("etaParticle", "", 60, -3, 3);
  TH1F * h_yParticle = new TH1F("yParticle", "", 60, -3, 3);
  TH2F * h_jetPull_2d = new TH2F("jetPull_2d","Jet pT vs Eta for pT>20",30,-3,+3,60,-3.5,3.5);
  TH2D * h_jetPull_3d = new TH2D("jetPull_3d","phi vs y vs pT" ,30,-3,+3,60,-3.5,3.5);
  //TH1F * pTCutoffs = new TH1F("pTcutoffs", "Jet Pull Angle -- with pThat Cutoffs", 60, 0.01, 500);
  //TH1F * pullAngle = new TH1F("pullAngle", "Jet Pull Angle vs pT", 60, -3.5, 3.5); 

  //DEFINE TREE VARIABLES +++++++++++++++++++++++++++++++++++
  Float_t pt[1000];  
  Float_t eta[1000];  
  Float_t phi[1000];  
  Float_t y[1000];
  Float_t genpt[1000];
  Float_t genphi[1000];
  Float_t geny[1000];
  Float_t geneta[1000];

  TVector2 jetPull;
  
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
  Float_t dot = 0;
  Double_t mag = 0;
  Double_t thrust_temp = 0;
  Double_t thrust_max = 0;
  Double_t dot_maj = 0;
  Double_t dot_min = 0;
  Double_t min_temp = 0;
  Double_t maj_temp = 0;
  Double_t thrust_maj_max =0;
  Double_t thrust_min_max = 0;
  TVector3 max_thrust_axis;
  TVector3 p3Norm;
  TVector2 r; 
  Float_t max_eta = 0;   Float_t temp_eta = 0;  
  Float_t max_phi = 0;   Float_t temp_phi = 0;
  Int_t max_nref;
  Int_t jetCount = 0;
  Int_t eventCount = 0;
  Float_t jetPullMag  = 0;
  Float_t pullAngle = 0; 

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
    
    string s = "";
    string str = s.append(filename[ifile]);
    file = TFile::Open(str.c_str());
    string w = Form("weights/weights_pp_%d.root", ifile+1);
    cout << w << endl; 
    weight_file = TFile::Open(w.c_str());
    //if(weighted) weight_file = TFile::Open(w.c_str());

    if (debug) cout << "\n **** =========================== New File ================================= **** \n ";
    cout << "File Name: " << filename[ifile] << endl;
    cout << "File Number: " << ifile << "/" << startfile << " to " << endfile << endl;

    //define trees and file +++++++++++++++++
    t = (TTree*)file->Get(Form("ak%dPFJetAnalyzer/t", radius));
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
      
      jetCount = 0;
      if(TMath::Abs(vz) > 15 || halo == 0) {
	//isGoodEvt = 0;
	continue;
      }
      //isGoodEvt = 1;
      
      // get the pt vector for each event which passes you jet selections based on pT and eta
      vector <float> pt_v;
      vector <float> eta_v;
      vector <float> phi_v;

      for(int ij = 0; ij<nref; ++ij){
	if(pt[ij] > ptCutLead && fabs(eta[ij]) < etaCut){
	  isGoodEvt = 1; 
	}
	if(pt[ij] > ptCut && fabs(eta[ij]) < etaCut){
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
      if(!isGoodEvt){continue;}

      if(debug) cout<< " \n ******* New Event ******** " << endl;
      if(debug) cout<< " ******* " << nentry << " ******** " << endl;
      
      //reset maximum values
      eventCount++;
      evt++;
      thrust_max = 0;

      vector <double> px;
      vector <double> py;
      vector <double> pz;
      
      for(Long64_t naxis = 0; naxis < NJets_Sel; ++naxis){
	float axis_jet_pt = pt_v[naxis];
	float axis_jet_eta = eta_v[naxis];
	float axis_jet_phi = phi_v[naxis];
      	px.push_back((double)axis_jet_pt * TMath::Cos(axis_jet_phi));
	py.push_back((double)axis_jet_pt * TMath::Sin(axis_jet_phi));
	pz.push_back((double)axis_jet_pt * TMath::SinH(axis_jet_eta));
      }

      //EVENT LOOP  ====================================================
      for(Long64_t njet = 0; njet < NJets_Sel; ++njet){	
	
	r.Set(0.0,0.0);
	jetPull.Set(0.0,0.0);
	if(debug) cout << "NEW JET: "<<  ngen << " *******************************" << endl;

	//PARTICLE LOOP STARTS HERE ===================================
	for(Long64_t npart = 0; npart < ngen; ++npart){

	  //Delta R^2 cut
	  /*
	  if((Double_t)(TMath::Power(genphi[npart],2) + TMath::Power(geneta[npart],2)) < 0.2){
	    continue;
	  }
	  */
	  if(debug) cout << "NEW PART: " << npart << " ================" << endl;
	  if(debug) cout << "y: " << geny[npart] << endl;
	  if(debug) cout << "azimuthal: " << genphi[npart] << endl;
	  
	  TVector2 partAxis = TVector2(geny[npart], genphi[npart]);
	  TVector2 jetAxis = TVector2(y[njet], phi[njet]);
	  TVector2 newJetPull; 
	  r.Set((partAxis.X() - jetAxis.X()),(partAxis.Y() - jetAxis.Y()));
	  if(debug) r.Print();
	  if(debug) cout << genpt[ngen] << endl;

	  jetPullMag = (r.X()*r.X() + r.Y()*r.Y())* genpt[ngen]/pt[njet];
	  newJetPull = jetPullMag*r;
	  jetPull.Set(jetPull.X() + newJetPull.X(), jetPull.Y() + newJetPull.Y());
	  if(debug) cout << "JET PULL MAG: " << jetPullMag << endl;

	  h_yParticle->Fill(geny[npart], pThat_weight);
	  h_etaParticle->Fill(geneta[npart], pThat_weight);
	  
	}//PARTICLE LOOP ENDS HERE ===========

	if(jetPullMag != 0){
	  h_jetPull->Fill(jetPullMag, pThat_weight);
	  h_jetPull_2d->Fill(jetPull.X(), jetPull.Y());
	  Double_t weightNum = pt[njet]*pThat_weight;
	  h_jetPull_3d->Fill(jetPull.X(), jetPull.Y(), pt[njet]);
	}
      }//end jet loop =======================
      
      //fill all the maximum values before finishing
      h_thrust->Fill(thrust_max);
      h_min->Fill(thrust_min_max);
      h_maj->Fill(thrust_maj_max);
      
      //h_jetPull->Fill(jetPull, pThat_weight);
      h_thrust->Fill(thrust_max, pThat_weight);
      h_min->Fill(thrust_min_max, pThat_weight);
      h_maj->Fill(thrust_maj_max, pThat_weight);
      
      h_nref->Fill(nref);
      h_jetCount->Fill(NJets_Sel);
      h_weight->Fill(pThat,pThat_weight);
      h_unweight->Fill(pThat);     
      
      pt_v.clear();
      eta_v.clear();
      phi_v.clear();
      
      px.clear();
      py.clear();
      pz.clear();
      
      //thrust_tree->Fill(); 
      
    }//end of event loop
    
    gROOT->GetListOfFiles()->Remove(file);
    gROOT->GetListOfFiles()->Remove(weight);
    
    cout << "Events Selected: " << eventCount << endl;
    cout << "File Finished" << endl; 
    
  }//end file loop
  
  //SCALE HISTOGRAMS +++++++++++++++++++++++++++++++++++
  //h_jetPull->Scale(1./h_jetPull->Integral());
  h_thrust->Scale(1./h_thrust->Integral());
  h_maj->Scale(1./h_maj->Integral());
  h_min->Scale(1./h_min->Integral());
  
  //CREATE PLOTS VS. dN/dT +++++++++++++++++++++++
  TH1F * h_jetPullScaled = DivideByBinWidth(h_jetPull, "jetPull_scaled");
  TH1F * h_T = DivideByBinWidth(h_thrust, "thrust_scaled");
  TH1F * h_Tmaj = DivideByBinWidth(h_maj, "thrust_maj_scaled");
  TH1F * h_Tmin = DivideByBinWidth(h_min, "thrust_min_scaled");
  
  save_File->cd();
  /* 
  h_jetPullScaled->Print("base");
  h_jetPull->Print("base"); 
  h_T->Print("base");
  h_Tmaj->Print("base");
  h_Tmin->Print("base");
  h_pT->Print("base");
  h_pTcut->Print("base");
  h_nref->Print("base");
  h_jetCount->Print("base");
  h_eta->Print("base");
  h_phi->Print("base");
  h_weight->Print("base");
  h_yParticle->Print("base");
  h_etaParticle->Print("base");
  */
  h_jetPullScaled->Write();
  h_jetPull->Write();
  h_T->Write();
  h_Tmaj->Write();
  h_Tmin->Write();
  h_thrust->Write();
  h_maj->Write();
  h_min->Write();
  h_pT->Write();
  h_pTcut->Write();
  h_nref->Write();
  h_jetCount->Write();
  h_eta->Write();
  h_phi->Write();
  h_weight->Write();
  h_unweight->Write();
  h_yParticle->Write();
  h_etaParticle->Write();
  h_jetPull_2d->Write();
  h_jetPull_3d->Write();
  
  save_File->Write();
  save_File->Close();
  gROOT->GetListOfFiles()->Remove(save_File);

  timer.Stop();
  cout<<"Macro finished: "<<endl;
  
  cout<<"CPU time (min)  = "<<(Float_t)timer.CpuTime()/60<<endl;
  cout<<"Real time (min) = "<<(Float_t)timer.RealTime()/60<<endl;
    
}//end of plot thrust
