// Jennifer Coulter
// June 9th 201
// Rutgers University, jennifer.coulter@cern.ch
//
// Plotting macro for sample thrust variables. 
//
//

#include <iostream>
#include <stdio.h>
#include <fstream>
#include <sstream>
#include <string>
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

using namespace std;

//plot jetPull
void jetPull_draw(){

  gStyle->SetOptStat(0);

  bool debug = true; 
  Int_t jobNum = 3;
  
  //define trees and file
  TFile * fin = TFile::Open(Form("pythia_jetPull_%d.root",jobNum));

  TH1F * h_jetPull = (TH1F*)fin->Get("jetPull");
  TH2F * h_jetPull_2d = (TH2F*)fin->Get("jetPull_2d");
  TH1F * h_unweighted = (TH1F*)fin->Get("unweighted");
  TH1F * h_weighted = (TH1F*)fin->Get("weighted");
  TH2F * h_jetPull_3d = (TH2F*)fin->Get("jetPull_3d");
  TH1F * h_eta = (TH1F*)fin->Get("etaParticle");
  TH1F * h_y = (TH1F*)fin->Get("yParticle");
  TH1F * h_pullAngle = (TH1F*)fin->Get("pullAngle");
  TH1F * h_30 = (TH1F*)fin->Get("pTcutoff30");  
  TH1F * h_60 = (TH1F*)fin->Get("pTcutoff60");
  TH1F * h_90 = (TH1F*)fin->Get("pTcutoff90");  
  TH1F * h_120 = (TH1F*)fin->Get("pTcutoff120");
  TH1F * h_leadingpT = (TH1F*)fin->Get("leading_pT");
  TH1F * h_leadingEta = (TH1F*)fin->Get("leadingJeteta");
  TH1F * h_leadingPhi = (TH1F*)fin->Get("leadingJetphi");
    
  TCanvas * canvas = new TCanvas("c","jetPull Test", 1100, 1100);
  canvas->Divide(2,3);

  if (debug) {
    cout << "JetPull Entries: " << endl;
    h_jetPull->Print("base");
    h_jetPull_3d->Print("base");
    h_90->Print("base");
  }
  
  //plot with scaling, error bars, and marker styles
  canvas->cd(1)->SetLogy();
  /*
  h_jetPull->SetTitle("JetPull  Magnitude vs. Counts"); 
  h_jetPull->SetXTitle("JetPull");
  h_jetPull->SetYTitle("Counts");
  h_jetPull->SetAxisRange(-5.5,5.5,"X");
  h_jetPull->GetXaxis()->CenterTitle();
  h_jetPull->GetYaxis()->CenterTitle();
  h_jetPull->Draw();
  */
  
  h_pullAngle->SetTitle("Pull Angle vs. Counts"); 
  h_pullAngle->SetXTitle("Pull Angle");
  h_pullAngle->SetYTitle("Counts");
  h_pullAngle->SetAxisRange(-0.1,5.5,"X");
  h_pullAngle->GetXaxis()->CenterTitle();
  h_pullAngle->GetYaxis()->CenterTitle();
  h_pullAngle->Draw();

  canvas->cd(2);
  h_jetPull_2d->SetTitle("JetPull"); 
  h_jetPull_2d->SetXTitle("y");
  h_jetPull_2d->SetYTitle("phi");
  h_jetPull_2d->SetAxisRange(-5.5,5.5,"X");
  h_jetPull_2d->GetXaxis()->CenterTitle();
  h_jetPull_2d->GetYaxis()->CenterTitle();
  h_jetPull_2d->Draw("colz");

  canvas->cd(3)->SetLogy();
  h_leadingpT->SetTitle("Leading Jet pT");
  h_leadingpT->SetXTitle("pT");
  h_leadingpT->SetYTitle("Counts");
  h_leadingpT->SetAxisRange(1,500,"X");
  h_leadingpT->GetXaxis()->CenterTitle();
  h_leadingpT->GetYaxis()->CenterTitle();
  h_leadingpT->Draw();
  /*
  h_unweighted->SetTitle("Unweighted PThat");
  h_unweighted->SetXTitle("PThat");
  h_unweighted->SetYTitle("Counts");
  h_unweighted->SetAxisRange(1,1200,"X");
  h_unweighted->GetXaxis()->CenterTitle();
  h_unweighted->GetYaxis()->CenterTitle();
  h_unweighted->Draw();
  */
  
  canvas->cd(4)->SetLogy();
  h_leadingEta->SetTitle("Eta Difference");
  h_leadingEta->SetXTitle("Delta Eta");
  h_leadingEta->SetYTitle("Counts");
  h_leadingEta->GetXaxis()->CenterTitle();
  h_leadingEta->GetYaxis()->CenterTitle();
  h_leadingEta->Draw();
  /*
  h_weighted->SetXTitle("PThat");
  h_weighted->SetTitle("Weighted PThat");
  h_weighted->SetYTitle("Counts");
  h_weighted->SetAxisRange(1,1200,"X");
  h_weighted->GetXaxis()->CenterTitle();
  h_weighted->GetYaxis()->CenterTitle();
  h_weighted->Draw();
  */

  canvas->cd(5)->SetLogy();
  h_leadingPhi->SetTitle("Leading Jet Phi");
  h_leadingPhi->SetXTitle("Phi");
  h_leadingPhi->SetYTitle("Counts");
  h_leadingPhi->GetXaxis()->CenterTitle();
  h_leadingPhi->GetYaxis()->CenterTitle();
  h_leadingPhi->Draw();
  /*
  h_jetPull_3d->SetTitle("JetPull vs. pT"); 
  h_jetPull_3d->SetXTitle("y");
  h_jetPull_3d->SetYTitle("phi");
  h_jetPull_3d->SetZTitle("pT");
  h_jetPull_3d->SetAxisRange(-3.5,3.5,"X");
  h_jetPull_3d->SetAxisRange(-3.5,3.5,"Y");
  h_jetPull_3d->SetAxisRange(0.00000001,100000000,"Z");
  h_jetPull_3d->GetXaxis()->CenterTitle();
  h_jetPull_3d->GetYaxis()->CenterTitle();
  h_jetPull_3d->GetZaxis()->CenterTitle();
  h_jetPull_3d->Draw("LEGO");

  h_leadingEta->SetTitle("Leading Jet Eta");
  h_leadingEta->SetXTitle("Eta");
  h_leadingEta->SetYTitle("Counts");
  h_leadingEta->GetXaxis()->CenterTitle();
  h_leadingEta->GetYaxis()->CenterTitle();
  h_leadingEta->Draw();
  */

  canvas->cd(6)->SetLogy();
  h_30->SetTitle("PT Chunks");
  h_30->SetXTitle("jetPull magnitude");
  h_30->SetYTitle("Counts");
  h_30->SetAxisRange(0.00000000000001,1,"Y");
  h_30->SetAxisRange(-5.5,5.5,"X");
  h_30->GetXaxis()->CenterTitle();
  h_30->GetYaxis()->CenterTitle();
  h_30->SetMarkerColor(6);
  h_30->SetMarkerStyle(22);
  h_60->SetMarkerStyle(23); 
  TLegend * a = new TLegend(0.7,.70,.9,.85);
  a->AddEntry(h_30, "30", "p");
  a->AddEntry(h_60, "60", "p");
  a->AddEntry(h_90, "90", "p");
  a->AddEntry(h_120, "120", "p");
  h_30->Draw();
  h_60->Draw("same");
  h_90->Draw("same");
  h_120->Draw("same");
  a->Draw("same");

    /*
  h_jetPull_3d->SetTitle("JetPull vs pT"); 
  h_jetPull_3d->Draw("psrlego2");
  h_jetPull_3d->SetAxisRange(0.1,10000,"Z");

  h_eta->SetTitle("Pseduorapidity, Rapidity Comparison");
  h_eta->SetXTitle("Eta, y");
  h_eta->SetYTitle("Counts");
  h_eta->SetAxisRange(-3.5,3.5,"X");
  h_eta->GetXaxis()->CenterTitle();
  h_eta->GetYaxis()->CenterTitle();
  h_eta->SetMarkerColor(6);
  h_eta->SetMarkerStyle(22);
  h_y->SetMarkerStyle(23); 
  TLegend * a = new TLegend(0.7,.70,.9,.85);
  a->AddEntry(h_eta, "pseudorapidity", "p");
  a->AddEntry(h_y, "rapidity", "p");
  h_eta->Draw();
  h_y->Draw("same");
  a->Draw("same");
    */
			   
  canvas->Print(Form("test_pp_jetPull_%d.pdf", jobNum)); 

}//end of plot thrust
