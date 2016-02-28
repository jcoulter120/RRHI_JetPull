// Jennifer Coulter
// June 9th 2015
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
  Int_t jobNum = 60;
  
  //define trees and file
  TFile * fin = TFile::Open(Form("pythia_jetPull_%d.root",jobNum));

  TH1F * h_jetPull = (TH1F*)fin->Get("jetPull");
  TH2F * h_jetPull_2d = (TH2F*)fin->Get("jetPull_2d");
  TH1F * h_unweighted = (TH1F*)fin->Get("unweighted");
  TH1F * h_weighted = (TH1F*)fin->Get("weighted");
  TH3F * h_jetPull_3d = (TH3F*)fin->Get("jetPull_3d");
  TH1F * h_eta = (TH1F*)fin->Get("etaParticle");
  TH1F * h_y = (TH1F*)fin->Get("yParticle");
    
  TCanvas * canvas = new TCanvas("c","jetPull Test", 1100, 1100);
  canvas->Divide(2,3);

  if (debug) {
    cout << "JetPull Entries: " << endl;
    h_jetPull->Print("base");
    h_jetPull_3d->Print("base");
  }
  
  //plot with scaling, error bars, and marker styles
  canvas->cd(1)->SetLogy();
  h_jetPull->SetTitle("JetPull Magnitude vs. Counts"); 
  h_jetPull->SetXTitle("JetPull");
  h_jetPull->SetYTitle("Counts");
  h_jetPull->SetAxisRange(-5.5,5.5,"X");
  h_jetPull->GetXaxis()->CenterTitle();
  h_jetPull->GetYaxis()->CenterTitle();
  h_jetPull->Draw();

  canvas->cd(2);
  h_jetPull_2d->SetTitle("JetPull"); 
  h_jetPull_2d->SetXTitle("y");
  h_jetPull_2d->SetYTitle("phi");
  h_jetPull_2d->SetAxisRange(-5.5,5.5,"X");
  h_jetPull_2d->GetXaxis()->CenterTitle();
  h_jetPull_2d->GetYaxis()->CenterTitle();
  h_jetPull_2d->Draw("colz");

  canvas->cd(3)->SetLogy();
  h_unweighted->SetXTitle("Unweighted PThat");
  h_unweighted->SetYTitle("Counts");
  h_unweighted->SetAxisRange(1,1200,"X");
  h_unweighted->GetXaxis()->CenterTitle();
  h_unweighted->GetYaxis()->CenterTitle();
  h_unweighted->Draw();
  
  canvas->cd(4)->SetLogy();
  h_weighted->SetXTitle("Weighted PThat");
  h_weighted->SetYTitle("Counts");
  h_weighted->SetAxisRange(1,1200,"X");
  h_weighted->GetXaxis()->CenterTitle();
  h_weighted->GetYaxis()->CenterTitle();
  h_weighted->Draw();

  canvas->cd(5);
  h_jetPull_3d->SetTitle("JetPull"); 
  h_jetPull_3d->SetXTitle("y");
  h_jetPull_3d->SetYTitle("phi");
  h_jetPull_3d->SetZTitle("pT");
  h_jetPull_3d->SetAxisRange(-100.5,100.5,"X");
  h_jetPull_3d->SetAxisRange(-100.5,100.5,"Y");
  h_jetPull_3d->SetAxisRange(-100.5,100.5,"Z");
  h_jetPull_3d->GetXaxis()->CenterTitle();
  h_jetPull_3d->GetYaxis()->CenterTitle();
  h_jetPull_3d->GetZaxis()->CenterTitle();
  h_jetPull_3d->Draw("LEGO1");
  

  canvas->cd(6)->SetLogy();
  h_eta->SetXTitle("Pseduorapidity, Rapidity");
  h_eta->SetYTitle("Counts");
  h_eta->SetAxisRange(1,1200,"X");
  h_eta->GetXaxis()->CenterTitle();
  h_eta->GetYaxis()->CenterTitle();
  h_eta->SetMarkerColor(3);
  TLegend * a = new TLegend(0.2,.70,.4,.85);
  a->AddEntry(h_eta, "pseudorapidity", "p");
  
  a->AddEntry(h_y, "rapidity", "p");
  h_eta->Draw();
  h_y->Draw("same");
  a->Draw("same");
  
  canvas->Print(Form("test_pp_jetPull_%d.pdf", jobNum)); 

}//end of plot thrust
