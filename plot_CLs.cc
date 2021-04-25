#include<iostream>
#include<fstream>
#include<sstream>
#include<cmath>
#include "stdlib.h"
using namespace std;

#include<map>

#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TLine.h"
#include "TMath.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "THStack.h"

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TBranch.h"

#include "TRandom3.h"
#include "TGaxis.h"
#include "TStyle.h"
#include "TPaveText.h"
#include "TText.h"
#include "TLatex.h"

#include "TCanvas.h"
#include "TVirtualPad.h"
#include "TPad.h"
#include "TLegend.h"
#include "TString.h"
#include "TColor.h"

#include "TPrincipal.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"

void plot_CLs()
{
  
  //////////////////////////////////////////////////////////////////////////////////////// Draw style
  
  gStyle->SetOptStat(0);
  //gStyle->SetPalette(kBird);

  double snWidth = 2;

  // use medium bold lines and thick markers
  gStyle->SetLineWidth(snWidth);
  gStyle->SetFrameLineWidth(snWidth);
  gStyle->SetHistLineWidth(snWidth);
  gStyle->SetFuncWidth(snWidth);
  gStyle->SetGridWidth(snWidth);
  gStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes
  gStyle->SetMarkerStyle(20);
  gStyle->SetMarkerSize(1.2);
  gStyle->SetEndErrorSize(4);
  gStyle->SetEndErrorSize(0);

  ////////////////////////////////////////////////////////////////////////////////////////

  TString roostr = "";
  
  int bins_theta = 100;
  int bins_dm2   = 100;

  ////// X: sin22t14, 1e-2 -> 1   ---> "log10()" ---> -2 -> 0
  ////// Y: m41^2,    1e-1 -> 20  ---> "log10()" ---> -1 -> 1.30103      
  TH2D *h2_space = new TH2D("h2_space", "h2_space", bins_theta, -2, 0, bins_dm2, -1, 1.30103);

  roostr = "out_sum.txt";
  ifstream InputFile(roostr, ios::in);
  if(!InputFile) { cerr<<" No input-list"<<endl; exit(1); }

  for(int idx=1; idx<=bins_theta*bins_dm2; idx++) {
    int xbin(0), ybin(0);
    double d4v(0), d3v(0), dd(0), CL(0);
    InputFile>>xbin>>ybin>>d4v>>d3v>>dd>>CL;
    h2_space->SetBinContent(xbin, ybin, CL);
  }

     
    ///////
    const int Ncontour = 3;
    double contours[Ncontour] = {0};
    contours[0] = 0.9;
    contours[1] = 0.95;
    contours[2] = 0.9973;
    
    int colors[Ncontour] = {kBlue, kGreen+3, kRed};

    TH2D *h2_counter_mc = (TH2D*)h2_space->Clone("h2_counter_mc");
    
    ///////
    TCanvas *canv_h2_counter_mc = new TCanvas("canv_h2_counter_mc", "canv_h2_counter_mc", 800, 600);
    h2_counter_mc->SetStats(0);
    h2_counter_mc->SetContour(Ncontour, contours);
    h2_counter_mc->Draw("cont z list");
    canv_h2_counter_mc->Update(); // Needed to force the plotting and retrieve the contours in TGraphs
     
    // Get Contours
    TObjArray *conts_mc = (TObjArray*)gROOT->GetListOfSpecials()->FindObject("contours");
    TList* contLevel_mc = NULL;
    TGraph* gh_curv_mc[10] = {0};

    Int_t nGraphs_mc    = 0;
    Int_t TotalConts_mc = 0;

    if (conts_mc == NULL){
      printf("*** No Contours Were Extracted!\n");
      TotalConts_mc = 0;
      exit(1);
    } else {
      TotalConts_mc = conts_mc->GetSize();
    }

    printf("TotalConts_mc = %d\n", TotalConts_mc);

    for(int i = 0; i < TotalConts_mc; i++){
      contLevel_mc = (TList*)conts_mc->At(i);
      printf("Contour %d has %d Graphs\n", i, contLevel_mc->GetSize());
      nGraphs_mc += contLevel_mc->GetSize();
    }

    nGraphs_mc = 0;
    for(int i = 0; i < TotalConts_mc; i++){
      contLevel_mc = (TList*)conts_mc->At(i);

      // Get first graph from list on curves on this level
      gh_curv_mc[i] = (TGraph*)contLevel_mc->First();     
    }

    ///////
    TGraph *gh_CLs_mc[Ncontour];
    for(int ic=0; ic<Ncontour; ic++) {
      gh_CLs_mc[ic] = new TGraph();
      roostr = TString::Format("gh_CLs_mc_%d", ic);
      gh_CLs_mc[ic]->SetName(roostr);
    }
  
    for(int ic=0; ic<Ncontour; ic++) {
      int np_curv = gh_curv_mc[ic]->GetN();
      for(int idx=0; idx<np_curv; idx++) {
	double t14 = 0;
	double m41 = 0;
	gh_curv_mc[ic]->GetPoint(idx, t14, m41);
	gh_CLs_mc[ic]->SetPoint(idx, pow(10., t14), pow(10., m41) );      
      }    
    }

    ////////////////////////////////////////////////////////////
    
    double xxlow = h2_space->GetXaxis()->GetBinLowEdge(1);
    double xxhgh = h2_space->GetXaxis()->GetBinUpEdge(bins_theta);
    
    double yylow = h2_space->GetYaxis()->GetBinLowEdge(1);
    double yyhgh = h2_space->GetYaxis()->GetBinUpEdge(bins_dm2);

    roostr = "h2_basic_CLs_mc";
    TH2D *h2_basic_CLs_mc = new TH2D(roostr, "", bins_theta, pow(10, xxlow), pow(10, xxhgh), bins_dm2, pow(10, yylow), pow(10, yyhgh));

    roostr = "canv_h2_basic_CLs_mc";
    TCanvas *canv_h2_basic_CLs_mc = new TCanvas(roostr, roostr, 800, 600);
    canv_h2_basic_CLs_mc->SetLeftMargin(0.15);
    canv_h2_basic_CLs_mc->SetRightMargin(0.1);
    canv_h2_basic_CLs_mc->SetBottomMargin(0.18);
    canv_h2_basic_CLs_mc->SetLogx();
    canv_h2_basic_CLs_mc->SetLogy();
    h2_basic_CLs_mc->Draw();
    h2_basic_CLs_mc->SetXTitle("sin^{2}2#theta_{14}");
    h2_basic_CLs_mc->SetYTitle("#Deltam^{2}_{41} [eV^{2}]");
    h2_basic_CLs_mc->GetXaxis()->CenterTitle(1);
    h2_basic_CLs_mc->GetYaxis()->CenterTitle(1);
    h2_basic_CLs_mc->GetXaxis()->SetTitleOffset(1.25);
    h2_basic_CLs_mc->GetYaxis()->SetTitleOffset(1.1);

    h2_basic_CLs_mc->GetYaxis()->SetLabelSize(0.045);
    h2_basic_CLs_mc->GetYaxis()->SetTitleSize(0.045);    
    h2_basic_CLs_mc->GetXaxis()->SetLabelSize(0.045);
    h2_basic_CLs_mc->GetXaxis()->SetTitleSize(0.045);
    h2_basic_CLs_mc->GetYaxis()->SetTitleOffset(1.4);
    h2_basic_CLs_mc->GetXaxis()->SetTitleOffset(1.4);
    
    for(int ic=0; ic<3; ic++) {
      gh_CLs_mc[ic]->Draw("same l");
      gh_CLs_mc[ic]->SetLineColor( colors[ic] );
    }

    TGraph *gh_n4 = new TGraph();
    gh_n4->SetPoint(0, 0.26, 7.25);
    gh_n4->SetMarkerStyle(22);
    gh_n4->SetMarkerSize(1.8);
    gh_n4->Draw("same p");

    TGraph *gh_RAA = new TGraph();
    gh_RAA->SetPoint(0, 0.165, 2.39);
    gh_RAA->SetMarkerStyle(20);
    gh_RAA->SetMarkerSize(1.8);
    gh_RAA->Draw("same p");

    TLegend *lg_mc = new TLegend(0.18-0.02,0.22,0.48-0.02,0.45);
    lg_mc->SetFillStyle(0);
    lg_mc->SetBorderSize(0);
    lg_mc->AddEntry(gh_CLs_mc[0], "CL_{s} exclusion, 90% CL", "l");
    lg_mc->AddEntry(gh_CLs_mc[1], "CL_{s} exclusion, 95% CL", "l");
    lg_mc->AddEntry(gh_CLs_mc[2], "CL_{s} exclusion, 99.73% CL", "l");
    lg_mc->AddEntry( gh_n4, "Neutrino-4 best-fit", "p" );
    lg_mc->AddEntry( gh_RAA, "SBL+Gallium Anomaly (RAA)", "p" );
    lg_mc->Draw();
    lg_mc->SetTextSize(0.045);
    lg_mc->SetTextFont(132);
    lg_mc->Draw();

    canv_h2_basic_CLs_mc->SaveAs("canv_h2_basic_CLs_mc.png");
  
}
