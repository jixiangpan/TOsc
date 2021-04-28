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

void func_get_contours(TH2D *h2_CL, TGraph *gh_CL[3], int index)
{
  TString roostr = "";
  
  const int Ncontour = 3;
  double contours[Ncontour] = {0};
  contours[0] = 0.9;
  contours[1] = 0.95;
  contours[2] = 0.9973;
  
  ///////
  roostr = TString::Format("canv_h2_CL_%d", index);
  TCanvas *canv_h2_CL = new TCanvas(roostr, roostr, 800, 600);
  h2_CL->SetStats(0);
  h2_CL->SetContour(Ncontour, contours);
  h2_CL->Draw("cont z list");
  canv_h2_CL->Update(); // Needed to force the plotting and retrieve the contours in TGraphs
     
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
  for(int ic=0; ic<Ncontour; ic++) {
    int np_curv = gh_curv_mc[ic]->GetN();
    for(int idx=0; idx<np_curv; idx++) {
      double t14 = 0;
      double m41 = 0;
      gh_curv_mc[ic]->GetPoint(idx, t14, m41);
      gh_CL[ic]->SetPoint(idx, pow(10., t14), pow(10., m41) );      
    }    
  }

  delete canv_h2_CL;
}

//////////////////////////////////////////////////////////////////
///////////////////////// MAIN ///////////////////////////////////
//////////////////////////////////////////////////////////////////

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
       
  ///////
  const int Ncontour = 3;
  //int colors[Ncontour] = {kBlue, kGreen+3, kRed};
  int colors[Ncontour] = {kBlue, kBlack, kRed};
  
  ////////////////////////////////////////////////////////////////////////////////////////
  
  int bins_theta = 10;
  int bins_dm2   = 10;

  
  ////// X: sin22t14, 1e-2 -> 1   ---> "log10()" ---> -2 -> 0
  ////// Y: m41^2,    1e-1 -> 20  ---> "log10()" ---> -1 -> 1.30103

  
  TH2D *h2_space_data = new TH2D("h2_space_data", "h2_space_data", bins_theta, -2, 0, bins_dm2, -1, 1.30103);
  
  TH2D *h2_space_sens = new TH2D("h2_space_sens", "h2_space_sens", bins_theta, -2, 0, bins_dm2, -1, 1.30103);

  roostr = "h2_space_sens_xx_plus_1sigma";
  TH2D *h2_space_sens_xx_plus_1sigma = new TH2D(roostr, roostr, bins_theta, -2, 0, bins_dm2, -1, 1.30103);

  roostr = "h2_space_sens_xx_minus_1sigma";
  TH2D *h2_space_sens_xx_minus_1sigma = new TH2D(roostr, roostr, bins_theta, -2, 0, bins_dm2, -1, 1.30103);

  
  //roostr = "out_sum.txt";
  //roostr = "out_sum_40by40.txt";
  roostr = "out_result_0001_0001.txt";
  
  ifstream InputFile(roostr, ios::in);
  if(!InputFile) { cerr<<" No input-list"<<endl; exit(1); }

  for(int idx=1; idx<=bins_theta*bins_dm2; idx++) {
    int xbin(0), ybin(0);
    double d4v(0), d3v(0), dd(0), CL(0);
    InputFile>>xbin>>ybin>>d4v>>d3v>>dd>>CL;

    ///
    h2_space_data->SetBinContent(xbin, ybin, CL);

    ///
    double CL_sens = 1 + TMath::Erf( (d4v-d3v)/sqrt(8*fabs(d4v)) );
    h2_space_sens->SetBinContent(xbin, ybin, 1 - CL_sens);
    
    ///
    double CL_xx_plus_1sigma = 1 + TMath::Erf( (d4v+2*sqrt(fabs(d4v)) -d3v-2*sqrt(fabs(d3v)))/sqrt(8*fabs(d4v+2*sqrt(fabs(d4v)))) );
    h2_space_sens_xx_plus_1sigma->SetBinContent(xbin, ybin, 1 - CL_xx_plus_1sigma);
 
    ///
    double CL_xx_minus_1sigma = 1 + TMath::Erf( (d4v-2*sqrt(fabs(d4v)) -d3v+2*sqrt(fabs(d3v)))/sqrt(8*fabs(d4v-2*sqrt(fabs(d4v)))) );
    h2_space_sens_xx_minus_1sigma->SetBinContent(xbin, ybin, 1 - CL_xx_minus_1sigma);
 
  }

  int index = 0;
  
  //////
  index++;
  TGraph *gh_CLs_data[Ncontour];
  for(int idx=0; idx<Ncontour; idx++) {
    roostr = TString::Format("gh_CLs_data_%d", idx);
    gh_CLs_data[idx] = new TGraph();
    gh_CLs_data[idx]->SetName(roostr);
    gh_CLs_data[idx]->SetLineColor( colors[idx] );
  }
  func_get_contours( h2_space_data, gh_CLs_data, index);

  //////
  index++;
  TGraph *gh_CLs_sens[Ncontour];
  for(int idx=0; idx<Ncontour; idx++) {
    roostr = TString::Format("gh_CLs_sens_%d", idx);
    gh_CLs_sens[idx] = new TGraph();
    gh_CLs_sens[idx]->SetName(roostr);
    gh_CLs_sens[idx]->SetLineColor( colors[idx] );
    gh_CLs_sens[idx]->SetLineStyle( 7 );
  }
  func_get_contours( h2_space_sens, gh_CLs_sens, index);
  
  //////
  index++;
  TGraph *gh_CLs_sens_xx_plus_1sigma[Ncontour];
  for(int idx=0; idx<Ncontour; idx++) {
    roostr = TString::Format("gh_CLs_sens_xx_plus_1sigma_%d", idx);
    gh_CLs_sens_xx_plus_1sigma[idx] = new TGraph();
    gh_CLs_sens_xx_plus_1sigma[idx]->SetName(roostr);
    gh_CLs_sens_xx_plus_1sigma[idx]->SetLineColor( colors[idx] );
    gh_CLs_sens_xx_plus_1sigma[idx]->SetLineStyle( 7 );
    gh_CLs_sens_xx_plus_1sigma[idx]->SetLineWidth( 4 );
  }
  func_get_contours( h2_space_sens_xx_plus_1sigma, gh_CLs_sens_xx_plus_1sigma, index);

  //////
  index++;
  TGraph *gh_CLs_sens_xx_minus_1sigma[Ncontour];
  for(int idx=0; idx<Ncontour; idx++) {
    roostr = TString::Format("gh_CLs_sens_xx_minus_1sigma_%d", idx);
    gh_CLs_sens_xx_minus_1sigma[idx] = new TGraph();
    gh_CLs_sens_xx_minus_1sigma[idx]->SetName(roostr);
    gh_CLs_sens_xx_minus_1sigma[idx]->SetLineColor( colors[idx] );
    gh_CLs_sens_xx_minus_1sigma[idx]->SetLineStyle( 7 );
    gh_CLs_sens_xx_minus_1sigma[idx]->SetLineWidth( 4 );
  }
  func_get_contours( h2_space_sens_xx_minus_1sigma, gh_CLs_sens_xx_minus_1sigma, index);

  ////////////////////////////////////////////////////////////
    
  double xxlow = h2_space_data->GetXaxis()->GetBinLowEdge(1);
  double xxhgh = h2_space_data->GetXaxis()->GetBinUpEdge(bins_theta);
    
  double yylow = h2_space_data->GetYaxis()->GetBinLowEdge(1);
  double yyhgh = h2_space_data->GetYaxis()->GetBinUpEdge(bins_dm2);

  roostr = "h2_basic_CLs_data";
  TH2D *h2_basic_CLs_data = new TH2D(roostr, "",
				     bins_theta, pow(10, xxlow), pow(10, xxhgh),
				     bins_dm2, pow(10, yylow), pow(10, yyhgh));

  roostr = "canv_h2_basic_CLs_data";
  TCanvas *canv_h2_basic_CLs_data = new TCanvas(roostr, roostr, 800, 600);
  canv_h2_basic_CLs_data->SetLeftMargin(0.15);
  canv_h2_basic_CLs_data->SetRightMargin(0.1);
  canv_h2_basic_CLs_data->SetBottomMargin(0.18);
  canv_h2_basic_CLs_data->SetLogx();
  canv_h2_basic_CLs_data->SetLogy();
  h2_basic_CLs_data->Draw();
  h2_basic_CLs_data->SetXTitle("sin^{2}2#theta_{14}");
  h2_basic_CLs_data->SetYTitle("#Deltam^{2}_{41} [eV^{2}]");
  h2_basic_CLs_data->GetXaxis()->CenterTitle(1);
  h2_basic_CLs_data->GetYaxis()->CenterTitle(1);
  h2_basic_CLs_data->GetXaxis()->SetTitleOffset(1.25);
  h2_basic_CLs_data->GetYaxis()->SetTitleOffset(1.1);

  h2_basic_CLs_data->GetYaxis()->SetLabelSize(0.045);
  h2_basic_CLs_data->GetYaxis()->SetTitleSize(0.045);    
  h2_basic_CLs_data->GetXaxis()->SetLabelSize(0.045);
  h2_basic_CLs_data->GetXaxis()->SetTitleSize(0.045);
  h2_basic_CLs_data->GetYaxis()->SetTitleOffset(1.4);
  h2_basic_CLs_data->GetXaxis()->SetTitleOffset(1.4);
   
  //gh_CLs_sens_xx_plus_1sigma[1]->Draw("same l");
  //gh_CLs_sens_xx_minus_1sigma[1]->Draw("same l");
  TGraph *gh_CLs_sens_95_pm_1sigma = new TGraph();
  for(int idx=0; idx<gh_CLs_sens_xx_plus_1sigma[1]->GetN(); idx++) {
    double xx(0), yy(0);
    gh_CLs_sens_xx_plus_1sigma[1]->GetPoint(idx,xx,yy);
    gh_CLs_sens_95_pm_1sigma->SetPoint( gh_CLs_sens_95_pm_1sigma->GetN(), xx, yy );
  }
  for(int idx=gh_CLs_sens_xx_minus_1sigma[1]->GetN()-1; idx>=0; idx--) {
    double xx(0), yy(0);
    gh_CLs_sens_xx_minus_1sigma[1]->GetPoint(idx,xx,yy);
    gh_CLs_sens_95_pm_1sigma->SetPoint( gh_CLs_sens_95_pm_1sigma->GetN(), xx, yy );
  }
  gh_CLs_sens_95_pm_1sigma->Draw("same f");
  gh_CLs_sens_95_pm_1sigma->SetFillColor(kGreen);
   
  for(int ic=0; ic<3; ic++) {
    if( ic==0 ) continue;
    if( ic==2 ) continue;
    gh_CLs_data[ic]->Draw("same l");
    gh_CLs_sens[ic]->Draw("same l");
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

  TLegend *lg_data = new TLegend(0.18-0.02,0.22,0.48-0.02,0.45);
  lg_data->SetFillStyle(0);
  lg_data->SetBorderSize(0);
  //lg_data->AddEntry(gh_CLs_data[0], "CL_{s} exclusion, 90% CL", "l");
  lg_data->AddEntry(gh_CLs_data[1], "CL_{s} exclusion, 95% CL", "l");
  //lg_data->AddEntry(gh_CLs_data[2], "CL_{s} exclusion, 99.73% CL", "l");
  lg_data->AddEntry( gh_n4, "Neutrino-4 best-fit", "p" );
  lg_data->AddEntry( gh_RAA, "SBL+Gallium Anomaly (RAA)", "p" );
  lg_data->Draw();
  lg_data->SetTextSize(0.045);
  lg_data->SetTextFont(132);
  lg_data->Draw();

  canv_h2_basic_CLs_data->SaveAs("canv_h2_basic_CLs_data.png");

  
}
