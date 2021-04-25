#include<iostream>
#include<fstream>
#include<sstream>
#include<cmath>
#include "stdlib.h"
using namespace std;

#include<map>

#include "WCPLEEANA/TOsc.h"

#include "WCPLEEANA/Configure_Osc.h"

#include "TApplication.h"

/////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////// MAIN //////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////

/*  
  usage:
  
  make clean
  make  
  ./read_TLee_v20 -f 1 -p 1
  
  ---> README:
  ---> Makefile: comment the line "ROOTSYS=/home/xji/data0/software/root_build", if you have your own "ROOTSYS"
  ---> minuit2 is in the ROOT
*/

int main(int argc, char** argv)
{
  TString roostr = "";
  
  double scaleF_POT = 1;
  int ifile = 1;

  int space_xbin = 1;
  int space_ybin = 1;
  
  for(int i=1; i<argc; i++) {
    
    if( strcmp(argv[i],"-p")==0 ) {
      stringstream convert( argv[i+1] );
      if(  !( convert>>scaleF_POT ) ) { cerr<<" ---> Error scaleF_POT !"<<endl; exit(1); }
    }

    if( strcmp(argv[i],"-f")==0 ) {
      stringstream convert( argv[i+1] );
      if(  !( convert>>ifile ) ) { cerr<<" ---> Error ifile !"<<endl; exit(1); }
    }

    if( strcmp(argv[i],"-x")==0 ) {
      stringstream convert( argv[i+1] );
      if(  !( convert>>space_xbin ) ) { cerr<<" ---> Error space_xbin !"<<endl; exit(1); }
    }

    if( strcmp(argv[i],"-y")==0 ) {
      stringstream convert( argv[i+1] );
      if(  !( convert>>space_ybin ) ) { cerr<<" ---> Error space_ybin !"<<endl; exit(1); }
    }
    
  }

  cout<<endl<<TString::Format(" ---> check, scaleF_POT %6.4f, ifile %d, space xbin/ybin %d %d",
			      scaleF_POT, ifile, space_xbin, space_ybin
			      )<<endl<<endl;
  
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

  if( !config_Osc::flag_display_graphics ) {
    gROOT->SetBatch( 1 );
  }

  ////////////////////////////////////////////////////////////////////////////////////////
  
  TApplication theApp("theApp",&argc,argv);

  ////////////////////////////////////////////////////////////////////////////////////////

  TOsc *Osc_test = new TOsc();

  ////////// just do it one time in the whole procedure

  Osc_test->channels_observation = config_Osc::channels_observation;

  Osc_test->Set_Spectra_MatrixCov(config_Osc::eventlist_dir,
				  config_Osc::event_summation_afterscale_file,
				  config_Osc::centralvalue_noosc_file,
				  config_Osc::syst_result_dir, "d", "e");
  
  Osc_test->scaleF_POT = scaleF_POT;  
  Osc_test->Apply_POT_scaled();

  ////////// can do any times
  
  Osc_test->flag_syst_flux       = config_Osc::flag_syst_flux;
  Osc_test->flag_syst_geant      = config_Osc::flag_syst_geant;
  Osc_test->flag_syst_Xs         = config_Osc::flag_syst_Xs;
  Osc_test->flag_syst_detector   = config_Osc::flag_syst_detector;
  Osc_test->flag_syst_additional = config_Osc::flag_syst_additional;
  Osc_test->flag_syst_MCstat     = config_Osc::flag_syst_MCstat;
  
  //////////

  // cout<<" test AA"<<endl;
  // Osc_test->Apply_Oscillation();

  // cout<<" test AB"<<endl;
  // Osc_test->Set_OscPars(0.5, 5);
  // Osc_test->Set_OscPars(0.26, 7.25);
  // Osc_test->Apply_Oscillation();

  // Osc_test->Set_OscPars(0.1, 2);
  // Osc_test->Set_Collapse();
  // Osc_test->Set_Asimov2dataFIT();
  // Osc_test->Set_data2dataFIT();  
  // Osc_test->Minimization_OscPars_FullCov(0.1, 1.8, 1);
  
  ////////////////////////////////////////////////////////
  
  if( 0 ) {
    double test_s22theta = 0.26;
    double test_dm2      = 7.25;
    
    //////    
    double dchi2_4vAsimov = 0;
    double dchi2_3vAsimov = 0;
    
    double chi2_4v_on_data = 0;
    double chi2_3v_on_data = 0;
    double dchi2_data      = 0;
    
    /////////////////////////////////// 4v Asimov  
    Osc_test->Set_OscPars(test_s22theta, test_dm2);
    Osc_test->Set_Collapse();
    Osc_test->Set_Asimov2dataFIT();
    //Osc_test->Minimization_OscPars_FullCov(test_s22theta, test_dm2, 1); cout<<"  ---> min chi2: " <<Osc_test->minimization_chi2<<endl;
    Osc_test->Minimization_OscPars_FullCov(0, 7.25, 1);
    dchi2_4vAsimov = (-Osc_test->minimization_chi2);
    
    /////////////////////////////////// 3v Asimov    
    Osc_test->Set_OscPars(0, 1);
    Osc_test->Set_Collapse();
    Osc_test->Set_Asimov2dataFIT();
    //Osc_test->Minimization_OscPars_FullCov(0, 1, 1); cout<<"  ---> min chi2: " <<Osc_test->minimization_chi2<<endl;
    Osc_test->Minimization_OscPars_FullCov(test_s22theta, test_dm2, 1);
    dchi2_3vAsimov = Osc_test->minimization_chi2;

    /////////////////////////////////// data  
    Osc_test->Set_data2dataFIT();
    
    Osc_test->Minimization_OscPars_FullCov(0, 1, 1);
    chi2_3v_on_data = Osc_test->minimization_chi2;
  
    Osc_test->Minimization_OscPars_FullCov(test_s22theta, test_dm2, 1);
    chi2_4v_on_data = Osc_test->minimization_chi2;
    
    dchi2_data = chi2_4v_on_data - chi2_3v_on_data;

    ///////    
    double delta_4v = dchi2_4vAsimov;
    double delta_3v = dchi2_3vAsimov;
    double delta_dd = dchi2_data;
    
    double numerator_erf = TMath::Erf( (delta_4v-delta_dd)/sqrt(8.*fabs(delta_4v)) );
    double denominator_erf = TMath::Erf( (delta_3v-delta_dd)/sqrt(8.*fabs(delta_3v)) );
    
    double CLs = (1+numerator_erf)/(1+denominator_erf);
    double result = 1 - CLs;  
    cout<<endl<<" ---> test result: "<<result<<endl<<endl;
  }

  ////////////////////////////////////////////////////////

  if( 1 ) {

    int bins_theta = 10;
    int bins_dm2   = 10;

    ////// X: sin22t14, 1e-2 -> 1   ---> "log10()" ---> -2 -> 0
    ////// Y: m41^2,    1e-1 -> 20  ---> "log10()" ---> -1 -> 1.30103      
    TH2D *h2_space = new TH2D("h2_space", "h2_space", bins_theta, -2, 0, bins_dm2, -1, 1.30103);
        
    double test_s22theta = pow( 10, h2_space->GetXaxis()->GetBinCenter(space_xbin) );
    double test_dm2      = pow( 10, h2_space->GetYaxis()->GetBinCenter(space_ybin) );
    
    //////    
    double dchi2_4vAsimov = 0;
    double dchi2_3vAsimov = 0;
    
    double chi2_4v_on_data = 0;
    double chi2_3v_on_data = 0;
    double dchi2_data      = 0;
    
    /////////////////////////////////// 4v Asimov  
    Osc_test->Set_OscPars(test_s22theta, test_dm2);
    Osc_test->Set_Collapse();
    Osc_test->Set_Asimov2dataFIT();
    //Osc_test->Minimization_OscPars_FullCov(test_s22theta, test_dm2, 1); cout<<"  ---> min chi2: " <<Osc_test->minimization_chi2<<endl;
    Osc_test->Minimization_OscPars_FullCov(0, 7.25, 1);
    dchi2_4vAsimov = (-Osc_test->minimization_chi2);
    
    /////////////////////////////////// 3v Asimov    
    Osc_test->Set_OscPars(0, 1);
    Osc_test->Set_Collapse();
    Osc_test->Set_Asimov2dataFIT();
    //Osc_test->Minimization_OscPars_FullCov(0, 1, 1); cout<<"  ---> min chi2: " <<Osc_test->minimization_chi2<<endl;
    Osc_test->Minimization_OscPars_FullCov(test_s22theta, test_dm2, 1);
    dchi2_3vAsimov = Osc_test->minimization_chi2;

    /////////////////////////////////// data  
    Osc_test->Set_data2dataFIT();
    
    Osc_test->Minimization_OscPars_FullCov(0, 1, 1);
    chi2_3v_on_data = Osc_test->minimization_chi2;
  
    Osc_test->Minimization_OscPars_FullCov(test_s22theta, test_dm2, 1);
    chi2_4v_on_data = Osc_test->minimization_chi2;
    
    dchi2_data = chi2_4v_on_data - chi2_3v_on_data;

    ///////    
    double delta_4v = dchi2_4vAsimov;
    double delta_3v = dchi2_3vAsimov;
    double delta_dd = dchi2_data;
    
    double numerator_erf = TMath::Erf( (delta_4v-delta_dd)/sqrt(8.*fabs(delta_4v)) );
    double denominator_erf = TMath::Erf( (delta_3v-delta_dd)/sqrt(8.*fabs(delta_3v)) );
    
    double CLs = (1+numerator_erf)/(1+denominator_erf);
    double CL = 1 - CLs;
    //cout<<" ---> "<<CL<<endl;

    roostr = TString::Format("out_result_%04d_%04d.txt", space_xbin, space_ybin);
    ofstream ListWrite(roostr, ios::out|ios::trunc);
    ListWrite<<TString::Format("%4d %4d %12.6f %12.6f %12.6f %18.15f",
			       space_xbin, space_ybin,
			       delta_4v, delta_3v, delta_dd, CL
			       )<<endl;
    ListWrite.close();
  }
  
  ////////////////////////////////////////////////////////
  
  if( 0 ) {

    int bins_theta = 10;
    int bins_dm2   = 10;

    ////// X: sin22t14, 1e-2 -> 1   ---> "log10()" ---> -2 -> 0
    ////// Y: m41^2,    1e-1 -> 20  ---> "log10()" ---> -1 -> 1.30103      
    TH2D *h2_space = new TH2D("h2_space", "h2_space", bins_theta, -2, 0, bins_dm2, -1, 1.30103);
    
    for(int ibin=1; ibin<=bins_theta; ibin++) {
      
      cout<<TString::Format(" ---> processing %4d/%4d", ibin, bins_theta)<<endl;
      
      for(int jbin=1; jbin<=bins_dm2; jbin++) {

	cout<<TString::Format(" ---> processing sub %4d - %4d", ibin, jbin)<<endl;
	
	double xcenter = h2_space->GetXaxis()->GetBinCenter(ibin);
	double ycenter = h2_space->GetYaxis()->GetBinCenter(jbin);

	double test_s22theta = pow( 10, xcenter );
	double test_dm2      = pow( 10, ycenter );

	//////
  
	double dchi2_4vAsimov = 0;
	double dchi2_3vAsimov = 0;

	double chi2_4v_on_data = 0;
	double chi2_3v_on_data = 0;
	double dchi2_data      = 0;

	/////////////////////////////////// 4v Asimov
  
	Osc_test->Set_OscPars(test_s22theta, test_dm2);
	Osc_test->Set_Collapse();
	Osc_test->Set_Asimov2dataFIT();
	//Osc_test->Minimization_OscPars_FullCov(test_s22theta, test_dm2, 1); cout<<"  ---> min chi2: " <<Osc_test->minimization_chi2<<endl;
	Osc_test->Minimization_OscPars_FullCov(0, 7.25, 1);
	dchi2_4vAsimov = (-Osc_test->minimization_chi2);
	//cout<<" chi2_4vAsimov "<<dchi2_4vAsimov<<endl;

	/////////////////////////////////// 3v Asimov
  
	Osc_test->Set_OscPars(0, 1);
	Osc_test->Set_Collapse();
	Osc_test->Set_Asimov2dataFIT();
	//Osc_test->Minimization_OscPars_FullCov(0, 1, 1); cout<<"  ---> min chi2: " <<Osc_test->minimization_chi2<<endl;
	Osc_test->Minimization_OscPars_FullCov(test_s22theta, test_dm2, 1);
	dchi2_3vAsimov = Osc_test->minimization_chi2;
	//cout<<" chi2_3vAsimov "<<dchi2_3vAsimov<<endl;

	/////////////////////////////////// data
  
	Osc_test->Set_data2dataFIT();
  
	Osc_test->Minimization_OscPars_FullCov(0, 1, 1);
	chi2_3v_on_data = Osc_test->minimization_chi2;
	//cout<<" chi2_3v_on_data "<<chi2_3v_on_data<<endl;
  
	Osc_test->Minimization_OscPars_FullCov(test_s22theta, test_dm2, 1);
	chi2_4v_on_data = Osc_test->minimization_chi2;
	//cout<<" chi2_4v_on_data "<<chi2_4v_on_data<<endl;
  
	dchi2_data = chi2_4v_on_data - chi2_3v_on_data;

	///////
	///////

	double delta_4v = dchi2_4vAsimov;
	double delta_3v = dchi2_3vAsimov;
	double delta_dd = dchi2_data;
  
	double numerator_erf = TMath::Erf( (delta_4v-delta_dd)/sqrt(8.*fabs(delta_4v)) );
	double denominator_erf = TMath::Erf( (delta_3v-delta_dd)/sqrt(8.*fabs(delta_3v)) );
  
	double CLs = (1+numerator_erf)/(1+denominator_erf);
	double result = 1 - CLs;  
	//cout<<" ---> result: "<<result<<endl;

	h2_space->SetBinContent(ibin, jbin, result);
	
      }// jbin
    }// ibin
     
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
      return 1;
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

    int bins_xx = 100;
    int bins_yy = 100;

    double xxlow = h2_space->GetXaxis()->GetBinLowEdge(1);
    double xxhgh = h2_space->GetXaxis()->GetBinUpEdge(bins_theta);
    
    double yylow = h2_space->GetYaxis()->GetBinLowEdge(1);
    double yyhgh = h2_space->GetYaxis()->GetBinUpEdge(bins_dm2);

    roostr = "h2_basic_CLs_mc";
    TH2D *h2_basic_CLs_mc = new TH2D(roostr, "", bins_xx, pow(10, xxlow), pow(10, xxhgh), bins_yy, pow(10, yylow), pow(10, yyhgh));

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

  ///////////////////////////////////////////////////////////////
  
  if( config_Osc::flag_display_graphics ) {
    cout<<endl<<" Enter Ctrl+c to end the program"<<endl;
    cout<<" Enter Ctrl+c to end the program"<<endl<<endl;  
    theApp.Run();
  }
    
  return 0;
}
