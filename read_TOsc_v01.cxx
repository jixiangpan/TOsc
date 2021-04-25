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
  
  for(int i=1; i<argc; i++) {    
    if( strcmp(argv[i],"-p")==0 ) {
      stringstream convert( argv[i+1] );
      if(  !( convert>>scaleF_POT ) ) { cerr<<" ---> Error scaleF_POT !"<<endl; exit(1); }
    }

    if( strcmp(argv[i],"-f")==0 ) {
      stringstream convert( argv[i+1] );
      if(  !( convert>>ifile ) ) { cerr<<" ---> Error ifile !"<<endl; exit(1); }
    }
  }
  
  cout<<endl<<" ---> check, scaleF_POT "<<scaleF_POT<<", ifile "<<ifile<<endl<<endl;

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

  // if( !config_Lee::flag_display_graphics ) {
  //   gROOT->SetBatch( 1 );
  // }

  ////////////////////////////////////////////////////////////////////////////////////////
  
  TApplication theApp("theApp",&argc,argv);

  ////////////////////////////////////////////////////////////////////////////////////////

  TOsc *Osc_test = new TOsc();

  ////////// just do it one time in the whole procedure

  Osc_test->scaleF_POT = scaleF_POT;
  
  Osc_test->channels_observation = config_Osc::channels_observation;;
  
  Osc_test->Set_Spectra_MatrixCov(config_Osc::eventlist_dir,
				  config_Osc::event_summation_afterscale_file,
				  config_Osc::centralvalue_noosc_file,
				  config_Osc::syst_result_dir, "d", "e");
  
  Osc_test->Apply_POT_scaled();

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
  
  //////////

  int bins_theta = 10;
  int bins_dm2   = 10;

  ////// X: sin22t14, 1e-2 -> 1   ---> "log10()" ---> -2 -> 0
  ////// Y: m41^2,    1e-1 -> 20  ---> "log10()" ---> -1 -> 1.30103  

  TH2D *h2_space = new TH2D("h2_space", "h2_space", bins_theta, -2, 0, bins_dm2, -1, 1.30103);

  
  ////////////////////// 

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
  cout<<endl<<" ---> test result: "<<result<<endl<<endl;


  if( 1 ) {

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

  
  
  cout<<endl<<" Enter Ctrl+c to end the program"<<endl;
  cout<<" Enter Ctrl+c to end the program"<<endl<<endl;
  
  theApp.Run();
  
  
  ////////////////////////////////////////////////////////////////////////////////////////
  /*
  TLee *Lee_test = new TLee();
    
  ////////// just do it one time in the whole procedure

  Lee_test->channels_observation   = config_Lee::channels_observation;
  Lee_test->syst_cov_flux_Xs_begin = config_Lee::syst_cov_flux_Xs_begin;
  Lee_test->syst_cov_flux_Xs_end   = config_Lee::syst_cov_flux_Xs_end;
  Lee_test->syst_cov_mc_stat_begin = config_Lee::syst_cov_mc_stat_begin;
  Lee_test->syst_cov_mc_stat_end   = config_Lee::syst_cov_mc_stat_end;  
  
  Lee_test->scaleF_POT = scaleF_POT;
  Lee_test->Set_config_file_directory(config_Lee::spectra_file, config_Lee::flux_Xs_directory,
                                      config_Lee::detector_directory, config_Lee::mc_directory);
  Lee_test->Set_Spectra_MatrixCov();
  Lee_test->Set_POT_implement();
  Lee_test->Set_TransformMatrix();

  ////////// can do any times
  
  Lee_test->flag_syst_flux_Xs    = config_Lee::flag_syst_flux_Xs;
  Lee_test->flag_syst_detector   = config_Lee::flag_syst_detector;
  Lee_test->flag_syst_additional = config_Lee::flag_syst_additional;
  Lee_test->flag_syst_mc_stat    = config_Lee::flag_syst_mc_stat;
  
  Lee_test->scaleF_Lee = config_Lee::Lee_strength_for_outputfile_covariance_matrix;
  Lee_test->scaleF_Lee = config_Lee::Lee_strength_for_GoF;
  Lee_test->Set_Collapse();

  Lee_test->flag_Lee_minimization_after_constraint = config_Lee::flag_Lee_minimization_after_constraint;
  
  //////////
  
  TFile *file_collapsed_covariance_matrix = new TFile("file_collapsed_covariance_matrix.root", "recreate");
  
  TTree *tree_config = new TTree("tree", "configure information");
  int flag_syst_flux_Xs = config_Lee::flag_syst_flux_Xs;
  int flag_syst_detector = config_Lee::flag_syst_detector;
  int flag_syst_additional = config_Lee::flag_syst_additional;
  int flag_syst_mc_stat = config_Lee::flag_syst_mc_stat;
  double user_Lee_strength_for_output_covariance_matrix = config_Lee::Lee_strength_for_outputfile_covariance_matrix;
  double user_scaleF_POT = scaleF_POT;
  vector<double>vc_val_GOF;
  vector<int>vc_val_GOF_NDF;
  tree_config->Branch("flag_syst_flux_Xs", &flag_syst_flux_Xs, "flag_syst_flux_Xs/I" );
  tree_config->Branch("flag_syst_detector", &flag_syst_detector, "flag_syst_detector/I" );
  tree_config->Branch("flag_syst_additional", &flag_syst_additional, "flag_syst_additional/I" );
  tree_config->Branch("flag_syst_mc_stat", &flag_syst_mc_stat, "flag_syst_mc_stat/I" );
  tree_config->Branch("user_Lee_strength_for_output_covariance_matrix", &user_Lee_strength_for_output_covariance_matrix,
                      "user_Lee_strength_for_output_covariance_matrix/D" );
  tree_config->Branch("user_scaleF_POT", &user_scaleF_POT, "user_scaleF_POT/D" );
  tree_config->Branch("vc_val_GOF", &vc_val_GOF);
  tree_config->Branch("vc_val_GOF_NDF", &vc_val_GOF_NDF);
  file_collapsed_covariance_matrix->cd();

  Lee_test->matrix_absolute_cov_newworld.Write("matrix_absolute_cov_newworld");// (bins, bins)
  Lee_test->matrix_absolute_flux_cov_newworld.Write("matrix_absolute_flux_cov_newworld");
  Lee_test->matrix_absolute_Xs_cov_newworld.Write("matrix_absolute_Xs_cov_newworld");
  Lee_test->matrix_absolute_detector_cov_newworld.Write("matrix_absolute_detector_cov_newworld");
  Lee_test->matrix_absolute_mc_stat_cov_newworld.Write("matrix_absolute_mc_stat_cov_newworld");
  Lee_test->matrix_absolute_additional_cov_newworld.Write("matrix_absolute_additional_cov_newworld");
                 
  for(auto it=Lee_test->matrix_input_cov_detector_sub.begin(); it!=Lee_test->matrix_input_cov_detector_sub.end(); it++) {
    int idx = it->first;
    roostr = TString::Format("matrix_absolute_detector_sub_cov_newworld_%02d", idx);
    Lee_test->matrix_absolute_detector_sub_cov_newworld[idx].Write(roostr);
  }
     
  Lee_test->matrix_pred_newworld.Write("matrix_pred_newworld");// (1, bins)
  Lee_test->matrix_data_newworld.Write("matrix_data_newworld");// (1, bins)  
  //file_collapsed_covariance_matrix->Close();
  
  //////////

  if( config_Lee::flag_plotting_systematics ) Lee_test->Plotting_systematics();
  
  //////////////////////////////////////////////////////////////////////////////////////// Goodness of fit
  
  //Lee_test->scaleF_Lee = config_Lee::Lee_strength_for_GoF;
  //Lee_test->Set_Collapse();

  if( config_Lee::flag_GoF_output2file_default_0 ) {
    file_collapsed_covariance_matrix->cd();

    for(auto it=Lee_test->map_data_spectrum_ch_bin.begin(); it!=Lee_test->map_data_spectrum_ch_bin.end(); it++) {
      int val_ch = it->first;
      int size_map = it->second.size();
      int size_before = 0;
      for(int idx=1; idx<val_ch; idx++) {
	int size_current = Lee_test->map_data_spectrum_ch_bin[idx].size();
	size_before += size_current;
      }
      
      vector<int>vc_target_chs;
      for(int ibin=1; ibin<size_map; ibin++) {
	vc_target_chs.push_back( size_before + ibin -1 );
      }
      
      vector<int>vc_support_chs;

      Lee_test->Exe_Goodness_of_fit_detailed( vc_target_chs, vc_support_chs, 100 + val_ch );

      vc_val_GOF.push_back( Lee_test->val_GOF_noConstrain );
      vc_val_GOF_NDF.push_back( Lee_test->val_GOF_NDF );
      //cout<<" ---> "<<Lee_test->val_GOF_NDF<<"\t"<<Lee_test->val_GOF_noConstrain<<endl;
    }
    
    tree_config->Fill();
    tree_config->Write();
    file_collapsed_covariance_matrix->Close();
  }
  
  bool flag_both_numuCC            = config_Lee::flag_both_numuCC;// 1
  bool flag_CCpi0_FC_by_numuCC     = config_Lee::flag_CCpi0_FC_by_numuCC;// 2
  bool flag_CCpi0_PC_by_numuCC     = config_Lee::flag_CCpi0_PC_by_numuCC;// 3
  bool flag_NCpi0_by_numuCC        = config_Lee::flag_NCpi0_by_numuCC;// 4
  bool flag_nueCC_PC_by_numuCC_pi0 = config_Lee::flag_nueCC_PC_by_numuCC_pi0;// 5
  bool flag_nueCC_HghE_FC_by_numuCC_pi0_nueFC = config_Lee::flag_nueCC_HghE_FC_by_numuCC_pi0_nueFC;// 6, HghE>800 MeV
  bool flag_nueCC_LowE_FC_by_all   = config_Lee::flag_nueCC_LowE_FC_by_all;// 7
  bool flag_nueCC_FC_by_all        = config_Lee::flag_nueCC_FC_by_all;// 8
  
  ///////////////////////// gof
  
  if( flag_nueCC_FC_by_all ) {
    vector<int>vc_target_chs;
    vc_target_chs.push_back( 1 );
    
    vector<int>vc_support_chs;
    vc_support_chs.push_back( 2 );
    vc_support_chs.push_back( 3 );
    vc_support_chs.push_back( 4 );
    vc_support_chs.push_back( 5 );
    vc_support_chs.push_back( 6 );
    vc_support_chs.push_back( 7 );

    Lee_test->Exe_Goodness_of_fit( vc_target_chs, vc_support_chs, 8 );
  }

  ///////////////////////// gof
  
  if( flag_nueCC_LowE_FC_by_all ) {
    TMatrixD matrix_gof_trans( Lee_test->bins_newworld, 26*4 + 11*3 );// oldworld, newworld
    for( int ibin=1; ibin<=26*4 + 11*3; ibin++) matrix_gof_trans(ibin-1, ibin-1) = 1;
    
    TMatrixD matrix_gof_trans_T( matrix_gof_trans.GetNcols(), matrix_gof_trans.GetNrows() );
    matrix_gof_trans_T.Transpose( matrix_gof_trans );

    TMatrixD matrix_gof_pred = Lee_test->matrix_pred_newworld * matrix_gof_trans;
    TMatrixD matrix_gof_data = Lee_test->matrix_data_newworld * matrix_gof_trans;
    TMatrixD matrix_gof_syst = matrix_gof_trans_T * (Lee_test->matrix_absolute_cov_newworld) * matrix_gof_trans;

    Lee_test->Exe_Goodness_of_fit( 8, matrix_gof_trans.GetNcols()-8, matrix_gof_pred, matrix_gof_data, matrix_gof_syst, 7);
  }

  
  if( 0 ) {
    TMatrixD matrix_gof_trans( Lee_test->bins_newworld, 26*6 + 11*3 );// oldworld, newworld
    for( int ibin=1; ibin<=26*6 + 11*3; ibin++) matrix_gof_trans(ibin-1, ibin-1) = 1;
    
    TMatrixD matrix_gof_trans_T( matrix_gof_trans.GetNcols(), matrix_gof_trans.GetNrows() );
    matrix_gof_trans_T.Transpose( matrix_gof_trans );

    TMatrixD matrix_gof_pred = Lee_test->matrix_pred_newworld * matrix_gof_trans;
    TMatrixD matrix_gof_data = Lee_test->matrix_data_newworld * matrix_gof_trans;
    TMatrixD matrix_gof_syst = matrix_gof_trans_T * (Lee_test->matrix_absolute_cov_newworld) * matrix_gof_trans;

    Lee_test->Exe_Goodness_of_fit( 8, matrix_gof_trans.GetNcols()-8, matrix_gof_pred, matrix_gof_data, matrix_gof_syst, 7);
  }

  
  if( 0 ) {// first 6 bins--> 1 bin, constrained by others
    int nbins_first = 6;
    TMatrixD matrix_gof_trans( Lee_test->bins_newworld, 1 + (26-nbins_first) + 26*3 + 11*3 );// oldworld, newworld
    for( int ibin=1; ibin<=nbins_first; ibin++) matrix_gof_trans(ibin-1, 0) = 1;
    for( int ibin=1; ibin<=26*4+11*3-nbins_first; ibin++) matrix_gof_trans(nbins_first+ibin-1, ibin) = 1;
    
    
    TMatrixD matrix_gof_trans_T( matrix_gof_trans.GetNcols(), matrix_gof_trans.GetNrows() );
    matrix_gof_trans_T.Transpose( matrix_gof_trans );

    TMatrixD matrix_gof_pred = Lee_test->matrix_pred_newworld * matrix_gof_trans;
    TMatrixD matrix_gof_data = Lee_test->matrix_data_newworld * matrix_gof_trans;
    TMatrixD matrix_gof_syst = matrix_gof_trans_T * (Lee_test->matrix_absolute_cov_newworld) * matrix_gof_trans;

    Lee_test->Exe_Goodness_of_fit( 1, matrix_gof_trans.GetNcols()-1, matrix_gof_pred, matrix_gof_data, matrix_gof_syst, 201);
  }  
  
  if( 0 ) {// first 6 bins, constrained by others
    int nbins_first = 6;
    TMatrixD matrix_gof_trans( Lee_test->bins_newworld, 26*4 + 11*3 );// oldworld, newworld
    for( int ibin=1; ibin<=26*4 + 11*3; ibin++) matrix_gof_trans(ibin-1, ibin-1) = 1;
    
    TMatrixD matrix_gof_trans_T( matrix_gof_trans.GetNcols(), matrix_gof_trans.GetNrows() );
    matrix_gof_trans_T.Transpose( matrix_gof_trans );

    TMatrixD matrix_gof_pred = Lee_test->matrix_pred_newworld * matrix_gof_trans;
    TMatrixD matrix_gof_data = Lee_test->matrix_data_newworld * matrix_gof_trans;
    TMatrixD matrix_gof_syst = matrix_gof_trans_T * (Lee_test->matrix_absolute_cov_newworld) * matrix_gof_trans;

    Lee_test->Exe_Goodness_of_fit( nbins_first, matrix_gof_trans.GetNcols()-nbins_first, matrix_gof_pred, matrix_gof_data, matrix_gof_syst, 202);
  }

  ///////////////////////// gof
  
  if( flag_nueCC_HghE_FC_by_numuCC_pi0_nueFC ) {
    TMatrixD matrix_gof_trans( Lee_test->bins_newworld, (26-8) + 26*3 + 11*3 );// oldworld, newworld
    for( int ibin=1; ibin<=(26-8) + 26*3 + 11*3; ibin++) matrix_gof_trans(8+ibin-1, ibin-1) = 1;
    
    TMatrixD matrix_gof_trans_T( matrix_gof_trans.GetNcols(), matrix_gof_trans.GetNrows() );
    matrix_gof_trans_T.Transpose( matrix_gof_trans );

    TMatrixD matrix_gof_pred = Lee_test->matrix_pred_newworld * matrix_gof_trans;
    TMatrixD matrix_gof_data = Lee_test->matrix_data_newworld * matrix_gof_trans;
    TMatrixD matrix_gof_syst = matrix_gof_trans_T * (Lee_test->matrix_absolute_cov_newworld) * matrix_gof_trans;

    Lee_test->Exe_Goodness_of_fit( (26-8), matrix_gof_trans.GetNcols()-(26-8), matrix_gof_pred, matrix_gof_data, matrix_gof_syst, 6);    
  }

  ///////////////////////// gof
  
  if( flag_nueCC_PC_by_numuCC_pi0) {    
    vector<int>vc_target_chs;
    vc_target_chs.push_back( 2 );
    
    vector<int>vc_support_chs;
    vc_support_chs.push_back( 3 );
    vc_support_chs.push_back( 4 );
    vc_support_chs.push_back( 5 );
    vc_support_chs.push_back( 6 );
    vc_support_chs.push_back( 7 );

    Lee_test->Exe_Goodness_of_fit( vc_target_chs, vc_support_chs, 5 );
  }
  
  ///////////////////////// gof
  
  if( flag_both_numuCC ) {
    vector<int>vc_target_chs;
    vc_target_chs.push_back( 3 );
    vc_target_chs.push_back( 4 );
    
    vector<int>vc_support_chs;
    
    Lee_test->Exe_Goodness_of_fit( vc_target_chs, vc_support_chs, 1 );
  }

  ///////////////////////// gof
  
  if( flag_CCpi0_FC_by_numuCC ) { 
    vector<int>vc_target_chs;
    vc_target_chs.push_back( 5 );
    
    vector<int>vc_support_chs;
    vc_support_chs.push_back( 3 );
    vc_support_chs.push_back( 4 );

    Lee_test->Exe_Goodness_of_fit( vc_target_chs, vc_support_chs, 2 );
  }
  
  ///////////////////////// gof
  
  if( flag_CCpi0_PC_by_numuCC ) {
    vector<int>vc_target_chs;
    vc_target_chs.push_back( 6 );
    
    vector<int>vc_support_chs;
    vc_support_chs.push_back( 3 );
    vc_support_chs.push_back( 4 );

    Lee_test->Exe_Goodness_of_fit( vc_target_chs, vc_support_chs, 3 );
  }
  
  ///////////////////////// gof
  
  if( flag_NCpi0_by_numuCC ) {
    vector<int>vc_target_chs;
    vc_target_chs.push_back( 7 );
    
    vector<int>vc_support_chs;
    vc_support_chs.push_back( 3 );
    vc_support_chs.push_back( 4 );

    Lee_test->Exe_Goodness_of_fit( vc_target_chs, vc_support_chs, 4 );
  }

  ////////////////////////////////////////////////////////////////////

  if( 0 ) {
    vector<int>vc_target_chs;
    vc_target_chs.push_back( 1 );
    vc_target_chs.push_back( 2 );
    
    vector<int>vc_support_chs;
    vc_support_chs.push_back( 3 );
    vc_support_chs.push_back( 4 );
    vc_support_chs.push_back( 5 );
    vc_support_chs.push_back( 6 );
    vc_support_chs.push_back( 7 );

    Lee_test->Exe_Goodness_of_fit( vc_target_chs, vc_support_chs, 101 );
  }

  if( 0 ) {
    vector<int>vc_target_chs;
    //vc_target_chs.push_back( 1 );
    vc_target_chs.push_back( 2 );
    
    vector<int>vc_support_chs;
    vc_support_chs.push_back( 3 );
    vc_support_chs.push_back( 4 );
    vc_support_chs.push_back( 5 );
    vc_support_chs.push_back( 6 );
    vc_support_chs.push_back( 7 );

    Lee_test->Exe_Goodness_of_fit( vc_target_chs, vc_support_chs, 102 );
  }

  
  if( 0 ) {
    vector<int>vc_target_chs;
    for(int idx=1; idx<=8; idx++) vc_target_chs.push_back( idx-1 );
    
    vector<int>vc_support_chs;
    for(int idx=9; idx<=Lee_test->bins_newworld; idx++) {
      if( idx>=26+1 && idx<=26+8 ) continue;
      vc_support_chs.push_back( idx-1 );
    }

    Lee_test->Exe_Goodness_of_fit_detailed( vc_target_chs, vc_support_chs, 101 );
  }

  

  
  ////////////////////////////////////////////////////////////////////

  bool flag_publicnote = 0;
  
  if( flag_publicnote ) {
    
    ///////////////////////// gof
    
    if( 1 ) {
      vector<int>vc_target_chs;
      vc_target_chs.push_back( 3 );
      vc_target_chs.push_back( 4 ); 
      vector<int>vc_support_chs;      
      Lee_test->Exe_Goodness_of_fit( vc_target_chs, vc_support_chs, 1 );
    }

    if( 1 ) {
      vector<int>vc_target_chs;
      vc_target_chs.push_back( 5 );
      vector<int>vc_support_chs;
      vc_support_chs.push_back( 3 );
      vc_support_chs.push_back( 4 );      
      Lee_test->Exe_Goodness_of_fit( vc_target_chs, vc_support_chs, 2 );
    }

    if( 1 ) {
      vector<int>vc_target_chs;
      vc_target_chs.push_back( 6 );
      vector<int>vc_support_chs;
      vc_support_chs.push_back( 3 );
      vc_support_chs.push_back( 4 );      
      Lee_test->Exe_Goodness_of_fit( vc_target_chs, vc_support_chs, 3 );
    }

    if( 1 ) {
      vector<int>vc_target_chs;
      vc_target_chs.push_back( 7 );
      vector<int>vc_support_chs;
      vc_support_chs.push_back( 3 );
      vc_support_chs.push_back( 4 );      
      Lee_test->Exe_Goodness_of_fit( vc_target_chs, vc_support_chs, 4 );
    }

    if( 0 ) {    
      vector<int>vc_target_chs;
      vc_target_chs.push_back( 2 );
      
      vector<int>vc_support_chs;
      vc_support_chs.push_back( 3 );
      vc_support_chs.push_back( 4 );
      vc_support_chs.push_back( 5 );
      vc_support_chs.push_back( 6 );
      vc_support_chs.push_back( 7 );
      
      Lee_test->Exe_Goodness_of_fit( vc_target_chs, vc_support_chs, 5 );
    }
    
    if( 0 ) {
      TMatrixD matrix_gof_trans( Lee_test->bins_newworld, (26-8) + 26*3 + 11*3 );// oldworld, newworld
      for( int ibin=1; ibin<=(26-8) + 26*3 + 11*3; ibin++) matrix_gof_trans(8+ibin-1, ibin-1) = 1;
    
      TMatrixD matrix_gof_trans_T( matrix_gof_trans.GetNcols(), matrix_gof_trans.GetNrows() );
      matrix_gof_trans_T.Transpose( matrix_gof_trans );

      TMatrixD matrix_gof_pred = Lee_test->matrix_pred_newworld * matrix_gof_trans;
      TMatrixD matrix_gof_data = Lee_test->matrix_data_newworld * matrix_gof_trans;
      TMatrixD matrix_gof_syst = matrix_gof_trans_T * (Lee_test->matrix_absolute_cov_newworld) * matrix_gof_trans;

      Lee_test->Exe_Goodness_of_fit( (26-8), matrix_gof_trans.GetNcols()-(26-8), matrix_gof_pred, matrix_gof_data, matrix_gof_syst, 6);    
    }

    if( 0 ) {
      TMatrixD matrix_gof_trans( Lee_test->bins_newworld, 26*4 + 11*3 );// oldworld, newworld
      for( int ibin=1; ibin<=26*4 + 11*3; ibin++) matrix_gof_trans(ibin-1, ibin-1) = 1;
    
      TMatrixD matrix_gof_trans_T( matrix_gof_trans.GetNcols(), matrix_gof_trans.GetNrows() );
      matrix_gof_trans_T.Transpose( matrix_gof_trans );

      TMatrixD matrix_gof_pred = Lee_test->matrix_pred_newworld * matrix_gof_trans;
      TMatrixD matrix_gof_data = Lee_test->matrix_data_newworld * matrix_gof_trans;
      TMatrixD matrix_gof_syst = matrix_gof_trans_T * (Lee_test->matrix_absolute_cov_newworld) * matrix_gof_trans;

      Lee_test->Exe_Goodness_of_fit( 8, matrix_gof_trans.GetNcols()-8, matrix_gof_pred, matrix_gof_data, matrix_gof_syst, 7);
    }

  }
  
  
  //////////////////////////////////////////////////////////////////////////////////////// LEE strength fitting

  if( config_Lee::flag_Lee_strength_data ) {

    Lee_test->Set_measured_data();// use the measured data as the input data for the fitting

    Lee_test->Minimization_Lee_strength_FullCov(2, 0);// (initial value, fix or not)

    // cout<<endl<<TString::Format(" ---> Best fit of Lee strength: chi2 %6.3f, %5.3f +/- %5.3f",
    //                             Lee_test->minimization_chi2,
    //                             Lee_test->minimization_Lee_strength_val,
    //                             Lee_test->minimization_Lee_strength_err
    //                             )<<endl<<endl;

    cout<<endl<<TString::Format(" ---> Best fit of Lee strength: %6.4f,  chi2 %6.3f",                          
                                Lee_test->minimization_Lee_strength_val,
				Lee_test->minimization_chi2
                                )<<endl<<endl;    

    cout<<endl<<TString::Format(" ---> Best fit of Lee strength: %6.4f +/- %6.4f,  chi2 %6.3f",                          
                                Lee_test->minimization_Lee_strength_val,
				Lee_test->minimization_Lee_strength_err,
				Lee_test->minimization_chi2
                                )<<endl<<endl;    

    /////////////////////////////////////////
    
    double gmin = Lee_test->minimization_chi2;
    TGraph *gh_scan = new TGraph();
    double slow = 0;
    double shgh = 3;
    int nscan = 100;
    double val_max_dchi2 = 0;
    double step = (shgh-slow)/nscan;
    for(int idx=1; idx<=nscan; idx++) {
      if( idx%(max(1, nscan/10))==0 ) cout<<Form(" ---> scan %4.2f, %3d", idx*1./nscan, idx)<<endl;
      double val_s = slow + (idx-1)*step;
      Lee_test->Minimization_Lee_strength_FullCov(val_s, 1);// (initial value, fix or not)
      double val_chi2 = Lee_test->minimization_chi2;
      gh_scan->SetPoint( gh_scan->GetN(), val_s, val_chi2 - gmin);
      if( val_max_dchi2<val_chi2 - gmin ) val_max_dchi2 = val_chi2 - gmin;
    }
    
    double val_dchi2at0 = gh_scan->Eval(0);
    double val_dchi2at1 = gh_scan->Eval(1);
    if( fabs(val_dchi2at0)<1e-6 ) val_dchi2at0 = 0;
    
    cout<<endl<<Form(" ---> dchi2 at LEE 0/1: %7.4f %7.4f", val_dchi2at0, val_dchi2at1 )<<endl<<endl;
    
    TCanvas *canv_gh_scan = new TCanvas("canv_gh_scan", "canv_gh_scan", 900, 650);
    canv_gh_scan->SetLeftMargin(0.15); canv_gh_scan->SetRightMargin(0.1);
    canv_gh_scan->SetTopMargin(0.1); canv_gh_scan->SetBottomMargin(0.15);    
    gh_scan->Draw("al");
    gh_scan->GetXaxis()->SetTitle("LEE strength"); gh_scan->GetYaxis()->SetTitle("#Delta#chi^{2}");    
    gh_scan->GetXaxis()->SetLabelSize(0.05); gh_scan->GetXaxis()->SetTitleSize(0.05);
    gh_scan->GetYaxis()->SetLabelSize(0.05); gh_scan->GetYaxis()->SetTitleSize(0.05);
    gh_scan->GetXaxis()->CenterTitle(); gh_scan->GetYaxis()->CenterTitle();
    gh_scan->GetXaxis()->SetTitleOffset(1.2);
    gh_scan->GetYaxis()->SetRangeUser(0, val_max_dchi2*1.1);
   
    TLine *lineA_dchi2at1 = new TLine(1, 0, 1, val_dchi2at1);    
    lineA_dchi2at1->Draw("same");
    lineA_dchi2at1->SetLineWidth(2);
    lineA_dchi2at1->SetLineColor(kBlue);
    lineA_dchi2at1->SetLineStyle(7);
    TLine *lineB_dchi2at1 = new TLine(0, val_dchi2at1, 1, val_dchi2at1);    
    lineB_dchi2at1->Draw("same");
    lineB_dchi2at1->SetLineWidth(2);
    lineB_dchi2at1->SetLineColor(kBlue);
    lineB_dchi2at1->SetLineStyle(7);
    auto *tt_text_data = new TLatex( 0.2, val_dchi2at1*1.1, Form("#Delta#chi^{2} = %4.3f", val_dchi2at1) );
    tt_text_data->SetTextAlign(11); tt_text_data->SetTextSize(0.05); tt_text_data->SetTextAngle(0);
    tt_text_data->SetTextFont(42);  tt_text_data->Draw(); tt_text_data->SetTextColor(kBlue);

    canv_gh_scan->SaveAs("canv_gh_scan.png");
    
  }
  
  ////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////// MicroBooNE suggested /////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////

  if( config_Lee::flag_chi2_data_H0 ) {

    Lee_test->Set_measured_data();// use the measured data as the input data for the fitting    
    
    Lee_test->scaleF_Lee = 0;
    Lee_test->Set_Collapse();
    
    Lee_test->Minimization_Lee_strength_FullCov(0, 1);
    
    double chi2 = Lee_test->minimization_chi2;
    int ndf = Lee_test->bins_newworld;
    double p_value = TMath::Prob( chi2, ndf );
    
    cout<<Form(" ---> flag_chi2_data_H0, chi2/ndf %8.2f %3d %8.4f, p-value %f", chi2,ndf, chi2/ndf, p_value)<<endl;
  }

  
  if( config_Lee::flag_dchi2_H0toH1 ) {
    
    Lee_test->Set_measured_data();// use the measured data as the input data for the fitting
    
    Lee_test->scaleF_Lee = 0;
    Lee_test->Set_Collapse();

    Lee_test->Minimization_Lee_strength_FullCov(0, 1);    
    double chi2_H0 = Lee_test->minimization_chi2;
    
    Lee_test->Minimization_Lee_strength_FullCov(1, 1);    
    double chi2_H1 = Lee_test->minimization_chi2;

    double dchi2 = chi2_H0 - chi2_H1;
    
    int ndf = 1;
    double p_value = TMath::Prob( dchi2, ndf );    
    cout<<Form(" ---> flag_dchi2_H0toH1, chi2/ndf %8.2f %3d %8.4f, p-value %f", dchi2, ndf, dchi2/ndf, p_value)<<endl;

    p_value = TMath::Prob( -dchi2, ndf );
    cout<<Form(" ---> flag_dchi2_H1toH0, chi2/ndf %8.2f %3d %8.4f, p-value %f", -dchi2, ndf, -dchi2/ndf, p_value)<<endl;
    
  }

  
  ////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////// Advanced Statistics Analysis /////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////
  //
  // (*) Hypothesis test always reject the null
  //
  // [*] the realization of the functionality is a combination of the tools:
  //
  //         Lee_test->Set_measured_data();
  //
  //         Lee_test->scaleF_Lee = #;
  //         Lee_test->Set_Collapse();
  //
  //         Lee_test->Set_toy_Asimov();
  //
  //         Lee_test->Set_Variations( # );
  //         Lee_test->Set_toy_Variation( # );
  //
  //         Lee_test->Minimization_Lee_strength_FullCov(#, #);
  //
  
  /////////////////////////////////////////////////////// example: do fitting on Asimov sample

  if( 0 ) {
    Lee_test->scaleF_Lee = 1;
    Lee_test->Set_Collapse();
  
    Lee_test->Set_toy_Asimov();// use the Asimov sample as the input data for the fitting

    Lee_test->Minimization_Lee_strength_FullCov(2, 0);// (initial value, fix or not)

    cout<<endl<<TString::Format(" ---> Best fit of Lee strength: chi2 %6.2f, %5.2f +/- %5.2f",
				Lee_test->minimization_chi2,
				Lee_test->minimization_Lee_strength_val,
				Lee_test->minimization_Lee_strength_err
				)<<endl<<endl;
  }

  ////////////////////////////////////////////////////// example: do fitting on variation sample

  if( 0 ) {
    Lee_test->scaleF_Lee = 1;
    Lee_test->Set_Collapse();

    Lee_test->Set_Variations( 10 );// generate 10 variation samples
    Lee_test->Set_toy_Variation( 4 );// use the 4th sample as the input data for the fitting
    
    Lee_test->Minimization_Lee_strength_FullCov(2, 0);// (initial value, fix or not)

    cout<<endl<<TString::Format(" ---> Best fit of Lee strength: chi2 %6.2f, %5.2f +/- %5.2f",
				Lee_test->minimization_chi2,
				Lee_test->minimization_Lee_strength_val,
				Lee_test->minimization_Lee_strength_err
				)<<endl<<endl;
  }
  
  //////////////////////////////////////////////////////////////////////////////////////// example: simple versus simple likelihood ratio test

  if( 0 ) {    
    Lee_test->Set_measured_data();// use the measured data as the input data for the fitting

    Lee_test->Minimization_Lee_strength_FullCov(1, 1);// (initial value, fix or not)
    double val_chi2_Lee = Lee_test->minimization_chi2;

    Lee_test->Minimization_Lee_strength_FullCov(0, 1);// (initial value, fix or not)
    double val_chi2_sm = Lee_test->minimization_chi2;

    double val_dchi2 = val_chi2_Lee - val_chi2_sm;
    
    cout<<endl<<TString::Format(" ---> dchi2 = Lee - sm: %7.4f, LEE %7.4f, sm %7.4f", val_dchi2, val_chi2_Lee, val_chi2_sm)<<endl<<endl;
  }

  ////////////////////////////////////////////////////////// sensitivity calcualtion by FC
  
  if( 0 ) {
    double chi2_null_null8sm_true8sm  = 0;
    double chi2_gmin_null8sm_true8sm  = 0;
    double chi2_null_null8Lee_true8Lee = 0;
    double chi2_gmin_null8Lee_true8Lee = 0;
    
    TFile *file_out = new TFile(TString::Format("file_out_%03d.root", ifile), "recreate");
    TTree *tree = new TTree("tree", "tree");
    tree->Branch("chi2_null_null8sm_true8sm", &chi2_null_null8sm_true8sm, "chi2_null_null8sm_true8sm/D" );
    tree->Branch("chi2_gmin_null8sm_true8sm", &chi2_gmin_null8sm_true8sm, "chi2_gmin_null8sm_true8sm/D" );
    tree->Branch("chi2_null_null8Lee_true8Lee", &chi2_null_null8Lee_true8Lee, "chi2_null_null8Lee_true8Lee/D" );
    tree->Branch("chi2_gmin_null8Lee_true8Lee", &chi2_gmin_null8Lee_true8Lee, "chi2_gmin_null8Lee_true8Lee/D" );

    int N_toy = 500;
        
    for(int itoy=1; itoy<=N_toy; itoy++) {
            
      if( itoy%max(N_toy/10,1)==0 ) {
	cout<<TString::Format(" ---> processing toy ( total cov ): %4.2f, %6d", itoy*1./N_toy, itoy)<<endl;
      }
      cout<<Form(" running %6d", itoy)<<endl;
      
      int status_fit = 0;
          
      /////////////////////////////////// null8sm, true8sm
      
      Lee_test->scaleF_Lee = 0;
      Lee_test->Set_Collapse();    
      Lee_test->Set_Variations(1);
      Lee_test->Set_toy_Variation(1);
    
      Lee_test->Minimization_Lee_strength_FullCov(0, 1);
      chi2_null_null8sm_true8sm = Lee_test->minimization_chi2;

      Lee_test->Minimization_Lee_strength_FullCov(1, 0);
      chi2_gmin_null8sm_true8sm = Lee_test->minimization_chi2;
      status_fit += Lee_test->minimization_status;
      
      /////////////////////////////////// null8Lee, true8Lee
      
      Lee_test->scaleF_Lee = 1;
      Lee_test->Set_Collapse();    
      Lee_test->Set_Variations(1);
      Lee_test->Set_toy_Variation(1);
    
      Lee_test->Minimization_Lee_strength_FullCov(1, 1);
      chi2_null_null8Lee_true8Lee = Lee_test->minimization_chi2;

      Lee_test->Minimization_Lee_strength_FullCov(1, 0);
      chi2_gmin_null8Lee_true8Lee = Lee_test->minimization_chi2;
      status_fit += Lee_test->minimization_status;
      
      ///////////////////////////////////
      
      if( status_fit!=0 ) continue;
      tree->Fill();
    }
    
    file_out->cd();
    tree->Write();
    file_out->Close();
    
  }

  //////////////////////////////////////////////// Sensitivity by Asimov sample

  if( 0 ) {

    ///////////////////////// reject SM
    
    Lee_test->scaleF_Lee = 1;
    Lee_test->Set_Collapse();
    
    Lee_test->Set_toy_Asimov();// use the Asimov sample as the input data for the fitting
    Lee_test->Minimization_Lee_strength_FullCov(0, 1);// (initial value, fix or not)

    double sigma_SM = sqrt( Lee_test->minimization_chi2 );
    cout<<TString::Format(" ---> Excluding  SM: %5.2f sigma", sigma_SM)<<endl;
    
    ///////////////////////// reject 1*LEE
    
    Lee_test->scaleF_Lee = 0;
    Lee_test->Set_Collapse();
    
    Lee_test->Set_toy_Asimov();// use the Asimov sample as the input data for the fitting
    Lee_test->Minimization_Lee_strength_FullCov(1, 1);// (initial value, fix or not)

    double sigma_Lee = sqrt( Lee_test->minimization_chi2 );
    cout<<TString::Format(" ---> Excluding LEE: %5.2f sigma", sigma_Lee)<<endl<<endl;;
    
  }

  ////////////////////////////////////////////////  Feldman-Cousins approach --> heavy computation cost

  if( 0 ) {
    
    /////////////// range: [low, hgh] with step
    
    double Lee_true_low = 0;
    double Lee_true_hgh = 3;
    double Lee_step     = 0.02;
    
    /////////////// dchi2 distribution 
    
    // int num_toy = 2;    
    // Lee_test->Exe_Feldman_Cousins(Lee_true_low, Lee_true_hgh, Lee_step, num_toy, ifile);

    /////////////// dchi2 of Asimov sample
    
    Lee_test->Exe_Fledman_Cousins_Asimov(Lee_true_low, Lee_true_hgh, Lee_step);

    /////////////// dchi2 of measured data
    
    Lee_test->Set_measured_data();    
    TMatrixD matrix_data_input_fc = Lee_test->matrix_data_newworld;    
    Lee_test->Exe_Fiedman_Cousins_Data( matrix_data_input_fc, Lee_true_low, Lee_true_hgh, Lee_step );
    
  }
  
  ////////////////////////////////////////////////////////////////////////////////////////

  cout<<endl;
  cout<<" ---> flag_syst_flux_Xs    "<<Lee_test->flag_syst_flux_Xs<<endl;
  cout<<" ---> flag_syst_detector   "<<Lee_test->flag_syst_detector<<endl;
  cout<<" ---> flag_syst_additional "<<Lee_test->flag_syst_additional<<endl;
  cout<<" ---> flag_syst_mc_stat    "<<Lee_test->flag_syst_mc_stat<<endl;  
  
  cout<<endl<<endl;
  cout<<" ---> Finish all the program"<<endl;
  cout<<endl<<endl;
  
  if( config_Lee::flag_display_graphics ) {
    cout<<endl<<" Enter Ctrl+c to end the program"<<endl;
    cout<<" Enter Ctrl+c to end the program"<<endl<<endl;
    
    theApp.Run();
  }
*/  
  return 0;
}
