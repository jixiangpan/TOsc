#ifndef TOsc_ana
#define TOsc_ana

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

/// minuit2
#include "Math/Functor.h"
#include "Minuit2/Minuit2Minimizer.h"

///////////////////////////////////////////////////////////////////////////////////////////////////////// TOsc

struct EventInfo {
  bool flag_FC;
  int pdg;
  double weight;
  double nueEtrue;
  double nueEreco;
  double baseline;
};

class TOsc {
 public:
  TOsc() {
    cout<<endl<<" ---> Hello TOscillation"<<endl<<endl;
    
    rand = new TRandom3(0);

    Osc_delta_m2_eV2  = 7.25;// neutrino-4, arXiv:2003.03199
    Osc_sin22theta_ee = 0;
    
    scaleF_POT = 1;
  }

  /////////////////////////////////////////////////////// data memeber
  /// all the "map_bin" index are from zero
  /// all the "map_ch" index are corresponding to the merge.root

  TRandom3 *rand;

  double Osc_delta_m2_eV2;
  double Osc_sin22theta_ee;

  // manually set the following two map
  map<int, int>ch_nue_from_intrinsic_sample;// must be only two terms, first is FC, second is PC
  map<int, int>zeroout_ch;
  
  vector<TString>eventlist_beforescale_file; // intrinsic nue at each run
  vector< vector<EventInfo> >eventlist_beforescale_runs;
  vector<TString>event_afterscale_file; // merge.root only for intrinsic nue MC at each run
  vector<double>event_scaleF_runs; // calcualted by events_afterscale_file/events_list_file
  // validation: sum (events_list_file * events_scaleF) = events_summation_afterscale_file

  map<int, TH1F*>map_h1f_basic_nue_binning;
  
  TMatrixD matrix_transform;
  int bins_oldworld;
  int bins_newworld;

  int channels_observation;
  map<int, map<int, double> >map_data_spectrum_ch_bin;
  map<int, double>map_data_spectrum_newworld_bin;

  map<int, TString>map_input_spectrum_ch_str;
  map<int, map<int, double> >map_input_spectrum_ch_bin;// this is the Pother
  //map<int, double>map_input_spectrum_oldworld_bin;
  map<int, map<int, double> >map_nue_intrinsic_noosc_spectrum_ch_bin;// this is the Pnue noosc, scaled by POT in Apply_POT_scaled()
  map<int, map<int, double> >map_nue_intrinsic_wiosc_spectrum_ch_bin;// this is the Pnue wiosc, scaled by POT in Apply_Oscillation()
  map<int, map<int, double> >map_pred_wiosc_ch_bin;// ---> this would be Pnue + Pother, set value in Apply_Oscillation()
  map<int, double>map_pred_wiosc_oldworld_bin;     // ---> this would be Pnue + Pother, set value in Apply_Oscillation()
    
  double scaleF_POT;

  TMatrixD matrix_data_newworld;
  TMatrixD matrix_pred_newworld;
  TMatrixD matrix_pred_oldworld;
  
  ////////////////////

  TH2D *h2_space;
  
  map<int, map<int, TMatrixD> >matrix_syst_frac_flux_before;
  map<int, map<int, TMatrixD> >matrix_syst_frac_geant_before;
  map<int, map<int, TMatrixD> >matrix_syst_frac_Xs_before;
  map<int, map<int, TMatrixD> >matrix_syst_frac_detector_before;
  map<int, map<int, TMatrixD> >matrix_syst_frac_additional_before;

  map<int, map<int, TMatrixD> >matrix_syst_frac_MCstat_after;

  TMatrixD matrix_syst_abs_flux_newworld;
  TMatrixD matrix_syst_abs_geant_newworld;
  TMatrixD matrix_syst_abs_Xs_newworld;
  TMatrixD matrix_syst_abs_detector_newworld;
  TMatrixD matrix_syst_abs_additional_newworld;
  TMatrixD matrix_syst_abs_MCstat_newworld;
  TMatrixD matrix_syst_abs_total_newworld;
  
  bool flag_syst_flux;
  bool flag_syst_geant;
  bool flag_syst_Xs;
  bool flag_syst_detector;
  bool flag_syst_additional;
  bool flag_syst_MCstat;

  ////////////////////

  TMatrixD matrix_dataFIT_newworld;

  int minimization_status;
  double minimization_chi2;
  double minimization_s22theta_val;
  double minimization_s22theta_err;
  double minimization_dm2_val;
  double minimization_dm2_err;
  
  /////////////////////////////////////////////////////// function memeber

  void Set_OscPars(double sin22theta_ee, double delta_m2_eV2)
  {    
    Osc_sin22theta_ee = sin22theta_ee;
    Osc_delta_m2_eV2 = delta_m2_eV2;
  }

  double ProbOsc(double nueEtrue, double baseline)
  {
    double weight = 1;
    weight = 1 - Osc_sin22theta_ee * pow(TMath::Sin(1.267 * Osc_delta_m2_eV2 * baseline/nueEtrue), 2);
    return weight;
  }

  void Set_Spectra_MatrixCov(TString eventlist_dir, TString event_summation_afterscale_file, TString centralvalue_noosc_file, TString flux_geant_Xs_file_dir, TString detector_file_dir, TString MCstat_file_dir);
  
  void Apply_POT_scaled();

  void Apply_Oscillation();
  // output: matrix_pred_oldworld, map_pred_wiosc_oldworld_bin, map_pred_wiosc_ch_bin, map_nue_intrinsic_wiosc_spectrum_ch_bin

  void Set_Collapse();

  void Set_Asimov2dataFIT() { matrix_dataFIT_newworld = matrix_pred_newworld; }

  void Set_data2dataFIT() { matrix_dataFIT_newworld = matrix_data_newworld; }
  
  void Minimization_OscPars_FullCov(double init_s22theta, double init_dm2, bool flag_fixed);
};


///////////////////////////////////////////////////////////////////////////////////////////////////////// TLee

class TLee {
public:
  TLee() {
    rand = new TRandom3(0);
    flag_individual_cov_newworld = true;
    flag_Lee_minimization_after_constraint = false;
  }

  /////////////////////////////////////////////////////// data memeber

  TRandom3 *rand;
  
  double scaleF_POT;
  double scaleF_Lee;

  bool flag_syst_flux_Xs;
  bool flag_syst_detector;
  bool flag_syst_additional;
  bool flag_syst_mc_stat;

  TString spectra_file;
  TString flux_Xs_directory;
  TString detector_directory;
  TString mc_directory;

  int channels_observation;
  int syst_cov_flux_Xs_begin;
  int syst_cov_flux_Xs_end;
  int syst_cov_mc_stat_begin;
  int syst_cov_mc_stat_end;
  
  /// bin index starts from 0, channel index from 1
  map<int, map<int, double> >map_input_spectrum_ch_bin;
  map<int, TString>map_input_spectrum_ch_str;
  map<int, double>map_input_spectrum_oldworld_bin;
  int bins_oldworld;
  int bins_newworld;

  map<int, int>map_Lee_ch;
  map<int, int>map_Lee_oldworld;

  map<int, map<int, double> >map_data_spectrum_ch_bin;
  map<int, double>map_data_spectrum_newworld_bin;
  TMatrixD matrix_data_newworld;
  
  TMatrixD matrix_input_cov_flux_Xs;  
  TMatrixD matrix_input_cov_detector;
  TMatrixD matrix_input_cov_additional;

  TMatrixD matrix_input_cov_flux;
  TMatrixD matrix_input_cov_Xs;
  map<int, TMatrixD>matrix_input_cov_detector_sub;
  
  map<int, TGraph*>gh_mc_stat_bin;

  TMatrixD matrix_transform;

  map<int, double>map_pred_spectrum_newworld_bin;
  TMatrixD matrix_pred_newworld;
  
  TMatrixD matrix_absolute_cov_oldworld;
  TMatrixD matrix_absolute_cov_newworld;

  bool flag_individual_cov_newworld;
  TMatrixD matrix_absolute_flux_cov_newworld;
  TMatrixD matrix_absolute_Xs_cov_newworld;
  TMatrixD matrix_absolute_detector_cov_newworld;
  map<int, TMatrixD>matrix_absolute_detector_sub_cov_newworld;
  TMatrixD matrix_absolute_mc_stat_cov_newworld;
  TMatrixD matrix_absolute_additional_cov_newworld;

  map<int, map<int, double> >map_toy_variation;
  
  map<int, double>map_fake_data;

  double val_GOF_noConstrain = 0;
  double val_GOF_wiConstrain = 0;
  int val_GOF_NDF = 0;
  
  int minimization_status;
  double minimization_chi2;
  double minimization_Lee_strength_val;
  double minimization_Lee_strength_err;
  int minimization_NDF;
  bool flag_Lee_minimization_after_constraint;
  
  /////////////////////////////////////////////////////// function member

  void Set_config_file_directory(TString spectra_file_, TString flux_Xs_directory_, TString detector_directory_, TString mc_directory_);

  //
  void Set_Spectra_MatrixCov();
  void Set_POT_implement();
  void Set_TransformMatrix();

  //
  void Set_Collapse();

  //
  void Plotting_systematics();

  // Y constrained by X
  int Exe_Goodness_of_fit(int num_Y, int num_X, TMatrixD matrix_pred, TMatrixD matrix_data, TMatrixD matrix_syst, int index);

  // target channels are constrained by support channels
  int Exe_Goodness_of_fit(vector<int>vc_target_chs, vector<int>vc_support_chs, int index);

  // target detailed channels are constrained by supported detailed channels
  int Exe_Goodness_of_fit_detailed(vector<int>vc_target_detailed_chs, vector<int>vc_support_detailed_chs, int index);
  
  // produceing pseudo-experiments
  void Set_Variations(int num_toy);

  // set input sample in to the minimization procedure
  void Set_toy_Asimov();
  void Set_toy_Variation(int itoy);
  void Set_measured_data();
  void Set_fakedata(TMatrixD matrix_fakedata);// 1, bins
  
  // minimization
  void Minimization_Lee_strength_FullCov(double Lee_initial_value, bool flag_fixed);
  
  // Feldman-Cousins approach
  void Exe_Feldman_Cousins(double Lee_true_low, double Lee_true_hgh, double step, int num_toy, int ifile);
  void Exe_Fledman_Cousins_Asimov(double Lee_true_low, double Lee_true_hgh, double step);
  void Exe_Fiedman_Cousins_Data(TMatrixD matrix_fakedata, double Lee_true_low, double Lee_true_hgh, double step);
};

#endif
