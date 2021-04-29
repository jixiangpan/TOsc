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

  /// manually set the following two map
  map<int, int>ch_nue_from_intrinsic_sample;// must be only two terms, first is FC, second is PC
  map<int, int>zeroout_ch;

  ///
  vector<TString>eventlist_beforescale_file; // intrinsic nue at each run
  vector< vector<EventInfo> >eventlist_beforescale_runs;
  vector<TString>event_afterscale_file; // merge.root only for intrinsic nue MC at each run
  vector<double>event_scaleF_runs; // calcualted by events_afterscale_file/events_list_file
  // validation: sum (events_list_file * events_scaleF) = events_summation_afterscale_file

  map<int, TH1F*>map_h1f_basic_nue_binning;

  TString g_flux_geant_Xs_file_dir;
  
  TMatrixD matrix_transform;
  int bins_oldworld;
  int bins_newworld;

  int channels_observation;
  map<int, map<int, double> >map_data_spectrum_ch_bin;
  map<int, double>map_data_spectrum_newworld_bin;

  map<int, TString>map_input_spectrum_ch_str;
  map<int, map<int, double> >map_input_spectrum_ch_bin;// this is the Pother
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
  
  bool flag_FIT_after_constraint;
    
  TMatrixD matrix_data_realdforFIT;
  TMatrixD matrix_pred_asimvforFIT;

  map<int, map<int,double> >map_toy_dataforFIT;
  
  TMatrixD matrix_dataFIT;
  TMatrixD matrix_predFIT;
  TMatrixD matrix_systFIT;

  int    minimization_status;
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

  void Set_Spectra_MatrixCov(TString eventlist_dir, TString event_summation_afterscale_file, TString centralvalue_noosc_file, TString flux_geant_Xs_file_dir);
  
  void Apply_POT_scaled();

  void Apply_Oscillation();
  // output: matrix_pred_oldworld, map_pred_wiosc_oldworld_bin, map_pred_wiosc_ch_bin, map_nue_intrinsic_wiosc_spectrum_ch_bin

  void Set_Collapse();

  void Set_EffectiveInput2FIT();

  void Set_Asimov2dataFIT() {
    matrix_dataFIT.Clear(); matrix_dataFIT.ResizeTo( 1, matrix_pred_asimvforFIT.GetNcols() );
    matrix_dataFIT = matrix_pred_asimvforFIT;
  }

  void Set_Toy2dataFIT(int itoy) {
    matrix_dataFIT.Clear(); matrix_dataFIT.ResizeTo( 1, matrix_pred_asimvforFIT.GetNcols() );
    for(int idx=0; idx<matrix_pred_asimvforFIT.GetNcols(); idx++) {
      matrix_dataFIT(0, idx) = map_toy_dataforFIT[itoy][idx];
    }
  }

  void Set_data2dataFIT() {
    matrix_dataFIT.Clear(); matrix_dataFIT.ResizeTo( 1, matrix_data_realdforFIT.GetNcols() );
    matrix_dataFIT = matrix_data_realdforFIT;
  }

  void Produce_Variations(int num_toy);
  
  void Minimization_OscPars_FullCov(double init_s22theta, double init_dm2, bool flag_fixed_all, bool flag_fixed_dm2);
};

#endif
