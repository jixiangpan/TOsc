#include "WCPLEEANA/TOsc.h"

#include "draw.icc"

/////////////////////////////////////////////////////////////////////////////////////////////////////////

namespace DataBase {
  double x1[101]={0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
                  11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
                  21, 22, 23, 24, 25, 26, 27, 28, 29, 30,
                  31, 32, 33, 34, 35, 36, 37, 38, 39, 40,
                  41, 42, 43, 44, 45, 46, 47, 48, 49, 50,
                  51, 52, 53, 54, 55, 56, 57, 58, 59, 60,
                  61, 62, 63, 64, 65, 66, 67, 68, 69, 70,
                  71, 72, 73, 74, 75, 76, 77, 78, 79, 80,
                  81, 82, 83, 84, 85, 86, 87, 88, 89, 90,
                  91, 92, 93, 94, 95, 96, 97, 98, 99, 100};
  
  double yl[101]={0, 0, 0, 0.856632, 1.70317, 2.51005, 3.32075, 4.14046, 4.9693, 5.80646, 6.65117,
                  7.5025, 8.35978, 9.22237, 10.0898, 10.9615, 11.8372, 12.7165, 13.5992, 14.4849, 15.3734,
                  16.2646, 17.1583, 18.0543, 18.9524, 19.8526, 20.7547, 21.6586, 22.5642, 23.4715, 24.3803,
                  25.2906, 26.2023, 27.1153, 28.0297, 28.9452, 29.8619, 30.7797, 31.6987, 32.6187, 33.5396,
                  34.4616, 35.3845, 36.3083, 37.2329, 38.1584, 39.0847, 40.0118, 40.9396, 41.8682, 42.7975,
                  43.7275, 44.6581, 45.5895, 46.5215, 47.454, 48.3873, 49.321, 50.2554, 51.1903, 52.1257,
                  53.0617, 53.9982, 54.9352, 55.8727, 56.8107, 57.7491, 58.6881, 59.6274, 60.5673, 61.5075,
                  62.4482, 63.3892, 64.3307, 65.2725, 66.2148, 67.1575, 68.1005, 69.0438, 69.9876, 70.9317,
                  71.8761, 72.8209, 73.766, 74.7114, 75.6572, 76.6033, 77.5497, 78.4964, 79.4434, 80.3907,
                  81.3383, 82.2862, 83.2342, 84.1827, 85.1314, 86.0804, 87.0296, 87.9791, 88.9288, 89.8788};

  double yh[101]={1.1478, 2.35971, 3.51917, 4.72422, 5.98186, 7.21064, 8.41858, 9.61053, 10.7896, 11.9582, 13.1179,
                  14.27, 15.4155, 16.5552, 17.6898, 18.8197, 19.9454, 21.0673, 22.1858, 23.3011, 24.4133,
                  25.5229, 26.6299, 27.7346, 28.837, 29.9374, 31.0358, 32.1322, 33.2271, 34.3201, 35.4117,
                  36.5017, 37.5904, 38.6776, 39.7635, 40.8483, 41.9318, 43.0141, 44.0955, 45.1757, 46.2549,
                  47.3331, 48.4104, 49.4868, 50.5623, 51.637, 52.7108, 53.7839, 54.8561, 55.9277, 56.9985,
                  58.0686, 59.1381, 60.2068, 61.275, 62.3425, 63.4094, 64.4757, 65.5415, 66.6066, 67.6713,
                  68.7354, 69.7989, 70.862, 71.9246, 72.9866, 74.0483, 75.1094, 76.1701, 77.2304, 78.2902,
                  79.3496, 80.4085, 81.4672, 82.5253, 83.5831, 84.6406, 85.6976, 86.7542, 87.8105, 88.8665,
                  89.9221, 90.9774, 92.0323, 93.0869, 94.1411, 95.1951, 96.2488, 97.3021, 98.3552, 99.4079,
                  100.46, 101.513, 102.564, 103.616, 104.667, 105.718, 106.769, 107.82, 108.87, 109.92};
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//// ccc
void TOsc::Minimization_OscPars_FullCov(double init_s22theta, double init_dm2, bool flag_fixed)
{
  TString roostr = "";
  
  minimization_status = (-100);
  minimization_chi2   = (-100);
    
  if( !flag_fixed ) {
    ROOT::Minuit2::Minuit2Minimizer min_Osc( ROOT::Minuit2::kMigrad );
    min_Osc.SetPrintLevel(0);
    min_Osc.SetStrategy(1); //0- cursory, 1- default, 2- thorough yet no more successful
    min_Osc.SetMaxFunctionCalls(500000);
    min_Osc.SetMaxIterations(500000);
    min_Osc.SetTolerance(1e-4); // tolerance*2e-3 = edm precision
    min_Osc.SetPrecision(1e-18); //precision in the target function

    /// set fitting parameters
    ROOT::Math::Functor Chi2Functor_Osc( [&](const double *par) {// FCN
	TString roostr = "";
	double chi2 = 0;
	double s22theta = par[0];
	double dm2 = par[1];

	Set_OscPars(s22theta, dm2);
	Set_Collapse();      

	/////////
	TMatrixD matrix_data = matrix_dataFIT_newworld;          
	TMatrixD matrix_pred = matrix_pred_newworld;      
	TMatrixD matrix_cov_syst = matrix_syst_abs_total_newworld;;
	int bins_fit = matrix_cov_syst.GetNrows();
      
	for(int idx=0; idx<bins_fit; idx++) {
	  double val_stat_cov = 0;        
	  double val_data = matrix_data(0, idx);
	  double val_pred = matrix_pred(0, idx);
        
	  if( val_data==0 ) { val_stat_cov = val_pred/2; }
	  else {
	    if( val_pred!=0 ) val_stat_cov = 3./( 1./val_data + 2./val_pred );
	    else val_stat_cov = val_data;
	  }
	
	  matrix_cov_syst(idx, idx) += val_stat_cov;
	  if( matrix_cov_syst(idx, idx)==0 ) matrix_cov_syst(idx, idx) = 1e-6;
	}

	/////////
	TMatrixD matrix_cov_total = matrix_cov_syst;
	TMatrixD matrix_cov_total_inv = matrix_cov_total; matrix_cov_total_inv.Invert();
      
	TMatrixD matrix_delta = matrix_pred - matrix_data;
	TMatrixD matrix_delta_T = matrix_delta.T(); matrix_delta.T();
   
	TMatrixD matrix_chi2 = matrix_delta * matrix_cov_total_inv *matrix_delta_T;
	chi2 = matrix_chi2(0,0);           
      
	return chi2;
      
      },// end of FCN
      2 // number of fitting parameters
      );
  
    min_Osc.SetFunction(Chi2Functor_Osc);
  
    min_Osc.SetVariable( 0, "sin_2_2theta", init_s22theta, 1e-3);
    min_Osc.SetVariable( 1, "delta_m2", init_dm2, 1e-3);
  
    //min_Osc.SetLowerLimitedVariable(0, "sin_2_2theta", init_s22theta, 1e-3, 0);
    //min_Osc.SetUpperLimitedVariable(0, "sin_2_2theta", init_s22theta, 1e-3, 1);
    //min_Osc.SetVariableLowerLimit(0, 0);
    //min_Osc.SetVariableUpperLimit(0, 1);

    min_Osc.SetLimitedVariable(0, "sin_2_2theta", init_s22theta, 1e-3, 0, 1);
    min_Osc.SetLowerLimitedVariable(1, "delta_m2", init_dm2, 1e-3, 0);  
  
    // if( flag_fixed ) {
    //   min_Osc.SetFixedVariable( 0, "sin_2_2theta", init_s22theta );
    //   min_Osc.SetFixedVariable( 1, "delta_m2", init_dm2 );
    // }
  
    /// do the minimization
    min_Osc.Minimize();
    int status_Osc = min_Osc.Status();
    const double *par_Osc = min_Osc.X();
    const double *par_Osc_err = min_Osc.Errors();

    if( status_Osc!=0 ) {
      cerr<<" -----------> Oscillation pars fitting failed "<<endl;
      minimization_status = status_Osc;
    }

    minimization_status = status_Osc;
    minimization_chi2 = min_Osc.MinValue();
    minimization_s22theta_val = par_Osc[0];
    minimization_s22theta_err = par_Osc_err[0];
    minimization_dm2_val = par_Osc[1];
    minimization_dm2_err = par_Osc_err[1];
  }

  /////////////////////////////////////////////
  /////////////////////////////////////////////
  /////////////////////////////////////////////
  
  if( flag_fixed ) {
    
    double chi2 = 0;
    double s22theta = init_s22theta;
    double dm2      = init_dm2;

    Set_OscPars(s22theta, dm2);
    Set_Collapse();      

    /////////    
    TMatrixD matrix_data = matrix_dataFIT_newworld;          
    TMatrixD matrix_pred = matrix_pred_newworld;      
    TMatrixD matrix_cov_syst = matrix_syst_abs_total_newworld;;
    int bins_fit = matrix_cov_syst.GetNrows();
      
    for(int idx=0; idx<bins_fit; idx++) {
      double val_stat_cov = 0;        
      double val_data = matrix_data(0, idx);
      double val_pred = matrix_pred(0, idx);
        
      if( val_data==0 ) { val_stat_cov = val_pred/2; }
      else {
	if( val_pred!=0 ) val_stat_cov = 3./( 1./val_data + 2./val_pred );
	else val_stat_cov = val_data;
      }
	
      matrix_cov_syst(idx, idx) += val_stat_cov;
      if( matrix_cov_syst(idx, idx)==0 ) matrix_cov_syst(idx, idx) = 1e-6;
    }

    /////////
    TMatrixD matrix_cov_total = matrix_cov_syst;
    TMatrixD matrix_cov_total_inv = matrix_cov_total; matrix_cov_total_inv.Invert();
      
    TMatrixD matrix_delta = matrix_pred - matrix_data;
    TMatrixD matrix_delta_T = matrix_delta.T(); matrix_delta.T();
   
    TMatrixD matrix_chi2 = matrix_delta * matrix_cov_total_inv *matrix_delta_T;
    chi2 = matrix_chi2(0,0);   

    minimization_chi2 = chi2;
	
  }// flag_fixed
 
}

//// ccc
void TOsc::Set_Collapse()
{
  Apply_Oscillation();
  
  /////////////////////////////////////// central value

  matrix_pred_newworld.ResizeTo(1, bins_newworld);
  matrix_pred_newworld.Clear();
  matrix_pred_newworld.ResizeTo(1, bins_newworld);
  matrix_pred_newworld = matrix_pred_oldworld * matrix_transform;

  /////////////////////////////////////// covariance matrix for the systematics

  int global_bin = h2_space->FindBin( log10(Osc_sin22theta_ee), log10(Osc_delta_m2_eV2) );
  int xbin(0), ybin(0), zbin(0);
  h2_space->GetBinXYZ(global_bin, xbin, ybin, zbin);
  if( xbin<1  ) xbin = 1;
  if( xbin>h2_space->GetNbinsX() ) xbin = h2_space->GetNbinsX();
  if( ybin<1  ) ybin = 1;
  if( ybin>h2_space->GetNbinsY() ) ybin = h2_space->GetNbinsY();

  if( fabs(Osc_sin22theta_ee)<1e-2 ) {// no oscillation case
    xbin = 0;
    ybin = 0;
    //cout<<" Using the syst covariance assuming no oscillation"<<endl;
  }
  
  ////////////////

  TMatrixD matrix_syst_abs_flux_before = matrix_syst_frac_flux_before[xbin][ybin];
  TMatrixD matrix_syst_abs_geant_before = matrix_syst_frac_geant_before[xbin][ybin];
  TMatrixD matrix_syst_abs_Xs_before = matrix_syst_frac_Xs_before[xbin][ybin];
  TMatrixD matrix_syst_abs_detector_before = matrix_syst_frac_detector_before[xbin][ybin];
  TMatrixD matrix_syst_abs_additional_before = matrix_syst_frac_additional_before[xbin][ybin];
  
  for(int idx=0; idx<bins_oldworld; idx++) {
    for(int jdx=0; jdx<bins_oldworld; jdx++) {
      double cv_i = matrix_pred_oldworld(0, idx);
      double cv_j = matrix_pred_oldworld(0, jdx);
      double cov_ij_frac = 0;

      /// flux
      cov_ij_frac = matrix_syst_abs_flux_before(idx, jdx);
      matrix_syst_abs_flux_before(idx, jdx) = cov_ij_frac * cv_i * cv_j;
   
      /// geant
      cov_ij_frac = matrix_syst_abs_geant_before(idx, jdx);
      matrix_syst_abs_geant_before(idx, jdx) = cov_ij_frac * cv_i * cv_j;

      /// Xs
      cov_ij_frac = matrix_syst_abs_Xs_before(idx, jdx);
      matrix_syst_abs_Xs_before(idx, jdx) = cov_ij_frac * cv_i * cv_j;

      /// detector
      cov_ij_frac = matrix_syst_abs_detector_before(idx, jdx);
      matrix_syst_abs_detector_before(idx, jdx) = cov_ij_frac * cv_i * cv_j;

      /// additional
      cov_ij_frac = matrix_syst_abs_additional_before(idx, jdx);
      matrix_syst_abs_additional_before(idx, jdx) = cov_ij_frac * cv_i * cv_j;               
    }// jdx
  }// idx

  
  TMatrixD matrix_transform_T = matrix_transform.T(); matrix_transform.T();

  /// flux
  matrix_syst_abs_flux_newworld.Clear();
  matrix_syst_abs_flux_newworld.ResizeTo(bins_newworld, bins_newworld);
  matrix_syst_abs_flux_newworld = matrix_transform_T * matrix_syst_abs_flux_before * matrix_transform;
    
  /// geant
  matrix_syst_abs_geant_newworld.Clear();
  matrix_syst_abs_geant_newworld.ResizeTo(bins_newworld, bins_newworld);
  matrix_syst_abs_geant_newworld = matrix_transform_T * matrix_syst_abs_geant_before * matrix_transform;
 
  /// Xs
  matrix_syst_abs_Xs_newworld.Clear();
  matrix_syst_abs_Xs_newworld.ResizeTo(bins_newworld, bins_newworld);
  matrix_syst_abs_Xs_newworld = matrix_transform_T * matrix_syst_abs_Xs_before * matrix_transform;
  
  /// detector
  matrix_syst_abs_detector_newworld.Clear();
  matrix_syst_abs_detector_newworld.ResizeTo(bins_newworld, bins_newworld);
  matrix_syst_abs_detector_newworld = matrix_transform_T * matrix_syst_abs_detector_before * matrix_transform;
    
  /// additional
  matrix_syst_abs_additional_newworld.Clear();
  matrix_syst_abs_additional_newworld.ResizeTo(bins_newworld, bins_newworld);
  matrix_syst_abs_additional_newworld = matrix_transform_T * matrix_syst_abs_additional_before * matrix_transform;
         
  /// MCstat
  matrix_syst_abs_MCstat_newworld.Clear();
  matrix_syst_abs_MCstat_newworld.ResizeTo(bins_newworld, bins_newworld);
  matrix_syst_abs_MCstat_newworld = matrix_syst_frac_MCstat_after[xbin][ybin];

  /// total
  matrix_syst_abs_total_newworld.Clear();
  matrix_syst_abs_total_newworld.ResizeTo(bins_newworld, bins_newworld);
  
  if(flag_syst_flux)       matrix_syst_abs_total_newworld += matrix_syst_abs_flux_newworld;
  if(flag_syst_geant)      matrix_syst_abs_total_newworld += matrix_syst_abs_geant_newworld;
  if(flag_syst_Xs)         matrix_syst_abs_total_newworld += matrix_syst_abs_Xs_newworld;
  if(flag_syst_detector)   matrix_syst_abs_total_newworld += matrix_syst_abs_detector_newworld;
  if(flag_syst_additional) matrix_syst_abs_total_newworld += matrix_syst_abs_additional_newworld;
  if(flag_syst_MCstat)     matrix_syst_abs_total_newworld += matrix_syst_abs_MCstat_newworld;

}

//// ccc
void TOsc::Apply_Oscillation()
{  
  if( (int)ch_nue_from_intrinsic_sample.size()!=2 ) {cerr<<" ch_nue_from_intrinsic_sample.size()!=2"<<endl; exit(1); }

  ///////////////////////////  
  int temp_FC = 0;
  int temp_PC = 0;
  for(auto it_int_map=ch_nue_from_intrinsic_sample.begin(); it_int_map!=ch_nue_from_intrinsic_sample.end(); it_int_map++) {
    int ich = it_int_map->first;
    if( it_int_map==ch_nue_from_intrinsic_sample.begin() ) temp_FC = ich;
    else temp_PC = ich;

    /// initialize the nue intrinsic with Oscillation
    map_h1f_basic_nue_binning[ich]->Reset();
    map_nue_intrinsic_wiosc_spectrum_ch_bin.clear();
  }
  
  int num_eventlist_file = eventlist_beforescale_runs.size();
  for(int idx=0; idx<num_eventlist_file; idx++) {
    int num_events = eventlist_beforescale_runs.at(idx).size();
    double scaleF_eachrun = event_scaleF_runs.at(idx);
    
    for(int ievent=0; ievent<num_events; ievent++) {
      double nueEtrue = eventlist_beforescale_runs.at(idx).at(ievent).nueEtrue;      
      double baseline = eventlist_beforescale_runs.at(idx).at(ievent).baseline;
      double nueEreco = eventlist_beforescale_runs.at(idx).at(ievent).nueEreco;
      double weight = eventlist_beforescale_runs.at(idx).at(ievent).weight;
      bool flag_FC = eventlist_beforescale_runs.at(idx).at(ievent).flag_FC;

      double prob_val = ProbOsc(nueEtrue, baseline);
      double weight_total = prob_val * weight * scaleF_eachrun * scaleF_POT;
      
      int ich = 0;
      if( flag_FC ) ich = temp_FC;
      else ich = temp_PC;

      map_h1f_basic_nue_binning[ich]->Fill( nueEreco, weight_total );
      
    }// for(int ievent=0; ievent<num_events; ievent++)    
  }// for(int idx=0; idx<num_eventlist_file; idx++)
  
  for(auto it_int_map=ch_nue_from_intrinsic_sample.begin(); it_int_map!=ch_nue_from_intrinsic_sample.end(); it_int_map++) {
    int ich = it_int_map->first;
    int bins = map_h1f_basic_nue_binning[ich]->GetNbinsX()+1;

    for(int ibin=1; ibin<=bins; ibin++) {
      double content = map_h1f_basic_nue_binning[ich]->GetBinContent(ibin);
      map_nue_intrinsic_wiosc_spectrum_ch_bin[ich][ibin-1] = content;
    }
  }

  ////////////////////////////
  map_pred_wiosc_ch_bin.clear();  

  for(auto it_it_map=map_input_spectrum_ch_bin.begin(); it_it_map!=map_input_spectrum_ch_bin.end(); it_it_map++) {
    int ich = it_it_map->first; int bins = it_it_map->second.size();
    for(int idx=0; idx<bins; idx++) {
      double val_Pother = map_input_spectrum_ch_bin[ich][idx];
      double val_Pnue = 0;
      if( ch_nue_from_intrinsic_sample.find(ich)!=ch_nue_from_intrinsic_sample.end() )
	{ val_Pnue = map_nue_intrinsic_wiosc_spectrum_ch_bin[ich][idx]; }
      map_pred_wiosc_ch_bin[ich][idx] = val_Pother + val_Pnue;      
    }
  }

  map_pred_wiosc_oldworld_bin.clear();
  matrix_pred_oldworld.Clear();
  matrix_pred_oldworld.ResizeTo(1, bins_oldworld);
  
  int line_temp_global = -1;
  for(auto it_it_map=map_input_spectrum_ch_bin.begin(); it_it_map!=map_input_spectrum_ch_bin.end(); it_it_map++) {
    int ich = it_it_map->first; int bins = it_it_map->second.size();
    for(int idx=0; idx<bins; idx++) {
      line_temp_global++;
      double content = map_pred_wiosc_ch_bin[ich][idx];
      map_pred_wiosc_oldworld_bin[line_temp_global] = content;
      matrix_pred_oldworld(0, line_temp_global) = content;
    }
  }

  //////////////////////////// validation during the development
  // if( 1 ) {
  //   //TFile *file_test = new TFile("./data_nue_beforeafter_scale/merge_dm2_7d25_s22t_0d26.root", "read");
  //   //TFile *file_test = new TFile("./data_nue_beforeafter_scale/merge_noosc.root", "read");
  //   for(auto it_it_map=map_input_spectrum_ch_bin.begin(); it_it_map!=map_input_spectrum_ch_bin.end(); it_it_map++) {
  //     int ich = it_it_map->first; int bins = it_it_map->second.size();
  //     if( zeroout_ch.find(ich)!=zeroout_ch.end() ) continue;      
  //     TH1F *h1f_temp = (TH1F*)file_test->Get(TString::Format("histo_%d", ich));
  //     for(int idx=0; idx<bins; idx++) {
  // 	double my_result = map_pred_wiosc_ch_bin[ich][idx];
  // 	double wcpframe_result = h1f_temp->GetBinContent(idx+1);
  // 	if( fabs(my_result-wcpframe_result)>1e-4 ) { cerr<<" fabs(my_result-wcpframe_result)>1e-4"<<endl; exit(1); }
  //     }
  //   }    
  // }  
  //////////////////////////// end the validation during the development
  
}

//// ccc
void TOsc::Apply_POT_scaled()
{
  cout<<endl<<" ---> Apply_POT_scaled"<<endl<<endl;

  //////////////// data
  for( auto it_it_map=map_data_spectrum_ch_bin.begin(); it_it_map!=map_data_spectrum_ch_bin.end(); it_it_map++) {
    int ich = it_it_map->first; int bins = it_it_map->second.size();
    for(int idx=0; idx<bins; idx++) map_data_spectrum_ch_bin[ich][idx] *= scaleF_POT;
  }

  matrix_data_newworld.Clear();
  matrix_data_newworld.ResizeTo(1, bins_newworld);
  for( auto it_double_map=map_data_spectrum_newworld_bin.begin(); it_double_map!= map_data_spectrum_newworld_bin.end(); it_double_map++) {
    int idx = it_double_map->first;
    map_data_spectrum_newworld_bin[idx] *= scaleF_POT;
    matrix_data_newworld(0, idx) = map_data_spectrum_newworld_bin[idx];
  }

  //////////////// pred
  for( auto it_it_map=map_input_spectrum_ch_bin.begin(); it_it_map!=map_input_spectrum_ch_bin.end(); it_it_map++ ) {
    int ich = it_it_map->first; int bins = it_it_map->second.size();
    for(int idx=0; idx<bins; idx++) map_input_spectrum_ch_bin[ich][idx] *= scaleF_POT;
  }

  for( auto it_it_map=map_nue_intrinsic_noosc_spectrum_ch_bin.begin(); it_it_map!=map_nue_intrinsic_noosc_spectrum_ch_bin.end(); it_it_map++) {
    int ich = it_it_map->first; int bins = it_it_map->second.size();
    for(int idx=0; idx<bins; idx++) map_nue_intrinsic_noosc_spectrum_ch_bin[ich][idx] *= scaleF_POT;
  }

  //////////////// MCstat
  double scaleF_POT2 = scaleF_POT * scaleF_POT;  
  for(int xbin=0; xbin<=h2_space->GetNbinsX(); xbin++) {
    for(int ybin=0; ybin<=h2_space->GetNbinsY(); ybin++) {
      if( ybin==0 && xbin!=0 ) continue;
      if( xbin==0 && ybin!=0 ) continue;
      
      for(int idx=0; idx<bins_newworld; idx++) {
	matrix_syst_frac_MCstat_after[xbin][ybin] *= scaleF_POT2;
      }
    }
  }
  
}

//// ccc
void TOsc::Set_Spectra_MatrixCov(TString eventlist_dir, TString event_summation_afterscale_file, TString centralvalue_noosc_file, TString flux_geant_Xs_file_dir, TString detector_file_dir, TString MCstat_file_dir)
{
  cout<<endl<<" ---> Set_Spectra_MatrixCov"<<endl<<endl;

  TString roostr = "";
  
  ////////////////////////////////////////
  
  ch_nue_from_intrinsic_sample[1] = 1;// FC
  ch_nue_from_intrinsic_sample[2] = 1;// PC
  
  cout<<" Read eventlist"<<endl;

  TString str_read_eventlist_filename = eventlist_dir+"eventlist_filename.txt";
  ifstream read_eventlist_filename(str_read_eventlist_filename, ios::in);
  if(!read_eventlist_filename) { cerr<<" No eventlist_filename.txt"<<endl; exit(1); }
  while( !read_eventlist_filename.eof() ) {
    TString str_before, str_after;
    read_eventlist_filename>>str_before>>str_after;
    if( str_before=="end" ) break;
    cout<<" ---> eventlist input file: "<<str_before<<"\t"<<str_after<<endl;
    eventlist_beforescale_file.push_back( eventlist_dir+str_before );
    event_afterscale_file.push_back( eventlist_dir+str_after );
  }

  //////
  
  int num_eventlist_file = eventlist_beforescale_file.size();
  cout<<" ---> eventlist input file num: "<<num_eventlist_file<<endl;

  for(int idx=0; idx<num_eventlist_file; idx++) {

    TFile *file_afterscale = new TFile( event_afterscale_file.at(idx), "read" );
    roostr = TString::Format("histo_%d", ch_nue_from_intrinsic_sample.begin()->first );
    TH1F *h1_nue_FC_afterscale = (TH1F*)file_afterscale->Get(roostr);
    int temp_bins = h1_nue_FC_afterscale->GetNbinsX();
    double temp_xlow = h1_nue_FC_afterscale->GetXaxis()->GetBinLowEdge(1);
    double temp_xhgh = h1_nue_FC_afterscale->GetXaxis()->GetBinUpEdge(temp_bins);
    double integral_afterscale = h1_nue_FC_afterscale->Integral();
    delete h1_nue_FC_afterscale;
    delete file_afterscale;

    TH1F *h1_FC_beforescale = new TH1F("h1_FC_beforescale", "", temp_bins, temp_xlow, temp_xhgh);    
    
    //////
    vector<EventInfo>eventlist_eachrun;
    
    TString eventlist_file = eventlist_beforescale_file.at(idx);
    ifstream read_eventlist_file_aa(eventlist_file, ios::in);
    string str_read_eventlist_file_aa = "";
    int line_read_eventlist_file_aa = 0;
    while( getline(read_eventlist_file_aa, str_read_eventlist_file_aa) ) line_read_eventlist_file_aa++;
    cout<<" processing file "<<eventlist_file<<"\t"<<line_read_eventlist_file_aa<<endl;
    
    ifstream read_eventlist_file_ab(eventlist_file, ios::in);
    for(int idx=1; idx<=line_read_eventlist_file_aa; idx++) {
      int pdg_val(0); double weight_val(0), nueEtrue(0), nueEreco(0), baseline(0); TString ch_name = "";
      read_eventlist_file_ab>>pdg_val>>weight_val>>nueEtrue>>nueEreco>>baseline>>ch_name;

      bool flag_FC = false;
      if( ch_name.Contains("_FC_") ) flag_FC = true;
      
      EventInfo eventinfo;
      eventinfo.flag_FC  = flag_FC;
      eventinfo.pdg      = pdg_val;
      eventinfo.weight   = weight_val;
      eventinfo.nueEtrue = nueEtrue;
      eventinfo.nueEreco = nueEreco;
      eventinfo.baseline = baseline;
      eventlist_eachrun.push_back( eventinfo );

      if( flag_FC ) h1_FC_beforescale->Fill( nueEreco, weight_val );      
    }// for(int idx=1; idx<=line_read_eventlist_file_aa; idx++)

    double integral_beforescale = h1_FC_beforescale->Integral();
    if( integral_beforescale==0 ) { cerr<<" integral_beforescale=0"<<endl; exit(1);  };
    double temp_scaleF = integral_afterscale/integral_beforescale;
    cout<<" processing file scaleF: "<<TString::Format("%8.6f", temp_scaleF)<<endl;
    delete h1_FC_beforescale;

    eventlist_beforescale_runs.push_back( eventlist_eachrun );
    event_scaleF_runs.push_back( temp_scaleF );    
  }// for(int idx=0; idx<num_eventlist_file; idx++)
  
  //////
  
  cout<<" ---> check eventlist after scale"<<endl;
    
  TFile *file_event_summation_afterscale = new TFile(event_summation_afterscale_file, "read");

  TH1F *h1_nue_FC_summation;
  TH1F *h1_nue_PC_summation;  
  for(auto it_double_map=ch_nue_from_intrinsic_sample.begin(); it_double_map!=ch_nue_from_intrinsic_sample.end(); it_double_map++) {
    int ich = it_double_map->first;

    if( it_double_map==ch_nue_from_intrinsic_sample.begin() ) {
      h1_nue_FC_summation = (TH1F*)file_event_summation_afterscale->Get( TString::Format("histo_%d", ich) );
      map_h1f_basic_nue_binning[ich] = (TH1F*)h1_nue_FC_summation->Clone(TString::Format("map_h1f_basic_nue_binning_%d", ich));
    }
    else {
      h1_nue_PC_summation = (TH1F*)file_event_summation_afterscale->Get( TString::Format("histo_%d", ich) );
      map_h1f_basic_nue_binning[ich] = (TH1F*)h1_nue_PC_summation->Clone(TString::Format("map_h1f_basic_nue_binning_%d", ich));
    }
    
    TH1F* h1f_temp = (TH1F*)file_event_summation_afterscale->Get( TString::Format("histo_%d", ich) );
    for(int ibin=1; ibin<=h1f_temp->GetNbinsX()+1; ibin++) {
      double content = h1f_temp->GetBinContent(ibin);
      map_nue_intrinsic_noosc_spectrum_ch_bin[ich][ibin-1] = content;
    }
  }
  
  ///
  TH1F *h1_nue_FC_subadd = (TH1F*)h1_nue_FC_summation->Clone("h1_nue_FC_subadd");
  h1_nue_FC_subadd->Reset();
  TH1F *h1_nue_PC_subadd = (TH1F*)h1_nue_PC_summation->Clone("h1_nue_PC_subadd");
  h1_nue_PC_subadd->Reset();

  for(int idx=0; idx<num_eventlist_file; idx++) {
    int num_events = eventlist_beforescale_runs.at(idx).size();
    double scaleF_eachrun = event_scaleF_runs.at(idx);
    
    for(int ievent=0; ievent<num_events; ievent++) {
      //double nueEtrue = eventlist_beforescale_runs.at(idx).at(ievent).nueEtrue;      
      //double baseline = eventlist_beforescale_runs.at(idx).at(ievent).baseline;
      double nueEreco = eventlist_beforescale_runs.at(idx).at(ievent).nueEreco;
      double weight = eventlist_beforescale_runs.at(idx).at(ievent).weight;
      bool flag_FC = eventlist_beforescale_runs.at(idx).at(ievent).flag_FC;

      if( flag_FC ) h1_nue_FC_subadd->Fill( nueEreco, weight*scaleF_eachrun );
      else h1_nue_PC_subadd->Fill( nueEreco, weight*scaleF_eachrun );
    }// for(int ievent=0; ievent<num_events; ievent++)
    
  }// for(int idx=0; idx<num_eventlist_file; idx++)

  for(int ibin=1; ibin<=h1_nue_FC_subadd->GetNbinsX(); ibin++) {
    double content_summation_FC = h1_nue_FC_summation->GetBinContent(ibin);
    double content_subadd_FC = h1_nue_FC_subadd->GetBinContent(ibin);
    double content_summation_PC = h1_nue_PC_summation->GetBinContent(ibin);
    double content_subadd_PC = h1_nue_PC_subadd->GetBinContent(ibin);
    if( fabs(content_summation_FC-content_subadd_FC)>1e-4 ) 
      { cerr<<" content_summation_FC != content_subadd_FC"<<endl; exit(1); }
    if( fabs(content_summation_PC-content_subadd_PC)>1e-4 ) 
      { cerr<<" content_summation_PC != content_subadd_PC"<<endl; exit(1); }   
  }
  
  ////////////////////////////////////////
  
  cout<<endl<<" Read centralvalue_noosc_file"<<endl;

  TFile *file_centralvalue_noosc = new TFile(centralvalue_noosc_file, "read");

  ////////
  TMatrixD *mat_collapse = (TMatrixD*)file_centralvalue_noosc->Get("mat_collapse");

  bins_oldworld = mat_collapse->GetNrows();
  bins_newworld = mat_collapse->GetNcols();
  
  matrix_transform.Clear();
  matrix_transform.ResizeTo( mat_collapse->GetNrows(), mat_collapse->GetNcols() );
  matrix_transform = (*mat_collapse);  

  ///////
  
  cout<<endl<<" Observations:"<<endl;

  int line_global_channels_observation = -1;
  for(int ich=1; ich<=channels_observation; ich++) {
    roostr = TString::Format("hdata_obsch_%d", ich);
    TH1F *h1_spectrum = (TH1F*)file_centralvalue_noosc->Get(roostr);
    cout<<Form(" %2d  %-20s   bin-num %2d", ich, roostr.Data(), h1_spectrum->GetNbinsX()+1)<<endl;
       
    for(int ibin=1; ibin<=h1_spectrum->GetNbinsX()+1; ibin++) {
      line_global_channels_observation++;
      double content = h1_spectrum->GetBinContent(ibin);
      map_data_spectrum_ch_bin[ich][ibin-1] = content;
      map_data_spectrum_newworld_bin[line_global_channels_observation] = map_data_spectrum_ch_bin[ich][ibin-1];
    }// ibin
  }// ich

  if( line_global_channels_observation+1!=bins_newworld )
    { cerr<<" line_global_channels_observation+1!=bins_newworld"<<endl; exit(1); }
  
  ///////
  
  cout<<endl<<" Predictions:"<<endl;
   
  map_input_spectrum_ch_str[1] = "nueCC_FC_norm";
  map_input_spectrum_ch_str[2] = "nueCC_PC_norm";
  map_input_spectrum_ch_str[3] = "numuCC_FC_norm";
  map_input_spectrum_ch_str[4] = "numuCC_PC_norm";
  map_input_spectrum_ch_str[5] = "CCpi0_FC_norm";
  map_input_spectrum_ch_str[6] = "CCpi0_PC_norm";
  map_input_spectrum_ch_str[7] = "NCpi0_norm";
  map_input_spectrum_ch_str[8] = "Lee_FC";
  map_input_spectrum_ch_str[9] = "Lee_PC"; 
  map_input_spectrum_ch_str[10]= "nueCC_FC_ext";
  map_input_spectrum_ch_str[11]= "nueCC_PC_ext";
  map_input_spectrum_ch_str[12]= "numuCC_FC_ext";
  map_input_spectrum_ch_str[13]= "numuCC_PC_ext";
  map_input_spectrum_ch_str[14]= "CCpi0_FC_ext";
  map_input_spectrum_ch_str[15]= "CCpi0_PC_ext";
  map_input_spectrum_ch_str[16]= "NCpi0_ext";
  
  zeroout_ch[8] = 1;
  zeroout_ch[9] = 1;
  for(auto it_zero=zeroout_ch.begin(); it_zero!=zeroout_ch.end(); it_zero++) 
    cout<<" ---> zeroout pred-ch: "<< it_zero->first <<endl;

  for(int ich=1; ich<=(int)map_input_spectrum_ch_str.size(); ich++) {
    roostr = TString::Format("histo_%d", ich);
    TH1F *h1_spectrum = (TH1F*)file_centralvalue_noosc->Get(roostr);
    int bins = h1_spectrum->GetNbinsX() + 1;    
    cout<<Form(" %2d  %-20s   bin-num %2d", ich, map_input_spectrum_ch_str[ich].Data(), bins)<<endl;

    for(int ibin=1; ibin<=bins; ibin++) {
      double content = h1_spectrum->GetBinContent(ibin);
      ////// case for no useful ch
      if( zeroout_ch.find(ich)!=zeroout_ch.end() ) { content = 0; }      
      map_input_spectrum_ch_bin[ich][ibin-1] = content;
    }    
  }

  ////// minus Pnue to get Pother
  for(auto it_map_input=map_input_spectrum_ch_bin.begin(); it_map_input!=map_input_spectrum_ch_bin.end(); it_map_input++) {
    int ich = it_map_input->first;
    int temp_size = it_map_input->second.size();
    if( ch_nue_from_intrinsic_sample.find(ich)!=ch_nue_from_intrinsic_sample.end() ) {
      for(int idx=0; idx<temp_size; idx++) {
	map_input_spectrum_ch_bin[ich][idx] -= map_nue_intrinsic_noosc_spectrum_ch_bin[ich][idx];	
      }
    }
  }
      
  //////
  int line_global_pred_input = -1;
  for(auto it_map_input=map_input_spectrum_ch_bin.begin(); it_map_input!=map_input_spectrum_ch_bin.end(); it_map_input++) {
    for(int idx=0; idx<(int)it_map_input->second.size(); idx++) line_global_pred_input++;
  }  
  if( line_global_pred_input+1!=bins_oldworld ) { cerr<<" line_global_pred_input+1!=bins_oldworld"<<endl; exit(1); }

  //////////////////////////////////////
  
  cout<<endl<<" ---> set the inputs of systematics covariance matrix"<<endl<<endl;
  
  ////// X: sin22t14, 1e-2 -> 1   ---> "log10()" ---> -2 -> 0
  ////// Y: m41^2,    1e-1 -> 20  ---> "log10()" ---> -1 -> 1.30103  
  
  const int bins_theta = 40;
  const int bins_dm2   = 40;
  h2_space = new TH2D("h2_space", "h2_space", bins_theta, -2, 0, bins_dm2, -1, 1.30103);

  for(int xbin=0; xbin<=bins_theta; xbin++) {
    cout<<TString::Format(" ---> initializing xbin %3d / %3d", xbin, bins_theta)<<endl;
    
    for(int ybin=0; ybin<=bins_dm2; ybin++) {      
      if( ybin==0 && xbin!=0 ) continue;
      if( xbin==0 && ybin!=0 ) continue;
      
      matrix_syst_frac_flux_before[xbin][ybin].ResizeTo(bins_oldworld, bins_oldworld);
      matrix_syst_frac_flux_before[xbin][ybin].Clear();
      matrix_syst_frac_flux_before[xbin][ybin].ResizeTo(bins_oldworld, bins_oldworld);
           
      matrix_syst_frac_geant_before[xbin][ybin].ResizeTo(bins_oldworld, bins_oldworld);
      matrix_syst_frac_geant_before[xbin][ybin].Clear();
      matrix_syst_frac_geant_before[xbin][ybin].ResizeTo(bins_oldworld, bins_oldworld);
          
      matrix_syst_frac_Xs_before[xbin][ybin].ResizeTo(bins_oldworld, bins_oldworld);
      matrix_syst_frac_Xs_before[xbin][ybin].Clear();
      matrix_syst_frac_Xs_before[xbin][ybin].ResizeTo(bins_oldworld, bins_oldworld);
          
      matrix_syst_frac_detector_before[xbin][ybin].ResizeTo(bins_oldworld, bins_oldworld);
      matrix_syst_frac_detector_before[xbin][ybin].Clear();
      matrix_syst_frac_detector_before[xbin][ybin].ResizeTo(bins_oldworld, bins_oldworld);
           
      matrix_syst_frac_additional_before[xbin][ybin].ResizeTo(bins_oldworld, bins_oldworld);
      matrix_syst_frac_additional_before[xbin][ybin].Clear();
      matrix_syst_frac_additional_before[xbin][ybin].ResizeTo(bins_oldworld, bins_oldworld);

      /////////////////////////////////////////// flux_geant_Xs
      
      for( int idx=1; idx<=17; idx++ ) {
	
	roostr = flux_geant_Xs_file_dir+TString::Format("result_syst_%d_%d/XsFlux/cov_%d.root", xbin, ybin, idx);
	
	TFile *file_temp = new TFile(roostr, "read");	
	TMatrixD *matrix_temp = (TMatrixD*)file_temp->Get( TString::Format("frac_cov_xf_mat_%d",idx) );
	if(idx<=13) matrix_syst_frac_flux_before[xbin][ybin] += (*matrix_temp);
	else if(idx<=16) matrix_syst_frac_geant_before[xbin][ybin] += (*matrix_temp);
	else matrix_syst_frac_Xs_before[xbin][ybin] += (*matrix_temp);
	delete matrix_temp;
	delete file_temp;	
      }

      /////////////////////////////////////////// MCstat
      
      roostr = flux_geant_Xs_file_dir+TString::Format("result_syst_%d_%d/mc_stat/0.log", xbin, ybin);
      
      int count_check_MCstat = 0;
      string line_check_MCstat;    
      ifstream file_check(roostr);
      while( getline(file_check, line_check_MCstat) ) count_check_MCstat++;      
      if( (count_check_MCstat-1)!=bins_newworld ) {
	cerr<<" Error (count_check_MCstat-1)!=bins_newworld: "<<roostr<<endl; exit(1);
      }
        
      matrix_syst_frac_MCstat_after[xbin][ybin].ResizeTo(bins_newworld, bins_newworld);
      matrix_syst_frac_MCstat_after[xbin][ybin].Clear();
      matrix_syst_frac_MCstat_after[xbin][ybin].ResizeTo(bins_newworld, bins_newworld);

      ifstream cppfile_temp(roostr, ios::in);
      int line_cppfile = -1;
      for(int idx=1; idx<=bins_newworld+1; idx++) {
	double Lee(0), run(0);
	int gbin = 0; int lbin = 0; double val_pred = 0; double mc_stat = 0; double nn_stat = 0;
	if(idx==1) { cppfile_temp>>Lee>>run; }
	else {
	  line_cppfile++;
	  cppfile_temp>>gbin>>lbin>>val_pred>>mc_stat>>nn_stat;
	  matrix_syst_frac_MCstat_after[xbin][ybin](line_cppfile, line_cppfile) = mc_stat;
	}
      }// idx
      cppfile_temp.close();
 
    }// ybin
  }// xbin
    
  //////////////////////////////////////
  
  matrix_dataFIT_newworld.Clear();
  matrix_dataFIT_newworld.ResizeTo(1, bins_newworld);
  
  //////////////////////////////////////
  
  cout<<endl;
  cout<<" ---------------------------------------"<<endl;
  cout<<" Make sure the \"zeroout_ch\" is correct"<<endl;
  for(auto it_int_map=zeroout_ch.begin(); it_int_map!=zeroout_ch.end(); it_int_map++) {
    cout<<" zeroout channel: "<<it_int_map->first<<endl;
  }
  
  cout<<endl<<" ---> Finish Set_Spectra_MatrixCov"<<endl<<endl;
  
}
