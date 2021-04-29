namespace config_Osc
{
  ////////// input files for spectra and covariance matrixes

  TString eventlist_dir = "./data_nue_beforeafter_scale/";
  TString event_summation_afterscale_file = "./data_nue_beforeafter_scale/zz_merge_summation_allruns.root";

  /// TOsc.cxx, manually set
  /// map<int, int>ch_nueCC_from_intrinsic_sample
  /// map<int, int>zeroout_ch;
  
  TString centralvalue_noosc_file = "./data_basic/merge_noosc.root";
  TString syst_result_dir = "./data_basic/";
  
  int channels_observation = 7;

  ////////// display graphics flag

  bool flag_display_graphics = 0;
  
  ////////// systematics flag
  
  bool flag_syst_flux       = 1;
  bool flag_syst_geant      = 1;
  bool flag_syst_Xs         = 1;
  bool flag_syst_detector   = 0;
  bool flag_syst_additional = 0;
  bool flag_syst_MCstat     = 1;

  bool flag_FIT_after_constraint = 1;
  
}
