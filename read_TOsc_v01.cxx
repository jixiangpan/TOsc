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

  int xlow = 1;
  int xhgh = 1;
  
  for(int i=1; i<argc; i++) {
    
    if( strcmp(argv[i],"-p")==0 ) {
      stringstream convert( argv[i+1] );
      if(  !( convert>>scaleF_POT ) ) { cerr<<" ---> Error scaleF_POT !"<<endl; exit(1); }
    }

    if( strcmp(argv[i],"-f")==0 ) {
      stringstream convert( argv[i+1] );
      if(  !( convert>>ifile ) ) { cerr<<" ---> Error ifile !"<<endl; exit(1); }
    }

    ///
    if( strcmp(argv[i],"-x")==0 ) {
      stringstream convert( argv[i+1] );
      if(  !( convert>>space_xbin ) ) { cerr<<" ---> Error space_xbin !"<<endl; exit(1); }
    }

    if( strcmp(argv[i],"-y")==0 ) {
      stringstream convert( argv[i+1] );
      if(  !( convert>>space_ybin ) ) { cerr<<" ---> Error space_ybin !"<<endl; exit(1); }
    }

    ///
    if( strcmp(argv[i],"-xl")==0 ) {
      stringstream convert( argv[i+1] );
      if(  !( convert>>xlow ) ) { cerr<<" ---> Error xlow !"<<endl; exit(1); }
    }

    if( strcmp(argv[i],"-xh")==0 ) {
      stringstream convert( argv[i+1] );
      if(  !( convert>>xhgh ) ) { cerr<<" ---> Error xhgh !"<<endl; exit(1); }
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
				  config_Osc::syst_result_dir);
  
  Osc_test->scaleF_POT = scaleF_POT;  
  Osc_test->Apply_POT_scaled();

  ////////// can do any times
  
  Osc_test->flag_syst_flux       = config_Osc::flag_syst_flux;
  Osc_test->flag_syst_geant      = config_Osc::flag_syst_geant;
  Osc_test->flag_syst_Xs         = config_Osc::flag_syst_Xs;
  Osc_test->flag_syst_detector   = config_Osc::flag_syst_detector;
  Osc_test->flag_syst_additional = config_Osc::flag_syst_additional;
  Osc_test->flag_syst_MCstat     = config_Osc::flag_syst_MCstat;

  Osc_test->flag_FIT_after_constraint = config_Osc::flag_FIT_after_constraint;
  
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
  // Osc_test->Minimization_OscPars_FullCov(0.1, 1.8, 1, 0);

  if( 0 ) {// best-fit value
    Osc_test->Set_OscPars(0.1, 1);
    Osc_test->Set_Collapse();
    Osc_test->Set_data2dataFIT();
    //Osc_test->Set_Asimov2dataFIT();
    Osc_test->Minimization_OscPars_FullCov(0.001, 0.8, 0, 0);
  }

  if( 0 ) {
    Osc_test->Set_OscPars(0, 1);
    Osc_test->Set_Collapse();
    Osc_test->Produce_Variations(1);
    Osc_test->Set_Toy2dataFIT(1);

    int flag_fit_good       = 0;
    double out_fit_chi2     = 1e6;
    double out_fit_s22theta = 0;
    double out_fit_dm2      = 0;        
        
    for(int ybin=1; ybin<=Osc_test->h2_space->GetNbinsY(); ybin++) {

      double val_dm2 = pow( 10, Osc_test->h2_space->GetYaxis()->GetBinCenter(ybin) );

      cout<<endl;
      
      Osc_test->Minimization_OscPars_FullCov(0, val_dm2, 1, 1);
      double fix_fit_chi2 = Osc_test->minimization_chi2;
      cout<<TString::Format(" ---> %2d fitting results: status %4d, min_chi2 %8.3f, s22t %8.6f, dm2 %9.6f",
			    ybin,
			    Osc_test->minimization_status,
			    Osc_test->minimization_chi2,
			    Osc_test->minimization_s22theta_val,
			    Osc_test->minimization_dm2_val
			    )<<endl;

      Osc_test->Minimization_OscPars_FullCov(0.1, val_dm2, 0, 1);
      int free_fit_status  = Osc_test->minimization_status;
      double free_fit_chi2 = Osc_test->minimization_chi2;
      cout<<TString::Format(" ---> %2d fitting results: status %4d, min_chi2 %8.3f, s22t %8.6f, dm2 %9.6f",
			    ybin,
			    Osc_test->minimization_status,
			    Osc_test->minimization_chi2,
			    Osc_test->minimization_s22theta_val,
			    Osc_test->minimization_dm2_val
			    )<<endl;

      if( free_fit_status==0 ) {

	if( free_fit_chi2>fix_fit_chi2 ) {
	  cout<<" ooo> free_fit_chi2 > fix_fit_chi2 "<<endl;
	  
	  Osc_test->Minimization_OscPars_FullCov(0.0001, val_dm2, 0, 1);
	  cout<<TString::Format(" oaa> %2d fitting results: status %4d, min_chi2 %8.3f, s22t %8.6f, dm2 %9.6f",
				ybin,
				Osc_test->minimization_status,
				Osc_test->minimization_chi2,
				Osc_test->minimization_s22theta_val,
				Osc_test->minimization_dm2_val
				)<<endl;

	  if( Osc_test->minimization_status==0 ) {
	    Osc_test->Minimization_OscPars_FullCov(Osc_test->minimization_s22theta_val, val_dm2, 0, 0);
	    cout<<TString::Format(" obb> %2d fitting results: status %4d, min_chi2 %8.3f, s22t %8.6f, dm2 %9.6f",
				  ybin,
				  Osc_test->minimization_status,
				  Osc_test->minimization_chi2,
				  Osc_test->minimization_s22theta_val,
				  Osc_test->minimization_dm2_val
				  )<<endl;

	    if( Osc_test->minimization_status==0 ) {
	      flag_fit_good = 1;
	      if( out_fit_chi2>Osc_test->minimization_chi2 ) {
		out_fit_chi2 = Osc_test->minimization_chi2;
		out_fit_s22theta = Osc_test->minimization_s22theta_val;
		out_fit_dm2 = Osc_test->minimization_dm2_val;
	      }
	    }
	    
	  }
	}// if( free_fit_chi2>fix_fit_chi2 )
	else {
	  cout<<" ppp> free_fit_chi2 < fix_fit_chi2 "<<endl;
	  
	  Osc_test->Minimization_OscPars_FullCov(Osc_test->minimization_s22theta_val, val_dm2, 0, 0);
	  cout<<TString::Format(" paa> %2d fitting results: status %4d, min_chi2 %8.3f, s22t %8.6f, dm2 %9.6f",
				ybin,
				Osc_test->minimization_status,
				Osc_test->minimization_chi2,
				Osc_test->minimization_s22theta_val,
				Osc_test->minimization_dm2_val
				)<<endl;

	  if( Osc_test->minimization_status==0 ) {
	    flag_fit_good = 1;
	    if( out_fit_chi2>Osc_test->minimization_chi2 ) {
	      out_fit_chi2 = Osc_test->minimization_chi2;
	      out_fit_s22theta = Osc_test->minimization_s22theta_val;
	      out_fit_dm2 = Osc_test->minimization_dm2_val;
	    }
	  }
	  
	}// else
	
      }// if( free_fit_status==0 )
      
      cout<<endl;
    }// for(int ybin=1; ybin<=Osc_test->h2_space->GetNbinsY(); ybin++)

    
    roostr = TString::Format("out_result_%05d.txt", ifile);
    ofstream ListWrite(roostr, ios::out|ios::trunc);
    ListWrite<<TString::Format("%d %14.6f %10.6f %10.6f", flag_fit_good, out_fit_chi2, out_fit_s22theta, out_fit_dm2);
    
    for(int idx=0; idx<(Osc_test->matrix_dataFIT.GetNcols()); idx++) {
      ListWrite<<TString::Format(" %5.1f", Osc_test->matrix_dataFIT(0, idx));
    }
    ListWrite<<endl;
      
    ListWrite.close();
    
  }
  
  //////////////////////////////////////////////////////// for one point
  
  if( 0 ) {

    cout<<endl<<" ---> test Oscillation parameters"<<endl<<endl;
    
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
    Osc_test->Minimization_OscPars_FullCov(0, 7.25, 1, 0);
    dchi2_4vAsimov = (-Osc_test->minimization_chi2);
    
    /////////////////////////////////// 3v Asimov    
    Osc_test->Set_OscPars(0, 1);
    Osc_test->Set_Collapse();
    Osc_test->Set_Asimov2dataFIT();
    Osc_test->Minimization_OscPars_FullCov(test_s22theta, test_dm2, 1, 0);
    dchi2_3vAsimov = Osc_test->minimization_chi2;

    /////////////////////////////////// data  
    Osc_test->Set_data2dataFIT();
  
    Osc_test->Minimization_OscPars_FullCov(0, 1, 1, 0);
    chi2_3v_on_data = Osc_test->minimization_chi2;
  
    Osc_test->Minimization_OscPars_FullCov(test_s22theta, test_dm2, 1, 0);
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

    double d4v = delta_4v;
    double d3v = delta_3v;          
    double CL_sens = 1 + TMath::Erf( (d4v-d3v)/sqrt(8*fabs(d4v)) );
    cout<<" sens "<< 1-CL_sens<<endl;

    cout<<endl<<d4v<<"\t"<<d3v<<"\t"<<delta_dd<<endl<<endl;
    cout<<" 4v/3v: "<<chi2_4v_on_data<<"\t"<<chi2_3v_on_data<<endl;

    
    // Osc_test->Set_OscPars(0, 1);
    // Osc_test->Set_Collapse();
    // for(int idx=0; idx<(Osc_test->bins_newworld); idx++) {
    //   double val_data = Osc_test->matrix_data_newworld(0, idx);
    //   double val_pred = Osc_test->matrix_pred_newworld(0, idx);
    //   cout<<TString::Format(" ---> %3d  %9.3f  %9.3f", idx+1, val_data, val_pred)<<endl;
    // }
    
  }

  //////////////////////////////////////////////////////// for jobs

  if( 0 ) {
    int bins_theta = 40;
    int bins_dm2   = 40;

    ////// X: sin22t14, 1e-2 -> 1   ---> "log10()" ---> -2 -> 0
    ////// Y: m41^2,    1e-1 -> 20  ---> "log10()" ---> -1 -> 1.30103      
    TH2D *h2_space = new TH2D("h2_space_jobs", "h2_space_jobs", bins_theta, -2, 0, bins_dm2, -1, 1.30103);
        
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
    Osc_test->Minimization_OscPars_FullCov(0, 7.25, 1, 0);
    dchi2_4vAsimov = (-Osc_test->minimization_chi2);
    
    /////////////////////////////////// 3v Asimov    
    Osc_test->Set_OscPars(0, 1);
    Osc_test->Set_Collapse();
    Osc_test->Set_Asimov2dataFIT();
    Osc_test->Minimization_OscPars_FullCov(test_s22theta, test_dm2, 1, 0);
    dchi2_3vAsimov = Osc_test->minimization_chi2;

    /////////////////////////////////// data  
    Osc_test->Set_data2dataFIT();
    
    Osc_test->Minimization_OscPars_FullCov(0, 1, 1, 0);
    chi2_3v_on_data = Osc_test->minimization_chi2;
  
    Osc_test->Minimization_OscPars_FullCov(test_s22theta, test_dm2, 1, 0);
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

    roostr = TString::Format("out_result_%04d_%04d.txt", space_xbin, space_ybin);
    ofstream ListWrite(roostr, ios::out|ios::trunc);
    ListWrite<<TString::Format("%4d %4d %20.10f %20.10f %20.10f %18.15f",
			       space_xbin, space_ybin,
			       delta_4v, delta_3v, delta_dd, CL
			       )<<endl;
    ListWrite.close();
  }
  
  //////////////////////////////////////////////////////// for whole scan at one time
  
  if( 0 ) {

    cout<<endl<<" ---> for whole scan at one time"<<endl<<endl;
    
    int bins_theta = 40;
    int bins_dm2   = 40;

    ////// X: sin22t14, 1e-2 -> 1   ---> "log10()" ---> -2 -> 0
    ////// Y: m41^2,    1e-1 -> 20  ---> "log10()" ---> -1 -> 1.30103      
    TH2D *h2_space = new TH2D("h2_space_whole", "h2_space_whole", bins_theta, -2, 0, bins_dm2, -1, 1.30103);
    
    for(int ibin=xlow; ibin<=xhgh; ibin++) {      
      cout<<TString::Format(" ---> processing %4d/%4d", ibin, bins_theta)<<endl;
      
      for(int jbin=1; jbin<=bins_dm2; jbin++) {
	//cout<<TString::Format(" ---> processing sub %4d - %4d", ibin, jbin)<<endl;
	
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
	Osc_test->Minimization_OscPars_FullCov(0, 7.25, 1, 0);
	dchi2_4vAsimov = (-Osc_test->minimization_chi2);
	double chi2_3v_on_4vAsimov = Osc_test->minimization_chi2;
	
	/////////////////////////////////// 3v Asimov  
	Osc_test->Set_OscPars(0, 1);
	Osc_test->Set_Collapse();
	Osc_test->Set_Asimov2dataFIT();
	Osc_test->Minimization_OscPars_FullCov(test_s22theta, test_dm2, 1, 0);
	dchi2_3vAsimov = Osc_test->minimization_chi2;
	double chi2_4v_on_3vAsimov = Osc_test->minimization_chi2;
	
	/////////////////////////////////// data
  
	Osc_test->Set_data2dataFIT();
  
	Osc_test->Minimization_OscPars_FullCov(0, 1, 1, 0);
	chi2_3v_on_data = Osc_test->minimization_chi2;
  
	Osc_test->Minimization_OscPars_FullCov(test_s22theta, test_dm2, 1, 0);
	chi2_4v_on_data = Osc_test->minimization_chi2;
  
	dchi2_data = chi2_4v_on_data - chi2_3v_on_data;

	///////
	///////

	double delta_4v = dchi2_4vAsimov;
	double delta_3v = dchi2_3vAsimov;
	double delta_dd = dchi2_data;
  
	double numerator_erf = TMath::Erf( (delta_4v-delta_dd)/sqrt(8.*fabs(delta_4v)) );
	double denominator_erf = TMath::Erf( (delta_3v-delta_dd)/sqrt(8.*fabs(delta_3v)) );
  
	double CLs = (1+numerator_erf)/(1+denominator_erf);
	double CL = 1 - CLs;  
	
	roostr = TString::Format("out_result_%04d_%04d.txt", ibin, jbin);
	ofstream ListWrite(roostr, ios::out|ios::trunc);
	ListWrite<<TString::Format("%4d %4d %20.10f %20.10f %20.10f %18.15f",
				   ibin, jbin,
				   delta_4v, delta_3v, delta_dd, CL
				   )<<endl;
	ListWrite.close();

	///////
	///////
	
	roostr = TString::Format("qfo_result_%04d_%04d.txt", ibin, jbin);
	ofstream ListWrite_qfo(roostr, ios::out|ios::trunc);
	ListWrite_qfo<<TString::Format("%4d %4d %20.12f %20.12f %20.12f %20.12f %10.6f %10.6ff",
				       ibin, jbin,
				       chi2_4v_on_data,
				       chi2_3v_on_data,				       
				       chi2_4v_on_3vAsimov,
				       chi2_3v_on_4vAsimov,
				       test_s22theta, test_dm2
				       )<<endl;
	ListWrite.close();
	
      }// jbin
    }// ibin
          
  }

  ///////////////////////////////////////////////////////////////

  cout<<endl;
  cout<<" flag_syst_flux       "<<Osc_test->flag_syst_flux<<endl;
  cout<<" flag_syst_geant      "<<Osc_test->flag_syst_geant<<endl;
  cout<<" flag_syst_Xs         "<<Osc_test->flag_syst_Xs<<endl;
  cout<<" flag_syst_detector   "<<Osc_test->flag_syst_detector<<endl;
  cout<<" flag_syst_additional "<<Osc_test->flag_syst_additional<<endl;
  cout<<" flag_syst_MCstat     "<<Osc_test->flag_syst_MCstat<<endl;
  cout<<endl;
  
  ///////////////////////////////////////////////////////////////
  
  if( config_Osc::flag_display_graphics ) {
    cout<<endl<<" Enter Ctrl+c to end the program"<<endl;
    cout<<" Enter Ctrl+c to end the program"<<endl<<endl;  
    theApp.Run();
  }
    
  return 0;
}
