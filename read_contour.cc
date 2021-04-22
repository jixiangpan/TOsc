void read_contour()
{
  gROOT->ProcessLine(".x ./lhcbStyle_edit.C");
  TString roostr = "";
  
  // Y: m41^2, eV^2
  // X: sin22t14
  
  TString roostr_4vfit3v = "log_results_4vfit3v";
  TString roostr_3vfit4v = "log_results_3vfit4v";
  
  const int max_Num = 500000;
  
  ///////
  const int Ncontour = 3;
  double contours[Ncontour] = {0};
  contours[0] = 0.9;
  contours[1] = 0.95;
  contours[2] = 0.9973;

  int colors[Ncontour] = {kBlue, kGreen+3, kRed};
  
  ///////
  double val_t14[max_Num] = {0};
  double val_m41[max_Num] = {0};

  int bin_t14[max_Num] = {0};
  int bin_m41[max_Num] = {0};
  
  double val_chi2_4vfit3v[max_Num] = {0};
  double val_chi2_3vfit4v[max_Num] = {0};
  
  ///////
  ///////
  for(int index=0; index<max_Num; index++) {
    val_t14[index] = -1;
    val_m41[index] = -1;
    bin_t14[index] = -1;
    bin_m41[index] = -1;

    val_chi2_4vfit3v[index] = -1;
    val_chi2_3vfit4v[index] = -1;
  }

  ///////
  //TH2D *h2_scan = new TH2D("h2_scan", "h2_scan", 60, -2, 0, 60, -1, 1.4);
  TH2D *h2_scan = new TH2D("h2_scan", "h2_scan", 80, -2, 0, 80, -1, 1.4);
  int xbins_scan = h2_scan->GetNbinsX();
  int ybins_scan = h2_scan->GetNbinsY();

  int line_scan = 0;
  for(int xbin=1; xbin<=xbins_scan; xbin++) {
    for(int ybin=1; ybin<=ybins_scan; ybin++) {
      line_scan++;
      double center_x = h2_scan->GetXaxis()->GetBinCenter(xbin);
      double center_y = h2_scan->GetYaxis()->GetBinCenter(ybin);
      
      double val_s_2_2_14 = pow(10, center_x);
      double val_delta_m41_2 = pow(10, center_y);
      
      int index = line_scan;
      val_t14[index] = val_s_2_2_14;
      val_m41[index] = val_delta_m41_2;
      
      bin_t14[index] = xbin;
      bin_m41[index] = ybin;
    }
  }
  
  ///////
  int max_index = 0;

  ///////
  ifstream read_4vfit3v(roostr_4vfit3v, ios::in);
  if( !read_4vfit3v ) {cerr<<" ---> read_4vfit3v doesn't exist!"<<endl; exit(0);}
  while( !read_4vfit3v.eof() ) {
    int index = 0;
    double chi2 = -1;
    double t14 = -1;
    double m41 = -1;
    read_4vfit3v>>index>>chi2>>t14>>m41;

    if( index==0 ) continue;
    val_chi2_4vfit3v[index] = chi2;
    if( index > max_index ) max_index = index;
  }

  ///////
  ifstream read_3vfit4v(roostr_3vfit4v, ios::in);
  if( !read_3vfit4v ) {cerr<<" ---> read_3vfit4v doesn't exist!"<<endl; exit(0);}
  while( !read_3vfit4v.eof() ) {
    int index = 0;
    double chi2 = -1;
    double t14 = -1;
    double m41 = -1;
    read_3vfit4v>>index>>chi2>>t14>>m41;

    if( index==0 ) continue;
    val_chi2_3vfit4v[index] = chi2;
    if( index > max_index ) max_index = index;
  }
 
  ///////
  /////////////////////////////////////////////////////////////////////////////////////////////////////// mc sensitivity
  /////////////////////////////////////////////////////////////////////////////////////////////////////// mc sensitivity
  
  TH2D *h2_CLs_mc = (TH2D*)h2_scan->Clone("h2_CLs_mc");
  TH2D *h2_counter_mc = (TH2D*)h2_scan->Clone("h2_counter_mc");

  for(int index=1; index<max_index; index++) {
    if( val_chi2_3vfit4v[index]>=0 && val_chi2_4vfit3v[index]>=0) {
 
      double delta_4v = -val_chi2_3vfit4v[index];
      double delta_3v = val_chi2_4vfit3v[index];
      double delta_dd = delta_3v;
      
      double numerator_erf = TMath::Erf( (delta_4v-delta_dd)/sqrt(8.*fabs(delta_4v)) );
      double denominator_erf = TMath::Erf( (delta_3v-delta_dd)/sqrt(8.*fabs(delta_3v)) );

      double CLs = (1+numerator_erf)/(1+denominator_erf);
	
      h2_CLs_mc->SetBinContent( bin_t14[index], bin_m41[index], CLs );
      h2_counter_mc->SetBinContent( bin_t14[index], bin_m41[index], 1-CLs );
    }
  }

  roostr = "canv_h2_CLs_mc";
  TCanvas *canv_h2_CLs_mc = new TCanvas(roostr, roostr, 800, 600);
  h2_CLs_mc->Draw("colz");

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
    return;
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

  ///////
  roostr = "h2_basic_CLs_mc";
  TH2D *h2_basic_CLs_mc = new TH2D(roostr, roostr, 100, 4.5e-3, 1, 100, 4.e-2, 25);

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
  
  int ic = 1;
  gh_CLs_mc[ic]->Draw("same l");
  gh_CLs_mc[ic]->SetLineColor( colors[ic] );

  //gh_CLs_mc[ic]->SaveAs("gh_CLs_mc_95_bin2bin10.root");

  TFile *gh_CLs_mc_95_bin2bin10 = new TFile("gh_CLs_mc_95_bin2bin10.root", "read");
  TGraph *gh_bin2bin10 = (TGraph*)gh_CLs_mc_95_bin2bin10->Get("gh_CLs_mc_1");
  gh_bin2bin10->Draw("same l");

  TFile *gh_CLs_mc_95_bin2binFree = new TFile("gh_CLs_mc_95_bin2binFree.root", "read");
  TGraph *gh_bin2binFree = (TGraph*)gh_CLs_mc_95_bin2binFree->Get("gh_CLs_mc_1");
  gh_bin2binFree->Draw("same l");
  
  gh_CLs_mc[1]->SetLineColor( colors[0] );
  gh_bin2bin10->SetLineColor( colors[1] );
  gh_bin2binFree->SetLineColor( colors[2] );
  

  // for(int ic=0; ic<Ncontour; ic++) {
  //   gh_CLs_mc[ic]->Draw("same l");
  //   gh_CLs_mc[ic]->SetLineColor( colors[ic] );
  // }

  TLegend *lg_mc = new TLegend(0.18,0.22,0.48,0.45);
  lg_mc->SetFillStyle(0);
  lg_mc->AddEntry("", "95% CL_{s} with 33 days", "");
  lg_mc->AddEntry(gh_CLs_mc[1], "Bin2Bin error 2%", "l");
  lg_mc->AddEntry(gh_bin2bin10, "Bin2Bin error 10%", "l");
  lg_mc->AddEntry(gh_bin2binFree, "Bin2Bin free", "l");
  lg_mc->Draw();
  lg_mc->SetTextSize(0.045);
  lg_mc->SetTextFont(132);
  lg_mc->Draw();

  canv_h2_basic_CLs_mc->SaveAs("canv_h2_basic_CLs_mc.png");
  canv_h2_basic_CLs_mc->SaveAs("canv_h2_basic_CLs_mc.pdf");
 
}
