//#include<iostream>
//#include<fstream>

void read_scale()
{
  TString roostr = "";

  ///////////////////////////////////////////////////////////

  TFile *roofile = new TFile("za_nueOsc_afterScale_merge.root", "read");
  
  TH1F *h1_nue_ff_FC = (TH1F*)roofile->Get("histo_1");
  TH1F *h1_nue_ff_PC = (TH1F*)roofile->Get("histo_2");

  ///////////////////////////////////////////////////////////

  roostr = "h1_nue_FC_mm"; TH1D *h1_nue_FC_mm = new TH1D(roostr, roostr, 30, 0, 3000);
  roostr = "h1_nue_PC_mm"; TH1D *h1_nue_PC_mm = new TH1D(roostr, roostr, 30, 0, 3000);  
  
  TString file_eventlist = "za_nueOsc_eventlist.txt";
  ifstream inputfile_aa(file_eventlist, ios::in);
  if(!inputfile_aa) {
    cerr<<" No inputfile of nueOsc before scaling"<<endl;
    exit(1);
  }
  
  string str_inputfile_aa = "";
  int line_inputfile_aa = 0;
  while( std::getline(inputfile_aa, str_inputfile_aa) ) line_inputfile_aa++;
  cout<<endl<<" ---> line_inputfile nueOsc: "<<line_inputfile_aa<<endl<<endl;

  ifstream inputfile_ab(file_eventlist, ios::in);
  for(int idx=1; idx<=line_inputfile_aa; idx++) {
    int pdg_val(0);
    double weight_val(0), nueEtrue(0), nueEreco(0), baseline(0);
    TString ch_name = "";
    inputfile_ab>>pdg_val>>weight_val>>nueEtrue>>nueEreco>>baseline>>ch_name;

    int flag_FC = 0;
    if( ch_name.Contains("_FC_") ) flag_FC = 1;
    //cout<<idx<<"\t"<<flag_FC<<"\t"<<ch_name<<endl;

    if( flag_FC ) h1_nue_FC_mm->Fill( nueEreco, weight_val );
    else h1_nue_PC_mm->Fill( nueEreco, weight_val );
    
  }

  double ff_num_FC = h1_nue_ff_FC->Integral(1, 30);
  double ff_num_PC = h1_nue_ff_PC->Integral(1, 30);

  double mm_num_FC = h1_nue_FC_mm->Integral(1, 30);
  double mm_num_PC = h1_nue_PC_mm->Integral(1, 30);

  double scale_FC = ff_num_FC/mm_num_FC;
  double scale_PC = ff_num_PC/mm_num_PC;
 
  h1_nue_FC_mm->Scale( scale_FC );
  h1_nue_PC_mm->Scale( scale_PC );

  ///////////////////////////////
  ///////////////////////////////

  roostr = "canv_nue";
  TCanvas *canv_nue = new TCanvas(roostr, roostr, 900, 650);
  
  h1_nue_ff_FC->Draw();
  h1_nue_ff_FC->SetStats(0);
  h1_nue_ff_FC->SetTitle("");
  h1_nue_ff_FC->SetLineColor(kBlue);
  h1_nue_ff_FC->SetMarkerColor(kBlue);
  h1_nue_ff_FC->SetMarkerStyle(20);

  h1_nue_ff_PC->Draw("same");
  h1_nue_ff_PC->SetLineColor(kRed);
  h1_nue_ff_PC->SetMarkerColor(kRed);
  h1_nue_ff_PC->SetMarkerStyle(20);
  
  h1_nue_FC_mm->Draw("same hist");
  h1_nue_FC_mm->SetLineColor(kBlue);

  h1_nue_PC_mm->Draw("same hist");
  h1_nue_PC_mm->SetLineColor(kRed);

  
  
}
