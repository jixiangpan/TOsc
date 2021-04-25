void gen_list_farm()
{
  int from = 1001;
  int to   = 5000;
  int line = 0;
  
  TString str_exe = "./read_TOsc_v01 ";
  TString File = "";
  TString roostr = "";
  
  cout<<endl;
  cout<<endl;

  ofstream ListWrite("cmd_list.csh", ios::out|ios::trunc);
    
  //for(int idx=from; idx<=to; idx++)
  for(int idx=1; idx<=100; idx++) {
  for(int jdx=1; jdx<=100; jdx++)
    {
      line++;
     
      if( line%100==0 ) 
      cout<<TString::Format(" ---> processing %5d", line)<<endl;
      
      ofstream WriteFile(TString::Format("./cmd_%06d.csh",line),ios::out|ios::trunc);
      WriteFile<<"#!/bin/tcsh"<<endl;
      WriteFile<<"source /afs/ihep.ac.cn/soft/dayabay/NuWa-64/opt/external/ROOT/5.30.06_python2.7/x86_64-slc5-gcc41-opt/root/bin/thisroot.csh"<<endl;
      WriteFile<<"cd /dybfs2/users/jixp/work_ss/TOsc"<<endl;

      roostr = TString::Format(" -x %d -y %d ", idx, jdx);
      
      WriteFile<<str_exe<<roostr<<endl;
      WriteFile<<endl;
      WriteFile.close();

      ListWrite<<Form("echo %9d", line)<<endl;
      ListWrite<<"usleep 10"<<endl;// microseconds
      ListWrite<<TString::Format("hep_sub cmd_%06d.csh",line)<<endl;
    }
}

  ListWrite.close();

  cout<<endl;
  cout<<endl;


}

