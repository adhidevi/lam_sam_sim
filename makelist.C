void makelist(TString maindir="/volatile/halla/parity/adhidevi/remoll_rootfiles/develop_br/", TString subdir="LH2_beam_V23", int Nfiles=5000, TString outfilename = "mylist.txt"){

  ofstream outfile;
  outfile.open(outfilename);

  int Nstart=1001;
  for(int ii=Nstart; ii<Nfiles+Nstart;ii++){
    outfile<<maindir.Data()<<subdir.Data()<<"/"<<subdir.Data()<<"_"<<ii<<".root"<<endl;
  }
}
