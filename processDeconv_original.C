TFile *fout;
vector<vector<TH1D*>> hAsym_e1,hAsym_eP1;
int separateW=0;

const int nSp=5;
const string spTit[nSp]={"e/#pi","e/#pi E>1","#gamma","neutron","primary e E>1"};
const string spH[nSp]={"e","e1","g","n","eP1"};
map<int,int> spM {{11,1},{211,1},{22,3},{2112,4}};

const int nDet=1;
const string detH[nDet]={"det28"};
map<int,int> dtM {{28,1}};

const int nSec=3;
const string secH[nSec]={"closed","trans","open"};

string fin;
int nFiles(0);
long currentEvNr(0);

const double pi = acos(-1);

void initHisto(int);
long processOne(string);
void process();
int findDetector(int &sector, double phi, double r);
void writeOutput();

const int nRings = 7;
double rMin[nRings][nSec],rMax[nRings][nSec];
int setSegmentation();

void processDeconv(const string& finName = "./remollout.root", int wRegions=0,int tiling=1){
  if(!setSegmentation())
    return;
  separateW = wRegions;
  fin = finName;
  initHisto(tiling);
  process();
  writeOutput();
}

void process(){
  if(fin==""){
    cout<<"\t did not find input file. Quitting!"<<endl;
    return 2;
  }

  long nTotEv(0);
  if( fin.find(".root") < fin.size() ){
    cout<<"Processing single file:\n\t"<<fin<<endl;
    nTotEv+=processOne(fin);
    nFiles=1;
  }else{
    cout<<"Attempting to process list of output from\n\t"<<fin<<endl;
    ifstream ifile(fin.c_str());
    string data;
    while(ifile>>data){
      cout<<" processing: "<<data<<endl;
      nTotEv+=processOne(data);
      nFiles++;
    }
  }
  cout<<"\nFinished processing a total of "<<nTotEv<<endl;
}

void initHisto(int tileConf){
  string foutNm = Form("%s_tileConf%d_procDeconvV1.root",fin.substr(0,fin.find(".")).c_str(),tileConf);
  fout = new TFile(foutNm.c_str(),"RECREATE");
  fout->mkdir("deconvolution","histos for deconvolution analysis");
  fout->cd("deconvolution");

  const string secNm[3]={"closed","transition","open"};
  for(int i=0;i<6;i++){
    vector<TH1D*> dt1,dt2;
    for(int j=0;j<3;j++){
      dt1.push_back(new TH1D(Form("hAsym_e1_R%d_S%d",i+1,j),
			     Form("rate weighted Asyms for Ring %d Sector %s;asymmetry [ppb]",i+1,secNm[j].c_str()),
			     100,-1000000,1000000));
      dt2.push_back(new TH1D(Form("hAsym_eP1_R%d_S%d",i+1,j),
			     Form("rate weighted Asyms for Ring %d Sector %s;asymmetry [ppb]",i+1,secNm[j].c_str()),
			     100,-1000000,1000000));
    }
    hAsym_e1.push_back(dt1);
    hAsym_eP1.push_back(dt2);
  } 
}

long processOne(string fnm){
  TFile *fin=TFile::Open(fnm.c_str(),"READ");
  TTree *t=(TTree*)fin->Get("T");
  if (t == 0) return 0;
  Double_t rate=0;
  std::vector<remollGenericDetectorHit_t> *hit=0;
  t->SetBranchAddress("rate", &rate);
  t->SetBranchAddress("hit", &hit);

  long nEntries = t->GetEntries();
  cout<<"\tTotal events: "<<nEntries<<endl;
  float currentProc=1,procStep=60;
  vector<int> procID;
  int sector(-1);
  double pi = acos(-1);

  //for (Long64_t event = 0; event < 5; t->GetEntry(event++)) {
  for (Long64_t event = 0; event < nEntries; t->GetEntry(event++)) {
    currentEvNr++;
    if( float(event+1)/nEntries*100 > currentProc){
      cout<<"at tree entry\t"<<event<<"\t"<< float(event+1)/nEntries*100<<endl;
      currentProc+=procStep;
    }
  }
  fin->Close();
  delete fin;
  return nEntries;
};

void writeOutput(){
  fout->cd("deconvolution");
  for(int i=0;i<6;i++)
    for(int j=0;j<3;j++){
      hAsym_e1[i][j]->Scale(1./nFiles);
      hAsym_e1[i][j]->Write();
      hAsym_eP1[i][j]->Scale(1./nFiles);
      hAsym_eP1[i][j]->Write();
    }
  fout->Close();
}

int findDetector(int &sector, double phi, double r){

  const double secPhi = fmod(phi, 2*pi/7);

  //0,1,2 == closed, transition, open
  if( secPhi < pi/28 )
    sector = 0;
  else if( secPhi < 3*pi/28 )
    sector = 1;
  else if( secPhi < 5*pi/28 )
    sector = 2;
  else if( secPhi < 7*pi/28 )
    sector = 1;
  else if( secPhi < 8*pi/28 )
    sector = 0;

  // this will just pick up the first hit and ignore the rest if there is overlap in tiling
  for(int i=0;i<nRings;i++)
    if(r >= rMin[i][sector] && r <= rMax[i][sector])
      return i;

  return -1;
}


int setSegmentation(){

    double rm3[nRings][nSec]={
      { 650.0,  650.0,  650.0},
      { 690.0,  690.0,  690.0},
      { 735.0,  735.0,  735.0},
      { 790.0,  790.0,  790.0},
      { 900.0,  900.0,  900.0},
      {1060.0, 1060.0, 1060.0},
      {1200.0, 1200.0, 1200.0}
    };
    double rM3[nRings][nSec]={
      { 690.0,  690.0,  690.0},
      { 735.0,  735.0,  735.0},
      { 790.0,  790.0,  790.0},
      { 900.0,  900.0,  900.0},
      {1060.0, 1060.0, 1060.0},
      {1160.0, 1160.0, 1160.0},
      {1500.0, 1500.0, 1500.0}
    };
    cout<<"using tiling configuration from 2021-05-01 for parallel segmented magnets; uniform R5=16cm"<<endl;
    for(int i=0;i<nRings;i++)
      for(int j=0;j<nSec;j++){
	rMin[i][j]=rm3[i][j];
	rMax[i][j]=rM3[i][j];
      }
    return 1;
}
