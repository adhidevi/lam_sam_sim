TFile *fout;
vector<vector<TH1D*>> hR_e1,hR_eP1;
vector<vector<TH2D*>> hXY_e1,hXY_eP1;

const int nSp = 5;
const string spTit[nSp]={"e/#pi","e/#pi E>1","#gamma","neutron","primary e E>1"};
const string spH[nSp]={"e","e1","g","n","eP1"};
map<int,int> spM{{11,1},{211,1},{22,3},{2112,4}};

const int nDet=1;
const string detH[nDet]={"det28"};
map<int,int> dtM{{28,1}};

const int nSecDet = 21;//7 rings x 3 sectors

const int nSec = 3;
const string secH[nSec]={"closed","trans","open"};

string fin;
int nFiles(0);
long nTotEv(0);
long currentEvNr(0);

const double pi = acos(-1);

void initHisto();
long processOne(string);
void process();
int findDetector(int &sector, double phi, double r);
void writeOutput();

const int nRings = 7;
double rMin[nRings][nSec], rMax[nRings][nSec];
int setSegmentation();

string tgt_gen_config = "LH2_beam_V27";

void processPlots(const string& finName = "./remollout.root"){
  if(!setSegmentation())
    return;
  fin = finName;
  initHisto();
  process();
  writeOutput();
}

void process(){
  if(fin==""){
    cout<<"\t did not find input file. Quitting!"<<endl;
    return 2;
  }

  if(fin.find(".root") < fin.size()){
    cout<<"Proessing single file:\n\t"<<fin<<endl;
    nTotEv+=processOne(fin);
    nFiles=1;
  }else{
    cout<<"Attempting to process list of output from\n\t"<<fin<<endl;
    ifstream ifile(fin.c_str());
    string data;
    while(ifile>>data){
      cout<<"processing: "<<data<<endl;
      nTotEv+=processOne(data);
      nFiles++;
    }
  }
  cout<<"\nFinished processing a total of "<<nTotEv<<endl;
}

void initHisto(){
  string foutNm = Form("./rootfiles/%s_SAM_radiation_Study.root",tgt_gen_config.c_str());
  fout = new TFile(foutNm.c_str(),"RECREATE");
  fout->mkdir("det28","histos for main detector analysis");
  fout->cd("det28");

  const string secNm[3]={"closed","transition","open"};
  for(int i=0;i<6;i++){
    vector<TH1D*> hr1,hr2;
    vector<TH2D*> hxy1,hxy2;
    for(int j=0;j<3;j++){
      hr1.push_back(new TH1D(Form("hR_e1_R%d_S%d",i+1,j),Form("rate weighted Radial dist. for Ring %d Sector %s;r [mm]",i+1,secNm[j].c_str()),1000,0,2000));
      hr2.push_back(new TH1D(Form("hR_eP1_R%d_S%d",i+1,j),Form("rate weighted Radial dist. for Ring %d Sector %s;r [mm]",i+1,secNm[j].c_str()),1000,0,2000));
      hxy1.push_back(new TH2D(Form("hXY_e1_R%d_S%d",i+1,j),Form("rate weighted XY dist. for Ring %d Sector %s;x [m]; y [mm]",i+1,secNm[j].c_str()),2000,-2000,2000,2000,-2000,2000));
      hxy2.push_back(new TH2D(Form("hXY_eP1_R%d_S%d",i+1,j),Form("rate weighted XY dist. for Ring %d Sector %s;x [mm]; y [mm]",i+1,secNm[j].c_str()),2000,-2000,2000,2000,-2000,2000));
    }
    hR_e1.push_back(hr1);
    hR_eP1.push_back(hr2);
    hXY_e1.push_back(hxy1);
    hXY_eP1.push_back(hxy2);
  }
}

long processOne(string fnm){
  TFile *fin = TFile::Open(fnm.c_str(),"READ");
//  if(!fin->IsOpen() || fin->IsZombie()){
//    cout<<"Problem: can't find file: "<<fnm<<endl;
//    fin->Close();
//    delete fin;
//    return 0;
//  }
  TTree *t = (TTree*)fin->Get("T");
  if(t==0) return 0;
  Double_t rate = 0;
  std::vector<remollGenericDetectorHit_t> *hit =0;
  
  t->SetBranchAddress("rate", &rate);
  t->SetBranchAddress("hit", &hit);

  long nEntries = t->GetEntries();
  cout<<"\tTotal events: "<<nEntries<<endl;
  float currentProc = 1, procStep = 60;
  vector<int> procID;
  int sector(-1);
  double pi = acos(-1);
  
  for(Long64_t event=0;event<nEntries;t->GetEntry(event++)){
    currentEvNr++;
    if(float(event+1)/nEntries*100 > currentProc){
      cout<<"at tree entry\t"<<event<<"\t"<<float(event+1)/nEntries*100<<endl;
      currentProc+=procStep;
    }

    for(int j=0;j<hit->size();j++){
      if(std::isnan(rate) || std::isinf(rate)) continue;
      if(rate==0) {rate = 1;};
      int sp = spM[int(abs(hit->at(j).pid))]-1;
      if(sp==-1) continue;
      int dt = dtM[int(hit->at(j).det)]-1;
      if(dt==-1) continue;
      double phi = atan2(hit->at(j).y,hit->at(j).x);
      if(phi<0) phi+=2*pi;
      int foundRing = findDetector(sector, phi, hit->at(j).r);
      if(foundRing==-1) continue;

      if(hit->at(j).e>1 && (abs(hit->at(j).pid)==11 || abs(hit->at(j).pid)==211)){
        if(dt==0){
          if(foundRing<6){
            hR_e1[foundRing][sector]->Fill(hit->at(j).r,rate);
            hXY_e1[foundRing][sector]->Fill(hit->at(j).x,hit->at(j).y,rate);
          }
        }
      }

      if(hit->at(j).trid==1 || hit->at(j).trid==2){
        if(dt==0){
          if(foundRing<6){
            hR_eP1[foundRing][sector]->Fill(hit->at(j).r,rate);
            hXY_eP1[foundRing][sector]->Fill(hit->at(j).x,hit->at(j).y,rate);
          }
        }
      }
    }
  }
  fin->Close();
  delete fin;
  return nEntries;
  };

void writeOutput(){
  fout->cd("det28");
  for(int i=0;i<6;i++){
    for(int j=0;j<3;j++){
      hR_e1[i][j]->Scale(1./nTotEv);
      hR_e1[i][j]->Write();
      hR_eP1[i][j]->Scale(1./nTotEv);
      hR_eP1[i][j]->Write();
      hXY_e1[i][j]->Scale(1./nTotEv);
      hXY_e1[i][j]->Write();
      hXY_eP1[i][j]->Scale(1./nTotEv);
      hXY_eP1[i][j]->Write();
    }
  }

  fout->Close();
}



int findDetector(int &sector, double phi, double r){
  const double secPhi = fmod(phi, 2*pi/7);
  //0, 1, 2 == closed, transition, open
  if(secPhi<pi/28)
    sector = 0;
  else if(secPhi<3*pi/28)
    sector = 1;
  else if(secPhi<5*pi/28)
    sector = 2;
  else if(secPhi<7*pi/28)
    sector = 1;
  else if(secPhi<8*pi/28)
    sector = 0;

//this will just pick up the first hit and ignore the rest if there is overlap in tiling
  for(int i=0;i<nRings;i++)
    if(r>=rMin[i][sector] && r<=rMax[i][sector])
      return i;
  
  return -1;
}

int setSegmentation(){
  double rmin[nRings][nSec]={
    { 660,  660,  660},
    { 690,  690,  690},
    { 750,  750,  750},
    { 810,  810,  810},
    { 930,  930,  930},
    {1070, 1070, 1070},
    {1170, 1170, 1170}
  };
  double rmax[nRings][nSec]={
    { 690,  690,  690},
    { 750,  750,  750},
    { 810,  810,  810},
    { 930,  930,  930},
    {1070, 1070, 1070},
    {1170, 1170, 1170},
    {1500, 1500, 1500}
  };
  for(int i=0;i<nRings;i++){
    for(int j=0;j<nSec;j++){
      rMin[i][j]=rmin[i][j];
      rMax[i][j]=rmax[i][j];
    }
  }
  return 1;
}
