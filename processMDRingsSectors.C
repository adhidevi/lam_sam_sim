const string spTit[] = {"e-/#pi-","e+/#pi+","#gamma","neutron","e/#pi E>1 MeV","e- trid==1 E>1 MeV"};
const int nSp = sizeof(spTit)/sizeof(*spTit);
const string spH[nSp] = {"epiM","epiP","g","n","epiMP1","eTrIdCut1"};
map<int,int> spM{{11,1},{-211,1},{-11,2},{211,2},{22,3},{2112,4}};
const string detH[]={"det28"};
const int nDet=sizeof(detH)/sizeof(*detH);
map<int,int> dtM{{28,1}};

const int nSec = 3;
const string secH[nSec] = {"closed","transition","open"};

TFile *fout;
vector<vector<TH1D*>> hR_epiM,hR_epiP,hR_g,hR_n,hR_epiMP1,hR_eP1;
vector<vector<TH2D*>> hXY_epiM,hXY_epiP,hXY_g,hXY_n,hXY_epiMP1,hXY_eP1;
const int nSecDet = 21; // 7 (ring, including pmts) x 3 (sectors)
TH1D *hR_rate[nSp];
string fin;
int nFiles(0);
long nTotEv(0);

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
int beamGen(1);

void processMDRingsSectors(const string& finName = "./remollout.root"){
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
  string foutNm = Form("./rootfiles/%s_det28_rings_sectors_rate_study.root",tgt_gen_config.c_str());
  fout = new TFile(foutNm.c_str(),"RECREATE");
  fout->mkdir("det28","histos for main detector analysis");
  fout->cd("det28");

  const string secNm[3]={"closed","transition","open"};
  for(int i=0;i<7;i++){
    vector<TH1D*> hr1,hr2,hr3,hr4,hr5,hr6;
    vector<TH2D*> hxy1,hxy2,hxy3,hxy4,hxy5,hxy6;
    for(int j=0;j<3;j++){
      hr1.push_back(new TH1D(Form("hR_epiM_R%d_S%d",i+1,j),Form("e-/#pi- rate weighted Radial dist. for Ring %d Sector %s;r [mm]",i+1,secNm[j].c_str()),1000,0,2000));
      hr2.push_back(new TH1D(Form("hR_epiP_R%d_S%d",i+1,j),Form("e+/#pi+ rate weighted Radial dist. for Ring %d Sector %s;r [mm]",i+1,secNm[j].c_str()),1000,0,2000));
      hr3.push_back(new TH1D(Form("hR_g_R%d_S%d",i+1,j),Form("#gamma rate weighted Radial dist. for Ring %d Sector %s;r [mm]",i+1,secNm[j].c_str()),1000,0,2000));
      hr4.push_back(new TH1D(Form("hR_n_R%d_S%d",i+1,j),Form("neutron rate weighted Radial dist. for Ring %d Sector %s;r [mm]",i+1,secNm[j].c_str()),1000,0,2000));
      hr5.push_back(new TH1D(Form("hR_epiMP1_R%d_S%d",i+1,j),Form("e/#pi E>1 MeV rate weighted Radial dist. for Ring %d Sector %s;r [mm]",i+1,secNm[j].c_str()),1000,0,2000));
      hr6.push_back(new TH1D(Form("hR_eP1_R%d_S%d",i+1,j),Form("e- trid==1 E>1 MeV rate weighted Radial dist. for Ring %d Sector %s;r [mm]",i+1,secNm[j].c_str()),1000,0,2000));
      hxy1.push_back(new TH2D(Form("hXY_epiM_R%d_S%d",i+1,j),Form("e-/#pi- rate weighted XY dist. for Ring %d Sector %s;x [m]; y [mm]",i+1,secNm[j].c_str()),2000,-2000,2000,2000,-2000,2000));
      hxy2.push_back(new TH2D(Form("hXY_epiP_R%d_S%d",i+1,j),Form("e+/#pi+ rate weighted XY dist. for Ring %d Sector %s;x [m]; y [mm]",i+1,secNm[j].c_str()),2000,-2000,2000,2000,-2000,2000));
      hxy3.push_back(new TH2D(Form("hXY_g_R%d_S%d",i+1,j),Form("#gamma rate weighted XY dist. for Ring %d Sector %s;x [m]; y [mm]",i+1,secNm[j].c_str()),2000,-2000,2000,2000,-2000,2000));
      hxy4.push_back(new TH2D(Form("hXY_n_R%d_S%d",i+1,j),Form("neutron rate weighted XY dist. for Ring %d Sector %s;x [m]; y [mm]",i+1,secNm[j].c_str()),2000,-2000,2000,2000,-2000,2000));
      hxy5.push_back(new TH2D(Form("hXY_epiMP1_R%d_S%d",i+1,j),Form("e/#pi E>1 MeV rate weighted XY dist. for Ring %d Sector %s;x [m]; y [mm]",i+1,secNm[j].c_str()),2000,-2000,2000,2000,-2000,2000));
      hxy6.push_back(new TH2D(Form("hXY_eP1_R%d_S%d",i+1,j),Form("e- trid==1 E>1 MeV rate weighted XY dist. for Ring %d Sector %s;x [mm]; y [mm]",i+1,secNm[j].c_str()),2000,-2000,2000,2000,-2000,2000));
    }
    hR_epiM.push_back(hr1);
    hR_epiP.push_back(hr2);
    hR_g.push_back(hr3);
    hR_n.push_back(hr4);
    hR_epiMP1.push_back(hr5);
    hR_eP1.push_back(hr6);
    hXY_epiM.push_back(hxy1);
    hXY_epiP.push_back(hxy2);
    hXY_g.push_back(hxy3);
    hXY_n.push_back(hxy4);
    hXY_epiMP1.push_back(hxy5);
    hXY_eP1.push_back(hxy6);
  }
  for(int iSp=0;iSp<nSp;iSp++){
    hR_rate[iSp] = new TH1D(Form("hR_rate_%s",spH[iSp].c_str()),Form("Sum for all rings and sectors %s",spTit[iSp].c_str()),nSecDet,0,nSecDet);
    for(int k=1;k<=nSecDet;k++){
      int ring = (k-1-(k-1)%3)/3+1;
      int sector = (k-1)%3;
      hR_rate[iSp]->GetXaxis()->SetBinLabel(k,Form("R%d_%s",ring,secNm[sector].c_str()));
    }
  }
}

long processOne(string fnm){
  TFile *fin = TFile::Open(fnm.c_str(),"READ");
  if(!fin->IsOpen() || fin->IsZombie()){
    cout<<"Problem: can't find file: "<<fnm<<endl;
    fin->Close();
    delete fin;
    return 0;
  }
  TTree *t = (TTree*)fin->Get("T");
  if(t==0) return 0;
  Double_t rate = 0;
  std::vector<remollGenericDetectorHit_t> *hit =0;
  
  t->SetBranchAddress("rate", &rate);
  t->SetBranchAddress("hit", &hit);

  long nEntries = t->GetEntries();
  cout<<"\tTotal events: "<<nEntries<<endl;
  vector<int> procID;
  int sector(-1);
  double pi = acos(-1);
  
  for(Long64_t event=0;event<nEntries;t->GetEntry(event++)){

    for(int j=0;j<hit->size();j++){
      if(std::isnan(rate) || std::isinf(rate)) continue;
      if(rate==0) {rate = 1;}
      if(beamGen) {rate = 1;}

      int dt = dtM[int(hit->at(j).det)]-1;
      if(dt==-1) continue;
      double phi = atan2(hit->at(j).y,hit->at(j).x);
      if(phi<0) phi+=2*pi;
      int foundRing = findDetector(sector, phi, hit->at(j).r);
      if(foundRing==-1) continue;

      if(hit->at(j).pid==11 || hit->at(j).pid==-211){
        hR_epiM[foundRing][sector]->Fill(hit->at(j).r,rate);
        hXY_epiM[foundRing][sector]->Fill(hit->at(j).x,hit->at(j).y,rate);
        hR_rate[0]->SetBinContent(foundRing*3+sector+1,rate + hR_rate[0]->GetBinContent(foundRing*3+sector+1));
      }
      if(hit->at(j).pid==-11 || hit->at(j).pid==211){
        hR_epiP[foundRing][sector]->Fill(hit->at(j).r,rate);
        hXY_epiP[foundRing][sector]->Fill(hit->at(j).x,hit->at(j).y,rate);
        hR_rate[1]->SetBinContent(foundRing*3+sector+1,rate + hR_rate[1]->GetBinContent(foundRing*3+sector+1));
      }
      if(hit->at(j).pid==22){
        hR_g[foundRing][sector]->Fill(hit->at(j).r,rate);
        hXY_g[foundRing][sector]->Fill(hit->at(j).x,hit->at(j).y,rate);
        hR_rate[2]->SetBinContent(foundRing*3+sector+1,rate + hR_rate[2]->GetBinContent(foundRing*3+sector+1));
      }
      if(hit->at(j).pid==2112){
        hR_n[foundRing][sector]->Fill(hit->at(j).r,rate);
        hXY_n[foundRing][sector]->Fill(hit->at(j).x,hit->at(j).y,rate);
        hR_rate[3]->SetBinContent(foundRing*3+sector+1,rate + hR_rate[3]->GetBinContent(foundRing*3+sector+1));
      }
      if(hit->at(j).e>1 && (abs(hit->at(j).pid)==11 || abs(hit->at(j).pid)==211)){
        hR_epiMP1[foundRing][sector]->Fill(hit->at(j).r,rate);
        hXY_epiMP1[foundRing][sector]->Fill(hit->at(j).x,hit->at(j).y,rate);
        hR_rate[4]->SetBinContent(foundRing*3+sector+1,rate + hR_rate[4]->GetBinContent(foundRing*3+sector+1));
      }

      if(hit->at(j).e>1 && hit->at(j).trid==1 && hit->at(j).pid==11){
        hR_eP1[foundRing][sector]->Fill(hit->at(j).r,rate);
        hXY_eP1[foundRing][sector]->Fill(hit->at(j).x,hit->at(j).y,rate);
        hR_rate[5]->SetBinContent(foundRing*3+sector+1,rate + hR_rate[5]->GetBinContent(foundRing*3+sector+1));
      }
    }
  }
  fin->Close();
  delete fin;
  return nEntries;
  };

void writeOutput(){
  fout->cd("det28");
cout<<nTotEv<<endl;
  for(int i=0;i<7;i++){
    for(int j=0;j<3;j++){
      hR_epiM[i][j]->Scale(1./nTotEv);
      hR_epiP[i][j]->Scale(1./nTotEv);
      hR_g[i][j]->Scale(1./nTotEv);
      hR_n[i][j]->Scale(1./nTotEv);
      hR_epiMP1[i][j]->Scale(1./nTotEv);
      hR_eP1[i][j]->Scale(1./nTotEv);
      hXY_epiM[i][j]->Scale(1./nTotEv);
      hXY_epiP[i][j]->Scale(1./nTotEv);
      hXY_g[i][j]->Scale(1./nTotEv);
      hXY_n[i][j]->Scale(1./nTotEv);
      hXY_epiMP1[i][j]->Scale(1./nTotEv);
      hXY_eP1[i][j]->Scale(1./nTotEv);
      hR_epiM[i][j]->Write();
      hR_epiP[i][j]->Write();
      hR_g[i][j]->Write();
      hR_n[i][j]->Write();
      hR_epiMP1[i][j]->Write();
      hR_eP1[i][j]->Write();
      hXY_epiM[i][j]->Write();
      hXY_epiP[i][j]->Write();
      hXY_g[i][j]->Write();
      hXY_n[i][j]->Write();
      hXY_epiMP1[i][j]->Write();
      hXY_eP1[i][j]->Write();
    }
  }
  for(int iSp=0;iSp<nSp;iSp++){
    hR_rate[iSp]->Scale(1.0/nTotEv);
    hR_rate[iSp]->Write();
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
