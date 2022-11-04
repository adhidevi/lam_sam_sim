TFile *outfile;
string fin;
int nFiles(0);
long nTotEv(0);
long currentEvNr(0);
void initHisto();
long processOne(string);
void process();
void writeOutput();

vector<vector<TH1D*>> hR, hR_pzG0, hR_pzL0;
vector<vector<TH2D*>> hXY, hXY_pzG0, hXY_pzL0;

const string spTit[] = {"e-/#pi-","e+/#pi+","#gamma","neutron","e-/e+","e- trid==1"};
const int nSp = sizeof(spTit)/sizeof(*spTit);
const string spH[nSp] = {"epiM","epiP","g","n","ee","eTrIdCut"};
map<int,int> spM {{11,1},{-211,1},{-11,2},{211,2},{22,3},{2112,4}};
const int DetNo = 176;
map<int,int> dtM {{DetNo,1}};
const double Ecut[] {0,1};//MeV
const int nEcut = sizeof(Ecut)/sizeof(*Ecut);
int beamGen(1);

string tgt_gen_config = "LH2_beam_V30";
string histR_title;
string histXY_title;

void process_rad_xy(const string& finName = "./remollout.root"){
  fin = finName;
  initHisto();
  process();
  writeOutput();
}

void process(){
  if(fin==""){
    cout<<"\t did not find input file. Quitting!"<<endl;
    exit(0);
  }
  if(fin.find(".root") < fin.size()){
    cout<<"Processing single file:\n\t"<<fin<<endl;
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
  string foutNm = Form("./rootfiles/%s_rad_xy.root",tgt_gen_config.c_str());
  outfile = new TFile(foutNm.c_str(),"RECREATE");
  outfile->mkdir(Form("det%d",DetNo),"histos for SAM analysis");
  outfile->cd(Form("det%d",DetNo));

  for(int iSp=0;iSp<nSp;iSp++){
    vector<TH1D*> hr, hr_pzG0, hr_pzL0;
    vector<TH2D*> hxy, hxy_pzG0, hxy_pzL0;
    for(int iEcut=0;iEcut<nEcut;iEcut++){
      histR_title = Form("%s (E>%.0f MeV) Radial dist. on %s (%s);Radius (mm);hits/%sthrownEvents",spTit[iSp].c_str(),Ecut[iEcut],Form("det%d",DetNo),tgt_gen_config.c_str(),"#");
      histXY_title = Form("%s (E>%.0f MeV) XY dist. on %s (%s);x (mm); y (mm);hits/%sthrownEvents",spTit[iSp].c_str(),Ecut[iEcut],Form("det%d",DetNo),tgt_gen_config.c_str(),"#");

      hr.push_back(new TH1D(Form("%s_R_det%d_EG%.0fMeV",spH[iSp].c_str(),DetNo,Ecut[iEcut]),histR_title.c_str(),700,0,700));
      hr_pzG0.push_back(new TH1D(Form("%s_RpzG0_det%d_EG%.0fMeV",spH[iSp].c_str(),DetNo,Ecut[iEcut]),histR_title.c_str(),700,0,700));
      hr_pzL0.push_back(new TH1D(Form("%s_RpzL0_det%d_EG%.0fMeV",spH[iSp].c_str(),DetNo,Ecut[iEcut]),histR_title.c_str(),700,0,700));
      hxy.push_back(new TH2D(Form("%s_XY_det%d_EG%.0fMeV",spH[iSp].c_str(),DetNo,Ecut[iEcut]),histXY_title.c_str(),1400,-700,700,1400,-700,700));
      hxy_pzG0.push_back(new TH2D(Form("%s_XYpzG0_det%d_EG%.0fMeV",spH[iSp].c_str(),DetNo,Ecut[iEcut]),histXY_title.c_str(),1400,-700,700,1400,-700,700));
      hxy_pzL0.push_back(new TH2D(Form("%s_XYpzL0_det%d_EG%.0fMeV",spH[iSp].c_str(),DetNo,Ecut[iEcut]),histXY_title.c_str(),1400,-700,700,1400,-700,700));
    }
    hR.push_back(hr);
    hR_pzG0.push_back(hr_pzG0);
    hR_pzL0.push_back(hr_pzL0);
    hXY.push_back(hxy);
    hXY_pzG0.push_back(hxy_pzG0);
    hXY_pzL0.push_back(hxy_pzL0);
  }
}

long processOne(string fnm){
  TFile *fin = TFile::Open(fnm.c_str(),"READ");
  TTree *t = (TTree*)fin->Get("T");
  if(t==0) return 0;
  Double_t rate = 0;
  std::vector<remollGenericDetectorHit_t> *hit = 0;
  t->SetBranchAddress("rate", &rate);
  t->SetBranchAddress("hit", &hit);

  long nEntries = t->GetEntries();
  cout<<"\tTotal events: "<<nEntries<<endl;
  float currentProc = 1, procStep = 60;
  
  for(Long64_t event=0;event<nEntries;t->GetEntry(event++)){
    currentEvNr++;
    if(float(event+1)/nEntries*100 > currentProc){
      cout<<"at tree entry\t"<<event<<"\t"<<float(event+1)/nEntries*100<<"%"<<endl;
      currentProc+=procStep;
    }

    for(int j=0;j<hit->size();j++){
      if(std::isnan(rate) || std::isinf(rate)) continue;
      if(rate==0) rate = 1;
      if(beamGen) rate = 1;
      int sp = spM[int(hit->at(j).pid)]-1;
      if(sp==-1) continue;
      int dt = dtM[int(hit->at(j).det)]-1;
      if(dt==-1) continue;
      for(int iEcut=0;iEcut<nEcut;iEcut++){
        if(hit->at(j).e>Ecut[iEcut]){
          hR[sp][iEcut]->Fill(hit->at(j).r,rate);
          hXY[sp][iEcut]->Fill(hit->at(j).x,hit->at(j).y,rate);
          if(hit->at(j).pz>=0){
            hR_pzG0[sp][iEcut]->Fill(hit->at(j).r,rate);
            hXY_pzG0[sp][iEcut]->Fill(hit->at(j).x,hit->at(j).y,rate);
          }else{
            hR_pzL0[sp][iEcut]->Fill(hit->at(j).r,rate);
            hXY_pzL0[sp][iEcut]->Fill(hit->at(j).x,hit->at(j).y,rate);
          }
          if(hit->at(j).pid==11 || hit->at(j).pid==-11){
            hR[4][iEcut]->Fill(hit->at(j).r,rate);
            hXY[4][iEcut]->Fill(hit->at(j).x,hit->at(j).y,rate);
            if(hit->at(j).pz>=0){
              hR_pzG0[4][iEcut]->Fill(hit->at(j).r,rate);
              hXY_pzG0[4][iEcut]->Fill(hit->at(j).x,hit->at(j).y,rate);
            }else{
              hR_pzL0[4][iEcut]->Fill(hit->at(j).r,rate);
              hXY_pzL0[4][iEcut]->Fill(hit->at(j).x,hit->at(j).y,rate);
            }
          }
          if(hit->at(j).trid==1 && hit->at(j).pid==11){
            hR[5][iEcut]->Fill(hit->at(j).r,rate);
            hXY[5][iEcut]->Fill(hit->at(j).x,hit->at(j).y,rate);
            if(hit->at(j).pz>=0){
              hR_pzG0[5][iEcut]->Fill(hit->at(j).r,rate);
              hXY_pzG0[5][iEcut]->Fill(hit->at(j).x,hit->at(j).y,rate);
            }else{
              hR_pzL0[5][iEcut]->Fill(hit->at(j).r,rate);
              hXY_pzL0[5][iEcut]->Fill(hit->at(j).x,hit->at(j).y,rate);
            }
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
  outfile->cd(Form("det%d",DetNo));
  cout<<nTotEv<<endl;
  for(int iSp=0;iSp<nSp;iSp++){
    for(int iEcut=0;iEcut<nEcut;iEcut++){
      hR[iSp][iEcut]->Scale(1./nTotEv);
      hXY[iSp][iEcut]->Scale(1./nTotEv);
      hR[iSp][iEcut]->Write();
      hXY[iSp][iEcut]->Write();
    }
  }
  outfile->Close();
}
