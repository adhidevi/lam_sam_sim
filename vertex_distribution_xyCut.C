//This macro produces vertex distribution for various virtual detectors
#include "remolltypes.hh"
#include <sstream>
#include <iostream>
#include <fstream>
void vertex_distribution_xyCut(){
  gROOT->Reset();
  gStyle->SetOptStat(0);
  gStyle->SetTitleYOffset(1.3);
  gStyle->SetPadGridX(1);
  gStyle->SetPadGridY(1);
  TGaxis::SetMaxDigits(3);

  double kinEcut = 0;//MeV
  string spTit[] = {"e-/#pi-","e+/#pi+","#gamma","neutron"};
  const int nSp = sizeof(spTit)/sizeof(*spTit);
  for(int iSp=0;iSp<nSp;iSp++){
    spTit[iSp] = Form("%s (KE>=%.0f MeV) ",spTit[iSp].c_str(),kinEcut);
  }
  const string spH[nSp] = {"epiM","epiP","g","n"};
  map<int,int> spM {{11,1},{-211,1},{-11,2},{211,2},{22,3},{2112,4}};

///Change the following lines for which detectors you want to include////
  string detH[] = {"det1751"};
//LAM
  const int nDet = sizeof(detH)/sizeof(*detH);
  const int Det[nDet] = {1751};
  map<int,int> dtM {{1751,1}};
////////////////////////////////////////////////////////////////////////

  const double x_min = -3000;
  const double x_max = 3000;
  const double z_min = -5000;
  const double z_max = 30000;
  const double x_minXm[nDet] = {862.21-66.55};
  const double x_maxXm[nDet] = {862.21};
  const double y_minXm[nDet] = {-2108.28};
  const double y_maxXm[nDet] = {-2108.28+85.0};
  const double x_minYm[nDet] = {-56.39/2.0};
  const double x_maxYm[nDet] = {56.39/2.0};
  const double y_minYm[nDet] = {-2517.07};
  const double y_maxYm[nDet] = {-2517.07+54.10};
  const int nbin = 1000;
  const string motor[] = {"xMotor","yMotor"};
  const int nMot = sizeof(motor)/sizeof(*motor);
  TH2F* h_VxVz[nSp][nDet][nMot];
  TH2F* h_VxVzPzG0[nSp][nDet][nMot];
  TH2F* h_VxVzPzL0[nSp][nDet][nMot];
  TH2F* h_VyVz[nSp][nDet][nMot];
  TH2F* h_VyVzPzG0[nSp][nDet][nMot];
  TH2F* h_VyVzPzL0[nSp][nDet][nMot];
  TObjArray Hlist[nDet];

///Change the following lines as needed////
  const string geometry = "develop_new";
  const string tgt_gen_config = "LH2_beam_V22";
  const string plotType = Form("vertex_distribution_xyCut_EG%.0f",kinEcut);
  int beamGen(1);

///Change this line for appropriate rootfile directory////
  TString rootfile_dir = "/volatile/halla/parity/adhidevi/remoll_rootfiles/develop_br";
//////////////////////////////////////////////////////////

  for(int iSp=0;iSp<nSp;iSp++){
   for(int iDet=0;iDet<nDet;iDet++){
    for(int iMot=0;iMot<nMot;iMot++){
      string titleXZ = Form("%s vertex dist. on %s plane (%s);vz (mm);vx (mm);hits/#thrownEvents",spTit[iSp].c_str(),detH[iDet].c_str(),tgt_gen_config.c_str());
      string titleYZ = Form("%s vertex dist. on %s plane (%s);vz (mm);vy (mm);hits/#thrownEvents",spTit[iSp].c_str(),detH[iDet].c_str(),tgt_gen_config.c_str());
      h_VxVz[iSp][iDet][iMot] = new TH2F(Form("%s_VxVz_%s_%s",detH[iDet].c_str(),spH[iSp].c_str(),motor[iMot].c_str()),titleXZ.c_str(),nbin,z_min,z_max,nbin,x_min,x_max);
      h_VxVzPzG0[iSp][iDet][iMot] = new TH2F(Form("%s_VxVzPzG0_%s_%s",detH[iDet].c_str(),spH[iSp].c_str(),motor[iMot].c_str()),titleXZ.c_str(),nbin,z_min,z_max,nbin,x_min,x_max);
      h_VxVzPzL0[iSp][iDet][iMot] = new TH2F(Form("%s_VxVzPzL0_%s_%s",detH[iDet].c_str(),spH[iSp].c_str(),motor[iMot].c_str()),titleXZ.c_str(),nbin,z_min,z_max,nbin,x_min,x_max);
      h_VyVz[iSp][iDet][iMot] = new TH2F(Form("%s_VyVz_%s_%s",detH[iDet].c_str(),spH[iSp].c_str(),motor[iMot].c_str()),titleYZ.c_str(),nbin,z_min,z_max,nbin,x_min,x_max);
      h_VyVzPzG0[iSp][iDet][iMot] = new TH2F(Form("%s_VyVzPzG0_%s_%s",detH[iDet].c_str(),spH[iSp].c_str(),motor[iMot].c_str()),titleYZ.c_str(),nbin,z_min,z_max,nbin,x_min,x_max);
      h_VyVzPzL0[iSp][iDet][iMot] = new TH2F(Form("%s_VyVzPzL0_%s_%s",detH[iDet].c_str(),spH[iSp].c_str(),motor[iMot].c_str()),titleYZ.c_str(),nbin,z_min,z_max,nbin,x_min,x_max);
    }
   }
  }

  int nfile=0;
  Long64_t nentry=0;
  long nTotEv=0;
  for(int ifile=1001;ifile<=6000;ifile++){
///Change this line for appropriate rootfiles////
    string infile = Form("%s/%s/%s_%d.root",rootfile_dir.Data(),tgt_gen_config.c_str(),tgt_gen_config.c_str(),ifile);
/////////////////////////////////////////////
    TFile *fin = TFile::Open(infile.c_str(),"READ");
    if(!fin || fin->IsZombie()){
      delete fin;
      continue;
      cout<<Form("Skipping %s. Recovered file.",infile.c_str())<<endl;
    }
    if(fin->TestBit(TFile::kRecovered)){
      continue;
      cout<<Form("Skipping %s. Recovered file.",infile.c_str())<<endl;
    }
    nfile++;
   
    TTree *T = (TTree*)fin->Get("T");
    if(T==0) return 0;

    nentry = T->GetEntries();
    cout<<Form("Found %lld entries in filesplit %d",nentry,nfile)<<endl;
    nTotEv+=nentry;

    std::vector<remollGenericDetectorHit_t> *hit =0;
    std::vector<remollEventParticle_t> *part = 0;
    remollBeamTarget_t *bm = 0;
    remollEvent_t *ev = 0;
    Double_t rate = 0;
    T->SetBranchAddress("hit", & hit);
    T->SetBranchAddress("part", & part);
    T->SetBranchAddress("bm", & bm);
    T->SetBranchAddress("ev",& ev);
    T->SetBranchAddress("rate", & rate);
   
    for(Long64_t ientry=0;ientry<nentry;ientry++){
      T->GetEntry(ientry);
      for(int j=0;j<hit->size();j++){
        if(std::isnan(rate) || std::isinf(rate)) continue;
        if(beamGen) rate = 1.0;

        int sp = spM[int(hit->at(j).pid)]-1;
        if(sp==-1) continue;
        int dt = dtM[int(hit->at(j).det)]-1;
        if(dt==-1) continue;
        if(hit->at(j).k<kinEcut) continue;

        if((hit->at(j).x>=x_minXm[dt] && hit->at(j).x<=x_maxXm[dt]) && (hit->at(j).y>=y_minXm[dt] && hit->at(j).y<=y_maxXm[dt])){
          h_VxVz[sp][dt][0]->Fill(hit->at(j).vz,hit->at(j).vx,rate);
          h_VyVz[sp][dt][0]->Fill(hit->at(j).vz,hit->at(j).vy,rate);
          if(hit->at(j).pz>=0){
            h_VxVzPzG0[sp][dt][0]->Fill(hit->at(j).vz,hit->at(j).vx,rate);
            h_VyVzPzG0[sp][dt][0]->Fill(hit->at(j).vz,hit->at(j).vy,rate);
          }else{
            h_VxVzPzL0[sp][dt][0]->Fill(hit->at(j).vz,hit->at(j).vx,rate);
            h_VyVzPzL0[sp][dt][0]->Fill(hit->at(j).vz,hit->at(j).vy,rate);
          }
        }
        if((hit->at(j).x>=x_minYm[dt] && hit->at(j).x<=x_maxYm[dt]) && (hit->at(j).y>=y_minYm[dt] && hit->at(j).y<=y_maxYm[dt])){
          h_VxVz[sp][dt][1]->Fill(hit->at(j).vz,hit->at(j).vx,rate);
          h_VyVz[sp][dt][1]->Fill(hit->at(j).vz,hit->at(j).vy,rate);
          if(hit->at(j).pz>=0){
            h_VxVzPzG0[sp][dt][1]->Fill(hit->at(j).vz,hit->at(j).vx,rate);
            h_VyVzPzG0[sp][dt][1]->Fill(hit->at(j).vz,hit->at(j).vy,rate);
          }else{
            h_VxVzPzL0[sp][dt][1]->Fill(hit->at(j).vz,hit->at(j).vx,rate);
            h_VyVzPzL0[sp][dt][1]->Fill(hit->at(j).vz,hit->at(j).vy,rate);
          }
        }
      }
    }
    delete T;
    delete fin;
  }

  cout<<Form("Total number of file splits: %d",nfile)<<endl;
  cout<<Form("Total number of entries: %ld",nTotEv)<<endl;
   
//////////////////////////////////////////
  TFile* outfile = new TFile(Form("./rootfiles/%s_%s_%s_xyCut.root",geometry.c_str(),tgt_gen_config.c_str(),plotType.c_str()),"recreate");
  for(int iDet=0;iDet<nDet;iDet++){
    outfile->mkdir(Form("%s",detH[iDet].c_str()));
    outfile->cd(Form("%s",detH[iDet].c_str()));
    for(int iSp=0;iSp<nSp;iSp++){
      for(int iMot=0;iMot<nMot;iMot++){
       if(beamGen){
         h_VxVz[iSp][iDet][iMot]->Scale(1.0/nTotEv);
         h_VxVzPzG0[iSp][iDet][iMot]->Scale(1.0/nTotEv);
         h_VxVzPzL0[iSp][iDet][iMot]->Scale(1.0/nTotEv);
         h_VyVz[iSp][iDet][iMot]->Scale(1.0/nTotEv);
         h_VyVzPzG0[iSp][iDet][iMot]->Scale(1.0/nTotEv);
         h_VyVzPzL0[iSp][iDet][iMot]->Scale(1.0/nTotEv);
       }else{
         h_VxVz[iSp][iDet][iMot]->Scale(1.0e-9/nfile);
         h_VxVzPzG0[iSp][iDet][iMot]->Scale(1.0e-9/nfile);
         h_VxVzPzL0[iSp][iDet][iMot]->Scale(1.0e-9/nfile);
         h_VyVz[iSp][iDet][iMot]->Scale(1.0e-9/nfile);
         h_VyVzPzG0[iSp][iDet][iMot]->Scale(1.0e-9/nfile);
         h_VyVzPzL0[iSp][iDet][iMot]->Scale(1.0e-9/nfile);
       }
         h_VxVz[iSp][iDet][iMot]->Write();
         h_VxVzPzG0[iSp][iDet][iMot]->Write();
         h_VxVzPzL0[iSp][iDet][iMot]->Write();
         h_VyVz[iSp][iDet][iMot]->Write();
         h_VyVzPzG0[iSp][iDet][iMot]->Write();
         h_VyVzPzL0[iSp][iDet][iMot]->Write();
      }
    }
  }
  outfile->Close();
}
