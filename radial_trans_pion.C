#include "remolltypes.hh"
#include <sstream>
#include <iostream>
#include <fstream>
const double pi = TMath::Pi();
void radial_trans_pion(){
  gROOT->Reset();
  gStyle->SetOptStat(0);
  gStyle->SetTitleYOffset(1.3);
  gStyle->SetPadGridX(1);
  gStyle->SetPadGridY(1);
  TGaxis::SetMaxDigits(3);
  
  double kinEcut = 1;//MeV
  const string spTit[] = {"#pi-","#pi+","e-","e+"};
  const int nSp = sizeof(spTit)/sizeof(*spTit);
  const string spH[nSp] = {"pi_minus","pi_plus","e_minus","e_plus"};
  map<int,int> spM {{-211,1},{211,2},{11,3},{-11,4}};

///Change the following lines for which detectors you want to include////
  string detH[] = {"det79"};
  const int nDet = sizeof(detH)/sizeof(*detH);
  const int Det[nDet] = {79};
  map<int,int> dtM {{79,1}};
////////////////////////////////////////////////////////////////////////

///Change the following lines as needed////
  const string geometry = "develop_new_pion";
  const string tgt_gen_config = "LH2_beam_V21b";
  const string plotType = Form("radial_trans_th_p_EG%.0f",kinEcut);
  int beamGen(1);

///Change this line for appropriate rootfile directory////
  TString rootfile_dir = "/volatile/halla/parity/adhidevi/remoll_rootfiles/develop_br";
//////////////////////////////////////////////////////////

  double x_min = 0;
  double x_max = 900;
  const int nbin = 200;
  double bin_width = (x_max-x_min)/nbin;
  const string weight = "rate";
  string weight_unit;
  if(beamGen)
    weight_unit = "hits/#thrownEvents";
  else
    weight_unit = "GHz";

  TH1F* h_r[nSp][nDet];
  TH1F* h_rPzG0[nSp][nDet];
  TH1F* h_rPzL0[nSp][nDet];
  TH2F* h_xy[nSp][nDet];
  TH2F* h_xyPzG0[nSp][nDet];
  TH2F* h_xyPzL0[nSp][nDet];
  TH1F* h_p[nSp][nDet];
  TH1F* h_pPzG0[nSp][nDet];
  TH1F* h_pPzL0[nSp][nDet];
  TH1F* h_th[nSp][nDet];
  TH1F* h_thPzG0[nSp][nDet];
  TH1F* h_thPzL0[nSp][nDet];
  TH2F* h_pTh[nSp][nDet];
  TH2F* h_pThPzG0[nSp][nDet];
  TH2F* h_pThPzL0[nSp][nDet];

  for(int iSp=0;iSp<nSp;iSp++){
   for(int iDet=0;iDet<nDet;iDet++){
    string titleR = Form("%s (KE>=%.0f MeV) Radial dist. on %s plane (%s);Radius (mm);%s/%.1fmm",spTit[iSp].c_str(),kinEcut,detH[iDet].c_str(),tgt_gen_config.c_str(),weight_unit.c_str(),bin_width);
    h_r[iSp][iDet] = new TH1F(Form("%s_r_%s_%s",detH[iDet].c_str(),spH[iSp].c_str(),weight.c_str()),titleR.c_str(),nbin,x_min,x_max);
    h_rPzG0[iSp][iDet] = new TH1F(Form("%s_rPzG0_%s_%s",detH[iDet].c_str(),spH[iSp].c_str(),weight.c_str()),titleR.c_str(),nbin,x_min,x_max);
    h_rPzL0[iSp][iDet] = new TH1F(Form("%s_rPzL0_%s_%s",detH[iDet].c_str(),spH[iSp].c_str(),weight.c_str()),titleR.c_str(),nbin,x_min,x_max);

    string titleXY = Form("%s (KE>=%.0f MeV) XY dist. on %s plane (%s);x (mm);y (mm)",spTit[iSp].c_str(),kinEcut,detH[iDet].c_str(),tgt_gen_config.c_str());
    h_xy[iSp][iDet] = new TH2F(Form("%s_xy_%s_%s",detH[iDet].c_str(),spH[iSp].c_str(),weight.c_str()),titleXY.c_str(),nbin,-x_max,x_max,nbin,-x_max,x_max);
    h_xyPzG0[iSp][iDet] = new TH2F(Form("%s_xyPzG0_%s_%s",detH[iDet].c_str(),spH[iSp].c_str(),weight.c_str()),titleXY.c_str(),nbin,-x_max,x_max,nbin,-x_max,x_max);
    h_xyPzL0[iSp][iDet] = new TH2F(Form("%s_xyPzL0_%s_%s",detH[iDet].c_str(),spH[iSp].c_str(),weight.c_str()),titleXY.c_str(),nbin,-x_max,x_max,nbin,-x_max,x_max);

   string titleTh = Form("%s (KE>=%.0f MeV) Theta dist. on %s plane (%s);Theta (deg);hits/#thrownEvents",spTit[iSp].c_str(),kinEcut,detH[iDet].c_str(),tgt_gen_config.c_str());
   h_th[iSp][iDet] = new TH1F(Form("%s_th_%s_%s",detH[iDet].c_str(),spH[iSp].c_str(),weight.c_str()),titleTh.c_str(),200,0,90);
   h_thPzG0[iSp][iDet] = new TH1F(Form("%s_thPzG0_%s_%s",detH[iDet].c_str(),spH[iSp].c_str(),weight.c_str()),titleTh.c_str(),200,0,90);
   h_thPzL0[iSp][iDet] = new TH1F(Form("%s_thPzL0_%s_%s",detH[iDet].c_str(),spH[iSp].c_str(),weight.c_str()),titleTh.c_str(),200,0,90);

   string titleP = Form("%s (KE>=%.0f MeV) Momentum dist. on %s plane (%s);Momentum (MeV);hits/#thrownEvents",spTit[iSp].c_str(),kinEcut,detH[iDet].c_str(),tgt_gen_config.c_str());
   h_p[iSp][iDet] = new TH1F(Form("%s_p_%s_%s",detH[iDet].c_str(),spH[iSp].c_str(),weight.c_str()),titleP.c_str(),1200,0,12000);
   h_pPzG0[iSp][iDet] = new TH1F(Form("%s_pPzG0_%s_%s",detH[iDet].c_str(),spH[iSp].c_str(),weight.c_str()),titleP.c_str(),1200,0,12000);
   h_pPzL0[iSp][iDet] = new TH1F(Form("%s_pPzL0_%s_%s",detH[iDet].c_str(),spH[iSp].c_str(),weight.c_str()),titleP.c_str(),1200,0,12000);

    string titlePTh = Form("%s (KE>=%.0f MeV) P vs Th dist. on %s plane (%s);Theta (deg);Momentum (MeV)",spTit[iSp].c_str(),kinEcut,detH[iDet].c_str(),tgt_gen_config.c_str());
    h_pTh[iSp][iDet] = new TH2F(Form("%s_pTh_%s_%s",detH[iDet].c_str(),spH[iSp].c_str(),weight.c_str()),titlePTh.c_str(),200,0,90,1200,0,12000);
    h_pThPzG0[iSp][iDet] = new TH2F(Form("%s_pThPzG0_%s_%s",detH[iDet].c_str(),spH[iSp].c_str(),weight.c_str()),titlePTh.c_str(),200,0,90,1200,0,12000);
    h_pThPzL0[iSp][iDet] = new TH2F(Form("%s_pThPzL0_%s_%s",detH[iDet].c_str(),spH[iSp].c_str(),weight.c_str()),titlePTh.c_str(),200,0,90,1200,0,12000);

    h_r[iSp][iDet]->Sumw2();
    h_rPzG0[iSp][iDet]->Sumw2();
    h_rPzL0[iSp][iDet]->Sumw2();
    h_p[iSp][iDet]->Sumw2();
    h_pPzG0[iSp][iDet]->Sumw2();
    h_pPzL0[iSp][iDet]->Sumw2();
    h_th[iSp][iDet]->Sumw2();
    h_thPzG0[iSp][iDet]->Sumw2();
    h_thPzL0[iSp][iDet]->Sumw2();
   }
  }

  int nfile=0;
  Long64_t nentry=0;
  long nTotEv=0;
  for(int ifile=1001;ifile<=6000;ifile++){
///Change this line for appropriate rootfiles////
    string infile = Form("%s/%s/%s_%d.root",rootfile_dir.Data(),tgt_gen_config.c_str(),tgt_gen_config.c_str(),ifile);
//////////////////////////////////////////////
    ifstream inf(infile.c_str());
    if(!inf){
      cout<<Form("Skipping %s. File doesn't exist.",infile.c_str())<<endl;
      continue;
    }
    TFile *fin = TFile::Open(infile.c_str(),"READ");
    if(fin->TestBit(TFile::kRecovered)){
      cout<<Form("Skipping %s. Recovered file.",infile.c_str())<<endl;
      continue;
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

        double px = hit->at(j).px;
        double py = hit->at(j).py;
        double pz = hit->at(j).pz;
        double theta = (180.0/pi)*atan2(sqrt(px*px+py*py),pz);

        h_r[sp][dt]->Fill(hit->at(j).r,rate);
        h_xy[sp][dt]->Fill(hit->at(j).x,hit->at(j).y,rate);
        h_p[sp][dt]->Fill(hit->at(j).p,rate);
        h_th[sp][dt]->Fill(theta,rate);
        h_pTh[sp][dt]->Fill(theta,hit->at(j).p,rate);

        if(hit->at(j).pz>=0){
          h_rPzG0[sp][dt]->Fill(hit->at(j).r,rate);
          h_xyPzG0[sp][dt]->Fill(hit->at(j).x,hit->at(j).y,rate);
          h_pPzG0[sp][dt]->Fill(hit->at(j).p,rate);
          h_thPzG0[sp][dt]->Fill(theta,rate);
          h_pThPzG0[sp][dt]->Fill(theta,hit->at(j).p,rate);
        }else{
          h_rPzL0[sp][dt]->Fill(hit->at(j).r,rate);
          h_xyPzL0[sp][dt]->Fill(hit->at(j).x,hit->at(j).y,rate);
          h_pPzL0[sp][dt]->Fill(hit->at(j).p,rate);
          h_thPzL0[sp][dt]->Fill(theta,rate);
          h_pThPzL0[sp][dt]->Fill(theta,hit->at(j).p,rate);
        }

      }
    }
    delete T;
    delete fin;
  }

  cout<<Form("Total number of file splits: %d",nfile)<<endl;
  cout<<Form("Total number of entries: %ld",nTotEv)<<endl;
   
//////////////////////////////////////////
  TFile* outfile = new TFile(Form("./rootfiles/%s_%s_%s.root",geometry.c_str(),tgt_gen_config.c_str(),plotType.c_str()),"recreate");

  for(int iDet=0;iDet<nDet;iDet++){
    outfile->mkdir(Form("%s",detH[iDet].c_str()));
    outfile->cd(Form("%s",detH[iDet].c_str()));
    for(int iSp=0;iSp<nSp;iSp++){
     if(beamGen){ 
      h_r[iSp][iDet]->Scale(1.0/nTotEv);
      h_rPzG0[iSp][iDet]->Scale(1.0/nTotEv);
      h_rPzL0[iSp][iDet]->Scale(1.0/nTotEv);
      h_xy[iSp][iDet]->Scale(1.0/nTotEv);
      h_xyPzG0[iSp][iDet]->Scale(1.0/nTotEv);
      h_xyPzL0[iSp][iDet]->Scale(1.0/nTotEv);
      h_p[iSp][iDet]->Scale(1.0/nTotEv);
      h_pPzG0[iSp][iDet]->Scale(1.0/nTotEv);
      h_pPzL0[iSp][iDet]->Scale(1.0/nTotEv);
      h_th[iSp][iDet]->Scale(1.0/nTotEv);
      h_thPzG0[iSp][iDet]->Scale(1.0/nTotEv);
      h_thPzL0[iSp][iDet]->Scale(1.0/nTotEv);
      h_pTh[iSp][iDet]->Scale(1.0/nTotEv);
      h_pThPzG0[iSp][iDet]->Scale(1.0/nTotEv);
      h_pThPzL0[iSp][iDet]->Scale(1.0/nTotEv);
     }else{
      h_r[iSp][iDet]->Scale(1.0e-9/nfile);
      h_rPzG0[iSp][iDet]->Scale(1.0e-9/nfile);
      h_rPzL0[iSp][iDet]->Scale(1.0e-9/nfile);
      h_xy[iSp][iDet]->Scale(1.0e-9/nfile);
      h_xyPzG0[iSp][iDet]->Scale(1.0e-9/nfile);
      h_xyPzL0[iSp][iDet]->Scale(1.0e-9/nfile);
      h_p[iSp][iDet]->Scale(1.0e-9/nfile);
      h_pPzG0[iSp][iDet]->Scale(1.0e-9/nfile);
      h_pPzL0[iSp][iDet]->Scale(1.0e-9/nfile);
      h_th[iSp][iDet]->Scale(1.0e-9/nfile);
      h_thPzG0[iSp][iDet]->Scale(1.0e-9/nfile);
      h_thPzL0[iSp][iDet]->Scale(1.0e-9/nfile);
      h_pTh[iSp][iDet]->Scale(1.0e-9/nfile);
      h_pThPzG0[iSp][iDet]->Scale(1.0e-9/nfile);
      h_pThPzL0[iSp][iDet]->Scale(1.0e-9/nfile);
     }
      h_r[iSp][iDet]->Write();
      h_rPzG0[iSp][iDet]->Write();
      h_rPzL0[iSp][iDet]->Write();
      h_xy[iSp][iDet]->Write();
      h_xyPzG0[iSp][iDet]->Write();
      h_xyPzL0[iSp][iDet]->Write();
      h_p[iSp][iDet]->Write();
      h_pPzG0[iSp][iDet]->Write();
      h_pPzL0[iSp][iDet]->Write();
      h_th[iSp][iDet]->Write();
      h_thPzG0[iSp][iDet]->Write();
      h_thPzL0[iSp][iDet]->Write();
      h_pTh[iSp][iDet]->Write();
      h_pThPzG0[iSp][iDet]->Write();
      h_pThPzL0[iSp][iDet]->Write();
    }
  }
}
