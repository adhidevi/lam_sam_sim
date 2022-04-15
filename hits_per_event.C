#include "remolltypes.hh"
#include <sstream>
#include <iostream>
#include <fstream>
const double pi = TMath::Pi();
const double lam_length = 360.0;//azimuthal length of LAM quartz
const double lam_rin = 1010.0;//inner radius of LAM quartz
double lam_angle = atan(lam_length/lam_rin);
double sep_mid = 2*pi/14.0;

void hits_per_event(){
  gROOT->Reset();
  gStyle->SetOptStat(0);
  gStyle->SetTitleYOffset(1.3);
  gStyle->SetPadGridX(1);
  gStyle->SetPadGridY(1);
  TGaxis::SetMaxDigits(3);

  const int Det = 174;
  const double inRcut = 1010;
  const double outRcut = 1130;
  const double kinEcut = 0;//MeV
  const string spTit[] = {"e-/#pi-","e+/#pi+","#gamma","neutron"};
  const int nSp = sizeof(spTit)/sizeof(*spTit);
  const string spH[nSp] = {"epiM","epiP","g","n"};
  map<int,int> spM {{11,1},{-211,1},{-11,2},{211,2},{22,3},{2112,4}};

  TString rootfile_dir = "/volatile/halla/parity/adhidevi/remoll_rootfiles/develop_br";
  const string geometry = "develop";
  const string tgt_gen_config = "LH2_beam_V16";
  const string plotType = Form("hits_dist_KEG%.0fMeV_split1",kinEcut);
  int beamGen(1);

  const int nbin_hits = 100;
  const double x_min = 0;
  const double x_max = 100;

  const int nbin_r = 500;
  const double r_min = 0;
  const double r_max = 1900;
  double bin_width = (r_max - r_min)/nbin_r;

  TH1D* h_hits[nSp];
  TH1F* h_rate[nSp];
  TH1F* h_ratePzG0[nSp];
  TH1F* h_ratePzL0[nSp];
  TH2F* h_xy[nSp];
  TH2F* h_xyPzG0[nSp];
  TH2F* h_xyPzL0[nSp];

  for(int iSp=0;iSp<nSp;iSp++){
    h_hits[iSp] = new TH1D(Form("hits_%s",spH[iSp].c_str()),Form("%s (KE>%.0f MeV) nhits at det%d plane (%s);nhits;Events",spTit[iSp].c_str(),kinEcut,Det,tgt_gen_config.c_str()),nbin_hits,x_min,x_max);
    h_hits[iSp]->Sumw2();

    string title1D = Form("%s Radial dist. on det%d plane (%s);Radius (mm);hits/#thrownEvents/%.1fmm",spTit[iSp].c_str(),Det,tgt_gen_config.c_str(),bin_width);

    h_rate[iSp] = new TH1F(Form("det%d_r_%s_rate",Det,spH[iSp].c_str()),title1D.c_str(),nbin_r,r_min,r_max);
    h_ratePzG0[iSp] = new TH1F(Form("det%d_rPzG0_%s_rate",Det,spH[iSp].c_str()),title1D.c_str(),nbin_r,r_min,r_max);
    h_ratePzL0[iSp] = new TH1F(Form("det%d_rPzL0_%s_rate",Det,spH[iSp].c_str()),title1D.c_str(),nbin_r,r_min,r_max);

    string title2D = Form("%s XY dist. on det%d plane (%s);x (mm);y (mm)",spTit[iSp].c_str(),Det,tgt_gen_config.c_str());

    h_xy[iSp] = new TH2F(Form("det%d_xy_%s_rate",Det,spH[iSp].c_str()),title2D.c_str(),nbin_r,-r_max,r_max,nbin_r,-r_max,r_max);
    h_xyPzG0[iSp] = new TH2F(Form("det%d_xyPzG0_%s_rate",Det,spH[iSp].c_str()),title2D.c_str(),nbin_r,-r_max,r_max,nbin_r,-r_max,r_max);
    h_xyPzL0[iSp] = new TH2F(Form("det%d_xyPzL0_%s_rate",Det,spH[iSp].c_str()),title2D.c_str(),nbin_r,-r_max,r_max,nbin_r,-r_max,r_max);

    h_rate[iSp]->Sumw2();
    h_ratePzG0[iSp]->Sumw2();
    h_ratePzL0[iSp]->Sumw2();
  }

  int nfile = 0;
  Long64_t nentry = 0;
  long nTotEv = 0;

  for(int ifile=1001;ifile<=3500;ifile++){
    string infile = Form("%s/%s/%s_%d.root",rootfile_dir.Data(),tgt_gen_config.c_str(),tgt_gen_config.c_str(),ifile);

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

    std::vector<remollGenericDetectorHit_t> *hit = 0;
    Double_t rate = 0;
    T->SetBranchAddress("hit", &hit);
    T->SetBranchAddress("rate", &rate);

    for(int ientry=0;ientry<nentry;ientry++){
      T->GetEntry(ientry);
      int eMcounter = 0;
      int ePcounter = 0;
      int gcounter = 0;
      int ncounter = 0;
      for(int j=0;j<hit->size();j++){
        if(std::isnan(rate) || std::isinf(rate)) continue;
        if(beamGen) rate = 1.0;

        double phi = hit->at(j).ph;
        if(!((abs(phi)>=1*sep_mid-lam_angle/2.0 && abs(phi)<=1*sep_mid+lam_angle/2.0) ||
            (abs(phi)>=3*sep_mid-lam_angle/2.0 && abs(phi)<=3*sep_mid+lam_angle/2.0) ||
            (abs(phi)>=5*sep_mid-lam_angle/2.0 && abs(phi)<=5*sep_mid+lam_angle/2.0) ||
            (abs(phi)>=7*sep_mid-lam_angle/2.0 && abs(phi)<=7*sep_mid+lam_angle/2.0) ||
            (abs(phi)>=9*sep_mid-lam_angle/2.0 && abs(phi)<=9*sep_mid+lam_angle/2.0) ||
            (abs(phi)>=11*sep_mid-lam_angle/2.0 && abs(phi)<=11*sep_mid+lam_angle/2.0) ||
            (abs(phi)>=13*sep_mid-lam_angle/2.0 && abs(phi)<=13*sep_mid+lam_angle/2.0))
          ) continue;
        
        if(hit->at(j).k<kinEcut) continue;
        if(hit->at(j).det!=Det) continue;
        if(hit->at(j).r<inRcut || hit->at(j).r>outRcut) continue;

        if(hit->at(j).pid==11 || hit->at(j).pid==-211){
          eMcounter++;
          h_rate[0]->Fill(hit->at(j).r,rate);
          h_xy[0]->Fill(hit->at(j).x,hit->at(j).y,rate);
          if(hit->at(j).pz>=0){
            h_ratePzG0[0]->Fill(hit->at(j).r,rate);
            h_xyPzG0[0]->Fill(hit->at(j).x,hit->at(j).y,rate);
          }else{
            h_ratePzL0[0]->Fill(hit->at(j).r,rate);
            h_xyPzL0[0]->Fill(hit->at(j).x,hit->at(j).y,rate);
          }
        }else if(hit->at(j).pid==-11 || hit->at(j).pid==211){
          ePcounter++;
          h_rate[1]->Fill(hit->at(j).r,rate);
          h_xy[1]->Fill(hit->at(j).x,hit->at(j).y,rate);
          if(hit->at(j).pz>=0){
            h_ratePzG0[1]->Fill(hit->at(j).r,rate);
            h_xyPzG0[1]->Fill(hit->at(j).x,hit->at(j).y,rate);
          }else{
            h_ratePzL0[1]->Fill(hit->at(j).r,rate);
            h_xyPzL0[1]->Fill(hit->at(j).x,hit->at(j).y,rate);
          }
        }else if(hit->at(j).pid==22){
          gcounter++;
          h_rate[2]->Fill(hit->at(j).r,rate);
          h_xy[2]->Fill(hit->at(j).x,hit->at(j).y,rate);
          if(hit->at(j).pz>=0){
            h_ratePzG0[2]->Fill(hit->at(j).r,rate);
            h_xyPzG0[2]->Fill(hit->at(j).x,hit->at(j).y,rate);
          }else{
            h_ratePzL0[2]->Fill(hit->at(j).r,rate);
            h_xyPzL0[2]->Fill(hit->at(j).x,hit->at(j).y,rate);
          }
        }else if(hit->at(j).pid==2112){
          ncounter++;
          h_rate[3]->Fill(hit->at(j).r,rate);
          h_xy[3]->Fill(hit->at(j).x,hit->at(j).y,rate);
          if(hit->at(j).pz>=0){
            h_ratePzG0[3]->Fill(hit->at(j).r,rate);
            h_xyPzG0[3]->Fill(hit->at(j).x,hit->at(j).y,rate);
          }else{
            h_ratePzL0[3]->Fill(hit->at(j).r,rate);
            h_xyPzL0[3]->Fill(hit->at(j).x,hit->at(j).y,rate);
          }
        }
      }
//        cout<<"Event: "<<ientry<<" nhits: "<<eMcounter<<"\t"<<ePcounter<<"\t"<<gcounter<<"\t"<<ncounter<<endl;
      if(eMcounter!=0)
        h_hits[0]->Fill(eMcounter);
      if(ePcounter!=0)
        h_hits[1]->Fill(ePcounter);
      if(gcounter!=0)
        h_hits[2]->Fill(gcounter);
      if(ncounter!=0)
        h_hits[3]->Fill(ncounter);

      eMcounter=0;
      ePcounter=0;
      gcounter=0;
      ncounter=0;
    }
    delete fin;
  }
  cout<<Form("Total number of file splits: %d",nfile)<<endl;
  cout<<Form("Total number of entries: %ld",nTotEv)<<endl;

  TFile* outfile = new TFile(Form("./rootfiles/%s_%s_%s.root",geometry.c_str(),tgt_gen_config.c_str(),plotType.c_str()),"recreate");
  for(int iSp=0;iSp<nSp;iSp++){
    h_hits[iSp]->Write();
    h_rate[iSp]->Scale(1.0/nTotEv);
    h_ratePzG0[iSp]->Scale(1.0/nTotEv);
    h_ratePzL0[iSp]->Scale(1.0/nTotEv);
    h_xy[iSp]->Scale(1.0/nTotEv);
    h_xyPzG0[iSp]->Scale(1.0/nTotEv);
    h_xyPzL0[iSp]->Scale(1.0/nTotEv);
    h_rate[iSp]->Write();
    h_ratePzG0[iSp]->Write();
    h_ratePzL0[iSp]->Write();
    h_xy[iSp]->Write();
    h_xyPzG0[iSp]->Write();
    h_xyPzL0[iSp]->Write();
  }

  TTree T("T","T");
  TBranch* nthrownEv = T.Branch("nthrownEv",&nTotEv,"nthrownEv/I");

  nthrownEv->Fill();
  T.SetEntries();
  outfile->cd();
  T.Write();
  outfile->Close();
}
