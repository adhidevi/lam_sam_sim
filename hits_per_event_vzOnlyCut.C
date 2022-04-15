#include "remolltypes.hh"
#include <sstream>
#include <iostream>
#include <fstream>
const double pi = TMath::Pi();
const double lam_length = 360.0;//azimuthal length of LAM quartz
const double lam_rin = 1010.0;//inner radius of LAM quartz
double lam_angle = atan(lam_length/lam_rin);
double sep_mid = 2*pi/14.0;

void hits_per_event_vzOnlyCut(){
  gROOT->Reset();
  gStyle->SetOptStat(0);
  gStyle->SetTitleYOffset(1.3);
  gStyle->SetPadGridX(1);
  gStyle->SetPadGridY(1);
  TGaxis::SetMaxDigits(3);

  const int Det = 174;
  const double inRcut = 1010;
  const double outRcut = 1130;
  const double kinEcut = 6;//MeV
  const string spTit[] = {"e-/#pi-","e+/#pi+","#gamma","neutron"};
  const int nSp = sizeof(spTit)/sizeof(*spTit);
  const string spH[nSp] = {"epiM","epiP","g","n"};
  map<int,int> spM {{11,1},{-211,1},{-11,2},{211,2},{22,3},{2112,4}};

  TString rootfile_dir = "/volatile/halla/parity/adhidevi/remoll_rootfiles/develop_br";
  const string geometry = "develop";
  const string tgt_gen_config = "LH2_beam_V18";
  const string plotType = Form("hits_dist_vzOnlyCut_KEG%.0fMeV_split2",kinEcut);
  int beamGen(1);

  const int nbin_hits = 100;
  const double x_min = 0;
  const double x_max = 100;

  const int nbin_r = 500;
  const double r_min = 0;
  const double r_max = 1900;
  double bin_width = (r_max - r_min)/nbin_r;
  const string vzCut[] = {"no_vzCut","target_vzCut","flange_vzCut","collar2OR_vzCut"};
  const int nCut = sizeof(vzCut)/sizeof(*vzCut);

  TH1D* h_hits[nSp][nCut];
  TH1F* h_rate[nSp][nCut];
  TH1F* h_ratePzG0[nSp][nCut];
  TH1F* h_ratePzL0[nSp][nCut];
  TH2F* h_xy[nSp][nCut];
  TH2F* h_xyPzG0[nSp][nCut];
  TH2F* h_xyPzL0[nSp][nCut];

  for(int iSp=0;iSp<nSp;iSp++){
   for(int iCut=0;iCut<nCut;iCut++){
    h_hits[iSp][iCut] = new TH1D(Form("nhits_%s_%s",spH[iSp].c_str(),vzCut[iCut].c_str()),Form("%s (KE>%.0f MeV) nhits at det%d plane (%s);nhits;Events",spTit[iSp].c_str(),kinEcut,Det,tgt_gen_config.c_str()),nbin_hits,x_min,x_max);

    string title1D = Form("%s Radial dist. on det%d plane (%s);Radius (mm);hits/#thrownEvents/%.1fmm",spTit[iSp].c_str(),Det,tgt_gen_config.c_str(),bin_width);

    h_rate[iSp][iCut] = new TH1F(Form("det%d_r_%s_%s",Det,spH[iSp].c_str(),vzCut[iCut].c_str()),title1D.c_str(),nbin_r,r_min,r_max);
    h_ratePzG0[iSp][iCut] = new TH1F(Form("det%d_rPzG0_%s_%s",Det,spH[iSp].c_str(),vzCut[iCut].c_str()),title1D.c_str(),nbin_r,r_min,r_max);
    h_ratePzL0[iSp][iCut] = new TH1F(Form("det%d_rPzL0_%s_%s",Det,spH[iSp].c_str(),vzCut[iCut].c_str()),title1D.c_str(),nbin_r,r_min,r_max);

    string title2D = Form("%s XY dist. on det%d plane (%s);x (mm);y (mm)",spTit[iSp].c_str(),Det,tgt_gen_config.c_str());

    h_xy[iSp][iCut] = new TH2F(Form("det%d_xy_%s_%s",Det,spH[iSp].c_str(),vzCut[iCut].c_str()),title2D.c_str(),nbin_r,-r_max,r_max,nbin_r,-r_max,r_max);
    h_xyPzG0[iSp][iCut] = new TH2F(Form("det%d_xyPzG0_%s_%s",Det,spH[iSp].c_str(),vzCut[iCut].c_str()),title2D.c_str(),nbin_r,-r_max,r_max,nbin_r,-r_max,r_max);
    h_xyPzL0[iSp][iCut] = new TH2F(Form("det%d_xyPzL0_%s_%s",Det,spH[iSp].c_str(),vzCut[iCut].c_str()),title2D.c_str(),nbin_r,-r_max,r_max,nbin_r,-r_max,r_max);

    h_hits[iSp][iCut]->Sumw2();
    h_rate[iSp][iCut]->Sumw2();
    h_ratePzG0[iSp][iCut]->Sumw2();
    h_ratePzL0[iSp][iCut]->Sumw2();
   }
  }

  int nfile = 0;
  Long64_t nentry = 0;
  long nTotEv = 0;

  for(int ifile=3501;ifile<=6000;ifile++){
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

    if(T->GetEntries()<50000){
      cout<<Form("Skipping %s. Bad file.",infile.c_str())<<endl;
      continue;
    }

    nentry = T->GetEntries();
    cout<<Form("Found %lld entries in filesplit %d",nentry,nfile)<<endl;
    nTotEv+=nentry;

    std::vector<remollGenericDetectorHit_t> *hit = 0;
    Double_t rate = 0;
    T->SetBranchAddress("hit", &hit);
    T->SetBranchAddress("rate", &rate);

    for(int ientry=0;ientry<nentry;ientry++){
      T->GetEntry(ientry);
      int spCounter[nSp] = {0};
      int spCounterT[nSp] = {0};
      int spCounterF[nSp] = {0};
      int spCounterC[nSp] = {0};
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
        
        int sp = spM[int(hit->at(j).pid)]-1;
        if(sp==-1) continue;
        if(hit->at(j).k<kinEcut) continue;
        if(hit->at(j).det!=Det) continue;
        if(hit->at(j).r<=inRcut || hit->at(j).r>=outRcut) continue;

        spCounter[sp]++;
        h_rate[sp][0]->Fill(hit->at(j).r,rate);
        h_xy[sp][0]->Fill(hit->at(j).x,hit->at(j).y,rate);
        if(hit->at(j).pz>=0){
          h_ratePzG0[sp][0]->Fill(hit->at(j).r,rate);
          h_xyPzG0[sp][0]->Fill(hit->at(j).x,hit->at(j).y,rate);
        }else{
          h_ratePzL0[sp][0]->Fill(hit->at(j).r,rate);
          h_xyPzL0[sp][0]->Fill(hit->at(j).x,hit->at(j).y,rate);
        }
        
        if(hit->at(j).vz<=-3875){
          spCounterT[sp]++;
          h_rate[sp][1]->Fill(hit->at(j).r,rate);
          h_xy[sp][1]->Fill(hit->at(j).x,hit->at(j).y,rate);
          if(hit->at(j).pz>=0){
            h_ratePzG0[sp][1]->Fill(hit->at(j).r,rate);
            h_xyPzG0[sp][1]->Fill(hit->at(j).x,hit->at(j).y,rate);
          }else{
            h_ratePzL0[sp][1]->Fill(hit->at(j).r,rate);
            h_xyPzL0[sp][1]->Fill(hit->at(j).x,hit->at(j).y,rate);
          }
        }
 
        if(hit->at(j).vz>18850&&hit->at(j).vz<18925){
          spCounterF[sp]++;
          h_rate[sp][2]->Fill(hit->at(j).r,rate);
          h_xy[sp][2]->Fill(hit->at(j).x,hit->at(j).y,rate);
          if(hit->at(j).pz>=0){
            h_ratePzG0[sp][2]->Fill(hit->at(j).r,rate);
            h_xyPzG0[sp][2]->Fill(hit->at(j).x,hit->at(j).y,rate);
          }else{
            h_ratePzL0[sp][2]->Fill(hit->at(j).r,rate);
            h_xyPzL0[sp][2]->Fill(hit->at(j).x,hit->at(j).y,rate);
          }
        }
        if(hit->at(j).vz>18925&&hit->at(j).vz<19090){
          spCounterC[sp]++;
          h_rate[sp][3]->Fill(hit->at(j).r,rate);
          h_xy[sp][3]->Fill(hit->at(j).x,hit->at(j).y,rate);
          if(hit->at(j).pz>=0){
            h_ratePzG0[sp][3]->Fill(hit->at(j).r,rate);
            h_xyPzG0[sp][3]->Fill(hit->at(j).x,hit->at(j).y,rate);
          }else{
            h_ratePzL0[sp][3]->Fill(hit->at(j).r,rate);
            h_xyPzL0[sp][3]->Fill(hit->at(j).x,hit->at(j).y,rate);
          }
        }

      }

      for(int iSp=0;iSp<nSp;iSp++){
        if(spCounter[iSp]!=0){
          h_hits[iSp][0]->Fill(spCounter[iSp]);
          spCounter[iSp]=0;
        }
      }
      for(int iSp=0;iSp<nSp;iSp++){
        if(spCounterT[iSp]!=0){
          h_hits[iSp][1]->Fill(spCounterT[iSp]);
          spCounterT[iSp]=0;
        }
      }
      for(int iSp=0;iSp<nSp;iSp++){
        if(spCounterF[iSp]!=0){
          h_hits[iSp][2]->Fill(spCounterF[iSp]);
          spCounterF[iSp]=0;
        }
      }
      for(int iSp=0;iSp<nSp;iSp++){
        if(spCounterC[iSp]!=0){
          h_hits[iSp][3]->Fill(spCounterC[iSp]);
          spCounterC[iSp]=0;
        }
      }
    }
    delete T;
    delete fin;
  }
  cout<<Form("Total number of file splits: %d",nfile)<<endl;
  cout<<Form("Total number of entries: %ld",nTotEv)<<endl;

  TFile* outfile = new TFile(Form("./rootfiles/%s_%s_%s.root",geometry.c_str(),tgt_gen_config.c_str(),plotType.c_str()),"recreate");
  for(int iSp=0;iSp<nSp;iSp++){
   for(int iCut=0;iCut<nCut;iCut++){
    h_hits[iSp][iCut]->Write();
    h_rate[iSp][iCut]->Scale(1.0/nTotEv);
    h_ratePzG0[iSp][iCut]->Scale(1.0/nTotEv);
    h_ratePzL0[iSp][iCut]->Scale(1.0/nTotEv);
    h_xy[iSp][iCut]->Scale(1.0/nTotEv);
    h_xyPzG0[iSp][iCut]->Scale(1.0/nTotEv);
    h_xyPzL0[iSp][iCut]->Scale(1.0/nTotEv);
    h_rate[iSp][iCut]->Write();
    h_ratePzG0[iSp][iCut]->Write();
    h_ratePzL0[iSp][iCut]->Write();
    h_xy[iSp][iCut]->Write();
    h_xyPzG0[iSp][iCut]->Write();
    h_xyPzL0[iSp][iCut]->Write();
   }
  }

  TTree Tree("T","T");
  TBranch* nthrownEv = Tree.Branch("nthrownEv",&nTotEv,"nthrownEv/I");

  nthrownEv->Fill();
  Tree.SetEntries();
  outfile->cd();
  Tree.Write();
  outfile->Close();
}
