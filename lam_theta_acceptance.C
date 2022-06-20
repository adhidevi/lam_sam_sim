#include "remolltypes.hh"
#include <sstream>
#include <iostream>
#include <fstream>
const double pi = TMath::Pi();
const double lam_length = 360;//azimuthal length of LAM quartz
const double lam_rin = 1010.0;//inner radius of LAM quartz
double lam_phi_angle = atan(lam_length/lam_rin);
double sep_mid = 2*pi/14.0;

void lam_theta_acceptance(){
  gROOT->Reset();
  gStyle->SetOptStat(0);
  gStyle->SetTitleYOffset(1.3);
  gStyle->SetPadGridX(1);
  gStyle->SetPadGridY(1);
  TGaxis::SetMaxDigits(3);

  double kinEcut = 0;//MeV
  string spTit[] = {"e-/#pi-","e+/#pi+","#gamma","neutron","e-/#pi- vz<=-3875","e- trid==1"};
  const int nSp = sizeof(spTit)/sizeof(*spTit);
  for(int iSp=0;iSp<nSp;iSp++){
    spTit[iSp] = Form("%s (KE>=%.0f MeV) ",spTit[iSp].c_str(),kinEcut);
  }
  const string spH[nSp] = {"epiM","epiP","g","n","epiMvzCut","eTrIdCut"};
  map<int,int> spM {{11,1},{-211,1},{-11,2},{211,2},{22,3},{2112,4}};

  string detH[] = {"det174","det2174"};
  const int nDet = sizeof(detH)/sizeof(*detH);
  const int Det[nDet] = {174,2174};
  map<int,int> dtM {{174,1},{2174,2}};

  const double x_min = 0.0;
  const double x_max = 1900;
  const double r_min[nDet] = {1010,0};
  const double r_max[nDet] = {1130,1900};
  const int nbin = 500;
  const string weight = "rate";
  TH2F* h_xy[nSp][nDet];
  TH2F* h_xyPzG0[nSp][nDet];
  TH2F* h_xyPzL0[nSp][nDet];
  TH1F* h_r[nSp][nDet];
  TH1F* h_rPzG0[nSp][nDet];
  TH1F* h_rPzL0[nSp][nDet];
  TH1F* h_th[nSp][nDet];
  TH1F* h_thPzG0[nSp][nDet];
  TH1F* h_thPzL0[nSp][nDet];

  const string geometry = "develop";
  const string tgt_gen_config = "LH2_beam_V18";
  const string plotType = Form("theta_distribution_EG%.0f",kinEcut);
  int beamGen(1);

  TString rootfile_dir = "/volatile/halla/parity/adhidevi/remoll_rootfiles/develop_br";

  for(int iSp=0;iSp<nSp;iSp++){
   for(int iDet=0;iDet<nDet;iDet++){
     string titleXY = Form("%s Transverse dist. on %s plane (%s),%.0f<r<%.0f;x (mm);y (mm);hits/#thrownEvents",spTit[iSp].c_str(),detH[iDet].c_str(),tgt_gen_config.c_str(),r_min[iDet],r_max[iDet]);
     h_xy[iSp][iDet] = new TH2F(Form("%s_xy_%s_%s",detH[iDet].c_str(),spH[iSp].c_str(),weight.c_str()),titleXY.c_str(),nbin,-x_max,x_max,nbin,-x_max,x_max);
     h_xyPzG0[iSp][iDet] = new TH2F(Form("%s_xyPzG0_%s_%s",detH[iDet].c_str(),spH[iSp].c_str(),weight.c_str()),titleXY.c_str(),nbin,-x_max,x_max,nbin,-x_max,x_max);
     h_xyPzL0[iSp][iDet] = new TH2F(Form("%s_xyPzL0_%s_%s",detH[iDet].c_str(),spH[iSp].c_str(),weight.c_str()),titleXY.c_str(),nbin,-x_max,x_max,nbin,-x_max,x_max);

     string titleR = Form("%s Radial dist. on %s plane (%s),%.0f<r<%.0f;r (mm);hits/#thrownEvents",spTit[iSp].c_str(),detH[iDet].c_str(),tgt_gen_config.c_str(),r_min[iDet],r_max[iDet]);
     h_r[iSp][iDet] = new TH1F(Form("%s_r_%s_%s",detH[iDet].c_str(),spH[iSp].c_str(),weight.c_str()),titleR.c_str(),nbin,x_min,x_max);
     h_rPzG0[iSp][iDet] = new TH1F(Form("%s_rPzG0_%s_%s",detH[iDet].c_str(),spH[iSp].c_str(),weight.c_str()),titleR.c_str(),nbin,x_min,x_max);
     h_rPzL0[iSp][iDet] = new TH1F(Form("%s_rPzL0_%s_%s",detH[iDet].c_str(),spH[iSp].c_str(),weight.c_str()),titleR.c_str(),nbin,x_min,x_max);

     string titleTh = Form("%s Theta dist. on %s plane (%s),%.0f<r<%.0f;Theta (deg);hits/#thrownEvents",spTit[iSp].c_str(),detH[iDet].c_str(),tgt_gen_config.c_str(),r_min[iDet],r_max[iDet]);
     h_th[iSp][iDet] = new TH1F(Form("%s_th_%s_%s",detH[iDet].c_str(),spH[iSp].c_str(),weight.c_str()),titleTh.c_str(),200,0,15);
     h_thPzG0[iSp][iDet] = new TH1F(Form("%s_thPzG0_%s_%s",detH[iDet].c_str(),spH[iSp].c_str(),weight.c_str()),titleTh.c_str(),200,0,15);
     h_thPzL0[iSp][iDet] = new TH1F(Form("%s_thPzL0_%s_%s",detH[iDet].c_str(),spH[iSp].c_str(),weight.c_str()),titleTh.c_str(),200,0,15);
   }
  }

  int nfile=0;
  Long64_t nentry=0;
  long nTotEv=0;
  for(int ifile=1001;ifile<=6000;ifile++){
    string infile = Form("%s/%s/%s_%d.root",rootfile_dir.Data(),tgt_gen_config.c_str(),tgt_gen_config.c_str(),ifile);
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
        if(hit->at(j).r<r_min[dt] || hit->at(j).r>r_max[dt]) continue;
        if(hit->at(j).k<kinEcut) continue;

        double phi = hit->at(j).ph;
        if(!((abs(phi)>=1*sep_mid-lam_phi_angle/2.0 && abs(phi)<=1*sep_mid+lam_phi_angle/2.0) ||
            (abs(phi)>=3*sep_mid-lam_phi_angle/2.0 && abs(phi)<=3*sep_mid+lam_phi_angle/2.0) ||
            (abs(phi)>=5*sep_mid-lam_phi_angle/2.0 && abs(phi)<=5*sep_mid+lam_phi_angle/2.0) ||
            (abs(phi)>=7*sep_mid-lam_phi_angle/2.0 && abs(phi)<=7*sep_mid+lam_phi_angle/2.0) ||
            (abs(phi)>=9*sep_mid-lam_phi_angle/2.0 && abs(phi)<=9*sep_mid+lam_phi_angle/2.0) ||
            (abs(phi)>=11*sep_mid-lam_phi_angle/2.0 && abs(phi)<=11*sep_mid+lam_phi_angle/2.0) ||
            (abs(phi)>=13*sep_mid-lam_phi_angle/2.0 && abs(phi)<=13*sep_mid+lam_phi_angle/2.0))
          ) continue;

        double px = hit->at(j).px;
        double py = hit->at(j).py;
        double pz = hit->at(j).pz;
        double theta = (180.0/pi)*atan2(sqrt(px*px+py*py),pz);

        h_xy[sp][dt]->Fill(hit->at(j).x,hit->at(j).y,rate);
        h_r[sp][dt]->Fill(hit->at(j).r,rate);
        h_th[sp][dt]->Fill(theta,rate);
        if(hit->at(j).pz>=0){
          h_xyPzG0[sp][dt]->Fill(hit->at(j).x,hit->at(j).y,rate);
          h_rPzG0[sp][dt]->Fill(hit->at(j).r,rate);
          h_thPzG0[sp][dt]->Fill(theta,rate);
        }else{
          h_xyPzL0[sp][dt]->Fill(hit->at(j).x,hit->at(j).y,rate);
          h_rPzL0[sp][dt]->Fill(hit->at(j).r,rate);
          h_thPzL0[sp][dt]->Fill(theta,rate);
        }
        if(hit->at(j).vz<=-3875 && (hit->at(j).pid==11 || hit->at(j).pid==-211)){
          h_xy[4][dt]->Fill(hit->at(j).x,hit->at(j).y,rate);
          h_r[4][dt]->Fill(hit->at(j).r,rate);
          h_th[4][dt]->Fill(theta,rate);
          if(hit->at(j).pz>=0){
            h_xyPzG0[4][dt]->Fill(hit->at(j).x,hit->at(j).y,rate);
            h_rPzG0[4][dt]->Fill(hit->at(j).r,rate);
            h_thPzG0[4][dt]->Fill(theta,rate);
          }else{
            h_xyPzL0[4][dt]->Fill(hit->at(j).x,hit->at(j).y,rate);
            h_rPzL0[4][dt]->Fill(hit->at(j).r,rate);
            h_thPzL0[4][dt]->Fill(theta,rate);
          }
        }
        if(hit->at(j).trid==1 && hit->at(j).pid==11){
          h_xy[5][dt]->Fill(hit->at(j).x,hit->at(j).y,rate);
          h_r[5][dt]->Fill(hit->at(j).r,rate);
          h_th[5][dt]->Fill(theta,rate);
          if(hit->at(j).pz>=0){
            h_xyPzG0[5][dt]->Fill(hit->at(j).x,hit->at(j).y,rate);
            h_rPzG0[5][dt]->Fill(hit->at(j).r,rate);
            h_thPzG0[5][dt]->Fill(theta,rate);
          }else{
            h_xyPzG0[5][dt]->Fill(hit->at(j).x,hit->at(j).y,rate);
            h_rPzG0[5][dt]->Fill(hit->at(j).r,rate);
            h_thPzL0[5][dt]->Fill(theta,rate);
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
  TFile* outfile = new TFile(Form("./rootfiles/%s_%s_%s_rCut.root",geometry.c_str(),tgt_gen_config.c_str(),plotType.c_str()),"recreate");
  for(int iDet=0;iDet<nDet;iDet++){
    outfile->mkdir(Form("%s",detH[iDet].c_str()));
    outfile->cd(Form("%s",detH[iDet].c_str()));
    for(int iSp=0;iSp<nSp;iSp++){
       h_xy[iSp][iDet]->Scale(1.0/nTotEv);
       h_xyPzG0[iSp][iDet]->Scale(1.0/nTotEv);
       h_xyPzL0[iSp][iDet]->Scale(1.0/nTotEv);
       h_r[iSp][iDet]->Scale(1.0/nTotEv);
       h_rPzG0[iSp][iDet]->Scale(1.0/nTotEv);
       h_rPzL0[iSp][iDet]->Scale(1.0/nTotEv);
       h_th[iSp][iDet]->Scale(1.0/nTotEv);
       h_thPzG0[iSp][iDet]->Scale(1.0/nTotEv);
       h_thPzL0[iSp][iDet]->Scale(1.0/nTotEv);
       h_xy[iSp][iDet]->Write();
       h_xyPzG0[iSp][iDet]->Write();
       h_xyPzL0[iSp][iDet]->Write();
       h_r[iSp][iDet]->Write();
       h_rPzG0[iSp][iDet]->Write();
       h_rPzL0[iSp][iDet]->Write();
       h_th[iSp][iDet]->Write();
       h_thPzG0[iSp][iDet]->Write();
       h_thPzL0[iSp][iDet]->Write();
    }
  }
  outfile->Close();
}
