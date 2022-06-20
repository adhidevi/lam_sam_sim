#include "remolltypes.hh"
#include <sstream>
#include <iostream>
#include <fstream>
const double pi = TMath::Pi();
const double lam_length = 360.0;//azimuthal length of LAM quartz
const double lam_rin = 1010.0;//inner radius of LAM quartz
double lam_angle = atan(lam_length/lam_rin);
double sep_mid = 2*pi/14.0;

void radial_trans_radialCut_LAMopenSec(){
  gROOT->Reset();
  gStyle->SetOptStat(0);
  gStyle->SetTitleYOffset(1.3);
  gStyle->SetPadGridX(1);
  gStyle->SetPadGridY(1);
  TGaxis::SetMaxDigits(3);
 
  double kinEcut = 0;//MeV 
  const string spTit[] = {"e-/#pi-","e+/#pi+","#gamma","neutron","e-/e+","e-/#pi- vz<=-3875","e- trid==1"};
  const int nSp = sizeof(spTit)/sizeof(*spTit);
  const string spH[nSp] = {"epiM","epiP","g","n","ee","epiMvzCut","eTrIdCut"};
  map<int,int> spM {{11,1},{-211,1},{-11,2},{211,2},{22,3},{2112,4}};

///Change the following lines for which detectors you want to include////
  string detH[] = {"det174","det178","det28"};
  const int nDet = sizeof(detH)/sizeof(*detH);
  const int Det[nDet] = {174,178,28};
  map<int,int> dtM {{174,1},{178,2},{28,3}};
////////////////////////////////////////////////////////////////////////

  double x_min = 0;
  double x_max = 1900;
  const int nbin = 500;
  double bin_width = (x_max-x_min)/nbin;
  const Int_t SectorTit[] = {5,6,7,1,2,3,4};
  const int nSector = 7;
// if beam generator use the following
  const string unit ="hits/#thrownEvents";

  TH1F* h_rate[nSp][nDet][nSector];
  TH1F* h_ratePzG0[nSp][nDet][nSector];
  TH1F* h_ratePzL0[nSp][nDet][nSector];
  TH2F* h_xy[nSp][nDet][nSector];
  TH2F* h_xyPzG0[nSp][nDet][nSector];
  TH2F* h_xyPzL0[nSp][nDet][nSector];

///Change the following lines as needed////
  const string geometry = "develop_new";
  const string tgt_gen_config = "LH2_beam_V1";
  const string plotType = Form("radial_trans_EG%.0f_openSector",kinEcut);
  int beamGen(1);
//////////////////////////////////////////

  TFile* outfile = new TFile(Form("./rootfiles/%s_%s_%s.root",geometry.c_str(),tgt_gen_config.c_str(),plotType.c_str()),"recreate");
///Change this line for appropriate rootfile directory////
  TString rootfile_dir = "/volatile/halla/parity/adhidevi/remoll_rootfiles/develop_br";
//////////////////////////////////////////////////////////

  for(int iSp=0;iSp<nSp;iSp++){
   for(int iDet=0;iDet<nDet;iDet++){
    for(int iSector=0;iSector<nSector;iSector++){
      string title1D = Form("%s (KE>=%.0f MeV) Radial dist. on %s plane sector %d (%s);Radius (mm);%s/%.1fmm",spTit[iSp].c_str(),kinEcut,detH[iDet].c_str(),SectorTit[iSector],tgt_gen_config.c_str(),unit.c_str(),bin_width);
      h_rate[iSp][iDet][iSector] = new TH1F(Form("%s_r_%s_rate%d",detH[iDet].c_str(),spH[iSp].c_str(),SectorTit[iSector]),title1D.c_str(),nbin,x_min,x_max);
      h_ratePzG0[iSp][iDet][iSector] = new TH1F(Form("%s_rPzG0_%s_rate%d",detH[iDet].c_str(),spH[iSp].c_str(),SectorTit[iSector]),title1D.c_str(),nbin,x_min,x_max);
      h_ratePzL0[iSp][iDet][iSector] = new TH1F(Form("%s_rPzL0_%s_rate%d",detH[iDet].c_str(),spH[iSp].c_str(),SectorTit[iSector]),title1D.c_str(),nbin,x_min,x_max);

      string title2D = Form("%s (KE>=%.0f MeV) XY dist. on %s plane sector %d (%s);x (mm);y (mm)",spTit[iSp].c_str(),kinEcut,detH[iDet].c_str(),SectorTit[iSector],tgt_gen_config.c_str());
      h_xy[iSp][iDet][iSector] = new TH2F(Form("%s_xy_%s_rate%d",detH[iDet].c_str(),spH[iSp].c_str(),SectorTit[iSector]),title2D.c_str(),nbin,-x_max,x_max,nbin,-x_max,x_max);
      h_xyPzG0[iSp][iDet][iSector] = new TH2F(Form("%s_xyPzG0_%s_rate%d",detH[iDet].c_str(),spH[iSp].c_str(),SectorTit[iSector]),title2D.c_str(),nbin,-x_max,x_max,nbin,-x_max,x_max);
      h_xyPzL0[iSp][iDet][iSector] = new TH2F(Form("%s_xyPzL0_%s_rate%d",detH[iDet].c_str(),spH[iSp].c_str(),SectorTit[iSector]),title2D.c_str(),nbin,-x_max,x_max,nbin,-x_max,x_max);

      h_rate[iSp][iDet][iSector]->Sumw2();
      h_ratePzG0[iSp][iDet][iSector]->Sumw2();
      h_ratePzL0[iSp][iDet][iSector]->Sumw2();
    }
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
//comment following line if want to plot all r
//        if(hit->at(j).r>100) continue;
        if(hit->at(j).k<kinEcut) continue;
        
        double phi = hit->at(j).ph;
        if(phi<0) phi+=2.0*pi;
        for(int iSector=0;iSector<nSector;iSector++){
          if(!((phi>=(2*iSector+1)*sep_mid-lam_angle/2.0 && phi<=(2*iSector+1)*sep_mid+lam_angle/2.0))) continue;

          h_rate[sp][dt][iSector]->Fill(hit->at(j).r,rate);
          h_xy[sp][dt][iSector]->Fill(hit->at(j).x,hit->at(j).y,rate);

          if(hit->at(j).pz>=0){
            h_ratePzG0[sp][dt][iSector]->Fill(hit->at(j).r,rate);
            h_xyPzG0[sp][dt][iSector]->Fill(hit->at(j).x,hit->at(j).y,rate);
          }else{
            h_ratePzL0[sp][dt][iSector]->Fill(hit->at(j).r,rate);
            h_xyPzL0[sp][dt][iSector]->Fill(hit->at(j).x,hit->at(j).y,rate);
          }

          if(hit->at(j).pid==11 || hit->at(j).pid==-11){
            h_rate[4][dt][iSector]->Fill(hit->at(j).r,rate);
            h_xy[4][dt][iSector]->Fill(hit->at(j).x,hit->at(j).y,rate);
            if(hit->at(j).pz>=0){
              h_ratePzG0[4][dt][iSector]->Fill(hit->at(j).r,rate);
              h_xyPzG0[4][dt][iSector]->Fill(hit->at(j).x,hit->at(j).y,rate);
            }else{
              h_ratePzL0[4][dt][iSector]->Fill(hit->at(j).r,rate);
              h_xyPzL0[4][dt][iSector]->Fill(hit->at(j).x,hit->at(j).y,rate);
            }
          }

          if(hit->at(j).vz<=-3875 && (hit->at(j).pid==11 || hit->at(j).pid==-211)){
            h_rate[5][dt][iSector]->Fill(hit->at(j).r,rate);
            h_xy[5][dt][iSector]->Fill(hit->at(j).x,hit->at(j).y,rate);
            if(hit->at(j).pz>=0){
              h_ratePzG0[5][dt][iSector]->Fill(hit->at(j).r,rate);
              h_xyPzG0[5][dt][iSector]->Fill(hit->at(j).x,hit->at(j).y,rate);
            }else{
              h_ratePzL0[5][dt][iSector]->Fill(hit->at(j).r,rate);
              h_xyPzG0[5][dt][iSector]->Fill(hit->at(j).x,hit->at(j).y,rate);
            }
          }

          if(hit->at(j).trid==1 && hit->at(j).pid==11){
            h_rate[6][dt][iSector]->Fill(hit->at(j).r,rate);
            h_xy[6][dt][iSector]->Fill(hit->at(j).x,hit->at(j).y,rate);
            if(hit->at(j).pz>=0){
              h_ratePzG0[6][dt][iSector]->Fill(hit->at(j).r,rate);
              h_xyPzG0[6][dt][iSector]->Fill(hit->at(j).x,hit->at(j).y,rate);
            }else{
              h_ratePzL0[6][dt][iSector]->Fill(hit->at(j).r,rate);
              h_xyPzG0[6][dt][iSector]->Fill(hit->at(j).x,hit->at(j).y,rate);
            }
          }
        }

      }
    }
    delete T;
    delete fin;
  }

  cout<<Form("Total number of file splits: %d",nfile)<<endl;
  cout<<Form("Total number of entries: %ld",nTotEv)<<endl;
   
  for(int iDet=0;iDet<nDet;iDet++){
    outfile->mkdir(Form("%s",detH[iDet].c_str()));
    outfile->cd(Form("%s",detH[iDet].c_str()));
    for(int iSp=0;iSp<nSp;iSp++){
      for(int iSector=0;iSector<nSector;iSector++){
       if(beamGen){
         h_rate[iSp][iDet][iSector]->Scale(1.0/nTotEv);
         h_ratePzG0[iSp][iDet][iSector]->Scale(1.0/nTotEv);
         h_ratePzL0[iSp][iDet][iSector]->Scale(1.0/nTotEv);
         h_xy[iSp][iDet][iSector]->Scale(1.0/nTotEv);
         h_xyPzG0[iSp][iDet][iSector]->Scale(1.0/nTotEv);
         h_xyPzL0[iSp][iDet][iSector]->Scale(1.0/nTotEv);
       }
       h_rate[iSp][iDet][iSector]->Write();
       h_ratePzG0[iSp][iDet][iSector]->Write();
       h_ratePzL0[iSp][iDet][iSector]->Write();
       h_xy[iSp][iDet][iSector]->Write();
       h_xyPzG0[iSp][iDet][iSector]->Write();
       h_xyPzL0[iSp][iDet][iSector]->Write();
      }
    }
  }
}
