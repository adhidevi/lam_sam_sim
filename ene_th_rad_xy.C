#include "remolltypes.hh"
#include <sstream>
#include <iostream>
#include <fstream>

const double pi = TMath::Pi();
TString rootfile_dir = "/volatile/halla/parity/adhidevi/remoll_rootfiles/develop_br";
const string spTit[] = {"e-/#pi-","e+/#pi+","#gamma","neutron","e-/e+","e- trid==1"};
const int nSp = sizeof(spTit)/sizeof(*spTit);
const string spH[nSp] ={"epiM","epiP","g","n","ee","eTrIdCut"};
map<int,int> spM {{11,1},{-211,1},{-11,2},{211,2},{22,3},{2112,4}};
const string detH[] = {"det178","det28","det176"};
const int nDet = sizeof(detH)/sizeof(*detH);
map<int,int> dtM {{178,1},{28,2},{176,3}};
int Det[nDet] = {178,28,176};

const double kinEcut[] = {0,1};//MeV
const int nkinEcut = sizeof(kinEcut)/sizeof(*kinEcut);
const double rInCut[] = {1400,1500};
const int nRcut = sizeof(rInCut)/sizeof(*rInCut);
const double rOutCut[] = {1500,1700};

TH1D* hE[nSp][nDet][nRcut];
TH1D* hE_pzG0[nSp][nDet][nRcut];
TH1D* hE_pzL0[nSp][nDet][nRcut];
TH1D* hTh[nSp][nDet][nRcut];
TH1D* hTh_pzG0[nSp][nDet][nRcut];
TH1D* hTh_pzL0[nSp][nDet][nRcut];

TH1D* hR[nSp][nDet][nkinEcut];
TH1D* hR_pzG0[nSp][nDet][nkinEcut];
TH1D* hR_pzL0[nSp][nDet][nkinEcut];
TH2D* hXY[nSp][nDet][nkinEcut];
TH2D* hXY_pzG0[nSp][nDet][nkinEcut];
TH2D* hXY_pzL0[nSp][nDet][nkinEcut];

void niceLogBins(TH1*);
TFile* outfile;
int beamGen(1);
const string tgt_gen_config = "LH2_beam_V24";
string histE_title;
string histR_title;
string histXY_title;
string histTh_title;

void ene_th_rad_xy(){
   gStyle->SetOptStat(0);
   for(int iSp=0;iSp<nSp;iSp++){
    for(int iDet=0;iDet<nDet;iDet++){
     for(int ikinEcut=0;ikinEcut<nkinEcut;ikinEcut++){
      if(beamGen){
      histR_title = Form("%s (KE>%.0f MeV) Radial dist. on %s (%s);Radius (mm);hits/%sthrownEvents",spTit[iSp].c_str(),kinEcut[ikinEcut],detH[iDet].c_str(),tgt_gen_config.c_str(),"#");
      histXY_title = Form("%s (KE>%.0f MeV) XY dist. on %s (%s);x (mm);y (mm);hits/%sthrownEvents",spTit[iSp].c_str(),kinEcut[ikinEcut],detH[iDet].c_str(),tgt_gen_config.c_str(),"#");
      }else{
      histR_title = Form("%s (KE>%.0f MeV) Radial dist. on %s (%s);Radius (mm);rate (GHz)",spTit[iSp].c_str(),kinEcut[ikinEcut],detH[iDet].c_str(),tgt_gen_config.c_str());
      histXY_title = Form("%s (KE>%.0f MeV) XY dist. on %s (%s);x (mm);y (mm);rate (GHz)",spTit[iSp].c_str(),kinEcut[ikinEcut],detH[iDet].c_str(),tgt_gen_config.c_str());
      }
      hR[iSp][iDet][ikinEcut] = new TH1D(Form("%s_R_%s_%.0fMeV",detH[iDet].c_str(),spH[iSp].c_str(),kinEcut[ikinEcut]),histR_title.c_str(),500,0,1900);
      hR_pzG0[iSp][iDet][ikinEcut] = new TH1D(Form("%s_RpzG0_%s_%.0fMeV",detH[iDet].c_str(),spH[iSp].c_str(),kinEcut[ikinEcut]),histR_title.c_str(),500,0,1900);
      hR_pzL0[iSp][iDet][ikinEcut] = new TH1D(Form("%s_RpzL0_%s_%.0fMeV",detH[iDet].c_str(),spH[iSp].c_str(),kinEcut[ikinEcut]),histR_title.c_str(),500,0,1900);
      hXY[iSp][iDet][ikinEcut] = new TH2D(Form("%s_XY_%s_%.0fMeV",detH[iDet].c_str(),spH[iSp].c_str(),kinEcut[ikinEcut]),histXY_title.c_str(),1000,-1900,1900,1000,-1900,1900);
      hXY_pzG0[iSp][iDet][ikinEcut] = new TH2D(Form("%s_XYpzG0_%s_%.0fMeV",detH[iDet].c_str(),spH[iSp].c_str(),kinEcut[ikinEcut]),histXY_title.c_str(),1000,-1900,1900,1000,-1900,1900);
      hXY_pzL0[iSp][iDet][ikinEcut] = new TH2D(Form("%s_XYpzL0_%s_%.0fMeV",detH[iDet].c_str(),spH[iSp].c_str(),kinEcut[ikinEcut]),histXY_title.c_str(),1000,-1900,1900,1000,-1900,1900);

      hR[iSp][iDet][ikinEcut]->Sumw2();
      hR_pzG0[iSp][iDet][ikinEcut]->Sumw2();
      hR_pzL0[iSp][iDet][ikinEcut]->Sumw2();
     }
     for(int iRcut=0;iRcut<nRcut;iRcut++){
      if(beamGen){
      histE_title = Form("%s (%.0f<=r<=%.0f) KE dist. on %s (%s);E (MeV);hits/%sthrownEvents",spTit[iSp].c_str(),rInCut[iRcut],rOutCut[iRcut],detH[iDet].c_str(),tgt_gen_config.c_str(),"#");
      histTh_title = Form("%s (%.0f<=r<=%.0f, KE>1 MeV) Theta dist. on %s (%s);Theta (deg);hits/%sthrownEvents",spTit[iSp].c_str(),rInCut[iRcut],rOutCut[iRcut],detH[iDet].c_str(),tgt_gen_config.c_str(),"#");
      }else{
      histE_title = Form("%s (%.0f<=r<=%.0f) KE dist. on %s (%s);E (MeV);rate (GHz)",spTit[iSp].c_str(),rInCut[iRcut],rOutCut[iRcut],detH[iDet].c_str(),tgt_gen_config.c_str());
      histTh_title = Form("%s (%.0f<=r<=%.0f, KE>1 MeV) Theta dist. on %s (%s);Theta (deg);rate (GHz)",spTit[iSp].c_str(),rInCut[iRcut],rOutCut[iRcut],detH[iDet].c_str(),tgt_gen_config.c_str());
      }
      hE[iSp][iDet][iRcut] = new TH1D(Form("%s_E_%s_r%d",detH[iDet].c_str(),spH[iSp].c_str(),iRcut),histE_title.c_str(),121,-8,5);
      hE_pzG0[iSp][iDet][iRcut] = new TH1D(Form("%s_EpzG0_%s_r%d",detH[iDet].c_str(),spH[iSp].c_str(),iRcut),histE_title.c_str(),121,-8,5);
      hE_pzL0[iSp][iDet][iRcut] = new TH1D(Form("%s_EpzL0_%s_r%d",detH[iDet].c_str(),spH[iSp].c_str(),iRcut),histE_title.c_str(),121,-8,5);
      hTh[iSp][iDet][iRcut] = new TH1D(Form("%s_Th_%s_r%d",detH[iDet].c_str(),spH[iSp].c_str(),iRcut),histTh_title.c_str(),200,0,90);
      hTh_pzG0[iSp][iDet][iRcut] = new TH1D(Form("%s_ThpzG0_%s_r%d",detH[iDet].c_str(),spH[iSp].c_str(),iRcut),histTh_title.c_str(),200,0,90);
      hTh_pzL0[iSp][iDet][iRcut] = new TH1D(Form("%s_ThpzL0_%s_r%d",detH[iDet].c_str(),spH[iSp].c_str(),iRcut),histTh_title.c_str(),200,0,90);

      hE[iSp][iDet][iRcut]->Sumw2();
      hE_pzG0[iSp][iDet][iRcut]->Sumw2();
      hE_pzL0[iSp][iDet][iRcut]->Sumw2();
      niceLogBins(hE[iSp][iDet][iRcut]);
      niceLogBins(hE_pzG0[iSp][iDet][iRcut]);
      niceLogBins(hE_pzL0[iSp][iDet][iRcut]);
      hTh[iSp][iDet][iRcut]->Sumw2();
      hTh_pzG0[iSp][iDet][iRcut]->Sumw2();
      hTh_pzL0[iSp][iDet][iRcut]->Sumw2();
     }
    }
   }

   int nfile=0;
   Long64_t nentry=0;
   long nTotEv=0;
   for(int ifile=3001;ifile<=4000;ifile++){
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

      std::vector<remollGenericDetectorHit_t> *hit =0;
      std::vector<remollEventParticle_t> *part =0;
      remollBeamTarget_t *bm =0;
      remollEvent_t *ev =0;
      Double_t rate =0;
      Double_t phi = 0;
      Double_t modphi = 0;
      T->SetBranchAddress("hit", & hit);
      T->SetBranchAddress("part", & part);
      T->SetBranchAddress("bm", & bm);
      T->SetBranchAddress("ev", & ev);
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
          phi = hit->at(j).ph;
          if(phi<0) phi +=2.0*TMath::Pi();
          modphi = fmod(phi,2.0*TMath::Pi()/7.0);
//comment following line if want to plot full azimuth
//          if(modphi<3.0*TMath::Pi()/28.0 || modphi>5.0*TMath::Pi()/28.0) continue;

          double px = hit->at(j).px;
          double py = hit->at(j).py;
          double pz = hit->at(j).pz;
          double theta = (180.0/pi)*atan2(sqrt(px*px+py*py),pz);

          for(int ikinEcut=0;ikinEcut<nkinEcut;ikinEcut++){
           if(hit->at(j).e>kinEcut[ikinEcut]){
            hR[sp][dt][ikinEcut]->Fill(hit->at(j).r,rate);
            hXY[sp][dt][ikinEcut]->Fill(hit->at(j).x,hit->at(j).y,rate);
            if(hit->at(j).pz>=0){
             hR_pzG0[sp][dt][ikinEcut]->Fill(hit->at(j).r,rate);
             hXY_pzG0[sp][dt][ikinEcut]->Fill(hit->at(j).x,hit->at(j).y,rate);
            }else{
             hR_pzL0[sp][dt][ikinEcut]->Fill(hit->at(j).r,rate);
             hXY_pzL0[sp][dt][ikinEcut]->Fill(hit->at(j).x,hit->at(j).y,rate);
            }
            if(hit->at(j).pid==11 || hit->at(j).pid==-11){
             hR[4][dt][ikinEcut]->Fill(hit->at(j).r,rate);
             hXY[4][dt][ikinEcut]->Fill(hit->at(j).x,hit->at(j).y,rate);
             if(hit->at(j).pz>=0){
              hR_pzG0[4][dt][ikinEcut]->Fill(hit->at(j).r,rate);
              hXY_pzG0[4][dt][ikinEcut]->Fill(hit->at(j).x,hit->at(j).y,rate);
             }else{
              hR_pzL0[4][dt][ikinEcut]->Fill(hit->at(j).r,rate);
              hXY_pzL0[4][dt][ikinEcut]->Fill(hit->at(j).x,hit->at(j).y,rate);
             }
            }
            if(hit->at(j).trid==1 && hit->at(j).pid==11){
              hR[5][dt][ikinEcut]->Fill(hit->at(j).r,rate);
              hXY[5][dt][ikinEcut]->Fill(hit->at(j).x,hit->at(j).y,rate);
              if(hit->at(j).pz>=0){
               hR_pzG0[5][dt][ikinEcut]->Fill(hit->at(j).r,rate);
               hXY_pzG0[5][dt][ikinEcut]->Fill(hit->at(j).x,hit->at(j).y,rate);
              }else{
               hR_pzL0[5][dt][ikinEcut]->Fill(hit->at(j).r,rate);
               hXY_pzL0[5][dt][ikinEcut]->Fill(hit->at(j).x,hit->at(j).y,rate);
              }
            }
           }
          }

          for(int iRcut=0;iRcut<nRcut;iRcut++){
           if(hit->at(j).r>=rInCut[iRcut] && hit->at(j).r<=rOutCut[iRcut]){
            hE[sp][dt][iRcut]->Fill(hit->at(j).k,rate);
            if(hit->at(j).pz>=0){
             hE_pzG0[sp][dt][iRcut]->Fill(hit->at(j).k,rate);
            }else{
             hE_pzL0[sp][dt][iRcut]->Fill(hit->at(j).k,rate);
            }
            if(hit->at(j).pid==11 || hit->at(j).pid==-11){
             hE[4][dt][iRcut]->Fill(hit->at(j).k,rate);
             if(hit->at(j).pz>=0){
              hE_pzG0[4][dt][iRcut]->Fill(hit->at(j).k,rate);
             }else{
              hE_pzL0[4][dt][iRcut]->Fill(hit->at(j).k,rate);
             }
            }
            if(hit->at(j).trid==1 && hit->at(j).pid==11){
              hE[5][dt][iRcut]->Fill(hit->at(j).k,rate);
              if(hit->at(j).pz>=0){
               hE_pzG0[5][dt][iRcut]->Fill(hit->at(j).k,rate);
              }else{
               hE_pzL0[5][dt][iRcut]->Fill(hit->at(j).k,rate);
              }
            }
            if(hit->at(j).k>1){
            hTh[sp][dt][iRcut]->Fill(theta,rate);
            if(hit->at(j).pz>=0){
             hTh_pzG0[sp][dt][iRcut]->Fill(theta,rate);
            }else{
             hTh_pzL0[sp][dt][iRcut]->Fill(theta,rate);
            }
            if(hit->at(j).pid==11 || hit->at(j).pid==-11){
             hTh[4][dt][iRcut]->Fill(theta,rate);
             if(hit->at(j).pz>=0){
              hTh_pzG0[4][dt][iRcut]->Fill(theta,rate);
             }else{
              hTh_pzL0[4][dt][iRcut]->Fill(theta,rate);
             }
            }
            if(hit->at(j).trid==1 && hit->at(j).pid==11){
              hTh[5][dt][iRcut]->Fill(theta,rate);
              if(hit->at(j).pz>=0){
               hTh_pzG0[5][dt][iRcut]->Fill(theta,rate);
              }else{
               hTh_pzL0[5][dt][iRcut]->Fill(theta,rate);
              }
             }
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
   outfile = TFile::Open(Form("./rootfiles/ene_th_rad_xy_DBM_%s_split3.root",tgt_gen_config.c_str()),"RECREATE");
   outfile->cd();
   for(int iDet=0;iDet<nDet;iDet++){
     outfile->mkdir(Form("%s",detH[iDet].c_str()));
     outfile->cd(Form("%s",detH[iDet].c_str()));
     for(int iSp=0;iSp<nSp;iSp++){
      for(int ikinEcut=0;ikinEcut<nkinEcut;ikinEcut++){
       if(beamGen){
         hR[iSp][iDet][ikinEcut]->Scale(1.0/nTotEv);
         hR_pzG0[iSp][iDet][ikinEcut]->Scale(1.0/nTotEv);
         hR_pzL0[iSp][iDet][ikinEcut]->Scale(1.0/nTotEv);
         hXY[iSp][iDet][ikinEcut]->Scale(1.0/nTotEv);
         hXY_pzG0[iSp][iDet][ikinEcut]->Scale(1.0/nTotEv);
         hXY_pzL0[iSp][iDet][ikinEcut]->Scale(1.0/nTotEv);
       }else{
         hR[iSp][iDet][ikinEcut]->Scale(1.0e-9/nfile);
         hR_pzG0[iSp][iDet][ikinEcut]->Scale(1.0e-9/nfile);
         hR_pzL0[iSp][iDet][ikinEcut]->Scale(1.0e-9/nfile);
         hXY[iSp][iDet][ikinEcut]->Scale(1.0e-9/nfile);
         hXY_pzG0[iSp][iDet][ikinEcut]->Scale(1.0e-9/nfile);
         hXY_pzL0[iSp][iDet][ikinEcut]->Scale(1.0e-9/nfile);
       }
       hR[iSp][iDet][ikinEcut]->Write();
       hR_pzG0[iSp][iDet][ikinEcut]->Write();
       hR_pzL0[iSp][iDet][ikinEcut]->Write();
       hXY[iSp][iDet][ikinEcut]->Write();
       hXY_pzG0[iSp][iDet][ikinEcut]->Write();
       hXY_pzL0[iSp][iDet][ikinEcut]->Write();
     }
      for(int iRcut=0;iRcut<nRcut;iRcut++){
       if(beamGen){
         hE[iSp][iDet][iRcut]->Scale(1.0/nTotEv);
         hE_pzG0[iSp][iDet][iRcut]->Scale(1.0/nTotEv);
         hE_pzL0[iSp][iDet][iRcut]->Scale(1.0/nTotEv);
         hTh[iSp][iDet][iRcut]->Scale(1.0/nTotEv);
         hTh_pzG0[iSp][iDet][iRcut]->Scale(1.0/nTotEv);
         hTh_pzL0[iSp][iDet][iRcut]->Scale(1.0/nTotEv);
       }else{
         hE[iSp][iDet][iRcut]->Scale(1.0e-9/nfile);
         hE_pzG0[iSp][iDet][iRcut]->Scale(1.0e-9/nfile);
         hE_pzL0[iSp][iDet][iRcut]->Scale(1.0e-9/nfile);
         hTh[iSp][iDet][iRcut]->Scale(1.0e-9/nfile);
         hTh_pzG0[iSp][iDet][iRcut]->Scale(1.0e-9/nfile);
         hTh_pzL0[iSp][iDet][iRcut]->Scale(1.0e-9/nfile);
       }
       hE[iSp][iDet][iRcut]->Write();
       hE_pzG0[iSp][iDet][iRcut]->Write();
       hE_pzL0[iSp][iDet][iRcut]->Write();
       hTh[iSp][iDet][iRcut]->Write();
       hTh_pzG0[iSp][iDet][iRcut]->Write();
       hTh_pzL0[iSp][iDet][iRcut]->Write();
     }
    }
   }
}

void niceLogBins(TH1*h)
{
  TAxis *axis = h->GetXaxis();
  int bins = axis->GetNbins();

  double from = axis->GetXmin();
  double to = axis->GetXmax();
  double width = (to - from) / bins;
  double *new_bins = new double[bins + 1];

  for (int i = 0; i <= bins; i++) {
    new_bins[i] = pow(10, from + i * width);
  }
  axis->Set(bins, new_bins);
  delete[] new_bins;
}
