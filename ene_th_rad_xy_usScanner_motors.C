#include "remolltypes.hh"
#include <sstream>
#include <iostream>
#include <fstream>

const double pi = TMath::Pi();
TString rootfile_dir = "/volatile/halla/parity/adhidevi/remoll_rootfiles/develop_br";
const string spTit[] = {"e-/#pi-","e+/#pi+","#gamma","neutron","e-/#pi- vz<=-3875","e- trid==1"};
const int nSp = sizeof(spTit)/sizeof(*spTit);
const string spH[nSp] ={"epiM","epiP","g","n","epiMvzCut","eTrIdCut"};
map<int,int> spM {{11,1},{-211,1},{-11,2},{211,2},{22,3},{2112,4}};
const string detH[] = {"det175","det1751"};
const int nDet = sizeof(detH)/sizeof(*detH);
map<int,int> dtM {{175,1},{1751,2}};
int Det[nDet] = {175,1751};
const Int_t nBand = 3;
double rIn[nBand] = {2174.1060,2462.97,2174.1060};
double rOut[nBand] = {2277.7732,2517.2279,2517.2279};
TH1D* hE[nSp][nDet][nBand];
TH1D* hE_pzG0[nSp][nDet][nBand];
TH1D* hE_pzL0[nSp][nDet][nBand];
TH1D* hR[nSp][nDet][nBand];
TH1D* hR_pzG0[nSp][nDet][nBand];
TH1D* hR_pzL0[nSp][nDet][nBand];
TH2D* hXY[nSp][nDet][nBand];
TH2D* hXY_pzG0[nSp][nDet][nBand];
TH2D* hXY_pzL0[nSp][nDet][nBand];
TH1D* hTh[nSp][nDet][nBand];
TH1D* hTh_pzG0[nSp][nDet][nBand];
TH1D* hTh_pzL0[nSp][nDet][nBand];

void niceLogBins(TH1*);
TFile* outfile;
int beamGen(1);
const string tgt_gen_config = "LH2_beam_V23";
string histE_title;
string histR_title;
string histXY_title;
string histTh_title;
const double kinEcut = 0;//MeV

void ene_th_rad_xy_usScanner_motors(){
   gStyle->SetOptStat(0);
   for(int iSp=0;iSp<nSp;iSp++){
    for(int iDet=0;iDet<nDet;iDet++){
     for(int iBand=0;iBand<nBand;iBand++){
      if(beamGen){
      histE_title = Form("%s (KE>%.0f MeV) KE dist. on %s (%s),%.4f<=r<=%.4f;E (MeV);hits/%sthrownEvents",spTit[iSp].c_str(),kinEcut,detH[iDet].c_str(),tgt_gen_config.c_str(),rIn[iBand],rOut[iBand],"#");
      histR_title = Form("%s (KE>%.0f MeV) Radial dist. on %s (%s),%.4f<=r<=%.4f;Radius (mm);hits/%sthrownEvents",spTit[iSp].c_str(),kinEcut,detH[iDet].c_str(),tgt_gen_config.c_str(),rIn[iBand],rOut[iBand],"#");
      histXY_title = Form("%s (KE>%.0f MeV) XY dist. on %s (%s),%.4f<=r<=%.4f;x (mm);y (mm);hits/%sthrownEvents",spTit[iSp].c_str(),kinEcut,detH[iDet].c_str(),tgt_gen_config.c_str(),rIn[iBand],rOut[iBand],"#");
      histTh_title = Form("%s (KE>%.0f MeV) Theta dist. on %s (%s),%.4f<=r<=%.4f;Theta (deg);hits/%sthrownEvents",spTit[iSp].c_str(),kinEcut,detH[iDet].c_str(),tgt_gen_config.c_str(),rIn[iBand],rOut[iBand],"#");
      }else{
      histE_title = Form("%s (KE>%.0f MeV) KE dist. on %s (%s),%.4f<=r<=%.4f;E (MeV);rate (GHz)",spTit[iSp].c_str(),kinEcut,detH[iDet].c_str(),tgt_gen_config.c_str(),rIn[iBand],rOut[iBand]);
      histR_title = Form("%s (KE>%.0f MeV) Radial dist. on %s (%s),%.4f<=r<=%.4f;Radius (mm);rate (GHz)",spTit[iSp].c_str(),kinEcut,detH[iDet].c_str(),tgt_gen_config.c_str(),rIn[iBand],rOut[iBand]);
      histXY_title = Form("%s (KE>%.0f MeV) XY dist. on %s (%s),%.4f<=r<=%.4f;x (mm);y (mm);rate (GHz)",spTit[iSp].c_str(),kinEcut,detH[iDet].c_str(),tgt_gen_config.c_str(),rIn[iBand],rOut[iBand]);
      histTh_title = Form("%s (KE>%.0f MeV) Theta dist. on %s (%s),%.4f<=r<=%.4f;Theta (deg);rate (GHz)",spTit[iSp].c_str(),kinEcut,detH[iDet].c_str(),tgt_gen_config.c_str(),rIn[iBand],rOut[iBand]);
      }
      hE[iSp][iDet][iBand] = new TH1D(Form("%s_E_%s_%d",detH[iDet].c_str(),spH[iSp].c_str(),iBand),histE_title.c_str(),121,-8,5);
      hE_pzG0[iSp][iDet][iBand] = new TH1D(Form("%s_EpzG0_%s_%d",detH[iDet].c_str(),spH[iSp].c_str(),iBand),histE_title.c_str(),121,-8,5);
      hE_pzL0[iSp][iDet][iBand] = new TH1D(Form("%s_EpzL0_%s_%d",detH[iDet].c_str(),spH[iSp].c_str(),iBand),histE_title.c_str(),121,-8,5);
      hR[iSp][iDet][iBand] = new TH1D(Form("%s_R_%s_%d",detH[iDet].c_str(),spH[iSp].c_str(),iBand),histR_title.c_str(),1000,0,3000);
      hR_pzG0[iSp][iDet][iBand] = new TH1D(Form("%s_RpzG0_%s_%d",detH[iDet].c_str(),spH[iSp].c_str(),iBand),histR_title.c_str(),1000,0,3000);
      hR_pzL0[iSp][iDet][iBand] = new TH1D(Form("%s_RpzL0_%s_%d",detH[iDet].c_str(),spH[iSp].c_str(),iBand),histR_title.c_str(),1000,0,3000);
      hXY[iSp][iDet][iBand] = new TH2D(Form("%s_XY_%s_%d",detH[iDet].c_str(),spH[iSp].c_str(),iBand),histXY_title.c_str(),1000,-3000,3000,1000,-3000,3000);
      hXY_pzG0[iSp][iDet][iBand] = new TH2D(Form("%s_XYpzG0_%s_%d",detH[iDet].c_str(),spH[iSp].c_str(),iBand),histXY_title.c_str(),1000,-3000,3000,1000,-3000,3000);
      hXY_pzL0[iSp][iDet][iBand] = new TH2D(Form("%s_XYpzL0_%s_%d",detH[iDet].c_str(),spH[iSp].c_str(),iBand),histXY_title.c_str(),1000,-3000,3000,1000,-3000,3000);
      hTh[iSp][iDet][iBand] = new TH1D(Form("%s_Th_%s_%d",detH[iDet].c_str(),spH[iSp].c_str(),iBand),histTh_title.c_str(),200,0,90);
      hTh_pzG0[iSp][iDet][iBand] = new TH1D(Form("%s_ThpzG0_%s_%d",detH[iDet].c_str(),spH[iSp].c_str(),iBand),histTh_title.c_str(),200,0,90);
      hTh_pzL0[iSp][iDet][iBand] = new TH1D(Form("%s_ThpzL0_%s_%d",detH[iDet].c_str(),spH[iSp].c_str(),iBand),histTh_title.c_str(),200,0,90);

      hR[iSp][iDet][iBand]->Sumw2();
      hR_pzG0[iSp][iDet][iBand]->Sumw2();
      hR_pzL0[iSp][iDet][iBand]->Sumw2();
      hTh[iSp][iDet][iBand]->Sumw2();
      hTh_pzG0[iSp][iDet][iBand]->Sumw2();
      hTh_pzL0[iSp][iDet][iBand]->Sumw2();

      niceLogBins(hE[iSp][iDet][iBand]);
      niceLogBins(hE_pzG0[iSp][iDet][iBand]);
      niceLogBins(hE_pzL0[iSp][iDet][iBand]);
     }
    }
   }
    outfile = TFile::Open(Form("./rootfiles/ene_th_rad_xy_usScanner_motors_%s_EG%.0f.root",tgt_gen_config.c_str(),kinEcut),"RECREATE");
    int nfile=0;
    Long64_t nentry=0;
    long nTotEv=0;
    for(int ifile=1001;ifile<=6000;ifile++){
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
           if(hit->at(j).k<kinEcut) continue;

           double px = hit->at(j).px;
           double py = hit->at(j).py;
           double pz = hit->at(j).pz;
           double theta = (180.0/pi)*atan2(sqrt(px*px+py*py),pz);

           for(int iBand=0;iBand<nBand;iBand++){
            if(hit->at(j).r>=rIn[iBand] && hit->at(j).r<=rOut[iBand]){
             hE[sp][dt][iBand]->Fill(hit->at(j).k,rate);
             hR[sp][dt][iBand]->Fill(hit->at(j).r,rate);
             hXY[sp][dt][iBand]->Fill(hit->at(j).x,hit->at(j).y,rate);
             hTh[sp][dt][iBand]->Fill(theta,rate);
             if(hit->at(j).pz>=0){
              hE_pzG0[sp][dt][iBand]->Fill(hit->at(j).k,rate);
              hR_pzG0[sp][dt][iBand]->Fill(hit->at(j).r,rate);
              hXY_pzG0[sp][dt][iBand]->Fill(hit->at(j).x,hit->at(j).y,rate);
              hTh_pzG0[sp][dt][iBand]->Fill(theta,rate);
             }else{
              hE_pzL0[sp][dt][iBand]->Fill(hit->at(j).k,rate);
              hR_pzL0[sp][dt][iBand]->Fill(hit->at(j).r,rate);
              hXY_pzL0[sp][dt][iBand]->Fill(hit->at(j).x,hit->at(j).y,rate);
              hTh_pzL0[sp][dt][iBand]->Fill(theta,rate);
             }
             if(hit->at(j).vz<=-3875 && (hit->at(j).pid==11 || hit->at(j).pid==-211)){
              hE[4][dt][iBand]->Fill(hit->at(j).k,rate);
              hR[4][dt][iBand]->Fill(hit->at(j).r,rate);
              hXY[4][dt][iBand]->Fill(hit->at(j).x,hit->at(j).y,rate);
              hTh[4][dt][iBand]->Fill(theta,rate);
              if(hit->at(j).pz>=0){
               hE_pzG0[4][dt][iBand]->Fill(hit->at(j).k,rate);
               hR_pzG0[4][dt][iBand]->Fill(hit->at(j).r,rate);
               hXY_pzG0[4][dt][iBand]->Fill(hit->at(j).x,hit->at(j).y,rate);
               hTh_pzG0[4][dt][iBand]->Fill(theta,rate);
              }else{
               hE_pzL0[4][dt][iBand]->Fill(hit->at(j).k,rate);
               hR_pzL0[4][dt][iBand]->Fill(hit->at(j).r,rate);
               hXY_pzL0[4][dt][iBand]->Fill(hit->at(j).x,hit->at(j).y,rate);
               hTh_pzL0[4][dt][iBand]->Fill(theta,rate);
              }
             }
             if(hit->at(j).trid==1 && hit->at(j).pid==11){
               hE[5][dt][iBand]->Fill(hit->at(j).k,rate);
               hR[5][dt][iBand]->Fill(hit->at(j).r,rate);
               hXY[5][dt][iBand]->Fill(hit->at(j).x,hit->at(j).y,rate);
               hTh[5][dt][iBand]->Fill(theta,rate);
               if(hit->at(j).pz>=0){
                hE_pzG0[5][dt][iBand]->Fill(hit->at(j).k,rate);
                hR_pzG0[5][dt][iBand]->Fill(hit->at(j).r,rate);
                hXY_pzG0[5][dt][iBand]->Fill(hit->at(j).x,hit->at(j).y,rate);
                hTh_pzG0[5][dt][iBand]->Fill(theta,rate);
               }else{
                hE_pzL0[5][dt][iBand]->Fill(hit->at(j).k,rate);
                hR_pzL0[5][dt][iBand]->Fill(hit->at(j).r,rate);
                hXY_pzL0[5][dt][iBand]->Fill(hit->at(j).x,hit->at(j).y,rate);
                hTh_pzL0[5][dt][iBand]->Fill(theta,rate);
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
    outfile->cd();
    for(int iDet=0;iDet<nDet;iDet++){
      outfile->mkdir(Form("%s",detH[iDet].c_str()));
      outfile->cd(Form("%s",detH[iDet].c_str()));
      for(int iSp=0;iSp<nSp;iSp++){
       for(int iBand=0;iBand<nBand;iBand++){
        if(beamGen){
          hE[iSp][iDet][iBand]->Scale(1.0/nTotEv);
          hE_pzG0[iSp][iDet][iBand]->Scale(1.0/nTotEv);
          hE_pzL0[iSp][iDet][iBand]->Scale(1.0/nTotEv);
          hR[iSp][iDet][iBand]->Scale(1.0/nTotEv);
          hR_pzG0[iSp][iDet][iBand]->Scale(1.0/nTotEv);
          hR_pzL0[iSp][iDet][iBand]->Scale(1.0/nTotEv);
          hXY[iSp][iDet][iBand]->Scale(1.0/nTotEv);
          hXY_pzG0[iSp][iDet][iBand]->Scale(1.0/nTotEv);
          hXY_pzL0[iSp][iDet][iBand]->Scale(1.0/nTotEv);
          hTh[iSp][iDet][iBand]->Scale(1.0/nTotEv);
          hTh_pzG0[iSp][iDet][iBand]->Scale(1.0/nTotEv);
          hTh_pzL0[iSp][iDet][iBand]->Scale(1.0/nTotEv);
        }else{
          hE[iSp][iDet][iBand]->Scale(1.0e-9/nfile);
          hE_pzG0[iSp][iDet][iBand]->Scale(1.0e-9/nfile);
          hE_pzL0[iSp][iDet][iBand]->Scale(1.0e-9/nfile);
          hR[iSp][iDet][iBand]->Scale(1.0e-9/nfile);
          hR_pzG0[iSp][iDet][iBand]->Scale(1.0e-9/nfile);
          hR_pzL0[iSp][iDet][iBand]->Scale(1.0e-9/nfile);
          hXY[iSp][iDet][iBand]->Scale(1.0e-9/nfile);
          hXY_pzG0[iSp][iDet][iBand]->Scale(1.0e-9/nfile);
          hXY_pzL0[iSp][iDet][iBand]->Scale(1.0e-9/nfile);
          hTh[iSp][iDet][iBand]->Scale(1.0e-9/nfile);
          hTh_pzG0[iSp][iDet][iBand]->Scale(1.0e-9/nfile);
          hTh_pzL0[iSp][iDet][iBand]->Scale(1.0e-9/nfile);
        }
        hE[iSp][iDet][iBand]->Write();
        hE_pzG0[iSp][iDet][iBand]->Write();
        hE_pzL0[iSp][iDet][iBand]->Write();
        hR[iSp][iDet][iBand]->Write();
        hR_pzG0[iSp][iDet][iBand]->Write();
        hR_pzL0[iSp][iDet][iBand]->Write();
        hXY[iSp][iDet][iBand]->Write();
        hXY_pzG0[iSp][iDet][iBand]->Write();
        hXY_pzL0[iSp][iDet][iBand]->Write();
        hTh[iSp][iDet][iBand]->Write();
        hTh_pzG0[iSp][iDet][iBand]->Write();
        hTh_pzL0[iSp][iDet][iBand]->Write();
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
