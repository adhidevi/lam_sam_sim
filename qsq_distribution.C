#include "remolltypes.hh"
#include <sstream>
#include <iostream>
#include <fstream>
void qsq_distribution(){
  gROOT->Reset();
  gStyle->SetOptStat(0);
  gStyle->SetTitleYOffset(1.3);
  gStyle->SetPadGridX(1);
  gStyle->SetPadGridY(1);
  TGaxis::SetMaxDigits(3);
  
  const TString Ecut[] = {"all E", "KE>1 MeV"};
  const int nEcut = sizeof(Ecut)/sizeof(*Ecut);
  const string spTit[] = {"e-/#pi-","e+/#pi+","#gamma","neutron"};
  const int nSp = sizeof(spTit)/sizeof(*spTit);
  const string spH[nSp] = {"epiM","epiP","g","n"};
  map<int,int> spM {{11,1},{-211,1},{-11,2},{211,2},{22,3},{2112,4}};

////Change the following line for which detectors you want to include////
  string detH[] = {"det174"};
  const int nDet = sizeof(detH)/sizeof(detH);
  const int Det = {174};
  map<int,int> dtM {{174,1}};
  const double inRcut[] = {1010};
  const double outRcut[] = {1130};
  const int nRcut = sizeof(inRcut)/sizeof(*inRcut);
  const string weight[] = {"rate","rateE"};
  const int nWt = sizeof(weight)/sizeof(*weight);
  const string weight_unit[nWt] = {"hits/#thrownEvents"};


  TString rootfile_dir = "/volatile/halla/parity/adhidevi/remoll_rootfiles/develop_br";
  const string geometry = "develop";
  const string tgt_gen_config = "LH2_beam_V12";
  const string plotType = "qsq_dist";
  int beamGen(1);

  const int nbin = 500;
  const double qsq_min = 0.0;
  const double qsq_max = 0.06;
  const double r_min = 0.0;
  const double r_max = 1900;
  double qsq_bin_width = (qsq_max - qsq_min)/nbin;
  double r_bin_width = (r_max - r_min)/nbin;
  TString hst;
  TH1D* h_q2[nSp][nDet][nWt][nEcut];
  TH1D* h_q2pri[nSp][nDet][nWt][nEcut];
  TH1D* h_q2sec[nSp][nDet][nWt][nEcut];
  TH1D* h_R[nSp][nDet][nWt][nEcut];
  TH1D* h_Rpri[nSp][nDet][nWt][nEcut];
  TH1D* h_Rsec[nSp][nDet][nWt][nEcut];
  
  for(int iSp=0;iSp<nSp;iSp++){
   for(int iDet=0;iDet<nDet;iDet++){
    for(int iWt=0;iWt<nWt;iWt++){
     for(int iEcut=0;iEcut<nEcut;iEcut++){
       hst = Form("%s_%s_%s_%d",spH[iSp].c_str(),detH[iDet].c_str(),weight[iWt].c_str(),iEcut);
       string titleQ2 = Form("%s (%s) Qsq dist. on %s plane (%s);Qsq (GeV/c)^{2};%s/%.1e (GeV/c)^{2}",spTit[iSp].c_str(),Ecut[iEcut].Data(),detH[iDet].c_str(),tgt_gen_config.c_str(),weight_unit[iWt].c_str(),qsq_bin_width);
       h_q2[iSp][iDet][iWt][iEcut] = new TH1D(Form("qsq_%s",hst.Data()),titleQ2.c_str(),nbin,qsq_min,qsq_max);
       h_q2pri[iSp][iDet][iWt][iEcut] = new TH1D(Form("qsq_%s_pri",hst.Data()),titleQ2.c_str(),nbin,qsq_min,qsq_max);
       h_q2sec[iSp][iDet][iWt][iEcut] = new TH1D(Form("qsq_%s_sec",hst.Data()),titleQ2.c_str(),nbin,qsq_min,qsq_max);
       string titleR = Form("%s (%s) Radial dist. on %s plane (%s);Radius (mm);%s/%.1f (mm)",spTit[iSp].c_str(),Ecut[iEcut].Data(),detH[iDet].c_str(),tgt_gen_config.c_str(),weight_unit[iWt].c_str(),r_bin_width);
       h_R[iSp][iDet][iWt][iEcut] = new TH1D(Form("r_%s",hst.Data()),titleR.c_str(),nbin,r_min,r_max);
       h_Rpri[iSp][iDet][iWt][iEcut] = new TH1D(Form("r_%s_pri",hst.Data()),titleR.c_str(),nbin,r_min,r_max);
       h_Rsec[iSp][iDet][iWt][iEcut] = new TH1D(Form("r_%s_sec",hst.Data()),titleR.c_str(),nbin,r_min,r_max);
       h_R[iSp][iDet][iWt][iEcut]->Sumw2();
       h_Rpri[iSp][iDet][iWt][iEcut]->Sumw2();
       h_Rsec[iSp][iDet][iWt][iEcut]->Sumw2();
     }
    }
   }
  }
  int nfile = 0;
  Long64_t nentry = 0;
  long nTotEv = 0;
  for(int ifile=1001;ifile<=1010;ifile++){
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
    remollEvent_t *ev = 0;
    Double_t rate = 0;
    T->SetBranchAddress("hit", &hit);
    T->SetBranchAddress("ev", &ev);
    T->SetBranchAddress("rate", &rate);

    for(int ientry=0;ientry<nentry;ientry++){
      T->GetEntry(ientry);
      for(int j=0;j<hit->size();j++){
        if(std::isnan(rate) || std::isinf(rate)) continue;
        if(beamGen) rate = 1.0;

        int sp = spM[int(hit->at(j).pid)]-1;
        if(sp==-1) continue;
        int dt = dtM[int(hit->at(j).det)]-1;
        if(dt==-1) continue;
        if(hit->at(j).r<inRcut[0] || hit->at(j).r>outRcut[0]) continue;
        
        h_q2[sp][dt][0][0]->Fill(ev->Q2*1.e-6,rate);
        h_q2[sp][dt][1][0]->Fill(ev->Q2*1.e-6,rate*hit->at(j).e);
        h_R[sp][dt][0][0]->Fill(hit->at(j).r,rate);
        h_R[sp][dt][1][0]->Fill(hit->at(j).r,rate*hit->at(j).e);
        if(hit->at(j).vz<=-3875){
          h_q2pri[sp][dt][0][0]->Fill(ev->Q2*1.e-6,rate);
          h_q2pri[sp][dt][1][0]->Fill(ev->Q2*1.e-6,rate*hit->at(j).e);
          h_Rpri[sp][dt][0][0]->Fill(hit->at(j).r,rate);
          h_Rpri[sp][dt][1][0]->Fill(hit->at(j).r,rate*hit->at(j).e);
        }else{
          h_q2sec[sp][dt][0][0]->Fill(ev->Q2*1.e-6,rate);
          h_q2sec[sp][dt][1][0]->Fill(ev->Q2*1.e-6,rate*hit->at(j).e);
          h_Rsec[sp][dt][0][0]->Fill(hit->at(j).r,rate);
          h_Rsec[sp][dt][1][0]->Fill(hit->at(j).r,rate*hit->at(j).e);
        }

        if(hit->at(j).k>=1){
          h_q2[sp][dt][0][1]->Fill(ev->Q2*1.e-6,rate);
          h_q2[sp][dt][1][1]->Fill(ev->Q2*1.e-6,rate*hit->at(j).e);
          h_R[sp][dt][0][1]->Fill(hit->at(j).r,rate);
          h_R[sp][dt][1][1]->Fill(hit->at(j).r,rate*hit->at(j).e);
          if(hit->at(j).vz<=-3875){
            h_q2pri[sp][dt][0][1]->Fill(ev->Q2*1.e-6,rate);
            h_q2pri[sp][dt][1][1]->Fill(ev->Q2*1.e-6,rate*hit->at(j).e);
            h_Rpri[sp][dt][0][1]->Fill(hit->at(j).r,rate);
            h_Rpri[sp][dt][1][1]->Fill(hit->at(j).r,rate*hit->at(j).e);
          }else{
            h_q2sec[sp][dt][0][1]->Fill(ev->Q2*1.e-6,rate);
            h_q2sec[sp][dt][1][1]->Fill(ev->Q2*1.e-6,rate*hit->at(j).e);
            h_Rsec[sp][dt][0][1]->Fill(hit->at(j).r,rate);
            h_Rsec[sp][dt][1][1]->Fill(hit->at(j).r,rate*hit->at(j).e);
          }
        }
      }
    }
    delete fin;
  }
  cout<<Form("Total number of file splits: %d",nfile)<<endl;
  cout<<Form("Total number of entries: %ld",nTotEv)<<endl;

  TFile* outfile = new TFile(Form("./rootfiles/%s_%s_%s.root",geometry.c_str(),tgt_gen_config.c_str(),plotType.c_str()),"recreate");

  for(int iDet=0;iDet<nDet;iDet++){
    outfile->mkdir(Form("%s",detH[iDet].c_str()));
    outfile->cd(Form("%s",detH[iDet].c_str()));
    for(int iSp=0;iSp<nSp;iSp++){
      for(int iWt=0;iWt<nWt;iWt++){
        for(int iEcut=0;iEcut<nEcut;iEcut++){
          if(beamGen){
            h_q2[iSp][iDet][iWt][iEcut]->Scale(1.0/nTotEv);
            h_q2pri[iSp][iDet][iWt][iEcut]->Scale(1.0/nTotEv);
            h_q2sec[iSp][iDet][iWt][iEcut]->Scale(1.0/nTotEv);
            h_R[iSp][iDet][iWt][iEcut]->Scale(1.0/nTotEv);
            h_Rpri[iSp][iDet][iWt][iEcut]->Scale(1.0/nTotEv);
            h_Rsec[iSp][iDet][iWt][iEcut]->Scale(1.0/nTotEv);
          }else{
            h_q2[iSp][iDet][iWt][iEcut]->Scale(1.0e-9/nfile);
            h_q2pri[iSp][iDet][iWt][iEcut]->Scale(1.0e-9/nfile);
            h_q2sec[iSp][iDet][iWt][iEcut]->Scale(1.0e-9/nfile);
            h_R[iSp][iDet][iWt][iEcut]->Scale(1.0e-9/nfile);
            h_Rpri[iSp][iDet][iWt][iEcut]->Scale(1.0e-9/nfile);
            h_Rsec[iSp][iDet][iWt][iEcut]->Scale(1.0e-9/nfile);
          }
        h_q2[iSp][iDet][iWt][iEcut]->Write();
        h_q2pri[iSp][iDet][iWt][iEcut]->Write();
        h_q2sec[iSp][iDet][iWt][iEcut]->Write();
        h_R[iSp][iDet][iWt][iEcut]->Write();
        h_Rpri[iSp][iDet][iWt][iEcut]->Write();
        h_Rsec[iSp][iDet][iWt][iEcut]->Write();
        }
      }
    }
  }

}
