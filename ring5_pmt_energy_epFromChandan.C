#include "../remolltypes.hh"
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>

void niceLogXBins(TH1*h){
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

void ring5_pmt_energy_epFromChandan(){
   gROOT->Reset();
   gStyle->SetOptStat(0);
   gStyle->SetTitleYOffset(1.3);
   gStyle->SetPadGridX(1);
   gStyle->SetPadGridY(1);
   TGaxis::SetMaxDigits(3);

   TString rootfile_dir = "/lustre19/expphy/volatile/halla/moller12gev/chandan/sim_out/Col6_design_ep_Oct18_sym/";
   
   int color[] = {1,2,3,4,6,7};
   const int nhist = sizeof(color)/sizeof(*color);
   TString particle_type[nhist] = {"gamma","e-/pi-","e+/pi+","neutron","e-/e+ E>1 MeV","primary E>1 MeV"};
   TH1F* h_E[nhist];
   for(int ihist=0;ihist<nhist;ihist++){
      h_E[ihist] = new TH1F(Form("h_E%d",ihist),Form("Energy dist. on Ring5 PMT region (beam gen);Energy (MeV);Rate (Hz/65uA)"),121,-8,4.1);
      h_E[ihist]->SetLineColor(color[ihist]);
      h_E[ihist]->Sumw2();
      niceLogXBins(h_E[ihist]);
   }
   TChain* T = new TChain("T");
   int nfileSplit=0;
   for(int fileSplit=1;fileSplit<=4400;fileSplit++){
       ifstream fileName(rootfile_dir+Form("remollout_ep%d.root",fileSplit));
       if(!fileName){
         cout<<Form("File %sremollout_ep%d.root doesn't exist. Escaping this...!!",rootfile_dir.Data(),fileSplit)<<endl;
         continue;
       }
       nfileSplit++;
       T->Add(rootfile_dir+Form("remollout_ep%d.root",fileSplit));
   }
   cout<<Form("Found %d number of file splits!",nfileSplit)<<endl;
   Long64_t nentry = T->GetEntries();
   std::vector<remollGenericDetectorHit_t> *fHit =0;
   Double_t fRate = 0.0;

   T->SetBranchAddress("hit", &fHit);
   T->SetBranchAddress("rate", &fRate);

   Double_t energy, hitr, asym, rate, trid, kinE;
   Int_t detector, pid;
   for(int ientry=0;ientry<nentry;ientry++){
      if(ientry%(nentry/10)==0)
        cout<<"analyzed "<<ientry<<" events!!"<<endl;
      T->GetEntry(ientry);
      for(size_t pk=0;pk<fHit->size();pk++){
         pid = (Int_t)TMath::Abs(fHit->at(pk).pid);
         detector = fHit->at(pk).det;
         energy = fHit->at(pk).e;
         kinE = fHit->at(pk).k;
         hitr = fHit->at(pk).r;
         trid = fHit->at(pk).trid;
         rate = fRate;
        if(detector==28 && hitr>1200 && hitr<1500){
          if(pid==22){
            h_E[0]->Fill(kinE,rate);
          }
          if(pid==11 || pid==-211){
            h_E[1]->Fill(kinE,rate);
          }
          if(pid==-11 || pid==211){
            h_E[2]->Fill(kinE,rate);
          }
          if(pid==2112){
            h_E[3]->Fill(kinE,rate);
          cout<<energy<<endl;
          }
          if((pid==11 || pid==-11) && kinE>1){
            h_E[4]->Fill(kinE,rate);
          }
          if(pid==11 && trid==1 && kinE>1){
            h_E[5]->Fill(kinE,rate);
          }
        }
      }
   }
   for(int ihist=0;ihist<nhist;ihist++){
      h_E[ihist]->Scale(1.0/nfileSplit);
   }
   TLatex latex;
   latex.SetNDC(1);
   latex.SetTextSize(0.04);
   TCanvas* c_E = new TCanvas("c_E");
   gPad->SetLogy();
   gPad->SetLogx();
   
   for(int ihist=0;ihist<nhist;ihist++){
      h_E[ihist]->SetOption("E");
      h_E[ihist]->Draw("hist same");
   }
   for(int ihist=0;ihist<nhist;ihist++){
      h_E[ihist]->SetTitle("Energy Distribution at Ring5 PMT Region (beam generator)");
      h_E[ihist]->GetYaxis()->SetRangeUser(1.e6,1.e10);
      latex.SetTextColor(color[ihist]);
      if(ihist<3)
        latex.DrawLatex(0.55,0.85-0.05*ihist,Form("%s",particle_type[ihist].Data()));
      else
        latex.DrawLatex(0.70,0.85-0.05*(ihist-3),Form("%s",particle_type[ihist].Data()));
   }
   c_E->SaveAs("./plots/LH2_epFromChandan_ring5_pmt_mom.pdf");

}
