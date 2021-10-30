//
#include "remolltypes.hh"
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>

void ds_scanner_radial_dist_elasticC12_noiseRateCut(){
   gROOT->Reset();
   gStyle->SetOptStat(0);
   gStyle->SetTitleYOffset(1.3);
   gStyle->SetPadGridX(1);
   gStyle->SetPadGridY(1);
   TGaxis::SetMaxDigits(3);

   Double_t rmin = 500;
   Double_t rmax = 700;
   int beam_curr = 65;

   TString rootfile_dir = "$VOLATILE/remoll_rootfiles/PMTShielding/";
   
   double x_min = 500;
   double x_max = 700;
   int nbin=x_max-x_min;
   double bin_width = (x_max-x_min)/nbin;
   TH1F* h_rate;
   
   h_rate = new TH1F(Form("h_rate"),Form("%s dist. on ds-scanner plane (rate weighted);Radius (mm);rate (GHz)/(%duA*%.1fmm)","Radial",beam_curr,bin_width),nbin,x_min,x_max);
   h_rate->SetLineColor(1);
   h_rate->Sumw2();

   TChain* T = new TChain("T");
   int nfileSplit=0;
   for(int fileSplit=1001;fileSplit<=1100;fileSplit++){
       nfileSplit++;
       T->Add(rootfile_dir+Form("PMTSh_elC12_magOFF/PMTSh_elC12_magOFF_%d.root",fileSplit));
   }
   cout<<Form("Found %d number of file splits!",nfileSplit)<<endl;
   Long64_t nentry = T->GetEntries();
   std::vector<remollGenericDetectorHit_t> *fHit =0;
   remollEvent_t *fEv =0;
   Double_t fRate=0.;
   T->SetBranchAddress("hit", &fHit);
   T->SetBranchAddress("ev", &fEv);
   T->SetBranchAddress("rate", &fRate);
   
   Double_t energy, hitr, asym, rate;
   Int_t detector, pid;
   for(int ientry=0;ientry<nentry;ientry++){
      if(ientry%(nentry/10)==0)
        cout<<"analyzed "<<ientry<<" events!!"<<endl;
      T->GetEntry(ientry);
      for(size_t pk=0;pk<fHit->size();pk++){
         pid = (Int_t)TMath::Abs(fHit->at(pk).pid);
         detector = fHit->at(pk).det;
         energy = fHit->at(pk).e;
         hitr = fHit->at(pk).r;
         asym = -1*fEv->A;
         rate = fRate/1.e9/nfileSplit;//Convert to GHz, divide by number of file split
        if(detector==176 && energy>1 && hitr>100 && pid==11 && (rate<0.001)){
          h_rate->Fill(hitr,rate);
        }
      }
   }
   
   TH1F* h_rateQ;
   h_rateQ = (TH1F*)h_rate->Clone(Form("h_rateQ"));
   h_rateQ->GetXaxis()->SetRangeUser(rmin,rmax);
   h_rateQ->SetLineColor(2);

   TLatex latex;
   latex.SetNDC(1);
   latex.SetTextSize(0.04);
   TCanvas* c_rate_linear = new TCanvas("c_rate_linear");
   h_rate->SetTitle("Radial Distribution of Electrons at Downstream Linear Scanner location with Magnet Off");
   h_rate->Draw("hist");
   h_rateQ->Draw("hist same");
   latex.SetTextColor(1);
   latex.DrawLatex(0.55,0.85,"elasticC12 magnet off");
   latex.SetTextColor(2);
   latex.DrawLatex(0.55,0.80,"ds scanner acceptance");
   c_rate_linear->SaveAs("./temp/ds_scanner_linear.pdf");

   TCanvas* c_rate_log = new TCanvas("c_rate_log");
   gPad->SetLogy();
   h_rate->Draw("hist");
   h_rateQ->Draw("hist same");
   latex.SetTextColor(1);
   latex.DrawLatex(0.55,0.86,"elasticC12 magnet off");
   latex.SetTextColor(2);
   latex.DrawLatex(0.55,0.82,"ds scanner acceptance");
   c_rate_log->SaveAs("./temp/ds_scanner_log.pdf");

//Now combine all pdf files saved in ./temp/ directory and save a single pdf file in ./plots/ directory
   gSystem->Exec(Form("pdfunite ./temp/ds_scanner_*.pdf ./plots/PMTSh_elC12_magOFF_electrons_ds_scanner.pdf"));
   gSystem->Exec(Form("rm -rf ./temp/ds_scanner_*.pdf"));
}
