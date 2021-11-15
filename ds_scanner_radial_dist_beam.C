//
#include "remolltypes.hh"
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>

void ds_scanner_radial_dist_beam(){
   gROOT->Reset();
   gStyle->SetOptStat(0);
   gStyle->SetTitleYOffset(1.3);
   gStyle->SetPadGridX(1);
   gStyle->SetPadGridY(1);
   TGaxis::SetMaxDigits(3);

   Double_t rmin = 500;
   Double_t rmax = 700;
   int beam_curr = 65;

   TString rootfile_dir = "$VOLATILE/remoll_rootfiles/default-geo/";
   
   double x_min = 200;
   double x_max = 1000;
   int nbin=200;
   double bin_width = (x_max-x_min)/nbin;
   TH1F* h_rate = new TH1F(Form("h_rate"),Form("Radial dist. on ds scanner plane (beam gen);Radius (mm);hits/#thrownEvents/%.1fmm",bin_width),nbin,x_min,x_max);
   h_rate->SetLineColor(1);
   h_rate->Sumw2();

   TChain* T = new TChain("T");
   int nfileSplit=0;
   for(int fileSplit=1001;fileSplit<=6000;fileSplit++){
       nfileSplit++;
       T->Add(rootfile_dir+Form("Optics1_beam_magOFF_V5/Optics1_beam_magOFF_V5_%d.root",fileSplit));
   }
   cout<<Form("Found %d number of file splits!",nfileSplit)<<endl;
   Long64_t nentry = T->GetEntries();
   std::vector<remollGenericDetectorHit_t> *fHit =0;
   T->SetBranchAddress("hit", &fHit);
   
   Double_t energy, hitr, asym, rate;
   Int_t detector, pid;
   for(int ientry=0;ientry<nentry;ientry++){
      if(ientry%(nentry/10)==0)
        cout<<"analyzed "<<ientry<<" events!!"<<endl;
      T->GetEntry(ientry);
      for(size_t pk=0;pk<fHit->size();pk++){
         pid = (Int_t)TMath::Abs(fHit->at(pk).pid);
         detector = fHit->at(pk).det;
         energy = fHit->at(pk).k;
         hitr = fHit->at(pk).r;
         rate = 1;
        if(detector==176 && energy>1 && (hitr>200 && hitr<1000) && pid==11){
          h_rate->Fill(hitr,rate);
        }
      }
   }
   
   TH1F* h_rateQ;
   h_rateQ = (TH1F*)h_rate->Clone(Form("h_rateQ"));
   h_rateQ->GetXaxis()->SetRangeUser(rmin,rmax);
   h_rateQ->SetLineColor(2);
   h_rate->Scale(1.0/nentry);
   h_rateQ->Scale(1.0/nentry);

   TLatex latex;
   latex.SetNDC(1);
   latex.SetTextSize(0.04);
   TCanvas* c_rate_linear = new TCanvas("c_rate_linear");
   h_rate->SetTitle("Radial Distribution of Electrons at Downstream Linear Scanner location with Magnet Off");
   h_rate->Draw("hist");
   h_rateQ->Draw("hist same");
   latex.SetTextColor(1);
   latex.DrawLatex(0.45,0.85,"beam generator (Optics1 DSC target)");
   latex.SetTextColor(2);
   latex.DrawLatex(0.45,0.80,"ds scanner acceptance");
   c_rate_linear->SaveAs("./temp/ds_scanner_linear.pdf");

   TCanvas* c_rate_log = new TCanvas("c_rate_log");
   gPad->SetLogy();
   h_rate->Draw("hist");
   h_rateQ->Draw("hist same");
   latex.SetTextColor(1);
   latex.DrawLatex(0.45,0.86,"beam generator (Optics1 DSC target)");
   latex.SetTextColor(2);
   latex.DrawLatex(0.45,0.82,"ds scanner acceptance");
   c_rate_log->SaveAs("./temp/ds_scanner_log.pdf");

//Now combine all pdf files saved in ./temp/ directory and save a single pdf file in ./plots/ directory
   gSystem->Exec(Form("pdfunite ./temp/ds_scanner_*.pdf ./plots/defaultGeo_electrons_ds_scanner_Optics1_beam_magOFF_V5.pdf"));
   gSystem->Exec(Form("rm -rf ./temp/ds_scanner_*.pdf"));
}
