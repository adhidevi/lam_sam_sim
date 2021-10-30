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

void ring5_pmt_energy_dist_simulate(){
   gROOT->Reset();
   gStyle->SetOptStat(0);
   gStyle->SetTitleYOffset(1.3);
   gStyle->SetPadGridX(1);
   gStyle->SetPadGridY(1);
   TGaxis::SetMaxDigits(3);

   TH1F* h_E = new TH1F("h_E","Energy dist. on Ring5 PMT region (simulation);Energy (MeV);Events",121,-8,5);
   niceLogXBins(h_E);
   h_E->SetLineColor(1);
   h_E->Sumw2();
   TRandom2* rand = new TRandom2(0);
   double nentry = 1000000;

   for(int ientry=0;ientry<nentry;ientry++){
      double r = rand->Uniform(1e-8,1e5);
      h_E->Fill(r);
   }

   h_E->Scale(1.0/nentry);
   TCanvas* c_E = new TCanvas("c_E");
   gPad->SetLogy();
   gPad->SetLogx();
   h_E->Draw("hist");
}
