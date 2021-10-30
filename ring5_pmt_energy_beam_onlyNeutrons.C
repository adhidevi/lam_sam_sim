#include "remolltypes.hh"
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

void ring5_pmt_energy_beam_onlyNeutrons(){

   TString rootfile_dir = "/volatile/halla/parity/adhidevi/remoll_rootfiles/PMTShielding/";
   TChain* T = new TChain("T");
   int nfile = 0;
   for(int ifile=1001;ifile<=2000;ifile++){
      nfile++;
      T->Add(Form("%sPMTSh_beam/PMTSh_beam_%d.root",rootfile_dir.Data(),ifile));
   }
   cout<<Form("Found %d files!!!",nfile)<<endl;

   TH1D* hist1 = new TH1D("hist1","Neutron Distribution (1200<=r<=1500) ring 5",121,-8,5); 
   TH1D* hist2 = new TH1D("hist2","e-/pi- Distribution (1200<=r<=1500) ring 5",121,-8,5); 
   TH1D* hist3 = new TH1D("hist3","e+/pi+ Distribution (1200<=r<=1500) ring 5",121,-8,5); 
   TH1D* hist4 = new TH1D("hist4","Gamma Distribution (1200<=r<=1500) ring 5",121,-8,5); 
   TH1D* hist5 = new TH1D("hist5","e-/e+ E>1 Distribution (1200<=r<=1500) ring 5",121,-8,5); 
   TH1D* hist6 = new TH1D("hist6","Primary E?1 Distribution (1200<=r<=1500) ring 5",121,-8,5); 
   niceLogXBins(hist1);
   niceLogXBins(hist2);
   niceLogXBins(hist3);
   niceLogXBins(hist4);
   niceLogXBins(hist5);
   niceLogXBins(hist6);
   hist1->SetLineColor(1);
   hist2->SetLineColor(2);
   hist3->SetLineColor(3);
   hist4->SetLineColor(4);
   hist5->SetLineColor(6);
   hist6->SetLineColor(7);
   Long64_t nentry = T->GetEntries();
   T->Draw("hit.k>>hist1","hit.det==28 && hit.pid==2112 && (hit.r>=1200 && hit.r<=1500)","goff");
   T->Draw("hit.k>>hist2","hit.det==28 && (hit.pid==11 || hit.pid==-211) && (hit.r>=1200 && hit.r<=1500)","goff");
   T->Draw("hit.k>>hist3","hit.det==28 && (hit.pid==-11 || hit.pid==211) && (hit.r>=1200 && hit.r<=1500)","goff");
   T->Draw("hit.k>>hist4","hit.det==28 && hit.pid==22 && (hit.r>=1200 && hit.r<=1500)","goff");
   T->Draw("hit.k>>hist5","hit.det==28 && (hit.pid==11 || hit.pid==-11) && hit.k>1 && (hit.r>=1200 && hit.r<=1500)","goff");
   T->Draw("hit.k>>hist6","hit.det==28 && hit.pid==11 && hit.trid==1 && hit.k>1 && (hit.r>=1200 && hit.r<=1500)","goff");
   TCanvas* c1 = new TCanvas();
   c1->SetLogx();
   c1->SetLogy();
   hist1->Scale(1./nentry);
   hist2->Scale(1./nentry);
   hist3->Scale(1./nentry);
   hist4->Scale(1./nentry);
   hist5->Scale(1./nentry);
   hist6->Scale(1./nentry);
   hist1->Draw("hist");
   hist1->GetYaxis()->SetRangeUser(1.e-9,1.e-2);
   hist2->Draw("hist same");
   hist3->Draw("hist same");
   hist4->Draw("hist same");
   hist5->Draw("hist same");
   hist6->Draw("hist same");

   TLatex latex;
   latex.SetNDC();
   latex.SetTextSize(0.04);
   latex.SetTextColor(1);
   latex.DrawLatex(0.15,0.85,"Neutron");
   latex.SetTextColor(2);
   latex.DrawLatex(0.15,0.80,"e-/#pi-");
   latex.SetTextColor(3);
   latex.DrawLatex(0.15,0.75,"e+/#pi+");
   latex.SetTextColor(4);
   latex.DrawLatex(0.30,0.85,"#gamma");
   latex.SetTextColor(6);
   latex.DrawLatex(0.30,0.80,"e-/e+ E>1");
   latex.SetTextColor(7);
   latex.DrawLatex(0.30,0.75,"e- primary E>1");
   c1->SaveAs("./plots/test_test_test.pdf");
}
