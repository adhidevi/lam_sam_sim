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

void loop_vs_tdraw(int particleID = 2112){
   TString particle;
   if(particleID==11) particle = "e-";
   if(particleID==22) particle = "#gamma";
   if(particleID==211) particle = "pi+";
   if(particleID==2112) particle = "neutron";

   TString rootfile_dir = "/volatile/halla/parity/adhidevi/remoll_rootfiles/PMTShielding/";
   TChain* T = new TChain("T");
   TChain* TT = new TChain("T");
   int nfile = 0;
   for(int ifile=1001;ifile<=1100;ifile++){
      nfile++;
      T->Add(Form("%sPMTSh_beam/PMTSh_beam_%d.root",rootfile_dir.Data(),ifile));
      TT->Add(Form("%sPMTSh_beam/PMTSh_beam_%d.root",rootfile_dir.Data(),ifile));
   }
   cout<<Form("Found %d files!!!",nfile)<<endl;

   TH1D* hist1 = new TH1D(Form("hist1_%d",particleID),"Hit Distribution using loop (1200<=r<=1500) ring 5",121,-8,5); 
   niceLogXBins(hist1);
   hist1->SetLineColor(1);

   TH1D* hist2 = new TH1D(Form("hist2_%d",particleID),"Hit Distribution using t-draw(1200<=r<=1500) ring 5",121,-8,5); 
   niceLogXBins(hist2);
   hist2->SetLineColor(2);
   Long64_t nentry = T->GetEntries();
   std::vector<remollGenericDetectorHit_t> *hit=0;
   T->SetBranchAddress("hit", &hit);
   
   for(int ientry=0;ientry<nentry;T->GetEntry(ientry++)){
      if(ientry%(nentry/10)==0)
        cout<<"analyzed "<<ientry<<" events!!"<<endl;
      for(int j=0;j<hit->size();j++){
         if(hit->at(j).pid!=particleID) continue;
         if(hit->at(j).det!=28) continue;
         if(hit->at(j).r<1200) continue;
         if(hit->at(j).r>1500) continue;
         hist1->Fill(hit->at(j).k);
      }
   }
   TT->Draw(Form("hit.k>>hist2_%d",particleID),Form("hit.det==28 && hit.pid==%d && (hit.r>=1200 && hit.r<=1500)",particleID),"goff");
   TCanvas* c1 = new TCanvas();
   c1->SetLogx();
   c1->SetLogy();
   hist1->Scale(1./nentry);
   hist2->Scale(1./nentry);
   hist1->Draw("hist");
   hist1->GetYaxis()->SetRangeUser(1.e-9,1.e-2);
   hist2->Draw("hist sames");
   
   gPad->Update();
   TPaveStats* st1 = (TPaveStats*)hist1->FindObject("stats");
   TPaveStats* st2 = (TPaveStats*)hist2->FindObject("stats");
   st1->SetY1NDC(0.90);
   st1->SetY2NDC(0.75);
   st2->SetY1NDC(0.75);
   st2->SetY2NDC(0.60);
   st1->SetTextColor(1);
   st2->SetTextColor(2);
   gPad->Modified();

   TLatex latex;
   latex.SetNDC();
   latex.SetTextSize(0.04);
   latex.SetTextColor(1);
   latex.DrawLatex(0.15,0.85,Form("%s (using loop)",particle.Data()));
   latex.SetTextColor(2);
   latex.DrawLatex(0.15,0.80,"using t-draw");
   c1->SaveAs("./plots/loop_vs_tdraw.pdf");
}
