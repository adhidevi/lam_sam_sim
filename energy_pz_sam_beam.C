//plot energy distribution sam plane
//
#include "remolltypes.hh"
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>

void energy_pz_sam_beam(){
   gROOT->Reset();
   gStyle->SetOptStat(111111);
   gStyle->SetTitleYOffset(1.3);
   gStyle->SetPadGridX(1);
   gStyle->SetPadGridY(1);
   TGaxis::SetMaxDigits(3);

   int beam_curr = 65;

   TString rootfile_dir = "$VOLATILE/remoll_rootfiles/PhotonBlocker";
   int x_min = 0;
   int x_max = 12000; 
   int nbin = 200;
   TH1F* h_energy = new TH1F("h_energy",Form("%s dist. on %s det plane (beam generator);Energy (MeV);hits/%sthrownEvents","Energy","SAM","#"),nbin,x_min,x_max);
   TH1F* h_energy_cut = new TH1F("h_energy_cut",Form("%s dist. on %s det plane (beam generator);Energy (MeV);hits/%sthrownEvents","Energy","SAM","#"),nbin,x_min,x_max);
   TH1F* h_pz = new TH1F("h_pz",Form("%s dist. on %s det plane (beam generator);pz (MeV/c);hits/%sthrownEvents","pz","SAM","#"),nbin,-1000,x_max);

   h_energy->SetLineColor(1);
   h_energy_cut->SetLineColor(2);
   h_pz->SetLineColor(1);
   h_energy->Sumw2();
   h_energy_cut->Sumw2();
   h_pz->Sumw2();

   TChain* T = new TChain("T");
   int nfile = 0;
   for(int ifile = 1001;ifile<=2000;ifile++){
      nfile++;
      T->Add(Form("%s/pB_beam/pB_beam_%d.root",rootfile_dir.Data(),ifile));
   }
   cout<<Form("Found %d file splits!!!",nfile)<<endl;

   Long64_t nentry = T->GetEntries();
   std::vector<remollGenericDetectorHit_t> *fHit =0;
   remollEvent_t *fEv =0;
   Double_t fRate=0.;
   T->SetBranchAddress("hit", &fHit);
   T->SetBranchAddress("ev", &fEv);
   T->SetBranchAddress("rate", &fRate);
   
   Double_t energy, hitr, hitx, hity, rate, hitpz;
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
         hitx = fHit->at(pk).x;
         hity = fHit->at(pk).y;
	 hitpz = fHit->at(pk).pz;
         rate = 1;//Rate is taken 1 for beam generator
         if(detector==176 && energy>1 && (hitr>45 && hitr<60) && pid==11){
            h_energy->Fill(energy,rate);
            h_pz->Fill(hitpz,rate);
         }
         if(detector==176 && energy>10500 && (hitr>45 && hitr<60) && pid==11){
            h_energy_cut->Fill(energy,rate);
         }
      }
   }
   
   h_energy->Scale(1./nentry);
   h_energy_cut->Scale(1./nentry);
   h_pz->Scale(1./nentry);

   TLatex latex;
   latex.SetNDC(1);
   latex.SetTextSize(0.04);

   TCanvas* c_rate_linear = new TCanvas("c_rate_linear");
   h_energy->Draw("hist");
   h_energy_cut->Draw("hist same");
   latex.SetTextColor(1);
   latex.DrawLatex(0.15,0.85,"Cut::E>1 && r>45 && r<60 && PID==11");
   latex.DrawLatex(0.15,0.80,Form("Integral = %.2e",h_energy->Integral()));
   latex.SetTextColor(2);
   latex.DrawLatex(0.15,0.75,"Cut::E>10500 && r>45 && r<60 && PID==11");
   latex.DrawLatex(0.15,0.70,Form("Integral = %.2e",h_energy_cut->Integral()));
   gPad->Update();
   TPaveStats* st = (TPaveStats*)h_energy->FindObject("stats");
   st->SetFillStyle(0);
   c_rate_linear->SaveAs("./temp/energy_linear.pdf");

   TCanvas* c_rate_log = new TCanvas("c_rate_log");
   gPad->SetLogy();
   h_energy->Draw("hist");
   h_energy_cut->Draw("hist same");
   latex.SetTextColor(1);
   latex.DrawLatex(0.15,0.85,"Cut::E>1 && r>45 && r<60 && PID==11");
   latex.DrawLatex(0.15,0.80,Form("Integral = %.2e",h_energy->Integral()));
   latex.SetTextColor(2);
   latex.DrawLatex(0.15,0.75,"Cut::E>10500 && r>45 && r<60 && PID==11");
   latex.DrawLatex(0.15,0.70,Form("Integral = %.2e",h_energy_cut->Integral()));
   gPad->Update();
   TPaveStats* st_ = (TPaveStats*)h_energy->FindObject("stats");
   st_->SetFillStyle(0);
   c_rate_log->SaveAs("./temp/energy_log.pdf");

//Now combine all pdf files saved in ./temp/ directory and save a single pdf file in ./plots/ directory
   gSystem->Exec(Form("pdfunite ./temp/*.pdf ./plots/energy_dist_electron_beamGen.pdf"));
   gSystem->Exec(Form("rm -rf ./temp/*.pdf"));
}
