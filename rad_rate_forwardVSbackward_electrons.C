//plot radial and transverse hit distribution on ring5, lam, usscanner, and dsscanner plane
//
#include "remolltypes.hh"
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>

void rad_rate_forwardVSbackward_electrons(){
   gROOT->Reset();
   gStyle->SetOptStat(0);
   gStyle->SetTitleYOffset(1.3);
   gStyle->SetPadGridX(1);
   gStyle->SetPadGridY(1);
   TGaxis::SetMaxDigits(3);

   int x_min = 0;
   int x_max = 1000;
   int nbin = 500;
   int beam_curr = 65;

   TString rootfile_dir = "$VOLATILE/remoll_rootfiles/PhotonBlocker";
   
   TH1F* h_pzAll = new TH1F("h_pzAll",Form("%s dist. on %s det plane;Radius (mm);hits/%sthrownEvents","Radial","SAM","#"),nbin,x_min,x_max);
   TH1F* h_pzG0 = new TH1F("h_pzG0",Form("%s dist. on %s det plane;Radius (mm);hits/%sthrownEvents","Radial","SAM","#"),nbin,x_min,x_max);
   TH1F* h_pzL0 = new TH1F("h_pzL0",Form("%s dist. on %s det plane;Radius (mm);hits/%sthrownEvents","Radial","SAM","#"),nbin,x_min,x_max);
   TH1F* h_pzAll_sam = new TH1F("h_pzAll_sam",Form("%s dist. on %s det plane;Radius (mm);hits/%sthrownEvents","Radial","SAM","#"),nbin,x_min,x_max);
   TH1F* h_pzG0_sam = new TH1F("h_pzG0_sam",Form("%s dist. on %s det plane;Radius (mm);hits/%sthrownEvents","Radial","SAM","#"),nbin,x_min,x_max);
   TH1F* h_pzL0_sam = new TH1F("h_pzL0_sam",Form("%s dist. on %s det plane;Radius (mm);hits/%sthrownEvents","Radial","SAM","#"),nbin,x_min,x_max);

   h_pzAll->SetLineColor(1);
   h_pzAll_sam->SetLineColor(1);
   h_pzG0->SetLineColor(2);
   h_pzG0_sam->SetLineColor(2);
   h_pzL0->SetLineColor(4);
   h_pzL0_sam->SetLineColor(4);
   h_pzAll->Sumw2();
   h_pzAll_sam->Sumw2();
   h_pzG0->Sumw2();
   h_pzG0_sam->Sumw2();
   h_pzL0->Sumw2();
   h_pzL0_sam->Sumw2();
   TChain* T = new TChain("T");
   int nfile = 0;
   for(int ifile = 1001;ifile<=2000;ifile++){
      nfile++;
      T->Add(Form("%s/pB_beam/pB_beam_%d.root",rootfile_dir.Data(),ifile));
   }
   cout<<Form("Found %d file splits!!!",nfile)<<endl;

   Long64_t nentry = T->GetEntries();
   std::vector<remollGenericDetectorHit_t> *fHit =0;
   T->SetBranchAddress("hit", &fHit);
   
   Double_t energy, hitr, rate, hitpz;
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
	 hitpz = fHit->at(pk).pz;
         rate = 1;//Rate is taken 1 for beam generator
         if(detector==176 && energy>1 && hitr>0 && pid==11){
           h_pzAll->Fill(hitr,rate);
  	   if(hitr>45 && hitr<60)
             h_pzAll_sam->Fill(hitr,rate);
           if(hitpz>0){
             h_pzG0->Fill(hitr,rate);
  	   if(hitr>45 && hitr<60)
             h_pzG0_sam->Fill(hitr,rate);
	   }else{
	     h_pzL0->Fill(hitr,rate);
  	   if(hitr>45 && hitr<60)
             h_pzL0_sam->Fill(hitr,rate);
           }
	 }
       }
   }
   
   TLatex latex;
   latex.SetNDC(1);
   latex.SetTextSize(0.04);
   h_pzAll->Scale(1./nentry);
   h_pzG0->Scale(1./nentry);
   h_pzL0->Scale(1./nentry);
   h_pzAll_sam->Scale(1./nentry);
   h_pzG0_sam->Scale(1./nentry);
   h_pzL0_sam->Scale(1./nentry);
 
   TString intRange = Form("%.0f<r<%.0f",45.0,60.0);
   double integral_pzAll = h_pzAll_sam->Integral();
   double integral_pzG0 = h_pzG0_sam->Integral();
   double integral_pzL0 = h_pzL0_sam->Integral();

   TCanvas* c_rate_linear= new TCanvas("c_rate_linear");
   h_pzAll->Draw("hist");
   h_pzG0->Draw("hist same");
   h_pzL0->Draw("hist same");
  latex.SetTextColor(1);
   latex.DrawLatex(0.45,0.85,Form("%s(%s):%.3e","Integral pzAll",intRange.Data(),integral_pzAll));
   latex.SetTextColor(2);
   latex.DrawLatex(0.45,0.80,Form("%s(%s):%.3e","Integral pzG0",intRange.Data(),integral_pzG0));
   latex.SetTextColor(4);
   latex.DrawLatex(0.45,0.75,Form("%s(%s):%.3e","Integral pzL0",intRange.Data(),integral_pzL0));
   c_rate_linear->SaveAs("./temp/p_rate_linear.pdf");

   TCanvas* c_rate_log= new TCanvas("c_rate_log");
   gPad->SetLogy();
   h_pzAll->Draw("hist");
   h_pzG0->Draw("hist same");
   h_pzL0->Draw("hist same");
  latex.SetTextColor(1);
   latex.DrawLatex(0.45,0.85,Form("%s(%s):%.3e","Integral pzAll",intRange.Data(),integral_pzAll));
   latex.SetTextColor(2);
   latex.DrawLatex(0.45,0.80,Form("%s(%s):%.3e","Integral pzG0",intRange.Data(),integral_pzG0));
   latex.SetTextColor(4);
   latex.DrawLatex(0.45,0.75,Form("%s(%s):%.3e","Integral pzL0",intRange.Data(),integral_pzL0));
   c_rate_log->SaveAs("./temp/p_rate_log.pdf");

//Now combine all pdf files saved in ./temp/ directory and save a single pdf file in ./plots/ directory
   gSystem->Exec(Form("pdfunite ./temp/*.pdf ./plots/rad_rate_forwardVSbackward_electrons.pdf"));
   gSystem->Exec(Form("rm -rf ./temp/*.pdf"));
}
