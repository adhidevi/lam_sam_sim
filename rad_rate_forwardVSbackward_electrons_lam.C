//plot radial and transverse hit distribution on ring5, lam, usscanner, and dsscanner plane
//
#include "../remolltypes.hh"
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>

void rad_rate_forwardVSbackward_electrons_lam(){
   gROOT->Reset();
   gStyle->SetOptStat(0);
   gStyle->SetTitleYOffset(1.3);
   gStyle->SetPadGridX(1);
   gStyle->SetPadGridY(1);
   TGaxis::SetMaxDigits(3);

   int x_min = 0;
   int x_max = 4500;
   int nbin = 1000;
   int beam_curr = 65;

   TString rootfile_dir = "$VOLATILE/remoll_rootfiles/PhotonBlocker";
   
   TH1F* h_pzAll = new TH1F("h_pzAll",Form("%s dist. on %s det plane;Radius (mm);hits/%sthrownEvents","Radial","SAM","#"),nbin,x_min,x_max);
   TH1F* h_pzG0 = new TH1F("h_pzG0",Form("%s dist. on %s det plane;Radius (mm);hits/%sthrownEvents","Radial","SAM","#"),nbin,x_min,x_max);
   TH1F* h_pzL0 = new TH1F("h_pzL0",Form("%s dist. on %s det plane;Radius (mm);hits/%sthrownEvents","Radial","SAM","#"),nbin,x_min,x_max);

   h_pzAll->SetLineColor(1);
   h_pzG0->SetLineColor(2);
   h_pzL0->SetLineColor(4);
   h_pzAll->Sumw2();
   h_pzG0->Sumw2();
   h_pzL0->Sumw2();
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
         if(detector==174 && energy>1 && hitr>952 && pid==11){
             h_pzAll->Fill(hitr,rate);
             if(hitpz>0)
             h_pzG0->Fill(hitr,rate);
	     else
	     h_pzL0->Fill(hitr,rate);
           }
       }
   }
   
   TLatex latex;
   latex.SetNDC(1);
   latex.SetTextSize(0.04);
   h_pzAll->Scale(1./nentry);
   h_pzG0->Scale(1./nentry);
   h_pzL0->Scale(1./nentry);
 
   TString intRange = Form("%.0f<=r<=%.0f",1013.0,1133.0);
   double integral_pzAll = h_pzAll->Integral(h_pzAll->FindBin(1013),h_pzAll->FindBin(1133));
   double integral_pzG0 = h_pzG0->Integral(h_pzG0->FindBin(1013),h_pzG0->FindBin(1133));
   double integral_pzL0 = h_pzL0->Integral(h_pzL0->FindBin(1013),h_pzL0->FindBin(1133));

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
   gSystem->Exec(Form("pdfunite ./temp/*.pdf ./plots/lam_rate_forwardVSbackward_electrons.pdf"));
   gSystem->Exec(Form("rm -rf ./temp/*.pdf"));
}
