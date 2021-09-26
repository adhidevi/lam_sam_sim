//plot radial hit distribution on ring5 and showermax plane
//
#include "remolltypes.hh"
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>

void radial_distribution_electron(){
   gROOT->Reset();
   gStyle->SetOptStat(0);
   gStyle->SetTitleYOffset(1.3);
   gStyle->SetPadGridX(1);
   gStyle->SetPadGridY(1);
   TGaxis::SetMaxDigits(3);

   Double_t ring5_rmin = 920;
   Double_t ring5_rmax = 1060;
   Double_t lam_rmin = 1013;
   Double_t lam_rmax = lam_rmin+120;
   int beam_curr = 65;

   int colorsim[] = {1,2,4};//for total events of ee, ep-el, and ep-inel
   int colorsimQ[] = {3,6,7};//for quartz accepted ee, ep-el, and ep-inel

   TString rootfile_dir = "$VOLATILE/remoll_rootfiles/PhotonBlocker/";
   TString file[] = {"sim_ee_lam/sim_ee_lam_*.root","sim_epel_lam/sim_epel_lam_*.root","sim_epinel_lam/sim_epinel_lam_*.root"};
        int nfile = sizeof(file)/sizeof(*file);
   TString sim[] = {"EE", "EP_el", "EP_inel"};
   
   int nbin = 500;
   double x_min = 500;
   double x_max = 1500;
   double bin_width = (x_max-x_min)/nbin;
   TH1F* h_ring5[nfile];
   TH1F* h_lam[nfile];
   TH1F* h_ring5rate[nfile];
   TH1F* h_lamrate[nfile];

   for(int ifile=0;ifile<nfile;ifile++){
      h_ring5[ifile] = new TH1F(Form("h_fing5[%d]",ifile),Form("%s distribution on ring 5 (rate*Asym weighted);hit.r (mm);rate*Asym (GHz*ppb)/(%duA*%.1fmm)","Radial",beam_curr,bin_width),nbin,x_min,x_max);
      h_ring5rate[ifile] = new TH1F(Form("h_fing5Qrate[%d]",ifile),Form("%s distribution on ring 5 (rate weighted);hit.r (mm);rate (GHz)/(%duA*%.1fmm)","Radial",beam_curr,bin_width),nbin,x_min,x_max);
      h_lam[ifile] = new TH1F(Form("h_lam[%d]",ifile),Form("%s distribution on lam det (rate*Asym weighted);hit.r (mm);rate*Asym (GHz*ppb)/(%duA*%.1fmm)","Radial",beam_curr,bin_width),nbin,x_min,x_max);
      h_lamrate[ifile] = new TH1F(Form("h_lamrate[%d]",ifile),Form("%s distribution on lam det (rate weighted);hit.r (mm);rate (GHz)/(%duA*%.1fmm)","Radial",beam_curr,bin_width),nbin,x_min,x_max);
      h_ring5[ifile]->SetLineColor(colorsim[ifile]);
      h_ring5rate[ifile]->SetLineColor(colorsim[ifile]);
      h_lam[ifile]->SetLineColor(colorsim[ifile]);
      h_lamrate[ifile]->SetLineColor(colorsim[ifile]);
      h_ring5[ifile]->Sumw2();
      h_ring5rate[ifile]->Sumw2();
      h_lam[ifile]->Sumw2();
      h_lamrate[ifile]->Sumw2();

      TChain* T = new TChain("T");
      T->Add(rootfile_dir+Form("%s",file[ifile].Data()));
      Long64_t nentry = T->GetEntries();
      std::vector<remollGenericDetectorHit_t> *fHit =0;
      remollEvent_t *fEv =0;
      Double_t fRate=0.;
      T->SetBranchAddress("hit", &fHit);
      T->SetBranchAddress("ev", &fEv);
      T->SetBranchAddress("rate", &fRate);
   
      Float_t energy(-1.e-12), hitr(-1.e-12), detector(-1.e-12), asym(-1.e-12), phi(-1.e-12), modphi(-1.e-12), rate(-1.e-12);
      Int_t pid(0), trid(0), mtrid(0);
   
      for(int ientry=0;ientry<nentry;ientry++){
         if(ientry%(nentry/10)==0)
           cout<<"analyzed "<<ientry<<" events!!"<<endl;
         T->GetEntry(ientry);
         for(size_t pk=0;pk<fHit->size();pk++){
            pid = (Int_t)TMath::Abs(fHit->at(pk).pid);
            detector = fHit->at(pk).det;
            energy = fHit->at(pk).e;
            hitr = fHit->at(pk).r;
	    trid = fHit->at(pk).trid;
	    mtrid = fHit->at(pk).mtrid;
            asym = -1*fEv->A;
            rate = fRate/1.e9/20.0;//Convert to GHz, divide by number of file split
            phi = fHit->at(pk).ph;
            if(detector==28 && energy>1000 && hitr>500 && pid==11){
              h_ring5[ifile]->Fill(hitr,asym*rate);
              h_ring5rate[ifile]->Fill(hitr,rate);
            }
            if(detector==174 && energy>1 && hitr>500 && pid==11){
              h_lam[ifile]->Fill(hitr,asym*rate);
              h_lamrate[ifile]->Fill(hitr,rate);
            }
         }
      }
   }
   TH1F* h_ring5Q[nfile];
   TH1F* h_lamQ[nfile];
   TH1F* h_ring5Qrate[nfile];
   TH1F* h_lamQrate[nfile];
   for(int ifile=0;ifile<nfile;ifile++){
      h_ring5Q[ifile] = (TH1F*)h_ring5[ifile]->Clone(Form("h_ring5Q[%d]",ifile));
      h_ring5Q[ifile]->GetXaxis()->SetRangeUser(ring5_rmin,ring5_rmax);
      h_ring5Q[ifile]->SetLineColor(colorsimQ[ifile]);
      h_ring5Qrate[ifile] = (TH1F*)h_ring5rate[ifile]->Clone(Form("h_ring5Qrate%d]",ifile));
      h_ring5Qrate[ifile]->GetXaxis()->SetRangeUser(ring5_rmin,ring5_rmax);
      h_ring5Qrate[ifile]->SetLineColor(colorsimQ[ifile]);
      h_lamQ[ifile] = (TH1F*)h_lam[ifile]->Clone(Form("h_lamQ[%d]",ifile));
      h_lamQ[ifile]->GetXaxis()->SetRangeUser(lam_rmin,lam_rmax);
      h_lamQ[ifile]->SetLineColor(colorsimQ[ifile]);
      h_lamQrate[ifile] = (TH1F*)h_lamrate[ifile]->Clone(Form("h_lamQrate[%d]",ifile));
      h_lamQrate[ifile]->GetXaxis()->SetRangeUser(lam_rmin,lam_rmax);
      h_lamQrate[ifile]->SetLineColor(colorsimQ[ifile]);
    }
//Let's plot the rate*asym weighted ring5 radial distribution with linear scale in y
   TCanvas* ring5_1_linear = new TCanvas("ring5_1_linear");
   h_ring5[1]->Draw("hist");
   h_ring5[0]->Draw("same hist");
   h_ring5[2]->Draw("same hist");
   h_ring5Q[0]->Draw("same hist");
   h_ring5Q[1]->Draw("same hist");
   h_ring5Q[2]->Draw("same hist");

   TLine* line[2];
   line[0] = new TLine(ring5_rmin,0.0,ring5_rmin,h_ring5[0]->GetMaximum()/2.0);
   line[1] = new TLine(ring5_rmax,0.0,ring5_rmax,h_ring5[0]->GetMaximum()/2.0);
   line[0]->SetLineWidth(2);
   line[1]->SetLineWidth(2);
   line[0]->SetLineColor(kRed-2);
   line[1]->SetLineColor(kRed-2);
   line[0]->Draw();
   line[1]->Draw();

   TLatex latex;
   latex.SetNDC(1);
   latex.SetTextSize(0.04);
   for(int ifile=0;ifile<nfile;ifile++){
   latex.SetTextColor(colorsim[ifile]);
   latex.DrawLatex(0.57,0.85-0.05*ifile,Form("%s",sim[ifile].Data()));
   latex.SetTextColor(colorsimQ[ifile]);
   latex.DrawLatex(0.67,0.85-0.05*ifile,Form("%s:%.1fGHz*ppb","accept",h_ring5Q[ifile]->Integral()));
//   latex.DrawLatex(0.12,0.85-0.05*ifile,Form("Int=%.1fGHz*ppb",h_ring5Q[ifile]->Integral()));
   }
   ring5_1_linear->SaveAs(Form("./temp/radial1_linear_ring5.pdf"));

//Let's plot the rate*asym weighted ring5 radial distribution with log scale in y
   TCanvas* ring5_1_log = new TCanvas("ring5_1_log");
   gPad->SetLogy(1);
   h_ring5[1]->Draw("hist");
   h_ring5[0]->Draw("hist same");
   h_ring5[2]->Draw("hist same");
   h_ring5Q[0]->Draw("hist same");
   h_ring5Q[1]->Draw("hist same");
   h_ring5Q[2]->Draw("hist same");
   for(int ifile=0;ifile<nfile;ifile++){
   latex.SetTextColor(colorsim[ifile]);
   latex.DrawLatex(0.57,0.85-0.05*ifile,Form("%s",sim[ifile].Data()));
   latex.SetTextColor(colorsimQ[ifile]);
   latex.DrawLatex(0.67,0.85-0.05*ifile,Form("%s:%.1fGHz*ppb","accept",h_ring5Q[ifile]->Integral()));
//   latex.DrawLatex(0.12,0.85-0.05*ifile,Form("Int=%.1fGHz*ppb",h_ring5Q[ifile]->Integral()));
   }
   line[0]->Draw();
   line[1]->Draw();
   ring5_1_log->SaveAs(Form("./temp/radial1_log_ring5.pdf"));

//Let's plot the rate weighted ring5 radial distribution with log scale in y
   TCanvas* ring5_2_linear = new TCanvas("ring5_2_linear");
   h_ring5rate[0]->Draw("hist");
   h_ring5rate[1]->Draw("hist same");
   h_ring5rate[2]->Draw("hist same");
   h_ring5Qrate[0]->Draw("hist same");
   h_ring5Qrate[1]->Draw("hist same");
   h_ring5Qrate[2]->Draw("hist same");
   for(int ifile=0;ifile<nfile;ifile++){
   latex.SetTextColor(colorsim[ifile]);
   latex.DrawLatex(0.57,0.85-0.05*ifile,Form("%s",sim[ifile].Data()));
   latex.SetTextColor(colorsimQ[ifile]);
   latex.DrawLatex(0.67,0.85-0.05*ifile,Form("%s:%.1fGHz","accept",h_ring5Qrate[ifile]->Integral()));
//   latex.DrawLatex(0.12,0.85-0.05*ifile,Form("Int=%.1fGHz",h_ring5Qrate[ifile]->Integral()));
   }
   line[0] = new TLine(ring5_rmin,0.0,ring5_rmin,h_ring5rate[0]->GetMaximum()/2.0);
   line[1] = new TLine(ring5_rmax,0.0,ring5_rmax,h_ring5rate[0]->GetMaximum()/2.0);
   line[0]->SetLineWidth(2);
   line[1]->SetLineWidth(2);
   line[0]->SetLineColor(kRed-2);
   line[1]->SetLineColor(kRed-2);
   line[0]->Draw();
   line[1]->Draw();
   ring5_2_linear->SaveAs(Form("./temp/radial2_linear_ring5.pdf"));

//Let's plot the rate weighted ring5 radial distribution with log scale in y
   TCanvas* ring5_2_log = new TCanvas("ring5_2_log");
   gPad->SetLogy(1);
   h_ring5rate[0]->Draw("hist");
   h_ring5rate[1]->Draw("hist same");
   h_ring5rate[2]->Draw("hist same");
   h_ring5Qrate[0]->Draw("hist same");
   h_ring5Qrate[1]->Draw("hist same");
   h_ring5Qrate[2]->Draw("hist same");
   for(int ifile=0;ifile<nfile;ifile++){
   latex.SetTextColor(colorsim[ifile]);
   latex.DrawLatex(0.57,0.85-0.05*ifile,Form("%s",sim[ifile].Data()));
   latex.SetTextColor(colorsimQ[ifile]);
   latex.DrawLatex(0.67,0.85-0.05*ifile,Form("%s:%.1fGHz","accept",h_ring5Qrate[ifile]->Integral()));
//   latex.DrawLatex(0.12,0.85-0.05*ifile,Form("Int=%.1fGHz",h_ring5Qrate[ifile]->Integral()));
   }
   line[0]->Draw();
   line[1]->Draw();
   ring5_2_log->SaveAs(Form("./temp/radial2_log_ring5.pdf"));

//Let's plot the rate*asym weighted lam radial distribution with linear scale in y
   TCanvas* lam_1_linear = new TCanvas("lam_1_linear");
   h_lam[1]->Draw("hist");
   h_lam[0]->Draw("same hist");
   h_lam[2]->Draw("same hist");
   h_lamQ[0]->Draw("same hist");
   h_lamQ[1]->Draw("same hist");
   h_lamQ[2]->Draw("same hist");

   line[0] = new TLine(lam_rmin,0.0,lam_rmin,h_lam[0]->GetMaximum()/2.0);
   line[1] = new TLine(lam_rmax,0.0,lam_rmax,h_lam[0]->GetMaximum()/2.0);
   line[0]->SetLineWidth(2);
   line[1]->SetLineWidth(2);
   line[0]->SetLineColor(kRed-2);
   line[1]->SetLineColor(kRed-2);
   line[0]->SetLineWidth(2);
   line[1]->SetLineWidth(2);
   line[0]->SetLineColor(kRed-2);
   line[1]->SetLineColor(kRed-2);
   line[0]->Draw();
   line[1]->Draw();

   for(int ifile=0;ifile<nfile;ifile++){
   latex.SetTextColor(colorsim[ifile]);
   latex.DrawLatex(0.57,0.85-0.05*ifile,Form("%s",sim[ifile].Data()));
   latex.SetTextColor(colorsimQ[ifile]);
   latex.DrawLatex(0.67,0.85-0.05*ifile,Form("%s:%.1fGHz*ppb","accept",h_lamQ[ifile]->Integral()));
//   latex.DrawLatex(0.12,0.85-0.05*ifile,Form("Int=%.1fGHz*ppb",h_lamQ[ifile]->Integral()));
   }
   lam_1_linear->SaveAs(Form("./temp/radial1_linear_lam.pdf"));

//Let's plot the rate*asym weighted lam radial distribution with log scale in y
   TCanvas* lam_1_log = new TCanvas("lam_1_log");
   gPad->SetLogy(1);
   h_lam[1]->Draw("hist");
   h_lam[0]->Draw("hist same");
   h_lam[2]->Draw("hist same");
   h_lamQ[0]->Draw("hist same");
   h_lamQ[1]->Draw("hist same");
   h_lamQ[2]->Draw("hist same");
   for(int ifile=0;ifile<nfile;ifile++){
   latex.SetTextColor(colorsim[ifile]);
   latex.DrawLatex(0.57,0.85-0.05*ifile,Form("%s",sim[ifile].Data()));
   latex.SetTextColor(colorsimQ[ifile]);
   latex.DrawLatex(0.67,0.85-0.05*ifile,Form("%s:%.1fGHz*ppb","accept",h_lamQ[ifile]->Integral()));
//   latex.DrawLatex(0.12,0.85-0.05*ifile,Form("Int=%.1fGHz*ppb",h_lamQ[ifile]->Integral()));
   }
   line[0]->Draw();
   line[1]->Draw();
   lam_1_log->SaveAs(Form("./temp/radial1_log_lam.pdf"));

//Let's plot the rate weighted lam radial distribution with log scale in y
   TCanvas* lam_2_linear = new TCanvas("lam_2_linear");
   h_lamrate[0]->Draw("hist");
   h_lamrate[1]->Draw("hist same");
   h_lamrate[2]->Draw("hist same");
   h_lamQrate[0]->Draw("hist same");
   h_lamQrate[1]->Draw("hist same");
   h_lamQrate[2]->Draw("hist same");
   for(int ifile=0;ifile<nfile;ifile++){
   latex.SetTextColor(colorsim[ifile]);
   latex.DrawLatex(0.57,0.85-0.05*ifile,Form("%s",sim[ifile].Data()));
   latex.SetTextColor(colorsimQ[ifile]);
   latex.DrawLatex(0.67,0.85-0.05*ifile,Form("%s:%.1fGHz","accept",h_lamQrate[ifile]->Integral()));
//   latex.DrawLatex(0.12,0.85-0.05*ifile,Form("Int=%.1fGHz",h_lamQrate[ifile]->Integral()));
   }
   line[0] = new TLine(lam_rmin,0.0,lam_rmin,h_lamrate[0]->GetMaximum()/2.0);
   line[1] = new TLine(lam_rmax,0.0,lam_rmax,h_lamrate[0]->GetMaximum()/2.0);
   line[0]->SetLineWidth(2);
   line[1]->SetLineWidth(2);
   line[0]->SetLineColor(kRed-2);
   line[1]->SetLineColor(kRed-2);
   line[0]->Draw();
   line[1]->Draw();
   lam_2_linear->SaveAs(Form("./temp/radial2_linear_lam.pdf"));

//Let's plot the rate weighted lam radial distribution with log scale in y
   TCanvas* lam_2_log = new TCanvas("lam_2_log");
   gPad->SetLogy(1);
   h_lamrate[0]->Draw("hist");
   h_lamrate[1]->Draw("hist same");
   h_lamrate[2]->Draw("hist same");
   h_lamQrate[0]->Draw("hist same");
   h_lamQrate[1]->Draw("hist same");
   h_lamQrate[2]->Draw("hist same");
   for(int ifile=0;ifile<nfile;ifile++){
   latex.SetTextColor(colorsim[ifile]);
   latex.DrawLatex(0.57,0.85-0.05*ifile,Form("%s",sim[ifile].Data()));
   latex.SetTextColor(colorsimQ[ifile]);
   latex.DrawLatex(0.67,0.85-0.05*ifile,Form("%s:%.1fGHz","accept",h_lamQrate[ifile]->Integral()));
//   latex.DrawLatex(0.12,0.85-0.05*ifile,Form("Int=%.1fGHz",h_lamQrate[ifile]->Integral()));
   }
   line[0]->Draw();
   line[1]->Draw();
   lam_2_log->SaveAs(Form("./temp/radial2_log_lam.pdf"));
//Now combine all pdf files saved in ./temp/ directory and save a single pdf file in ./plots/ directory
   gSystem->Exec(Form("pdfunite ./temp/*.pdf ./plots/radial_distribution_electrons.pdf"));
   gSystem->Exec(Form("rm -rf ./temp/radial*.pdf"));
}
