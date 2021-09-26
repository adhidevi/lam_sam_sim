//plot radial hit distribution on ring5, lam, usscanner, and dsscanner plane
//
#include "remolltypes.hh"
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>

void ds_scanner_radial_dist_magON(){
   gROOT->Reset();
   gStyle->SetOptStat(0);
   gStyle->SetTitleYOffset(1.3);
   gStyle->SetPadGridX(1);
   gStyle->SetPadGridY(1);
   TGaxis::SetMaxDigits(3);

   Double_t rmin = 500;
   Double_t rmax = 700;
   int beam_curr = 65;

   int colorsim[] = {1,2};//for total events of ee, ep-el
   int colorsimQ[] = {3,6};//for quartz accepted ee, ep-el

   TString rootfile_dir = "$VOLATILE/remoll_rootfiles/PhotonBlocker/";
   TString file[] = {"pB_ee10/pB_ee10","thcommin002_epel/thcommin002_epel"};
   int nfile = sizeof(file)/sizeof(*file);
   TString sim[] = {"EE mag on", "EP_el mag on"};
   
   double x_min = 200;
   double x_max = 1000;
   int nbin=x_max-x_min;
   double bin_width = (x_max-x_min)/nbin;
   TH1F* h_rate[nfile];
   
   for(int ifile=0;ifile<nfile;ifile++){
      h_rate[ifile] = new TH1F(Form("h_rate[%d]",ifile),Form("Radial Distribution Downstream Linear Scanner location with Magnet On;Radius (mm);rate (GHz)/(%duA*%.1fmm)",beam_curr,bin_width),nbin,x_min,x_max);
      h_rate[ifile]->SetLineColor(colorsim[ifile]);
      h_rate[ifile]->Sumw2();

      TChain* T = new TChain("T");
      int nfileSplit=0;
      for(int fileSplit=1001;fileSplit<=1020;fileSplit++){
        if(fileSplit==1016) continue;
          nfileSplit++;
          T->Add(rootfile_dir+Form("%s_%d.root",file[ifile].Data(),fileSplit));
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
	   if(detector==176 && energy>1){
             h_rate[ifile]->Fill(hitr,rate);
	   }
         }
      }
   }
   
   TH1F* h_rateQ[nfile];
   for(int ifile=0;ifile<nfile;ifile++){
      h_rateQ[ifile] = (TH1F*)h_rate[ifile]->Clone(Form("h_rateQ[%d]",ifile));
      h_rateQ[ifile]->GetXaxis()->SetRangeUser(rmin,rmax);
      h_rateQ[ifile]->SetLineColor(colorsimQ[ifile]);
    }

   TLine* line[2];
   TLatex latex;
   latex.SetNDC(1);
   latex.SetTextSize(0.04);
   TCanvas* c_rate_linear = new TCanvas("c_rate_linear");
   h_rate[1]->Draw("hist");
   h_rateQ[1]->Draw("same hist");
   h_rate[0]->Draw("same hist");
   h_rateQ[0]->Draw("same hist");
   line[0] = new TLine(rmin,0.0,rmin,h_rate[1]->GetMaximum()/2.0);
   line[1] = new TLine(rmax,0.0,rmax,h_rate[1]->GetMaximum()/2.0);
   for(int ii=0;ii<2;ii++){
      line[ii]->SetLineWidth(2);
      line[ii]->SetLineColor(kRed-2);
      line[ii]->Draw();
   }
   for(int ifile=0;ifile<nfile;ifile++){
      latex.SetTextColor(colorsim[ifile]);
      latex.DrawLatex(0.47,0.85-0.05*ifile,Form("%s",sim[ifile].Data()));
      latex.SetTextColor(colorsimQ[ifile]);
      latex.DrawLatex(0.67,0.85-0.05*ifile,Form("%s:%.1fGHz","accept",h_rateQ[ifile]->Integral()));
   }
   c_rate_linear->SaveAs("./temp/p_rate_linear.pdf");

   TCanvas* c_rate_log = new TCanvas("c_rate_log");
   gPad->SetLogy();
   h_rate[1]->Draw("hist");
   h_rateQ[1]->Draw("same hist");
   h_rate[0]->Draw("same hist");
   h_rateQ[0]->Draw("same hist");
   line[0] = new TLine(rmin,0.0,rmin,h_rate[1]->GetMaximum()/2.0);
   line[1] = new TLine(rmax,0.0,rmax,h_rate[1]->GetMaximum()/2.0);
   for(int ii=0;ii<2;ii++){
      line[ii]->SetLineWidth(2);
      line[ii]->SetLineColor(kRed-2);
      line[ii]->Draw();
   }
   for(int ifile=0;ifile<nfile;ifile++){
      latex.SetTextColor(colorsim[ifile]);
      latex.DrawLatex(0.47,0.85-0.05*ifile,Form("%s",sim[ifile].Data()));
      latex.SetTextColor(colorsimQ[ifile]);
      latex.DrawLatex(0.67,0.85-0.05*ifile,Form("%s:%.1fGHz","accept",h_rateQ[ifile]->Integral()));
   }
      c_rate_log->SaveAs("./temp/p_rate_log.pdf");
//Now combine all pdf files saved in ./temp/ directory and save a single pdf file in ./plots/ directory
   gSystem->Exec(Form("pdfunite ./temp/*.pdf ./plots/allParticles_dist_ds_scanner_magON_09222021.pdf"));
   gSystem->Exec(Form("rm -rf ./temp/*.pdf"));
}
