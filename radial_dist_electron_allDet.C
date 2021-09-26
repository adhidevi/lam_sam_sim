//plot radial hit distribution on ring5, lam, usscanner, and dsscanner plane
//
#include "remolltypes.hh"
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>

void radial_dist_electron_allDet(){
   gROOT->Reset();
   gStyle->SetOptStat(0);
   gStyle->SetTitleYOffset(1.3);
   gStyle->SetPadGridX(1);
   gStyle->SetPadGridY(1);
   TGaxis::SetMaxDigits(3);

   TString DetName[] = {"ring5","lam","sam"};
   const int nDet = sizeof(DetName)/sizeof(*DetName);
   int DetNo[nDet] = {28,174,176};
   Double_t rmin[nDet] = {920,1013,45};
   Double_t rmax[nDet] = {1060,1133,60};
   Double_t eCut[nDet] = {1000,1,1};
   Double_t rCut[nDet] = {500,500,100};
   int beam_curr = 65;

   int colorsim[] = {1,2,4};//for total events of ee, ep-el, and ep-inel
   int colorsimQ[] = {3,6,7};//for quartz accepted ee, ep-el, and ep-inel

   TString rootfile_dir = "$VOLATILE/remoll_rootfiles/PhotonBlocker/";
   TString file[] = {"sim_ee_lam_usscan_sam/sim_ee_lam_usscan_sam_*.root","sim_epel_lam_usscan_sam/sim_epel_lam_usscan_sam_*.root","sim_epinel_lam_usscan_sam/sim_epinel_lam_usscan_sam_*.root"};
   int nfile = sizeof(file)/sizeof(*file);
   TString sim[] = {"EE", "EP_el", "EP_inel"};
   
   int nbin = 500;
   double x_min[nDet] = {500,500,30};
   double x_max[nDet] = {1500,1500,100};
   TH1F* h_asym[nfile][nDet];
   TH1F* h_energy[nfile][nDet];
   TH1F* h_rate[nfile][nDet];
   
   for(int ifile=0;ifile<nfile;ifile++){
      for(int iDet=0;iDet<nDet;iDet++){
         double bin_width = (x_max[iDet]-x_min[iDet])/nbin;
         h_asym[ifile][iDet] = new TH1F(Form("h_asym[%d][%d]",ifile,iDet),Form("%s dist. on %s det plane (rate*Asym weighted);hit.r (mm);rate*Asym (GHz*ppb)/(%duA*%.1fmm)","Radial",DetName[iDet].Data(),beam_curr,bin_width),nbin,x_min[iDet],x_max[iDet]);
         h_energy[ifile][iDet] = new TH1F(Form("h_energy[%d][%d]",ifile,iDet),Form("%s dist. on %s det plane (rate*E weighted);hit.r (mm);rate*E (GHz*GeV)/(%duA*%.1fmm)","Radial",DetName[iDet].Data(),beam_curr,bin_width),nbin,x_min[iDet],x_max[iDet]);
         h_rate[ifile][iDet] = new TH1F(Form("h_rate[%d][%d]",ifile,iDet),Form("%s dist. on %s det plane (rate weighted);hit.r (mm);rate (GHz)/(%duA*%.1fmm)","Radial",DetName[iDet].Data(),beam_curr,bin_width),nbin,x_min[iDet],x_max[iDet]);
         h_asym[ifile][iDet]->SetLineColor(colorsim[ifile]);
         h_energy[ifile][iDet]->SetLineColor(colorsim[ifile]);
         h_rate[ifile][iDet]->SetLineColor(colorsim[ifile]);
         h_asym[ifile][iDet]->Sumw2();
         h_energy[ifile][iDet]->Sumw2();
         h_rate[ifile][iDet]->Sumw2();
      }
      TChain* T = new TChain("T");
      T->Add(rootfile_dir+Form("%s",file[ifile].Data()));
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
            rate = fRate/1.e9/20.0;//Convert to GHz, divide by number of file split
           if(detector==28 && energy>1000 && hitr>500 && pid==11){
             h_asym[ifile][0]->Fill(hitr,asym*rate);
             h_energy[ifile][0]->Fill(hitr,energy*1e-3*rate);
             h_rate[ifile][0]->Fill(hitr,rate);
           }
	   if(detector==174 && energy>1 && hitr>500 && pid==11){
             h_asym[ifile][1]->Fill(hitr,asym*rate);
             h_energy[ifile][1]->Fill(hitr,energy*1e-3*rate);
             h_rate[ifile][1]->Fill(hitr,rate);
	   }
	   if(detector==176 && energy>1 && hitr>30 && pid==11){
             h_asym[ifile][2]->Fill(hitr,asym*rate);
             h_energy[ifile][2]->Fill(hitr,energy*1e-3*rate);
             h_rate[ifile][2]->Fill(hitr,rate);
	   }
         }
      }
   }
   
   TH1F* h_asymQ[nfile][nDet];
   TH1F* h_energyQ[nfile][nDet];
   TH1F* h_rateQ[nfile][nDet];
   for(int ifile=0;ifile<nfile;ifile++){
   for(int iDet=0;iDet<nDet;iDet++){
      h_asymQ[ifile][iDet] = (TH1F*)h_asym[ifile][iDet]->Clone(Form("h_asymQ[%d][%d]",ifile,iDet));
      h_asymQ[ifile][iDet]->GetXaxis()->SetRangeUser(rmin[iDet],rmax[iDet]);
      h_asymQ[ifile][iDet]->SetLineColor(colorsimQ[ifile]);
      h_energyQ[ifile][iDet] = (TH1F*)h_energy[ifile][iDet]->Clone(Form("h_energyQ[%d][%d]",ifile,iDet));
      h_energyQ[ifile][iDet]->GetXaxis()->SetRangeUser(rmin[iDet],rmax[iDet]);
      h_energyQ[ifile][iDet]->SetLineColor(colorsimQ[ifile]);
      h_rateQ[ifile][iDet] = (TH1F*)h_rate[ifile][iDet]->Clone(Form("h_rateQ[%d][%d]",ifile,iDet));
      h_rateQ[ifile][iDet]->GetXaxis()->SetRangeUser(rmin[iDet],rmax[iDet]);
      h_rateQ[ifile][iDet]->SetLineColor(colorsimQ[ifile]);
    }
    }

   TLine* line[2];
   TLatex latex;
   latex.SetNDC(1);
   latex.SetTextSize(0.04);
   TCanvas* c_asym_linear[nDet];
   for(int iDet=0;iDet<nDet;iDet++){
      c_asym_linear[iDet] = new TCanvas(Form("c_asym_linear_%d",iDet));
      h_asym[1][iDet]->Draw("hist");
      h_asymQ[1][iDet]->Draw("same hist");
      h_asym[0][iDet]->Draw("same hist");
      h_asymQ[0][iDet]->Draw("same hist");
      h_asym[2][iDet]->Draw("same hist");
      h_asymQ[2][iDet]->Draw("same hist");
      line[0] = new TLine(rmin[iDet],0.0,rmin[iDet],h_asym[1][iDet]->GetMaximum()/2.0);
      line[1] = new TLine(rmax[iDet],0.0,rmax[iDet],h_asym[1][iDet]->GetMaximum()/2.0);
      for(int ii=0;ii<2;ii++){
         line[ii]->SetLineWidth(2);
         line[ii]->SetLineColor(kRed-2);
         line[ii]->Draw();
      }
      for(int ifile=0;ifile<nfile;ifile++){
         latex.SetTextColor(colorsim[ifile]);
         latex.DrawLatex(0.57,0.85-0.05*ifile,Form("%s",sim[ifile].Data()));
         latex.SetTextColor(colorsimQ[ifile]);
         latex.DrawLatex(0.67,0.85-0.05*ifile,Form("%s:%.1fGHz*ppb","accept",h_asymQ[ifile][iDet]->Integral()));
      }
      c_asym_linear[iDet]->SaveAs(Form("./temp/p_asym_linear_%d.pdf",iDet));
   }

   TCanvas* c_asym_log[nDet];
   for(int iDet=0;iDet<nDet;iDet++){
      c_asym_log[iDet] = new TCanvas(Form("c_asym_log_%d",iDet));
      gPad->SetLogy();
      h_asym[1][iDet]->Draw("hist");
      h_asymQ[1][iDet]->Draw("same hist");
      h_asym[0][iDet]->Draw("same hist");
      h_asymQ[0][iDet]->Draw("same hist");
      h_asym[2][iDet]->Draw("same hist");
      h_asymQ[2][iDet]->Draw("same hist");
      line[0] = new TLine(rmin[iDet],0.0,rmin[iDet],h_asym[1][iDet]->GetMaximum()/2.0);
      line[1] = new TLine(rmax[iDet],0.0,rmax[iDet],h_asym[1][iDet]->GetMaximum()/2.0);
      for(int ii=0;ii<2;ii++){
         line[ii]->SetLineWidth(2);
         line[ii]->SetLineColor(kRed-2);
         line[ii]->Draw();
      }
      for(int ifile=0;ifile<nfile;ifile++){
         latex.SetTextColor(colorsim[ifile]);
         latex.DrawLatex(0.57,0.85-0.05*ifile,Form("%s",sim[ifile].Data()));
         latex.SetTextColor(colorsimQ[ifile]);
         latex.DrawLatex(0.67,0.85-0.05*ifile,Form("%s:%.1fGHz*ppb","accept",h_asymQ[ifile][iDet]->Integral()));
      }
      c_asym_log[iDet]->SaveAs(Form("./temp/p_asym_log_%d.pdf",iDet));
   }

   TCanvas* c_energy_linear[nDet];
   for(int iDet=0;iDet<nDet;iDet++){
      c_energy_linear[iDet] = new TCanvas(Form("c_energy_linear_%d",iDet));
      h_energy[1][iDet]->Draw("hist");
      h_energyQ[1][iDet]->Draw("same hist");
      h_energy[0][iDet]->Draw("same hist");
      h_energyQ[0][iDet]->Draw("same hist");
      h_energy[2][iDet]->Draw("same hist");
      h_energyQ[2][iDet]->Draw("same hist");
      line[0] = new TLine(rmin[iDet],0.0,rmin[iDet],h_energy[1][iDet]->GetMaximum()/2.0);
      line[1] = new TLine(rmax[iDet],0.0,rmax[iDet],h_energy[1][iDet]->GetMaximum()/2.0);
      for(int ii=0;ii<2;ii++){
         line[ii]->SetLineWidth(2);
         line[ii]->SetLineColor(kRed-2);
         line[ii]->Draw();
      }
      for(int ifile=0;ifile<nfile;ifile++){
         latex.SetTextColor(colorsim[ifile]);
         latex.DrawLatex(0.57,0.85-0.05*ifile,Form("%s",sim[ifile].Data()));
         latex.SetTextColor(colorsimQ[ifile]);
         latex.DrawLatex(0.67,0.85-0.05*ifile,Form("%s:%.1fGHz*GeV","accept",h_energyQ[ifile][iDet]->Integral()));
      }
      c_energy_linear[iDet]->SaveAs(Form("./temp/p_energy_linear_%d.pdf",iDet));
   }

   TCanvas* c_energy_log[nDet];
   for(int iDet=0;iDet<nDet;iDet++){
      c_energy_log[iDet] = new TCanvas(Form("c_energy_log_%d",iDet));
      gPad->SetLogy();
      h_energy[1][iDet]->Draw("hist");
      h_energyQ[1][iDet]->Draw("same hist");
      h_energy[0][iDet]->Draw("same hist");
      h_energyQ[0][iDet]->Draw("same hist");
      h_energy[2][iDet]->Draw("same hist");
      h_energyQ[2][iDet]->Draw("same hist");
      line[0] = new TLine(rmin[iDet],0.0,rmin[iDet],h_energy[1][iDet]->GetMaximum()/2.0);
      line[1] = new TLine(rmax[iDet],0.0,rmax[iDet],h_energy[1][iDet]->GetMaximum()/2.0);
      for(int ii=0;ii<2;ii++){
         line[ii]->SetLineWidth(2);
         line[ii]->SetLineColor(kRed-2);
         line[ii]->Draw();
      }
      for(int ifile=0;ifile<nfile;ifile++){
         latex.SetTextColor(colorsim[ifile]);
         latex.DrawLatex(0.57,0.85-0.05*ifile,Form("%s",sim[ifile].Data()));
         latex.SetTextColor(colorsimQ[ifile]);
         latex.DrawLatex(0.67,0.85-0.05*ifile,Form("%s:%.1fGHz*GeV","accept",h_energyQ[ifile][iDet]->Integral()));
      }
      c_energy_log[iDet]->SaveAs(Form("./temp/p_energy_log_%d.pdf",iDet));
   }

   TCanvas* c_rate_linear[nDet];
   for(int iDet=0;iDet<nDet;iDet++){
      c_rate_linear[iDet] = new TCanvas(Form("c_rate_linear_%d",iDet));
      if(iDet==2){
      h_rate[1][iDet]->Draw("hist");
      h_rateQ[1][iDet]->Draw("same hist");
      h_rate[0][iDet]->Draw("same hist");
      h_rateQ[0][iDet]->Draw("same hist");
      h_rate[2][iDet]->Draw("same hist");
      h_rateQ[2][iDet]->Draw("same hist");
      }else{
      h_rate[0][iDet]->Draw("hist");
      h_rateQ[0][iDet]->Draw("same hist");
      h_rate[1][iDet]->Draw("same hist");
      h_rateQ[1][iDet]->Draw("same hist");
      h_rate[2][iDet]->Draw("same hist");
      h_rateQ[2][iDet]->Draw("same hist");
      }
      line[0] = new TLine(rmin[iDet],0.0,rmin[iDet],h_rate[1][iDet]->GetMaximum()/2.0);
      line[1] = new TLine(rmax[iDet],0.0,rmax[iDet],h_rate[1][iDet]->GetMaximum()/2.0);
      for(int ii=0;ii<2;ii++){
         line[ii]->SetLineWidth(2);
         line[ii]->SetLineColor(kRed-2);
         line[ii]->Draw();
      }
      for(int ifile=0;ifile<nfile;ifile++){
         latex.SetTextColor(colorsim[ifile]);
         latex.DrawLatex(0.57,0.85-0.05*ifile,Form("%s",sim[ifile].Data()));
         latex.SetTextColor(colorsimQ[ifile]);
         latex.DrawLatex(0.67,0.85-0.05*ifile,Form("%s:%.1fGHz","accept",h_rateQ[ifile][iDet]->Integral()));
      }
      c_rate_linear[iDet]->SaveAs(Form("./temp/p_rate_linear_%d.pdf",iDet));
   }

   TCanvas* c_rate_log[nDet];
   for(int iDet=0;iDet<nDet;iDet++){
      c_rate_log[iDet] = new TCanvas(Form("c_rate_log_%d",iDet));
      gPad->SetLogy();
      if(iDet==2){
      h_rate[1][iDet]->Draw("hist");
      h_rateQ[1][iDet]->Draw("same hist");
      h_rate[0][iDet]->Draw("same hist");
      h_rateQ[0][iDet]->Draw("same hist");
      h_rate[2][iDet]->Draw("same hist");
      h_rateQ[2][iDet]->Draw("same hist");
      }else{
      h_rate[0][iDet]->Draw("hist");
      h_rateQ[0][iDet]->Draw("same hist");
      h_rate[1][iDet]->Draw("same hist");
      h_rateQ[1][iDet]->Draw("same hist");
      h_rate[2][iDet]->Draw("same hist");
      h_rateQ[2][iDet]->Draw("same hist");
      }
      line[0] = new TLine(rmin[iDet],0.0,rmin[iDet],h_rate[1][iDet]->GetMaximum()/2.0);
      line[1] = new TLine(rmax[iDet],0.0,rmax[iDet],h_rate[1][iDet]->GetMaximum()/2.0);
      for(int ii=0;ii<2;ii++){
         line[ii]->SetLineWidth(2);
         line[ii]->SetLineColor(kRed-2);
         line[ii]->Draw();
      }
      for(int ifile=0;ifile<nfile;ifile++){
         latex.SetTextColor(colorsim[ifile]);
         latex.DrawLatex(0.57,0.85-0.05*ifile,Form("%s",sim[ifile].Data()));
         latex.SetTextColor(colorsimQ[ifile]);
         latex.DrawLatex(0.67,0.85-0.05*ifile,Form("%s:%.1fGHz","accept",h_rateQ[ifile][iDet]->Integral()));
      }
      c_rate_log[iDet]->SaveAs(Form("./temp/p_rate_log_%d.pdf",iDet));
   }
//Now combine all pdf files saved in ./temp/ directory and save a single pdf file in ./plots/ directory
   gSystem->Exec(Form("pdfunite ./temp/*.pdf ./plots/radial_dist_electron_allDet.pdf"));
   gSystem->Exec(Form("rm -rf ./temp/radial*.pdf"));
}
