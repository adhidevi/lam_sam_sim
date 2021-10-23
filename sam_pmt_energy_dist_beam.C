//
#include "../remolltypes.hh"
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>

void sam_pmt_energy_dist_beam(){
   gROOT->Reset();
   gStyle->SetOptStat(0);
   gStyle->SetTitleYOffset(1.3);
   gStyle->SetPadGridX(1);
   gStyle->SetPadGridY(1);
   TGaxis::SetMaxDigits(3);

   Double_t rmin = 1200;
   Double_t rmax = 1500;
   int beam_curr = 65;

   TString rootfile_dir = "$VOLATILE/remoll_rootfiles/PMTShielding/";
   
   double x_minR = 1000;
   double x_maxR = 1700;
   double nbinR = 200;
   double bin_widthR = (x_maxR-x_minR)/nbinR;
   TH1F* h_R_g = new TH1F(Form("h_R_g"),Form("Radial dist. on SAM PMT region (beam gen);Radius (mm);hits/#thrownEvents/%.1fmm",bin_widthR),nbinR,x_minR,x_maxR);//gamma
   TH1F* h_R_epi = new TH1F(Form("h_R_epi"),Form("Radial dist. on SAM PMT region (beam gen);Radius (mm);hits/#thrownEvents/%.1fmm",bin_widthR),nbinR,x_minR,x_maxR);//electron and pi-
   TH1F* h_R_ppi = new TH1F(Form("h_R_ppi"),Form("Radial dist. on SAM PMT region (beam gen);Radius (mm);hits/#thrownEvents/%.1fmm",bin_widthR),nbinR,x_minR,x_maxR);//positron and pi+
   TH1F* h_R_n = new TH1F(Form("h_R_n"),Form("Radial dist. on SAM PMT region (beam gen);Radius (mm);hits/#thrownEvents/%.1fmm",bin_widthR),nbinR,x_minR,x_maxR);//neutron
   TH1F* h_R_epG1 = new TH1F(Form("h_R_epG1"),Form("Radial dist. on SAM PMT region (beam gen);Radius (mm);hits/#thrownEvents/%.1fmm",bin_widthR),nbinR,x_minR,x_maxR);//electron and positron E>1MeV
   TH1F* h_R_priG1 = new TH1F(Form("h_R_priG1"),Form("Radial dist. on SAM PMT region (beam gen);Radius (mm);hits/#thrownEvents/%.1fmm",bin_widthR),nbinR,x_minR,x_maxR);//primary E>1MeV
   double x_minE = 1.e-8;
   double x_maxE = 12000;
   double nbinE = 1000;
   double bin_widthE = (x_maxE-x_minE)/nbinE;
   TH1F* h_E_g = new TH1F(Form("h_E_g"),Form("Energy dist. on SAM PMT region (beam gen);Energy (MeV);hits/#thrownEvents/%.1fmm",bin_widthE),nbinE,x_minE,x_maxE);
   TH1F* h_E_epi = new TH1F(Form("h_E_epi"),Form("Energy dist. on SAM PMT region (beam gen);Energy (MeV);hits/#thrownEvents/%.1fmm",bin_widthE),nbinE,x_minE,x_maxE);
   TH1F* h_E_ppi = new TH1F(Form("h_E_ppi"),Form("Energy dist. on SAM PMT region (beam gen);Energy (MeV);hits/#thrownEvents/%.1fmm",bin_widthE),nbinE,x_minE,x_maxE);
   TH1F* h_E_n = new TH1F(Form("h_E_n"),Form("Energy dist. on SAM PMT region (beam gen);Energy (MeV);hits/#thrownEvents/%.1fmm",bin_widthE),nbinE,x_minE,x_maxE);
   TH1F* h_E_epG1 = new TH1F(Form("h_E_epG1"),Form("Energy dist. on SAM PMT region (beam gen);Energy (MeV);hits/#thrownEvents/%.1fmm",bin_widthE),nbinE,x_minE,x_maxE);
   TH1F* h_E_priG1 = new TH1F(Form("h_E_priG1"),Form("Energy dist. on SAM PMT region (beam gen);Energy (MeV);hits/#thrownEvents/%.1fmm",bin_widthE),nbinE,x_minE,x_maxE);
   h_R_g->SetLineColor(1);
   h_R_epi->SetLineColor(2);
   h_R_ppi->SetLineColor(3);
   h_R_n->SetLineColor(4);
   h_R_epG1->SetLineColor(6);
   h_R_priG1->SetLineColor(7);
   h_R_g->Sumw2();
   h_R_epi->Sumw2();
   h_R_ppi->Sumw2();
   h_R_n->Sumw2();
   h_R_epG1->Sumw2();
   h_R_priG1->Sumw2();
   h_E_g->SetLineColor(1);
   h_E_epi->SetLineColor(2);
   h_E_ppi->SetLineColor(3);
   h_E_n->SetLineColor(4);
   h_E_epG1->SetLineColor(6);
   h_E_priG1->SetLineColor(1);
   h_E_g->Sumw2();
   h_E_epi->Sumw2();
   h_E_ppi->Sumw2();
   h_E_n->Sumw2();
   h_E_epG1->Sumw2();
   h_E_priG1->Sumw2();

   TChain* T = new TChain("T");
   int nfileSplit=0;
   for(int fileSplit=1001;fileSplit<=1010;fileSplit++){
       nfileSplit++;
       T->Add(rootfile_dir+Form("PMTSh_beam/PMTSh_beam_%d.root",fileSplit));
   }
   cout<<Form("Found %d number of file splits!",nfileSplit)<<endl;
   Long64_t nentry = T->GetEntries();
   std::vector<remollGenericDetectorHit_t> *fHit =0;
   T->SetBranchAddress("hit", &fHit);
   
   Double_t energy, hitr, asym, rate, trid;
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
         trid = fHit->at(pk).trid;
         rate = 1;
        if(detector==28 && hitr>rmin && hitr<rmax){
          if(pid==22){
            h_R_g->Fill(hitr,rate);
            h_E_g->Fill(energy,rate);
          }else if(pid==11 || pid==-211){
            h_R_epi->Fill(hitr,rate);
            h_E_epi->Fill(energy,rate);
          }else if(pid==-11 || pid==211){
            h_R_ppi->Fill(hitr,rate);
            h_E_ppi->Fill(energy,rate);
          }else if(pid==2112){
            h_R_n->Fill(hitr,rate);
            h_E_n->Fill(energy,rate);
          }else if((pid==11 || pid==-11) && energy>1){
            h_R_epG1->Fill(hitr,rate);
            h_E_epG1->Fill(energy,rate);
          }else if(trid==1 && energy>1){
            h_R_priG1->Fill(hitr,rate);
            h_E_priG1->Fill(energy,rate);
          }
        }
      }
   }
   
   h_R_g->Scale(1.0/nentry);
   h_R_epi->Scale(1.0/nentry);
   h_R_ppi->Scale(1.0/nentry);
   h_R_n->Scale(1.0/nentry);
   h_R_epG1->Scale(1.0/nentry);
   h_R_priG1->Scale(1.0/nentry);
   h_E_g->Scale(1.0/nentry);
   h_E_epi->Scale(1.0/nentry);
   h_E_ppi->Scale(1.0/nentry);
   h_E_n->Scale(1.0/nentry);
   h_E_epG1->Scale(1.0/nentry);
   h_E_priG1->Scale(1.0/nentry);

   TLatex latex;
   latex.SetNDC(1);
   latex.SetTextSize(0.04);
   TCanvas* c_R = new TCanvas("c_R");
   gPad->SetLogy();
   h_R_g->SetTitle("Radial Distribution at SAM PMT Region (beam generator)");
   h_R_g->Draw("hist");
   h_R_epi->Draw("hist same");
   h_R_ppi->Draw("hist same");
   h_R_n->Draw("hist same");
   h_R_epG1->Draw("hist same");
   h_R_priG1->Draw("hist same");
   latex.SetTextColor(1);
   latex.DrawLatex(0.55,0.85,"gamma");
   latex.SetTextColor(2);
   latex.DrawLatex(0.55,0.80,"e-/pi-");
   latex.SetTextColor(3);
   latex.DrawLatex(0.55,0.75,"e+/pi+");
   latex.SetTextColor(4);
   latex.DrawLatex(0.70,0.85,"neutron");
   latex.SetTextColor(6);
   latex.DrawLatex(0.70,0.80,"e-/e+ E>1 MeV");
   latex.SetTextColor(7);
   latex.DrawLatex(0.70,0.75,"primary E>1 MeV");
   c_R->SaveAs("./temp/sam_pmt_radial.pdf");

   TCanvas* c_E = new TCanvas("c_E");
   gPad->SetLogy();
   gPad->SetLogx();
   h_E_g->Draw("prof");
   h_E_epi->Draw("prof same");
   h_E_ppi->Draw("prof same");
   h_E_n->Draw("prof same");
   h_E_epG1->Draw("prof same");
   h_E_priG1->Draw("prof same");
   latex.SetTextColor(1);
   latex.DrawLatex(0.55,0.85,"gamma");
   latex.SetTextColor(2);
   latex.DrawLatex(0.55,0.80,"e-/pi-");
   latex.SetTextColor(3);
   latex.DrawLatex(0.55,0.75,"e+/pi+");
   latex.SetTextColor(4);
   latex.DrawLatex(0.70,0.85,"neutron");
   latex.SetTextColor(6);
   latex.DrawLatex(0.70,0.80,"e-/e+ E>1 MeV");
   latex.SetTextColor(7);
   latex.DrawLatex(0.70,0.75,"primary E>1 MeV");
   c_E->SaveAs("./temp/sam_pmt_energy.pdf");

   gSystem->Exec(Form("pdfunite ./temp/sam_pmt_*.pdf ./plots/PMTSh_beam_sam_pmt_bkg.pdf"));
   gSystem->Exec(Form("rm -rf ./temp/sam_pmt_*.pdf"));
}
