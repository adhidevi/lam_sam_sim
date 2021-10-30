//
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
  for (int i = 0; i <= bins; i++){
    new_bins[i] = pow(10, from + i * width);
  }
  axis->Set(bins, new_bins);
  delete[] new_bins;
}

void sam_pmt_energy_dist_beam(){
   gROOT->Reset();
   gStyle->SetOptStat(0);
   gStyle->SetTitleYOffset(1.3);
   gStyle->SetPadGridX(1);
   gStyle->SetPadGridY(1);
   TGaxis::SetMaxDigits(3);

   Double_t rmin = 290;
   Double_t rmax = 480;
   int beam_curr = 65;

   TString rootfile_dir = "/volatile/halla/parity/adhidevi/remoll_rootfiles/PMTShielding/";
   
   double x_minR = 270;
   double x_maxR = 500;
   double nbinR = 121;
   double bin_widthR = (x_maxR-x_minR)/nbinR;
   int color[] = {1,2,3,4,6,7};
   const int nhist = sizeof(color)/sizeof(*color);
   TString particle_type[nhist] = {"gamma","e-/pi-","e+/pi+","neutron","e-/e+ E>1 MeV","primary E>1 MeV"};
   TH1F* h_R[nhist];
   TH1F* h_E[nhist];

   for(int ihist=0;ihist<nhist;ihist++){
      h_R[ihist] = new TH1F(Form("h_R%d",ihist),Form("Radial dist. on SAM PMT region (beam gen);Radius (mm);hits/#thrownEvents/%.1fmm",bin_widthR),nbinR,x_minR,x_maxR);
      h_E[ihist] = new TH1F(Form("h_E%d",ihist),Form("Energy dist. on SAM PMT region (beam gen);Energy (MeV);hits/#thrownEvents"),121,-8,5);
      h_R[ihist]->SetLineColor(color[ihist]);
      h_E[ihist]->SetLineColor(color[ihist]);
      h_R[ihist]->Sumw2();
      h_E[ihist]->Sumw2();
      niceLogXBins(h_E[ihist]);
   }
   TChain* T = new TChain("T");
   int nfileSplit=0;
   for(int fileSplit=1001;fileSplit<=2000;fileSplit++){
       TString infile = Form("%sPMTSh_beam/PMTSh_beam_%d.root",rootfile_dir.Data(),fileSplit);
       ifstream fileName(infile);
       if(!fileName){
         cout<<Form("File %s doesn't exist. Escapping this...!!",infile.Data())<<endl;
         continue;
       }
       nfileSplit++;
       T->Add(infile);
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
         energy = fHit->at(pk).p;
         hitr = fHit->at(pk).r;
         trid = fHit->at(pk).trid;
         rate = 1;
        if(detector==176 && hitr>rmin && hitr<rmax){
          if(pid==22){
            h_R[0]->Fill(hitr,rate);
            h_E[0]->Fill(energy,rate);
          }
          if(pid==11 || pid==-211){
            h_R[1]->Fill(hitr,rate);
            h_E[1]->Fill(energy,rate);
          }
          if(pid==-11 || pid==211){
            h_R[2]->Fill(hitr,rate);
            h_E[2]->Fill(energy,rate);
          }
          if(pid==2112){
            h_R[3]->Fill(hitr,rate);
            h_E[3]->Fill(energy,rate);
          }
          if((pid==11 || pid==-11) && energy>1){
            h_R[4]->Fill(hitr,rate);
            h_E[4]->Fill(energy,rate);
          }
          if(trid==1 && energy>1){
            h_R[5]->Fill(hitr,rate);
            h_E[5]->Fill(energy,rate);
          }
        }
      }
   }
   for(int ihist=0;ihist<nhist;ihist++){
      h_R[ihist]->Scale(1.0/nentry);
      h_E[ihist]->Scale(1.0/nentry);
   }
   TLatex latex;
   latex.SetNDC(1);
   latex.SetTextSize(0.04);
   TCanvas* c_R = new TCanvas("c_R");
   gPad->SetLogy();
   for(int ihist=0;ihist<nhist;ihist++){
      h_R[ihist]->SetTitle("Radial Distribution at SAM PMT Region (beam generator)");
      h_R[ihist]->Draw("hist same");
      h_R[ihist]->GetYaxis()->SetRangeUser(1.e-9,1.e-2);
      latex.SetTextColor(color[ihist]);
      if(ihist<3)
        latex.DrawLatex(0.33,0.85-0.05*ihist,Form("%s: %.2e",particle_type[ihist].Data(),h_R[ihist]->Integral()));
      else
        latex.DrawLatex(0.56,0.85-0.05*(ihist-3),Form("%s: %.2e",particle_type[ihist].Data(),h_R[ihist]->Integral()));
   }
   latex.SetTextColor(kRed+3);
   latex.DrawLatex(0.12,0.85,"Radial Cut::");
   latex.DrawLatex(0.12,0.80,Form("(%.0f,%.0f)mm",rmin,rmax));
   c_R->SaveAs("./temp/sam_pmt_mom_radial.pdf");

   TCanvas* c_E = new TCanvas("c_E");
   gPad->SetLogy();
   gPad->SetLogx();
   for(int ihist=0;ihist<nhist;ihist++){
      h_E[ihist]->SetTitle("Energy Distribution at SAM PMT Region (beam generator)");
      h_E[ihist]->Draw("hist same");
      h_E[ihist]->GetYaxis()->SetRangeUser(1.e-9,1.e-2);
      latex.SetTextColor(color[ihist]);
      if(ihist<3)
        latex.DrawLatex(0.33,0.85-0.05*ihist,Form("%s: %.2e",particle_type[ihist].Data(),h_E[ihist]->Integral()));
      else
        latex.DrawLatex(0.56,0.85-0.05*(ihist-3),Form("%s: %.2e",particle_type[ihist].Data(),h_E[ihist]->Integral()));
   }
   latex.SetTextColor(kRed+3);
   latex.DrawLatex(0.12,0.85,"Radial Cut::");
   latex.DrawLatex(0.12,0.80,Form("(%.0f,%.0f)mm",rmin,rmax));
   c_E->SaveAs("./temp/sam_pmt_mom.pdf");

   gSystem->Exec(Form("pdfunite ./temp/sam_pmt_mom*.pdf ./plots/PMTSh_beam_sam_pmt_mom.pdf"));
   gSystem->Exec(Form("rm -rf ./temp/sam_pmt_mom*.pdf"));

}
