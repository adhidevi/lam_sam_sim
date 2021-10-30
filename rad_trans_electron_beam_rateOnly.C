//plot radial and transverse hit distribution on ring5, lam, usscanner, and dsscanner plane
//
#include "remolltypes.hh"
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>

void rad_trans_electron_beam_rateOnly(){
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
   int beam_curr = 65;

   TString rootfile_dir = "$VOLATILE/remoll_rootfiles/PhotonBlocker";
   
   double x_min1D[nDet] = {500,500,0};
   double x_max1D[nDet] = {1500,1500,100};
   double x_min2D[nDet] = {-2000,-2000,-175};
   double x_max2D[nDet] = {2000,2000,175};
   double eCut[nDet] = {1000,1,1};
   double rCut[nDet] = {500,500,0};

   TH1F* h_rate[nDet];
   TH2F* hxy_rate[nDet];
   TH2F* hxy_rateQ[nDet];


   for(int iDet=0;iDet<nDet;iDet++){
      int nbin1D = 500;
      int nbin2D = 400;
      h_rate[iDet] = new TH1F(Form("h_rate[%d]",iDet),Form("%s dist. on %s det plane;Radius (mm);hits/%sthrownEvents","Radial",DetName[iDet].Data(),"#"),nbin1D,x_min1D[iDet],x_max1D[iDet]);
      hxy_rate[iDet] = new TH2F(Form("hxy_rate[%d]",iDet),Form("%s dist. on %s det plane;x (mm);y (mm); hits/%sthrownEvents","Transverse",DetName[iDet].Data(),"#"),nbin2D,x_min2D[iDet],x_max2D[iDet],nbin2D,x_min2D[iDet],x_max2D[iDet]);
      hxy_rateQ[iDet] = new TH2F(Form("hxy_rateQ[%d]",iDet),Form("%s dist. on %s det plane;x (mm);y (mm); hits/%sthrownEvents","Transverse",DetName[iDet].Data(),"#"),nbin2D,x_min2D[iDet],x_max2D[iDet],nbin2D,x_min2D[iDet],x_max2D[iDet]);

      h_rate[iDet]->SetLineColor(1);
      h_rate[iDet]->Sumw2();
   }
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
   
   Double_t energy, hitr, rate, hitx, hity;
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
         rate = 1;//Rate is taken 1 for beam generator
         for(int iDet=0;iDet<nDet;iDet++){
           if(detector==DetNo[iDet] && energy>eCut[iDet] && hitr>rCut[iDet] && pid==11){
              h_rate[iDet]->Fill(hitr,rate);
              hxy_rate[iDet]->Fill(hitx,hity,rate);
	     if(hitr>rmin[iDet] && hitr<rmax[iDet]){
              hxy_rateQ[iDet]->Fill(hitx,hity,rate);
             }
           }
        }
      }
   }
   
   TH1F* h_rateQ[nDet];
   for(int iDet=0;iDet<nDet;iDet++){
      h_rate[iDet]->Scale(1./nentry);
      hxy_rate[iDet]->Scale(1./nentry);
      hxy_rateQ[iDet]->Scale(1./nentry);

      h_rateQ[iDet] = (TH1F*)h_rate[iDet]->Clone(Form("h_rateQ[%d]",iDet));
      h_rateQ[iDet]->GetXaxis()->SetRangeUser(rmin[iDet],rmax[iDet]);
      h_rateQ[iDet]->SetLineColor(4);
    }

   TLine* line[2];
   TLatex latex;
   latex.SetNDC(1);
   latex.SetTextSize(0.04);

   TCanvas* c_rate_linear[nDet];
   for(int iDet=0;iDet<nDet;iDet++){
      c_rate_linear[iDet] = new TCanvas(Form("c_rate_linear_%d",iDet));
      h_rate[iDet]->Draw("hist");
      h_rateQ[iDet]->Draw("same hist");
      line[0] = new TLine(rmin[iDet],0.0,rmin[iDet],h_rate[iDet]->GetMaximum()/2.0);
      line[1] = new TLine(rmax[iDet],0.0,rmax[iDet],h_rate[iDet]->GetMaximum()/2.0);
      for(int ii=0;ii<2;ii++){
         line[ii]->SetLineWidth(2);
         line[ii]->SetLineColor(kRed-2);
         line[ii]->Draw();
      }
      latex.SetTextColor(kRed-2);
      latex.DrawLatex(0.70,0.70,Form("%s",Form("R_{min}=%.0f mm",rmin[iDet])));
      latex.DrawLatex(0.70,0.65,Form("%s",Form("R_{max}=%.0f mm",rmax[iDet])));
      latex.SetTextColor(1);
      latex.DrawLatex(0.45,0.85,Form("%s","beam"));
      latex.SetTextColor(4);
      double integral = h_rateQ[iDet]->Integral();
      latex.DrawLatex(0.55,0.85,Form("%s:%.3e","accept",integral));
      if(iDet==2){
      latex.SetTextColor(1);
      latex.DrawLatex(0.45,0.80,Form("%s:%.3e","accept(60<=r<=75)",h_rate[iDet]->Integral(h_rate[iDet]->FindBin(60),h_rate[iDet]->FindBin(75))));
      }
      c_rate_linear[iDet]->SaveAs(Form("./temp/p_rate_linear_%d.pdf",iDet));
   }

   TCanvas* c_rate_log[nDet];
   for(int iDet=0;iDet<nDet;iDet++){
      c_rate_log[iDet] = new TCanvas(Form("c_rate_log_%d",iDet));
      gPad->SetLogy();
      h_rate[iDet]->Draw("hist");
      h_rateQ[iDet]->Draw("same hist");
      line[0] = new TLine(rmin[iDet],0.0,rmin[iDet],h_rate[iDet]->GetMaximum()/2.0);
      line[1] = new TLine(rmax[iDet],0.0,rmax[iDet],h_rate[iDet]->GetMaximum()/2.0);
      for(int ii=0;ii<2;ii++){
         line[ii]->SetLineWidth(2);
         line[ii]->SetLineColor(kRed-2);
         line[ii]->Draw();
      }
      latex.SetTextColor(kRed-2);
      latex.DrawLatex(0.70,0.70,Form("%s",Form("R_{min}=%.0f mm",rmin[iDet])));
      latex.DrawLatex(0.70,0.65,Form("%s",Form("R_{max}=%.0f mm",rmax[iDet])));
      latex.SetTextColor(1);
      latex.DrawLatex(0.45,0.85,Form("%s","beam"));
      latex.SetTextColor(4);
      double integral = h_rateQ[iDet]->Integral();
      latex.DrawLatex(0.55,0.85,Form("%s:%.3e","accept",integral));
      if(iDet==2){
      latex.SetTextColor(1);
      latex.DrawLatex(0.45,0.80,Form("%s:%.3e","accept(60<=r<=75)",h_rate[iDet]->Integral(h_rate[iDet]->FindBin(60),h_rate[iDet]->FindBin(75))));
      }
      c_rate_log[iDet]->SaveAs(Form("./temp/p_rate_log_%d.pdf",iDet));
   }

   TCanvas* cxy_rate[nDet];
   for(int iDet=0;iDet<nDet;iDet++){
      cxy_rate[iDet] = new TCanvas(Form("cxy_rate%d",iDet));
      if(gPad) gPad->SetRightMargin(0.15);
      hxy_rate[iDet]->Draw("colz");
      TArc* arc[2];
      arc[0] = new TArc(0.0,0.0,rmin[iDet],0,360);
      arc[1] = new TArc(0.0,0.0,rmax[iDet],0,360);
      for(int ii=0;ii<2;ii++){
         arc[ii]->SetLineWidth(2);
         arc[ii]->SetLineColor(kRed-2);
         arc[ii]->SetFillStyle(0);
         arc[ii]->SetLineStyle(7);
         arc[ii]->Draw();
      }
      latex.SetTextColor(kRed-2);
      latex.DrawLatex(0.20,0.85,Form("%s",Form("R_{min}=%.0f mm",rmin[iDet])));
      latex.DrawLatex(0.20,0.80,Form("%s",Form("R_{max}=%.0f mm",rmax[iDet])));
      latex.SetTextColor(1);
      latex.DrawLatex(0.45,0.85,Form("%s","beam"));
      latex.SetTextColor(4);
      double integral = hxy_rateQ[iDet]->Integral();
      latex.DrawLatex(0.55,0.85,Form("%s:%.3e","accept",integral));
      cxy_rate[iDet]->SaveAs(Form("./temp/p_rate_xy_%d.pdf",iDet));
   }
//Now combine all pdf files saved in ./temp/ directory and save a single pdf file in ./plots/ directory
   gSystem->Exec(Form("pdfunite ./temp/*.pdf ./plots/rad_trans_dist_electron_beam_rateOnly_test.pdf"));
   gSystem->Exec(Form("rm -rf ./temp/*.pdf"));
}
