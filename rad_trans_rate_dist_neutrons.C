//plot radial neutrons hit distribution on SAM plane
//
#include "remolltypes.hh"
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>

void rad_trans_rate_dist_neutrons(){
   gROOT->Reset();
   gStyle->SetOptStat(0);
   gStyle->SetTitleYOffset(1.3);
   gStyle->SetPadGridX(1);
   gStyle->SetPadGridY(1);
   TGaxis::SetMaxDigits(3);

   int DetNo[] = {176,176,176,176,176};
   const int nDet = sizeof(DetNo)/sizeof(*DetNo);
   TString DetName[] = {"Total","Scanner","PMT+Base"," ","SAM"};
   Double_t rmin[nDet] = {0,500,294,60,45};
   Double_t rmax[nDet] = {1500,700,472,75,60};
   int nbin = 200;
   int beam_curr = 65;

   TString rootfile_dir = "/volatile/halla/parity/adhidevi/remoll_rootfiles/PMTShielding";
   
   int colorF[nDet] = {1,7,3,4,6};
   int colorB[nDet] = {2,8,9,41,46};

   TH1F* h_forward[nDet];
   TH1F* h_backward[nDet];

   for(int iDet=0;iDet<nDet;iDet++){
      h_forward[iDet] = new TH1F(Form("h_forward[%d]",iDet),Form("%s dist. on %s det plane (beam gen);Radius (mm);hits/%sthrownEvents","Neutron","SAM","#"),nbin,0,1500);
      h_backward[iDet] = new TH1F(Form("h_backward[%d]",iDet),Form("%s dist. on %s det plane (beam gen);Radius (mm);hits/%sthrownEvents","Neutron","SAM","#"),nbin,0,1500);

      h_forward[iDet]->SetLineColor(colorF[iDet]);
      h_backward[iDet]->SetLineColor(colorB[iDet]);
      h_forward[iDet]->Sumw2();
      h_backward[iDet]->Sumw2();
   }
   TChain* T = new TChain("T");
   int nfile = 0;
   for(int ifile = 1001;ifile<=2000;ifile++){
      string infile = (Form("%s/PMTSh_beam_V3/PMTSh_beam_V3_%d.root",rootfile_dir.Data(),ifile));
      ifstream filename(infile);
      if(!filename) continue;
      nfile++;
      T->Add(Form("%s",infile.c_str()));
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
	 if(std::isnan(rate) || std::isinf(rate)) continue;
         pid = (Int_t)fHit->at(pk).pid;
         detector = fHit->at(pk).det;
         energy = fHit->at(pk).e;
         hitr = fHit->at(pk).r;
         hitpz = fHit->at(pk).pz;
         rate = 1;//Rate is taken 1 for beam generator
         for(int iDet=0;iDet<nDet;iDet++){
           if(detector==DetNo[iDet] && (hitr>=rmin[iDet] && hitr<=rmax[iDet]) && pid==2112){
	    if(hitpz>=0)
             h_forward[iDet]->Fill(hitr,rate);
	    else
             h_backward[iDet]->Fill(hitr,rate);
           }
        }
      }
   }
   
   TLatex latex;
   latex.SetNDC(1);
   latex.SetTextSize(0.04);

   TCanvas* c_rate_linear= new TCanvas("c_rate_linear");
   for(int iDet=0;iDet<nDet;iDet++){
      h_forward[iDet]->Scale(1./nentry);
      h_backward[iDet]->Scale(1./nentry);
   }
   h_backward[0]->Draw("hist");
   h_backward[0]->GetYaxis()->SetRangeUser(1.e-7,1.e-5);
   h_forward[0]->Draw("hist same");
   for(int iDet=1;iDet<nDet;iDet++){
//      h_forward[iDet]->Draw("hist same");
//      h_backward[iDet]->Draw("hist same");

      TString intRange = Form("%.0f<r<%.0f",rmin[iDet],rmax[iDet]);
      double integralF = h_forward[iDet]->Integral();
      double integralB = h_backward[iDet]->Integral();
      latex.SetTextColor(colorF[0]);
      latex.DrawLatex(0.10,0.90-0.05*iDet,Form("%s(%s):%.3e","Forward",intRange.Data(),integralF));
      latex.SetTextColor(colorB[0]);
      latex.DrawLatex(0.50,0.90-0.05*iDet,Form("%s(%s):%.3e","Backward",intRange.Data(),integralB));
   }
      c_rate_linear->SaveAs("./temp/neutrons_linear.pdf");

   TCanvas* c_rate_log= new TCanvas("c_rate_log");
   gPad->SetLogy();
   h_backward[0]->Draw("hist");
   h_forward[0]->Draw("hist same");
   for(int iDet=1;iDet<nDet;iDet++){
//      h_forward[iDet]->Draw("hist same");
//      h_backward[iDet]->Draw("hist same");

      TString intRange = Form("%.0f<r<%.0f",rmin[iDet],rmax[iDet]);
      double integralF = h_forward[iDet]->Integral();
      double integralB = h_backward[iDet]->Integral();
      latex.SetTextColor(colorF[0]);
      latex.DrawLatex(0.10,0.90-0.05*iDet,Form("%s(%s):%.3e","Forward",intRange.Data(),integralF));
      latex.SetTextColor(colorB[0]);
      latex.DrawLatex(0.50,0.90-0.05*iDet,Form("%s(%s):%.3e","Backward",intRange.Data(),integralB));
   }
      c_rate_log->SaveAs("./temp/neutrons_log.pdf");

//Now combine all pdf files saved in ./temp/ directory and save a single pdf file in ./plots/ directory
   gSystem->Exec(Form("pdfunite ./temp/neutrons_*.pdf ./plots/PMTSh_beam_V3_FB_neutrons_allE_det176.pdf"));
   gSystem->Exec(Form("rm -rf ./temp/neutrons_*.pdf"));
}
