//plot energy distribution on ring5, lam, and sam plane
//
#include "remolltypes.hh"
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>

void energy_electrons_beam(){
   gROOT->Reset();
   gStyle->SetOptStat(111111);
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
   
   double x_min1D[nDet] = {0,0,0};
   double x_max1D[nDet] = {8000,4000,12000};
   double eCut[nDet] = {0,0,0};
   TH1F* h_energy[nDet];
   
   for(int iDet=0;iDet<nDet;iDet++){
      int nbin1D = 200;
      h_energy[iDet] = new TH1F(Form("h_energy[%d]",iDet),Form("%s dist. on %s det plane (beam generator);Energy (MeV);hits/%sthrownEvents","Energy",DetName[iDet].Data(),"#"),nbin1D,x_min1D[iDet],x_max1D[iDet]);

      h_energy[iDet]->SetLineColor(1);
      h_energy[iDet]->Sumw2();
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
   remollEvent_t *fEv =0;
   Double_t fRate=0.;
   T->SetBranchAddress("hit", &fHit);
   T->SetBranchAddress("ev", &fEv);
   T->SetBranchAddress("rate", &fRate);
   
   Double_t energy, hitr, hitx, hity, rate;
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
           if(detector==DetNo[iDet] && energy>eCut[iDet] && (hitr>rmin[iDet] && hitr<rmax[iDet]) && pid==11){
              h_energy[iDet]->Fill(energy,rate);
           }
        }
      }
   }
   
   TH1F* h_energyQ[nDet];
   for(int iDet=0;iDet<nDet;iDet++){
      h_energy[iDet]->Scale(1./nentry);
    }

   TLatex latex;
   latex.SetNDC(1);
   latex.SetTextSize(0.04);

   TCanvas* c_rate_linear[nDet];
   for(int iDet=0;iDet<nDet;iDet++){
      c_rate_linear[iDet] = new TCanvas(Form("c_rate_linear_%d",iDet));
      h_energy[iDet]->Draw("hist");
      latex.SetTextColor(kRed-2);
      latex.DrawLatex(0.15,0.85,Form("Integral = %.2e",h_energy[iDet]->Integral()));
      latex.DrawLatex(0.15,0.80,Form("Cut::E>%.0f && r>%.0f && r<%.0f && PID==11",eCut[iDet],rmin[iDet],rmax[iDet]));
      gPad->Update();
      TPaveStats* st = (TPaveStats*)h_energy[iDet]->FindObject("stats");
      st->SetFillStyle(0);
      c_rate_linear[iDet]->SaveAs(Form("./temp/p_rate_linear_%d.pdf",iDet));
   }

   TCanvas* c_rate_log[nDet];
   for(int iDet=0;iDet<nDet;iDet++){
      c_rate_log[iDet] = new TCanvas(Form("c_rate_log_%d",iDet));
      gPad->SetLogy(1);
      h_energy[iDet]->Draw("hist");
      latex.SetTextColor(kRed-2);
      latex.DrawLatex(0.15,0.80,Form("Integral = %.2e",h_energy[iDet]->Integral()));
      latex.DrawLatex(0.15,0.85,Form("Cut::E>%.0f && r>%.0f && r<%.0f && PID==11",eCut[iDet],rmin[iDet],rmax[iDet]));
      gPad->Update();
      TPaveStats* st = (TPaveStats*)h_energy[iDet]->FindObject("stats");
      st->SetFillStyle(0);
      c_rate_log[iDet]->SaveAs(Form("./temp/p_rate_log_%d.pdf",iDet));
   }

//Now combine all pdf files saved in ./temp/ directory and save a single pdf file in ./plots/ directory
   gSystem->Exec(Form("pdfunite ./temp/*.pdf ./plots/energyG0_dist_electron_beam.pdf"));
   gSystem->Exec(Form("rm -rf ./temp/*.pdf"));
}
