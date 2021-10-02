//plot energy distribution on ring5, lam, and sam plane
//
#include "../remolltypes.hh"
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>

void energy_radius_electron_beam(){
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
   double e_min[nDet] = {0,0,0};
   double e_max[nDet] = {8000,4000,12000};
   double eCut[nDet] = {1000,1,1};

   int beam_curr = 65;

   TString rootfile_dir = "$VOLATILE/remoll_rootfiles/PhotonBlocker";
   
   TH2F* h_er[nDet];
   
   for(int iDet=0;iDet<nDet;iDet++){
      int nbinx = 500;
      int nbiny = 200;
      h_er[iDet] = new TH2F(Form("h_er[%d]",iDet),Form("%s electron dist. on %s det plane (beam generator);Radius (mm); Energy (MeV);hits/%sthrownEvents","Energy vs Radius",DetName[iDet].Data(),"#"),nbinx,rmin[iDet],rmax[iDet],nbiny,e_min[iDet],e_max[iDet]);
   }
   TChain* T = new TChain("T");
   int nfile = 0;
   for(int ifile = 1001;ifile<=1010;ifile++){
      nfile++;
      T->Add(Form("%s/pB_beam/pB_beam_%d.root",rootfile_dir.Data(),ifile));
   }
   cout<<Form("Found %d file splits!!!",nfile)<<endl;

   Long64_t nentry = T->GetEntries();
   std::vector<remollGenericDetectorHit_t> *fHit =0;
   T->SetBranchAddress("hit", &fHit);
   
   Double_t energy, hitr, rate;
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
         rate = 1;//Rate is taken 1 for beam generator
         for(int iDet=0;iDet<nDet;iDet++){
           if(detector==DetNo[iDet] && energy>eCut[iDet] && (hitr>rmin[iDet] && hitr<rmax[iDet]) && pid==11){
              h_er[iDet]->Fill(hitr,energy,rate);
           }
        }
      }
   }
   
   for(int iDet=0;iDet<nDet;iDet++){
      h_er[iDet]->Scale(1./nentry);
    }

   TCanvas* c_rate_linear[nDet];
   for(int iDet=0;iDet<nDet;iDet++){
      c_rate_linear[iDet] = new TCanvas(Form("c_rate_linear_%d",iDet));
      h_er[iDet]->Draw();
      c_rate_linear[iDet]->SaveAs(Form("./temp/p_rate_linear_%d.pdf",iDet));
   }

//Now combine all pdf files saved in ./temp/ directory and save a single pdf file in ./plots/ directory
   gSystem->Exec(Form("pdfunite ./temp/*.pdf ./plots/energyG1vsRadial_dist_electron_beam.pdf"));
   gSystem->Exec(Form("rm -rf ./temp/*.pdf"));
}
