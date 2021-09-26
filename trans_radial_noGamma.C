//plot y vs x (transverse) looking downstream from the scattering chamber hit distribution on  plane
//
#include "remolltypes.hh"
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>

void trans_radial_noGamma(){
   gROOT->Reset();
   gStyle->SetOptStat(0);
   gStyle->SetTitleYOffset(1.3);
   gStyle->SetPadGridX(1);
   gStyle->SetPadGridY(1);
   TGaxis::SetMaxDigits(3);

   TString rootgenerator_dir = "$VOLATILE/remoll_rootfiles/PhotonBlocker/";
   TString generator[] = {"sim_ee_lam","sim_epel_lam","sim_epinel_lam"};
   int ngenerator = sizeof(generator)/sizeof(*generator);

   TString sim[] = {"ee", "epel", "epinel"};

   int detNum[] = {174,28};
   int ndet = sizeof(detNum)/sizeof(*detNum);
   TFile* outfile = TFile::Open("outrootfiles/sim_lam_noGamma_LAMrG1013mm.root","RECREATE"); 
   for(int idet=0;idet<ndet;idet++)
   outfile->mkdir(Form("det%d",detNum[idet]));
   
   double beam_curr = 65.0;//uA
   TH2F* trans_r[ndet][ngenerator];
   TH2F* trans_rE[ndet][ngenerator];
   TH2F* trans_rA[ndet][ngenerator];
   TH1F* radial_r[ndet][ngenerator];
   TH1F* radial_rE[ndet][ngenerator];
   TH1F* radial_rA[ndet][ngenerator];
   for(int igenerator=0;igenerator<ngenerator;igenerator++){
     for(int idet=0;idet<ndet;idet++){
      trans_r[idet][igenerator] = new TH2F(Form("det%d_r_%s",detNum[idet],sim[igenerator].Data()),Form("%s on det%d (wt: rate/%.0fuA); x (mm);y (mm)",sim[igenerator].Data(),detNum[idet],beam_curr),400,-2000,2000,400,-2000,2000);
      trans_rE[idet][igenerator] = new TH2F(Form("det%d_rE_%s",detNum[idet],sim[igenerator].Data()),Form("%s on det%d (wt: rate*E/%.0fuA); x (mm);y (mm)",sim[igenerator].Data(),detNum[idet],beam_curr),400,-2000,2000,400,-2000,2000);
      trans_rA[idet][igenerator] = new TH2F(Form("det%d_rA_%s",detNum[idet],sim[igenerator].Data()),Form("%s on det%d (wt: rate*A/%.0fuA); x (mm);y (mm)",sim[igenerator].Data(),detNum[idet],beam_curr),400,-2000,2000,400,-2000,2000);
      radial_r[idet][igenerator] = new TH1F(Form("radial_det%d_r_%s",detNum[idet],sim[igenerator].Data()),Form("%s on det%d (wt: rate/%.0fuA); r (mm);rate (GHz/%.0fuA)",sim[igenerator].Data(),detNum[idet],beam_curr,beam_curr),500,0,4000);
      radial_rE[idet][igenerator] = new TH1F(Form("radial_det%d_rE_%s",detNum[idet],sim[igenerator].Data()),Form("%s on det%d (wt: rate*E/%.0fuA); r (mm);rate*E (GHz*GeV/%.0fuA)",sim[igenerator].Data(),detNum[idet],beam_curr,beam_curr),500,0,4000);
      radial_rA[idet][igenerator] = new TH1F(Form("radial_det%d_rA_%s",detNum[idet],sim[igenerator].Data()),Form("%s on det%d (wt: rate*A/%.0fuA); r (mm);rate*A (GHz*ppb/%.0fuA)",sim[igenerator].Data(),detNum[idet],beam_curr,beam_curr),500,0,4000);
      }
    }

   for(int igenerator=0;igenerator<ngenerator;igenerator++){
      TChain* T = new TChain("T");
      T->Add(rootgenerator_dir+Form("%s/%s_*.root",generator[igenerator].Data(),generator[igenerator].Data()));
      Long64_t nentry = T->GetEntries();
      std::vector<remollGenericDetectorHit_t> *fHit =0;
      remollEvent_t *fEv =0;
      Double_t fRate=0.;
      T->SetBranchAddress("hit", &fHit);
      T->SetBranchAddress("ev", &fEv);
      T->SetBranchAddress("rate", &fRate);
   
      Float_t energy(-1.e-12), hitr(-1.e-12), detector(-1.e-12), hitx(-1.e-12), hity(-1.e-12), phi(-1.e-12), rate(1.e-12), asym(-1.e-12);
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
            hitx = fHit->at(pk).x;
            hity = fHit->at(pk).y;
	    trid = fHit->at(pk).trid;
	    mtrid = fHit->at(pk).mtrid;
            phi = fHit->at(pk).ph;
	    asym = -1*fEv->A;
	    rate = fRate/1.e9/20.0;
            if(detector==detNum[0] && hitr>1013 && energy>1 && pid!=22){
              trans_r[0][igenerator]->Fill(hitx,hity,rate);
              trans_rE[0][igenerator]->Fill(hitx,hity,rate*energy/1.e3);
              trans_rA[0][igenerator]->Fill(hitx,hity,rate*asym);
              radial_r[0][igenerator]->Fill(hitr,rate);
              radial_rE[0][igenerator]->Fill(hitr,rate*energy/1.e3);
              radial_rA[0][igenerator]->Fill(hitr,rate*asym);
            }
            if(detector==detNum[1] && hitr>500 && energy>1000 && pid!=22){
              trans_r[1][igenerator]->Fill(hitx,hity,rate);
              trans_rE[1][igenerator]->Fill(hitx,hity,rate*energy/1.e3);
              trans_rA[1][igenerator]->Fill(hitx,hity,rate*asym);
              radial_r[1][igenerator]->Fill(hitr,rate);
              radial_rE[1][igenerator]->Fill(hitr,rate*energy/1.e3);
              radial_rA[1][igenerator]->Fill(hitr,rate*asym);
            }
         }
      }
  }
    for(int idet=0;idet<ndet;idet++){
       outfile->cd(Form("det%d",detNum[idet]));
       for(int igenerator=0;igenerator<ngenerator;igenerator++){
          trans_r[idet][igenerator]->Write();
          trans_rE[idet][igenerator]->Write();
          trans_rA[idet][igenerator]->Write();
          radial_r[idet][igenerator]->Write();
          radial_rE[idet][igenerator]->Write();
          radial_rA[idet][igenerator]->Write();
       }
    }   
}
