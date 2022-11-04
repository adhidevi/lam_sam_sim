#include "remolltypes.hh"
#include <sstream>
#include <iostream>
#include <fstream>

const double pi = TMath::Pi();
TString rootfile_dir = "/volatile/halla/parity/adhidevi/remoll_rootfiles/develop_br";
const string spTit[] = {"e-/#pi-","e+/#pi+","#gamma","neutron","e-/e+","e- trid==1"};
const int nSp = sizeof(spTit)/sizeof(*spTit);
const string spH[nSp] ={"epiM","epiP","g","n","ee","eTrIdCut"};
map<int,int> spM {{11,1},{-211,1},{-11,2},{211,2},{22,3},{2112,4}};
const string detH[] = {"det28","det1761","det1762","det1763","det1764"};
const int nDet = sizeof(detH)/sizeof(*detH);
map<int,int> dtM {{28,1},{1761,2},{1762,3},{1763,4},{1764,5}};
int Det[nDet] = {28,1761,1762,1763,1764};

const double kinEcut[] = {0,1};//MeV
const int nkinEcut = sizeof(kinEcut)/sizeof(*kinEcut);

TH1D* hR[nSp][nDet][nkinEcut];
TH1D* hR_pzG0[nSp][nDet][nkinEcut];
TH1D* hR_pzL0[nSp][nDet][nkinEcut];
TH2D* hXY[nSp][nDet][nkinEcut];
TH2D* hXY_pzG0[nSp][nDet][nkinEcut];
TH2D* hXY_pzL0[nSp][nDet][nkinEcut];

TH1D* hR_E[nSp][nDet][nkinEcut];//energy weighted
TH1D* hR_E_pzG0[nSp][nDet][nkinEcut];//energy weighted
TH1D* hR_E_pzL0[nSp][nDet][nkinEcut];//energy weighted
TH2D* hXY_E[nSp][nDet][nkinEcut];//energy weighted
TH2D* hXY_E_pzG0[nSp][nDet][nkinEcut];//energy weighted
TH2D* hXY_E_pzL0[nSp][nDet][nkinEcut];//energy weighted

TFile* outfile;
int beamGen(1);
const string tgt_gen_config = "LH2_beam_V29";
string histR_title;
string histXY_title;
string histR_E_title;
string histXY_E_title;

void sam_radiation_study_rad_xy(){
   gStyle->SetOptStat(0);

   for(int iSp=0;iSp<nSp;iSp++){
    for(int iDet=0;iDet<nDet;iDet++){
     for(int ikinEcut=0;ikinEcut<nkinEcut;ikinEcut++){
      if(beamGen){
      histR_title = Form("%s (KE>%.0f MeV) Radial dist. on %s (%s);Radius (mm);hits/%sthrownEvents",spTit[iSp].c_str(),kinEcut[ikinEcut],detH[iDet].c_str(),tgt_gen_config.c_str(),"#");
      histXY_title = Form("%s (KE>%.0f MeV) XY dist. on %s (%s);x (mm);y (mm);hits/%sthrownEvents",spTit[iSp].c_str(),kinEcut[ikinEcut],detH[iDet].c_str(),tgt_gen_config.c_str(),"#");
      histR_E_title = Form("%s (KE>%.0f MeV) Energy Weighted Radial dist. on %s (%s);Radius (mm);hits/%sthrownEvents",spTit[iSp].c_str(),kinEcut[ikinEcut],detH[iDet].c_str(),tgt_gen_config.c_str(),"#");
      histXY_E_title = Form("%s (KE>%.0f MeV) Energy Weighted XY dist. on %s (%s);x (mm);y (mm);hits/%sthrownEvents",spTit[iSp].c_str(),kinEcut[ikinEcut],detH[iDet].c_str(),tgt_gen_config.c_str(),"#");
      }else{
      histR_title = Form("%s (KE>%.0f MeV) Radial dist. on %s (%s);Radius (mm);rate (GHz)",spTit[iSp].c_str(),kinEcut[ikinEcut],detH[iDet].c_str(),tgt_gen_config.c_str());
      histXY_title = Form("%s (KE>%.0f MeV) XY dist. on %s (%s);x (mm);y (mm);rate (GHz)",spTit[iSp].c_str(),kinEcut[ikinEcut],detH[iDet].c_str(),tgt_gen_config.c_str());
      histR_E_title = Form("%s (KE>%.0f MeV) Energy Weighted Radial dist. on %s (%s);Radius (mm);rate (GHz)",spTit[iSp].c_str(),kinEcut[ikinEcut],detH[iDet].c_str(),tgt_gen_config.c_str());
      histXY_E_title = Form("%s (KE>%.0f MeV) Energy Weighted XY dist. on %s (%s);x (mm);y (mm);rate (GHz)",spTit[iSp].c_str(),kinEcut[ikinEcut],detH[iDet].c_str(),tgt_gen_config.c_str());
      }
      hR[iSp][iDet][ikinEcut] = new TH1D(Form("%s_R_%s_%.0fMeV",detH[iDet].c_str(),spH[iSp].c_str(),kinEcut[ikinEcut]),histR_title.c_str(),1000,0,2000);
      hR_pzG0[iSp][iDet][ikinEcut] = new TH1D(Form("%s_RpzG0_%s_%.0fMeV",detH[iDet].c_str(),spH[iSp].c_str(),kinEcut[ikinEcut]),histR_title.c_str(),1000,0,2000);
      hR_pzL0[iSp][iDet][ikinEcut] = new TH1D(Form("%s_RpzL0_%s_%.0fMeV",detH[iDet].c_str(),spH[iSp].c_str(),kinEcut[ikinEcut]),histR_title.c_str(),1000,0,2000);
      hXY[iSp][iDet][ikinEcut] = new TH2D(Form("%s_XY_%s_%.0fMeV",detH[iDet].c_str(),spH[iSp].c_str(),kinEcut[ikinEcut]),histXY_title.c_str(),2000,-2000,2000,2000,-2000,2000);
      hXY_pzG0[iSp][iDet][ikinEcut] = new TH2D(Form("%s_XYpzG0_%s_%.0fMeV",detH[iDet].c_str(),spH[iSp].c_str(),kinEcut[ikinEcut]),histXY_title.c_str(),2000,-2000,2000,2000,-2000,2000);
      hXY_pzL0[iSp][iDet][ikinEcut] = new TH2D(Form("%s_XYpzL0_%s_%.0fMeV",detH[iDet].c_str(),spH[iSp].c_str(),kinEcut[ikinEcut]),histXY_title.c_str(),2000,-2000,2000,2000,-2000,2000);

      hR_E[iSp][iDet][ikinEcut] = new TH1D(Form("%s_RE_%s_%.0fMeV",detH[iDet].c_str(),spH[iSp].c_str(),kinEcut[ikinEcut]),histR_E_title.c_str(),1000,0,2000);
      hR_E_pzG0[iSp][iDet][ikinEcut] = new TH1D(Form("%s_REpzG0_%s_%.0fMeV",detH[iDet].c_str(),spH[iSp].c_str(),kinEcut[ikinEcut]),histR_E_title.c_str(),1000,0,2000);
      hR_E_pzL0[iSp][iDet][ikinEcut] = new TH1D(Form("%s_REpzL0_%s_%.0fMeV",detH[iDet].c_str(),spH[iSp].c_str(),kinEcut[ikinEcut]),histR_E_title.c_str(),1000,0,2000);
      hXY_E[iSp][iDet][ikinEcut] = new TH2D(Form("%s_XYE_%s_%.0fMeV",detH[iDet].c_str(),spH[iSp].c_str(),kinEcut[ikinEcut]),histXY_E_title.c_str(),2000,-2000,2000,2000,-2000,2000);
      hXY_E_pzG0[iSp][iDet][ikinEcut] = new TH2D(Form("%s_XYEpzG0_%s_%.0fMeV",detH[iDet].c_str(),spH[iSp].c_str(),kinEcut[ikinEcut]),histXY_E_title.c_str(),2000,-2000,2000,2000,-2000,2000);
      hXY_E_pzL0[iSp][iDet][ikinEcut] = new TH2D(Form("%s_XYEpzL0_%s_%.0fMeV",detH[iDet].c_str(),spH[iSp].c_str(),kinEcut[ikinEcut]),histXY_E_title.c_str(),2000,-2000,2000,2000,-2000,2000);

      hR[iSp][iDet][ikinEcut]->Sumw2();
      hR_pzG0[iSp][iDet][ikinEcut]->Sumw2();
      hR_pzL0[iSp][iDet][ikinEcut]->Sumw2();
      hR_E[iSp][iDet][ikinEcut]->Sumw2();
      hR_E_pzG0[iSp][iDet][ikinEcut]->Sumw2();
      hR_E_pzL0[iSp][iDet][ikinEcut]->Sumw2();
     }
    }
   }

   int nfile=0;
   Long64_t nentry=0;
   long nTotEv=0;
   for(int ifile=1001;ifile<=3000;ifile++){
      string infile = Form("%s/%s/%s_%d.root",rootfile_dir.Data(),tgt_gen_config.c_str(),tgt_gen_config.c_str(),ifile);
      ifstream inf(infile.c_str());
      if(!inf){
        cout<<Form("Skipping %s. File doesn't exist.",infile.c_str())<<endl;
        continue;
      }
      TFile *fin = TFile::Open(infile.c_str(),"READ");
      if(fin->TestBit(TFile::kRecovered)){
        cout<<Form("Skipping %s. Recovered file.",infile.c_str())<<endl;
        continue;
      }
      nfile++;
      
      TTree *T = (TTree*)fin->Get("T");
      if(T==0) return 0;
    
      nentry = T->GetEntries();
      cout<<Form("Found %lld entries in filesplit %d",nentry,nfile)<<endl;
      nTotEv+=nentry;

      std::vector<remollGenericDetectorHit_t> *hit =0;
      std::vector<remollEventParticle_t> *part =0;
      remollBeamTarget_t *bm =0;
      remollEvent_t *ev =0;
      Double_t rate =0;
      Double_t phi = 0;
      Double_t modphi = 0;
      T->SetBranchAddress("hit", & hit);
      T->SetBranchAddress("part", & part);
      T->SetBranchAddress("bm", & bm);
      T->SetBranchAddress("ev", & ev);
      T->SetBranchAddress("rate", & rate);
     
      for(Long64_t ientry=0;ientry<nentry;ientry++){
        T->GetEntry(ientry);
        for(int j=0;j<hit->size();j++){
          if(std::isnan(rate) || std::isinf(rate)) continue;
          if(beamGen) rate = 1.0;
          
          int sp = spM[int(hit->at(j).pid)]-1;
          if(sp==-1) continue;
          int dt = dtM[int(hit->at(j).det)]-1;
          if(dt==-1) continue;
          phi = hit->at(j).ph;
          if(phi<0) phi +=2.0*TMath::Pi();
          modphi = fmod(phi,2.0*TMath::Pi()/7.0);

          double energy = hit->at(j).e;

          for(int ikinEcut=0;ikinEcut<nkinEcut;ikinEcut++){
           if(hit->at(j).e>kinEcut[ikinEcut]){
            hR[sp][dt][ikinEcut]->Fill(hit->at(j).r,rate);
            hXY[sp][dt][ikinEcut]->Fill(hit->at(j).x,hit->at(j).y,rate);
            hR_E[sp][dt][ikinEcut]->Fill(hit->at(j).r,rate*energy/1.0e3);//convert energy to GeV
            hXY_E[sp][dt][ikinEcut]->Fill(hit->at(j).x,hit->at(j).y,rate*energy/1.0e3);
            if(hit->at(j).pz>=0){
             hR_pzG0[sp][dt][ikinEcut]->Fill(hit->at(j).r,rate);
             hXY_pzG0[sp][dt][ikinEcut]->Fill(hit->at(j).x,hit->at(j).y,rate);
             hR_E_pzG0[sp][dt][ikinEcut]->Fill(hit->at(j).r,rate*energy/1.0e3);
             hXY_E_pzG0[sp][dt][ikinEcut]->Fill(hit->at(j).x,hit->at(j).y,rate*energy/1.0e3);
            }else{
             hR_pzL0[sp][dt][ikinEcut]->Fill(hit->at(j).r,rate);
             hXY_pzL0[sp][dt][ikinEcut]->Fill(hit->at(j).x,hit->at(j).y,rate);
             hR_E_pzL0[sp][dt][ikinEcut]->Fill(hit->at(j).r,rate*energy/1.0e3);
             hXY_E_pzL0[sp][dt][ikinEcut]->Fill(hit->at(j).x,hit->at(j).y,rate*energy/1.0e3);
            }
            if(hit->at(j).pid==11 || hit->at(j).pid==-11){
             hR[4][dt][ikinEcut]->Fill(hit->at(j).r,rate);
             hXY[4][dt][ikinEcut]->Fill(hit->at(j).x,hit->at(j).y,rate);
             hR_E[4][dt][ikinEcut]->Fill(hit->at(j).r,rate*energy/1.0e3);
             hXY_E[4][dt][ikinEcut]->Fill(hit->at(j).x,hit->at(j).y,rate*energy/1.0e3);
             if(hit->at(j).pz>=0){
              hR_pzG0[4][dt][ikinEcut]->Fill(hit->at(j).r,rate);
              hXY_pzG0[4][dt][ikinEcut]->Fill(hit->at(j).x,hit->at(j).y,rate);
              hR_E_pzG0[4][dt][ikinEcut]->Fill(hit->at(j).r,rate*energy/1.0e3);
              hXY_E_pzG0[4][dt][ikinEcut]->Fill(hit->at(j).x,hit->at(j).y,rate*energy/1.0e3);
             }else{
              hR_pzL0[4][dt][ikinEcut]->Fill(hit->at(j).r,rate);
              hXY_pzL0[4][dt][ikinEcut]->Fill(hit->at(j).x,hit->at(j).y,rate);
              hR_E_pzL0[4][dt][ikinEcut]->Fill(hit->at(j).r,rate*energy/1.0e3);
              hXY_E_pzL0[4][dt][ikinEcut]->Fill(hit->at(j).x,hit->at(j).y,rate*energy/1.0e3);
             }
            }
            if(hit->at(j).trid==1 && hit->at(j).pid==11){
              hR[5][dt][ikinEcut]->Fill(hit->at(j).r,rate);
              hXY[5][dt][ikinEcut]->Fill(hit->at(j).x,hit->at(j).y,rate);
              hR_E[5][dt][ikinEcut]->Fill(hit->at(j).r,rate*energy/1.0e3);
              hXY_E[5][dt][ikinEcut]->Fill(hit->at(j).x,hit->at(j).y,rate*energy/1.0e3);
              if(hit->at(j).pz>=0){
               hR_pzG0[5][dt][ikinEcut]->Fill(hit->at(j).r,rate);
               hXY_pzG0[5][dt][ikinEcut]->Fill(hit->at(j).x,hit->at(j).y,rate);
               hR_E_pzG0[5][dt][ikinEcut]->Fill(hit->at(j).r,rate*energy/1.0e3);
               hXY_E_pzG0[5][dt][ikinEcut]->Fill(hit->at(j).x,hit->at(j).y,rate*energy/1.0e3);
              }else{
               hR_pzL0[5][dt][ikinEcut]->Fill(hit->at(j).r,rate);
               hXY_pzL0[5][dt][ikinEcut]->Fill(hit->at(j).x,hit->at(j).y,rate);
               hR_E_pzL0[5][dt][ikinEcut]->Fill(hit->at(j).r,rate*energy/1.0e3);
               hXY_E_pzL0[5][dt][ikinEcut]->Fill(hit->at(j).x,hit->at(j).y,rate*energy/1.0e3);
              }
            }
           }
          }
        }
      }
      delete T;
      delete fin;
   }
   
   cout<<Form("Total number of file splits: %d",nfile)<<endl;
   cout<<Form("Total number of entries: %ld",nTotEv)<<endl;
   outfile = TFile::Open(Form("./rootfiles/sam_radiation_study_rad_xy_%s.root",tgt_gen_config.c_str()),"RECREATE");
   outfile->cd();
   for(int iDet=0;iDet<nDet;iDet++){
     outfile->mkdir(Form("%s",detH[iDet].c_str()));
     outfile->cd(Form("%s",detH[iDet].c_str()));
     for(int iSp=0;iSp<nSp;iSp++){
      for(int ikinEcut=0;ikinEcut<nkinEcut;ikinEcut++){
       if(beamGen){
         hR[iSp][iDet][ikinEcut]->Scale(1.0/nTotEv);
         hR_pzG0[iSp][iDet][ikinEcut]->Scale(1.0/nTotEv);
         hR_pzL0[iSp][iDet][ikinEcut]->Scale(1.0/nTotEv);
         hXY[iSp][iDet][ikinEcut]->Scale(1.0/nTotEv);
         hXY_pzG0[iSp][iDet][ikinEcut]->Scale(1.0/nTotEv);
         hXY_pzL0[iSp][iDet][ikinEcut]->Scale(1.0/nTotEv);
         hR_E[iSp][iDet][ikinEcut]->Scale(1.0/nTotEv);
         hR_E_pzG0[iSp][iDet][ikinEcut]->Scale(1.0/nTotEv);
         hR_E_pzL0[iSp][iDet][ikinEcut]->Scale(1.0/nTotEv);
         hXY_E[iSp][iDet][ikinEcut]->Scale(1.0/nTotEv);
         hXY_E_pzG0[iSp][iDet][ikinEcut]->Scale(1.0/nTotEv);
         hXY_E_pzL0[iSp][iDet][ikinEcut]->Scale(1.0/nTotEv);
       }else{
         hR[iSp][iDet][ikinEcut]->Scale(1.0e-9/nfile);
         hR_pzG0[iSp][iDet][ikinEcut]->Scale(1.0e-9/nfile);
         hR_pzL0[iSp][iDet][ikinEcut]->Scale(1.0e-9/nfile);
         hXY[iSp][iDet][ikinEcut]->Scale(1.0e-9/nfile);
         hXY_pzG0[iSp][iDet][ikinEcut]->Scale(1.0e-9/nfile);
         hXY_pzL0[iSp][iDet][ikinEcut]->Scale(1.0e-9/nfile);
         hR_E[iSp][iDet][ikinEcut]->Scale(1.0e-9/nfile);
         hR_E_pzG0[iSp][iDet][ikinEcut]->Scale(1.0e-9/nfile);
         hR_E_pzL0[iSp][iDet][ikinEcut]->Scale(1.0e-9/nfile);
         hXY_E[iSp][iDet][ikinEcut]->Scale(1.0e-9/nfile);
         hXY_E_pzG0[iSp][iDet][ikinEcut]->Scale(1.0e-9/nfile);
         hXY_E_pzL0[iSp][iDet][ikinEcut]->Scale(1.0e-9/nfile);
       }
       hR[iSp][iDet][ikinEcut]->Write();
       hR_pzG0[iSp][iDet][ikinEcut]->Write();
       hR_pzL0[iSp][iDet][ikinEcut]->Write();
       hXY[iSp][iDet][ikinEcut]->Write();
       hXY_pzG0[iSp][iDet][ikinEcut]->Write();
       hXY_pzL0[iSp][iDet][ikinEcut]->Write();
       hR_E[iSp][iDet][ikinEcut]->Write();
       hR_E_pzG0[iSp][iDet][ikinEcut]->Write();
       hR_E_pzL0[iSp][iDet][ikinEcut]->Write();
       hXY_E[iSp][iDet][ikinEcut]->Write();
       hXY_E_pzG0[iSp][iDet][ikinEcut]->Write();
       hXY_E_pzL0[iSp][iDet][ikinEcut]->Write();
     }
    }
   }
}
