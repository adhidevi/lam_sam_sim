//This macro produces radial distribution for a given detector and given species and integrates within the given radial extent.
//
#include "remolltypes.hh"
#include <sstream>
#include <iostream>
#include <fstream>
void Optics1_vs_LH2_compare(){
   gROOT->Reset();
   gStyle->SetOptStat(0);
   gStyle->SetTitleYOffset(1.3);
   gStyle->SetPadGridX(1);
   gStyle->SetPadGridY(1);
   TGaxis::SetMaxDigits(3);
   const string spTit[] = {"e-/#pi-","e+/#pi+","#gamma","neutron","e-/e+ E>1","primary e E>1"};
   const int nSp = sizeof(spTit)/sizeof(*spTit);
   const string spH[nSp] = {"epiM","epiP","g","n","ee1","pri1"};
   map<int,int> spM {{11,1},{-211,1},{-11,2},{211,2},{22,3},{2112,4}};

   const int Det = 28;
   map<int,int> dtM {{28,1}};
   const int x_min = 0;
   const int x_max = 1500;
   const int nbins = 500;
   double bin_width = (x_max-x_min)/nbins;
   const int r_min = 920;
   const int r_max = 1060;
   
   int beamGen(1);//Change this if using physics generators
   const string geometry = "defaultGeo";//defaultGeo or PMTSh
   const string sim[] = {"Optics1_beam","Optics1_beam_magON_V4"};
   const int nsim = sizeof(sim)/sizeof(*sim);
   TH1D* h_rate[nsim][nSp];
   TH1D* h_rateQ[nsim][nSp];

   TFile *outfile = new TFile(Form("./rootfiles/%s_%s_%s_det%d_radial.root",geometry.c_str(),sim[0].c_str(),sim[1].c_str(),Det),"recreate");
   
   string rootfile_dir = "/volatile/halla/parity/adhidevi/remoll_rootfiles/default-geo";
   
   Long64_t nentry[nsim]={0,0};
   int nfile[nsim] = {0,0};

   for(int isim=0;isim<nsim;isim++){
    for(int iSp=0;iSp<nSp;iSp++){
      h_rate[isim][iSp] = new TH1D(Form("det%d_%s_%s",Det,spH[iSp].c_str(),sim[isim].c_str()),Form("Radial Distribution of %s at det%d;Radius (mm); hits/%sthrownEvents)",spTit[iSp].c_str(),Det,"#"),nbins,x_min,x_max);
      h_rateQ[isim][iSp] = new TH1D(Form("det%dQ_%s_%s",Det,spH[iSp].c_str(),sim[isim].c_str()),Form("Radial Distribution of %s at det%d;Radius (mm); hits/%sthrownEvents)",spTit[iSp].c_str(),Det,"#"),nbins,x_min,x_max);
      h_rate[isim][iSp]->Sumw2();
      h_rateQ[isim][iSp]->Sumw2();
    }

      TChain* T = new TChain("T");
      for(int ifile=1001;ifile<=2000;ifile++){
         string infile = Form("%s/%s/%s_%d.root",rootfile_dir.c_str(),sim[isim].c_str(),sim[isim].c_str(),ifile);
         TFile *fin = new TFile(infile.c_str(),"READ");
         if(!fin->IsOpen() || fin->IsZombie()){
           cout<<"Problem: can't find file: "<<infile<<endl;
           continue;
         }else if(fin->TestBit(TFile::kRecovered)){
           cout<<"Problem: Recovered file: "<<infile<<endl;
           continue;
         }
         nfile[isim]++;
         T->Add(Form("%s",infile.c_str()));
      }
      cout<<Form("Found %d file splits in %s!!!",nfile[isim],sim[isim].c_str())<<endl;

      nentry[isim] = T->GetEntries();
      cout<<Form("Total number of entries to be analyzed: %lld in %s",nentry[isim],sim[isim].c_str())<<endl;

      std::vector<remollGenericDetectorHit_t> *hit = 0;
      std::vector<remollEventParticle_t> *part = 0;
      remollBeamTarget_t *bm = 0;
      remollEvent_t *ev = 0;
      Double_t rate = 0;

      T->SetBranchAddress("hit", &hit);
      T->SetBranchAddress("part", &part);
      T->SetBranchAddress("bm", &bm);
      T->SetBranchAddress("ev", &ev);
      T->SetBranchAddress("rate", &rate);

      for(Long64_t ientry=0;ientry<nentry[isim];ientry++){
         T->GetEntry(ientry);
         if(ientry%(nentry[isim]/10)==0)
            cout<<Form("Analyzed %lld events!!!",ientry)<<endl;
         for(size_t j=0;j<hit->size();j++){
            if(std::isnan(rate) || std::isinf(rate)) continue;
            if(beamGen) rate = 1.0;

            int sp = spM[int(hit->at(j).pid)]-1;
            if(sp==-1) continue;
            int dt = dtM[int(hit->at(j).det)]-1;
            if(dt==-1) continue;

            h_rate[isim][sp]->Fill(hit->at(j).r,rate);

            if(hit->at(j).r>=r_min && hit->at(j).r<=r_max)
              h_rateQ[isim][sp]->Fill(hit->at(j).r,rate);
            
            if(hit->at(j).k>1 && (hit->at(j).pid==11 || hit->at(j).pid==-11)){
              h_rate[isim][4]->Fill(hit->at(j).r,rate);
              if(hit->at(j).r>=r_min && hit->at(j).r<=r_max)
                h_rateQ[isim][4]->Fill(hit->at(j).r,rate);
            }
            if(hit->at(j).k>1 && hit->at(j).trid==1){
              h_rate[isim][5]->Fill(hit->at(j).r,rate);
              if(hit->at(j).r>=r_min && hit->at(j).r<=r_max)
                h_rateQ[isim][5]->Fill(hit->at(j).r,rate);
            }
         }
      }
   }

   for(int isim=0;isim<nsim;isim++){
      outfile->mkdir(Form("%s",sim[isim].c_str()));
      outfile->cd(Form("%s",sim[isim].c_str()));
      for(int iSp=0;iSp<nSp;iSp++){
       if(beamGen){
         h_rate[isim][iSp]->Scale(1./nentry[isim]);
         h_rateQ[isim][iSp]->Scale(1./nentry[isim]);
       }else{
         h_rate[isim][iSp]->Scale(1./nfile[isim]);
         h_rateQ[isim][iSp]->Scale(1./nfile[isim]);
       }   
       h_rate[isim][iSp]->Write();
       h_rateQ[isim][iSp]->Write();
      }
   }

   int color[nsim][nSp] = {{2,3,1,4,6,7},{8,9,40,41,42,43}};
   
   TLatex latex;
   latex.SetNDC();
   latex.SetTextSize(0.04);

   TH1D* h_rateCpy[nsim][nSp];
   int r_elastic_min = 690;
   int r_elastic_max = 780;

   for(int isim=0;isim<nsim;isim++){
     for(int iSp=0;iSp<nSp;iSp++){
        h_rate[isim][iSp]->SetLineColor(color[isim][iSp]);
        h_rateQ[isim][iSp]->SetLineColor(color[isim][iSp]);
       h_rateCpy[isim][iSp] = (TH1D*)h_rate[isim][iSp]->Clone(Form("det%dCpy_%s_%s",Det,spH[iSp].c_str(),sim[isim].c_str()));
       h_rateCpy[isim][iSp]->GetXaxis()->SetRangeUser(r_elastic_min,r_elastic_max);
        h_rateCpy[isim][iSp]->SetLineColor(color[isim][iSp]);
     }
   }

   TCanvas* c_r_lin[nSp];
   for(int iSp=0;iSp<nSp;iSp++){
     c_r_lin[iSp] = new TCanvas(Form("c_r_lin_%s",spH[iSp].c_str()));
     h_rate[0][iSp]->Draw("hist");
     h_rate[0][iSp]->GetYaxis()->SetRangeUser(1.5e-8,1.0);
     h_rate[1][iSp]->Draw("hist sames");
     h_rateQ[0][iSp]->Draw("hist sames");
     h_rateQ[1][iSp]->Draw("hist sames");
     h_rateCpy[0][iSp]->Draw("hist sames");
     h_rateCpy[1][iSp]->Draw("hist sames");

     latex.SetTextColor(color[0][iSp]);
     latex.DrawLatex(0.15,0.85,Form("%s",sim[0].c_str()));
     latex.SetTextColor(color[1][iSp]);
     latex.DrawLatex(0.15,0.80,Form("%s",sim[1].c_str()));
     latex.SetTextColor(color[0][iSp]);
     latex.DrawLatex(0.51,0.85,Form("intg[%d<=r<=%d]: %.3e",r_min,r_max,h_rateQ[0][iSp]->Integral()));
     latex.SetTextColor(color[1][iSp]);
     latex.DrawLatex(0.51,0.80,Form("intg[%d<=r<=%d]: %.3e",r_min,r_max,h_rateQ[1][iSp]->Integral()));

     latex.SetTextColor(color[0][iSp]);
     latex.DrawLatex(0.51,0.75,Form("intg[%d<=r<=%d]: %.3e",r_elastic_min,r_elastic_max,h_rateCpy[0][iSp]->Integral()));
     latex.SetTextColor(color[1][iSp]);
     latex.DrawLatex(0.51,0.70,Form("intg[%d<=r<=%d]: %.3e",r_elastic_min,r_elastic_max,h_rateCpy[1][iSp]->Integral()));

     c_r_lin[iSp]->SaveAs(Form("./temp/det%d_lin_%s.pdf",Det,spH[iSp].c_str()));
   }
   TCanvas* c_r_log[nSp];
   for(int iSp=0;iSp<nSp;iSp++){
     c_r_log[iSp] = new TCanvas(Form("c_r_log_%s",spH[iSp].c_str()));
     gPad->SetLogy();
     h_rate[0][iSp]->Draw("hist");
     h_rate[0][iSp]->GetYaxis()->SetRangeUser(1.5e-8,1.0);
     h_rate[1][iSp]->Draw("hist sames");
     h_rateQ[0][iSp]->Draw("hist sames");
     h_rateQ[1][iSp]->Draw("hist sames");
     h_rateCpy[0][iSp]->Draw("hist sames");
     h_rateCpy[1][iSp]->Draw("hist sames");

     latex.SetTextColor(color[0][iSp]);
     latex.DrawLatex(0.15,0.85,Form("%s",sim[0].c_str()));
     latex.SetTextColor(color[1][iSp]);
     latex.DrawLatex(0.15,0.80,Form("%s",sim[1].c_str()));
     latex.SetTextColor(color[0][iSp]);
     latex.DrawLatex(0.51,0.85,Form("intg[%d<=r<=%d]: %.3e",r_min,r_max,h_rateQ[0][iSp]->Integral()));
     latex.SetTextColor(color[1][iSp]);
     latex.DrawLatex(0.51,0.80,Form("intg[%d<=r<=%d]: %.3e",r_min,r_max,h_rateQ[1][iSp]->Integral()));

     latex.SetTextColor(color[0][iSp]);
     latex.DrawLatex(0.51,0.75,Form("intg[%d<=r<=%d]: %.3e",r_elastic_min,r_elastic_max,h_rateCpy[0][iSp]->Integral()));
     latex.SetTextColor(color[1][iSp]);
     latex.DrawLatex(0.51,0.70,Form("intg[%d<=r<=%d]: %.3e",r_elastic_min,r_elastic_max,h_rateCpy[1][iSp]->Integral()));

     c_r_log[iSp]->SaveAs(Form("./temp/det%d_log_%s.pdf",Det,spH[iSp].c_str()));
   }
   
   gSystem->Exec(Form("pdfunite ./temp/det%d_*.pdf ./plots/%s_%s_%s_det%d_radial.pdf",Det,geometry.c_str(),sim[0].c_str(),sim[1].c_str(),Det));
   gSystem->Exec(Form("rm -rf ./temp/det%d_*.pdf",Det));
}
