//This macro produces radial distribution for various virtual detectors
#include "remolltypes.hh"
#include <sstream>
#include <iostream>
#include <fstream>
void radialCut_manyDet_largeBinSize(){
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

////Change the following lines for which detectors you want to include////
   string detH[] = {"det173","det167","det168","det169","det170","det171","det172","det174","det175","det28","det30","det177","det176","det31"};
//just downstream of target, ring 5, just befor SAM (previous SAM VD location), SAM location in remoll, SAM location in CAD
   const int nDet = sizeof(detH)/sizeof(*detH);
   const int Det[nDet] = {173,167,168,169,170,171,172,174,175,28,30,177,176,31};
   map<int,int> dtM {{173,1},{167,2},{168,3},{169,4},{170,5},{171,6},{172,7},{174,8},{175,9},{28,10},{30,11},{177,12},{176,13},{31,14}};
/////////////////////////////////////////////////////////////////////////

   double x_min = 0;
   double x_max = 1000;
   const int nbin = 200;
   TH1F* h_rate[nSp][nDet];
   TH1F* h_ratePzG0[nSp][nDet];
   TH1F* h_ratePzL0[nSp][nDet];

////Change the following lines as needed////
   const string geometry = "defaultGeo";//defaultGeo or PMTSh
   const string target = "Optics1";
   const string generator = "beam";
   int beamGen(1);
   const string config = "magOFF_V4";
///////////////////////////////////////////

   TFile* outfile = new TFile(Form("./rootfiles/%s_%s_%s_%s_radialCut_largeBinSize.root",geometry.c_str(),
                            target.c_str(),generator.c_str(),config.c_str()),"recreate");
////Change this line for appropriate rootfile directory////
   TString rootfile_dir = "/volatile/halla/parity/adhidevi/remoll_rootfiles/default-geo";
///////////////////////////////////////////////////////////

   for(int iSp=0;iSp<nSp;iSp++){
    for(int iDet=0;iDet<nDet;iDet++){
      if(Det[iDet] == 167 || Det[iDet] == 168) x_min = 20;
      if(Det[iDet] == 169 || Det[iDet] == 170 || Det[iDet] == 171) x_min = 40;
      if(Det[iDet] == 172) x_min = 100;
      if(Det[iDet] == 174 || Det[iDet] == 175) x_min = 150;
      if(Det[iDet] == 176 || Det[iDet] == 177 || Det[iDet]<50) x_min = 200;
      double bin_width = (x_max-x_min)/nbin;

      h_rate[iSp][iDet] = new TH1F(Form("%s_r_%s",detH[iDet].c_str(),spH[iSp].c_str()),
                                  Form("Radial dist. on %s plane (%s gen, %s tgt, %s);Radius (mm);hits/#thrownEvents/%.1fmm",
                                  detH[iDet].c_str(),generator.c_str(),target.c_str(),
                                  config.c_str(),bin_width),nbin,x_min,x_max);
      h_ratePzG0[iSp][iDet] = new TH1F(Form("%s_rPzG0_%s",detH[iDet].c_str(),spH[iSp].c_str()),
                                  Form("Radial dist. on %s plane pz>=0 (%s gen, %s tgt, %s);Radius (mm);hits/#thrownEvents/%.1fmm",
                                  detH[iDet].c_str(),generator.c_str(),target.c_str(),
                                  config.c_str(),bin_width),nbin,x_min,x_max);
      h_ratePzL0[iSp][iDet] = new TH1F(Form("%s_rPzL0_%s",detH[iDet].c_str(),spH[iSp].c_str()),
                                  Form("Radial dist. on %s plane pz<0 (%s gen, %s tgt, %s);Radius (mm);hits/#thrownEvents/%.1fmm",
                                  detH[iDet].c_str(),generator.c_str(),target.c_str(),
                                  config.c_str(),bin_width),nbin,x_min,x_max);

      h_rate[iSp][iDet]->Sumw2();
      h_ratePzG0[iSp][iDet]->Sumw2();
      h_ratePzL0[iSp][iDet]->Sumw2();
    }
   }

   TChain* T = new TChain("T");
   int nfile=0;

   for(int ifile=1001;ifile<=6000;ifile++){
////Change this line for appropriate rootfiles////
       string infile = Form("%s/Optics1_beam_magOFF_V4/Optics1_beam_magOFF_V4_%d.root",
                        rootfile_dir.Data(),ifile);
//////////////////////////////////////////////////
       TFile *fin = new TFile(infile.c_str(),"READ");
       if(!fin->IsOpen() || fin->IsZombie()){
         cout<<"Problem: can't find file: "<<infile<<endl;
         continue;
       }else if(fin->TestBit(TFile::kRecovered)){
         cout<<"Problem: Recovered file: "<<infile<<endl;
         continue;
       }
       nfile++;
       T->Add(Form("%s",infile.c_str()));
   }
   cout<<Form("Found %d number of file splits!",nfile)<<endl;
   Long64_t nentry = T->GetEntries();
   cout<<Form("Total number of entries to be analyzed: %lld",nentry)<<endl;
 
   std::vector<remollGenericDetectorHit_t> *hit =0;
   std::vector<remollEventParticle_t> *part = 0;
   remollBeamTarget_t *bm = 0;
   remollEvent_t *ev = 0;
   Double_t rate = 0;
   T->SetBranchAddress("hit", & hit);
   T->SetBranchAddress("part", & part);
   T->SetBranchAddress("bm", & bm);
   T->SetBranchAddress("ev",& ev);
   T->SetBranchAddress("rate", & rate);
   
   for(Long64_t ientry=0;ientry<nentry;ientry++){
      T->GetEntry(ientry);
      if(ientry%(nentry/10)==0)
        cout<<Form("analyzed %lld events!!",ientry)<<endl;
      for(int j=0;j<hit->size();j++){
         if(std::isnan(rate) || std::isinf(rate)) continue;
         if(beamGen) rate = 1.0;

         int sp = spM[int(hit->at(j).pid)]-1;
         if(sp==-1) continue;
         int dt = dtM[int(hit->at(j).det)]-1;
         if(dt==-1) continue;

         h_rate[sp][dt]->Fill(hit->at(j).r,rate);
         if(hit->at(j).pz>=0)
           h_ratePzG0[sp][dt]->Fill(hit->at(j).r,rate);
         else
           h_ratePzL0[sp][dt]->Fill(hit->at(j).r,rate);
           
         if(hit->at(j).k>1 && (hit->at(j).pid==11 || hit->at(j).pid==-11)){
           h_rate[4][dt]->Fill(hit->at(j).r,rate);
           if(hit->at(j).pz>=0)
             h_ratePzG0[4][dt]->Fill(hit->at(j).r,rate);
           else
             h_ratePzL0[4][dt]->Fill(hit->at(j).r,rate);
         }
         if(hit->at(j).k>1 && hit->at(j).trid==1){
           h_rate[5][dt]->Fill(hit->at(j).r,rate);
           if(hit->at(j).pz>=0)
             h_ratePzG0[5][dt]->Fill(hit->at(j).r,rate);
           else
             h_ratePzL0[5][dt]->Fill(hit->at(j).r,rate);
         }
      }
   }
   
   for(int iDet=0;iDet<nDet;iDet++){
      outfile->mkdir(Form("%s",detH[iDet].c_str()));
      outfile->cd(Form("%s",detH[iDet].c_str()));
      for(int iSp=0;iSp<nSp;iSp++){
       if(beamGen){
         h_rate[iSp][iDet]->Scale(1.0/nentry);
         h_ratePzG0[iSp][iDet]->Scale(1.0/nentry);
         h_ratePzL0[iSp][iDet]->Scale(1.0/nentry);
       }else{
         h_rate[iSp][iDet]->Scale(1.0/nfile);
         h_ratePzG0[iSp][iDet]->Scale(1.0/nfile);
         h_ratePzL0[iSp][iDet]->Scale(1.0/nfile);
       }
         h_rate[iSp][iDet]->Write();
         h_ratePzG0[iSp][iDet]->Write();
         h_ratePzL0[iSp][iDet]->Write();
      }
   }

   int color[nSp] = {2,3,1,4,6,7};
   int colorPzG0[nSp] = {30,30,30,30,30,30};
   int colorPzL0[nSp] = {48,48,48,48,48,48};
   TLatex latex;
   latex.SetNDC(1);
   latex.SetTextSize(0.04);

   for(int iDet=0;iDet<nDet;iDet++){
      for(int iSp=0;iSp<nSp;iSp++){
         h_rate[iSp][iDet]->SetLineColor(color[iSp]);
         h_ratePzG0[iSp][iDet]->SetLineColor(colorPzG0[iSp]);
         h_ratePzL0[iSp][iDet]->SetLineColor(colorPzL0[iSp]);
      }
   }

   TCanvas* c_r_lin[nDet];
   for(int iDet=0;iDet<nDet;iDet++){
      c_r_lin[iDet] = new TCanvas(Form("c_r_lin_%s",detH[iDet].c_str()));
      for(int iSp=0;iSp<nSp;iSp++){
         if(iSp==0)
         h_rate[iSp][iDet]->Draw("hist");
         else
         h_rate[iSp][iDet]->Draw("hist same");

         latex.SetTextColor(color[iSp]);
         if(iSp<3)
           latex.DrawLatex(0.33,0.85-0.05*iSp,Form("%s",spTit[iSp].c_str()));
         else
           latex.DrawLatex(0.56,0.85-0.05*(iSp-3),Form("%s",spTit[iSp].c_str()));
      }
      c_r_lin[iDet]->SaveAs(Form("./temp/%s_%s_%s_%s_radialCut_largeBinSize_lin.pdf",geometry.c_str(),
                           target.c_str(),generator.c_str(),config.c_str()));
    }

   TCanvas* c_r_log[nDet];
   for(int iDet=0;iDet<nDet;iDet++){
      c_r_log[iDet] = new TCanvas(Form("c_r_log_%s",detH[iDet].c_str()));
      gPad->SetLogy(1);
      for(int iSp=0;iSp<nSp;iSp++){
         if(iSp==0)
         h_rate[iSp][iDet]->Draw("hist");
         else
         h_rate[iSp][iDet]->Draw("hist same");

         latex.SetTextColor(color[iSp]);
         if(iSp<3)
           latex.DrawLatex(0.33,0.85-0.05*iSp,Form("%s",spTit[iSp].c_str()));
         else
           latex.DrawLatex(0.56,0.85-0.05*(iSp-3),Form("%s",spTit[iSp].c_str()));
      }
      c_r_log[iDet]->SaveAs(Form("./temp/%s_%s_%s_%s_radialCut_largeBinSize_log.pdf",geometry.c_str(),
                           target.c_str(),generator.c_str(),config.c_str()));
   }

   TCanvas* c_rFB_lin[nSp][nDet];
   for(int iDet=0;iDet<nDet;iDet++){
      for(int iSp=0;iSp<nSp;iSp++){
         c_rFB_lin[iSp][iDet] = new TCanvas(Form("%s_rFB_lin_%s",spH[iSp].c_str(),detH[iDet].c_str()));
         h_rate[iSp][iDet]->Draw("hist");
         h_ratePzG0[iSp][iDet]->Draw("hist same");
         h_ratePzL0[iSp][iDet]->Draw("hist same");
         h_rate[3][iDet]->GetYaxis()->SetRangeUser(1.e-7,1.e-5);

         latex.SetTextColor(color[iSp]);
         latex.DrawLatex(0.25,0.85,Form("%s all Pz",spTit[iSp].c_str()));
         latex.SetTextColor(colorPzG0[iSp]);
         latex.DrawLatex(0.25,0.80,Form("%s Pz>=0",spTit[iSp].c_str()));
         latex.SetTextColor(colorPzL0[iSp]);
         latex.DrawLatex(0.25,0.75,Form("%s Pz<0",spTit[iSp].c_str()));
         c_rFB_lin[iSp][iDet]->SaveAs(Form("./temp/%s_%s_%s_%s_radialCut_largeBinSize_%s_%s_lin.pdf",geometry.c_str(),
                   target.c_str(),generator.c_str(),config.c_str(),spH[iSp].c_str(),detH[iDet].c_str()));
      }
    }

   TCanvas* c_rFB_log[nSp][nDet];
   for(int iDet=0;iDet<nDet;iDet++){
      for(int iSp=0;iSp<nSp;iSp++){
         c_rFB_log[iSp][iDet] = new TCanvas(Form("%s_rFB_log_%s",spH[iSp].c_str(),detH[iDet].c_str()));
         gPad->SetLogy(1);
         h_rate[iSp][iDet]->Draw("hist");
         h_ratePzG0[iSp][iDet]->Draw("hist same");
         h_ratePzL0[iSp][iDet]->Draw("hist same");

         latex.SetTextColor(color[iSp]);
         latex.DrawLatex(0.25,0.85,Form("%s all Pz",spTit[iSp].c_str()));
         latex.SetTextColor(colorPzG0[iSp]);
         latex.DrawLatex(0.25,0.80,Form("%s Pz>=0",spTit[iSp].c_str()));
         latex.SetTextColor(colorPzL0[iSp]);
         latex.DrawLatex(0.25,0.75,Form("%s Pz<0",spTit[iSp].c_str()));
         c_rFB_log[iSp][iDet]->SaveAs(Form("./temp/%s_%s_%s_%s_radialCut_largeBinSize_%s_%s_log.pdf",geometry.c_str(),
                    target.c_str(),generator.c_str(),config.c_str(),spH[iSp].c_str(),detH[iDet].c_str()));
      }
   }
//Now combine all pdf files saved in ./temp/ directory and save a single pdf file in ./plots/ directory
   gSystem->Exec(Form("pdfunite ./temp/%s_%s_%s_%s_radialCut_largeBinSize_* ./plots/%s_%s_%s_%s_radialCut_largeBinSize.pdf",
                geometry.c_str(),target.c_str(),generator.c_str(),config.c_str(),geometry.c_str(),
                target.c_str(),generator.c_str(),config.c_str()));
   gSystem->Exec(Form("rm -rf ./temp/%s_%s_%s_%s_radialCut_largeBinSize_*",geometry.c_str(),target.c_str(),
                generator.c_str(),config.c_str()));

}
