//This macro produces radial distribution for various virtual detectors
#include "remolltypes.hh"
#include <sstream>
#include <iostream>
#include <fstream>
void quick_radialCut_manyDet(){
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

///Change the following lines for which detectors you want to include////
  string detH[] = {"det173","det167","det168","det169","det170","det171","det172","det174","det175","det28","det30","det177","det176","det31"};
  const int nDet = sizeof(detH)/sizeof(*detH);
  const int Det[nDet] = {173,167,168,169,170,171,172,174,175,28,30,177,176,31};
  map<int,int> dtM {{173,1},{167,2},{168,3},{169,4},{170,5},{171,6},{172,7},{174,8},{175,9},{28,10},{30,11},{177,12},{176,13},{31,14}};
////////////////////////////////////////////////////////////////////////

  double x_min = 0;
  double x_max = 1500;
  const int nbin = 200;
  TH1F* h_rate[nSp][nDet];
  TH1F* h_ratePzG0[nSp][nDet];
  TH1F* h_ratePzL0[nSp][nDet];

///Change the following lines as needed////
  const string geometry = "defaultGeo";//defaultGeo or PMTSh
  const string tgt_gen_config = "Optics1_beam_magON_V4";
  const string plotType = "quick_radialCut_200bins";
  int beamGen(1);
//////////////////////////////////////////

  TFile* outfile = new TFile(Form("./rootfiles/%s_%s_%s.root",geometry.c_str(),tgt_gen_config.c_str(),plotType.c_str()),"recreate");
///Change this line for appropriate rootfile directory////
  TString rootfile_dir = "/volatile/halla/parity/adhidevi/remoll_rootfiles/default-geo";
//////////////////////////////////////////////////////////

  for(int iSp=0;iSp<nSp;iSp++){
   for(int iDet=0;iDet<nDet;iDet++){
     if(Det[iDet] == 167 || Det[iDet] == 168) x_min = 20;
     if(Det[iDet] == 169 || Det[iDet] == 170 || Det[iDet] == 171) x_min = 40;
     if(Det[iDet] == 172) x_min = 100;
     if(Det[iDet] == 173) x_min = 0;
     if(Det[iDet] == 174 || Det[iDet] == 175) x_min = 150;
     if(Det[iDet] == 176 || Det[iDet] == 177 || Det[iDet]<50) x_min = 200;
     double bin_width = (x_max-x_min)/nbin;

     h_rate[iSp][iDet] = new TH1F(Form("%s_r_%s",detH[iDet].c_str(),spH[iSp].c_str()),Form("Radial dist. on %s plane (%s);Radius (mm);hits/#thrownEvents/%.1fmm",detH[iDet].c_str(),tgt_gen_config.c_str(),bin_width),nbin,x_min,x_max);
     h_ratePzG0[iSp][iDet] = new TH1F(Form("%s_rPzG0_%s",detH[iDet].c_str(),spH[iSp].c_str()),Form("Radial dist. on %s plane pz>=0 (%s);Radius (mm);hits/#thrownEvents/%.1fmm",detH[iDet].c_str(),tgt_gen_config.c_str(),bin_width),nbin,x_min,x_max);
     h_ratePzL0[iSp][iDet] = new TH1F(Form("%s_rPzL0_%s",detH[iDet].c_str(),spH[iSp].c_str()),Form("Radial dist. on %s plane pz<0 (%s);Radius (mm);hits/#thrownEvents/%.1fmm",detH[iDet].c_str(),tgt_gen_config.c_str(),bin_width),nbin,x_min,x_max);

     h_rate[iSp][iDet]->Sumw2();
     h_ratePzG0[iSp][iDet]->Sumw2();
     h_ratePzL0[iSp][iDet]->Sumw2();
   }
  }

  int nfile=0;
  Long64_t nentry=0;
  long nTotEv=0;
  for(int ifile=1001;ifile<=2000;ifile++){
///Change this line for appropriate rootfiles////
    string infile = Form("%s/%s/%s_%d.root",rootfile_dir.Data(),tgt_gen_config.c_str(),tgt_gen_config.c_str(),ifile);
//////////////////////////////////////////////
    ifstream inf(infile.c_str());
    if(!inf){
      cout<<Form("Skipping %s. File doesn't exist.",infile.c_str())<<endl;
      continue;
    }
    TFile *fin = TFile::Open(infile.c_str(),"READ");
    if(fin->TestBit(TFile::kRecovered)){
      cout<<Form("Skipping %s. Recovered file.",infile.c_str())<<endl;
      fin->Close(); delete fin; return 0;
    }
    nfile++;
   
    TTree *T = (TTree*)fin->Get("T");
    if(T==0) return 0;

    nentry = T->GetEntries();
    cout<<Form("Found %lld entries in filesplit %d",nentry,nfile)<<endl;
    nTotEv+=nentry;

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
    delete fin;
  }
  cout<<Form("Total number of file splits: %d",nfile)<<endl;
  cout<<Form("Total number of entries: %ld",nTotEv)<<endl;
 
   
   for(int iDet=0;iDet<nDet;iDet++){
      outfile->mkdir(Form("%s",detH[iDet].c_str()));
      outfile->cd(Form("%s",detH[iDet].c_str()));
      for(int iSp=0;iSp<nSp;iSp++){
       if(beamGen){
         h_rate[iSp][iDet]->Scale(1.0/nTotEv);
         h_ratePzG0[iSp][iDet]->Scale(1.0/nTotEv);
         h_ratePzL0[iSp][iDet]->Scale(1.0/nTotEv);
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
      if(Det[iDet]==173)
      h_rate[5][iDet]->Draw("hist");
      else
      h_rate[2][iDet]->Draw("hist");
      for(int iSp=0;iSp<nSp;iSp++){
         h_rate[iSp][iDet]->Draw("hist same");
         latex.SetTextColor(color[iSp]);
         if(Det[iDet]>=167 && Det[iDet]<=173){
         if(iSp<3)
           latex.DrawLatex(0.60,0.85-0.05*iSp,Form("%s",spTit[iSp].c_str()));
         else
           latex.DrawLatex(0.75,0.85-0.05*(iSp-3),Form("%s",spTit[iSp].c_str()));
         }else{
         if(iSp<3)
           latex.DrawLatex(0.33,0.85-0.05*iSp,Form("%s",spTit[iSp].c_str()));
         else
           latex.DrawLatex(0.46,0.85-0.05*(iSp-3),Form("%s",spTit[iSp].c_str()));
         }
      }
      c_r_lin[iDet]->SaveAs(Form("./temp/%s_%s_%s_%s_lin.pdf",geometry.c_str(),tgt_gen_config.c_str(),plotType.c_str(),detH[iDet].c_str()));
    }

   TCanvas* c_r_log[nDet];
   for(int iDet=0;iDet<nDet;iDet++){
      c_r_log[iDet] = new TCanvas(Form("c_r_log_%s",detH[iDet].c_str()));
      gPad->SetLogy(1);
      if(Det[iDet]==173)
      h_rate[5][iDet]->Draw("hist");
      else
      h_rate[2][iDet]->Draw("hist");
      for(int iSp=0;iSp<nSp;iSp++){
         h_rate[iSp][iDet]->Draw("hist same");
         latex.SetTextColor(color[iSp]);
         if(Det[iDet]>=167 && Det[iDet]<=173){
         if(iSp<3)
           latex.DrawLatex(0.60,0.85-0.05*iSp,Form("%s",spTit[iSp].c_str()));
         else
           latex.DrawLatex(0.73,0.85-0.05*(iSp-3),Form("%s",spTit[iSp].c_str()));
         }else{
         if(iSp<3)
           latex.DrawLatex(0.33,0.85-0.05*iSp,Form("%s",spTit[iSp].c_str()));
         else
           latex.DrawLatex(0.46,0.85-0.05*(iSp-3),Form("%s",spTit[iSp].c_str()));
         }
      }
      c_r_log[iDet]->SaveAs(Form("./temp/%s_%s_%s_%s_log.pdf",geometry.c_str(),tgt_gen_config.c_str(),plotType.c_str(),detH[iDet].c_str()));
   }

   TCanvas* c_rFB_lin[nSp][nDet];
   for(int iDet=0;iDet<nDet;iDet++){
      for(int iSp=0;iSp<nSp;iSp++){
         c_rFB_lin[iSp][iDet] = new TCanvas(Form("%s_rFB_lin_%s",spH[iSp].c_str(),detH[iDet].c_str()));
         h_rate[iSp][iDet]->Draw("hist");
         h_ratePzG0[iSp][iDet]->Draw("hist same");
         h_ratePzL0[iSp][iDet]->Draw("hist same");

         latex.SetTextColor(color[iSp]);
         latex.DrawLatex(0.25,0.85,Form("%s all Pz",spTit[iSp].c_str()));
         latex.SetTextColor(colorPzG0[iSp]);
         latex.DrawLatex(0.25,0.80,Form("%s Pz>=0",spTit[iSp].c_str()));
         latex.SetTextColor(colorPzL0[iSp]);
         latex.DrawLatex(0.25,0.75,Form("%s Pz<0",spTit[iSp].c_str()));
         c_rFB_lin[iSp][iDet]->SaveAs(Form("./temp/%s_%s_%s_%s_%s_lin.pdf",geometry.c_str(),tgt_gen_config.c_str(),plotType.c_str(),spH[iSp].c_str(),detH[iDet].c_str()));
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
         c_rFB_log[iSp][iDet]->SaveAs(Form("./temp/%s_%s_%s_%s_%s_log.pdf",geometry.c_str(),tgt_gen_config.c_str(),plotType.c_str(),spH[iSp].c_str(),detH[iDet].c_str()));
      }
   }
//Now combine all pdf files saved in ./temp/ directory and save a single pdf file in ./plots/ directory
   gSystem->Exec(Form("pdfunite ./temp/%s_%s_%s_* ./plots/%s_%s_%s.pdf",geometry.c_str(),tgt_gen_config.c_str(),plotType.c_str(),geometry.c_str(),tgt_gen_config.c_str(),plotType.c_str()));
   gSystem->Exec(Form("rm -rf ./temp/%s_%s_%s_*",geometry.c_str(),tgt_gen_config.c_str(),plotType.c_str()));

}
