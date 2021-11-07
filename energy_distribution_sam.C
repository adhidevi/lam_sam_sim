//This macro produces the energy distribution of SAM VD at given two radial ranges
#include "remolltypes.hh"
#include <sstream>
#include <iostream>
#include <fstream>
    string rootfile_dir = "/volatile/halla/parity/adhidevi/remoll_rootfiles/PMTShielding";
    const string spTit[] = {"e-/#pi-","e+/#pi+","#gamma","neutron","e-/e+ E>1","primary e E>1"};
    const int nSp = sizeof(spTit)/sizeof(*spTit);
    const string spH[nSp] ={"epiM","epiP","g","n","ee1","pri1"};
    map<int,int> spM {{11,1},{-211,1},{-11,2},{211,2},{22,3},{2112,4}};
    int color[nSp] = {2,3,1,4,6,7};
    const string detH[] = {"det176"};
    map<int,int> dtM {{176,1}};
    const int nDet = sizeof(detH)/sizeof(*detH);
    int Det[nDet] = {176};
    double rmin1[nDet] = {875};
    double rmax1[nDet] = {955};
    double rmin2[nDet] = {1095};
    double rmax2[nDet] = {1295};
    TH1D* eRate1[nSp][nDet];
    TH1D* eRate2[nSp][nDet];
    void niceLogBins(TH1*);
    TFile* outfile; 
void energy_distribution_sam(){
   gStyle->SetOptStat(0);
   for(int iSp=0;iSp<nSp;iSp++){
   for(int iDet=0;iDet<nDet;iDet++){
      eRate1[iSp][iDet] = new TH1D(Form("%s_E_%s_r1",detH[iDet].c_str(),spH[iSp].c_str()),Form("Kinetic Energy for %s, %.0f<=r<=%.0f (beam generator);E (MeV);hits/%sthrownEvents",detH[iDet].c_str(),rmin1[iDet],rmax1[iDet],"#"),121,-8,5);
      eRate2[iSp][iDet] = new TH1D(Form("%s_E_%s_r2",detH[iDet].c_str(),spH[iSp].c_str()),Form("Kinetic Energy for %s, %.0f<=r<=%.0f (beam generator);E (MeV);hits/%sthrownEvents",detH[iDet].c_str(),rmin2[iDet],rmax2[iDet],"#"),121,-8,5);
      eRate1[iSp][iDet]->SetLineColor(color[iSp]);
      eRate2[iSp][iDet]->SetLineColor(color[iSp]);
      niceLogBins(eRate1[iSp][iDet]);
      niceLogBins(eRate2[iSp][iDet]);
   }
   }
    TChain* T = new TChain("T");
    int nfile = 0;
    const string sim = "PMTSh_beam_V3";
    for(int ifile=1001;ifile<=2000;ifile++){
       string infile = Form("%s/%s/%s_%d.root",rootfile_dir.c_str(),sim.c_str(),sim.c_str(),ifile);
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

    cout<<Form("Found %d files!!!",nfile)<<endl;
    Long64_t nentry = T->GetEntries();
    cout<<Form("Total number of entries to be analyzed: %lld",nentry)<<endl;

    outfile = new TFile(Form("./rootfiles/det176_energy_dist_%s_det176.root",sim.c_str()),"recreate");

    Double_t rate=1;
    std::vector<remollGenericDetectorHit_t> *hit=0;
    T->SetBranchAddress("hit", &hit);
   
    for(Long64_t ientry=0;ientry<nentry;ientry++){
       T->GetEntry(ientry);
       if(ientry%(nentry/10)==0)
         cout<<Form("analyzed %lld enents!!",ientry)<<endl;
       for(int j=0;j<hit->size();j++){
          if(std::isnan(rate) || std::isinf(rate)) continue;
          if(rate==0) rate=1;
          int sp = spM[int(hit->at(j).pid)]-1;
          if(sp==-1) continue;
          int dt = dtM[int(hit->at(j).det)]-1;
          if(dt==-1) continue;

          if(hit->at(j).r>=rmin1[dt] && hit->at(j).r<=rmax1[dt]){
            eRate1[sp][dt]->Fill(hit->at(j).k,rate);
            if(hit->at(j).k>1 && (hit->at(j).pid==11 || hit->at(j).pid==-11)){
              eRate1[4][dt]->Fill(hit->at(j).k,rate);
            }
            if(hit->at(j).k>1 && hit->at(j).trid==1){
              eRate1[5][dt]->Fill(hit->at(j).k,rate);
            }
          }

          if(hit->at(j).r>=rmin2[dt] && hit->at(j).r<=rmax2[dt]){
            eRate2[sp][dt]->Fill(hit->at(j).k,rate);
            if(hit->at(j).k>1 && (hit->at(j).pid==11 || hit->at(j).pid==-11)){
              eRate2[4][dt]->Fill(hit->at(j).k,rate);
            }
            if(hit->at(j).k>1 && hit->at(j).trid==1){
              eRate2[5][dt]->Fill(hit->at(j).k,rate);
            }
          }
       }
    }

    TCanvas* c_E1[nDet];
    TCanvas* c_E2[nDet];
    TLatex latex;
    latex.SetNDC(1);
    latex.SetTextSize(0.04);
    for(int iDet=0;iDet<nDet;iDet++){
       c_E1[iDet] = new TCanvas(Form("c_E1_%d",iDet));
       for(int iSp=0;iSp<nSp;iSp++){
          gPad->SetLogx();
          gPad->SetLogy();
          eRate1[iSp][iDet]->Scale(1./nentry);
          eRate1[iSp][iDet]->GetYaxis()->SetRangeUser(1.0e-9,1.0e-2);
          eRate1[iSp][iDet]->SetMarkerStyle(20);
          eRate1[iSp][iDet]->SetMarkerSize(0.5);
          eRate1[iSp][iDet]->Draw("hist same");
          eRate1[iSp][iDet]->Write();
          latex.SetTextColor(color[iSp]);
          if(iSp<3)
            latex.DrawLatex(0.33,0.85-0.05*iSp,Form("%s: %.2e",spTit[iSp].c_str(),eRate1[iSp][iDet]->Integral()));
          else
            latex.DrawLatex(0.56,0.85-0.05*(iSp-3),Form("%s: %.2e",spTit[iSp].c_str(),eRate1[iSp][iDet]->Integral()));
       }
    latex.SetTextColor(kRed+3);
    latex.DrawLatex(0.12,0.85,"Radial Cut::");
    latex.DrawLatex(0.12,0.80,Form("[%.0f,%.0f]mm",rmin1[iDet],rmax1[iDet]));
    c_E1[iDet]->SaveAs(Form("./temp/kinetic_energy1_det%d.pdf",Det[iDet]));
    }
    for(int iDet=0;iDet<nDet;iDet++){
       c_E2[iDet] = new TCanvas(Form("c_E2_%d",iDet));
       for(int iSp=0;iSp<nSp;iSp++){
          gPad->SetLogx();
          gPad->SetLogy();
          eRate2[iSp][iDet]->Scale(1./nentry);
          eRate2[iSp][iDet]->GetYaxis()->SetRangeUser(1.0e-9,1.0e-2);
          eRate2[iSp][iDet]->SetMarkerStyle(20);
          eRate2[iSp][iDet]->SetMarkerSize(0.5);
          eRate2[iSp][iDet]->Draw("hist same");
          eRate2[iSp][iDet]->Write();
          latex.SetTextColor(color[iSp]);
          if(iSp<3)
            latex.DrawLatex(0.33,0.85-0.05*iSp,Form("%s: %.2e",spTit[iSp].c_str(),eRate2[iSp][iDet]->Integral()));
          else
            latex.DrawLatex(0.56,0.85-0.05*(iSp-3),Form("%s: %.2e",spTit[iSp].c_str(),eRate2[iSp][iDet]->Integral()));
       }
    latex.SetTextColor(kRed+3);
    latex.DrawLatex(0.12,0.85,"Radial Cut::");
    latex.DrawLatex(0.12,0.80,Form("[%.0f,%.0f]mm",rmin2[iDet],rmax2[iDet]));
    c_E2[iDet]->SaveAs(Form("./temp/kinetic_energy2_det%d.pdf",Det[iDet]));
    }
    gSystem->Exec(Form("pdfunite ./temp/kinetic_energy*.pdf ./plots/PMTSh_beam_kinE_det%d_det%d.pdf",Det[0],Det[0]));
    gSystem->Exec(Form("rm -rf ./temp/kinetic_energy*.pdf"));

}

void niceLogBins(TH1*h)
{
  TAxis *axis = h->GetXaxis();
  int bins = axis->GetNbins();

  double from = axis->GetXmin();
  double to = axis->GetXmax();
  double width = (to - from) / bins;
  double *new_bins = new double[bins + 1];

  for (int i = 0; i <= bins; i++) {
    new_bins[i] = pow(10, from + i * width);
  }
  axis->Set(bins, new_bins);
  delete[] new_bins;
}
