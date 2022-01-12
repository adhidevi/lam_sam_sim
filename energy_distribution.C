#include "remolltypes.hh"
#include <sstream>
#include <iostream>
#include <fstream>
    TString rootfile_dir = "/volatile/halla/parity/adhidevi/remoll_rootfiles/PMTShielding/";
    const string spTit[] = {"e-/#pi-","e+/#pi+","#gamma","neutron","e-/e+ E>1","primary e E>1"};
    const int nSp = sizeof(spTit)/sizeof(*spTit);
    const string spH[nSp] ={"epiM","epiP","g","n","ee1","pri1"};
    map<int,int> spM {{11,1},{-211,1},{-11,2},{211,2},{22,3},{2112,4}};
    int color[nSp] = {2,3,1,4,6,7};
    const string detH[] = {"det28","det27"};
    map<int,int> dtM {{28,1},{27,2}};
    const int nDet = sizeof(detH)/sizeof(*detH);
    int Det[nDet] = {28,27};
    double rmin[nDet] = {1371.95,1371.95};
    double rmax[nDet] = {1651,1651};
    TH1D* eRate[nSp][nDet];
    void niceLogBins(TH1*);
    TFile* outfile; 
void energy_distribution(){
   gStyle->SetOptStat(0);
   for(int iSp=0;iSp<nSp;iSp++){
   for(int iDet=0;iDet<nDet;iDet++){
      eRate[iSp][iDet] = new TH1D(Form("%s_E_%s",detH[iDet].c_str(),spH[iSp].c_str()),Form("Kinetic Energy at PMT region for %s (beam generator);E (MeV);hits/%sthrownEvents",detH[iDet].c_str(),"#"),121,-8,5);
      eRate[iSp][iDet]->SetLineColor(color[iSp]);
      niceLogBins(eRate[iSp][iDet]);
   }
   }
    outfile = new TFile("./rootfiles/energy_dist_PMTSh_beam_V10.root","recreate");
    TChain* T = new TChain("T");
    int nfile = 0;
    for(int ifile=1001;ifile<=6000;ifile++){
       nfile++;
       T->Add(Form("%sPMTSh_beam_V10/PMTSh_beam_V10_%d.root",rootfile_dir.Data(),ifile));
    }
    cout<<Form("Found %d files!!!",nfile)<<endl;
    Long64_t nentry = T->GetEntries();
    cout<<Form("Total number of entries to be analyzed: %lld",nentry)<<endl;
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

          if(hit->at(j).r>=rmin[dt] && hit->at(j).r<=rmax[dt]){
            eRate[sp][dt]->Fill(hit->at(j).k,rate);
            if(hit->at(j).k>1 && (hit->at(j).pid==11 || hit->at(j).pid==-11)){
              eRate[4][dt]->Fill(hit->at(j).k,rate);
            }
            if(hit->at(j).k>1 && hit->at(j).trid==1){
              eRate[5][dt]->Fill(hit->at(j).k,rate);
            }
          }
       }
    }

    TCanvas* c_E[nDet];
    TLatex latex;
    latex.SetNDC(1);
    latex.SetTextSize(0.04);
    for(int iDet=0;iDet<nDet;iDet++){
       c_E[iDet] = new TCanvas(Form("c_E_%d",iDet));
       for(int iSp=0;iSp<nSp;iSp++){
          gPad->SetLogx();
          gPad->SetLogy();
          eRate[iSp][iDet]->Scale(1./nentry);
          eRate[iSp][iDet]->GetYaxis()->SetRangeUser(1.0e-9,1.0e-2);
          eRate[iSp][iDet]->SetMarkerStyle(20);
          eRate[iSp][iDet]->SetMarkerSize(0.5);
          eRate[iSp][iDet]->Draw("hist same");
          eRate[iSp][iDet]->Write();
          latex.SetTextColor(color[iSp]);
          if(iSp<3)
            latex.DrawLatex(0.33,0.85-0.05*iSp,Form("%s: %.2e",spTit[iSp].c_str(),eRate[iSp][iDet]->Integral()));
          else
            latex.DrawLatex(0.56,0.85-0.05*(iSp-3),Form("%s: %.2e",spTit[iSp].c_str(),eRate[iSp][iDet]->Integral()));
       }
    latex.SetTextColor(kRed+3);
    latex.DrawLatex(0.12,0.85,"Radial Cut::");
    latex.DrawLatex(0.12,0.80,Form("[%.0f,%.0f]mm",rmin[iDet],rmax[iDet]));
    c_E[iDet]->SaveAs(Form("./temp/kinetic_energy_det%d.pdf",Det[iDet]));
    }
    gSystem->Exec(Form("pdfunite ./temp/kinetic_energy_det*.pdf ./plots/PMTSh_beam_kinE_det%d_det%d.pdf",Det[0],Det[1]));
    gSystem->Exec(Form("rm -rf ./temp/kinetic_energy_det*.pdf"));

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
