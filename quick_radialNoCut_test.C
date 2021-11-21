//This macro produces radial distribution for various virtual detectors
#include "remolltypes.hh"
#include <sstream>
#include <iostream>
#include <fstream>
void quick_radialNoCut_test(){
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
  string detH[] = {"det28"};
  const int nDet = sizeof(detH)/sizeof(*detH);
  const int Det[nDet] = {28};
  map<int,int> dtM {{28,1}};
////////////////////////////////////////////////////////////////////////

  double x_min = 0;
  double x_max = 1500;
  const int nbin = 200;
  double bin_width = (x_max-x_min)/nbin;
  TH1F* h_rate[nSp][nDet];
  TH1F* h_ratePzG0[nSp][nDet];
  TH1F* h_ratePzL0[nSp][nDet];

///Change the following lines as needed////
  const string geometry = "PMTShielding_ep1";//defaultGeo or PMTSh
  const string tgt_gen_config = "remollout_ep";
  const string plotType = "radialNoCut_Ch";
  int beamGen(0);
//////////////////////////////////////////

  TFile* outfile = new TFile(Form("./rootfiles/%s_%s_%s.root",geometry.c_str(),tgt_gen_config.c_str(),plotType.c_str()),"recreate");
///Change this line for appropriate rootfile directory////
  TString rootfile_dir = "/volatile/halla/moller12gev/chandan/sim_out/";
//////////////////////////////////////////////////////////

  for(int iSp=0;iSp<nSp;iSp++){
   for(int iDet=0;iDet<nDet;iDet++){
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
  for(int ifile=1;ifile<=10;ifile++){
///Change this line for appropriate rootfiles////
    string infile = Form("%s/%s/%s%d.root",rootfile_dir.Data(),geometry.c_str(),tgt_gen_config.c_str(),ifile);
//////////////////////////////////////////////
    ifstream inf(infile.c_str());
    TFile *fin = TFile::Open(infile.c_str(),"READ");
    if(!inf){
      cout<<Form("Skipping %s. File diesn't exist.",infile.c_str())<<endl;
      continue;
    }
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
//        if(beamGen) rate = 1.0;

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
}
