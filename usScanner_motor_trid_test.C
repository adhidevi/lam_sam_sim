#include "remolltypes.hh"
#include <sstream>
#include <iostream>
#include <fstream>

TString rootfile_dir = "/volatile/halla/parity/adhidevi/remoll_rootfiles/develop_br";
  const string spTit[] = {"e-/#pi-","e+/#pi+","#gamma","neutron"};
  const int nSp = sizeof(spTit)/sizeof(*spTit);
  const string spH[nSp] = {"epiM","epiP","g","n"};
  map<int,int> spM {{11,1},{-211,1},{-11,2},{211,2},{22,3},{2112,4}};
  string detH[] = {"det1175","det2175"};
  map<int,int> dtM {{1175,1},{2175,2}};
  const int nDet = sizeof(detH)/sizeof(*detH);
  const int Det[nDet] = {1175,2175};
  TH1F* h_d_r[nSp][nDet];
  TH2F* h_xy[nSp][nDet];
  int beamGen(1);
  const string geometry = "develop";
  const string tgt_gen_config = "LH2_beam_V19";
  const string plotType = Form("usScanner_motor_trid_NovzCut_EG0");

void usScanner_motor_trid_test(){
  gStyle->SetOptStat(0);
  for(int iSp=0;iSp<nSp;iSp++){
   for(int iDet=0;iDet<nDet;iDet++){
     string titleR = Form("%s Radial dist. on %s plane (%s);Radius (mm);hits/%sthrownEvents/%.1fmm",spTit[iSp].c_str(),detH[iDet].c_str(),tgt_gen_config.c_str(),"#",(3000.0-0.0)/1000);
     h_d_r[iSp][iDet] = new TH1F(Form("%s_r_%s_rate",detH[iDet].c_str(),spH[iSp].c_str()),titleR.c_str(),1000,0,3000);

     string titleXY = Form("%s XY dist. on %s plane (%s);x (mm);y (mm)",spTit[iSp].c_str(),detH[iDet].c_str(),tgt_gen_config.c_str());
     h_xy[iSp][iDet] = new TH2F(Form("%s_xy_%s_rate",detH[iDet].c_str(),spH[iSp].c_str()),titleXY.c_str(),1000,-3000,3000,1000,-3000,3000);
     h_d_r[iSp][iDet]->Sumw2();
   }
  }

  int nfile=0;
  Long64_t nentry=0;
  long nTotEv=0;
  for(int ifile=1001;ifile<=1100;ifile++){
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
      continue;
    }
    nfile++;
   
    TTree *T = (TTree*)fin->Get("T");
    if(T==0) return 0;

    nentry = T->GetEntries();
    cout<<Form("Found %lld entries in filesplit %d",nentry,nfile)<<endl;
    nTotEv+=nentry;

    std::vector<remollGenericDetectorHit_t> *hit =0;
    Double_t rate = 0;
    T->SetBranchAddress("hit", & hit);
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

        h_d_r[sp][dt]->Fill(hit->at(j).r,rate);
        h_xy[sp][dt]->Fill(hit->at(j).x,hit->at(j).y,rate);
      }
    }
    delete T;
    delete fin;
  }

  cout<<Form("Total number of file splits: %d",nfile)<<endl;
  cout<<Form("Total number of entries: %ld",nTotEv)<<endl;
//////////////////////////////////////////
  TFile* outfile = new TFile(Form("./rootfiles/%s_%s_%s_NoPhiCut_test.root",geometry.c_str(),tgt_gen_config.c_str(),plotType.c_str()),"recreate");
//////////////////////////////////////////   
  for(int iDet=0;iDet<nDet;iDet++){
    outfile->mkdir(Form("%s",detH[iDet].c_str()));
    outfile->cd(Form("%s",detH[iDet].c_str()));
    for(int iSp=0;iSp<nSp;iSp++){
      h_d_r[iSp][iDet]->Scale(1.0/nTotEv);
      h_xy[iSp][iDet]->Scale(1.0/nTotEv);

      h_d_r[iSp][iDet]->Write("",TObject::kOverwrite);
      h_xy[iSp][iDet]->Write("",TObject::kOverwrite);
    }
  }
}
