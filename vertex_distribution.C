//This macro produces radial and transverse distribution for various virtual detectors
#include "remolltypes.hh"
#include <sstream>
#include <iostream>
#include <fstream>
void vertex_distribution(){
  gROOT->Reset();
  gStyle->SetOptStat(0);
  gStyle->SetTitleYOffset(1.3);
  gStyle->SetPadGridX(1);
  gStyle->SetPadGridY(1);
  TGaxis::SetMaxDigits(3);
  
  const string spTit[] = {"e-/#pi- all E","e+/#pi+ all E","#gamma all E","neutron all E","e-/e+ all KE","primary all E"};
  const int nSp = sizeof(spTit)/sizeof(*spTit);
  const string spH[nSp] = {"epiM","epiP","g","n","ee","pri"};
  map<int,int> spM {{11,1},{-211,1},{-11,2},{211,2},{22,3},{2112,4}};

///Change the following lines for which detectors you want to include////
  string detH[] = {"det28"};
  const int nDet = sizeof(detH)/sizeof(*detH);
  const int Det[nDet] = {28};
  map<int,int> dtM {{28,1}};
////////////////////////////////////////////////////////////////////////

  const double x_min = -1500;
  const double x_max = 1500;
  const double z_min = -10000;
  const double z_max = 40000;
  const int nbin = 500;
  const double rmin = 200;
  const double rmax = 600;
  const string weight[] = {"rate","rateE"};
  const int nWt = sizeof(weight)/sizeof(*weight);
  TH2F* h_VxVz[nSp][nDet][nWt];
  TH2F* h_VxVzPzG0[nSp][nDet][nWt];
  TH2F* h_VxVzPzL0[nSp][nDet][nWt];

///Change the following lines as needed////
  const string geometry = "defaultGeo";//defaultGeo or PMTSh
  const string tgt_gen_config = "LH2_beam";
  const string plotType = "vertex_distribution_allE";
  int beamGen(1);

//////////////////////////////////////////

  TFile* outfile = new TFile(Form("./rootfiles/%s_%s_%s.root",geometry.c_str(),tgt_gen_config.c_str(),plotType.c_str()),"recreate");
///Change this line for appropriate rootfile directory////
  TString rootfile_dir = "/volatile/halla/parity/adhidevi/remoll_rootfiles/default-geo";
//////////////////////////////////////////////////////////

  for(int iSp=0;iSp<nSp;iSp++){
   for(int iDet=0;iDet<nDet;iDet++){
    for(int iWt=0;iWt<nWt;iWt++){
      string title2D = Form("%s vertex dist. on %s plane (%s), %.0f<=r=<%.0f;vz (mm);vx (mm);hits/#thrownEvents",spTit[iSp].c_str(),detH[iDet].c_str(),tgt_gen_config.c_str(),rmin,rmax);
      h_VxVz[iSp][iDet][iWt] = new TH2F(Form("%s_VxVz_%s_%s",detH[iDet].c_str(),spH[iSp].c_str(),weight[iWt].c_str()),title2D.c_str(),nbin,z_min,z_max,nbin,x_min,x_max);
      h_VxVzPzG0[iSp][iDet][iWt] = new TH2F(Form("%s_VxVzPzG0_%s_%s",detH[iDet].c_str(),spH[iSp].c_str(),weight[iWt].c_str()),title2D.c_str(),nbin,z_min,z_max,nbin,x_min,x_max);
      h_VxVzPzL0[iSp][iDet][iWt] = new TH2F(Form("%s_VxVzPzL0_%s_%s",detH[iDet].c_str(),spH[iSp].c_str(),weight[iWt].c_str()),title2D.c_str(),nbin,z_min,z_max,nbin,x_min,x_max);
    }
   }
  }

  int nfile=0;
  Long64_t nentry=0;
  long nTotEv=0;
  for(int ifile=1001;ifile<=1927;ifile++){
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
//comment following line if want to plot all r
        if(hit->at(j).r<200 || hit->at(j).r>600) continue;

        h_VxVz[sp][dt][0]->Fill(hit->at(j).vz,hit->at(j).vx,rate);
        h_VxVz[sp][dt][1]->Fill(hit->at(j).vz,hit->at(j).vx,rate*hit->at(j).e);

        if(hit->at(j).pz>=0){
          h_VxVzPzG0[sp][dt][0]->Fill(hit->at(j).vz,hit->at(j).vx,rate);
          h_VxVzPzG0[sp][dt][1]->Fill(hit->at(j).vz,hit->at(j).vx,rate*hit->at(j).e);
        }else{
          h_VxVzPzL0[sp][dt][0]->Fill(hit->at(j).vz,hit->at(j).vx,rate);
          h_VxVzPzL0[sp][dt][1]->Fill(hit->at(j).vz,hit->at(j).vx,rate*hit->at(j).e);
        }

        if(hit->at(j).pid==11 || hit->at(j).pid==-11){
          h_VxVz[4][dt][0]->Fill(hit->at(j).vz,hit->at(j).vx,rate);
          h_VxVz[4][dt][1]->Fill(hit->at(j).vz,hit->at(j).vx,rate*hit->at(j).e);
          if(hit->at(j).pz>=0){
            h_VxVzPzG0[4][dt][0]->Fill(hit->at(j).vz,hit->at(j).vx,rate);
            h_VxVzPzG0[4][dt][1]->Fill(hit->at(j).vz,hit->at(j).vx,rate*hit->at(j).e);
          }else{
            h_VxVzPzL0[4][dt][0]->Fill(hit->at(j).vz,hit->at(j).vx,rate);
            h_VxVzPzL0[4][dt][1]->Fill(hit->at(j).vz,hit->at(j).vx,rate*hit->at(j).e);
          }
        }

        if(hit->at(j).trid==1){
          h_VxVz[5][dt][0]->Fill(hit->at(j).vz,hit->at(j).vx,rate);
          h_VxVz[5][dt][1]->Fill(hit->at(j).vz,hit->at(j).vx,rate*hit->at(j).e);
          if(hit->at(j).pz>=0){
            h_VxVzPzG0[5][dt][0]->Fill(hit->at(j).vz,hit->at(j).vx,rate);
            h_VxVzPzG0[5][dt][1]->Fill(hit->at(j).vz,hit->at(j).vx,rate*hit->at(j).e);
          }else{
            h_VxVzPzG0[5][dt][0]->Fill(hit->at(j).vz,hit->at(j).vx,rate);
            h_VxVzPzG0[5][dt][1]->Fill(hit->at(j).vz,hit->at(j).vx,rate*hit->at(j).e);
          }
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
      for(int iWt=0;iWt<nWt;iWt++){
       if(beamGen){
         h_VxVz[iSp][iDet][iWt]->Scale(1.0/nTotEv);
         h_VxVzPzG0[iSp][iDet][iWt]->Scale(1.0/nTotEv);
         h_VxVzPzL0[iSp][iDet][iWt]->Scale(1.0/nTotEv);
       }else{
         h_VxVz[iSp][iDet][iWt]->Scale(1.0e-9/nfile);
         h_VxVzPzG0[iSp][iDet][iWt]->Scale(1.0e-9/nfile);
         h_VxVzPzL0[iSp][iDet][iWt]->Scale(1.0e-9/nfile);
       }
         h_VxVz[iSp][iDet][iWt]->Write();
         h_VxVzPzG0[iSp][iDet][iWt]->Write();
         h_VxVzPzL0[iSp][iDet][iWt]->Write();
      }
    }
  }
}
