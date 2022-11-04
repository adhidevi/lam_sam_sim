//This macro produces vertex distribution for various virtual detectors
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

  double kinEcut = 1;//MeV
  string spTit[] = {"e-/#pi-","e+/#pi+","#gamma","neutron","e-/e+","e- trid==1"};
  const int nSp = sizeof(spTit)/sizeof(*spTit);
  for(int iSp=0;iSp<nSp;iSp++){
    spTit[iSp] = Form("%s (KE>=%.0f MeV) ",spTit[iSp].c_str(),kinEcut);
  }
  const string spH[nSp] = {"epiM","epiP","g","n","ee","eTrIdCut"};
  map<int,int> spM {{11,1},{-211,1},{-11,2},{211,2},{22,3},{2112,4}};

///Change the following lines for which detectors you want to include////
  string detH[] = {"det176"};
//LAM
  const int nDet = sizeof(detH)/sizeof(*detH);
  const int Det[nDet] = {176};
  map<int,int> dtM {{176,1}};
////////////////////////////////////////////////////////////////////////

  const double x_min = -100;
  const double x_max = 100;
  const double z_min = -8000;
  const double z_max = 30000;
  const double r_min[nDet] = {53};
  const double r_max[nDet] = {63};
  const int nbin = 500;
  const string weight[] = {"rate"/*,"rateE"*/};
  const int nWt = sizeof(weight)/sizeof(*weight);
  TH2F* h_VxVz[nSp][nDet][nWt];
  TH2F* h_VxVzPzG0[nSp][nDet][nWt];
  TH2F* h_VxVzPzL0[nSp][nDet][nWt];
  TH2F* h_VyVz[nSp][nDet][nWt];
  TH2F* h_VyVzPzG0[nSp][nDet][nWt];
  TH2F* h_VyVzPzL0[nSp][nDet][nWt];
  TH2F* h_VrVz[nSp][nDet][nWt];
  TH2F* h_VrVzPzG0[nSp][nDet][nWt];
  TH2F* h_VrVzPzL0[nSp][nDet][nWt];
  TH2F* h_xy[nSp][nDet][nWt];
  TH2F* h_xyPzG0[nSp][nDet][nWt];
  TH2F* h_xyPzL0[nSp][nDet][nWt];
  TH1F* h_r[nSp][nDet][nWt];
  TH1F* h_rPzG0[nSp][nDet][nWt];
  TH1F* h_rPzL0[nSp][nDet][nWt];

///Change the following lines as needed////
  const string geometry = "sam_study";
  const string tgt_gen_config = "sam_symm_study";
  const string plotType = Form("vertex_distribution_rmin53_rmax63_EG%.0f",kinEcut);
  int beamGen(1);

///Change this line for appropriate rootfile directory////
  TString rootfile_dir = "/volatile/halla/moller12gev/amgunsch/remoll_rootfiles/FieldMap_07_2022";
//////////////////////////////////////////////////////////

  for(int iSp=0;iSp<nSp;iSp++){
   for(int iDet=0;iDet<nDet;iDet++){
    for(int iWt=0;iWt<nWt;iWt++){
      string titleXZ = Form("%s vertex dist. on %s plane (%s),%.0f<r<%.0f;vz (mm);vx (mm);hits/#thrownEvents",spTit[iSp].c_str(),detH[iDet].c_str(),tgt_gen_config.c_str(),r_min[iDet],r_max[iDet]);
      string titleYZ = Form("%s vertex dist. on %s plane (%s),%.0f<r<%.0f;vz (mm);vy (mm);hits/#thrownEvents",spTit[iSp].c_str(),detH[iDet].c_str(),tgt_gen_config.c_str(),r_min[iDet],r_max[iDet]);
      string titleRZ = Form("%s vertex dist. on %s plane (%s),%.0f<r<%.0f;vz (mm);vr (mm);hits/#thrownEvents",spTit[iSp].c_str(),detH[iDet].c_str(),tgt_gen_config.c_str(),r_min[iDet],r_max[iDet]);
      string titleXY = Form("%s xy dist. on %s plane (%s),%.0f<r<%.0f;x (mm);y (mm);hits/#thrownEvents",spTit[iSp].c_str(),detH[iDet].c_str(),tgt_gen_config.c_str(),r_min[iDet],r_max[iDet]);
      string titleR = Form("%s radial dist. on %s plane (%s),%.0f<r<%.0f;r (mm);hits/#thrownEvents",spTit[iSp].c_str(),detH[iDet].c_str(),tgt_gen_config.c_str(),r_min[iDet],r_max[iDet]);

      h_VxVz[iSp][iDet][iWt] = new TH2F(Form("%s_VxVz_%s_%s",detH[iDet].c_str(),spH[iSp].c_str(),weight[iWt].c_str()),titleXZ.c_str(),nbin,z_min,z_max,nbin,x_min,x_max);
      h_VxVzPzG0[iSp][iDet][iWt] = new TH2F(Form("%s_VxVzPzG0_%s_%s",detH[iDet].c_str(),spH[iSp].c_str(),weight[iWt].c_str()),titleXZ.c_str(),nbin,z_min,z_max,nbin,x_min,x_max);
      h_VxVzPzL0[iSp][iDet][iWt] = new TH2F(Form("%s_VxVzPzL0_%s_%s",detH[iDet].c_str(),spH[iSp].c_str(),weight[iWt].c_str()),titleXZ.c_str(),nbin,z_min,z_max,nbin,x_min,x_max);
      h_VyVz[iSp][iDet][iWt] = new TH2F(Form("%s_VyVz_%s_%s",detH[iDet].c_str(),spH[iSp].c_str(),weight[iWt].c_str()),titleYZ.c_str(),nbin,z_min,z_max,nbin,x_min,x_max);
      h_VyVzPzG0[iSp][iDet][iWt] = new TH2F(Form("%s_VyVzPzG0_%s_%s",detH[iDet].c_str(),spH[iSp].c_str(),weight[iWt].c_str()),titleYZ.c_str(),nbin,z_min,z_max,nbin,x_min,x_max);
      h_VyVzPzL0[iSp][iDet][iWt] = new TH2F(Form("%s_VyVzPzL0_%s_%s",detH[iDet].c_str(),spH[iSp].c_str(),weight[iWt].c_str()),titleYZ.c_str(),nbin,z_min,z_max,nbin,x_min,x_max);
      h_VrVz[iSp][iDet][iWt] = new TH2F(Form("%s_VrVz_%s_%s",detH[iDet].c_str(),spH[iSp].c_str(),weight[iWt].c_str()),titleRZ.c_str(),nbin,z_min,z_max,nbin,x_min,x_max);
      h_VrVzPzG0[iSp][iDet][iWt] = new TH2F(Form("%s_VrVzPzG0_%s_%s",detH[iDet].c_str(),spH[iSp].c_str(),weight[iWt].c_str()),titleRZ.c_str(),nbin,z_min,z_max,nbin,x_min,x_max);
      h_VrVzPzL0[iSp][iDet][iWt] = new TH2F(Form("%s_VrVzPzL0_%s_%s",detH[iDet].c_str(),spH[iSp].c_str(),weight[iWt].c_str()),titleRZ.c_str(),nbin,z_min,z_max,nbin,x_min,x_max);

      h_xy[iSp][iDet][iWt] = new TH2F(Form("%s_xy_%s_%s",detH[iDet].c_str(),spH[iSp].c_str(),weight[iWt].c_str()),titleXY.c_str(),nbin,x_min,x_max,nbin,x_min,x_max);
      h_xyPzG0[iSp][iDet][iWt] = new TH2F(Form("%s_xyPzG0_%s_%s",detH[iDet].c_str(),spH[iSp].c_str(),weight[iWt].c_str()),titleXY.c_str(),nbin,x_min,x_max,nbin,x_min,x_max);
      h_xyPzL0[iSp][iDet][iWt] = new TH2F(Form("%s_xyPzL0_%s_%s",detH[iDet].c_str(),spH[iSp].c_str(),weight[iWt].c_str()),titleXY.c_str(),nbin,x_min,x_max,nbin,x_min,x_max);
      h_r[iSp][iDet][iWt] = new TH1F(Form("%s_r_%s_%s",detH[iDet].c_str(),spH[iSp].c_str(),weight[iWt].c_str()),titleR.c_str(),nbin,0,x_max);
      h_rPzG0[iSp][iDet][iWt] = new TH1F(Form("%s_rPzG0_%s_%s",detH[iDet].c_str(),spH[iSp].c_str(),weight[iWt].c_str()),titleR.c_str(),nbin,0,x_max);
      h_rPzL0[iSp][iDet][iWt] = new TH1F(Form("%s_rPzL0_%s_%s",detH[iDet].c_str(),spH[iSp].c_str(),weight[iWt].c_str()),titleR.c_str(),nbin,0,x_max);
    }
   }
  }

  int nfile=0;
  Long64_t nentry=0;
  long nTotEv=0;
  for(int ifile=1001;ifile<=6000;ifile++){
///Change this line for appropriate rootfiles////
    string infile = Form("%s/%s/%s_%d.root",rootfile_dir.Data(),tgt_gen_config.c_str(),tgt_gen_config.c_str(),ifile);
/////////////////////////////////////////////
    TFile *fin = TFile::Open(infile.c_str(),"READ");
    if(!fin || fin->IsZombie()){
      delete fin;
      continue;
      cout<<Form("Skipping %s. Recovered file.",infile.c_str())<<endl;
    }
    if(fin->TestBit(TFile::kRecovered)){
      continue;
      cout<<Form("Skipping %s. Recovered file.",infile.c_str())<<endl;
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
        if(hit->at(j).r<r_min[dt] || hit->at(j).r>r_max[dt]) continue;
        if(hit->at(j).k<kinEcut) continue;
        
//        double vr = TMath::Sqrt(pow(hit->at(j).vx,2)+pow(hit->at(j).vy,2));
        double vr = TMath::Sqrt(hit->at(j).vx*hit->at(j).vx+hit->at(j).vy*hit->at(j).vy);

        h_VxVz[sp][dt][0]->Fill(hit->at(j).vz,hit->at(j).vx,rate);
//        h_VxVz[sp][dt][1]->Fill(hit->at(j).vz,hit->at(j).vx,rate*hit->at(j).e);
        h_VyVz[sp][dt][0]->Fill(hit->at(j).vz,hit->at(j).vy,rate);
//        h_VyVz[sp][dt][1]->Fill(hit->at(j).vz,hit->at(j).vy,rate*hit->at(j).e);
        h_VrVz[sp][dt][0]->Fill(hit->at(j).vz,vr,rate);
//        h_VrVz[sp][dt][1]->Fill(hit->at(j).vz,vr,rate*hit->at(j).e);
        h_xy[sp][dt][0]->Fill(hit->at(j).x,hit->at(j).y,rate);
//        h_xy[sp][dt][1]->Fill(hit->at(j).x,hit->at(j).y,rate*hit->at(j).e);
        h_r[sp][dt][0]->Fill(hit->at(j).r,rate);
//        h_r[sp][dt][1]->Fill(hit->at(j).r,rate*hit->at(j).e);

        if(hit->at(j).pz>=0){
          h_VxVzPzG0[sp][dt][0]->Fill(hit->at(j).vz,hit->at(j).vx,rate);
//          h_VxVzPzG0[sp][dt][1]->Fill(hit->at(j).vz,hit->at(j).vx,rate*hit->at(j).e);
          h_VyVzPzG0[sp][dt][0]->Fill(hit->at(j).vz,hit->at(j).vy,rate);
//          h_VyVzPzG0[sp][dt][1]->Fill(hit->at(j).vz,hit->at(j).vy,rate*hit->at(j).e);
          h_VrVzPzG0[sp][dt][0]->Fill(hit->at(j).vz,vr,rate);
//          h_VrVzPzG0[sp][dt][1]->Fill(hit->at(j).vz,vr,rate*hit->at(j).e);
          h_xyPzG0[sp][dt][0]->Fill(hit->at(j).x,hit->at(j).y,rate);
//          h_xyPzG0[sp][dt][1]->Fill(hit->at(j).x,hit->at(j).y,rate*hit->at(j).e);
          h_rPzG0[sp][dt][0]->Fill(hit->at(j).r,rate);
//          h_rPzG0[sp][dt][1]->Fill(hit->at(j).r,rate*hit->at(j).e);
        }else{
          h_VxVzPzL0[sp][dt][0]->Fill(hit->at(j).vz,hit->at(j).vx,rate);
//          h_VxVzPzL0[sp][dt][1]->Fill(hit->at(j).vz,hit->at(j).vx,rate*hit->at(j).e);
          h_VyVzPzL0[sp][dt][0]->Fill(hit->at(j).vz,hit->at(j).vy,rate);
//          h_VyVzPzL0[sp][dt][1]->Fill(hit->at(j).vz,hit->at(j).vy,rate*hit->at(j).e);
          h_VrVzPzL0[sp][dt][0]->Fill(hit->at(j).vz,vr,rate);
//          h_VrVzPzL0[sp][dt][1]->Fill(hit->at(j).vz,vr,rate*hit->at(j).e);
          h_xyPzL0[sp][dt][0]->Fill(hit->at(j).x,hit->at(j).y,rate);
//          h_xyPzL0[sp][dt][1]->Fill(hit->at(j).x,hit->at(j).y,rate*hit->at(j).e);
          h_rPzL0[sp][dt][0]->Fill(hit->at(j).r,rate);
//          h_rPzL0[sp][dt][1]->Fill(hit->at(j).r,rate*hit->at(j).e);
        }

        if(hit->at(j).pid==11 || hit->at(j).pid==-11){
          h_VxVz[4][dt][0]->Fill(hit->at(j).vz,hit->at(j).vx,rate);
//          h_VxVz[4][dt][1]->Fill(hit->at(j).vz,hit->at(j).vx,rate*hit->at(j).e);
          h_VyVz[4][dt][0]->Fill(hit->at(j).vz,hit->at(j).vy,rate);
//          h_VyVz[4][dt][1]->Fill(hit->at(j).vz,hit->at(j).vy,rate*hit->at(j).e);
          h_VrVz[4][dt][0]->Fill(hit->at(j).vz,vr,rate);
//          h_VrVz[4][dt][1]->Fill(hit->at(j).vz,vr,rate*hit->at(j).e);
          h_xy[4][dt][0]->Fill(hit->at(j).x,hit->at(j).y,rate);
//          h_xy[4][dt][1]->Fill(hit->at(j).x,hit->at(j).y,rate*hit->at(j).e);
          h_r[4][dt][0]->Fill(hit->at(j).r,rate);
//          h_r[4][dt][1]->Fill(hit->at(j).r,rate*hit->at(j).e);
          if(hit->at(j).pz>=0){
            h_VxVzPzG0[4][dt][0]->Fill(hit->at(j).vz,hit->at(j).vx,rate);
//            h_VxVzPzG0[4][dt][1]->Fill(hit->at(j).vz,hit->at(j).vx,rate*hit->at(j).e);
            h_VyVzPzG0[4][dt][0]->Fill(hit->at(j).vz,hit->at(j).vy,rate);
//            h_VyVzPzG0[4][dt][1]->Fill(hit->at(j).vz,hit->at(j).vy,rate*hit->at(j).e);
            h_VrVzPzG0[4][dt][0]->Fill(hit->at(j).vz,vr,rate);
//            h_VrVzPzG0[4][dt][1]->Fill(hit->at(j).vz,vr,rate*hit->at(j).e);
            h_xyPzG0[4][dt][0]->Fill(hit->at(j).x,hit->at(j).y,rate);
//            h_xyPzG0[4][dt][1]->Fill(hit->at(j).x,hit->at(j).y,rate*hit->at(j).e);
            h_rPzG0[4][dt][0]->Fill(hit->at(j).r,rate);
//            h_rPzG0[4][dt][1]->Fill(hit->at(j).r,rate*hit->at(j).e);
          }else{
            h_VxVzPzL0[4][dt][0]->Fill(hit->at(j).vz,hit->at(j).vx,rate);
//            h_VxVzPzL0[4][dt][1]->Fill(hit->at(j).vz,hit->at(j).vx,rate*hit->at(j).e);
            h_VyVzPzL0[4][dt][0]->Fill(hit->at(j).vz,hit->at(j).vy,rate);
//            h_VyVzPzL0[4][dt][1]->Fill(hit->at(j).vz,hit->at(j).vy,rate*hit->at(j).e);
            h_VrVzPzL0[4][dt][0]->Fill(hit->at(j).vz,vr,rate);
//            h_VrVzPzL0[4][dt][1]->Fill(hit->at(j).vz,vr,rate*hit->at(j).e);
            h_xyPzL0[4][dt][0]->Fill(hit->at(j).x,hit->at(j).y,rate);
//            h_xyPzL0[4][dt][1]->Fill(hit->at(j).x,hit->at(j).y,rate*hit->at(j).e);
            h_rPzL0[4][dt][0]->Fill(hit->at(j).r,rate);
//            h_rPzL0[4][dt][1]->Fill(hit->at(j).r,rate*hit->at(j).e);
          }
        }

        if(hit->at(j).trid==1 && hit->at(j).pid==11){
          h_VxVz[5][dt][0]->Fill(hit->at(j).vz,hit->at(j).vx,rate);
//          h_VxVz[5][dt][1]->Fill(hit->at(j).vz,hit->at(j).vx,rate*hit->at(j).e);
          h_VyVz[5][dt][0]->Fill(hit->at(j).vz,hit->at(j).vy,rate);
//          h_VyVz[5][dt][1]->Fill(hit->at(j).vz,hit->at(j).vy,rate*hit->at(j).e);
          h_VrVz[5][dt][0]->Fill(hit->at(j).vz,vr,rate);
//          h_VrVz[5][dt][1]->Fill(hit->at(j).vz,vr,rate*hit->at(j).e);
          h_xy[5][dt][0]->Fill(hit->at(j).x,hit->at(j).y,rate);
//          h_xy[5][dt][1]->Fill(hit->at(j).x,hit->at(j).y,rate*hit->at(j).e);
          h_r[5][dt][0]->Fill(hit->at(j).r,rate);
//          h_r[5][dt][1]->Fill(hit->at(j).r,rate*hit->at(j).e);
          if(hit->at(j).pz>=0){
            h_VxVzPzG0[5][dt][0]->Fill(hit->at(j).vz,hit->at(j).vx,rate);
//            h_VxVzPzG0[5][dt][1]->Fill(hit->at(j).vz,hit->at(j).vx,rate*hit->at(j).e);
            h_VyVzPzG0[5][dt][0]->Fill(hit->at(j).vz,hit->at(j).vy,rate);
//            h_VyVzPzG0[5][dt][1]->Fill(hit->at(j).vz,hit->at(j).vy,rate*hit->at(j).e);
            h_VrVzPzG0[5][dt][0]->Fill(hit->at(j).vz,vr,rate);
//            h_VrVzPzG0[5][dt][1]->Fill(hit->at(j).vz,vr,rate*hit->at(j).e);
            h_xyPzG0[5][dt][0]->Fill(hit->at(j).x,hit->at(j).y,rate);
//            h_xyPzG0[5][dt][1]->Fill(hit->at(j).x,hit->at(j).y,rate*hit->at(j).e);
            h_rPzG0[5][dt][0]->Fill(hit->at(j).r,rate);
//            h_rPzG0[5][dt][1]->Fill(hit->at(j).r,rate*hit->at(j).e);
          }else{
            h_VxVzPzG0[5][dt][0]->Fill(hit->at(j).vz,hit->at(j).vx,rate);
//            h_VxVzPzG0[5][dt][1]->Fill(hit->at(j).vz,hit->at(j).vx,rate*hit->at(j).e);
            h_VyVzPzG0[5][dt][0]->Fill(hit->at(j).vz,hit->at(j).vy,rate);
//            h_VyVzPzG0[5][dt][1]->Fill(hit->at(j).vz,hit->at(j).vy,rate*hit->at(j).e);
            h_VrVzPzG0[5][dt][0]->Fill(hit->at(j).vz,vr,rate);
//            h_VrVzPzG0[5][dt][1]->Fill(hit->at(j).vz,vr,rate*hit->at(j).e);
            h_xyPzL0[5][dt][0]->Fill(hit->at(j).x,hit->at(j).y,rate);
//            h_xyPzL0[5][dt][1]->Fill(hit->at(j).x,hit->at(j).y,rate*hit->at(j).e);
            h_rPzL0[5][dt][0]->Fill(hit->at(j).r,rate);
//            h_rPzL0[5][dt][1]->Fill(hit->at(j).r,rate*hit->at(j).e);
          }
        }

      }
    }
    delete T;
    delete fin;
  }

  cout<<Form("Total number of file splits: %d",nfile)<<endl;
  cout<<Form("Total number of entries: %ld",nTotEv)<<endl;
   
//////////////////////////////////////////
  TFile* outfile = new TFile(Form("./rootfiles/%s_%s_%s_rCut.root",geometry.c_str(),tgt_gen_config.c_str(),plotType.c_str()),"recreate");
  for(int iDet=0;iDet<nDet;iDet++){
    outfile->mkdir(Form("%s",detH[iDet].c_str()));
    outfile->cd(Form("%s",detH[iDet].c_str()));
    for(int iSp=0;iSp<nSp;iSp++){
      for(int iWt=0;iWt<nWt;iWt++){
       if(beamGen){
         h_VxVz[iSp][iDet][iWt]->Scale(1.0/nTotEv);
         h_VxVzPzG0[iSp][iDet][iWt]->Scale(1.0/nTotEv);
         h_VxVzPzL0[iSp][iDet][iWt]->Scale(1.0/nTotEv);
         h_VyVz[iSp][iDet][iWt]->Scale(1.0/nTotEv);
         h_VyVzPzG0[iSp][iDet][iWt]->Scale(1.0/nTotEv);
         h_VyVzPzL0[iSp][iDet][iWt]->Scale(1.0/nTotEv);
         h_VrVz[iSp][iDet][iWt]->Scale(1.0/nTotEv);
         h_VrVzPzG0[iSp][iDet][iWt]->Scale(1.0/nTotEv);
         h_VrVzPzL0[iSp][iDet][iWt]->Scale(1.0/nTotEv);
         h_xy[iSp][iDet][iWt]->Scale(1.0/nTotEv);
         h_xyPzG0[iSp][iDet][iWt]->Scale(1.0/nTotEv);
         h_xyPzL0[iSp][iDet][iWt]->Scale(1.0/nTotEv);
         h_r[iSp][iDet][iWt]->Scale(1.0/nTotEv);
         h_rPzG0[iSp][iDet][iWt]->Scale(1.0/nTotEv);
         h_rPzL0[iSp][iDet][iWt]->Scale(1.0/nTotEv);
       }else{
         h_VxVz[iSp][iDet][iWt]->Scale(1.0e-9/nfile);
         h_VxVzPzG0[iSp][iDet][iWt]->Scale(1.0e-9/nfile);
         h_VxVzPzL0[iSp][iDet][iWt]->Scale(1.0e-9/nfile);
         h_VyVz[iSp][iDet][iWt]->Scale(1.0e-9/nfile);
         h_VyVzPzG0[iSp][iDet][iWt]->Scale(1.0e-9/nfile);
         h_VyVzPzL0[iSp][iDet][iWt]->Scale(1.0e-9/nfile);
         h_VrVz[iSp][iDet][iWt]->Scale(1.0e-9/nfile);
         h_VrVzPzG0[iSp][iDet][iWt]->Scale(1.0e-9/nfile);
         h_VrVzPzL0[iSp][iDet][iWt]->Scale(1.0e-9/nfile);
         h_xy[iSp][iDet][iWt]->Scale(1.0e-9/nfile);
         h_xyPzG0[iSp][iDet][iWt]->Scale(1.0e-9/nfile);
         h_xyPzL0[iSp][iDet][iWt]->Scale(1.0e-9/nfile);
         h_r[iSp][iDet][iWt]->Scale(1.0e-9/nfile);
         h_rPzG0[iSp][iDet][iWt]->Scale(1.0e-9/nfile);
         h_rPzL0[iSp][iDet][iWt]->Scale(1.0e-9/nfile);
       }
         h_VxVz[iSp][iDet][iWt]->Write();
         h_VxVzPzG0[iSp][iDet][iWt]->Write();
         h_VxVzPzL0[iSp][iDet][iWt]->Write();
         h_VyVz[iSp][iDet][iWt]->Write();
         h_VyVzPzG0[iSp][iDet][iWt]->Write();
         h_VyVzPzL0[iSp][iDet][iWt]->Write();
         h_VrVz[iSp][iDet][iWt]->Write();
         h_VrVzPzG0[iSp][iDet][iWt]->Write();
         h_VrVzPzL0[iSp][iDet][iWt]->Write();
         h_xy[iSp][iDet][iWt]->Write();
         h_xyPzG0[iSp][iDet][iWt]->Write();
         h_xyPzL0[iSp][iDet][iWt]->Write();
         h_r[iSp][iDet][iWt]->Write();
         h_rPzG0[iSp][iDet][iWt]->Write();
         h_rPzL0[iSp][iDet][iWt]->Write();
      }
    }
  }
  outfile->Close();
}
