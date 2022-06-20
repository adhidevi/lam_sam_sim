//This macro produces radial and transverse distribution for various virtual detectors
#include "remolltypes.hh"
#include <sstream>
#include <iostream>
#include <fstream>
void radial_trans_radialCut_EG1(){
  gROOT->Reset();
  gStyle->SetOptStat(0);
  gStyle->SetTitleYOffset(1.3);
  gStyle->SetPadGridX(1);
  gStyle->SetPadGridY(1);
  TGaxis::SetMaxDigits(3);
  
  double kinEcut = 0;//MeV
  const string spTit[] = {"e-/#pi-","e+/#pi+","#gamma","neutron","e-/e+","e-/#pi- vz<=-3875","e- trid==1"};
  const int nSp = sizeof(spTit)/sizeof(*spTit);
  const string spH[nSp] = {"epiM","epiP","g","n","ee","epiMvzCut","eTrIdCut"};
  map<int,int> spM {{11,1},{-211,1},{-11,2},{211,2},{22,3},{2112,4}};

///Change the following lines for which detectors you want to include////
  string detH[] = {"det176","det177"};
  const int nDet = sizeof(detH)/sizeof(*detH);
  const int Det[nDet] = {176,177};
  map<int,int> dtM {{176,1},{177,2}};
////////////////////////////////////////////////////////////////////////

  double x_min = 0;
  double x_max = 500;
  const int nbinR = 500;
  const int nbinXY = 1000;
  double binR_width = (x_max-x_min)/nbinR;
  const string weight[] = {"rate"/*,"rateE"*/};
  const int nWt = sizeof(weight)/sizeof(*weight);
// if physics generators use the following
//  const string weight_unit[nWt] ={"rate (GHz)","rate*E (GHz*MeV)"};
// if beam generator use the following
  const string weight_unit[nWt] ={"hits/#thrownEvents"/*,"E*hits/#thrownEvents"*/};

  TH1F* h_rate[nSp][nDet][nWt];
  TH1F* h_ratePzG0[nSp][nDet][nWt];
  TH1F* h_ratePzL0[nSp][nDet][nWt];
  TH2F* h_xy[nSp][nDet][nWt];
  TH2F* h_xyPzG0[nSp][nDet][nWt];
  TH2F* h_xyPzL0[nSp][nDet][nWt];

///Change the following lines as needed////
  const string geometry = "develop_new";
  const string tgt_gen_config = "LH2_beam_V1";
  const string plotType = Form("radial_trans_sam_EG%.0f",kinEcut);
  int beamGen(1);

///Change this line for appropriate rootfile directory////
  TString rootfile_dir = "/volatile/halla/parity/adhidevi/remoll_rootfiles/develop_br";
//////////////////////////////////////////////////////////

  for(int iSp=0;iSp<nSp;iSp++){
   for(int iDet=0;iDet<nDet;iDet++){
    for(int iWt=0;iWt<nWt;iWt++){
      string title1D = Form("%s (KE>=%.0f MeV) Radial dist. on %s plane (%s);Radius (mm);%s/%.1fmm",spTit[iSp].c_str(),kinEcut,detH[iDet].c_str(),tgt_gen_config.c_str(),weight_unit[iWt].c_str(),binR_width);
      h_rate[iSp][iDet][iWt] = new TH1F(Form("%s_r_%s_%s",detH[iDet].c_str(),spH[iSp].c_str(),weight[iWt].c_str()),title1D.c_str(),nbinR,x_min,x_max);
      h_ratePzG0[iSp][iDet][iWt] = new TH1F(Form("%s_rPzG0_%s_%s",detH[iDet].c_str(),spH[iSp].c_str(),weight[iWt].c_str()),title1D.c_str(),nbinR,x_min,x_max);
      h_ratePzL0[iSp][iDet][iWt] = new TH1F(Form("%s_rPzL0_%s_%s",detH[iDet].c_str(),spH[iSp].c_str(),weight[iWt].c_str()),title1D.c_str(),nbinR,x_min,x_max);

      string title2D = Form("%s (KE>=%.0f MeV) XY dist. on %s plane (%s);x (mm);y (mm)",spTit[iSp].c_str(),kinEcut,detH[iDet].c_str(),tgt_gen_config.c_str());
      h_xy[iSp][iDet][iWt] = new TH2F(Form("%s_xy_%s_%s",detH[iDet].c_str(),spH[iSp].c_str(),weight[iWt].c_str()),title2D.c_str(),nbinXY,-x_max,x_max,nbinXY,-x_max,x_max);
      h_xyPzG0[iSp][iDet][iWt] = new TH2F(Form("%s_xyPzG0_%s_%s",detH[iDet].c_str(),spH[iSp].c_str(),weight[iWt].c_str()),title2D.c_str(),nbinXY,-x_max,x_max,nbinXY,-x_max,x_max);
      h_xyPzL0[iSp][iDet][iWt] = new TH2F(Form("%s_xyPzL0_%s_%s",detH[iDet].c_str(),spH[iSp].c_str(),weight[iWt].c_str()),title2D.c_str(),nbinXY,-x_max,x_max,nbinXY,-x_max,x_max);

      h_rate[iSp][iDet][iWt]->Sumw2();
      h_ratePzG0[iSp][iDet][iWt]->Sumw2();
      h_ratePzL0[iSp][iDet][iWt]->Sumw2();
    }
   }
  }

  int nfile=0;
  Long64_t nentry=0;
  long nTotEv=0;
  for(int ifile=1001;ifile<=6000;ifile++){
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
//        if(hit->at(j).r>100) continue;
        if(hit->at(j).k<kinEcut) continue;

        h_rate[sp][dt][0]->Fill(hit->at(j).r,rate);
//        h_rate[sp][dt][1]->Fill(hit->at(j).r,rate*hit->at(j).e);
        h_xy[sp][dt][0]->Fill(hit->at(j).x,hit->at(j).y,rate);
//        h_xy[sp][dt][1]->Fill(hit->at(j).x,hit->at(j).y,rate*hit->at(j).e);

        if(hit->at(j).pz>=0){
          h_ratePzG0[sp][dt][0]->Fill(hit->at(j).r,rate);
//          h_ratePzG0[sp][dt][1]->Fill(hit->at(j).r,rate*hit->at(j).e);
          h_xyPzG0[sp][dt][0]->Fill(hit->at(j).x,hit->at(j).y,rate);
//          h_xyPzG0[sp][dt][1]->Fill(hit->at(j).x,hit->at(j).y,rate*hit->at(j).e);
        }else{
          h_ratePzL0[sp][dt][0]->Fill(hit->at(j).r,rate);
//          h_ratePzL0[sp][dt][1]->Fill(hit->at(j).r,rate*hit->at(j).e);
          h_xyPzL0[sp][dt][0]->Fill(hit->at(j).x,hit->at(j).y,rate);
//          h_xyPzL0[sp][dt][1]->Fill(hit->at(j).x,hit->at(j).y,rate*hit->at(j).e);
        }

        if(hit->at(j).pid==11 || hit->at(j).pid==-11){
          h_rate[4][dt][0]->Fill(hit->at(j).r,rate);
//          h_rate[4][dt][1]->Fill(hit->at(j).r,rate*hit->at(j).e);
          h_xy[4][dt][0]->Fill(hit->at(j).x,hit->at(j).y,rate);
//          h_xy[4][dt][1]->Fill(hit->at(j).x,hit->at(j).y,rate*hit->at(j).e);
          if(hit->at(j).pz>=0){
            h_ratePzG0[4][dt][0]->Fill(hit->at(j).r,rate);
//            h_ratePzG0[4][dt][1]->Fill(hit->at(j).r,rate*hit->at(j).e);
            h_xyPzG0[4][dt][0]->Fill(hit->at(j).x,hit->at(j).y,rate);
//            h_xyPzG0[4][dt][1]->Fill(hit->at(j).x,hit->at(j).y,rate*hit->at(j).e);
          }else{
            h_ratePzL0[4][dt][0]->Fill(hit->at(j).r,rate);
//            h_ratePzL0[4][dt][1]->Fill(hit->at(j).r,rate*hit->at(j).e);
            h_xyPzL0[4][dt][0]->Fill(hit->at(j).x,hit->at(j).y,rate);
//            h_xyPzL0[4][dt][1]->Fill(hit->at(j).x,hit->at(j).y,rate*hit->at(j).e);
          }
        }

        if(hit->at(j).vz<=-3875 && (hit->at(j).pid==11 || hit->at(j).pid==-211)){
          h_rate[5][dt][0]->Fill(hit->at(j).r,rate);
//          h_rate[5][dt][1]->Fill(hit->at(j).r,rate*hit->at(j).e);
          h_xy[5][dt][0]->Fill(hit->at(j).x,hit->at(j).y,rate);
//          h_xy[5][dt][1]->Fill(hit->at(j).x,hit->at(j).y,rate*hit->at(j).e);
          if(hit->at(j).pz>=0){
            h_ratePzG0[5][dt][0]->Fill(hit->at(j).r,rate);
//            h_ratePzG0[5][dt][1]->Fill(hit->at(j).r,rate*hit->at(j).e);
            h_xyPzG0[5][dt][0]->Fill(hit->at(j).x,hit->at(j).y,rate);
//            h_xyPzG0[5][dt][1]->Fill(hit->at(j).x,hit->at(j).y,rate*hit->at(j).e);
          }else{
            h_ratePzL0[5][dt][0]->Fill(hit->at(j).r,rate);
//            h_ratePzL0[5][dt][1]->Fill(hit->at(j).r,rate*hit->at(j).e);
            h_xyPzG0[5][dt][0]->Fill(hit->at(j).x,hit->at(j).y,rate);
//            h_xyPzG0[5][dt][1]->Fill(hit->at(j).x,hit->at(j).y,rate*hit->at(j).e);
          }
        }
        if(hit->at(j).trid==1 && hit->at(j).pid==11){
          h_rate[6][dt][0]->Fill(hit->at(j).r,rate);
//          h_rate[6][dt][1]->Fill(hit->at(j).r,rate*hit->at(j).e);
          h_xy[6][dt][0]->Fill(hit->at(j).x,hit->at(j).y,rate);
//          h_xy[6][dt][1]->Fill(hit->at(j).x,hit->at(j).y,rate*hit->at(j).e);
          if(hit->at(j).pz>=0){
            h_ratePzG0[6][dt][0]->Fill(hit->at(j).r,rate);
//            h_ratePzG0[6][dt][1]->Fill(hit->at(j).r,rate*hit->at(j).e);
            h_xyPzG0[6][dt][0]->Fill(hit->at(j).x,hit->at(j).y,rate);
//            h_xyPzG0[6][dt][1]->Fill(hit->at(j).x,hit->at(j).y,rate*hit->at(j).e);
          }else{
            h_ratePzL0[6][dt][0]->Fill(hit->at(j).r,rate);
//            h_ratePzL0[6][dt][1]->Fill(hit->at(j).r,rate*hit->at(j).e);
            h_xyPzG0[6][dt][0]->Fill(hit->at(j).x,hit->at(j).y,rate);
//            h_xyPzG0[6][dt][1]->Fill(hit->at(j).x,hit->at(j).y,rate*hit->at(j).e);
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
  TFile* outfile = new TFile(Form("./rootfiles/%s_%s_%s.root",geometry.c_str(),tgt_gen_config.c_str(),plotType.c_str()),"recreate");

  for(int iDet=0;iDet<nDet;iDet++){
    outfile->mkdir(Form("%s",detH[iDet].c_str()));
    outfile->cd(Form("%s",detH[iDet].c_str()));
    for(int iSp=0;iSp<nSp;iSp++){
      for(int iWt=0;iWt<nWt;iWt++){
       if(beamGen){
         h_rate[iSp][iDet][iWt]->Scale(1.0/nTotEv);
         h_ratePzG0[iSp][iDet][iWt]->Scale(1.0/nTotEv);
         h_ratePzL0[iSp][iDet][iWt]->Scale(1.0/nTotEv);
         h_xy[iSp][iDet][iWt]->Scale(1.0/nTotEv);
         h_xyPzG0[iSp][iDet][iWt]->Scale(1.0/nTotEv);
         h_xyPzL0[iSp][iDet][iWt]->Scale(1.0/nTotEv);
       }else{
         h_rate[iSp][iDet][iWt]->Scale(1.0e-9/nfile);//convert to GHz with 1.0e-9
         h_ratePzG0[iSp][iDet][iWt]->Scale(1.0e-9/nfile);
         h_ratePzL0[iSp][iDet][iWt]->Scale(1.0e-9/nfile);
         h_xy[iSp][iDet][iWt]->Scale(1.0e-9/nfile);
         h_xyPzG0[iSp][iDet][iWt]->Scale(1.0e-9/nfile);
         h_xyPzL0[iSp][iDet][iWt]->Scale(1.0e-9/nfile);
       }
         h_rate[iSp][iDet][iWt]->Write();
         h_ratePzG0[iSp][iDet][iWt]->Write();
         h_ratePzL0[iSp][iDet][iWt]->Write();
         h_xy[iSp][iDet][iWt]->Write();
         h_xyPzG0[iSp][iDet][iWt]->Write();
         h_xyPzL0[iSp][iDet][iWt]->Write();
      }
    }
  }
}
