//This macro produces radial and transverse distribution for various virtual detectors
#include "remolltypes.hh"
#include <sstream>
#include <iostream>
#include <fstream>
void radial_trans_radialCut_EG1_dsScanner(){
  gROOT->Reset();
  gStyle->SetOptStat(0);
  gStyle->SetTitleYOffset(1.3);
  gStyle->SetPadGridX(1);
  gStyle->SetPadGridY(1);
  TGaxis::SetMaxDigits(3);
  
//  const string spTit[] = {"e-/#pi- (KE>1 MeV)","e+/#pi+ (KE>1 MeV)","#gamma (KE>1 MeV)","neutron (KE>1 MeV)","e-/e+ (KE>1 MeV)","primary (KE>1 MeV)"};
  const string spTit[] = {"e-/#pi- all E","e+/#pi+ all E","#gamma all E","neutron all E","e-/e+ all KE","primary all E"};
  const int nSp = sizeof(spTit)/sizeof(*spTit);
  const string spH[nSp] = {"epiM","epiP","g","n","ee","pri"};
  map<int,int> spM {{11,1},{-211,1},{-11,2},{211,2},{22,3},{2112,4}};

///Change the following lines for which detectors you want to include////
<<<<<<< HEAD
  string detH[] = {"det28","det177"};
  const int nDet = sizeof(detH)/sizeof(*detH);
  const int Det[nDet] = {28,177};
  map<int,int> dtM {{28,1},{177,2}};
////////////////////////////////////////////////////////////////////////

  double x_min = 0;
  double x_max = 1000;
=======
  string detH[] = {"det174","det28","det30","det176"};
  const int nDet = sizeof(detH)/sizeof(*detH);
  const int Det[nDet] = {174,28,30,176};
  map<int,int> dtM {{174,1},{28,2},{30,3},{176,4}};
////////////////////////////////////////////////////////////////////////

  double x_min = 0;
  double x_max = 1500;
>>>>>>> 0c4366e67ad16fe5832ddc8a375d8e73cd22dbf0
  const int nbin = 200;
  double bin_width = (x_max-x_min)/nbin;
  const string weight[] = {"rate","rateE"};
  const int nWt = sizeof(weight)/sizeof(*weight);
<<<<<<< HEAD
//  if physics generators
//  const string weight_unit[nWt] ={"rate (GHz)","rate*E (GHz*MeV)","rate*A (GHz*ppb)"};
//if beam generator
  const string weight_unit[nWt] ={"hits/#thrownEvents","E*hits/#thrownEvents","A*hits/#thrownEvents"};
=======
  const string weight_unit[nWt] ={"hits/#thrownEvents","E*hits/#thrownEvents"};
>>>>>>> 0c4366e67ad16fe5832ddc8a375d8e73cd22dbf0
  TH1F* h_rate[nSp][nDet][nWt];
  TH1F* h_ratePzG0[nSp][nDet][nWt];
  TH1F* h_ratePzL0[nSp][nDet][nWt];
  TH2F* h_xy[nSp][nDet][nWt];
  TH2F* h_xyPzG0[nSp][nDet][nWt];
  TH2F* h_xyPzL0[nSp][nDet][nWt];

///Change the following lines as needed////
<<<<<<< HEAD
  const string geometry = "PMTSh";//defaultGeo or PMTSh
  const string tgt_gen_config = "PMTSh_beam_V6";
  const string plotType = "openSectors_radial_trans_rNoCut_allE";//rCut or rNoCut and EG1 or allE
=======
  const string geometry = "defaultGeo";//defaultGeo or PMTSh
  const string tgt_gen_config = "LH2_beam_magON_V5";
  const string plotType = "radial_trans_rNoCut_allE";//rCut or rNoCut and EG1 or allE
>>>>>>> 0c4366e67ad16fe5832ddc8a375d8e73cd22dbf0
  int beamGen(1);
//////////////////////////////////////////

  TFile* outfile = new TFile(Form("./rootfiles/%s_%s_%s.root",geometry.c_str(),tgt_gen_config.c_str(),plotType.c_str()),"recreate");
///Change this line for appropriate rootfile directory////
  TString rootfile_dir = "/volatile/halla/moller12gev/devi/remoll_rootfiles/default-geo";
//////////////////////////////////////////////////////////

  for(int iSp=0;iSp<nSp;iSp++){
   for(int iDet=0;iDet<nDet;iDet++){
    for(int iWt=0;iWt<nWt;iWt++){
      string title1D = Form("%s Radial dist. on %s plane (%s);Radius (mm);%s/%.1fmm",spTit[iSp].c_str(),detH[iDet].c_str(),tgt_gen_config.c_str(),weight_unit[iWt].c_str(),bin_width);
      h_rate[iSp][iDet][iWt] = new TH1F(Form("%s_r_%s_%s",detH[iDet].c_str(),spH[iSp].c_str(),weight[iWt].c_str()),title1D.c_str(),nbin,x_min,x_max);
      h_ratePzG0[iSp][iDet][iWt] = new TH1F(Form("%s_rPzG0_%s_%s",detH[iDet].c_str(),spH[iSp].c_str(),weight[iWt].c_str()),title1D.c_str(),nbin,x_min,x_max);
      h_ratePzL0[iSp][iDet][iWt] = new TH1F(Form("%s_rPzL0_%s_%s",detH[iDet].c_str(),spH[iSp].c_str(),weight[iWt].c_str()),title1D.c_str(),nbin,x_min,x_max);

      string title2D = Form("%s XY dist. on %s plane (%s);x (mm);y (mm)",spTit[iSp].c_str(),detH[iDet].c_str(),tgt_gen_config.c_str());
      h_xy[iSp][iDet][iWt] = new TH2F(Form("%s_xy_%s_%s",detH[iDet].c_str(),spH[iSp].c_str(),weight[iWt].c_str()),title2D.c_str(),nbin,-x_max,x_max,nbin,-x_max,x_max);
      h_xyPzG0[iSp][iDet][iWt] = new TH2F(Form("%s_xyPzG0_%s_%s",detH[iDet].c_str(),spH[iSp].c_str(),weight[iWt].c_str()),title2D.c_str(),nbin,-x_max,x_max,nbin,-x_max,x_max);
      h_xyPzL0[iSp][iDet][iWt] = new TH2F(Form("%s_xyPzL0_%s_%s",detH[iDet].c_str(),spH[iSp].c_str(),weight[iWt].c_str()),title2D.c_str(),nbin,-x_max,x_max,nbin,-x_max,x_max);

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
    Double_t phi = 0;
    Double_t modphi = 0;
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
        phi = hit->at(j).ph;
        if(phi<0) phi +=2.0*TMath::Pi();
        modphi = fmod(phi,2.0*TMath::Pi()/7.0);
//comment following line if want to plot full azimuth
//        if(modphi<3.0*TMath::Pi()/28.0 ||modphi>5.0*TMath::Pi()/28.0) continue;

//comment following line if want to plot all r
<<<<<<< HEAD
        if(hit->at(j).r>1000) continue;
=======
        if(hit->at(j).r>1500) continue;
>>>>>>> 0c4366e67ad16fe5832ddc8a375d8e73cd22dbf0
//comment following line if want to plot all E
//        if(hit->at(j).k<1) continue;

        h_rate[sp][dt][0]->Fill(hit->at(j).r,rate);
        h_rate[sp][dt][1]->Fill(hit->at(j).r,rate*hit->at(j).e);
        h_xy[sp][dt][0]->Fill(hit->at(j).x,hit->at(j).y,rate);
        h_xy[sp][dt][1]->Fill(hit->at(j).x,hit->at(j).y,rate*hit->at(j).e);

        if(hit->at(j).pz>=0){
          h_ratePzG0[sp][dt][0]->Fill(hit->at(j).r,rate);
          h_ratePzG0[sp][dt][1]->Fill(hit->at(j).r,rate*hit->at(j).e);
          h_xyPzG0[sp][dt][0]->Fill(hit->at(j).x,hit->at(j).y,rate);
          h_xyPzG0[sp][dt][1]->Fill(hit->at(j).x,hit->at(j).y,rate*hit->at(j).e);
        }else{
          h_ratePzL0[sp][dt][0]->Fill(hit->at(j).r,rate);
          h_ratePzL0[sp][dt][1]->Fill(hit->at(j).r,rate*hit->at(j).e);
          h_xyPzL0[sp][dt][0]->Fill(hit->at(j).x,hit->at(j).y,rate);
          h_xyPzL0[sp][dt][1]->Fill(hit->at(j).x,hit->at(j).y,rate*hit->at(j).e);
        }

        if(hit->at(j).pid==11 || hit->at(j).pid==-11){
          h_rate[4][dt][0]->Fill(hit->at(j).r,rate);
          h_rate[4][dt][1]->Fill(hit->at(j).r,rate*hit->at(j).e);
          h_xy[4][dt][0]->Fill(hit->at(j).x,hit->at(j).y,rate);
          h_xy[4][dt][1]->Fill(hit->at(j).x,hit->at(j).y,rate*hit->at(j).e);
          if(hit->at(j).pz>=0){
            h_ratePzG0[4][dt][0]->Fill(hit->at(j).r,rate);
            h_ratePzG0[4][dt][1]->Fill(hit->at(j).r,rate*hit->at(j).e);
            h_xyPzG0[4][dt][0]->Fill(hit->at(j).x,hit->at(j).y,rate);
            h_xyPzG0[4][dt][1]->Fill(hit->at(j).x,hit->at(j).y,rate*hit->at(j).e);
          }else{
            h_ratePzL0[4][dt][0]->Fill(hit->at(j).r,rate);
            h_ratePzL0[4][dt][1]->Fill(hit->at(j).r,rate*hit->at(j).e);
            h_xyPzL0[4][dt][0]->Fill(hit->at(j).x,hit->at(j).y,rate);
            h_xyPzL0[4][dt][1]->Fill(hit->at(j).x,hit->at(j).y,rate*hit->at(j).e);
          }
        }

        if(hit->at(j).trid==1){
          h_rate[5][dt][0]->Fill(hit->at(j).r,rate);
          h_rate[5][dt][1]->Fill(hit->at(j).r,rate*hit->at(j).e);
          h_xy[5][dt][0]->Fill(hit->at(j).x,hit->at(j).y,rate);
          h_xy[5][dt][1]->Fill(hit->at(j).x,hit->at(j).y,rate*hit->at(j).e);
          if(hit->at(j).pz>=0){
            h_ratePzG0[5][dt][0]->Fill(hit->at(j).r,rate);
            h_ratePzG0[5][dt][1]->Fill(hit->at(j).r,rate*hit->at(j).e);
            h_xyPzG0[5][dt][0]->Fill(hit->at(j).x,hit->at(j).y,rate);
            h_xyPzG0[5][dt][1]->Fill(hit->at(j).x,hit->at(j).y,rate*hit->at(j).e);
          }else{
            h_ratePzL0[5][dt][0]->Fill(hit->at(j).r,rate);
            h_ratePzL0[5][dt][1]->Fill(hit->at(j).r,rate*hit->at(j).e);
            h_xyPzG0[5][dt][0]->Fill(hit->at(j).x,hit->at(j).y,rate);
            h_xyPzG0[5][dt][1]->Fill(hit->at(j).x,hit->at(j).y,rate*hit->at(j).e);
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
         h_rate[iSp][iDet][iWt]->Scale(1.0/nTotEv);
         h_ratePzG0[iSp][iDet][iWt]->Scale(1.0/nTotEv);
         h_ratePzL0[iSp][iDet][iWt]->Scale(1.0/nTotEv);
         h_xy[iSp][iDet][iWt]->Scale(1.0/nTotEv);
         h_xyPzG0[iSp][iDet][iWt]->Scale(1.0/nTotEv);
         h_xyPzL0[iSp][iDet][iWt]->Scale(1.0/nTotEv);
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
