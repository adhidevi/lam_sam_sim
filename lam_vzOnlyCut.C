#include "remolltypes.hh"
#include <sstream>
#include <iostream>
#include <fstream>
const double pi = TMath::Pi();
const double lam_length = 360.0;//azimuthal length of LAM quartz
const double lam_rin = 1010.0;//inner radius of LAM quartz
double lam_angle = atan(lam_length/lam_rin);
double sep_mid = 2*pi/14.0;

void lam_vzOnlyCut(){
  gROOT->Reset();
  gStyle->SetOptStat(0);
  gStyle->SetTitleYOffset(1.3);
  gStyle->SetPadGridX(1);
  gStyle->SetPadGridY(1);
  TGaxis::SetMaxDigits(3);

  const double kinEcut = 1;//MeV
  const string spTit[] = {"e-/#pi-","e+/#pi+","#gamma","neutron","e- trid=1"};
  const int nSp = sizeof(spTit)/sizeof(*spTit);
  const string spH[nSp] = {"epiM","epiP","g","n","eTrIdCut"};
  map<int,int> spM {{11,1},{-211,1},{-11,2},{211,2},{22,3},{2112,4}};

///Change the following lines for which detectors you want to include////
  string detH[] = {"det2174","det3174","det174"};
  const int nDet = sizeof(detH)/sizeof(*detH);
  const int Det[nDet] = {2174,3174,174};
  map<int,int> dtM {{2174,1},{3174,2},{174,3}};
////////////////////////////////////////////////////////////////////////

  double rmin=0.0;
  double rmax=1900.0;
  int nbins=500;
  const string vzCut[] = {"no_vzCut","target_vzCut","flange_vzCut","collar2OR_vzCut"};
  const int nCut = sizeof(vzCut)/sizeof(*vzCut);
  const string vzCut_unit[nCut] ={"hits/#thrownEvents","hits/#thrownEvents","hits/#thrownEvents","hits/#thrownEvents"};

  TH1F* h_d_r[nSp][nDet][nCut];
  TH1F* h_d_rPzG0[nSp][nDet][nCut];
  TH1F* h_d_rPzL0[nSp][nDet][nCut];
  TH2F* h_xy[nSp][nDet][nCut];
  TH2F* h_xyPzG0[nSp][nDet][nCut];
  TH2F* h_xyPzL0[nSp][nDet][nCut];
///Change the following lines as needed////
  const string geometry = "develop";
  const string tgt_gen_config = "LH2_beam_V18";
  const string plotType = Form("lam_vzOnlyCut_EG%.0f_split2",kinEcut);//rCut or rNoCut and EG1 or allE
  int beamGen(1);
///Change this line for appropriate rootfile directory////
  TString rootfile_dir = "/volatile/halla/parity/adhidevi/remoll_rootfiles/develop_br";
//////////////////////////////////////////////////////////

  for(int iSp=0;iSp<nSp;iSp++){
   for(int iDet=0;iDet<nDet;iDet++){
    for(int iCut=0;iCut<nCut;iCut++){
      string titleR = Form("%s (KE>%.0f MeV) Radial dist. on %s plane (%s);Radius (mm);%s/%.1fmm",spTit[iSp].c_str(),kinEcut,detH[iDet].c_str(),tgt_gen_config.c_str(),vzCut_unit[iCut].c_str(),(rmax-rmin)/nbins);
      h_d_r[iSp][iDet][iCut] = new TH1F(Form("%s_r_%s_%s",detH[iDet].c_str(),spH[iSp].c_str(),vzCut[iCut].c_str()),titleR.c_str(),nbins,rmin,rmax);
      h_d_rPzG0[iSp][iDet][iCut] = new TH1F(Form("%s_rPzG0_%s_%s",detH[iDet].c_str(),spH[iSp].c_str(),vzCut[iCut].c_str()),titleR.c_str(),nbins,rmin,rmax);
      h_d_rPzL0[iSp][iDet][iCut] = new TH1F(Form("%s_rPzL0_%s_%s",detH[iDet].c_str(),spH[iSp].c_str(),vzCut[iCut].c_str()),titleR.c_str(),nbins,rmin,rmax);

      string titleXY = Form("%s (KE>%.0f MeV) XY dist. on %s plane (%s);x (mm);y (mm)",spTit[iSp].c_str(),kinEcut,detH[iDet].c_str(),tgt_gen_config.c_str());
      h_xy[iSp][iDet][iCut] = new TH2F(Form("%s_xy_%s_%s",detH[iDet].c_str(),spH[iSp].c_str(),vzCut[iCut].c_str()),titleXY.c_str(),nbins,-rmax,rmax,nbins,-rmax,rmax);
      h_xyPzG0[iSp][iDet][iCut] = new TH2F(Form("%s_xyPzG0_%s_%s",detH[iDet].c_str(),spH[iSp].c_str(),vzCut[iCut].c_str()),titleXY.c_str(),nbins,-rmax,rmax,nbins,-rmax,rmax);
      h_xyPzL0[iSp][iDet][iCut] = new TH2F(Form("%s_xyPzL0_%s_%s",detH[iDet].c_str(),spH[iSp].c_str(),vzCut[iCut].c_str()),titleXY.c_str(),nbins,-rmax,rmax,nbins,-rmax,rmax);
      h_d_r[iSp][iDet][iCut]->Sumw2();
      h_d_rPzG0[iSp][iDet][iCut]->Sumw2();
      h_d_rPzL0[iSp][iDet][iCut]->Sumw2();
    }
   }
  }

  int nfile=0;
  Long64_t nentry=0;
  long nTotEv=0;
  for(int ifile=3501;ifile<=6000;ifile++){
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
        if(hit->at(j).k<kinEcut) continue;

        h_d_r[sp][dt][0]->Fill(hit->at(j).r,rate);
        h_xy[sp][dt][0]->Fill(hit->at(j).x,hit->at(j).y,rate);
        if(hit->at(j).pz>=0){
          h_d_rPzG0[sp][dt][0]->Fill(hit->at(j).r,rate);
          h_xyPzG0[sp][dt][0]->Fill(hit->at(j).x,hit->at(j).y,rate);
        }else{
          h_d_rPzL0[sp][dt][0]->Fill(hit->at(j).r,rate);
          h_xyPzL0[sp][dt][0]->Fill(hit->at(j).x,hit->at(j).y,rate);
        }

        if(hit->at(j).trid==1 && hit->at(j).pid==11){
        h_d_r[4][dt][0]->Fill(hit->at(j).r,rate);
        h_xy[4][dt][0]->Fill(hit->at(j).x,hit->at(j).y,rate);
        if(hit->at(j).pz>=0){
          h_d_rPzG0[4][dt][0]->Fill(hit->at(j).r,rate);
          h_xyPzG0[4][dt][0]->Fill(hit->at(j).x,hit->at(j).y,rate);
        }else{
          h_d_rPzL0[4][dt][0]->Fill(hit->at(j).r,rate);
          h_xyPzL0[4][dt][0]->Fill(hit->at(j).x,hit->at(j).y,rate);
        }

        }
 
        if(hit->at(j).vz<=-3875){ 
          h_d_r[sp][dt][1]->Fill(hit->at(j).r,rate);
          h_xy[sp][dt][1]->Fill(hit->at(j).x,hit->at(j).y,rate);
          if(hit->at(j).pz>=0){
            h_d_rPzG0[sp][dt][1]->Fill(hit->at(j).r,rate);
            h_xyPzG0[sp][dt][1]->Fill(hit->at(j).x,hit->at(j).y,rate);
          }else{
            h_d_rPzL0[sp][dt][1]->Fill(hit->at(j).r,rate);
            h_xyPzL0[sp][dt][1]->Fill(hit->at(j).x,hit->at(j).y,rate);
          }
        }

        if(hit->at(j).vz>18850&&hit->at(j).vz<18925){ 
          h_d_r[sp][dt][2]->Fill(hit->at(j).r,rate);
          h_xy[sp][dt][2]->Fill(hit->at(j).x,hit->at(j).y,rate);
          if(hit->at(j).pz>=0){
            h_d_rPzG0[sp][dt][2]->Fill(hit->at(j).r,rate);
            h_xyPzG0[sp][dt][2]->Fill(hit->at(j).x,hit->at(j).y,rate);
          }else{
            h_d_rPzL0[sp][dt][2]->Fill(hit->at(j).r,rate);
            h_xyPzL0[sp][dt][2]->Fill(hit->at(j).x,hit->at(j).y,rate);
          }
        }

        if(hit->at(j).vz>18925&&hit->at(j).vz<19090){ 
          h_d_r[sp][dt][3]->Fill(hit->at(j).r,rate);
          h_xy[sp][dt][3]->Fill(hit->at(j).x,hit->at(j).y,rate);
          if(hit->at(j).pz>=0){
            h_d_rPzG0[sp][dt][3]->Fill(hit->at(j).r,rate);
            h_xyPzG0[sp][dt][3]->Fill(hit->at(j).x,hit->at(j).y,rate);
          }else{
            h_d_rPzL0[sp][dt][3]->Fill(hit->at(j).r,rate);
            h_xyPzL0[sp][dt][3]->Fill(hit->at(j).x,hit->at(j).y,rate);
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
  TFile* outfile = new TFile(Form("./rootfiles/%s_%s_%s_NoPhiCut.root",geometry.c_str(),tgt_gen_config.c_str(),plotType.c_str()),"recreate");
//////////////////////////////////////////   
  for(int iDet=0;iDet<nDet;iDet++){
    outfile->mkdir(Form("%s",detH[iDet].c_str()));
    outfile->cd(Form("%s",detH[iDet].c_str()));
    for(int iSp=0;iSp<nSp;iSp++){
      for(int iCut=0;iCut<nCut;iCut++){
         h_d_r[iSp][iDet][iCut]->Scale(1.0/nTotEv);
         h_d_rPzG0[iSp][iDet][iCut]->Scale(1.0/nTotEv);
         h_d_rPzL0[iSp][iDet][iCut]->Scale(1.0/nTotEv);
         h_xy[iSp][iDet][iCut]->Scale(1.0/nTotEv);
         h_xyPzG0[iSp][iDet][iCut]->Scale(1.0/nTotEv);
         h_xyPzL0[iSp][iDet][iCut]->Scale(1.0/nTotEv);

         h_d_r[iSp][iDet][iCut]->Write("",TObject::kOverwrite);
         h_d_rPzG0[iSp][iDet][iCut]->Write("",TObject::kOverwrite);
         h_d_rPzL0[iSp][iDet][iCut]->Write("",TObject::kOverwrite);
         h_xy[iSp][iDet][iCut]->Write("",TObject::kOverwrite);
         h_xyPzG0[iSp][iDet][iCut]->Write("",TObject::kOverwrite);
         h_xyPzL0[iSp][iDet][iCut]->Write("",TObject::kOverwrite);
      }
    }
  }
}
