#include "remolltypes.hh"
#include <sstream>
#include <iostream>
#include <fstream>
void isValid1(std::vector<remollGenericDetectorHit_t> *fHit, std::vector<int> &xMotor_trid, std::vector<int> &yMotor_trid);

void usScanner_motor_trid(){
  gROOT->Reset();
  gStyle->SetOptStat(0);
  gStyle->SetTitleYOffset(1.3);
  gStyle->SetPadGridX(1);
  gStyle->SetPadGridY(1);
  TGaxis::SetMaxDigits(3);

  const double kinEcut = 0;//MeV
  const string spTit[] = {"e-/#pi-","e+/#pi+","#gamma","neutron","e-/e+","e-/#pi- vz<=-3875","e- trid=1"};
  const int nSp = sizeof(spTit)/sizeof(*spTit);
  const string spH[nSp] = {"epiM","epiP","g","n","ee","epiMvzCut","eTrIdCut"};
  map<int,int> spM {{11,1},{-211,1},{-11,2},{211,2},{22,3},{2112,4}};

///Change the following lines for which detectors you want to include////
  string detH[] = {"det1175","det2175","det28"};
  const int nDet = sizeof(detH)/sizeof(*detH);
  const int Det[nDet] = {1175,2175,28};
  map<int,int> dtM {{1175,1},{2175,2},{28,3}};
////////////////////////////////////////////////////////////////////////

  double rmin=0.0;
  double rmax=3000.0;
  int nbins=1000;
  const string tridCut[] = {"xMotor_trid","yMotor_trid","rate"};
  const int nCut = sizeof(tridCut)/sizeof(*tridCut);
  const string tridCut_unit[nCut] ={"hits/#thrownEvents","hits/#thrownEvents","hits/#thrownEvents"};

  TH1F* h_d_r[nSp][nDet][nCut];
  TH1F* h_d_rPzG0[nSp][nDet][nCut];
  TH1F* h_d_rPzL0[nSp][nDet][nCut];
  TH2F* h_xy[nSp][nDet][nCut];
  TH2F* h_xyPzG0[nSp][nDet][nCut];
  TH2F* h_xyPzL0[nSp][nDet][nCut];
///Change the following lines as needed////
  const string geometry = "develop";
  const string tgt_gen_config = "LH2_beam_V19";
  const string plotType = Form("usScanner_motor_trid_NovzCut_EG%.0f",kinEcut);//rCut or rNoCut and EG1 or allE
  int beamGen(1);
///Change this line for appropriate rootfile directory////
  TString rootfile_dir = "/volatile/halla/parity/adhidevi/remoll_rootfiles/develop_br";
//////////////////////////////////////////////////////////

  for(int iSp=0;iSp<nSp;iSp++){
   for(int iDet=0;iDet<nDet;iDet++){
    for(int iCut=0;iCut<nCut;iCut++){
      string titleR = Form("%s Radial dist. on %s plane (%s);Radius (mm);%s/%.1fmm",spTit[iSp].c_str(),detH[iDet].c_str(),tgt_gen_config.c_str(),tridCut_unit[iCut].c_str(),(rmax-rmin)/nbins);
      h_d_r[iSp][iDet][iCut] = new TH1F(Form("%s_r_%s_%s",detH[iDet].c_str(),spH[iSp].c_str(),tridCut[iCut].c_str()),titleR.c_str(),nbins,rmin,rmax);
      h_d_rPzG0[iSp][iDet][iCut] = new TH1F(Form("%s_rPzG0_%s_%s",detH[iDet].c_str(),spH[iSp].c_str(),tridCut[iCut].c_str()),titleR.c_str(),nbins,rmin,rmax);
      h_d_rPzL0[iSp][iDet][iCut] = new TH1F(Form("%s_rPzL0_%s_%s",detH[iDet].c_str(),spH[iSp].c_str(),tridCut[iCut].c_str()),titleR.c_str(),nbins,rmin,rmax);

      string titleXY = Form("%s XY dist. on %s plane (%s);x (mm);y (mm)",spTit[iSp].c_str(),detH[iDet].c_str(),tgt_gen_config.c_str());
      h_xy[iSp][iDet][iCut] = new TH2F(Form("%s_xy_%s_%s",detH[iDet].c_str(),spH[iSp].c_str(),tridCut[iCut].c_str()),titleXY.c_str(),nbins,-rmax,rmax,nbins,-rmax,rmax);
      h_xyPzG0[iSp][iDet][iCut] = new TH2F(Form("%s_xyPzG0_%s_%s",detH[iDet].c_str(),spH[iSp].c_str(),tridCut[iCut].c_str()),titleXY.c_str(),nbins,-rmax,rmax,nbins,-rmax,rmax);
      h_xyPzL0[iSp][iDet][iCut] = new TH2F(Form("%s_xyPzL0_%s_%s",detH[iDet].c_str(),spH[iSp].c_str(),tridCut[iCut].c_str()),titleXY.c_str(),nbins,-rmax,rmax,nbins,-rmax,rmax);
      h_d_r[iSp][iDet][iCut]->Sumw2();
      h_d_rPzG0[iSp][iDet][iCut]->Sumw2();
      h_d_rPzL0[iSp][iDet][iCut]->Sumw2();
    }
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

      std::vector<int> xMotor_trid;
      std::vector<int> yMotor_trid;
      std::vector<int>::iterator xMotor_it;
      std::vector<int>::iterator yMotor_it;
      isValid1(hit,xMotor_trid,yMotor_trid);

      for(int j=0;j<hit->size();j++){
        if(std::isnan(rate) || std::isinf(rate)) continue;
        if(beamGen) rate = 1.0;

        xMotor_it = find(xMotor_trid.begin(),xMotor_trid.end(),hit->at(j).trid);
        yMotor_it = find(yMotor_trid.begin(),yMotor_trid.end(),hit->at(j).trid);

        int sp = spM[int(hit->at(j).pid)]-1;
        if(sp==-1) continue;
        int dt = dtM[int(hit->at(j).det)]-1;
        if(dt==-1) continue;
        if(hit->at(j).k<kinEcut) continue;

        if(xMotor_it!=xMotor_trid.end()){
         if(hit->at(j).trid==*xMotor_it){
           h_d_r[sp][dt][0]->Fill(hit->at(j).r,rate);
           h_xy[sp][dt][0]->Fill(hit->at(j).x,hit->at(j).y,rate);

           if(hit->at(j).pz>=0){
             h_d_rPzG0[sp][dt][0]->Fill(hit->at(j).r,rate);
             h_xyPzG0[sp][dt][0]->Fill(hit->at(j).x,hit->at(j).y,rate);
           }else{
             h_d_rPzL0[sp][dt][0]->Fill(hit->at(j).r,rate);
             h_xyPzL0[sp][dt][0]->Fill(hit->at(j).x,hit->at(j).y,rate);
           }

           if(hit->at(j).pid==11 || hit->at(j).pid==-11){
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

           if(hit->at(j).vz<=-3875 && (hit->at(j).pid==11 || hit->at(j).pid==-211)){
             h_d_r[5][dt][0]->Fill(hit->at(j).r,rate);
             h_xy[5][dt][0]->Fill(hit->at(j).x,hit->at(j).y,rate);
             if(hit->at(j).pz>=0){
               h_d_rPzG0[5][dt][0]->Fill(hit->at(j).r,rate);
               h_xyPzG0[5][dt][0]->Fill(hit->at(j).x,hit->at(j).y,rate);
             }else{
               h_d_rPzL0[5][dt][0]->Fill(hit->at(j).r,rate);
               h_xyPzG0[5][dt][0]->Fill(hit->at(j).x,hit->at(j).y,rate);
             }
           }

           if(hit->at(j).trid==1 && hit->at(j).pid==11){
             h_d_r[6][dt][0]->Fill(hit->at(j).r,rate);
             h_xy[6][dt][0]->Fill(hit->at(j).x,hit->at(j).y,rate);
             if(hit->at(j).pz>=0){
               h_d_rPzG0[6][dt][0]->Fill(hit->at(j).r,rate);
               h_xyPzG0[6][dt][0]->Fill(hit->at(j).x,hit->at(j).y,rate);
             }else{
               h_d_rPzL0[6][dt][0]->Fill(hit->at(j).r,rate);
               h_xyPzG0[6][dt][0]->Fill(hit->at(j).x,hit->at(j).y,rate);
             }
           }
         }
        }
        if(yMotor_it!=yMotor_trid.end()){
         if(hit->at(j).trid==*yMotor_it){
           h_d_r[sp][dt][1]->Fill(hit->at(j).r,rate);
           h_xy[sp][dt][1]->Fill(hit->at(j).x,hit->at(j).y,rate);

           if(hit->at(j).pz>=0){
             h_d_rPzG0[sp][dt][1]->Fill(hit->at(j).r,rate);
             h_xyPzG0[sp][dt][1]->Fill(hit->at(j).x,hit->at(j).y,rate);
           }else{
             h_d_rPzL0[sp][dt][1]->Fill(hit->at(j).r,rate);
             h_xyPzL0[sp][dt][1]->Fill(hit->at(j).x,hit->at(j).y,rate);
           }

           if(hit->at(j).pid==11 || hit->at(j).pid==-11){
             h_d_r[4][dt][1]->Fill(hit->at(j).r,rate);
             h_xy[4][dt][1]->Fill(hit->at(j).x,hit->at(j).y,rate);
             if(hit->at(j).pz>=0){
               h_d_rPzG0[4][dt][1]->Fill(hit->at(j).r,rate);
               h_xyPzG0[4][dt][1]->Fill(hit->at(j).x,hit->at(j).y,rate);
             }else{
               h_d_rPzL0[4][dt][1]->Fill(hit->at(j).r,rate);
               h_xyPzL0[4][dt][1]->Fill(hit->at(j).x,hit->at(j).y,rate);
             }
           }

           if(hit->at(j).vz<=-3875 && (hit->at(j).pid==11 || hit->at(j).pid==-211)){
             h_d_r[5][dt][1]->Fill(hit->at(j).r,rate);
             h_xy[5][dt][1]->Fill(hit->at(j).x,hit->at(j).y,rate);
             if(hit->at(j).pz>=0){
               h_d_rPzG0[5][dt][1]->Fill(hit->at(j).r,rate);
               h_xyPzG0[5][dt][1]->Fill(hit->at(j).x,hit->at(j).y,rate);
             }else{
               h_d_rPzL0[5][dt][1]->Fill(hit->at(j).r,rate);
               h_xyPzG0[5][dt][1]->Fill(hit->at(j).x,hit->at(j).y,rate);
             }
           }

           if(hit->at(j).trid==1 && hit->at(j).pid==11){
             h_d_r[6][dt][1]->Fill(hit->at(j).r,rate);
             h_xy[6][dt][1]->Fill(hit->at(j).x,hit->at(j).y,rate);
             if(hit->at(j).pz>=0){
               h_d_rPzG0[6][dt][1]->Fill(hit->at(j).r,rate);
               h_xyPzG0[6][dt][1]->Fill(hit->at(j).x,hit->at(j).y,rate);
             }else{
               h_d_rPzL0[6][dt][1]->Fill(hit->at(j).r,rate);
               h_xyPzG0[6][dt][1]->Fill(hit->at(j).x,hit->at(j).y,rate);
             }
           }
         }
        }
        h_d_r[sp][dt][2]->Fill(hit->at(j).r,rate);
        h_xy[sp][dt][2]->Fill(hit->at(j).x,hit->at(j).y,rate);
        if(hit->at(j).pz>=0){
          h_d_rPzG0[sp][dt][2]->Fill(hit->at(j).r,rate);
          h_xyPzG0[sp][dt][2]->Fill(hit->at(j).x,hit->at(j).y,rate);
        }else{
          h_d_rPzL0[sp][dt][2]->Fill(hit->at(j).r,rate);
          h_xyPzL0[sp][dt][2]->Fill(hit->at(j).x,hit->at(j).y,rate);
        }

        if(hit->at(j).pid==11 || hit->at(j).pid==-11){
          h_d_r[4][dt][2]->Fill(hit->at(j).r,rate);
          h_xy[4][dt][2]->Fill(hit->at(j).x,hit->at(j).y,rate);
          if(hit->at(j).pz>=0){
            h_d_rPzG0[4][dt][2]->Fill(hit->at(j).r,rate);
            h_xyPzG0[4][dt][2]->Fill(hit->at(j).x,hit->at(j).y,rate);
          }else{
            h_d_rPzL0[4][dt][2]->Fill(hit->at(j).r,rate);
            h_xyPzL0[4][dt][2]->Fill(hit->at(j).x,hit->at(j).y,rate);
          }
        }
        if(hit->at(j).vz<=-3875 && (hit->at(j).pid==11 || hit->at(j).pid==-211)){
          h_d_r[5][dt][2]->Fill(hit->at(j).r,rate);
          h_xy[5][dt][2]->Fill(hit->at(j).x,hit->at(j).y,rate);
          if(hit->at(j).pz>=0){
            h_d_rPzG0[5][dt][2]->Fill(hit->at(j).r,rate);
            h_xyPzG0[5][dt][2]->Fill(hit->at(j).x,hit->at(j).y,rate);
          }else{
            h_d_rPzL0[5][dt][2]->Fill(hit->at(j).r,rate);
            h_xyPzG0[5][dt][2]->Fill(hit->at(j).x,hit->at(j).y,rate);
          }
        }
        if(hit->at(j).trid==1 && hit->at(j).pid==11){
          h_d_r[6][dt][2]->Fill(hit->at(j).r,rate);
          h_xy[6][dt][2]->Fill(hit->at(j).x,hit->at(j).y,rate);
          if(hit->at(j).pz>=0){
            h_d_rPzG0[6][dt][2]->Fill(hit->at(j).r,rate);
            h_xyPzG0[6][dt][2]->Fill(hit->at(j).x,hit->at(j).y,rate);
          }else{
            h_d_rPzL0[6][dt][2]->Fill(hit->at(j).r,rate);
            h_xyPzG0[6][dt][2]->Fill(hit->at(j).x,hit->at(j).y,rate);
          }
        }
      }
    }
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

void isValid1(std::vector<remollGenericDetectorHit_t> *hit, std::vector<int> &xMotor_trid, std::vector<int> &yMotor_trid){
  for(size_t j=0;j<hit->size();j++){
     if(hit->at(j).det==1175/*&&hit->at(j).vz>=21356.51&&hit->at(j).vz<=21356.51+85.0*/){
       xMotor_trid.push_back(hit->at(j).trid);
     }
     if(hit->at(j).det==2175/*&&hit->at(j).vz>=21305.66&&hit->at(j).vz<=21305.66+56.39*/){
       yMotor_trid.push_back(hit->at(j).trid);
     }
  }
}
