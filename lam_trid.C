#include "remolltypes.hh"
#include <sstream>
#include <iostream>
#include <fstream>
void isValid1(std::vector<remollGenericDetectorHit_t> *fHit, std::vector<int> &flange_trid, std::vector<int> &collar2OR_trid);
const double pi = TMath::Pi();
const double lam_length = 360.0;//azimuthal length of LAM quartz
const double lam_rin = 1010.0;//inner radius of LAM quartz
double lam_angle = atan(lam_length/lam_rin);
double sep_mid = 2*pi/14.0;

void lam_trid(){
  gROOT->Reset();
  gStyle->SetOptStat(0);
  gStyle->SetTitleYOffset(1.3);
  gStyle->SetPadGridX(1);
  gStyle->SetPadGridY(1);
  TGaxis::SetMaxDigits(3);

  const string spTit[] = {"e-/#pi- (KE>1 MeV)","e+/#pi+ (KE>1 MeV)","#gamma (KE>1 MeV)","neutron (KE>1 MeV)","e-/e+ (KE>1 MeV)","e-/#pi- vz<=-3875 all KE","e- trid=1 (KE>1 MeV)"};
//  const string spTit[] = {"e-/#pi- all KE","e+/#pi+ all KE","#gamma all KE","neutron all KE","e-/e+ all KE","e-/#pi- vz<=-3875 all KE","e- trid=1 all KE"};
  const int nSp = sizeof(spTit)/sizeof(*spTit);
  const string spH[nSp] = {"epiM","epiP","g","n","ee","epiMvzCut","eTrIdCut"};
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
  const string tridCut[] = {"flange_trid","collar2OR_trid","rate"};
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
  const string tgt_gen_config = "LH2_beam_V18";
  const string plotType = "lam_trid_vzCut_EG1_split2";//rCut or rNoCut and EG1 or allE
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

      std::vector<int> flange_trid;
      std::vector<int> collar2OR_trid;
      std::vector<int>::iterator flange_it;
      std::vector<int>::iterator collar2OR_it;
      isValid1(hit,flange_trid,collar2OR_trid);

      for(int j=0;j<hit->size();j++){
        if(std::isnan(rate) || std::isinf(rate)) continue;
        if(beamGen) rate = 1.0;

        flange_it = find(flange_trid.begin(),flange_trid.end(),hit->at(j).trid);
        collar2OR_it = find(collar2OR_trid.begin(),collar2OR_trid.end(),hit->at(j).trid);

        int sp = spM[int(hit->at(j).pid)]-1;
        if(sp==-1) continue;
        int dt = dtM[int(hit->at(j).det)]-1;
        if(dt==-1) continue;
        if(hit->at(j).k<1) continue;

        if(flange_it!=flange_trid.end()){
         if(hit->at(j).trid==*flange_it){
           double theta = (180.0/pi)*atan2(sqrt(pow(hit->at(j).px,2)+pow(hit->at(j).py,2)),hit->at(j).pz);
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
        if(collar2OR_it!=collar2OR_trid.end()){
         if(hit->at(j).trid==*collar2OR_it){
           double theta = (180.0/pi)*atan2(sqrt(pow(hit->at(j).px,2)+pow(hit->at(j).py,2)),hit->at(j).pz);
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

void isValid1(std::vector<remollGenericDetectorHit_t> *hit, std::vector<int> &flange_trid, std::vector<int> &collar2OR_trid){
  for(size_t j=0;j<hit->size();j++){
//    double phi = hit->at(j).ph;
//    if(phi<0) phi +=2.0*pi;
//    double modphi = fmod(phi,2.0*pi/7.0);
//    if(modphi<3.0*pi/28.0 || modphi>5.0*pi/28.0) continue;
//    for(int sep=1;sep<=14;sep++){
//      if(sep%2==1 && (abs(phi)>=sep*sep_mid-lam_angle/2.0 && abs(phi)<=sep*sep_mid+lam_angle/2.0)){
        if(hit->at(j).det==3174&&hit->at(j).vz>18850&&hit->at(j).vz<18925){
          flange_trid.push_back(hit->at(j).trid);
        }
        if(hit->at(j).det==2174&&hit->at(j).vz>18925&&hit->at(j).vz<19090){
          collar2OR_trid.push_back(hit->at(j).trid);
        }
//      }
//    }
  }
}
