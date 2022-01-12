#include "remolltypes.hh"
#include <sstream>
#include <iostream>
#include <fstream>
<<<<<<< HEAD
    TString rootfile_dir = "/volatile/halla/parity/adhidevi/remoll_rootfiles/PMTShielding";
=======
    TString rootfile_dir = "/volatile/halla/parity/adhidevi/remoll_rootfiles/default-geo";
>>>>>>> 0c4366e67ad16fe5832ddc8a375d8e73cd22dbf0
    const string spTit[] = {"e-/#pi-","e+/#pi+","#gamma","neutron","e-/e+ E>1","primary E>1"};
    const int nSp = sizeof(spTit)/sizeof(*spTit);
    const string spH[nSp] ={"epiM","epiP","g","n","ee1","pri1"};
    map<int,int> spM {{11,1},{-211,1},{-11,2},{211,2},{22,3},{2112,4}};
    int color[nSp] = {2,3,1,4,6,7};
<<<<<<< HEAD
    const string detH[] = {"det28","det27"};
    map<int,int> dtM {{28,1},{27,2}};
    const int nDet = sizeof(detH)/sizeof(*detH);
    int Det[nDet] = {28,27};
    double rmin[nDet] = {1371.95,1371.95};
    double rmax[nDet] = {1651,1651};
=======
    const string detH[] = {"det28","det176"};
    map<int,int> dtM {{28,1},{176,2}};
    const int nDet = sizeof(detH)/sizeof(*detH);
    int Det[nDet] = {28,176};
    double rmin[nDet] = {200,200};
    double rmax[nDet] = {600,600};
>>>>>>> 0c4366e67ad16fe5832ddc8a375d8e73cd22dbf0
    TH1D* eRate[nSp][nDet];
    void niceLogBins(TH1*);
    TFile* outfile;
    int beamGen(1);
<<<<<<< HEAD
    const string tgt_gen_config = "PMTSh_beam_V10";
=======
    const string tgt_gen_config = "LH2_beam";
>>>>>>> 0c4366e67ad16fe5832ddc8a375d8e73cd22dbf0
void energy_distribution_smallR(){
   gStyle->SetOptStat(0);
   for(int iSp=0;iSp<nSp;iSp++){
   for(int iDet=0;iDet<nDet;iDet++){
      eRate[iSp][iDet] = new TH1D(Form("%s_E_%s",detH[iDet].c_str(),spH[iSp].c_str()),Form("Kinetic Energy on %s (%s), %.0f<=r<=%.0f;E (MeV);hits/%sthrownEvents",detH[iDet].c_str(),tgt_gen_config.c_str(),rmin[iDet],rmax[iDet],"#"),121,-8,5);
      eRate[iSp][iDet]->SetLineColor(color[iSp]);
      niceLogBins(eRate[iSp][iDet]);
   }
   }
<<<<<<< HEAD
    outfile = new TFile(Form("./rootfiles/DBM_kinE_dist_%s.root",tgt_gen_config.c_str()),"recreate");
    int nfile=0;
    Long64_t nentry=0;
    long nTotEv=0;
    for(int ifile=1001;ifile<=6000;ifile++){
=======
    outfile = new TFile(Form("./rootfiles/ring5_200to600mm_kinE_dist_%s.root",tgt_gen_config.c_str()),"recreate");
    int nfile=0;
    Long64_t nentry=0;
    long nTotEv=0;
    for(int ifile=1001;ifile<=1927;ifile++){
>>>>>>> 0c4366e67ad16fe5832ddc8a375d8e73cd22dbf0
       string infile = Form("%s/%s/%s_%d.root",rootfile_dir.Data(),tgt_gen_config.c_str(),tgt_gen_config.c_str(),ifile);
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
       std::vector<remollEventParticle_t> *part =0;
       remollBeamTarget_t *bm =0;
       remollEvent_t *ev =0;
       Double_t rate =0;
       Double_t phi = 0;
       Double_t modphi = 0;
       T->SetBranchAddress("hit", & hit);
       T->SetBranchAddress("part", & part);
       T->SetBranchAddress("bm", & bm);
       T->SetBranchAddress("ev", & ev);
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
//           if(modphi<3.0*TMath::Pi()/28.0 || modphi>5.0*TMath::Pi()/28.0) continue;

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
       delete fin;
    }
   
    cout<<Form("Total number of file splits: %d",nfile)<<endl;
    cout<<Form("Total number of entries: %ld",nTotEv)<<endl;

    outfile->cd();
    for(int iDet=0;iDet<nDet;iDet++){
       for(int iSp=0;iSp<nSp;iSp++){
          eRate[iSp][iDet]->Scale(1./nTotEv);
          eRate[iSp][iDet]->Write();
       }
    }
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
