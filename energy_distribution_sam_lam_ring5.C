#include "remolltypes.hh"
#include <sstream>
#include <iostream>
#include <fstream>

TString rootfile_dir = "/volatile/halla/parity/adhidevi/remoll_rootfiles/develop_br";
const string spTit[] = {"e-/#pi-","e+/#pi+","#gamma","neutron","e-/e+ (E>1MeV)","primary (KE>1MeV)"};
const int nSp = sizeof(spTit)/sizeof(*spTit);
const string spH[nSp] ={"epiM","epiP","g","n","ee1","pri1"};
map<int,int> spM {{11,1},{-211,1},{-11,2},{211,2},{22,3},{2112,4}};
int color[nSp] = {2,3,1,4,6,7};
const string detH[] = {"det174","det5719"};
map<int,int> dtM {{174,1},{5719,2}};
const int nDet = sizeof(detH)/sizeof(*detH);
int Det[nDet] = {174,5719};
double rmin[nDet] = {1010,1010};
double rmax[nDet] = {1130,1130};
TH1D* eRate[nSp][nDet];
TH1D* eRate_pzG0[nSp][nDet];
TH1D* eRate_pzL0[nSp][nDet];
void niceLogBins(TH1*);
TFile* outfile;
int beamGen(1);
const string tgt_gen_config = "LH2_beam_V12";
string hist_title;

void energy_distribution_sam_lam_ring5(){
   gStyle->SetOptStat(0);
   for(int iSp=0;iSp<nSp;iSp++){
    for(int iDet=0;iDet<nDet;iDet++){
      if(beamGen)
      hist_title = Form("Kinetic Energy on %s (%s), %.0fmm<=r<=%.0fmm;E (MeV);hits/%sthrownEvents",detH[iDet].c_str(),tgt_gen_config.c_str(),rmin[iDet],rmax[iDet],"#");
      else
      hist_title = Form("Kinetic Energy on %s (%s), %.0fmm<=r<=%.0fmm;E (MeV);rate (GHz)",detH[iDet].c_str(),tgt_gen_config.c_str(),rmin[iDet],rmax[iDet]);

      eRate[iSp][iDet] = new TH1D(Form("%s_E_%s",detH[iDet].c_str(),spH[iSp].c_str()),hist_title.c_str(),121,-8,5);
      eRate_pzG0[iSp][iDet] = new TH1D(Form("%s_EpzG0_%s",detH[iDet].c_str(),spH[iSp].c_str()),hist_title.c_str(),121,-8,5);
      eRate_pzL0[iSp][iDet] = new TH1D(Form("%s_EpzL0_%s",detH[iDet].c_str(),spH[iSp].c_str()),hist_title.c_str(),121,-8,5);

      eRate[iSp][iDet]->SetLineColor(color[iSp]);
      eRate_pzG0[iSp][iDet]->SetLineColor(color[iSp]);
      eRate_pzL0[iSp][iDet]->SetLineColor(color[iSp]);
      niceLogBins(eRate[iSp][iDet]);
      niceLogBins(eRate_pzG0[iSp][iDet]);
      niceLogBins(eRate_pzL0[iSp][iDet]);
    }
   }
    outfile = TFile::Open(Form("./rootfiles/lam_kinE_%s.root",tgt_gen_config.c_str()),"RECREATE");
    int nfile=0;
    Long64_t nentry=0;
    long nTotEv=0;
    for(int ifile=1001;ifile<=6000;ifile++){
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
           
           if(hit->at(j).r>=rmin[dt] && hit->at(j).r<=rmax[dt]){
             eRate[sp][dt]->Fill(hit->at(j).k,rate);
              if(hit->at(j).pz>=0){
               eRate_pzG0[sp][dt]->Fill(hit->at(j).k,rate);
              }else{
               eRate_pzL0[sp][dt]->Fill(hit->at(j).k,rate);
              }
             if(hit->at(j).k>1 && (hit->at(j).pid==11 || hit->at(j).pid==-11)){
               eRate[4][dt]->Fill(hit->at(j).k,rate);
              if(hit->at(j).pz>=0){
               eRate_pzG0[4][dt]->Fill(hit->at(j).k,rate);
              }else{
               eRate_pzL0[4][dt]->Fill(hit->at(j).k,rate);
              }
             }
             if(hit->at(j).k>1 && hit->at(j).vz<=-3875){
               eRate[5][dt]->Fill(hit->at(j).k,rate);
               if(hit->at(j).pz>=0){
                eRate_pzG0[5][dt]->Fill(hit->at(j).k,rate);
               }else{
                eRate_pzL0[5][dt]->Fill(hit->at(j).k,rate);
               }
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
        if(beamGen){
          eRate[iSp][iDet]->Scale(1.0/nTotEv);
          eRate_pzG0[iSp][iDet]->Scale(1.0/nTotEv);
          eRate_pzL0[iSp][iDet]->Scale(1.0/nTotEv);
        }else{
          eRate[iSp][iDet]->Scale(1.0e-9/nfile);
          eRate_pzG0[iSp][iDet]->Scale(1.0e-9/nfile);
          eRate_pzL0[iSp][iDet]->Scale(1.0e-9/nfile);
        }
        eRate[iSp][iDet]->Write();
        eRate_pzG0[iSp][iDet]->Write();
        eRate_pzL0[iSp][iDet]->Write();
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
