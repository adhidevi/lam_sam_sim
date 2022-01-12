#include "remolltypes.hh"
#include <sstream>
#include <iostream>
#include <fstream>
    TString rootfile_dir = "/volatile/halla/parity/adhidevi/remoll_rootfiles/PMTShielding";
    const string spTit[] = {"e-/#pi-","e+/#pi+","#gamma","neutron","e-/e+ E>1","primary E>1"};
    const int nSp = sizeof(spTit)/sizeof(*spTit);
    const string spH[nSp] ={"epiM","epiP","g","n","ee1","pri1"};
    map<int,int> spM {{11,1},{-211,1},{-11,2},{211,2},{22,3},{2112,4}};
    int color[nSp] = {2,3,1,4,6,7};
    const string detH[] = {"det186","det187","det188","det189","det286","det287","det288","det289","det386","det387","det388","det389","det486","det487","det488","det489","det586","det587","det588","det589","det686","det687","det688","det689","det786","det787","det788","det789","det886","det887","det888","det889"};
    map<int,int> dtM {{186,1},{187,2},{188,3},{189,4},{286,5},{287,6},{288,7},{289,8},{386,9},{387,10},{388,11},{389,12},{486,13},{487,14},{488,15},{489,16},{586,17},{587,18},{588,19},{589,20},{686,21},{687,22},{688,23},{689,24},{786,25},{787,26},{788,27},{789,28},{886,29},{887,30},{888,31},{889,32}};
    const int nDet = sizeof(detH)/sizeof(*detH);
    int Det[nDet] = {186,187,188,189,286,287,288,289,386,387,388,389,486,487,488,489,586,587,588,589,686,687,688,689,786,787,788,789,886,887,888,889};
    TH1D* eRate[nSp][nDet];
    TH1D* eRatePzG0[nSp][nDet];
    TH1D* eRatePzL0[nSp][nDet];
    void niceLogBins(TH1*);
    TFile* outfile;
    int beamGen(1);
    const string tgt_gen_config = "PMTSh_beam_V7";
void energy_SAMcathodePMT(){
   gStyle->SetOptStat(0);
   for(int iSp=0;iSp<nSp;iSp++){
   for(int iDet=0;iDet<nDet;iDet++){
      eRate[iSp][iDet] = new TH1D(Form("%s_E_%s",detH[iDet].c_str(),spH[iSp].c_str()),Form("Kinetic Energy on %s (%s);E (MeV);hits/%sthrownEvents",detH[iDet].c_str(),tgt_gen_config.c_str(),"#"),121,-8,5);
      eRatePzG0[iSp][iDet] = new TH1D(Form("%s_EPzG0_%s",detH[iDet].c_str(),spH[iSp].c_str()),Form("Kinetic Energy on %s (%s) pz>=0;E (MeV);hits/%sthrownEvents",detH[iDet].c_str(),tgt_gen_config.c_str(),"#"),121,-8,5);
      eRatePzL0[iSp][iDet] = new TH1D(Form("%s_EPzL0_%s",detH[iDet].c_str(),spH[iSp].c_str()),Form("Kinetic Energy on %s (%s) pz<0;E (MeV);hits/%sthrownEvents",detH[iDet].c_str(),tgt_gen_config.c_str(),"#"),121,-8,5);
      eRate[iSp][iDet]->SetLineColor(color[iSp]);
      niceLogBins(eRate[iSp][iDet]);
      niceLogBins(eRatePzG0[iSp][iDet]);
      niceLogBins(eRatePzL0[iSp][iDet]);
   }
   }
    outfile = new TFile(Form("./rootfiles/SAMcathodePMT_kinE_dist_%s.root",tgt_gen_config.c_str()),"recreate");
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

           eRate[sp][dt]->Fill(hit->at(j).k,rate);
           if(hit->at(j).pz>=0)
             eRatePzG0[sp][dt]->Fill(hit->at(j).k,rate);
           else
             eRatePzL0[sp][dt]->Fill(hit->at(j).k,rate);
           if(hit->at(j).k>1 && (hit->at(j).pid==11 || hit->at(j).pid==-11)){
             eRate[4][dt]->Fill(hit->at(j).k,rate);
             if(hit->at(j).pz>=0)
               eRatePzG0[4][dt]->Fill(hit->at(j).k,rate);
             else
               eRatePzL0[4][dt]->Fill(hit->at(j).k,rate);
           }
           if(hit->at(j).k>1 && hit->at(j).trid==1){
             eRate[5][dt]->Fill(hit->at(j).k,rate);
             if(hit->at(j).pz>=0)
               eRatePzG0[5][dt]->Fill(hit->at(j).k,rate);
             else
               eRatePzL0[5][dt]->Fill(hit->at(j).k,rate);
           }
         }
       }
       delete fin;
    }
   
    cout<<Form("Total number of file splits: %d",nfile)<<endl;
    cout<<Form("Total number of entries: %ld",nTotEv)<<endl;

    outfile->cd();
    for(int iDet=0;iDet<nDet;iDet++){
    outfile->mkdir(detH[iDet].c_str());
    outfile->cd(detH[iDet].c_str());
       for(int iSp=0;iSp<nSp;iSp++){
          eRate[iSp][iDet]->Scale(1./nTotEv);
          eRatePzG0[iSp][iDet]->Scale(1./nTotEv);
          eRatePzL0[iSp][iDet]->Scale(1./nTotEv);
          eRate[iSp][iDet]->Write();
          eRatePzG0[iSp][iDet]->Write();
          eRatePzL0[iSp][iDet]->Write();
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
