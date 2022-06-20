//void skimTree(string finNm, int testRun=0, int beamGen=0)
//if testRun=1 it only does 10 files max, if beamGen=1 it rescales things...
//manually type in detector number you want in macro as mydet below (i.e. mydet=28)
//manually adjust cut for what you want to find in processOne below (i.e hit.trid==1)


// Run commands in this order:
// 
// //Load in the script, and run it
// reroot 'skimTree("/volatile/halla/moller12gev/palatchi/remoll_ferrous_parallel_100M_20211211/remoll_ferrous_parallel_100M_20211211_100kEv_001/remollout.root",0,1)'
//
// OR point to a textfile list of root files with file names  
// reroot 'skimTree("mylist.txt",0,1)'
//

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include <vector>
using namespace std;

int mydet=1751;
TFile *ofile;
TTree *otree;
TBranch *b_rate,*b_hit;
Double_t newrate=0;
vector<remollGenericDetectorHit_t>  *newhit=0;
TString outputdir = "/volatile/halla/parity/adhidevi/remoll_rootfiles/develop_br/LH2_beam_V23/skimTree";
double scaleRate(1);

long processOne(string fnm);
long getEvents(string);

void skimTree(string finNm, int testRun=0, int beamGen=1){
  
  long nTotHits(0);
  int nFiles(0);

  if(finNm==""){
    cout<<"\t did not find input file. Quitting!"<<endl;
    return;
  }

  ofile = new TFile(Form("%s/o_remollSkimTree%d.root",outputdir.Data(),mydet),"RECREATE");

  otree     = new TTree("T", "skim tree");
  b_rate = otree->Branch("rate", &newrate);
  b_hit   = otree->Branch("hit", &newhit);

  if( finNm.find(".root") < finNm.size() ){
    cout<<"Processing single file:\n\t"<<finNm<<endl;
    if(beamGen){
      scaleRate = getEvents(finNm);
    }

    nTotHits+=processOne(finNm);
    nFiles=1;
  }else{
    cout<<"Attempting to process list of output from\n\t"<<finNm<<endl;
    ifstream ifile(finNm.c_str());
    string data;
    while(ifile>>data){
      nFiles++;
      if(beamGen)
	scaleRate += getEvents(data);

      if(testRun==1 && nFiles>10)
	break;
    }
    ifile.close();
    if(!beamGen)
      scaleRate = nFiles;
    nFiles=0;

    cout<<"\t scale rate by "<<scaleRate<<endl;

    ifile.open(finNm.c_str());
    while(ifile>>data){
      nFiles++;
      cout<<" processing: "<<data<<endl;
      nTotHits+=processOne(data);
      if(testRun==1 && nFiles>10)
	break;
    }
  }

  ofile->cd();
  otree->Write();
  ofile->Close();

  cout<<"\nFinished processing a total of "<<nTotHits<<" hits."<<endl;
}

long getEvents(string fnm){
  TFile* ifile = new TFile(fnm.c_str(),"READ");
  TTree* itree = (TTree*)ifile->Get("T");
  long evn = itree->GetEntries();
  ifile->Close();
  delete ifile;

  return evn;
}


long processOne(string fnm){

  TFile* ifile = new TFile(fnm.c_str(),"READ");
  TTree* itree = (TTree*)ifile->Get("T");

  long nEntries = itree->GetEntries();
  cout<<"\tTotal events: "<<nEntries<<endl;

  Double_t rate;
  vector<remollGenericDetectorHit_t>  *hit=0;
  itree->SetBranchAddress("rate", &rate);
  itree->SetBranchAddress("hit", &hit);

  long nHits(0);
  double tempnum;
  for (long i=0; i < nEntries;i++){
    itree->GetEntry(i);
    if(i%10000==0){
     cout<<"done with "<<i<<" entries!"<<endl;
    }
    tempnum=0;
    for(int j=0; j<hit->size(); j++){
      if(hit->at(j).det == mydet && hit->at(j).trid == 1 && (hit->at(j).r>=2174.1060 && hit->at(j).r<=2277.7732))
	{//make cut to avoid double counting ... actually can just do j=0 
	  //   	  if(hit->at(j).px!=hit->at(j-1).px&&hit->at(j).py!=hit->at(j-1).py&&hit->at(j).pz!=hit->at(j-1).pz){
	  if(tempnum==0){
	    tempnum++;//count one detector hit of trackid =1 per entry, otherwise double counting
	    newhit->push_back(hit->at(j));	
	    newrate = rate/scaleRate;
	  
	    if(rate == 0)
	      newrate = 1/scaleRate;
	    nHits++;
     	  }
	}
    }
    if(newhit->size()>0){
      otree->Fill();
      newhit->clear();
      newrate = 0;
    }
  }
  ifile->Close();
  delete ifile;

  return nHits;
}

