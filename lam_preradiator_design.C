#include "remolltypes.hh"
void set_plot_style();
void isValid1(std::vector<remollGenericDetectorHit_t> *fHit, std::vector<int> &col2trid, std::vector<int> &col2trid1);
void RotateXY(double &x,double &y);
double getAngle(double x, double y);
const double pi = acos(-1);
const double septant = (2*pi/7.);
const double midangle = 360./14.;
const double septantStart = 3*septant;
const double septantStop = septantStart+septant;

void lam_preradiator_design(){
TString source = "/volatile/halla/parity/adhidevi/remoll_rootfiles/develop_br/LH2_beam_V13/LH2_beam_V13_1001.root";
TString out = "./rootfiles/test_LH2_beam_V13.root";
TChain T("T");
T.Add(Form("%s",source.Data()));
Int_t nEvents = T.GetEntries();
std::cout<<"Number of entries  "<<nEvents<<std::endl;
TFile f(Form("%s",out.Data()),"RECREATE");
set_plot_style();
TTree skimTree("skimTree","skimTree");

Int_t det(0), pid(0), trid(0), mtrid(0), hitid(0), event(0);
Double_t hite(-1.e-12), hitr(-1.e-12), hitt(-1.e-12), hitz(-1.e-12), hitvz(-1.e-12);

TBranch* b_det_id = skimTree.Branch("det_id",&det,"det_id/I");
TBranch* b_pid = skimTree.Branch("pid",&pid,"pid/I");
TBranch* b_trid = skimTree.Branch("trid",&trid,"trid/I");
TBranch* b_mtrid = skimTree.Branch("mtrid",&mtrid,"mtrid/I");
TBranch* b_hit_id = skimTree.Branch("hit_id",&hitid,"hit_id/I");
TBranch* b_event_id = skimTree.Branch("event_id",&event,"event_id/I");
TBranch* b_hite = skimTree.Branch("hite",&hite,"hite/D");
TBranch* b_hitr = skimTree.Branch("hitr",&hitr,"hitr/D");
TBranch* b_hitt = skimTree.Branch("hitt",&hitt,"hitt/D");
TBranch* b_hitz = skimTree.Branch("hitz",&hitz,"hitz/D");
TBranch* b_hitvz = skimTree.Branch("hitvz",&hitvz,"hitvz/D");

std::map<int,string> snParticle  {{11,"electron"}};//,{-11,"positron"},{22,"photon"},{2112,"neutron"}};
std::map<int,string> snDet  {{65,"6BEnt"},{66,"6BExit"},{5719,"Coll2Ent"},{2174,"Coll2"},{1174,"PreRad"},{174,"LAM"}};

TString sDet[] = {"6BEnt","6BExit","Coll2Ent","Coll2","PreRad","LAM"};

const Int_t nDet=sizeof(sDet)/sizeof(*sDet);;

long ecounter =0;
long gcounter =0;
long ncounter =0;
Double_t fRate=0;
remollEvent_t *fEvent=0;
std::vector<remollGenericDetectorHit_t> *fHit=0;
std::vector<remollEventParticle_t> *fPart=0;

T.SetBranchAddress("hit", &fHit);
T.SetBranchAddress("rate", &fRate);

for(size_t j=0;j<nEvents;j++){
  T.GetEntry(j);
  event = j;

  std::vector<int> col2trid;
  std::vector<int>::iterator col2it;
  std::vector<int> col2trid1;
  std::vector<int>::iterator col2it1;

  if(fRate==0) fRate=1;
  isValid1(fHit,col2trid,col2trid1);

////////////////////////////////
  Int_t adet[fHit->size()];
  Int_t apid[fHit->size()];
  Int_t atrid[fHit->size()];
  Int_t amtrid[fHit->size()];
////////////////////////////////

  for(size_t i=0;i<fHit->size();i++){
    remollGenericDetectorHit_t hit=fHit->at(i); 

    det = hit.det;
    pid = abs(hit.pid);
    trid = hit.trid;
    mtrid = hit.mtrid;
    hitid = i;
    hite = hit.e;
    hitr = hit.r;
    hitt = hit.t;
    hitz = hit.z;
    hitvz = hit.vz;

////////////////////////////////////
    adet[i] = hit.pid;
    apid[i] = abs(hit.pid);
    atrid[i] = hit.trid;
    amtrid[i] = hit.mtrid;
    if(adet[i-1]==adet[i] && apid[i-1]==apid[i] && atrid[i-1]==atrid[i] && amtrid[i-1]==amtrid[i]) continue;
////////////////////////////////////

    if(!((det==65 || det==66 || det==5719 || det==174 || det==1174 || det==2174) && hit.pz>0)) continue;
    col2it = find(col2trid.begin(),col2trid.end(),trid);
    col2it1 = find(col2trid1.begin(),col2trid1.end(),trid);
    if(!(col2it!=col2trid.end())) continue;
    if((det==65 || det==66 || det==5719 || det==174 || det==1174 || det==2174)){
      if(pid==11)ecounter++;
      if(pid==22)gcounter++;
      if(pid==2112)ncounter++;
      b_det_id->Fill();
      b_pid->Fill();
      b_trid->Fill();
      b_mtrid->Fill();
      b_hit_id->Fill();
      b_event_id->Fill();
      b_hite->Fill();
      b_hitr->Fill();
      b_hitt->Fill();
      b_hitz->Fill();
      b_hitvz->Fill();
      std::cout<<"Event "<<event<<" hit# "<<hitid<<" det "<<det<<" pid "<<pid<<" time "<<hitt<<" r "<<hitr<<" E "<<hite<<" trid "<<trid<<" mtrid  "<<mtrid<<" hit.z "<< hitz<<"  hit.vz "<<hitvz<<std::endl;
    }
  }
}
cout<<"Electron Counter "<<ecounter<<"  Gamma Counter "<<gcounter<<"  Neutron Counter "<<ncounter<<endl;
skimTree.SetEntries();
f.cd();
skimTree.Write();
f.Close();
}

void set_plot_style(){
  const Int_t NRGBs = 5;
  const Int_t NCont = 255;
  Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  gStyle->SetNumberContours(NCont);
}

void isValid1(std::vector<remollGenericDetectorHit_t> *fHit, std::vector<int> &col2trid, std::vector<int> &col2trid1){
  for(size_t i=0;i<fHit->size();i++){
    remollGenericDetectorHit_t hit=fHit->at(i); 
    int det = hit.det; 
    int pid = hit.pid; 
    int hite = hit.e;
    double hitx = hit.x;
    double hity = hit.y;
    if(hity<0)hity=-hity;
    RotateXY(hitx,hity);
    double angle = getAngle(hitx,hity);
    if(det==1174 && hit.r>=1010 && hit.r<=1130 && hit.pz>0 && hit.vz<19087){//97.7 -33.7-35
      col2trid.push_back(hit.trid);
    }
    if(det==174 && hit.r>=1010 && hit.r<=1130 && hit.pz>0){
      col2trid1.push_back(hit.trid);
    }
  }
}

void RotateXY(double &x, double &y){
  double angle = getAngle(x,y);
  const double s1=sin(septant);
  const double c1=cos(septant);
  const double s2=sin(2*septant);
  const double c2=cos(2*septant);
  const double s3=sin(3*septant);
  const double c3=cos(3*septant);
  double tx,ty;
  if (angle>=0 && angle<septant){
    tx = x*c3-y*s3;
    ty = x*s3+y*c3;
    x = tx;
    y = (ty<0) ? -ty : ty;
  }
  if (angle>=septant && angle<2*septant){
    tx = x*c2-y*s2;
    ty = x*s2+y*c2;
    x = tx;
    y = (ty<0) ? -ty : ty;
  }
  if (angle>=2*septant && angle<3*septant){
    tx = x*c1-y*s1;
    ty = x*s1+y*c1;
    x=tx;
    y=(ty<0) ? -ty : ty;
    }
  return;
}

double getAngle(double x, double y){
  double angle=atan2(y,x);
  return (angle<0) ? (2*pi+angle) : angle;
}
