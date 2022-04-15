#include "remolltypes.hh"
void set_plot_style();
void isValid1(std::vector<remollGenericDetectorHit_t> *fHit, std::vector<int> &col2trid);
void RotateXY(double &x,double &y);
double getAngle(double x, double y);
const double pi = acos(-1);
const double septant = (2*pi/7.);
const double midangle = 360./14.;
const double septantStart = 3*septant;
const double septantStop = septantStart+septant;

int theta_dist(){
TString source = "/volatile/halla/moller12gev/devi/remoll_rootfiles/develop_br/LH2_beam_V7/LH2_beam_V7_*.root";
TString out = "./rootfiles/theta_dist_LH2_beam_V7.root";
TChain* T = new TChain("T");
T->Add(Form("%s",source.Data()));

Int_t nEvents = T->GetEntries();
std::cout<<"Number of entries  "<<nEvents<<std::endl;
Double_t weight = 1.;///nEvents;
TFile f(Form("%s",out.Data()),"RECREATE");
set_plot_style();

TString part;

TString sParticle[] = {"electron"};//,"positron","photon","neutron"};
const Int_t nParticle=sizeof(sParticle)/sizeof(*sParticle);
std::map<int,string> snParticle  {{11,"electron"}};//,{-11,"positron"},{22,"photon"},{2112,"neutron"}};
std::map<int,string> snDet  {{179,"Collar2Ent"}};
TString sDet[] = {"Collar2OuterRingEnt"};
const Int_t nDet = sizeof(sDet)/sizeof(*sDet);

TH2D* h_d_xy[nDet][nParticle];
TH2D* h_d_rz[nDet][nParticle];
TH1D* h_d_r[nDet][nParticle];
TH1D* h_d_p[nDet][nParticle];
TH1D* h_d_th[nDet][nParticle];

TH2D* h_d1_xy[nDet][nParticle];
TH2D* h_d1_rz[nDet][nParticle];
TH2D* h_d1_rph[nDet][nParticle];
TH1D* h_d1_r[nDet][nParticle];
TH1D* h_d1_p[nDet][nParticle];
TH1D* h_d_ph[nDet][nParticle];

for(Int_t k=0;k<nDet;k++){
 for(Int_t i=0;i<nParticle;i++){
  part=Form("h_%s_%s",sDet[k].Data(),sParticle[i].Data());
  h_d_p[k][i] = new TH1D(part+"_d_p",Form("Momemtum dist. at  %s for %s",sDet[k].Data(),sParticle[i].Data()),1200,0,12000);
  h_d_rz[k][i] = new TH2D(part+"_d_rz",Form("Vertex rz dist. at %s for %s",sDet[k].Data(),sParticle[i].Data()),3800,-6000,32000,300,-300,300);
  h_d_xy[k][i] = new TH2D(part+"_d_xy",Form("xy dist. at %s for %s",sDet[k].Data(),sParticle[i].Data()),300,-310.,-10.,125,-5,120);
  h_d_r[k][i] = new TH1D(part+"_d_r",Form("R dist. at  %s for %s",sDet[k].Data(),sParticle[i].Data()),60,0,300);
  h_d_th[k][i] = new TH1D(part+"_d_th",Form("Scatt. Angle at  %s for %s",sDet[k].Data(),sParticle[i].Data()),75,0,15.0);

  h_d1_p[k][i] = new TH1D(part+"_d1_p",Form("Momemtum dist. at  %s for %s",sDet[k].Data(),sParticle[i].Data()),1200,0,12000);
  h_d1_rz[k][i] = new TH2D(part+"_d1_rz",Form("Vertex rz dist. at %s for %s",sDet[k].Data(),sParticle[i].Data()),3800,-6000,32000,200,-200,200);
  h_d1_rph[k][i] = new TH2D(part+"_d1_rph",Form("r vs Ph at %s for %s",sDet[k].Data(),sParticle[i].Data()),300,0,300,100,-1.,1.0);
  h_d1_xy[k][i] = new TH2D(part+"_d1_xy",Form("xy dist. at %s for %s",sDet[k].Data(),sParticle[i].Data()),300,-300.,300.,300,-300,300);
  h_d1_r[k][i] = new TH1D(part+"_d1_r",Form("R dist. at  %s for %s",sDet[k].Data(),sParticle[i].Data()),60,0,300);
  h_d_ph[k][i] = new TH1D(part+"_d_ph",Form("Phi. Angle at  %s for %s",sDet[k].Data(),sParticle[i].Data()),100,-1.0,1.0);
}
}


Double_t fRate=1;
remollEvent_t *fEvent=0;
std::vector<remollGenericDetectorHit_t> *fHit=0;
std::vector<remollEventParticle_t> *fPart=0;

T->SetBranchAddress("hit", &fHit);
//T->SetBranchAddress("rate", &fRate);

for(size_t j=0;j<nEvents;j++){
  T->GetEntry(j);
  if(j%=5000)cout<<"Analyzed: "<<j<<" events!"<<endl;
  bool det1_th=false;
  bool det2_th=false;
  bool det1_ph=false;
  bool det2_ph=false;
  std::vector<int> col2trid;
  std::vector<int>::iterator col2it;
  std::vector<int> col2trid1;
  std::vector<int>::iterator col2it1;

  std::vector<int> col2trid2;//col2 exit
  std::vector<int>::iterator col2it2;
  std::vector<int> col2trid3;//coll2 exit
  std::vector<int>::iterator col2it3;

  if(fRate==0) fRate=1;
  isValid1(fHit,col2trid);
  for(size_t i=0;i<fHit->size();i++){
    remollGenericDetectorHit_t hit=fHit->at(i); 

    Int_t det = hit.det;
    Int_t pid = abs(hit.pid);
    Int_t trid = hit.trid;
    bool test_th = det1_th && det2_th;//these are implemented to avoid particles which are scattered multiple times 
    if(!(det==179 && hit.vz<=-3875 && hit.mtrid==0 && hit.pz>=0 && hit.p>500)) continue;
    col2it = find(col2trid.begin(),col2trid.end(),trid);

    if((det==179 ) && pid==11){
    double hitx = hit.x; 
    double hity = hit.y;
    double hite = hit.e;
    double px = hit.px; 
    double py = hit.py;
    double pz = hit.pz;
    //int trid = hit.trid;
    hite=1.;
    part=Form("h_%s_%s",snDet[det].c_str(),snParticle[pid].c_str());

    if(hity<0)hity=-hity;
    RotateXY(hitx,hity);
    if(col2it!=col2trid.end() && !test_th){
    if(trid==*col2it){
    h_d_xy[0][0]->Fill(hitx,hity,(fRate)*weight);
    h_d_rz[0][0]->Fill(hit.vz,sqrt(hit.vx*hit.vx+hit.vy*hit.vy),(fRate)*weight);
    h_d_r[0][0]->Fill(hit.r,(fRate)*weight);
    h_d_p[0][0]->Fill(hit.p,(fRate)*weight);
    double theta = (180./pi)*atan2(sqrt(px*px+py*py),pz);
    h_d_th[0][0]->Fill(theta,(fRate)*weight);}
    }    
    }
}

}


for(Int_t k=0;k<1;k++){
  for(Int_t i=0;i<nParticle;i++){
      part=Form("h_%s_%s",sDet[k].Data(),sParticle[i].Data());
      h_d_xy[k][i]->Write("", TObject::kOverwrite);
      h_d_p[k][i]->Write("", TObject::kOverwrite);
      h_d_rz[k][i]->Write("", TObject::kOverwrite);
      h_d_r[k][i]->Write("", TObject::kOverwrite);
      h_d_th[k][i]->Write("", TObject::kOverwrite);
  }
}


return 0;
}
void set_plot_style()
{
    const Int_t NRGBs = 5;
    const Int_t NCont = 255;
    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
    Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
    Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);
}

void isValid1(std::vector<remollGenericDetectorHit_t> *fHit, std::vector<int> &col2trid)
{
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
    if(det==179 && hit.r>=1010 && hit.r<=1133.65 && pid==11 && hit.vz<=-3875 && hit.pz>=0 && hit.p>500){//97.7 -33.7-35
    col2trid.push_back(hit.trid);
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
    x=tx;
    y=(ty<0) ? -ty : ty;
    }
if (angle>=septant && angle<2*septant){
    tx = x*c2-y*s2;
    ty = x*s2+y*c2;
    x=tx;
    y=(ty<0) ? -ty : ty;
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
