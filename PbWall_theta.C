/*
This program is written to see what outer radius of collimator 2 is really needed.
Chandan Ghosh April 2020

*/
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
int PbWall_theta(TString source, TString out)
{
TChain T("T");
T.Add(Form("%s",source.Data()));
/*const Int_t nfile=10;
TString add_file[nfile]={""};
for(int ij=1;ij<=nfile;ij++){
ostringstream string1;
string1<<Form("/lustre19/expphy/volatile/halla/moller12gev/chandan/sim_out/Collar_moller_NoShld/remollout_moller%d.root",ij);
//string1<<Form("remollout_beam%d.root",ij);
add_file[ij-1]=string1.str();
//cout<<string1.str()<<endl;
T.Add(add_file[ij-1]);
}*/
//T.Add("/lustre19/expphy/volatile/halla/parity/chandan/sim_out/ExistingColl/remollout_beam*.root");
Int_t nEvents = T.GetEntries();
std::cout<<"Number of entries  "<<nEvents<<std::endl;
Double_t weight = 1.;///nEvents;
//weight = 1.;
//TFile f("/lustre19/expphy/volatile/halla/parity/chandan/sim_out/DSConfig0Shield_ep/DSConfig0Shield_ep_disk_28.root","RECREATE");
TFile f(Form("%s",out.Data()),"RECREATE");
set_plot_style();

std::map<TString,TH2D*> h_d_xy;
std::map<TString,TH2D*> h_d_rz;
std::map<TString,TH1D*> h_d_r;
std::map<TString,TH1D*> h_d_p;
std::map<TString,TH1D*> h_d_th;

std::map<TString,TH2D*> h_d1_xy;
std::map<TString,TH2D*> h_d1_rz;
std::map<TString,TH2D*> h_d1_rph;
std::map<TString,TH1D*> h_d1_r;
std::map<TString,TH1D*> h_d1_p;
std::map<TString,TH1D*> h_d_ph;


TString part,part1,part2,part3,part4;
const Int_t nEnergy=1;
//0:Total; 1:hit.p<1; 2: hit.p>=1 && hit.p<10;  3: hit.p>=10 && hit.p<100; 4: hit.p>=100;
const Int_t nParticle=1;
const Int_t nDet=2;

string sParticle[nParticle] = {"electron"};//,"positron","photon","neutron"};
std::map<int,string> snParticle  {{11,"electron"}};//,{-11,"positron"},{22,"photon"},{2112,"neutron"}};
std::map<int,string> snDet  {{74,"PbWallEnt"},{75,"PbWallExit"}};
string sDet[nDet] = {"PbWallEnt","PbWallExit"};


for(Int_t k=0;k<nDet;k++){
 for(Int_t i=0;i<nParticle;i++){
  part=Form("h_%s_%s",sDet[k].c_str(),sParticle[i].c_str());
  h_d_p[part] = new TH1D(part+"_d_p",Form("Momemtum dist. at  %s for %s",sDet[k].c_str(),sParticle[i].c_str()),1200,0,12000);
  h_d_rz[part] = new TH2D(part+"_d_rz",Form("Vertex rz dist. at %s for %s",sDet[k].c_str(),sParticle[i].c_str()),3800,-6000,32000,300,-300,300);
  h_d_xy[part] = new TH2D(part+"_d_xy",Form("xy dist. at %s for %s",sDet[k].c_str(),sParticle[i].c_str()),300,-310.,-10.,125,-5,120);
  h_d_r[part] = new TH1D(part+"_d_r",Form("R dist. at  %s for %s",sDet[k].c_str(),sParticle[i].c_str()),60,0,300);
  h_d_th[part] = new TH1D(part+"_d_th",Form("Scatt. Angle at  %s for %s",sDet[k].c_str(),sParticle[i].c_str()),75,0,15.0);

  h_d1_p[part] = new TH1D(part+"_d1_p",Form("Momemtum dist. at  %s for %s",sDet[k].c_str(),sParticle[i].c_str()),1200,0,12000);
  h_d1_rz[part] = new TH2D(part+"_d1_rz",Form("Vertex rz dist. at %s for %s",sDet[k].c_str(),sParticle[i].c_str()),3800,-6000,32000,200,-200,200);
  h_d1_rph[part] = new TH2D(part+"_d1_rph",Form("r vs Ph at %s for %s",sDet[k].c_str(),sParticle[i].c_str()),300,0,300,100,-1.,1.0);
  h_d1_xy[part] = new TH2D(part+"_d1_xy",Form("xy dist. at %s for %s",sDet[k].c_str(),sParticle[i].c_str()),300,-300.,300.,300,-300,300);
  h_d1_r[part] = new TH1D(part+"_d1_r",Form("R dist. at  %s for %s",sDet[k].c_str(),sParticle[i].c_str()),60,0,300);
  h_d_ph[part] = new TH1D(part+"_d_ph",Form("Phi. Angle at  %s for %s",sDet[k].c_str(),sParticle[i].c_str()),100,-1.0,1.0);
}
}


Double_t fRate=0;
remollEvent_t *fEvent=0;
std::vector<remollGenericDetectorHit_t> *fHit=0;
std::vector<remollEventParticle_t> *fPart=0;

T.SetBranchAddress("hit", &fHit);
T.SetBranchAddress("rate", &fRate);

for(size_t j=0;j<nEvents;j++){
  T.GetEntry(j);
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
  /*for(it=det28trid.begin();it!=det28trid.end();++it)
        std::cout<<"Inside main event number "<<j<<" track id "<<*it<<std::endl;*/
  for(size_t i=0;i<fHit->size();i++){
    remollGenericDetectorHit_t hit=fHit->at(i); 

    Int_t det = hit.det;
    Int_t pid = abs(hit.pid);
    Int_t trid = hit.trid;
    bool test_th = det1_th && det2_th;//these are implemented to avoid particles which are scattered multiple times 
    bool test_ph = det1_ph && det2_ph; 
    //if(!((det==76 || det==77) && hit.vz<=-3875 && hit.pz>=0)) continue;
    if(!(det==74 && hit.vz<=-3875 && hit.mtrid==0 && hit.pz>=0 && hit.p>500)) continue;
    col2it = find(col2trid.begin(),col2trid.end(),trid);

    //if(!(col2it!=col2trid.end())) {/*std::cout<<"Inside Continue loop::Event "<<j<<" hit "<<i<<" det "<<det<<" det28trid size "<<det28trid.size()<<" col2 size "<<col2trid.size()<<" col4exit size "<<col4trid.size()<<std::endl;*/continue;}
    //std::cout<<"After Continue :: Event "<<j<<" hit "<<i<<" trid "<<trid<<" iterator find "<<std::endl;  
    if((det==74 ) && pid==11){
//	std::cout<<"Hi Event "<<j<<" hit "<<i<<" det "<<det<<" pid "<<pid<<" r "<<hit.r<<" trid "<<hit.trid<<" mtrid  "<<hit.mtrid<<" phi "<< (180./acos(-1))*atan(hit.y/hit.x)<<std::endl;
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
		//if(hit.r>98 || hit.r<95)
		//std::cout<<"Event # "<<j<<" hit #"<<i<<" detector # "<<hit.det<<" hit radius "<<hit.r<<" p "<<hit.p<<" pid "<<hit.pid<<" vz "<<hit.vz<<" trid "<<hit.trid<<" hit.mtrid "<<hit.mtrid<<" pz "<<hit.pz<<std::endl;
		h_d_xy[part]->Fill(hitx,hity,(fRate)*weight);
		h_d_rz[part]->Fill(hit.vz,sqrt(hit.vx*hit.vx+hit.vy*hit.vy),(fRate)*weight);
		h_d_r[part]->Fill(hit.r,(fRate)*weight);
		h_d_p[part]->Fill(hit.p,(fRate)*weight);
		double theta = (180./pi)*atan2(sqrt(px*px+py*py),pz);
		h_d_th[part]->Fill(theta,(fRate)*weight);}
		//std::cout<<std::endl;
		}	
		}
	}

}


for(Int_t k=0;k<1;k++){
  for(Int_t i=0;i<nParticle;i++){
      part=Form("h_%s_%s",sDet[k].c_str(),sParticle[i].c_str());
      h_d_xy[part]->Write("", TObject::kOverwrite);
      h_d_p[part]->Write("", TObject::kOverwrite);
      h_d_rz[part]->Write("", TObject::kOverwrite);
      h_d_r[part]->Write("", TObject::kOverwrite);
      h_d_th[part]->Write("", TObject::kOverwrite);

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
    //double angle = atan2(hity,hitx); 
    //if(det==68 && angle<2.888) std::cout<<"det "<<det<<" x "<<hitx<<" y "<<hity<<" phi "<<angle<<std::endl;
    /*if(det==68 && angle<2.888 && abs(pid)==11){
	//std:;cout<<"Inside isvalid1 det "<<det<<" r "<<hit.r<<"   track id "<<hit.trid<<std::endl;
	CutPlanetrid.push_back(hit.trid);
    }
    if(det==28 && hit.r>=640 &&hit.r<=1195 && pid==11){
	//std:;cout<<"Inside isvalid1 det "<<det<<" r "<<hit.r<<"   track id "<<hit.trid<<std::endl;
	det28trid.push_back(hit.trid);
    }*/
    if(det==74 && hit.r>=227 && hit.r<=235 && pid==11 && hit.vz<=-3875 && hit.pz>=0 && hit.p>500){//97.7 -33.7-35
	//std:;cout<<"Inside isvalid1 det "<<det<<" r "<<hit.r<<"   track id "<<hit.trid<<std::endl;
	col2trid.push_back(hit.trid);
    }
 }
 // std::cout<<"end of isvalid"<<std::endl;
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
