#include "remolltypes.hh"
//void set_plot_style();
void isValid1(std::vector<remollGenericDetectorHit_t> *fHit, std::vector<int> &collar2trid);
void RotateXY(double &x, double &y);
double getAngle(double x, double y);
const double pi = TMath::Pi();
const double septant = 2*pi/7.0;
//const double midangle = 360.0/14.0;//not in use
//const double septantStart = 3*septant;//not in use
//const double septantStop = septantStazt+septant;//not in use

void collar2_outerRing_taperAngle(){
  map<TString,TH2D*> h_d_xy;
  map<TString,TH2D*> h_d_rz;
  map<TString,TH1D*> h_d_r;
  map<TString,TH1D*> h_d_p;
  map<TString,TH1D*> h_d_th;

  TString part;
  string sParticle[] = {"electron"};
  const int nParticle = sizeof(sParticle)/sizeof(*sParticle);
  string sDet[] = {"Ent"};
  const int nDet = sizeof(sDet)/sizeof(*sDet);
  map<int,string> snParticle {{11,"electron"}};
  map<int,string> snDet {{179,"Ent"}};

  for(int iDet=0;iDet<nDet;iDet++){
   for(int iParticle=0;iParticle<nParticle;iParticle++){
    part=Form("h_%s_%s",sDet[iDet].c_str(),sParticle[iParticle].c_str());
    h_d_xy[part] = new TH2D(part+"_d_xy",Form("xy dist. at %s for %s",sDet[iDet].c_str(),sParticle[iParticle].c_str()),300,-1200,-900,300,-50,250);
    h_d_rz[part] = new TH2D(part+"_d_rz",Form("Vertex rz dist. at %s for %s",sDet[iDet].c_str(),sParticle[iParticle].c_str()),3800,-6000,32000,100,-1200,1200);
    h_d_r[part] = new TH1D(part+"_d_r",Form("R dist. at %s for %s",sDet[iDet].c_str(),sParticle[iParticle].c_str()),300,900,1200);
    h_d_p[part] = new TH1D(part+"_d_p",Form("Momentum dist. at %s for %s",sDet[iDet].c_str(),sParticle[iParticle].c_str()),1200,0,12000);
    h_d_th[part] = new TH1D(part+"_d_th",Form("Scatt. Angle dist. at %s for %s",sDet[iDet].c_str(),sParticle[iParticle].c_str()),100,0,15.0);
   }
  }

  int nfile=0;
  Long64_t nentry=0;
  long nTotEv=0;
  TString rootfile_dir = "/volatile/halla/moller12gev/devi/remoll_rootfiles/develop_br";
  string tgt_gen_config = "LH2_beam_V7";
  for(int ifile=1001;ifile<=6000;ifile++){
    string infile = Form("%s/%s/%s_%d.root",rootfile_dir.Data(),tgt_gen_config.c_str(),tgt_gen_config.c_str(),ifile);
    ifstream inf(infile.c_str());
    if(!inf){
      cout<<Form("Skipping %s. File doesn't exits.",infile.c_str())<<endl;
      continue;
    }
    TFile * fin = TFile::Open(infile.c_str(),"READ");
    if(fin->TestBit(TFile::kRecovered)){
      cout<<Form("Skippint %s. Recovered file.",infile.c_str())<<endl;
      continue;
    }
    nfile++;
    TTree* T = (TTree*)fin->Get("T");
    if(T==0) return 0;

    nentry = T->GetEntries();
    cout<<Form("Found %lld entries in filesplit %d",nentry,nfile)<<endl;
    nTotEv+=nentry;

    Double_t fRate=0;
    remollEvent_t *fEvent=0;
    std::vector<remollGenericDetectorHit_t> *fHit=0;
    T->SetBranchAddress("hit", &fHit);
    T->SetBranchAddress("rate", &fRate);

    for(size_t j=0;j<nentry;j++){
      T->GetEntry(j);
      std::vector<int> collar2trid;
      std::vector<int>::iterator collar2it;
      if(fRate==0) fRate=1;
      isValid1(fHit,collar2trid);
      for(size_t i=0;i<fHit->size();i++){
        remollGenericDetectorHit_t hit = fHit->at(i);
        Int_t det = hit.det;
        Int_t pid = hit.pid;
        Int_t trid = hit.trid;
        if(!(det==179 && hit.vz<=-3875 && hit.pz>=0)) continue;
        collar2it = find(collar2trid.begin(),collar2trid.end(),trid);
        if(det==179 && pid==11){
          double hitx = hit.x;
          double hity = hit.y;
          double hite = hit.e;
          double px = hit.px;
          double py = hit.py;
          double pz = hit.pz;
          part = Form("h_%s_%s",snDet[det].c_str(),snParticle[pid].c_str());

          if(hity<0)hity=-hity;
          RotateXY(hitx,hity);
          if(collar2it!=collar2trid.end()){
            if(trid==*collar2it){
              h_d_xy[part]->Fill(hitx,hity,fRate);
              h_d_rz[part]->Fill(hit.vz,sqrt(hit.vx*hit.vx+hit.vy*hit.vy),fRate);
              h_d_r[part]->Fill(hit.r,fRate);
              h_d_p[part]->Fill(hit.p,fRate);
              double theta = (180.0/pi)*atan2(sqrt(px*px+py*py),pz);
              h_d_th[part]->Fill(theta,fRate);
            }
          }
        }
      }
    }
    delete fin;
  }
  
  cout<<Form("Total number of file splits: %d",nfile)<<endl;
  cout<<Form("Total number of entries: %ld",nTotEv)<<endl;

  string outfile = Form("./rootfiles/%s_theta_1.root",tgt_gen_config.c_str());
  TFile* fout = TFile::Open(outfile.c_str(),"RECREATE");
//  set_plot_style();

  for(int iDet=0;iDet<nDet;iDet++){
    for(int iParticle=0;iParticle<nParticle;iParticle++){
      part = Form("h_%s_%s",sDet[iDet].c_str(),sParticle[iParticle].c_str());
      h_d_xy[part]->Write("",TObject::kOverwrite);
      h_d_rz[part]->Write("",TObject::kOverwrite);
      h_d_r[part]->Write("",TObject::kOverwrite);
      h_d_p[part]->Write("",TObject::kOverwrite);
      h_d_th[part]->Write("",TObject::kOverwrite);
    }
  }
}
/*
void set_plot_style(){
  const int NRGBs = 5;
  const int NCont = 255;
  Double_t stops[NRGBs] = {0.00,0.34,0.61,0.84,1.00};
  Double_t red[NRGBs] = {0.00,0.00,0.87,1.00,0.51};
  Double_t green[NRGBs] = {0.00,0.81,1.00,0.20,0.00};
  Double_t blue[NRGBs] = {0.51,1.00,0.12,0.00,0.00};
  TColor::CreateGradientColorTable(NRGBs,stops,red,green,blue,NCont);
  gStyle->SetNumberContours(NCont);
}
*/
void isValid1(std::vector<remollGenericDetectorHit_t> *fHit, std::vector<int> &collar2trid){
  for(size_t i=0;i<fHit->size();i++){
    remollGenericDetectorHit_t hit=fHit->at(i);
    int det = hit.det;
    int pid = hit.pid;
    double hite = hit.e;
    double hitx = hit.x;
    double hity = hit.y;
    if(hity<0)hity=-hity;
    RotateXY(hitx,hity);
    double angle = getAngle(hitx,hity);
    if(det==179 && hit.r>=1010 && hit.r<=1133.65 && pid==11 && hit.vz<=-3875 && hit.pz>=0 && hit.p>500){
      collar2trid.push_back(hit.trid);
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
  double tx, ty;
  if(angle>=0 && angle<septant){
    tx = x*c3-y*s3;
    ty = x*s3+y*c3;
    x = tx;
    y = (ty<0) ? -ty : ty;
  }
  if(angle>=septant && angle<2*septant){
    tx = x*c2-y*s2;
    ty = x*s2+y*c2;
    x = tx;
    y = (ty<0) ? -ty : ty;
  }
  if(angle>=2*septant && angle<3*septant){
    tx = x*c1-y*s1;
    ty = x*s1+y*c1;
    x = tx;
    y = (ty<0) ? -ty : ty;
  }
  return;
}

double getAngle(double x, double y){
  double angle=atan2(y,x);
  return (angle<0) ? (2*pi+angle) : angle;
}
