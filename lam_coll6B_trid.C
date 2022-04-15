#include "remolltypes.hh"
void isValid1(std::vector<remollGenericDetectorHit_t> *fHit, std::vector<int> &lamtrid);
const double pi = TMath::Pi();
//void RotateXY(double &x, double &y);
//double getAngle(double x, double y);

void lam_coll6B_trid(){
  map<TString,TH2D*> h_d_xy;
  map<TString,TH2D*> h_d_rz;
  map<TString,TH1D*> h_d_r;
  map<TString,TH1D*> h_d_p;
  map<TString,TH1D*> h_d_th;

  TString part;
  string sParticle[] = {"electron"};
  const int nParticle = sizeof(sParticle)/sizeof(*sParticle);
  string sDet[] = {"coll6BEnt"/*,"coll6BExit","collar2Ent"*/};
  const int nDet = sizeof(sDet)/sizeof(*sDet);
  map<int,string> snParticle {{11,"electron"}};
  map<int,string> snDet {{65,"coll6BEnt"}/*,{66,"coll6BExit"},{5719,"collar2Ent"}*/};
  double rmin, rmax;
  for(int iDet=0;iDet<nDet;iDet++){
   if(sDet[iDet]=="coll6BEnt"){rmin=300; rmax=800;}
   if(sDet[iDet]=="coll6BExit"){rmin=300; rmax=800;}
   if(sDet[iDet]=="collar2Ent"){rmin=800; rmax=1300;}
   for(int iParticle=0;iParticle<nParticle;iParticle++){
    part=Form("h_%s_%s",sDet[iDet].c_str(),sParticle[iParticle].c_str());
    h_d_xy[part] = new TH2D(part+"_d_xy",Form("xy dist. at %s for %s",sDet[iDet].c_str(),sParticle[iParticle].c_str()),400,-rmax,rmax,400,-rmax,rmax);
    h_d_rz[part] = new TH2D(part+"_d_rz",Form("Vertex rz dist. at %s for %s",sDet[iDet].c_str(),sParticle[iParticle].c_str()),3800,-6000,32000,400,-rmax,rmax);
    h_d_r[part] = new TH1D(part+"_d_r",Form("R dist. at %s for %s",sDet[iDet].c_str(),sParticle[iParticle].c_str()),500,rmin,rmax);
    h_d_p[part] = new TH1D(part+"_d_p",Form("Momentum dist. at %s for %s",sDet[iDet].c_str(),sParticle[iParticle].c_str()),1200,0,12000);
    h_d_th[part] = new TH1D(part+"_d_th",Form("Scatt. Angle dist. at %s for %s",sDet[iDet].c_str(),sParticle[iParticle].c_str()),100,0,15.0);
   }
  }

  int nfile=0;
  Long64_t nentry=0;
  long nTotEv=0;
  TString rootfile_dir = "/volatile/halla/parity/adhidevi/remoll_rootfiles/develop_br";
  string tgt_gen_config = "LH2_beam_V11";
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
      std::vector<int> lamtrid;
      std::vector<int>::iterator lamit;
      if(fRate==0) fRate=1;
      isValid1(fHit,lamtrid);
      for(size_t i=0;i<fHit->size();i++){
        remollGenericDetectorHit_t hit = fHit->at(i);
        Int_t det = hit.det;
        Int_t pid = hit.pid;
        Int_t trid = hit.trid;
        lamit = find(lamtrid.begin(),lamtrid.end(),trid);
        if(det==65 && pid==11 && hit.vz<=-3875 && hit.pz>=0){
          double hitx = hit.x;
          double hity = hit.y;
          double hite = hit.e;
          double px = hit.px;
          double py = hit.py;
          double pz = hit.pz;
          part = Form("h_%s_%s",snDet[det].c_str(),snParticle[pid].c_str());

          if(lamit!=lamtrid.end()){
            if(trid==*lamit){
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

  string outfile = Form("./rootfiles/%s_lam_%s.root",tgt_gen_config.c_str(),sDet[0].c_str());
  TFile* fout = TFile::Open(outfile.c_str(),"RECREATE");

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
void isValid1(std::vector<remollGenericDetectorHit_t> *fHit, std::vector<int> &lamtrid){
  for(size_t i=0;i<fHit->size();i++){
    remollGenericDetectorHit_t hit=fHit->at(i);
    int det = hit.det;
    int pid = hit.pid;
    double phi = hit.ph;
    if(phi<0) phi +=2.0*pi;
    double modphi = fmod(phi,2.0*pi/7.0);
    if(modphi<3.0*pi/28.0 || modphi>5.0*pi/28.0) continue;

    if(det==174 && pid==11 && hit.vz<=-3875 && hit.pz>=0 && (hit.r>=1010 && hit.r<=1130)){
      lamtrid.push_back(hit.trid);
    }
  }
}
/*
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
*/
