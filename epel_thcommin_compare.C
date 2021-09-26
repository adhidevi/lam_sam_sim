#include "remolltypes.hh"
void epel_thcommin_compare(){
    TString dir = "$VOLATILE/remoll_rootfiles/PhotonBlocker/";
    TString files[] = {"thcommin02_epel/thcommin02_epel_*.root","thcommin002_epel/thcommin002_epel_*.root","thcommin0002_epel/thcommin0002_epel_*.root"};
    const int nfiles = sizeof(files)/sizeof(*files);
    TH1D* h_radial[nfiles];
    double thcommin[nfiles] = {0.2,0.02,0.002};
    int color[nfiles] = {1,2,4};
    for(int ifiles=0;ifiles<nfiles;ifiles++){
    h_radial[ifiles] = new TH1D(Form("h_thcommin_%.3f",thcommin[ifiles]),Form("Different thcommin rate-weighted radial distribution (epel);r (mm);rate (GHz/65uA)"),500,500,1500);
    h_radial[ifiles]->SetLineColor(color[ifiles]);
    h_radial[ifiles]->Sumw2();

    TChain* T = new TChain("T");
    T->Add(dir+Form("%s",files[ifiles].Data()));
    Long64_t nentry = T->GetEntries();
    std::vector<remollGenericDetectorHit_t> *fHit = 0;
    remollEvent_t *fEv = 0;
    Double_t fRate = 0.;
    T->SetBranchAddress("hit", & fHit);
    T->SetBranchAddress("ev", & fEv);
    T->SetBranchAddress("rate", & fRate);
    
    Double_t energy, hitr, rate;
    Int_t detector, pid;

    for(int ientry=0;ientry<nentry;ientry++){
       if(ientry%(nentry/10)==0)
         cout<<"analyzed "<<ientry<<" events!!"<<endl;
       T->GetEntry(ientry);
       for(size_t pk=0;pk<fHit->size();pk++){
          pid = (Int_t)TMath::Abs(fHit->at(pk).pid);
          detector = fHit->at(pk).det;
          energy = fHit->at(pk).e;
          hitr = fHit->at(pk).r;
          rate = fRate/1.e9/20.0;
          if(detector==28 && energy>1 && pid==11){
            h_radial[ifiles]->Fill(hitr,rate);
          }
       }
    }
    }
    TCanvas* c1 = new TCanvas("c1","c1",900,700);
    h_radial[0]->Draw("hist");
    h_radial[1]->Draw("hist sames");
    h_radial[2]->Draw("hist sames");
    gPad->Update();
    TPaveStats* st[nfiles];
    for(int ifiles=0;ifiles<nfiles;ifiles++){
       st[ifiles] = (TPaveStats*)h_radial[ifiles]->FindObject("stats");
       st[ifiles]->SetTextColor(color[ifiles]);
       st[ifiles]->SetY1NDC(0.9-0.15*ifiles);
       st[ifiles]->SetY2NDC(0.75-0.15*ifiles);
       st[ifiles]->SetX1NDC(0.9);
       st[ifiles]->SetX1NDC(0.75);
    }
    gPad->Modified();
    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.04);
    latex.SetTextColor(1);
    latex.DrawLatex(0.15,0.85,"thcommin=0.2 deg");
    latex.SetTextColor(2);
    latex.DrawLatex(0.15,0.80,"thcommin=0.02 deg");
    latex.SetTextColor(4);
    latex.DrawLatex(0.15,0.75,"thcommin=0.002 deg");
    c1->SaveAs(Form("./plots/thcommin_compare_epel_ring5.pdf"));
}
