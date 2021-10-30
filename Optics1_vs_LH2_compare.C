#include "remolltypes.hh"
void Optics1_vs_LH2_compare(){
   gStyle->SetTitleYOffset(1.3);
   gStyle->SetOptStat(0);
   TString rootfile_dir = "/volatile/halla/parity/adhidevi/remoll_rootfiles/default-geo";
   TString file[] = {"LH2_beam","Optics1_beam"};
   const int nfile = sizeof(file)/sizeof(*file);
   TString sim[nfile] = {"LH2","Optics1"};
   TH1D* hist[nfile];
   TH1D* histR5[nfile];
   const int color[nfile] = {1,2};
   const int colorR5[nfile] = {4,6};
   const double rmin = 920;
   const double rmax = 1060;
   
   const int xmin = 0;
   const int xmax = 1500;
   const int nbins = 200;
   Long64_t nentry;

   for(int ifile=0;ifile<nfile;ifile++){
      hist[ifile] = new TH1D(Form("%s",sim[ifile].Data()),Form("Radial Distribution at Ring 5;Radius (mm); hits/%sthrownEvents)","#"),nbins,xmin,xmax);
      histR5[ifile] = new TH1D(Form("%s_R5",sim[ifile].Data()),Form("Radial Distribution at Ring 5;Radius (mm); hits/%sthrownEvents)","#"),nbins,xmin,xmax);
      hist[ifile]->SetLineColor(color[ifile]);
      histR5[ifile]->SetLineColor(colorR5[ifile]);
      hist[ifile]->Sumw2();
      histR5[ifile]->Sumw2();
      
      TChain* T = new TChain("T");
      int nsplit = 0;
      for(int isplit=1001;isplit<=1927;isplit++){
         nsplit++;
         T->Add(Form("%s/%s/%s_%d.root",rootfile_dir.Data(),file[ifile].Data(),file[ifile].Data(),isplit));
      }
      cout<<Form("Found %d file splits!!!",nsplit)<<endl;
      nentry = T->GetEntries();
      std::vector<remollGenericDetectorHit_t> *fHit = 0;
      remollEvent_t *fEv = 0;
      Double_t fRate = 0;

      T->SetBranchAddress("hit", &fHit);

      Double_t energy, hitr, rate, beamp, kinE;
      Int_t detector, pid;
      for(int ientry=0;ientry<nentry;ientry++){
         if(ientry%(nentry/10)==0)
            cout<<Form("Analyzed %d events!!!",ientry)<<endl;
         T->GetEntry(ientry);
         for(size_t pk=0;pk<fHit->size();pk++){
            pid = (Int_t)fHit->at(pk).pid;
            detector = fHit->at(pk).det;
            energy = fHit->at(pk).e;
            kinE = fHit->at(pk).k;
            hitr = fHit->at(pk).r;
            rate = 1;
            if(detector==28 && kinE>1 && hitr>0 && pid==11){
              hist[ifile]->Fill(hitr,rate);
              if(hitr>=rmin && hitr<=rmax){
                histR5[ifile]->Fill(hitr,rate);
              }
            }
         }
      }
   }
   
   TLatex latex;
   latex.SetNDC();
   latex.SetTextSize(0.04);
   TCanvas* c_linear = new TCanvas("c_linear");
   hist[0]->SetTitle("Radial Distribution of Electrons at Ring 5");
   hist[0]->Scale(1./nentry);
   hist[1]->Scale(1./nentry);
   histR5[0]->Scale(1./nentry);
   histR5[1]->Scale(1./nentry);
   hist[0]->Draw("hist");
   hist[0]->GetYaxis()->SetRangeUser(1.5e-8,1.0);
   hist[1]->Draw("hist sames");
   histR5[0]->Draw("hist sames");
   histR5[1]->Draw("hist sames");
   latex.SetTextColor(color[0]);
   latex.DrawLatex(0.40,0.85,Form("LH2"));
   latex.SetTextColor(color[1]);
   latex.DrawLatex(0.40,0.80,Form("Optics1"));
   latex.SetTextColor(colorR5[0]);
   latex.DrawLatex(0.50,0.85,Form("Ring 5 ingetral: %.3e",histR5[0]->Integral()));
   latex.SetTextColor(colorR5[1]);
   latex.DrawLatex(0.50,0.80,Form("Ring 5 ingetral: %.3e",histR5[1]->Integral()));
/*
   gPad->Update();
   TPaveStats* st1 = (TPaveStats*)hist[0]->FindObject("stats");
   TPaveStats* st2 = (TPaveStats*)hist[1]->FindObject("stats");
   st1->SetTextColor(color[0]);
   st2->SetTextColor(color[1]);
   st1->SetX1NDC(0.75);
   st1->SetX2NDC(0.90);
   st2->SetX1NDC(0.75);
   st2->SetX2NDC(0.90);
   st1->SetY1NDC(0.75);
   st1->SetY2NDC(0.90);
   st2->SetY1NDC(0.60);
   st2->SetY2NDC(0.75);
   st1->SetFillStyle(0);
   st2->SetFillStyle(0);
   gPad->Modified();
*/
   c_linear->SaveAs("./temp/ring5_linear.pdf");

   TCanvas* c_log = new TCanvas("c_log");
   gPad->SetLogy();
   hist[0]->Draw("hist");
   hist[1]->Draw("hist sames");
   histR5[0]->Draw("hist sames");
   histR5[1]->Draw("hist sames");
   latex.SetTextColor(color[0]);
   latex.DrawLatex(0.40,0.85,Form("LH2"));
   latex.SetTextColor(color[1]);
   latex.DrawLatex(0.40,0.80,Form("Optics1"));
   latex.SetTextColor(colorR5[0]);
   latex.DrawLatex(0.50,0.85,Form("Ring 5 ingetral: %.3e",histR5[0]->Integral()));
   latex.SetTextColor(colorR5[1]);
   latex.DrawLatex(0.50,0.80,Form("Ring 5 ingetral: %.3e",histR5[1]->Integral()));
   c_log->SaveAs("./temp/ring5_log.pdf");
   
   gSystem->Exec(Form("pdfunite ./temp/ring5_*.pdf ./plots/defaultGeo_LH2_Optics1_ring5.pdf"));
   gSystem->Exec(Form("rm -rf ./temp/ring5_*.pdf"));
}
