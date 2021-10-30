void skimTree(int fileNo){
    const string fileDir = "/volatile/halla/parity/adhidevi/remoll_rootfiles/PhotonBlocker/sim_ee_lam/";
    TFile* oldfile = new TFile(Form("%ssim_ee_lam_%d.root",fileDir.c_str(),fileNo),"READ");
    TTree* oldTree = (TTree*)oldfile->Get("T");
    oldTree->SetBranchStatus("*",0);
    oldTree->SetBranchStatus("hit.pid",1);
    oldTree->SetBranchStatus("hit.r",1);
    oldTree->SetBranchStatus("hit.e",1);
    oldTree->SetBranchStatus("hit.k",1);
    oldTree->SetBranchStatus("rate",1);
    oldTree->SetBranchStatus("hit.det",1);
    TFile* newfile = new TFile(Form("%sskim/skim_sim_ee_lam_%d.root",fileDir.c_str(),fileNo),"recreate");
    TTree* newTree = oldTree->CopyTree("");
    newTree->Print();
    newTree->Write();
    oldTree->SetBranchStatus("*",1);
    delete oldTree;
    delete newfile;
}
