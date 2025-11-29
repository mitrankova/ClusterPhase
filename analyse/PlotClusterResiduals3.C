// PlotClusterResiduals.C
// Usage in ROOT:
// root -l
// .L PlotClusterResiduals.C+
// PlotClusterResiduals();

#include <TFile.h>
#include <TTree.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TString.h>
#include <iostream>

void PlotClusterResiduals3()
{
  std::string inFileName = "/sphenix/tg/tg01/hf/mitrankova/ZigZagSimulation/LayerZigZag/output_phase/Pi_G4sPHENIX_10evts_*_phase.root";
//std::string inFileName = "/sphenix/user/mitrankova/ClusterPhase/output_phase/ZigZagG4sPHENIX10evts_0_phase.root";
//snprintf(infileresid1, sizeof(infileresid1), "/sphenix/tg/tg01/hf/mitrankova/AuAu_2025_ClusterPhase/75396/Au_Au_Resid_200evts_0skip_75396-*_phase.root");

  std::string outdir      = "/sphenix/user/mitrankova/ClusterPhase/analyse/plots/";
  std::string outfilename = "Sim_Truth_Pi_cluster_residual";

  const char* treeName = "phase_tree";
  // ========= cuts =========
  const float pt_min      = 4.0;
  const float pt_max      = 1e9;
  const int   quality_min = -1e9;
  const int   quality_max = 50;
  const int   ntpc_min    = 40;
  const int   ntpc_max    = 1e9;

  // layer regions: [7,7+16), [7+16,7+2*16), [7+2*16,7+3*16)
  const int startLayer  = 7;
  const int width       = 16;

  const int r1_min = startLayer;
  const int r1_max = startLayer + width;

  const int r2_min = startLayer + width;
  const int r2_max = startLayer + 2*width;

  const int r3_min = startLayer + 2*width;
  const int r3_max = startLayer + 3*width;

  // ========= open file & tree =========
  /*TFile* inFile = TFile::Open(inFileName.c_str(), "READ");
  if (!inFile || inFile->IsZombie())
  {
    std::cerr << "Error opening file " << inFileName << std::endl;
    return;
  }

  TTree* tree = dynamic_cast<TTree*>(inFile->Get(treeName));
  if (!tree)
  {
    std::cerr << "Error: tree " << treeName << " not found." << std::endl;
    inFile->Close();
    return;
  }
  */


  TChain *tree = new TChain(treeName,"tree");
  tree->Add(inFileName.c_str());


  // ========= branch variables =========
  float pt = 0.f;
  int   quality = 0;
  int   ntpc = 0;

  float cluster_phase_phi = 0.f;
  float cluster_residual_rphi = 0.f;
  int   cluster_layer = -1;
  int   cluster_size_phi = 0;
  std::vector<int> *hitside;

  tree->SetBranchAddress("pt", &pt);
  tree->SetBranchAddress("quality", &quality);
  tree->SetBranchAddress("ntpc", &ntpc);

  tree->SetBranchAddress("cluster_phase_phi",     &cluster_phase_phi);
  tree->SetBranchAddress("sim_cluster_residual_rphi", &cluster_residual_rphi);
  tree->SetBranchAddress("cluster_layer",         &cluster_layer);
  tree->SetBranchAddress("cluster_size_phi",      &cluster_size_phi);
  tree->SetBranchAddress("hit_side",      &hitside);

  // ========= book histograms =========
  int   nBinsPhi   = 120;
  float phiMin     = -1.2;
  float phiMax     =  1.2;

  int   nBinsRes   = 200;
  float resMin     = -0.2;//-0.06;
  float resMax     =  0.2;//0.06;

  

  TH2F* h_r1 = new TH2F("h_r1",
                        ";cluster phase;cluster residual r#phi [cm]",
                        nBinsPhi, phiMin, phiMax,
                        nBinsRes, resMin, resMax);
  h_r1->GetYaxis()->SetTitleOffset(1.5);

  TH2F* h_r2 = new TH2F("h_r2",
                        ";cluster phase;cluster residual r#phi [cm]",
                        nBinsPhi, phiMin, phiMax,
                        nBinsRes, resMin, resMax);
  h_r2->GetYaxis()->SetTitleOffset(1.5);
  TH2F* h_r3 = new TH2F("h_r3",
                        ";cluster phase;cluster residual r#phi [cm]",
                        nBinsPhi, phiMin, phiMax,
                        nBinsRes, resMin, resMax);
  h_r3->GetYaxis()->SetTitleOffset(1.5);
  // ========= loop over entries =========
  Long64_t nEntries = tree->GetEntries();
  std::cout << "Entries: " << nEntries << std::endl;

  for (Long64_t i = 0; i < nEntries; ++i)
  {
    tree->GetEntry(i);

    // event-level cuts
    if (pt < pt_min || pt > pt_max) continue;
    if (quality < quality_min || quality > quality_max) continue;
    if (ntpc < ntpc_min || ntpc > ntpc_max) continue;
    if (cluster_size_phi <= 1) continue;
    cout<<"hitside = "<<  hitside->at(0)<<" "<<hitside->back()<<endl;
    if ((hitside->at(0)!=1) ||( hitside->back()!=1)) continue;

    // choose region by cluster_layer
    if (cluster_layer >= r1_min && cluster_layer < r1_max)
    {
      h_r1->Fill(cluster_phase_phi, cluster_residual_rphi);
    }
    else if (cluster_layer >= r2_min && cluster_layer < r2_max)
    {
      h_r2->Fill(cluster_phase_phi, cluster_residual_rphi);
    }
    else if (cluster_layer >= r3_min && cluster_layer < r3_max)
    {
      h_r3->Fill(cluster_phase_phi, cluster_residual_rphi);
    }
  }

  // ========= build mean profiles (mean residual vs phi) =========
  TProfile* p_r1 = h_r1->ProfileX("p_r1");
  TProfile* p_r2 = h_r2->ProfileX("p_r2");
  TProfile* p_r3 = h_r3->ProfileX("p_r3");

  // cosmetics for profiles
  p_r1->SetLineColor(kBlack);
  p_r1->SetMarkerColor(kBlack);
  p_r1->SetMarkerStyle(20);
  p_r1->SetMarkerSize(0.3);

  p_r2->SetLineColor(kBlack);
  p_r2->SetMarkerColor(kBlack);
  p_r2->SetMarkerStyle(20);
  p_r2->SetMarkerSize(0.3);

  p_r3->SetLineColor(kBlack);
  p_r3->SetMarkerColor(kBlack);
  p_r3->SetMarkerStyle(20);
  p_r3->SetMarkerSize(0.3);


  p_r1->SetErrorOption("s");  // "s" = standard deviation (RMS)
  p_r2->SetErrorOption("s");
  p_r3->SetErrorOption("s");


  // ========= draw & save PNG =========


    TLatex latex;
  latex.SetNDC(true);
  latex.SetTextSize(0.05);
  latex.SetTextFont(42);
  latex.SetTextColor(kRed+2);

  TString txt_pt;
  txt_pt.Form("#it{p_{T} > %.1f GeV}", pt_min);

  TLatex latex_module[3];
  TString txt_module[3];
  for(int iregion = 0; iregion < 3; ++iregion)
  {
  latex_module[iregion].SetNDC(true);
  latex_module[iregion].SetTextSize(0.05);
  latex_module[iregion].SetTextFont(42);
  latex_module[iregion].SetTextColor(kBlack);

  
  txt_module[iregion].Form("R%d region", iregion+1);
  }


  gStyle->SetOptStat(0);
  TCanvas* c = new TCanvas("c_resid_vs_phase",
                           "cluster_residual_{r#phi} vs cluster_phase_{#phi}",
                           1500, 500);
  c->Divide(3,1);

  c->cd(1);
  h_r1->Draw("COLZ");
  p_r1->Draw("same");
  latex_module[0].DrawLatex(0.4, 0.91, txt_module[0]);

  c->cd(2);
  h_r2->Draw("COLZ");
  p_r2->Draw("same");
  latex_module[1].DrawLatex(0.4, 0.91, txt_module[1]);
  latex.DrawLatex(0.75, 0.96, txt_pt);
  latex.DrawLatex(0.00, 0.96,"#pi^{#minus} Simulation");
  latex.DrawLatex(0.3, 0.96,"Using truth positions");

  c->cd(3);
  h_r3->Draw("COLZ");
  p_r3->Draw("same");
  latex_module[2].DrawLatex(0.4, 0.91, txt_module[2]);

  c->SaveAs((outdir + outfilename + ".png").c_str());

  // ========= save hists & profiles to ROOT file =========
  TFile* outFile =
      new TFile((outdir + outfilename + ".root").c_str(),
                "RECREATE");
  c->Write();
  h_r1->Write();
  h_r2->Write();
  h_r3->Write();
  p_r1->Write();
  p_r2->Write();
  p_r3->Write();
  
  outFile->Close();

  //inFile->Close();

  std::cout << "Saved PNG and ROOT histograms/profiles." << std::endl;
}
