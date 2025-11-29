// PlotClusterPhase.C
// Usage in ROOT:
// root -l
// .L PlotClusterPhase.C+
// PlotClusterPhase();

#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TString.h>
#include <iostream>

void PlotClusterPhase2()
{
  // ========= input/output =========
 // std::string inFileName = "/sphenix/tg/tg01/hf/mitrankova/ZigZagSimulation/LayerZigZag/output_phase/LayerZigZagG4sPHENIX_truth_*_phase.root";
  std::string inFileName = "/sphenix/tg/tg01/hf/mitrankova/ZigZagSimulation/LayerZigZag/output_phase/Pi_G4sPHENIX_10evts_*_phase.root";
  std::string outdir      = "/sphenix/user/mitrankova/ClusterPhase/analyse/plots/";
  std::string outfilename = "Sim_Fit_Pi_cluster_phase";

  const char* treeName = "phase_tree";

  // ========= cuts (same logic as your phase macro; pt_min=3.0 here) =========
  const float pt_min      = 4.0;
  const float pt_max      = 1e9;
  const int   quality_min = -1e9;
  const int   quality_max = 100;
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

  // ========= open file & tree (TChain) =========
  TChain* tree = new TChain(treeName, "tree");
  tree->Add(inFileName.c_str());

  // ========= branch variables =========
  float pt     = 0.f;
  int   quality = 0;
  int   ntpc    = 0;

  float cluster_phase_phi   = 0.f;
  float cluster_phase_time  = 0.f;
  int   cluster_layer       = -1;
  int   cluster_size_phi    = 0;
  int   cluster_size_time   = 0;
  int cluster_adc         = 0.f;
  std::vector<int> *hitside;

  tree->SetBranchAddress("pt",              &pt);
  tree->SetBranchAddress("quality",         &quality);
  tree->SetBranchAddress("ntpc",            &ntpc);

  tree->SetBranchAddress("cluster_phase_phi",  &cluster_phase_phi);
  tree->SetBranchAddress("cluster_phase_time", &cluster_phase_time);
  tree->SetBranchAddress("cluster_layer",      &cluster_layer);
  tree->SetBranchAddress("cluster_size_phi",   &cluster_size_phi);
  tree->SetBranchAddress("cluster_size_time",  &cluster_size_time);
  tree->SetBranchAddress("cluster_adc",        &cluster_adc);
  tree->SetBranchAddress("hit_side",      &hitside);

  // ========= book histograms (same binning as your phase macro) =========
  double phase_bin_width      = 0.04;
  int    nbins_phase          = 120 - 1;
  double max_phase            = (phase_bin_width/2)*(nbins_phase/2);
  double min_phase            = -(phase_bin_width/2)*(nbins_phase/2);

  double phase_time_bin_width = 0.04;
  int    nbins_phase_time     = 120 - 1;
  double max_phase_time       = (phase_time_bin_width/2)*(nbins_phase_time/2);
  double min_phase_time       = -(phase_time_bin_width/2)*(nbins_phase_time/2);

  double max_nhits       = 5,  min_nhits       = 1;
  double max_nhits_time  = 7,  min_nhits_time  = 1;
  double max_adc         = 500, min_adc        = 0;

  int nbins_nhits      = (max_nhits      - min_nhits);
  int nbins_nhits_time = (max_nhits_time - min_nhits_time);
  int nbins_adc        = ((max_adc - min_adc)/10.);

  TH2F* h_size_phase_phi[3];
  TH2F* h_size_phase_time[3];
  TH2F* h_adc_phase_phi[3];
  TH2F* h_adc_phase_time[3];
  TH1D* phase_phi_hist[3][5];
  TH1D* phase_phi_hist_full[3];
  TH1D* phase_time_hist[3][(int)max_nhits_time];
  TH1D* phase_time_hist_full[3];

  for(int iregion = 0; iregion < 3; ++iregion)
  {
     h_size_phase_phi[iregion]  = new TH2F(Form("h_size_phase_phi_R%d", iregion+1),
    ";Phase in #phi;Cluster size in #phi",
    nbins_phase,     min_phase,     max_phase,
    nbins_nhits,     min_nhits,     max_nhits);
    
    h_size_phase_time[iregion] = new TH2F(Form("h_size_phase_time_R%d", iregion+1),
    ";Phase in time;Cluster size in time",
    nbins_phase_time, min_phase_time, max_phase_time,
    nbins_nhits_time, min_nhits_time, max_nhits_time);
    h_adc_phase_phi[iregion] = new TH2F(Form("h_adc_phase_phi_R%d", iregion+1),
    ";Phase;ADC",
    nbins_phase, min_phase, max_phase,
    nbins_adc,   min_adc,   max_adc);
    h_adc_phase_time[iregion] = new TH2F(Form("h_adc_phase_time_R%d", iregion+1),
    ";Phase;ADC",
    nbins_phase_time, min_phase_time, max_phase_time,
    nbins_adc,        min_adc,        max_adc);

    h_size_phase_phi[iregion]->GetYaxis()->SetTitleOffset(1.5);
    h_size_phase_time[iregion]->GetYaxis()->SetTitleOffset(1.5);
    h_adc_phase_phi[iregion]->GetYaxis()->SetTitleOffset(1.5);
    h_adc_phase_time[iregion]->GetYaxis()->SetTitleOffset(1.5);
    
  }

  int colors[] = {kRed, kBlue, kGreen+2, kMagenta, kOrange+7};

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
        cout<<"hitside = "<<  hitside->at(0)<<" "<<hitside->back()<<endl;
    if ((hitside->at(0)!=1) ||( hitside->back()!=1)) continue;

    // restrict to R2 region (same as cut_*_r2 in your Draw-based macro)
    int region = (int)((cluster_layer-7)/16);

    // size vs phase in phi (require cluster_size_phi > 1)
    if (cluster_size_phi > 1)
    {
      h_size_phase_phi[region]->Fill(cluster_phase_phi, cluster_size_phi);
    }

    // size vs phase in time (require cluster_size_time > 1)
    if (cluster_size_time > 1)
    {
      h_size_phase_time[region]->Fill(cluster_phase_time, cluster_size_time);
    }

    // ADC vs phase (no size cut)
    h_adc_phase_phi[region]->Fill(cluster_phase_phi,  cluster_adc);
    h_adc_phase_time[region]->Fill(cluster_phase_time, cluster_adc);
  }

  // ========= Canvas 1: cluster size vs phase + projections =========

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
  h_size_phase_phi[0]->Draw("COLZ");
  latex_module[0].DrawLatex(0.4, 0.91, txt_module[0]);

  c->cd(2);
  h_size_phase_phi[1]->Draw("COLZ");
  latex_module[1].DrawLatex(0.4, 0.91, txt_module[1]);
  latex.DrawLatex(0.7, 0.96, txt_pt);
  latex.DrawLatex(0.05, 0.96,"#pi^{#minus} Simulation");

  c->cd(3);
  h_size_phase_phi[2]->Draw("COLZ");
  latex_module[2].DrawLatex(0.4, 0.91, txt_module[2]);
  

  c->SaveAs((outdir + outfilename + "_Phase.png").c_str());

  TCanvas* c1 = new TCanvas("c1", "Cluster Phase Analysis", 1500, 500);
  c1->Divide(3, 1);
  
  for(int iregion=0; iregion<3; ++iregion)
  {

      c1->cd(iregion+1);
      gPad->SetRightMargin(0.01);
      auto leg_phi = new TLegend(0.5,0.79,0.98,0.88);
      leg_phi->SetNColumns(2);
      leg_phi->SetBorderSize(0);
      leg_phi->SetTextSize(0.03);        // try 0.025–0.03
      //leg_phi->SetEntrySeparation(0.1);  // smaller distance between rows
      //leg_phi->SetMargin(0.01);           // less horizontal padding before markers
      phase_phi_hist_full[iregion] =
      h_size_phase_phi[iregion]->ProjectionX("phase_phi_hist_size_full", 1+1, nbins_nhits);
      phase_phi_hist_full[iregion]->SetLineColor(kBlack);
      phase_phi_hist_full[iregion]->SetLineWidth(2);
      phase_phi_hist_full[iregion]->Draw("");
      leg_phi->AddEntry(phase_phi_hist_full[iregion], "All clusters", "l");

      
      for (int i = 0; i < (int)max_nhits; ++i)
      {
        int ybin = h_size_phase_phi[iregion]->GetYaxis()->FindBin(i+2); // sizes 2,3,4,5
        phase_phi_hist[iregion][i] =
          h_size_phase_phi[iregion]->ProjectionX(Form("phase_phi_hist_size_%d", i+2),
                                        ybin, ybin);
        phase_phi_hist[iregion][i]->SetLineColor(colors[i]);
        phase_phi_hist[iregion][i]->GetYaxis()->SetTitle("Counts");
        phase_phi_hist[iregion][i]->GetYaxis()->SetTitleOffset(0.5);
        phase_phi_hist[iregion][i]->Draw("same");
        leg_phi->AddEntry(phase_phi_hist[iregion][i],
                          Form("%d pads", i+2), "l");
      }
      leg_phi->Draw();
      latex_module[iregion].DrawLatex(0.4, 0.91, txt_module[iregion]);
      if(iregion==1) {
        latex.DrawLatex(0.7, 0.96, txt_pt);
        latex.DrawLatex(0.05, 0.96,"#pi^{#minus} Simulation");
      }

  }
  c1->Update();
  c1->SaveAs((outdir + outfilename + "_dNdPhase_phi.png").c_str());


  TCanvas* c2 = new TCanvas("c2", "Cluster Phase Analysis", 1500, 500);
  c2->Divide(3, 1);
  for(int iregion=0; iregion<3; ++iregion)
  {
    c2->cd(iregion+1);
    gPad->SetRightMargin(0.01);
    auto leg_time = new TLegend(0.5,0.79,0.98,0.88);
    leg_time->SetBorderSize(0);
    leg_time->SetNColumns(2);
    leg_time->SetTextSize(0.03);     
      phase_time_hist_full[iregion] =
      h_size_phase_time[iregion]->ProjectionX("phase_time_hist_size_full", 1+1, nbins_nhits_time);
      phase_time_hist_full[iregion]->SetLineColor(kBlack);
      phase_time_hist_full[iregion]->SetLineWidth(2);
      phase_time_hist_full[iregion]->Draw("");
      leg_time->AddEntry(phase_time_hist_full[iregion], "All clusters", "l");
  for (int i = 0; i < max_nhits_time; i++){
    int ybin = h_size_phase_time[iregion]->GetYaxis()->FindBin(i+2);
    phase_time_hist[iregion][i] = h_size_phase_time[iregion]->ProjectionX(Form("phase_time_hist_size_%d", i+2), ybin, ybin);
    phase_time_hist[iregion][i]->SetLineColor(colors[i]);
    phase_time_hist[iregion][i]->GetYaxis()->SetTitle("Counts");
    phase_time_hist[iregion][i]->GetYaxis()->SetTitleOffset(1.5);
    phase_time_hist[iregion][i]->Draw("same");
    leg_time->AddEntry(phase_time_hist[iregion][i], Form("%d timebins", i+2), "l");
  }
  leg_time->Draw();
  latex_module[iregion].DrawLatex(0.4, 0.91, txt_module[iregion]);
  if(iregion==1) {
        latex.DrawLatex(0.7, 0.96, txt_pt);
        latex.DrawLatex(0.05, 0.96,"#pi^{#minus} Simulation");
      }
  }
  c2->Update();
  c2->SaveAs((outdir + outfilename + "_dNdPhase_time.png").c_str());


  // ========= Print statistics =========
  std::cout << "\n========== Phase Distribution Statistics ==========" << std::endl;
  std::cout << "Total clusters (before cuts): " << nEntries << std::endl;
  std::cout << "pt in [" << pt_min << ", " << pt_max << "]" << std::endl;
  std::cout << "quality in [" << quality_min << ", " << quality_max << "]" << std::endl;
  std::cout << "ntpc in [" << ntpc_min << ", " << ntpc_max << "]" << std::endl;
  std::cout << "R2 layer range: [" << r2_min << ", " << r2_max << ")" << std::endl;
  std::cout << "===================================================\n" << std::endl;

  // ========= Save histograms to ROOT file =========
  TFile* outFile =
    new TFile((outdir + outfilename + ".root").c_str(), "RECREATE");
  for(int iregion=0; iregion<3; ++iregion)
  {
    h_size_phase_phi[iregion]->Write();
    h_size_phase_time[iregion]->Write();
    h_adc_phase_phi[iregion]->Write();
    h_adc_phase_time[iregion]->Write();

    phase_phi_hist_full[iregion]->Write(Form("phase_phi_prof_region_%d_npads_all", iregion)); 
    for (int i = 0; i < (int)max_nhits; ++i)
    {
      if (phase_phi_hist[iregion][i]) phase_phi_hist[iregion][i]->Write(Form("phase_phi_prof_region_%d_npads_%d", iregion, i+2));
    }
    phase_time_hist_full[iregion]->Write(Form("phase_time_prof_region_%d_npads_all", iregion));
    for (int i = 0; i < (int)max_nhits_time; ++i)
    {
      if (phase_time_hist[iregion][i]) phase_time_hist[iregion][i]->Write(Form("phase_time_prof_region_%d_npads_%d", iregion, i+2));
    }
  }

  c->Write();
  c1->Write();
  c2->Write();

  outFile->Close();

  std::cout << "Saved PNGs and ROOT histograms." << std::endl;
}
