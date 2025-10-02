#ifndef CLUSTERPHASEANALYSIS_H
#define CLUSTERPHASEANALYSIS_H

#include <fun4all/SubsysReco.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase/TpcDefs.h>
#include <trackbase/ActsGeometry.h>

#include <string>
#include <vector>

class PHCompositeNode;
class TFile;
class TTree;
class ActsGeometry;
class TrkrClusterContainer;
class TrkrHitSetContainer;
class TrkrClusterHitAssoc;
class SvtxTrack;
class SvtxTrackMap;


class ClusterPhaseAnalysis : public SubsysReco
{
 public:
  ClusterPhaseAnalysis(const std::string &name = "ClusterPhaseAnalysis",
                       const std::string &filename = "cluster_phase_output.root");

  ~ClusterPhaseAnalysis() override;

  int Init(PHCompositeNode *topNode) override;
  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;

  void set_cluster_container_name(const std::string &name) { m_clusterContainerName = name; }
  void setMinPt(double value){ m_minPt = value;}
  void setMaxQuality(double value){ m_maxQuality = value;}
  void setRequiredTpcClusters(int value){ m_minTpcClusters = value;}

 private:
  void reset_tree_vars();
  void processTracks();
  bool checkTrack(SvtxTrack *track);
  void FillClusters(SvtxTrack *track);
  int  getNodes(PHCompositeNode *topNode);
  //int get_phase_from_hit(TrkrDefs::hitkey hitkey, TrkrDefs::hitsetkey hitsetkey, TrkrHitSetContainer *hitmap);

    SvtxTrackMap *trackmap = nullptr;
    ActsGeometry *geometry = nullptr;
    TrkrClusterContainer *clustermap = nullptr;

  // Output file and tree
  std::string m_outputFileName;
  TFile *m_outputFile;
  TTree *m_tree;

  // Node names
  std::string m_clusterContainerName = "TRKR_CLUSTER";
  std::string m_trackMapName = "SvtxTrackMap";

  // Event-level variables
  int m_event;
  std::vector<int> m_firedTriggers;
  uint64_t m_gl1BunchCrossing = std::numeric_limits<uint64_t>::quiet_NaN();
  uint64_t m_bco = std::numeric_limits<uint64_t>::quiet_NaN();
  uint64_t m_bcotr = std::numeric_limits<uint64_t>::quiet_NaN();

  // Cluster-level variables (vectors to store multiple clusters per event)
  float m_cluster_x = std::numeric_limits<float>::quiet_NaN();
  float m_cluster_y = std::numeric_limits<float>::quiet_NaN();
  float m_cluster_z = std::numeric_limits<float>::quiet_NaN();
  float m_cluster_r = std::numeric_limits<float>::quiet_NaN();
  float m_cluster_phi = std::numeric_limits<float>::quiet_NaN();
  float m_cluster_eta = std::numeric_limits<float>::quiet_NaN();
  float m_cluster_e = std::numeric_limits<float>::quiet_NaN();
  float m_cluster_adc = std::numeric_limits<float>::quiet_NaN();
  int m_cluster_max_adc = std::numeric_limits<int>::quiet_NaN();
  unsigned int m_cluster_layer = std::numeric_limits<unsigned int>::quiet_NaN();
  unsigned int m_cluster_size = std::numeric_limits<unsigned int>::quiet_NaN();
  float m_cluster_phase = std::numeric_limits<float>::quiet_NaN();

  
  // Truth cluster information
  std::vector<float> m_truth_cluster_x;
  std::vector<float> m_truth_cluster_y;
  std::vector<float> m_truth_cluster_z;
  std::vector<float> m_truth_cluster_r;
  std::vector<float> m_truth_cluster_phi;

  // Phase information for each cluster

  
  // Hit-level information (nested vectors: outer = cluster index, inner = hits in that cluster)
  std::vector<std::vector<float>> m_hit_x;
  std::vector<std::vector<float>> m_hit_y;
  std::vector<std::vector<float>> m_hit_z;
  std::vector<std::vector<int>> m_hit_adc;
  std::vector<std::vector<int>> m_hit_phase;
  std::vector<std::vector<unsigned int>> m_hit_layer;

  // Phase histogram data (for dN/dPhase calculation)
  std::vector<int> m_phase_counts;  // Histogram bins for phases
  static const int m_n_phase_bins = 21;  // -10 to +10
  static const int m_phase_min = -10;
  static const int m_phase_max = 10;
  
  double m_minPt = 0.5;
  double m_maxQuality = 100.0;
  double m_maxResidual = 1.0;
  int m_minTpcClusters = 30;
  bool m_ifzerofield = false;



};

#endif  // CLUSTERPHASEANALYSIS_H