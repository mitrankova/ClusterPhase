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
class TrkrHitSetContainer;
//class PHG4TpcCylinderGeomContainer;
class PHG4TpcGeomContainer;
class PHG4CylinderGeomContainer;
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
  void setIsSimulation(int value){ isSimulation = value;}

 private:
  void reset_tree_vars();
  void processTracks();
  bool checkTrack(SvtxTrack *track);
  void FillClusters(SvtxTrack *track);
  void ComputeFitTruthAtLayer(SvtxTrack *track);
  void FindTruthClusters(uint64_t key, Acts::Vector3 glob);
  void FillHits(uint64_t ckey);
  void CalculatePhase(uint64_t ckey);
  int  getNodes(PHCompositeNode *topNode);
  //int get_phase_from_hit(TrkrDefs::hitkey hitkey, TrkrDefs::hitsetkey hitsetkey, TrkrHitSetContainer *hitmap);

    SvtxTrackMap *trackmap = nullptr;
    ActsGeometry *geometry = nullptr;
    TrkrClusterContainer *clustermap = nullptr;
    TrkrClusterContainer *truthClustersmap = nullptr;
    TrkrHitSetContainer *hitmap = nullptr;
    TrkrClusterHitAssoc *clusterhitassocmap = nullptr;
    //PHG4TpcCylinderGeomContainer *tpcGeom = nullptr;
    PHG4TpcGeomContainer *tpcGeom = nullptr;

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


  float m_vdrift = std::numeric_limits<float>::quiet_NaN();

  // Cluster-level variables
  float m_cluster_x = std::numeric_limits<float>::quiet_NaN();
  float m_cluster_y = std::numeric_limits<float>::quiet_NaN();
  float m_cluster_z = std::numeric_limits<float>::quiet_NaN();
  float m_cluster_time = std::numeric_limits<float>::quiet_NaN();
  float m_cluster_r = std::numeric_limits<float>::quiet_NaN();
  float m_cluster_phi = std::numeric_limits<float>::quiet_NaN();
  float m_cluster_eta = std::numeric_limits<float>::quiet_NaN();
  float m_cluster_e = std::numeric_limits<float>::quiet_NaN();
  int m_cluster_adc = std::numeric_limits<int>::quiet_NaN();
  int m_cluster_max_adc = std::numeric_limits<int>::quiet_NaN();
  int m_cluster_layer = std::numeric_limits<int>::quiet_NaN();
  int m_cluster_size = std::numeric_limits<int>::quiet_NaN();
  int m_cluster_size_phi = std::numeric_limits<int>::quiet_NaN();
  int m_cluster_size_time = std::numeric_limits<int>::quiet_NaN();

  // Phase-level information

  float m_cluster_phase_phi = std::numeric_limits<float>::quiet_NaN();
  float m_cluster_phase_time = std::numeric_limits<float>::quiet_NaN();

  // Truth cluster information
  float m_truth_cluster_x = std::numeric_limits<float>::quiet_NaN();
  float m_truth_cluster_y = std::numeric_limits<float>::quiet_NaN();
  float m_truth_cluster_z = std::numeric_limits<float>::quiet_NaN();
  float m_truth_cluster_r = std::numeric_limits<float>::quiet_NaN();
  float m_truth_cluster_phi = std::numeric_limits<float>::quiet_NaN();
  float m_truth_cluster_time = std::numeric_limits<float>::quiet_NaN();

  float m_cluster_residual_rphi = std::numeric_limits<float>::quiet_NaN();
  float m_cluster_residual_time = std::numeric_limits<float>::quiet_NaN();
  std::map<TrkrDefs::cluskey, float> residual_rphi_map;
  std::map<TrkrDefs::cluskey, float> residual_time_map;
    std::map<TrkrDefs::cluskey, float> m_truth_cluster_x_map;
  std::map<TrkrDefs::cluskey, float> m_truth_cluster_y_map;
  std::map<TrkrDefs::cluskey, float> m_truth_cluster_z_map;

  // Hit-level information ( vectors:   hits in that cluster)
  std::vector<uint32_t> m_hitkeys;
  std::vector<int> m_hit_pad;
  std::vector<int> m_hit_tbin;
  std::vector<int> m_hit_layer;
  std::vector<int> m_hit_region;
  std::vector<int> m_hit_side;
  std::vector<int> m_hit_adc;
  std::vector<float> m_hit_x;
  std::vector<float> m_hit_y;
  std::vector<float> m_hit_z;
  std::vector<float> m_hit_r;
  std::vector<float> m_hit_phi;
  std::vector<float> m_hit_time;




  
  double m_minPt = 0.5;
  double m_maxQuality = 100.0;
  double m_maxResidual = 1.0;
  int m_minTpcClusters = 30;
  bool m_ifzerofield = false;
  bool isSimulation = true;



};

#endif  // CLUSTERPHASEANALYSIS_H