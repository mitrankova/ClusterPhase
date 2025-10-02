#include "ClusterPhaseAnalysis.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase/TpcDefs.h>
#include <trackbase/TrkrHit.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainer.h>
#include <trackbase/TrkrClusterHitAssoc.h>
#include <trackbase/ActsGeometry.h>  
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/TrackSeedContainer.h>

#include <TFile.h>
#include <TTree.h>

#include <cmath>
#include <iostream>
#include <map>


namespace
{

    // square
    template <class T>
    inline constexpr T square(const T& x)
    {
        return x * x;
    }

    // radius
    template <class T>
    inline T r(const T& x, const T& y)
    {
        return std::sqrt(square(x) + square(y));
    }

    // get cluster keys from a given track
    std::vector<TrkrDefs::cluskey> get_cluster_keys(SvtxTrack* track)
    {
        std::vector<TrkrDefs::cluskey> out;
        for (const auto& seed : {track->get_silicon_seed(), track->get_tpc_seed()})
        {
            if (seed)
            {
                std::copy(seed->begin_cluster_keys(), seed->end_cluster_keys(), std::back_inserter(out));
            }
        }
        return out;
    }


}  // namespace

//___________________________________________________________________________________

ClusterPhaseAnalysis::ClusterPhaseAnalysis(const std::string &name,
                                           const std::string &filename)
  : SubsysReco(name)
  , m_outputFileName(filename)
  , m_outputFile(nullptr)
  , m_tree(nullptr)
  , m_clusterContainerName("TRKR_CLUSTER")
  , m_event(0)
{
}
//___________________________________________________________________________________
ClusterPhaseAnalysis::~ClusterPhaseAnalysis()
{
  if (m_outputFile)
  {
    m_outputFile->Close();
    delete m_outputFile;
  }
}
//___________________________________________________________________________________
int ClusterPhaseAnalysis::Init(PHCompositeNode* /*topNode*/)
{
  std::cout << "ClusterPhaseAnalysis::Init - Creating output file: " << m_outputFileName << std::endl;

  m_outputFile = new TFile(m_outputFileName.c_str(), "RECREATE");
  if (!m_outputFile || m_outputFile->IsZombie())
  {
    std::cout << "ClusterPhaseAnalysis::Init - Error: Cannot create output file" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_tree = new TTree("phase_tree", "Cluster Phase Analysis Tree");

  // Event-level branches
  m_tree->Branch("event", &m_event, "event/I");

  // Cluster-level branches
  m_tree->Branch("cluster_x", &m_cluster_x, "m_cluster_x/F");
  m_tree->Branch("cluster_y", &m_cluster_y, "m_cluster_y/F");
  m_tree->Branch("cluster_z", &m_cluster_z, "m_cluster_z/F");
  m_tree->Branch("cluster_r", &m_cluster_r, "m_cluster_r/F");
  m_tree->Branch("cluster_phi", &m_cluster_phi, "m_cluster_phi/F");
  m_tree->Branch("cluster_eta", &m_cluster_eta, "m_cluster_eta/F");
  m_tree->Branch("cluster_adc", &m_cluster_adc, "m_cluster_adc/I");
  m_tree->Branch("cluster_max_adc", &m_cluster_max_adc, "m_cluster_max_adc/I");
  m_tree->Branch("cluster_layer", &m_cluster_layer, "m_cluster_layer/I");
  m_tree->Branch("cluster_size", &m_cluster_size, "m_cluster_size/I");
  m_tree->Branch("cluster_phase", &m_cluster_phase, "m_cluster_phase/F");

  // Truth cluster branches
  m_tree->Branch("truth_cluster_x", &m_truth_cluster_x);
  m_tree->Branch("truth_cluster_y", &m_truth_cluster_y);
  m_tree->Branch("truth_cluster_z", &m_truth_cluster_z);
  m_tree->Branch("truth_cluster_r", &m_truth_cluster_r);
  m_tree->Branch("truth_cluster_phi", &m_truth_cluster_phi);

  // Hit-level branches (nested vectors)
  m_tree->Branch("hit_x", &m_hit_x);
  m_tree->Branch("hit_y", &m_hit_y);
  m_tree->Branch("hit_z", &m_hit_z);
  m_tree->Branch("hit_adc", &m_hit_adc);
  m_tree->Branch("hit_phase", &m_hit_phase);
  m_tree->Branch("hit_layer", &m_hit_layer);

  // Phase histogram
  m_tree->Branch("phase_counts", &m_phase_counts);

  m_event = 0;

  std::cout << "ClusterPhaseAnalysis::Init - Initialization complete" << std::endl;

  return Fun4AllReturnCodes::EVENT_OK;
}
//___________________________________________________________________________________
int ClusterPhaseAnalysis::InitRun(PHCompositeNode* topNode)
{
  std::cout << "ClusterPhaseAnalysis::InitRun - Starting run" << std::endl;
     
  if (getNodes(topNode) != Fun4AllReturnCodes::EVENT_OK)
  {
      std::cout << "ClusterPhaseAnalysis::InitRun - Abort Event:: getNodes - Event not ok " << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}
//___________________________________________________________________________________
int ClusterPhaseAnalysis::process_event(PHCompositeNode* /*topNode*/)
{
  if (Verbosity())
  {
    std::cout << "ClusterPhaseAnalysis::process_event - Processing event " << m_event << std::endl;
  }

  processTracks();

  return Fun4AllReturnCodes::EVENT_OK;
}




//___________________________________________________________________________________
void ClusterPhaseAnalysis::processTracks()
{
    
    if (Verbosity())
    {
        std::cout << "ClusterPhaseAnalysis::processTracks - proto track size " << trackmap->size() << std::endl;
    }

    for (const auto& [trackKey, track] : *trackmap)
    {
        if (!track)
        {
            continue;
        }
        if (Verbosity() > 1)
        {
            std::cout << "ClusterPhaseAnalysis::processTracks - Processing track key " << trackKey << std::endl;
        }

        
        if (checkTrack(track))
        {

            FillClusters(track);

        }
    }

 //   return Fun4AllReturnCodes::EVENT_OK;
}

//___________________________________________________________________________________
bool ClusterPhaseAnalysis::checkTrack(SvtxTrack* track)
{

    bool ifGoodTrack = true;
    if (track->get_pt() < m_minPt)
    {
        std::cout<<" Track pt is lower then "<<m_minPt<<std::endl;
        ifGoodTrack = false;
        return false;
    }

    if (m_ifzerofield == false && (track->get_quality() > m_maxQuality))
    {
        std::cout<<" Track quality is higher then "<<m_maxQuality<<std::endl;
        ifGoodTrack = false;
        return false;
    }

    int m_ntpc=0;
    for (const auto& ckey : get_cluster_keys(track))
    {
        const auto detId = TrkrDefs::getTrkrId(ckey);
        if (detId == TrkrDefs::tpcId)
        {
            m_ntpc++;
        }else{
            ifGoodTrack = false;
            return false;
        }
    }

    if (m_ntpc<m_minTpcClusters) 
    {
        ifGoodTrack = false;
        return false;
    }

    if (ifGoodTrack && Verbosity() > 2)
    {
        std::cout << "ClusterPhaseAnalysis::checkTrack - pt: " << track->get_pt() <<" ntpc: "<<m_ntpc<<" qualiuty: "<<track->get_quality()<<" If zero field: "<<m_ifzerofield<< std::endl;
    }

    return true;
}
//___________________________________________________________________________________
void ClusterPhaseAnalysis::FillClusters(SvtxTrack* track)
{

    for (const auto& cluskey : get_cluster_keys(track))
    {
      auto *cluster = clustermap->findCluster(cluskey);
      std::cout<<"        TPC cluster with key: "<<cluskey<<std::endl;
      if (!cluster)
      {
        continue;
      }

      Acts::Vector3 glob;
      glob = geometry->getGlobalPosition(cluskey, cluster);
      
      float x = glob.x();
      float y = glob.y();
      float z = glob.z();
      float r = sqrt(x * x + y * y);
      float phi = atan2(y, x);
      float theta = atan2(r, z);
      float eta = -log(tan(theta / 2.0));
    
      int adc = cluster->getAdc();
      int max_adc = cluster->getMaxAdc();
    
      unsigned int layer = TrkrDefs::getLayer(cluskey);

      m_cluster_x = x;
      m_cluster_y = y;
      m_cluster_z = z;
      m_cluster_r = r;
      m_cluster_phi = phi;
      m_cluster_eta = eta;
      m_cluster_adc = adc;
      m_cluster_max_adc = max_adc;
      m_cluster_layer = layer;

      
      std::cout<<"FILL TREE "<<std::endl;
      m_tree->Fill();
     
  }
}
//_____________________________________________________________________________________________

int ClusterPhaseAnalysis::End(PHCompositeNode* /*topNode*/)
{
  std::cout << "\n========================================" << std::endl;
  std::cout << "ClusterPhaseAnalysis::End" << std::endl;
  std::cout << "========================================" << std::endl;
  std::cout << "Processed " << m_event << " events" << std::endl;
  std::cout << "Output file: " << m_outputFileName << std::endl;
  std::cout << "========================================\n" << std::endl;

  if (m_outputFile)
  {
    m_outputFile->cd();
    m_tree->Write();
    m_outputFile->Close();
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//_______________________________________________________________________________
int ClusterPhaseAnalysis::getNodes(PHCompositeNode* topNode)
{
    trackmap = findNode::getClass<SvtxTrackMap>(topNode, m_trackMapName);
    clustermap = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
    geometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
    
    return Fun4AllReturnCodes::EVENT_OK;
}
//___________________________________________________________________________________
void ClusterPhaseAnalysis::reset_tree_vars()
{


  m_truth_cluster_x.clear();
  m_truth_cluster_y.clear();
  m_truth_cluster_z.clear();
  m_truth_cluster_r.clear();
  m_truth_cluster_phi.clear();

  m_hit_x.clear();
  m_hit_y.clear();
  m_hit_z.clear();
  m_hit_adc.clear();
  m_hit_phase.clear();
  m_hit_layer.clear();

  m_phase_counts.clear();
}