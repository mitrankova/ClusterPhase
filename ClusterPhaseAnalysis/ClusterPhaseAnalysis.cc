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

    //! get cluster keys from a given track
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
int ClusterPhaseAnalysis::InitRun(PHCompositeNode* /*topNode*/)
{
  std::cout << "ClusterPhaseAnalysis::InitRun - Starting run" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}
//___________________________________________________________________________________
int ClusterPhaseAnalysis::process_event(PHCompositeNode* /*topNode*/)
{
  if (Verbosity())
  {
    std::cout << "ClusterPhaseAnalysis::process_event - Processing event " << m_event << std::endl;
  }

  // Reset tree variables
 // reset_tree_vars();

  // Get node pointers
 // auto *geometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
//  auto *clustermap = findNode::getClass<TrkrClusterContainer>(topNode, m_clusterContainerName);
 // auto *hitmap = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
 // auto *clusterhitassocmap = findNode::getClass<TrkrClusterHitAssoc>(topNode, "TRKR_CLUSTERHITASSOC");
  //auto *truthClustersmap = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_TRUTHCLUSTERCONTAINER");
/*
  if (!clustermap)
  {
    std::cout << "ClusterPhaseAnalysis::process_event - Error: cluster container " 
              << m_clusterContainerName << " not found" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  if (!hitmap)
  {
    std::cout << "ClusterPhaseAnalysis::process_event - Error: TRKR_HITSET not found" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  if (!clusterhitassocmap)
  {
    std::cout << "ClusterPhaseAnalysis::process_event - Error: TRKR_CLUSTERHITASSOC not found" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
*/
  /*auto *gl1 = findNode::getClass<Gl1RawHit>(topNode, "GL1RAWHIT");
  if (gl1)
  {
    m_bco = gl1->get_bco();
    auto lbshift = m_bco << 24U;
    m_bcotr = lbshift >> 24U;
  }
  else
  {
    Gl1Packet* gl1PacketInfo = findNode::getClass<Gl1Packet>(topNode, "GL1RAWHIT");
    if (!gl1PacketInfo)
    {
      m_bco = std::numeric_limits<uint64_t>::quiet_NaN();
      m_bcotr = std::numeric_limits<uint64_t>::quiet_NaN();
    }
    m_firedTriggers.clear();

    if (gl1PacketInfo)
    {
      m_gl1BunchCrossing = gl1PacketInfo->getBunchNumber();
      uint64_t triggervec = gl1PacketInfo->getScaledVector();
      m_bco = gl1PacketInfo->getBCO();
      auto lbshift = m_bco << 24U;
      m_bcotr = lbshift >> 24U;
      for (int i = 0; i < 64; i++)
      {
        bool trig_decision = ((triggervec & 0x1U) == 0x1U);
        if (trig_decision)
        {
          m_firedTriggers.push_back(i);
        }
        triggervec = (triggervec >> 1U) & 0xffffffffU;
      }
    }
  }*/

  // Initialize phase histogram
   // m_phase_counts.resize(m_n_phase_bins, 0);
    processTracks();
  // Loop over all clusters
 // 
 /*

    // Second pass: calculate phases and store hit information
    int cluster_size = 0;
    int cluster_phase_value = 0;

    hit_range = clusterhitassocmap->getHits(cluskey);
    
    for (auto hit_iter = hit_range.first; hit_iter != hit_range.second; ++hit_iter)
    {
      TrkrDefs::hitkey hitkey = hit_iter->second;
      
      cluster_size++;

      // Get the hitset
      TrkrHitSetContainer::Iterator hitset_iter = hitmap->findOrAddHitSet(hitsetkey);
      if (hitset_iter == hitmap->getHitSets().second)
      {
        continue;
      }
      
      TrkrHitSet *hitset = hitset_iter->second;
      if (!hitset) continue;
      
      TrkrHit *hit = hitset->getHit(hitkey);
      if (!hit) continue;

      // Get hit ADC
      int hit_adc_val = hit->getAdc();
      
      // Calculate phase relative to maximum ADC time bin
      int current_tbin = TpcDefs::getTBin(hitkey);
      int phase = current_tbin - max_adc_tbin;

      // If this is the maximum ADC hit, set cluster phase
      if (hit_adc_val == max_adc_value && current_tbin == max_adc_tbin)
      {
        cluster_phase_value = 0;  // By definition, phase = 0 for max ADC
      }

      // Get hit layer
      unsigned int hit_layer = TrkrDefs::getLayer(hitsetkey);
      
      // Store hit information
      // Note: For hit position, we use cluster position as approximation
      // For more precise hit positions, you would need to decode the hitkey
      hit_x_vec.push_back(x);
      hit_y_vec.push_back(y);
      hit_z_vec.push_back(z);
      hit_adc_vec.push_back(hit_adc_val);
      hit_phase_vec.push_back(phase);
      hit_layer_vec.push_back(hit_layer);
    }

    m_cluster_size.push_back(cluster_size);
    m_cluster_phase.push_back(cluster_phase_value);

    // Update phase histogram
    int bin_index = cluster_phase_value - m_phase_min;
    if (bin_index >= 0 && bin_index < m_n_phase_bins)
    {
      m_phase_counts[bin_index]++;
    }

    // Store hit information for this cluster
    m_hit_x.push_back(hit_x_vec);
    m_hit_y.push_back(hit_y_vec);
    m_hit_z.push_back(hit_z_vec);
    m_hit_adc.push_back(hit_adc_vec);
    m_hit_phase.push_back(hit_phase_vec);
    m_hit_layer.push_back(hit_layer_vec);
  }
}

  // Process truth clusters if available
  if (truthClustersmap)
  {
    auto truth_range = truthClustersmap->getClusters();
    for (auto iter = truth_range.first; iter != truth_range.second; ++iter)
    {
      TrkrCluster *truth_cluster = iter->second;
      if (!truth_cluster) continue;

      float x = truth_cluster->getX();
      float y = truth_cluster->getY();
      float z = truth_cluster->getZ();
      float r = sqrt(x * x + y * y);
      float phi = atan2(y, x);

      m_truth_cluster_x.push_back(x);
      m_truth_cluster_y.push_back(y);
      m_truth_cluster_z.push_back(z);
      m_truth_cluster_r.push_back(r);
      m_truth_cluster_phi.push_back(phi);
    }
  }

  // Fill tree
  m_tree->Fill();
  m_event++;
*/
//return returnVal;
  return Fun4AllReturnCodes::EVENT_OK;
}

//___________________________________________________________________________________
/*
int ClusterPhaseAnalysis::get_phase_from_hit(TrkrDefs::hitkey hitkey, 
                                              TrkrDefs::hitsetkey hitsetkey,
                                              TrkrHitSetContainer *hitmap)
{
  // This function is kept for potential future use
  // The phase calculation is now done directly in process_event
  
  // Get the hitset
  TrkrHitSetContainer::Iterator hitset_iter = hitmap->findOrAddHitSet(hitsetkey);
  if (hitset_iter == hitmap->getHitSets().second)
  {
    return 0;
  }

  TrkrHitSet *hitset = hitset_iter->second;
  if (!hitset)
  {
    return 0;
  }
  
  // Find the time bin with maximum ADC
  int max_adc = -1;
  int max_tbin = -1;
  
  // Get the current hit's time bin
  int current_tbin = TpcDefs::getTBin(hitkey);
  
  // Scan through all hits in this hitset to find max ADC tbin
  TrkrHitSet::ConstRange hit_range = hitset->getHits();
  for (TrkrHitSet::ConstIterator hit_iter = hit_range.first; 
       hit_iter != hit_range.second; ++hit_iter)
  {
    TrkrHit *hit = hit_iter->second;
    if (!hit) continue;
    
    int adc = hit->getAdc();
    
    if (adc > max_adc)
    {
      max_adc = adc;
      max_tbin = TpcDefs::getTBin(hit_iter->first);
    }
  }
  
  // Calculate phase relative to max ADC time bin
  // Phase = 0: tbin with largest ADC
  // Phase = -1: one tbin before max
  // Phase = +1: one tbin after max
  int phase = current_tbin - max_tbin;
  
  return phase;
}
*/
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

        
        if (checkTrack(track)){

            FillClusters(track);

        }
    }

 //   return Fun4AllReturnCodes::EVENT_OK;
}

//___________________________________________________________________________________
bool ClusterPhaseAnalysis::checkTrack(SvtxTrack* track)
{


    if (track->get_pt() < m_minPt)
    {
      //  m_rejected_tracks_pt++;
        std::cout<<" Track pt is lower then "<<m_minPt<<std::endl;
        return false;
    }

    if (m_ifzerofield == false && (track->get_quality() > m_maxQuality))
    {
       // m_rejected_tracks_quality++;
        std::cout<<" Track quality is higher then "<<m_maxQuality<<std::endl;
        return false;
    }

    // make sure cluster is from TPC
    int m_ntpc=0;
    for (const auto& ckey : get_cluster_keys(track))
    {
        const auto detId = TrkrDefs::getTrkrId(ckey);
        if (detId == TrkrDefs::tpcId)
        {
            m_ntpc++;
        }else{
            return false;
        }
    }

    if (m_ntpc<m_minTpcClusters) 
    {
      //  m_rejected_tracks_nhits++;
        return false;
    }

    if (Verbosity() > 2)
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
      // Get cluster properties
      float x = glob.x();
      float y = glob.y();
      float z = glob.z();
      float r = sqrt(x * x + y * y);
      float phi = atan2(y, x);
      float theta = atan2(r, z);
      float eta = -log(tan(theta / 2.0));
      
      // Get cluster error/energy
      //float e = 0;//cluster->getError(0, 0);
      
      // Get cluster ADC
      int adc = cluster->getAdc();
      int max_adc = cluster->getMaxAdc();
      
      // Get layer
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

      // Get hitset key for this cluster
      //TrkrDefs::hitsetkey hitsetkey = TrkrDefs::getHitSetKeyFromClusKey(cluskey);
      /*
      // Get associated hits for this cluster
      std::vector<float> hit_x_vec;
      std::vector<float> hit_y_vec;
      std::vector<float> hit_z_vec;
      std::vector<int> hit_adc_vec;
      std::vector<int> hit_phase_vec;
      std::vector<unsigned int> hit_layer_vec;
     
      // First pass: find the hit with maximum ADC in this cluster
      int max_adc_value = -1;
      int max_adc_tbin = -1;

      auto hit_range = clusterhitassocmap->getHits(cluskey);
      
      // Find hit with maximum ADC
      for (auto hit_iter = hit_range.first; hit_iter != hit_range.second; ++hit_iter)
      {
        TrkrDefs::hitkey hitkey = hit_iter->second;
        
        // Get the hitset
        TrkrHitSetContainer::Iterator hitset_iter = hitmap->findOrAddHitSet(hitsetkey);
        if (hitset_iter == hitmap->getHitSets().second)
        {
          continue;
        }
        
        TrkrHitSet *hitset = hitset_iter->second;
        if (!hitset) continue;
        
        TrkrHit *hit = hitset->getHit(hitkey);
        if (!hit) continue;

        int hit_adc_val = hit->getAdc();
        
        if (hit_adc_val > max_adc_value)
        {
          max_adc_value = hit_adc_val;
          max_adc_tbin = TpcDefs::getTBin(hitkey);
        }
      }
      */

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