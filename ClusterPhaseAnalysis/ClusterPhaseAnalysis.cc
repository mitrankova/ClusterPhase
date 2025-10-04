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



#include <g4detectors/PHG4CylinderGeomContainer.h>
#include <g4detectors/PHG4TpcCylinderGeom.h>
#include <g4detectors/PHG4TpcCylinderGeomContainer.h>

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

    if(!track)
      {
	return out;
      }
    //{track->get_silicon_seed(),}
    for (const auto& seed :  {track->get_tpc_seed()} )
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
  m_tree->Branch("vdrift", &m_vdrift, "vdrift/F");

  // Cluster-level branches
  m_tree->Branch("cluster_x", &m_cluster_x, "m_cluster_x/F");
  m_tree->Branch("cluster_y", &m_cluster_y, "m_cluster_y/F");
  m_tree->Branch("cluster_z", &m_cluster_z, "m_cluster_z/F");
  m_tree->Branch("cluster_time", &m_cluster_time, "m_cluster_time/F");
  m_tree->Branch("cluster_r", &m_cluster_r, "m_cluster_r/F");
  m_tree->Branch("cluster_phi", &m_cluster_phi, "m_cluster_phi/F");
  m_tree->Branch("cluster_eta", &m_cluster_eta, "m_cluster_eta/F");
  m_tree->Branch("cluster_adc", &m_cluster_adc, "m_cluster_adc/I");
  m_tree->Branch("cluster_max_adc", &m_cluster_max_adc, "m_cluster_max_adc/I");
  m_tree->Branch("cluster_layer", &m_cluster_layer, "m_cluster_layer/I");
  m_tree->Branch("cluster_size", &m_cluster_size, "m_cluster_size/I");
  m_tree->Branch("cluster_size_phi", &m_cluster_size_phi, "m_cluster_size_phi/I");
  m_tree->Branch("cluster_size_time", &m_cluster_size_time, "m_cluster_size_time/I");

  m_tree->Branch("cluster_phase_phi", &m_cluster_phase_phi, "m_cluster_phase_phi/F");
  m_tree->Branch("cluster_phase_time", &m_cluster_phase_time, "m_cluster_phase_time/F");

  // Truth cluster branches
  m_tree->Branch("truth_cluster_x", &m_truth_cluster_x,"m_truth_cluster_x/F");
  m_tree->Branch("truth_cluster_y", &m_truth_cluster_y, "m_truth_cluster_y/F");
  m_tree->Branch("truth_cluster_z", &m_truth_cluster_z, "m_truth_cluster_z/F");
  m_tree->Branch("truth_cluster_r", &m_truth_cluster_r, "m_truth_cluster_r/F");
  m_tree->Branch("truth_cluster_phi", &m_truth_cluster_phi, "m_truth_cluster_phi/F");
  m_tree->Branch("truth_cluster_time", &m_truth_cluster_time, "m_truth_cluster_time/F"); 
  m_tree->Branch("cluster_residual_rphi", &m_cluster_residual_rphi, "m_cluster_residual_rphi/F");
  m_tree->Branch("cluster_residual_time", &m_cluster_residual_time, "m_cluster_residual_time/F");


  // Hit-level branches (nested vectors)
  m_tree->Branch("hit_keys", &m_hitkeys);
  m_tree->Branch("hit_side", &m_hit_side);
  m_tree->Branch("hit_layer", &m_hit_layer);
  m_tree->Branch("hit_region", &m_hit_region);
  m_tree->Branch("hit_pad", &m_hit_pad);
  m_tree->Branch("hit_tbin", &m_hit_tbin);
  m_tree->Branch("hit_adc", &m_hit_adc);
  m_tree->Branch("hit_x", &m_hit_x);
  m_tree->Branch("hit_y", &m_hit_y);
  m_tree->Branch("hit_z", &m_hit_z);
  m_tree->Branch("hit_r", &m_hit_r);
  m_tree->Branch("hit_phi", &m_hit_phi);
  m_tree->Branch("hit_time", &m_hit_time);


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
    if (Verbosity())
    {
        std::cout << "ClusterPhaseAnalysis::FillClusters for track " << std::endl;
    }
     m_vdrift = geometry->get_drift_velocity();

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
      float r = std::sqrt(square(x) + square(y));
      float phi = atan2(y, x);
      float theta = atan2(r, z);
      float eta = -log(tan(theta / 2.0));
    
      int adc = cluster->getAdc();
      int max_adc = cluster->getMaxAdc();
    
      unsigned int layer = TrkrDefs::getLayer(cluskey);

      m_cluster_x = x;
      m_cluster_y = y;
      m_cluster_z = z;
      m_cluster_time = cluster->getLocalY();
      m_cluster_r = r;
      m_cluster_phi = phi;
      m_cluster_eta = eta;
      m_cluster_adc = adc;
      m_cluster_max_adc = max_adc;
      m_cluster_layer = layer;

      if(isSimulation){
        FindTruthClusters(cluskey, glob);
      }
      

      reset_tree_vars();



      FillHits(cluskey);  
      CalculatePhase(cluskey);


      std::cout<<"FILL TREE"<<std::endl;
      m_tree->Fill();
  }
}
//___________________________________________________________________________________
void ClusterPhaseAnalysis::FindTruthClusters(uint64_t key, Acts::Vector3 glob)
{
    if (Verbosity())
    {
        std::cout << "ClusterPhaseAnalysis::FindTruthClusters from simulation " << std::endl;
    }
    float rx = glob.x();
    float ry = glob.y();
    float rz = glob.z();
    float r_reco = std::sqrt(square(rx) + square(ry));
    float phi_reco = atan2(ry, rx);
    //float theta = atan2(r_reco, rz);
    uint32_t hskey = TrkrDefs::getHitSetKeyFromClusKey(key);
    auto tr_range = truthClustersmap->getClusters(hskey); 
    float best_metric = std::numeric_limits<float>::max();
    TrkrCluster* best = nullptr;
    Acts::Vector3 best_glob = Acts::Vector3::Zero();

    for (auto it = tr_range.first; it != tr_range.second; ++it)
    {
        uint64_t tckey = it->first;  
        TrkrCluster* tclus = it->second;
         if (!tclus) continue;


        Acts::Vector3 tglob = geometry->getGlobalPosition(tckey, tclus);
        float tx = tglob.x();
        float ty = tglob.y();
        float tz = tglob.z();


        float phi_truth = std::atan2(ty, tx);
        float dphi = std::remainder(phi_reco - phi_truth, 2.f * static_cast<float>(M_PI));
        float metric = std::hypot(rz - tz, r_reco * dphi);
      
        if (metric < best_metric) 
        { 
          best_metric = metric; 
          best = tclus; 
          best_glob = tglob; 
        }
    }

    if (best)
    {
      m_truth_cluster_x = static_cast<float>(best_glob.x());
      m_truth_cluster_y = static_cast<float>(best_glob.y());
      m_truth_cluster_z = static_cast<float>(best_glob.z());
      m_truth_cluster_time = best->getLocalY();
      m_truth_cluster_r = std::sqrt(square(m_truth_cluster_x) + square(m_truth_cluster_y));
      m_truth_cluster_phi = std::atan2(m_truth_cluster_y, m_truth_cluster_x);
      m_cluster_residual_rphi = r_reco * std::remainder(phi_reco - m_truth_cluster_phi, 2.f * static_cast<float>(M_PI));
      m_cluster_residual_time = m_cluster_time - m_truth_cluster_time;

      //std::cout<<"FILL TREE: truth cluster "<<std::endl;
     // m_tree->Fill();
    }
}
//_____________________________________________________________________________________________
void ClusterPhaseAnalysis::FillHits(uint64_t ckey)
{
  if (Verbosity())
  {
      std::cout << "ClusterPhaseAnalysis::FillHits for cluster " <<ckey<< std::endl;
  }
    if (!clusterhitassocmap)
  {
    std::cerr << "ERROR: Can't find TRKR_CLUSTERHITASSOC node!" << std::endl;
    return;
  }

  /*

  size_t total_assoc_tpc = 0;
  for (unsigned int layer = 7; layer <= 54; ++layer)
  {
 
  auto crange = clustermap->getClusters(layer); // returns pair<ConstIter, ConstIter>
  for (auto cit = crange.first; cit != crange.second; ++cit)
  {
   
    const TrkrDefs::cluskey k = cit->first;
    uint64_t  ku = cit->first;
    auto hrange = clusterhitassocmap->getHits(k);
    const size_t nh = std::distance(hrange.first, hrange.second);
    total_assoc_tpc += nh;

    if (Verbosity() && nh > 0)
    {
      const auto hskey = TrkrDefs::getHitSetKeyFromClusKey(k);
      std::cout << "layer=" << layer
                << " hskey= " << hskey 
                << " ckey= "   << k  << " ckey(uint64_t)= "   << ku 
                << " assoc_hits=" << nh << " -- "<<(hrange.first)->second<<" "<<(hrange.second)->second<<"\n";
    }
  }
}
std::cout << "Total TPC associated hits counted via per-cluster queries: "
          << total_assoc_tpc << std::endl;
*/

  uint32_t hskey = TrkrDefs::getHitSetKeyFromClusKey(ckey);
  TrkrHitSet *hitset = hitmap->findHitSet(hskey);

  int hitlayer = TrkrDefs::getLayer(hskey);
  int hitside = TpcDefs::getSide(hskey);
  auto *geoLayer = tpcGeom->GetLayerCellGeom(hitlayer);

  auto hitrange = clusterhitassocmap->getHits(ckey);
 /* if(!isSimulation){
    const auto idx   = TrkrDefs::getClusIndex(ckey);
    const TrkrDefs::hitsetkey hs_layer_only = static_cast<TrkrDefs::hitsetkey>(hitlayer);
    const TrkrDefs::cluskey cluskey_layer = TrkrDefs::genClusKey(hs_layer_only, idx);
    hitrange = clusterhitassocmap->getHits(cluskey_layer);

  }*/
  /*if (Verbosity())
  {
      std::cout << "ClusterPhaseAnalysis::FillHits -- Number of hits: "<<hitmap->size()<<" Number of associated hits: " << std::distance(hitrange.first, hitrange.second)<<" ( "<<clusterhitassocmap->size()<<" )" << std::endl;
  }*/
  /*const auto det    = TrkrDefs::getTrkrId(hskey);
  if (det != TrkrDefs::tpcId) {
    std::cerr << "Non-TPC assoc encountered, skipping (det=" << int(det) << ")\n";
    return;
  }*/

  for (auto hit_iter = hitrange.first; hit_iter != hitrange.second; ++hit_iter)
  {
        // std::cout<<"clus key: "<<hit_iter->first<<" hit key "<<hit_iter->second<<" pad = "<<static_cast<uint16_t>((hit_iter->second >> 16) & 0xFFFF)<<" tbin = "<<static_cast<uint16_t>((hit_iter->second ) & 0xFFFF)<<std::endl;
       //auto *hit = hitr->second;
        uint32_t hitkey = hit_iter->second;
      
        //std::cout<<"clus key: "<<hitkey<<std::endl;
        

        auto hitpad = TpcDefs::getPad(hitkey);
        auto hittbin = TpcDefs::getTBin(hitkey);
        TrkrHit *hit = hitset->getHit(hitkey);
       // std::cout<<"   TPC hit with key: "<<hitkey<<" pad: "<<hitpad<<" tbin: "<<hittbin <<" adc: "<<hit->getAdc()<<std::endl;
        auto hitphi = geoLayer->get_phicenter(hitpad, hitside);
        auto hitradius = geoLayer->get_radius();
        float AdcClockPeriod = geoLayer->get_zstep();
        auto glob = geometry->getGlobalPositionTpc(hskey, hitkey, hitphi, hitradius, AdcClockPeriod);
        float hitgx = glob.x();
        float hitgy = glob.y();
        float hitgz = glob.z();
        //float hit_phi = phi; 

        m_hitkeys.push_back(hitkey);
        m_hit_side.push_back(hitside);
        m_hit_layer.push_back(hitlayer);
        m_hit_region.push_back((int)(hitlayer-7)/16);
        m_hit_pad.push_back((int)hitpad);
        m_hit_tbin.push_back((int)hittbin);
        m_hit_adc.push_back(hit->getAdc());
        m_hit_x.push_back(hitgx);
        m_hit_y.push_back(hitgy);
        m_hit_z.push_back(hitgz);
        m_hit_phi.push_back(hitphi);
        m_hit_r.push_back(hitradius);
        m_hit_time.push_back(hittbin*AdcClockPeriod);
       
  }
}

//_____________________________________________________________________________________________
void ClusterPhaseAnalysis::CalculatePhase(uint64_t ckey)
{
  if (Verbosity())
  {
      std::cout << "ClusterPhaseAnalysis::CalculatePhase" << std::endl;
  }

  uint32_t hskey = TrkrDefs::getHitSetKeyFromClusKey(ckey);
 // TrkrHitSet *hitset = hitmap->findHitSet(hskey);
  int hitlayer = TrkrDefs::getLayer(hskey);
  auto *geoLayer = tpcGeom->GetLayerCellGeom(hitlayer);
  float phiwidth = geoLayer->get_phistep();
  float twidth = geoLayer->get_zstep();
  const size_t N = m_hit_adc.size();

  if (Verbosity())
  {
      std::cout << "ClusterPhaseAnalysis::CalculatePhase -- Cluster phi size = "<< m_cluster_size_phi <<"  Cluster t size = "<<m_cluster_size_time<< std::endl;
  }
  if (N == 0) 
  {
    m_cluster_phase_phi = std::numeric_limits<float>::quiet_NaN();
    m_cluster_phase_time = std::numeric_limits<float>::quiet_NaN();
     std::cout << "ClusterPhaseAnalysis::CalculatePhase -- No associated hits found for this cluster "<< std::endl;
    return;
  }
  std::map<float, int> sum_by_tbin;
  std::map<float, int> sum_by_pad;
  //int max_adc_sum_by_pad_idx =0,max_adc_sum_by_tbin_idx =0;
  int max_adc_sum_by_pad = 0, max_adc_sum_by_tbin =0;
  //int  max_adc_sum_by_pad_pad = -1, max_adc_sum_by_tbin_tbin=-1;
  float max_adc_phi_bin = 0, max_adc_time_bin = 0;

  for (size_t i = 0; i < N; ++i) 
  {
    sum_by_tbin[m_hit_time[i]] += static_cast<double>(m_hit_adc[i]);
    sum_by_pad[m_hit_phi[i]] += static_cast<double>(m_hit_adc[i]);
  }


  for (const auto& entry : sum_by_pad) 
  {
    if (entry.second > max_adc_sum_by_pad) 
    {
        max_adc_sum_by_pad = entry.second;
        max_adc_phi_bin = entry.first; 
        //max_adc_sum_by_pad_pad = entry.first;
    }
  }
  float max_adc_phi = max_adc_phi_bin;
  float dNdPhase_phi = (float)(m_cluster_phi - max_adc_phi)/(phiwidth);


  for (const auto& entry : sum_by_tbin) 
  {
    if (entry.second > max_adc_sum_by_tbin) 
    {
        max_adc_sum_by_tbin = entry.second;
        max_adc_time_bin = entry.first;
        //max_adc_sum_by_tbin_tbin = entry.first;
    }
  }
  
  float max_adc_time = max_adc_time_bin;
  float dNdPhase_time = (float)(m_cluster_time- max_adc_time)/(twidth);

  m_cluster_phase_phi =dNdPhase_phi;
  m_cluster_phase_time = dNdPhase_time;
  m_cluster_size_phi = sum_by_pad.size();
  m_cluster_size_time = sum_by_tbin.size();
               // dNdPhase = std::round(dNdPhase * 100) / 100;
                //float phi_truth = atan2(gy_clust_truth, gx_clust_truth);
  
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
    if(isSimulation){
      truthClustersmap = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_TRUTHCLUSTERCONTAINER");
    }
    hitmap = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
    tpcGeom = findNode::getClass<PHG4TpcCylinderGeomContainer>(topNode, "CYLINDERCELLGEOM_SVTX");
    clusterhitassocmap = findNode::getClass<TrkrClusterHitAssoc>(topNode, "TRKR_CLUSTERHITASSOC");
    if (!clusterhitassocmap)
    {
      std::cerr << "ERROR: Can't find TRKR_CLUSTERHITASSOC node!" << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }
    
    return Fun4AllReturnCodes::EVENT_OK;
}
//___________________________________________________________________________________
void ClusterPhaseAnalysis::reset_tree_vars()
{

        m_hitkeys.clear();
        m_hit_side.clear();
        m_hit_layer.clear();
        m_hit_region.clear();
        m_hit_pad.clear();
        m_hit_tbin.clear();
        m_hit_adc.clear();
        m_hit_x.clear();
        m_hit_y.clear();
        m_hit_z.clear();
        m_hit_r.clear();
        m_hit_phi.clear();
        m_hit_time.clear();


}