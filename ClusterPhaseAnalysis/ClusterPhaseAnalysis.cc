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
#include <trackbase/TrkrClusterHitAssocv3.h>
#include <trackbase/ActsGeometry.h>  
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/TrackSeedContainer.h>



#include <g4detectors/PHG4CylinderGeomContainer.h>
//#include <g4detectors/PHG4TpcCylinderGeom.h>
//#include <g4detectors/PHG4TpcCylinderGeomContainer.h>

#include <g4detectors/PHG4TpcGeom.h>
#include <g4detectors/PHG4TpcGeomContainer.h>

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


  // ------ existing helpers (square, r, get_cluster_keys) stay here ------

  // Small epsilon for numeric safety
  inline double eps() { return 1e-12; }

  // Taubin circle fit in XY (y vs x)
  // Returns true on success and fills xc, yc, R
  bool fitCircleTaubinXY(const std::vector<std::pair<double,double> >& pts,
                         double& xc, double& yc, double& R)
  {
    const size_t n = pts.size();
    if (n < 3) return false;

    // 1) Compute centroid
    double meanx = 0.0, meany = 0.0;
    for (size_t i = 0; i < n; ++i) { meanx += pts[i].first; meany += pts[i].second; }
    meanx /= double(n);
    meany /= double(n);

    // 2) Centralized moments
    double Suu=0, Suv=0, Svv=0, Suuu=0, Svvv=0, Suvv=0, Svuu=0;
    for (size_t i = 0; i < n; ++i)
    {
      const double u = pts[i].first  - meanx;
      const double v = pts[i].second - meany;
      const double uu = u*u;
      const double vv = v*v;
      Suu  += uu;
      Svv  += vv;
      Suv  += u*v;
      Suuu += uu*u;
      Svvv += vv*v;
      Suvv += u*vv;
      Svuu += v*uu;
    }

    // 3) Solve linear system for (uc, vc) in centered coordinates
    const double A = Suu;
    const double B = Suv;
    const double C = Svv;
    const double D = 0.5*(Suuu + Suvv);
    const double E = 0.5*(Svvv + Svuu);

    const double det = A*C - B*B;
    if (std::fabs(det) < eps()) return false;

    const double uc = (D*C - B*E)/det;
    const double vc = (A*E - B*D)/det;

    xc = meanx + uc;
    yc = meany + vc;

    // 4) Radius
    double r2mean = 0.0;
    for (size_t i = 0; i < n; ++i)
    {
      const double dx = pts[i].first  - xc;
      const double dy = pts[i].second - yc;
      r2mean += dx*dx + dy*dy;
    }
    r2mean /= double(n);
    R = (r2mean > 0.0) ? std::sqrt(r2mean) : 0.0;

    return (R > 0.0);
  }

  // Intersect fitted circle (xc,yc,R) with origin-centered circle of radius rLayer.
  // Returns up to 2 points in outPts. True if intersection exists.
  bool intersectTwoCircles(double xc, double yc, double R, double rLayer,
                           std::pair<double,double>& p1,
                           std::pair<double,double>& p2,
                           bool& twoSolutions)
  {
    const double d = std::sqrt(xc*xc + yc*yc);
    // No intersection: too far or one circle contained in the other
    if (d < eps()) // concentric, treat specially
    {
      if (std::fabs(R - rLayer) < eps()) {
        // Infinite intersections; degenerate; return false to avoid ambiguity
        return false;
      } else {
        return false;
      }
    }
    if (d > R + rLayer + 1e-9) return false;
    if (d < std::fabs(R - rLayer) - 1e-9) return false;

    // a: distance from origin to the chord midpoint along c-hat
    const double a = (rLayer*rLayer + d*d - R*R) / (2.0*d);
    const double h2 = rLayer*rLayer - a*a;
    if (h2 < -1e-10) return false; // numerical
    const double h = (h2 <= 0.0) ? 0.0 : std::sqrt(std::max(0.0, h2));

    // Unit vector from origin to (xc,yc)
    const double ux = xc/d;
    const double uy = yc/d;

    // Base point on the line connecting origins
    const double x2 = a*ux;
    const double y2 = a*uy;

    // Perpendicular unit vector
    const double px = -uy;
    const double py =  ux;

    // Two intersection points
    p1 = std::make_pair(x2 + h*px, y2 + h*py);
    p2 = std::make_pair(x2 - h*px, y2 - h*py);
    twoSolutions = (h > eps());
    return true;
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
  if(clusterhitassocmap->size() == 0){
      std::cout << "WARNING:TRKR_CLUSTERHITASSOC is empty! Skip event" << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
  }else{
    std::cout << "TRKR_CLUSTERHITASSOC has size "<< clusterhitassocmap->size() << std::endl;
  }
 // clusterhitassocmap->identify();
  processTracks();

  return Fun4AllReturnCodes::EVENT_OK;
}




//___________________________________________________________________________________
void ClusterPhaseAnalysis::processTracks()
{
    std::cout<<"\n=========================================================="<<std::endl;
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
            if(!isSimulation) ComputeFitTruthAtLayer(track);
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
      std::cout<<"\n---------------------------------------------------------"<<std::endl;
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

      m_cluster_residual_rphi = residual_rphi_map.count(cluskey) ? residual_rphi_map[cluskey] : NAN;
      m_cluster_residual_time = residual_time_map.count(cluskey) ? residual_time_map[cluskey] : NAN;

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
  std::cout<<" hskey: "<<hskey<<std::endl;
  TrkrHitSet *hitset = hitmap->findHitSet(hskey);
  std::cout<<" hitset: "<<hitset<<std::endl;

  int hitlayer = TrkrDefs::getLayer(hskey);
  int hitside = TpcDefs::getSide(hskey);
  std::cout<<" hitlayer: "<<hitlayer<<" hitside: "<<hitside<<std::endl;




  
     std::cout <<"clusterhitassocmap size "<< clusterhitassocmap->size()<<std::endl;

  auto hitrange = clusterhitassocmap->getHits(ckey);


 /* if(!isSimulation){
    const auto idx   = TrkrDefs::getClusIndex(ckey);
    const TrkrDefs::hitsetkey hs_layer_only = static_cast<TrkrDefs::hitsetkey>(hitlayer);
    const TrkrDefs::cluskey cluskey_layer = TrkrDefs::genClusKey(hs_layer_only, idx);
    hitrange = clusterhitassocmap->getHits(cluskey_layer);

  }*/
  //if (Verbosity())
  //{
      std::cout << "ClusterPhaseAnalysis::FillHits -- Number of hits: "<<hitmap->size()<<" Number of associated hits: " << std::distance(hitrange.first, hitrange.second)<<" ( "<<clusterhitassocmap->size()<<" )" << std::endl;
  //}
  /*const auto det    = TrkrDefs::getTrkrId(hskey);
  if (det != TrkrDefs::tpcId) {
    std::cerr << "Non-TPC assoc encountered, skipping (det=" << int(det) << ")\n";
    return;
  }*/
    auto *geoLayer = tpcGeom->GetLayerCellGeom(hitlayer);
    geoLayer->identify();
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
        auto hitt =  geoLayer->get_zcenter(hittbin);
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
        m_hit_time.push_back(hitt);
       // m_hit_time.push_back(hittbin*AdcClockPeriod);
       
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
    std::cout<<"   hit phi "<<m_hit_phi[i]<<" time "<<m_hit_time[i]<<" adc "<<m_hit_adc[i]<<std::endl;
  }


  for (const auto& entry : sum_by_pad) 
  {
    std::cout<<"   phi "<<entry.first<<" adc sum "<<entry.second<<std::endl;
    if (entry.second > max_adc_sum_by_pad) 
    {
        max_adc_sum_by_pad = entry.second;
        max_adc_phi_bin = entry.first; 
      
        //max_adc_sum_by_pad_pad = entry.first;
    }
  }

  float dNdPhase_phi = (float)(m_cluster_phi - max_adc_phi_bin)/(phiwidth);


  for (const auto& entry : sum_by_tbin) 
  {
    std::cout<<"   time "<<entry.first<<" adc sum "<<entry.second<<std::endl;
    if (entry.second > max_adc_sum_by_tbin) 
    {
        max_adc_sum_by_tbin = entry.second;
        max_adc_time_bin = entry.first;
        //max_adc_sum_by_tbin_tbin = entry.first;
    }
  }

  float dNdPhase_time = (float)(m_cluster_time- max_adc_time_bin)/(twidth);

  m_cluster_size_phi = sum_by_pad.size();
  m_cluster_size_time = sum_by_tbin.size();
  m_cluster_phase_phi =dNdPhase_phi;
  m_cluster_phase_time = dNdPhase_time;
  if ((m_cluster_size_phi==2 && std::fabs(dNdPhase_phi)>0.5)||(m_cluster_size_phi==3 &&  std::fabs(dNdPhase_phi)>1)||(m_cluster_size_phi==4 &&  std::fabs(dNdPhase_phi)>1.5)||(m_cluster_size_phi==5 &&  std::fabs(dNdPhase_phi)>2))
  {
        m_cluster_phase_phi =0;
        m_cluster_phase_time =0;
        m_cluster_size_phi = 0;
        m_cluster_size_time = 0;
        std::cout<<"ClusterPhaseAnalysis::CalculatePhase -- ALARM! "<<std::endl;
        for(size_t i=0; i <  N ; ++i) 
        {
            std::cout<<"   hit pad "<<m_hit_pad[i]<<" time bin "<<m_hit_tbin[i]<<" adc "<<m_hit_adc[i]<<";   hit phi "<<m_hit_phi[i]<<" time "<<m_hit_time[i]<<std::endl;
        }
  }
 if (Verbosity())
  {
      std::cout << "ClusterPhaseAnalysis::CalculatePhase -- Cluster phi size = "<< m_cluster_size_phi <<"  Cluster t size = "<<m_cluster_size_time<< std::endl;

      std::cout<<"Cluster phi "<<m_cluster_phi<<" phase "<<dNdPhase_phi<<" max adc phi "<<max_adc_phi_bin<<" adc sum "<<max_adc_sum_by_pad<<" phiwidth "<<phiwidth<<std::endl;
      std::cout<<"Cluster time "<<m_cluster_time<<" phase time "<<dNdPhase_time<<" max adc time "<<max_adc_time_bin<<" adc sum "<<max_adc_sum_by_tbin<<" timebinwidth "<<twidth<<std::endl;
      std::cout<<"\n---------------------------------------------------------"<<std::endl;
  }
}
//_____________________________________________________________________________________________
void ClusterPhaseAnalysis::ComputeFitTruthAtLayer(SvtxTrack* track)
{
  const int minPointsForFit = 10;
  const int nLayerSkip = 1;

  // --- Step 1: Collect XY points per region (R1, R2, R3) ---
  std::array<std::vector<std::pair<double,double>>, 3> ptsXY;

  const auto clusterKeys = get_cluster_keys(track);
  ptsXY[0].reserve(64); ptsXY[1].reserve(64); ptsXY[2].reserve(64);

for (const auto& ckey : clusterKeys)
{
  const unsigned int layer = TrkrDefs::getLayer(ckey);
  const int region = (int)(layer - 7) / 16;
  if (region < 0 || region > 2) continue;

  const unsigned int regionFirst = 7 + static_cast<unsigned int>(region) * 16;
  const unsigned int regionLast  = regionFirst + 15;

  if (layer <= regionFirst + static_cast<unsigned int>(nLayerSkip) ||
      layer >= regionLast  - static_cast<unsigned int>(nLayerSkip))
    continue;

  TrkrCluster* c = clustermap->findCluster(ckey);
  if (!c) continue;

  const Acts::Vector3 g = geometry->getGlobalPosition(ckey, c);
  ptsXY[region].emplace_back(g.x(), g.y());
}


  // --- Step 2: Fit circles for all regions ---
  double Xc[3] = {0}, Yc[3] = {0}, R[3] = {0};
  for (int i = 0; i < 3; ++i)
  {
    if (ptsXY[i].size() < static_cast<size_t>(minPointsForFit))
      return; // abort if any region lacks enough clusters

    if (!fitCircleTaubinXY(ptsXY[i], Xc[i], Yc[i], R[i]))
      return;
  }

  // --- Step 3: Loop again to compute truth positions per cluster ---
  for (const auto& ckey : clusterKeys)
  {
    const unsigned int layer = TrkrDefs::getLayer(ckey);
    const int region = (layer - 7) / 16;
    if (region < 0 || region > 2) continue;

   auto* geoLayer = tpcGeom->GetLayerCellGeom(layer);
    if (!geoLayer) continue;
    const double rLayer = geoLayer->get_radius();

    TrkrCluster* c = clustermap->findCluster(ckey);
    if (!c) continue;

    const Acts::Vector3 g = geometry->getGlobalPosition(ckey, c);
    const double cx = g.x(), cy = g.y(), phiReco = std::atan2(cy, cx), rReco = std::hypot(cx, cy);

    std::pair<double,double> P1, P2;
    bool twoSol = false;
    if (!intersectTwoCircles(Xc[region], Yc[region], R[region], rLayer, P1, P2, twoSol))
      continue;

    // Pick intersection closest in azimuth
    const double phi1 = std::atan2(P1.second, P1.first);
    double pickx = P1.first, picky = P1.second;
    if (twoSol)
    {
      const double phi2 = std::atan2(P2.second, P2.first);
      if (std::fabs(std::remainder(phiReco - phi2, 2.0*M_PI)) <
          std::fabs(std::remainder(phiReco - phi1, 2.0*M_PI)))
      {
        pickx = P2.first; picky = P2.second;
      }
    }

    // --- Fill "fit-truth" info ---
    const double rFit = std::hypot(pickx, picky);
    const double phiFit = std::atan2(picky, pickx);

    m_truth_cluster_x    = float(pickx);
    m_truth_cluster_y    = float(picky);
    m_truth_cluster_z    = float(g.z());          // keep z same as cluster
    m_truth_cluster_r    = float(rFit);
    m_truth_cluster_phi  = float(phiFit);
    m_truth_cluster_time = m_cluster_time;

    // --- Residuals ---
  //  m_cluster_residual_rphi  = float(rReco * std::remainder(phiReco - phiFit, 2.0*M_PI));
    //m_cluster_residual_time  = std::numeric_limits<float>::quiet_NaN();
    residual_rphi_map[ckey] = float(rReco * std::remainder(phiReco - phiFit, 2.0*M_PI));
    residual_time_map[ckey] = std::numeric_limits<float>::quiet_NaN();
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
    if(isSimulation){
      truthClustersmap = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_TRUTHCLUSTERCONTAINER");
    }
    hitmap = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
   // tpcGeom = findNode::getClass<PHG4TpcCylinderGeomContainer>(topNode, "CYLINDERCELLGEOM_SVTX");
    tpcGeom = findNode::getClass<PHG4TpcGeomContainer>(topNode, "TPCGEOMCONTAINER");
    clusterhitassocmap = findNode::getClass<TrkrClusterHitAssoc>(topNode, "TRKR_CLUSTERHITASSOC");
    if (!clusterhitassocmap)
    {
      std::cerr << "ERROR: Can't find TRKR_CLUSTERHITASSOC node!" << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }
    std::cout<<" TRKR_CLUSTERHITASSOC size "<<clusterhitassocmap->size()<<std::endl;
    
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