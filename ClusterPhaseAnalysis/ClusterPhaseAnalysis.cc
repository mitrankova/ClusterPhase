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
#include <trackreco/ALICEKF.h>
#include <trackreco/GPUTPCTrackParam.h>






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
  struct SimpleKalmanState
{
  double x, y, z;     // position
  double tx, ty;      // slopes dx/dz, dy/dz
  double cov[5][5];   // covariance matrix
};

struct SimpleMeasurement
{
  double x, y, z;
  double sigma2; // measurement variance
};


  // Small epsilon for numeric safety
  inline double eps() { return 1e-12; }

  // Taubin circle fit in XY (y vs x)
  // Returns true on success and fills xc, yc, R

struct TrackCircleGuess {
  double xc = 0.;
  double yc = 0.;
  double R  = 0.;
};

TrackCircleGuess initialCircleFromTrack(const SvtxTrack* track, double Bz_T)
{
  TrackCircleGuess out;
  if (!track || std::fabs(Bz_T) < 1e-6) {
    out.R = 1e9;
    return out;
  }

  const double q = (track->get_charge() == 0 ? 1.0 : track->get_charge());
  const double pt = std::max(1e-6, double(track->get_pt())); // GeV/c
  const double phi0 = track->get_phi(); // rad
  const double x0 = track->get_x();     // cm
  const double y0 = track->get_y();     // cm

  // pT [GeV/c], B [T]  →  R [m];  multiply by 100 to get cm
  const double Rphys = 100.0 * pt / (0.3 * std::fabs(Bz_T)); // cm

  // charge–field orientation: curvature sign
  const double s = (q * Bz_T > 0) ? +1.0 : -1.0;

  // circle center is 90° rotated from momentum direction
  out.xc = x0 + s * Rphys * std::cos(phi0 + s * M_PI_2);
  out.yc = y0 + s * Rphys * std::sin(phi0 + s * M_PI_2);
  out.R  = Rphys;

  std::cout<<" initialCircleFromTrack: xc "<<out.xc<<" yc "<<out.yc<<" R "<<out.R<<std::endl;

  return out;
}


bool fitCircleTaubinXY(const std::vector<std::pair<double,double>>& pts,
                       double& xc, double& yc, double& R)
{
  const size_t n = pts.size();
  if (n < 3) return false;

  auto runWeighted = [](const std::vector<std::pair<double,double>>& p,
                        const std::vector<double>& w,
                        double& xc_out, double& yc_out, double& R_out) -> bool
  {
    double sumw = 0., meanx = 0., meany = 0.;
    for (size_t i=0;i<p.size();++i){
      sumw += w[i];
      meanx += p[i].first  * w[i];
      meany += p[i].second * w[i];
    }
    meanx /= sumw;  meany /= sumw;

    double Suu=0.,Svv=0.,Suv=0.,Suuu=0.,Svvv=0.,Suvv=0.,Svuu=0.;
    for (size_t i=0;i<p.size();++i){
      const double wi = w[i];
      const double u = p[i].first  - meanx;
      const double v = p[i].second - meany;
      Suu  += wi*u*u;
      Svv  += wi*v*v;
      Suv  += wi*u*v;
      Suuu += wi*u*u*u;
      Svvv += wi*v*v*v;
      Suvv += wi*u*v*v;
      Svuu += wi*v*u*u;
    }

    const double A = Suu;
    const double B = Suv;
    const double C = Svv;
    const double D = 0.5*(Suuu + Suvv);
    const double E = 0.5*(Svvv + Svuu);

    const double det = A*C - B*B;
    if (std::fabs(det) < 1e-14*(A*A + C*C)) return false;

    const double uc = (D*C - B*E) / det;
    const double vc = (A*E - B*D) / det;

    xc_out = meanx + uc;
    yc_out = meany + vc;

    double r2 = 0.;
    for (size_t i=0;i<p.size();++i){
      const double dx = p[i].first  - xc_out;
      const double dy = p[i].second - yc_out;
      r2 += w[i]*(dx*dx + dy*dy);
    }
    r2 /= sumw;
    R_out = std::sqrt(std::max(r2,0.0));
    return std::isfinite(R_out) && R_out>1e-6;
  };

  // --- initial unweighted ---
  std::vector<double> w(n,1.0);
  if(!runWeighted(pts,w,xc,yc,R)) return false;

  // --- iterative robust weighting ---
  const double huberK = 1.5; // cm cutoff
  for(int iter=0; iter<3; ++iter)
  {
    std::vector<double> res(n);
    for(size_t i=0;i<n;++i)
      res[i] = std::fabs(std::hypot(pts[i].first - xc, pts[i].second - yc) - R);

    for(size_t i=0;i<n;++i){
      double r = res[i];
      w[i] = (r<=huberK) ? 1.0 : huberK / r;
    }

    if(!runWeighted(pts,w,xc,yc,R)) break;
  }

  return true;
}


  // Intersect fitted circle (xc,yc,R) with origin-centered circle of radius rLayer.
  // Returns up to 2 points in outPts. True if intersection exists.
// ---------- helper: circle-circle intersection ----------
// Circle1: (xc,yc), R1 ; Circle2: (0,0), R2
bool intersectTwoCircles(
  double xc, double yc,
  double R1,          // track circle radius
  double R2,          // layer radius (around origin)
  std::pair<double,double>& P1,
  std::pair<double,double>& P2,
  bool& twoSol)
{
  const double dx = xc;
  const double dy = yc;
  const double d  = std::hypot(dx, dy);
  const double eps = 1e-9;

  // no intersection (too far apart, contained, or degenerate)
  if (d < eps || d > R1 + R2 + eps || d < std::fabs(R1 - R2) - eps)
  {
    twoSol = false;
    return false;
  }

  // *** THIS is the crucial line ***
  // distance from origin (layer center) to the "base point" along the line to (xc,yc)
  const double a = (R2*R2 - R1*R1 + d*d) / (2.0*d);
  // (your version used R1^2 - R2^2 + d^2, which is reversed)

  double h2 = R2*R2 - a*a;
  if (h2 < 0.0) h2 = 0.0;
  const double h = std::sqrt(h2);

  // base point on the line from origin to (xc,yc)
  const double px = (a/d) * dx;
  const double py = (a/d) * dy;

  // perpendicular offsets
  const double offx = (h/d) * (-dy);
  const double offy = (h/d) * ( dx);

  P1 = std::make_pair(px + offx, py + offy);
  P2 = std::make_pair(px - offx, py - offy);

  twoSol = (h > eps);
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
            if(!isSimulation) ComputeFitTruthAtLayer( track);
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
      m_truth_cluster_x = m_truth_cluster_x_map.count(cluskey) ? m_truth_cluster_x_map[cluskey] : NAN;
      m_truth_cluster_y = m_truth_cluster_y_map.count(cluskey) ? m_truth_cluster_y_map[cluskey] : NAN;
      m_truth_cluster_z = m_truth_cluster_z_map.count(cluskey) ? m_truth_cluster_z_map[cluskey] : NAN;

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
  if (!track || !geometry || !clustermap || !tpcGeom) return;

  // --- constants for region definition ---
  const unsigned int first_tpc_layer = 7;
  const unsigned int n_layers_region = 16;   // 3 regions of 16 layers each
  const unsigned int last_tpc_layer  = first_tpc_layer + 3*n_layers_region - 1; // 7+48-1 = 54

  // --- collect clusters for this track, grouped by region ---
  struct RegionData
  {
    std::vector<std::pair<double,double> > pts;      // (x,y)
    std::vector<TrkrDefs::cluskey>        keys;      // cluster keys in this region
  };

  RegionData regions[3];
  std::vector<std::pair<double,double> > allPts;

  const std::vector<TrkrDefs::cluskey> clusterKeys = get_cluster_keys(track);

  for (size_t i = 0; i < clusterKeys.size(); ++i)
  {
    const TrkrDefs::cluskey ckey = clusterKeys[i];
    const unsigned int layer = TrkrDefs::getLayer(ckey);

    // only TPC layers
    if (layer < first_tpc_layer || layer > last_tpc_layer) continue;

    const unsigned int idx = (layer - first_tpc_layer) / n_layers_region;
    if (idx > 2) continue;  // safety

    TrkrCluster* c = clustermap->findCluster(ckey);
    if (!c) continue;

    const Acts::Vector3 g = geometry->getGlobalPosition(ckey, c);
    const double x = g.x();
    const double y = g.y();

    regions[idx].pts.push_back(std::make_pair(x,y));
    regions[idx].keys.push_back(ckey);
    allPts.push_back(std::make_pair(x,y));
  }

  // nothing to do if we have no TPC clusters
  if (allPts.size() < 3) return;

  // --- Step 1: global circle seed from track parameters ---
  TrackCircleGuess init = initialCircleFromTrack(track, 1.4); // Bz ~ 1.4 T for sPHENIX

  double xc_global = init.xc;
  double yc_global = init.yc;
  double R_global  = init.R;

  // If track-based seed is bad, try a global Taubin fit on all points
  {
    double xc_tmp = xc_global;
    double yc_tmp = yc_global;
    double R_tmp  = R_global;

    if (fitCircleTaubinXY(allPts, xc_tmp, yc_tmp, R_tmp))
    {
      xc_global = xc_tmp;
      yc_global = yc_tmp;
      R_global  = R_tmp;
    }
    // else keep init from track
  }

  // --- Step 2: fit per region with robust Taubin, constrained by the global circle ---
  double xc_reg[3] = { xc_global, xc_global, xc_global };
  double yc_reg[3] = { yc_global, yc_global, yc_global };
  double R_reg[3]  = { R_global , R_global , R_global  };

  const double max_rel_R_deviation = 0.30; // 30% deviation from global R allowed

  for (int ireg = 0; ireg < 3; ++ireg)
  {
    std::vector<std::pair<double,double> >& pts = regions[ireg].pts;

    if (pts.size() < 3)
    {
      // not enough points in this region -> use global circle
      xc_reg[ireg] = xc_global;
      yc_reg[ireg] = yc_global;
      R_reg[ireg]  = R_global;
      continue;
    }

    // Optionally pre-filter using global circle to remove obvious outliers
    std::vector<std::pair<double,double> > use_pts;
    use_pts.reserve(pts.size());
    for (size_t ip = 0; ip < pts.size(); ++ip)
    {
      const double dx = pts[ip].first  - xc_global;
      const double dy = pts[ip].second - yc_global;
      const double r_to_center = std::sqrt(dx*dx + dy*dy);

      if (std::fabs(r_to_center - R_global) < 0.25 * R_global)
      {
        use_pts.push_back(pts[ip]);  // keep points reasonably close to global circle
      }
    }

    if (use_pts.size() < 3)
    {
      // filter was too strict → fall back to all points in region
      use_pts = pts;
    }

    double xc_fit = xc_global;
    double yc_fit = yc_global;
    double R_fit  = R_global;

    // robust Taubin fit on this region
    if (!fitCircleTaubinXY(use_pts, xc_fit, yc_fit, R_fit))
    {
      // if fit fails, use global
      xc_reg[ireg] = xc_global;
      yc_reg[ireg] = yc_global;
      R_reg[ireg]  = R_global;
      continue;
    }

    // if radius is wildly inconsistent with global, clamp back to global
    if (std::fabs(R_fit - R_global) / std::max(R_global, 1e-6) > max_rel_R_deviation)
    {
      xc_reg[ireg] = xc_global;
      yc_reg[ireg] = yc_global;
      R_reg[ireg]  = R_global;
    }
    else
    {
      xc_reg[ireg] = xc_fit;
      yc_reg[ireg] = yc_fit;
      R_reg[ireg]  = R_fit;
    }
  }

  // --- Step 3: for each cluster, compute intersection of its region circle with layer radius ---
  for (int ireg = 0; ireg < 3; ++ireg)
  {
    const std::vector<TrkrDefs::cluskey>& keys = regions[ireg].keys;
    if (keys.empty()) continue;

    for (size_t ik = 0; ik < keys.size(); ++ik)
    {
      const TrkrDefs::cluskey ckey = keys[ik];
      const unsigned int layer = TrkrDefs::getLayer(ckey);

      PHG4TpcGeom* geoLayer = tpcGeom->GetLayerCellGeom(layer);
      if (!geoLayer) continue;

      const double rLayer = geoLayer->get_radius();  // cylinder radius for this layer

      TrkrCluster* c = clustermap->findCluster(ckey);
      if (!c) continue;

      const Acts::Vector3 g = geometry->getGlobalPosition(ckey, c);
      const double xm = g.x();
      const double ym = g.y();
      const double zm = g.z();

      const double rReco = std::sqrt(xm*xm + ym*ym);

      // Intersection between fitted circle in this region and circle around beamline with radius rLayer
      std::pair<double,double> P1, P2;
      bool twoSol = false;
      if (!intersectTwoCircles(
             xc_reg[ireg], yc_reg[ireg], R_reg[ireg],
             rLayer,
             P1, P2, twoSol))
      {
        // no intersection → skip this cluster
        continue;
      }

      // choose solution closer in XY to measured cluster position
      const double d1 = std::sqrt( (P1.first - xm)*(P1.first - xm)
                                 + (P1.second - ym)*(P1.second - ym) );
      double d2 = std::numeric_limits<double>::infinity();
      if (twoSol)
      {
        d2 = std::sqrt( (P2.first - xm)*(P2.first - xm)
                      + (P2.second - ym)*(P2.second - ym) );
      }

      const double xfit = (d1 <= d2) ? P1.first  : P2.first;
      const double yfit = (d1 <= d2) ? P1.second : P2.second;

      const double phiReco = std::atan2(ym, xm);
      const double phiFit  = std::atan2(yfit, xfit);
      const double dphi    = std::remainder(phiReco - phiFit, 2.0 * M_PI);

      // store "fit truth" position in maps (z: keep measured for now)
      m_truth_cluster_x_map[ckey] = static_cast<float>(xfit);
      m_truth_cluster_y_map[ckey] = static_cast<float>(yfit);
      m_truth_cluster_z_map[ckey] = static_cast<float>(zm);

      std::cout<<"ClusterPhaseAnalysis::ComputeFitTruthAtLayer -- ckey "<<ckey
               <<" layer "<<layer
               <<" rLayer "<<rLayer
               <<" reco (x,y)=("<<xm<<","<<ym<<") rReco "<<rReco
               <<" fit (x,y)=("<<xfit<<","<<yfit<<")"
               <<" dphi "<<dphi
               <<std::endl;
      // residuals
      residual_rphi_map[ckey] = static_cast<float>(rReco * dphi);
      // time residual: we don't touch here, leave for other calibration
      residual_time_map[ckey] = std::numeric_limits<float>::quiet_NaN();
    }
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