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

  // std::cout<<" initialCircleFromTrack: xc "<<out.xc<<" yc "<<out.yc<<" R "<<out.R<<std::endl;

    return out;
  }
  // -------------------------------------------------------------
  // Solve 3x3 system A x = b (simple Gaussian elimination)
  // -------------------------------------------------------------
  bool solve3x3(const double A[3][3], const double b[3], double x[3])
  {
    // Copy so we can modify
    double M[3][3];
    double rhs[3];
    for (int i = 0; i < 3; ++i)
    {
      rhs[i] = b[i];
      for (int j = 0; j < 3; ++j) M[i][j] = A[i][j];
    }

    // Forward elimination
    for (int k = 0; k < 3; ++k)
    {
      // Find pivot
      int piv = k;
      double maxAbs = std::fabs(M[k][k]);
      for (int i = k + 1; i < 3; ++i)
      {
        double val = std::fabs(M[i][k]);
        if (val > maxAbs)
        {
          maxAbs = val;
          piv = i;
        }
      }

      if (maxAbs < 1e-20) return false;

      // Swap rows if needed
      if (piv != k)
      {
        for (int j = 0; j < 3; ++j) std::swap(M[k][j], M[piv][j]);
        std::swap(rhs[k], rhs[piv]);
      }

      // Eliminate below
      const double diag = M[k][k];
      for (int i = k + 1; i < 3; ++i)
      {
        const double f = M[i][k] / diag;
        rhs[i] -= f * rhs[k];
        for (int j = k; j < 3; ++j) M[i][j] -= f * M[k][j];
      }
    }

    // Back substitution
    for (int i = 2; i >= 0; --i)
    {
      double sum = rhs[i];
      for (int j = i + 1; j < 3; ++j) sum -= M[i][j] * x[j];
      if (std::fabs(M[i][i]) < 1e-20) return false;
      x[i] = sum / M[i][i];
    }
    return true;
  }

  // Sagitta model in C++: same as Python sagitta_model
  static void sagittaModelCpp(
      double S, double x0, double invR,
      const std::vector<double>& up,
      std::vector<double>& f_out)
  {
    const size_t n = up.size();
    f_out.resize(n);

    for (size_t i = 0; i < n; ++i)
    {
      const double dx  = up[i] - x0;
      const double dx2 = dx*dx;
      const double dx4 = dx2*dx2;
      const double dx6 = dx4*dx2;

      const double invR3 = invR*invR*invR;
      const double invR5 = invR3*invR*invR;

      f_out[i] =
        S
        - 0.5   * invR  * dx2
        - 0.125 * invR3 * dx4
        - 0.0625* invR5 * dx6;
    }
  }

  bool fitCircleSagittaXY(
    const std::vector<std::pair<double,double>>& pts,
    double& xc_out, double& yc_out, double& R_out)
  {
    const size_t n = pts.size();
    if (n < 3) return false;

    // ------------------------------------------------------------
    // 1) Straight line fit y = m x + b  (same as rotate_points_to_x_axis)
    // ------------------------------------------------------------
    double Sx = 0., Sy = 0., Sxx = 0., Sxy = 0.;
    for (size_t i = 0; i < n; ++i)
    {
      const double x = pts[i].first;
      const double y = pts[i].second;
      Sx  += x;
      Sy  += y;
      Sxx += x*x;
      Sxy += x*y;
    }

    const double denom = n*Sxx - Sx*Sx;
    double m = 0., b = 0.;
    if (std::fabs(denom) > 1e-20)
    {
      m = (n*Sxy - Sx*Sy) / denom;
      b = (Sy - m*Sx) / n;
    }
    else
    {
      m = 0.;
      b = (n > 0) ? (Sy / n) : 0.;
    }

    const double theta = std::atan(m);
    const double cth   = std::cos(theta);
    const double sth   = std::sin(theta);
    //std::cout<<" m = "<<m<<" theta = "<<theta<<" b = "<<b<<std::endl;

    // ------------------------------------------------------------
    // 2) Rotate to (x', y') where chord is horizontal
    // ------------------------------------------------------------
    std::vector<double> up(n), vp(n);
    for (size_t i = 0; i < n; ++i)
    {
      const double x = pts[i].first;
      const double y = pts[i].second;

      const double x_shift = x;
      const double y_shift = y - b; // line passes through origin

      // xr = cos θ * x + sin θ * (y-b)
      // yr = -sin θ * x + cos θ * (y-b)
      up[i] =  cth * x_shift + sth * y_shift;
      vp[i] = -sth * x_shift + cth * y_shift;
    // std::cout<<"  Horizonal line x = "<<up[i]<<" y = "<<vp[i]<<std::endl;
    }

    // ------------------------------------------------------------
    // 3) Initial guesses for S, x0, invR  (same as Python)
    // ------------------------------------------------------------
    double S0   = vp[0];
    double x0_0 = up[0];
    for (size_t i = 1; i < n; ++i)
    {
      if (vp[i] > S0)
      {
        S0   = vp[i];
        x0_0 = up[i];
      }
    }

    double S    = S0;          // sagitta height
    double x0   = x0_0;        // apex position in x'
    double invR = 1.0 / 100.0; // R ~ 100 cm
  // std::cout<<"Prelim Params 0 : S0  "<<S<<" x0 "<<x0<<" inv R "<<invR<<std::endl;
  // ------------------------------------------------------------
  // 4) Gauss–Newton with NUMERICAL Jacobian (like SciPy)
  // ------------------------------------------------------------
  const int    maxIter   = 50;
  const double tolDelta  = 1e-8;
  const double eps_param = 1e-6;

  std::vector<double> f(n), r(n);
  //std::cout << "Prelim Params 0 : S0 " << S << " x0 " << x0 << " invR " << invR << std::endl;

  //int it_used = 0;
  int nfev = 0;

  for (int iter = 0; iter < maxIter; ++iter)
  {
  // it_used = iter + 1;

    // baseline model and residuals
    sagittaModelCpp(S, x0, invR, up, f);
    nfev += 1;

    double cost = 0.0;
    for (size_t i = 0; i < n; ++i)
    {
      r[i] = f[i] - vp[i];
      cost += 0.5 * r[i] * r[i];
    }

    // optional: print per-iteration status
    // std::cout << "Iter " << iter << " cost = " << cost
    //           << "  S=" << S << " x0=" << x0 << " invR=" << invR << std::endl;

    double JTJ[3][3] = { {0.,0.,0.}, {0.,0.,0.}, {0.,0.,0.} };
    double JTr[3]    = { 0., 0., 0. };

    // --- compute perturbed models once per parameter ---
    std::vector<double> fS_plus(n), fS_minus(n);
    std::vector<double> fx0_plus(n), fx0_minus(n);
    std::vector<double> fR_plus(n), fR_minus(n);

    sagittaModelCpp(S + eps_param, x0, invR, up, fS_plus);
    sagittaModelCpp(S - eps_param, x0, invR, up, fS_minus);
    sagittaModelCpp(S, x0 + eps_param, invR, up, fx0_plus);
    sagittaModelCpp(S, x0 - eps_param, invR, up, fx0_minus);
    sagittaModelCpp(S, x0, invR + eps_param, up, fR_plus);
    sagittaModelCpp(S, x0, invR - eps_param, up, fR_minus);
    nfev += 6; // we did 6 extra evaluations

    for (size_t i = 0; i < n; ++i)
    {
      const double dfdS    = (fS_plus[i]  - fS_minus[i])  / (2.0 * eps_param);
      const double dfdx0   = (fx0_plus[i] - fx0_minus[i]) / (2.0 * eps_param);
      const double dfdinvR = (fR_plus[i]  - fR_minus[i])  / (2.0 * eps_param);

      const double J[3] = { dfdS, dfdx0, dfdinvR };
      const double ri   = r[i];

      for (int a = 0; a < 3; ++a)
      {
        JTr[a] += J[a] * ri;
        for (int j = 0; j < 3; ++j)
          JTJ[a][j] += J[a] * J[j];
      }
    }

    // Solve (J^T J) Δ = - J^T r
    for (int a = 0; a < 3; ++a) JTr[a] = -JTr[a];

    double delta[3] = {0.,0.,0.};
    if (!solve3x3(JTJ, JTr, delta)) break;

    const double alpha = 1.0; // SciPy uses line search, we keep full step
    S    += alpha * delta[0];
    x0   += alpha * delta[1];
    invR += alpha * delta[2];

    if (!std::isfinite(S) || !std::isfinite(x0) || !std::isfinite(invR))
      break;

    if (invR == 0.) invR = 1.0 / 100.0;

    const double maxDelta =
      std::max(std::fabs(delta[0]),
      std::max(std::fabs(delta[1]), std::fabs(delta[2])));

    if (maxDelta < tolDelta) break;
  }

  // after the loop:
  const double R = (invR != 0.0) ? (1.0 / invR) : 1e9;
  if (!std::isfinite(R) || R == 0.) return false;

  // final debug print (Python-style)
  double final_cost = 0.0;
  sagittaModelCpp(S, x0, invR, up, f);
  for (size_t i = 0; i < n; ++i)
  {
    const double ri = f[i] - vp[i];
    final_cost += 0.5 * ri * ri;
  }
  /*
  std::cout << "---------------- SAGITTA FIT RESULT ----------------\n"
            << "S       = " << S    << "\n"
            << "x0      = " << x0   << "\n"
            << "invR    = " << invR << "\n"
            << "R       = " << R    << "\n"
            << "Final cost (0.5 * sum r^2) = " << final_cost << "\n"
            << "Iterations = " << it_used << "\n"
            << "Function evaluations (approx) = " << nfev << "\n"
            << "----------------------------------------------------\n";

  */


        // signed radius
      // R = (invR != 0.0) ? (1.0 / invR) : 1e9;
        if (!std::isfinite(R)) return false;

        // geometry uses |R|
        const double Rmag = std::fabs(R);

        // 5) Back to (x,y): center in rotated frame, allow signed R
        const double Xc_prime = x0;
        const double Yc_prime = S - R;     // signed curvature

        const double xc_shift =  cth * Xc_prime - sth * Yc_prime;
        const double yc_shift =  sth * Xc_prime + cth * Yc_prime;

        const double xc = xc_shift;
        const double yc = yc_shift + b;

        if (!std::isfinite(xc) || !std::isfinite(yc))
            return false;

        xc_out = xc;
        yc_out = yc;
        R_out  = Rmag;     // return positive geometric radius
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


// wrap to (-pi, pi]
static inline float wrap_pi(float a)
{
  while (a >  M_PI) a -= 2.0f*M_PI;
  while (a <= -M_PI) a += 2.0f*M_PI;
  return a;
}

// check if phi is in (down, up] along the wrapped circle
static inline bool in_sector_interval(float phi, float down, float up)
{
  // normal case: interval does NOT cross the -pi/pi seam
  if (down <= up) return (phi > down && phi <= up);

  // seam-crossing case, e.g. down=+170deg, up=-170deg
  return (phi > down || phi <= up);
}

// returns 0..11, or -1 if something goes wrong
int sector_from_phi(float phi)
{
  const float SectorAngle = (float)(M_PI/6.0);     // 30 deg
  const float AngleStart  = (float)(M_PI - SectorAngle/2.0); // 165 deg

  const float p = wrap_pi(phi);

  for (int s = 0; s < 12; ++s)
  {
    float up   = wrap_pi(AngleStart - s*SectorAngle);
    float down = wrap_pi(up - SectorAngle);

    if (in_sector_interval(p, down, up)) return s;
  }

  return -1;
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
  m_tree->Branch("ntpc", &m_ntpc, "n_ntpc/I");
  m_tree->Branch("quality", &m_quality, "n_quality/I");
  m_tree->Branch("pt", &m_pt, "pt/F");
  m_tree->Branch("phi0", &m_phi0, "phi0/F");
  m_tree->Branch("eta0", &m_eta0, "eta0/F");



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
  m_tree->Branch("cluster_side", &m_cluster_side, "m_cluster_side/I");
  m_tree->Branch("cluster_sector", &m_cluster_sector, "m_cluster_sector/I");
  m_tree->Branch("cluster_e", &m_cluster_e, "m_cluster_e/F");
  m_tree->Branch("cluster_e_simple", &m_cluster_e_simple, "m_cluster_e_simple/F");
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

  // Simulation Truth cluster branches
  m_tree->Branch("sim_truth_cluster_x", &m_sim_truth_cluster_x,"m_sim_truth_cluster_x/F");
  m_tree->Branch("sim_truth_cluster_y", &m_sim_truth_cluster_y, "m_sim_truth_cluster_y/F");
  m_tree->Branch("sim_truth_cluster_z", &m_sim_truth_cluster_z, "m_sim_truth_cluster_z/F");
  m_tree->Branch("sim_truth_cluster_r", &m_sim_truth_cluster_r, "m_sim_truth_cluster_r/F");
  m_tree->Branch("sim_truth_cluster_phi", &m_sim_truth_cluster_phi, "m_sim_truth_cluster_phi/F");
  m_tree->Branch("sim_truth_cluster_time", &m_sim_truth_cluster_time, "m_sim_truth_cluster_time/F"); 
  m_tree->Branch("sim_cluster_residual_rphi", &m_sim_cluster_residual_rphi, "m_sim_cluster_residual_rphi/F");
  m_tree->Branch("sim_cluster_residual_time", &m_sim_cluster_residual_time, "m_sim_cluster_residual_time/F");
  m_tree->Branch("cluster_resolution_rphi", &m_cluster_resolution_rphi, "m_cluster_resolution_rphi/F");
  m_tree->Branch("cluster_resolution_phi", &m_cluster_resolution_phi, "m_cluster_resolution_phi/F");



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
    if((int)trackmap->size()>ntracksmax){
        std::cout << "WARNING: SvtxTrackMap size "<<trackmap->size()<<" exceeds maximum "<<ntracksmax<<"! Skip event"<<std::endl;
        return;
    } else if ((int)trackmap->size()<ntracksmin){
        std::cout << "WARNING: SvtxTrackMap is empty! Skip event"<<std::endl;
        return;
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
           // if(!isSimulation) 
            ComputeFitTruthAtLayer( track);
            ComputeLinearFitFromNeighbors(track);  
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
    if (m_ifzerofield == false)
    {
        m_pt = track->get_pt();
        m_phi0 = track->get_phi();
        m_eta0 = track->get_eta();
    }

    if (m_ifzerofield == false && (track->get_quality() > m_maxQuality))
    {
        m_quality = track->get_quality();
        std::cout<<" Track quality is higher then "<<m_maxQuality<<std::endl;
        ifGoodTrack = false;
        return false;
    }

     m_ntpc=0;
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
        std::cout << "ClusterPhaseAnalysis::checkTrack - pt: " << track->get_pt() <<" ntpc: "<<m_ntpc<<" qualiuty: "<<m_quality<<" If zero field: "<<m_ifzerofield<< std::endl;
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
      float phi_err = cluster->getPhiError();
     
    
      unsigned int layer = TrkrDefs::getLayer(cluskey);
      const int side   = (int)TpcDefs::getSide(cluskey);     // 0 or 1
      const int sector = sector_from_phi(phi);

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
      m_cluster_side = side;
      m_cluster_sector = sector;
      m_cluster_e = phi_err;


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
      std::cout<<"  !!!!!! cluster "<<cluskey<<" m_cluster_x "<<m_cluster_x<<" m_cluster_y "<<m_cluster_y<<" m_cluster_z "<<m_cluster_z<<std::endl;
      std::cout<<"  !!!!!! cluster "<<cluskey<<" m_truth_cluster_x "<<m_truth_cluster_x<<" m_truth_cluster_y "<<m_truth_cluster_y<<" m_truth_cluster_z "<<m_truth_cluster_z<<"   -- m_cluster_residual_rphi "<<m_cluster_residual_rphi<<std::endl;
      std::cout<<"FILL TREE"<<std::endl;


      m_cluster_resolution_rphi = resolution_rphi_map.count(cluskey) ? resolution_rphi_map[cluskey] : NAN;
      m_cluster_resolution_phi = resolution_phi_map.count(cluskey) ? resolution_phi_map[cluskey] : NAN;

      std::cout << "  !!!!!! cluster " << cluskey 
                << " resolution_rphi " << m_cluster_resolution_rphi 
                << " resolution_phi " << m_cluster_resolution_phi << std::endl;

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
      m_sim_truth_cluster_x = static_cast<float>(best_glob.x());
      m_sim_truth_cluster_y = static_cast<float>(best_glob.y());
      m_sim_truth_cluster_z = static_cast<float>(best_glob.z());
      m_sim_truth_cluster_time = best->getLocalY();
      m_sim_truth_cluster_r = std::sqrt(square(m_sim_truth_cluster_x) + square(m_sim_truth_cluster_y));
      m_sim_truth_cluster_phi = std::atan2(m_sim_truth_cluster_y, m_sim_truth_cluster_x);
      m_sim_cluster_residual_rphi = r_reco * std::remainder(phi_reco - m_sim_truth_cluster_phi, 2.f * static_cast<float>(M_PI));
      m_sim_cluster_residual_time = m_cluster_time - m_sim_truth_cluster_time;

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
  float radius = geoLayer->get_radius();
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
  int iphi2_sum = 0, iphi_sum=0;
    int adc_sum = 0;
  for (size_t i = 0; i < N; ++i) 
  {
    iphi_sum += m_hit_pad[i] * m_hit_adc[i];
    iphi2_sum += square(m_hit_pad[i]) * m_hit_adc[i];
    adc_sum += m_hit_adc[i];
    sum_by_tbin[m_hit_time[i]] += static_cast<double>(m_hit_adc[i]);
    sum_by_pad[m_hit_phi[i]] += static_cast<double>(m_hit_adc[i]);
  // std::cout<<"   hit phi "<<m_hit_phi[i]<<" time "<<m_hit_time[i]<<" adc "<<m_hit_adc[i]<<std::endl;
  }

  float  clusiphi = iphi_sum / adc_sum;
  double phi_cov = (iphi2_sum / adc_sum - square(clusiphi)) * pow(phiwidth, 2);
  double phi_err = radius*std::sqrt(phi_cov / (adc_sum * 0.14));
  m_cluster_e_simple = phi_err;

  for (const auto& entry : sum_by_pad) 
  {
   // std::cout<<"   phi "<<entry.first<<" adc sum "<<entry.second<<std::endl;
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
   // std::cout<<"   time "<<entry.first<<" adc sum "<<entry.second<<std::endl;
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
    if (layer == 7 || layer == 22 || layer == 23 || layer == 38 || layer == 39 || layer == 54) continue;
    if (layer < first_tpc_layer || layer > last_tpc_layer) continue;

    const unsigned int idx = (layer - first_tpc_layer) / n_layers_region;
    if (idx > 2) continue;  // safety

    TrkrCluster* c = clustermap->findCluster(ckey);
    if (!c) continue;

    const Acts::Vector3 g = geometry->getGlobalPosition(ckey, c);
    const double x = g.x();
    const double y = g.y();

    if (c->getEdge() != 0) continue;
    
      /*const int side   = (int)TpcDefs::getSide(ckey);          // 0 or 1
      const int sector = sector_from_phi(std::atan2(g.y(), g.x()));

      const std::vector<int>& allowed = good_sectors[side][(int)(layer - 7)/16];
      if (std::find(allowed.begin(), allowed.end(), sector) == allowed.end())
      {
        continue;  // not a “good” sector for this side/region
      }
*/
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

    if (fitCircleSagittaXY(allPts, xc_tmp, yc_tmp, R_tmp))
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

  //const double max_rel_R_deviation = 0.30; // 30% deviation from global R allowed

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
      /*
      const double dx = pts[ip].first  - xc_global;
      const double dy = pts[ip].second - yc_global;
      const double r_to_center = std::sqrt(dx*dx + dy*dy);

      if (std::fabs(r_to_center - R_global) < 0.25 * R_global)
      {
        use_pts.push_back(pts[ip]);  // keep points reasonably close to global circle
      }
      */
      use_pts.push_back(pts[ip]);
      std::cout<<"ClusterPhaseAnalysis::ComputeFitTruthAtLayer -- Region "<<ireg
               <<" using point x = "<<pts[ip].first<<" y = "<<pts[ip].second<<std::endl;
    }

    if (use_pts.size() < 3)
    {
      // filter was too strict → fall back to all points in region
      use_pts = pts;
    }

    double xc_fit = xc_global;
    double yc_fit = yc_global;
    double R_fit  = R_global;

    
    std::cout<<"ClusterPhaseAnalysis::ComputeFitTruthAtLayer -- Region "<<ireg
             <<" fitting "<<use_pts.size()<<" points. Global Fit xc = "<<xc_global<<" yc = "<<yc_global<<" R = "<<R_global<<std::endl;

    if (!fitCircleSagittaXY(use_pts, xc_fit, yc_fit, R_fit))
    {
      // if fit fails, use global
      xc_reg[ireg] = xc_global;
      yc_reg[ireg] = yc_global;
      R_reg[ireg]  = R_global;
      std::cout<<"ClusterPhaseAnalysis::ComputeFitTruthAtLayer -- Region "<<ireg
               <<" sagitta fit failed, using global circle for this region."<<std::endl;
      continue;
    }
      std::cout<<"ClusterPhaseAnalysis::ComputeFitTruthAtLayer -- Region "<<ireg
               <<" fitted weighted: xc = "<<xc_fit<<" yc = "<<yc_fit<<" R = "<<R_fit<<std::endl;


    // if radius is wildly inconsistent with global, clamp back to global
   /* if (std::fabs(R_fit - R_global) / std::max(R_global, 1e-6) > max_rel_R_deviation)
    {
      std::cout<<"ClusterPhaseAnalysis::ComputeFitTruthAtLayer -- Region "<<ireg
               <<" fitted R = "<<R_fit<<" deviates too much from global R = "<<R_global
               <<", using global circle for this region."<<std::endl;
      xc_reg[ireg] = xc_global;
      yc_reg[ireg] = yc_global;
      R_reg[ireg]  = R_global;
    }*/
   // else
    {

      xc_reg[ireg] = xc_fit;
      yc_reg[ireg] = yc_fit;
      R_reg[ireg]  = R_fit;
    }
          std::cout<<"ClusterPhaseAnalysis::ComputeFitTruthAtLayer -- Region "<<ireg
               <<" Final taken fit: xc = "<<xc_reg[ireg]<<" yc = "<<yc_reg[ireg]<<" R = "<<R_reg[ireg]<<std::endl;
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
             rReco,
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
void ClusterPhaseAnalysis::ComputeLinearFitFromNeighbors(SvtxTrack* track)
{
  if (!track || !geometry || !clustermap || !tpcGeom) return;

  // --- constants for region definition ---
  const unsigned int first_tpc_layer = 7;
  const unsigned int last_tpc_layer  = 54;

  // --- collect all clusters for this track, organized by layer ---
  std::map<unsigned int, TrkrDefs::cluskey> clustersByLayer;
  
  const std::vector<TrkrDefs::cluskey> clusterKeys = get_cluster_keys(track);

  for (size_t i = 0; i < clusterKeys.size(); ++i)
  {
    const TrkrDefs::cluskey ckey = clusterKeys[i];
    const unsigned int layer = TrkrDefs::getLayer(ckey);

    // only TPC layers
    if (layer == 7 || layer == 22 || layer == 23 || layer == 38 || layer == 39 || layer == 54) continue;
    if (layer < first_tpc_layer || layer > last_tpc_layer) continue;

    TrkrCluster* c = clustermap->findCluster(ckey);
    if (!c) continue;
    
    if (c->getEdge() != 0) continue;

    clustersByLayer[layer] = ckey;  // Only one cluster per layer
  }

  if (clustersByLayer.empty()) return;

  // --- For each cluster, fit line through neighbors in layer-1 and layer+1 ---
  for (auto& layerPair : clustersByLayer)
  {
    const unsigned int currentLayer = layerPair.first;
    const TrkrDefs::cluskey ckey = layerPair.second;

    // Get neighboring layers
    const unsigned int layerMinus1 = currentLayer - 1;
    const unsigned int layerPlus1  = currentLayer + 1;

    // Check if neighboring layers exist and have clusters
    auto itMinus1 = clustersByLayer.find(layerMinus1);
    auto itPlus1  = clustersByLayer.find(layerPlus1);

    if (itMinus1 == clustersByLayer.end() || itPlus1 == clustersByLayer.end())
    {
      std::cout << "ClusterPhaseAnalysis::ComputeLinearFitFromNeighbors -- Layer " << currentLayer
                << " missing neighbors, skipping" << std::endl;
      continue;
    }

    // Get the current cluster position
    TrkrCluster* c = clustermap->findCluster(ckey);
    if (!c) continue;

    const Acts::Vector3 gCurrent = geometry->getGlobalPosition(ckey, c);
    const double xCurrent = gCurrent.x();
    const double yCurrent = gCurrent.y();
   // const double zCurrent = gCurrent.z();
    const double rCurrent = std::sqrt(xCurrent*xCurrent + yCurrent*yCurrent);
    const double phiCurrent = std::atan2(yCurrent, xCurrent);

    // Get neighbor cluster from layer-1
    TrkrCluster* cMinus1 = clustermap->findCluster(itMinus1->second);
    if (!cMinus1) continue;
    
    const Acts::Vector3 gMinus1 = geometry->getGlobalPosition(itMinus1->second, cMinus1);
    const double xMinus1 = gMinus1.x();
    const double yMinus1 = gMinus1.y();
    const double rMinus1 = std::sqrt(xMinus1*xMinus1 + yMinus1*yMinus1);
    double phiMinus1 = std::atan2(yMinus1, xMinus1);

    // Get neighbor cluster from layer+1
    TrkrCluster* cPlus1 = clustermap->findCluster(itPlus1->second);
    if (!cPlus1) continue;
    
    const Acts::Vector3 gPlus1 = geometry->getGlobalPosition(itPlus1->second, cPlus1);
    const double xPlus1 = gPlus1.x();
    const double yPlus1 = gPlus1.y();
    const double rPlus1 = std::sqrt(xPlus1*xPlus1 + yPlus1*yPlus1);
    double phiPlus1 = std::atan2(yPlus1, xPlus1);

    // Handle phi wrapping around 2π for consistency
    phiMinus1 = phiCurrent + std::remainder(phiMinus1 - phiCurrent, 2.0 * M_PI);
    phiPlus1  = phiCurrent + std::remainder(phiPlus1 - phiCurrent, 2.0 * M_PI);

    // Linear fit with just 2 points: phi = a + b*r
    // slope b = (phi2 - phi1) / (r2 - r1)
    // intercept a = phi1 - b*r1
    double denom = rPlus1 - rMinus1;
    if (std::fabs(denom) < 1e-12)
    {
      std::cout << "ClusterPhaseAnalysis::ComputeLinearFitFromNeighbors -- ckey " << ckey
                << " layer " << currentLayer << " zero denominator in linear fit" << std::endl;
      continue;
    }

    double slope = (phiPlus1 - phiMinus1) / denom;
    double intercept = phiMinus1 - slope * rMinus1;

    // Calculate fitted phi at current r
    double phiFit = intercept + slope * rCurrent;

    // Calculate residual
    double dphi = std::remainder(phiCurrent - phiFit, 2.0 * M_PI);
    double residual_rphi = rCurrent * dphi;

    // Store results
    //m_truth_cluster_x_map[ckey] = static_cast<float>(rCurrent * std::cos(phiFit));
    //m_truth_cluster_y_map[ckey] = static_cast<float>(rCurrent * std::sin(phiFit));
    //m_truth_cluster_z_map[ckey] = static_cast<float>(zCurrent);
    
    resolution_rphi_map[ckey] = static_cast<float>(residual_rphi);
    resolution_phi_map[ckey] = static_cast<float>(dphi);
   // resolution_time_map[ckey] = std::numeric_limits<float>::quiet_NaN();

    std::cout << "ClusterPhaseAnalysis::ComputeLinearFitFromNeighbors -- ckey " << ckey
              << " layer " << currentLayer
              << " rCurrent " << rCurrent << " phiCurrent " << phiCurrent
              << " layer-1: r=" << rMinus1 << " phi=" << phiMinus1
              << " layer+1: r=" << rPlus1 << " phi=" << phiPlus1
              << " phiFit " << phiFit
              << " residual_rphi " << residual_rphi
              << std::endl;
  }
}
//_______________________________________________________________________________

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