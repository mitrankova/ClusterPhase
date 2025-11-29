#ifndef MACRO_FUN4ALLG4SPHENIX_C
#define MACRO_FUN4ALLG4SPHENIX_C

#include <GlobalVariables.C>

#include <DisplayOn.C>
#include <G4Setup_sPHENIX.C>
#include <G4_Mbd.C>
#include <G4_CaloTrigger.C>
#include <G4_Centrality.C>
#include <G4_DSTReader.C>
#include <G4_Global.C>
#include <G4_HIJetReco.C>
#include <G4_Input.C>
#include <G4_Jets.C>
#include <G4_KFParticle.C>
#include <G4_ParticleFlow.C>
#include <G4_Production.C>
#include <G4_TopoClusterReco.C>

#include <Trkr_RecoInit.C>
#include <Trkr_Clustering.C>
#include <Trkr_LaserClustering.C>
#include <Trkr_Reco.C>
#include <Trkr_Eval.C>
#include <Trkr_QA.C>

#include <Trkr_Diagnostics.C>
#include <G4_User.C>
#include <QA.C>

#include <ffamodules/FlagHandler.h>
#include <ffamodules/HeadReco.h>
#include <ffamodules/SyncReco.h>
#include <ffamodules/CDBInterface.h>

#include <fun4all/Fun4AllDstOutputManager.h>
#include <fun4all/Fun4AllOutputManager.h>
#include <fun4all/Fun4AllServer.h>

#include <Trkr_Clustering.C>
#include <Trkr_Reco.C>
#include <Trkr_RecoInit.C>
#include <Trkr_TpcReadoutInit.C>

#include <trackingdiagnostics/TrackResiduals.h>
#include </sphenix/user/mitrankova/ClusterPhase/ClusterPhaseAnalysis/ClusterPhaseAnalysis.h>

#include <phool/PHRandomSeed.h>
#include <phool/recoConsts.h>

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libffamodules.so)
R__LOAD_LIBRARY(libTrackingDiagnostics.so)
R__LOAD_LIBRARY(/sphenix/user/mitrankova/ClusterPhase/ClusterPhaseAnalysis/install/lib/libClusterPhaseAnalysis.so)


int Fun4All_G4_sPHENIX(
    const int nEvents = 2,
    const int processId = 0,
    const std::string outdir  = "/sphenix/user/mitrankova/ClusterPhase/",
    const string &outputFile = "Pi_G4sPHENIX_"
    )
{

    const int skip = 0;
 // std::string outdir  = "/sphenix/tg/tg01/hf/mitrankova/ZigZagSimulation/NoZigZag/";
  TString outfile = outdir + "output/"+ outputFile + to_string(nEvents) + "evts_" + to_string(processId);
  std::string theOutfile = outfile.Data();
  TString outfile_resid = outdir + "output_resid/"+ outputFile + to_string(nEvents) + "evts_" + to_string(processId);
  std::string theOutfile_resid = outfile_resid.Data();
  TString outfile_phase = outdir + "output_phase/"+ outputFile + to_string(nEvents) + "evts_" + to_string(processId);
  std::string theOutfile_phase = outfile_phase.Data();


  Fun4AllServer *se = Fun4AllServer::instance();
  se->Verbosity(2);


  PHRandomSeed::Verbosity(1);

  recoConsts *rc = recoConsts::instance();


  //===============
  // Input options
  //===============

  Input::VERBOSITY = 0;


   Input::SIMPLE = true;

   Input::UPSILON = false;
   Input::UPSILON_NUMBER = 1; // if you need 3 of them
   Input::UPSILON_VERBOSITY = 0;

  InputInit();

    if (Input::SIMPLE)
    {
      INPUTGENERATOR::SimpleEventGenerator[0]->add_particles("pi-", 1);

      INPUTGENERATOR::SimpleEventGenerator[0]->set_vertex_distribution_function(PHG4SimpleEventGenerator::Gaus,
                                                                              PHG4SimpleEventGenerator::Gaus,
                                                                              PHG4SimpleEventGenerator::Gaus);
      INPUTGENERATOR::SimpleEventGenerator[0]->set_vertex_distribution_mean(0., 0., 0.);
      INPUTGENERATOR::SimpleEventGenerator[0]->set_vertex_distribution_width(0.01, 0.01, 5.);
      
      INPUTGENERATOR::SimpleEventGenerator[0]->set_eta_range(0, 1);
      INPUTGENERATOR::SimpleEventGenerator[0]->set_phi_range(-M_PI, M_PI);
      INPUTGENERATOR::SimpleEventGenerator[0]->set_pt_range(3, 10);
      INPUTGENERATOR::SimpleEventGenerator[0]->set_power_law_n(8);
    }



    if (Input::UPSILON)
    {
      INPUTGENERATOR::VectorMesonGenerator[0]->add_decay_particles("e", 0);
      INPUTGENERATOR::VectorMesonGenerator[0]->set_rapidity_range(-1, 1);
      INPUTGENERATOR::VectorMesonGenerator[0]->set_pt_range(4., 10.);
      // Y species - select only one, last one wins
      INPUTGENERATOR::VectorMesonGenerator[0]->set_upsilon_1s();
    }

    InputRegister();


    rc->set_IntFlag("RUNNUMBER",1);

    SyncReco *sync = new SyncReco();
    se->registerSubsystem(sync);

    HeadReco *head = new HeadReco();
    se->registerSubsystem(head);
  

    FlagHandler *flag = new FlagHandler();
    se->registerSubsystem(flag);

 
  //======================
  // What to run
  //======================

  // QA, main switch
  Enable::QA = false;

 
  // central tracking
  Enable::MVTX = true;
  Enable::MVTX_CELL = Enable::MVTX && true;
  Enable::MVTX_CLUSTER = Enable::MVTX_CELL && true;


  Enable::INTT = true;
//  Enable::INTT_ABSORBER = true; // enables layerwise support structure readout
//  Enable::INTT_SUPPORT = true; // enable global support structure readout
  Enable::INTT_CELL = Enable::INTT && true;
  Enable::INTT_CLUSTER = Enable::INTT_CELL && true;

  Enable::TPC = true;
  Enable::TPC_ABSORBER = true;
  Enable::TPC_CELL = Enable::TPC && true;
  Enable::TPC_CLUSTER = Enable::TPC_CELL && true;


  Enable::MICROMEGAS = true;
  Enable::MICROMEGAS_CELL = Enable::MICROMEGAS && true;
  Enable::MICROMEGAS_CLUSTER = Enable::MICROMEGAS_CELL && true;


  Enable::TRACKING_TRACK = (Enable::MICROMEGAS_CLUSTER && Enable::TPC_CLUSTER && Enable::INTT_CLUSTER && Enable::MVTX_CLUSTER) && true;
  Enable::GLOBAL_RECO = (Enable::MBDFAKE || Enable::MBDRECO || Enable::TRACKING_TRACK) && true;
  Enable::TRACKING_EVAL = Enable::TRACKING_TRACK && Enable::GLOBAL_RECO && true;
  Enable::TRACKING_QA = Enable::TRACKING_TRACK && Enable::QA && true;

  // only do track matching if TRACKINGTRACK is also used
  Enable::TRACK_MATCHING = Enable::TRACKING_TRACK && true;
  Enable::TRACK_MATCHING_TREE = Enable::TRACK_MATCHING && true;
  Enable::TRACK_MATCHING_TREE_CLUSTERS = Enable::TRACK_MATCHING_TREE && true;


   TRACKING::pp_mode = true;
  // TRACKING::pp_extended_readout_time = 20000;

  // set flags to simulate and correct TPC distortions, specify distortion and correction files
  //G4TPC::ENABLE_STATIC_DISTORTIONS = true;
  //G4TPC::static_distortion_filename = std::string("/sphenix/user/rcorliss/distortion_maps/2023.02/Summary_hist_mdc2_UseFieldMaps_AA_event_0_bX180961051_0.distortion_map.hist.root");
  //G4TPC::ENABLE_STATIC_CORRECTIONS = true;
  //G4TPC::static_correction_filename = std::string("/sphenix/user/rcorliss/distortion_maps/2023.02/Summary_hist_mdc2_UseFieldMaps_AA_smoothed_average.correction_map.hist.root");
  //G4TPC::ENABLE_AVERAGE_CORRECTIONS = false;

  

  Enable::MAGNET = true;
  Enable::MAGNET_ABSORBER = true;




  Enable::BLACKHOLE = true;

  //===============
  // conditions DB flags
  //===============
  Enable::CDB = true;
  // global tag
  rc->set_StringFlag("CDB_GLOBALTAG",CDB::global_tag);
  // 64 bit timestamp
  rc->set_uint64Flag("TIMESTAMP",CDB::timestamp);

  //Enable::MVTX_APPLYMISALIGNMENT = true;
  //ACTSGEOM::mvtx_applymisalignment = Enable::MVTX_APPLYMISALIGNMENT; 

  G4Init();

  //---------------------
  // GEANT4 Detector description
  //---------------------

  G4Setup();
  

  //------------------
  // Detector Division
  //------------------


  if (Enable::MVTX_CELL) Mvtx_Cells();
  if (Enable::INTT_CELL) Intt_Cells();
  if (Enable::TPC_CELL) TPC_Cells();
  if (Enable::MICROMEGAS_CELL) Micromegas_Cells();


  //--------------
  // SVTX tracking
  //--------------
G4TPC::DO_HIT_ASSOCIATION = true; 
  TrackingInit();
    
  if (Enable::MVTX_CLUSTER) Mvtx_Clustering();
  if (Enable::INTT_CLUSTER) Intt_Clustering();
  if (Enable::TPC_CLUSTER)
    {
	  TPC_Clustering();
	}
    
  if (Enable::MICROMEGAS_CLUSTER) Micromegas_Clustering();


  
  Tracking_Reco();
  



  TString residoutfile = theOutfile_resid + "_resid.root";
  std::string residstring(residoutfile.Data());

  auto resid = new TrackResiduals("TrackResiduals");
  resid->outfileName(residstring);
  resid->alignment(false);
  resid->hitTree();
  resid->clusterTree();
  //resid->EventTree();
  resid->convertSeeds(G4TRACKING::convert_seeds_to_svtxtracks);
  resid->Verbosity(0);
  //se->registerSubsystem(resid);

std::cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<std::endl;
std::cout<<"G4TPC::TPC_GAS_MIXTURE = "<<G4TPC::TPC_GAS_MIXTURE<<std::endl;
std::cout<<"G4TPC::ArCF4Isobutane_drift_velocity = "<<G4TPC::ArCF4Isobutane_drift_velocity<<std::endl;
std::cout<<"diffusion_long = "<<G4TPC::ArCF4Isobutane_diffusion_long<<std::endl;
std::cout<<"diffusion_trans = "<<G4TPC::ArCF4Isobutane_diffusion_trans<<std::endl;



  TString phaseoutfile = theOutfile_phase + "_phase.root";
  std::string phasestring(phaseoutfile.Data());
  ClusterPhaseAnalysis *ana = new ClusterPhaseAnalysis("ClusterPhaseAnalysis", phasestring);
  ana->Verbosity(10);
  ana->setIsSimulation(true);
  //ana->set_cluster_container_name("TRKR_CLUSTER");  // Set cluster container name if different from default
  se->registerSubsystem(ana);




  //--------------
  // Set up Input Managers
  //--------------

  InputManagers();


  //-----------------
  // Event processing
  //-----------------
  if (Enable::DISPLAY)
  {
    DisplayOn();

    gROOT->ProcessLine("Fun4AllServer *se = Fun4AllServer::instance();");
    gROOT->ProcessLine("PHG4Reco *g4 = (PHG4Reco *) se->getSubsysReco(\"PHG4RECO\");");

    cout << "-------------------------------------------------" << endl;
    cout << "You are in event display mode. Run one event with" << endl;
    cout << "se->run(1)" << endl;
    cout << "Run Geant4 command with following examples" << endl;
    gROOT->ProcessLine("displaycmd()");

    return 0;
  }


  if (nEvents <= 0)
  {
    return 0;
  }


  se->skip(skip);
  se->run(nEvents);
  se->PrintTimer();



  //-----
  // Exit
  //-----

  CDBInterface::instance()->Print(); // print used DB files
  se->End();
  std::cout << "All done" << std::endl;
  delete se;
  if (Enable::PRODUCTION)
  {
    Production_MoveOutput();
  }

  gSystem->Exit(0);
  return 0;
}
#endif