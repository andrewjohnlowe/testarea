// the first two come for free when using AthAlgorithm
//#include "GaudiKernel/MsgStream.h"
//#include "StoreGate/StoreGateSvc.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IToolSvc.h"

#include "JetEvent/JetCollection.h"
#include "JetEvent/Jet.h"
#include "JetEvent/JetConstituentIterator.h"
#include "JetUtils/JetCollectionHelper.h"
#include "JetUtils/JetSorters.h"
#include "JetTagInfo/TruthInfo.h"
#include "JetTagInfo/BaseTagInfo.h"

// #include "NavFourMom/IParticleContainer.h"
// #include "NavFourMom/INavigable4MomentumCollection.h"

#include "GaudiKernel/ITHistSvc.h"
#include "GaudiKernel/Property.h"
#include "AthenaKernel/errorcheck.h"
#include "TTree.h"

#include "UserAnalysis/MyAnalysis.h"
#include "UserAnalysis/Kludges.h"
#include "UserAnalysis/Numeric.h"
#include "UserAnalysis/JetFlags.h" // Common jet reco flag access

#include "UserAnalysis/JetShapes.h"
#include "UserAnalysis/EventShapes.h"
#include "UserAnalysis/SphericitySpherocity.h"
//#include "UserAnalysis/AngularCorrelation.h"
#include "UserAnalysis/FoxWolfram.h"
#include "UserAnalysis/RadiationVariables.h"

#include "McParticleEvent/TruthParticleContainer.h"
#include "VxVertex/VxContainer.h"

// Event Info
#include "EventInfo/EventInfo.h"
#include "EventInfo/EventID.h"
#include "EventInfo/EventType.h"

// #include <algorithm>
// #include <math.h>
// #include <functional>

#include <string>
#include <tr1/unordered_map>

//#include <fenv.h> // Turn off FPE crashes in execute method

//#include "valgrind/callgrind.h"

const int bcfVersion = 0; // 1: Euclidian distance, 2: zero, 0: infinite line
const int maxTermsFW = 1000; // Max number of terms in Fox-Wolfram summations

//#define VERBOSE WARNING

// // debugging macros so we can pin down message provenance at a glance
#include <iostream>
// #define DEBUG(x)							\
//   std::cout << "DEBUG:\t" << "\t" << #x << "\t" << x << "\t"		\
//   << __FUNCTION__ << "\t" << __FILE__ << "\t" << __LINE__ << std::endl;
//#define DEBUG(x)
#define DEBUG(x)							\
  ATH_MSG_DEBUG("DEBUG:\t" << "\t" << #x << "\t" << x << "\t"		\
		<< __FUNCTION__ << "\t" << __FILE__ << "\t" << __LINE__);
#define DEBUG(x)
//#define OUTPUT_CONSTITUENTS

using namespace Analysis;

MyAnalysis::MyAnalysis(const std::string& name, ISvcLocator* pSvcLocator):
AthAlgorithm(name, pSvcLocator),
m_jvfTool("JetVertexAssociationTool/JetVertexAssociationTool", this),
m_jetcollname("JetCollection"),
m_cutIP3DSV1(1.55),
m_vxCandidatesName("VxPrimaryCandidate") {
  // switches to control the analysis through job options
  declareProperty("JVFTool", m_jvfTool);
  declareProperty("JetCollection", m_jetcollname);
  declareProperty("JetCollections", m_jetcollnames);
  declareProperty("BjetWeightCut", m_cutIP3DSV1);
  declareProperty("VxContainerName", m_vxCandidatesName);
  declareProperty("McParticleContainer", m_truthParticleContainerName = "SpclMC");
}
MyAnalysis::~MyAnalysis() {}

StatusCode MyAnalysis::beginRun() {

  ATH_MSG_INFO("Initializing MyAnalysis (before eventloop)");

  if(m_jvfTool.retrieve().isFailure()) {
    ATH_MSG_ERROR("Creation of algTool JetVertexAssociationTool FAILED!");
    return StatusCode::FAILURE;
  }
  
  ATH_MSG_DEBUG("Will use JetCollections: ");
  for(unsigned int i = 0; i < m_jetcollnames.size(); ++i) {
    ATH_MSG_DEBUG(m_jetcollnames[i]); 
  }
  
  return StatusCode::SUCCESS;
} 

StatusCode MyAnalysis::initialize() {

  ATH_MSG_INFO("Initializing MyAnalysis");

  if(m_jvfTool.retrieve().isFailure()) {
    ATH_MSG_ERROR("Creation of algTool JetVertexAssociationTool FAILED!");
    return StatusCode::FAILURE;
  }
  
  ATH_MSG_DEBUG("Will use JetCollections: ");
  for(unsigned int i = 0; i < m_jetcollnames.size(); ++i) {
    ATH_MSG_DEBUG(m_jetcollnames[i]); 
  }
  
  /** get a handle on the NTuple and histogramming service */
  StatusCode sc = service("THistSvc", m_thistSvc);
  if(sc.isFailure()) {
    ATH_MSG_ERROR("Unable to retrieve pointer to THistSvc");
    return sc;
  }

  _jStr_numPVJ = new std::vector<int>;
  _jStr_PVxJ = new std::vector<double>;
  _jStr_PVyJ = new std::vector<double>;
  _jStr_PVzJ = new std::vector<double>;
  _jStr_PVrJ = new std::vector<double>;
  _jStr_PVchiSqrJ = new std::vector<double>;
  _jStr_PVnDoFJ = new std::vector<double>;
  _jStr_PVfitJ = new std::vector<double>;
  _jStr_PVnTrkJ = new std::vector<int>;
  _jStr_PVtypeJ = new std::vector<int>;

  // EventShape, whole event, jets as input:
  _jStr_evtEJ_etaJ1 = new std::vector<double>;
  _jStr_evtEJ_etaJ2 = new std::vector<double>;
  _jStr_evtEJ_etaJ3 = new std::vector<double>;
  _jStr_evtEJ_etaJ4 = new std::vector<double>;
  _jStr_evtEJ_pTJ1 = new std::vector<double>;
  _jStr_evtEJ_pTJ2 = new std::vector<double>;
  _jStr_evtEJ_pTJ3 = new std::vector<double>;
  _jStr_evtEJ_pTJ4 = new std::vector<double>;

  _jStr_evtEJ_fE1 = new std::vector<double>;
  _jStr_evtEJ_fE2 = new std::vector<double>;
  _jStr_evtEJ_fE3 = new std::vector<double>;
  _jStr_evtEJ_fET1 = new std::vector<double>;
  _jStr_evtEJ_fET2 = new std::vector<double>;
  _jStr_evtEJ_fET3 = new std::vector<double>;
  _jStr_evtEJ_mct = new std::vector<double>;
  _jStr_evtEJ_q = new std::vector<double>;
  _jStr_evtEJ_mjj = new std::vector<double>;
  _jStr_evtEJ_mTjj = new std::vector<double>;
  _jStr_evtEJ_mjjj = new std::vector<double>;
  _jStr_evtEJ_mTjjj = new std::vector<double>;
  _jStr_evtEJ_mjjjj = new std::vector<double>;
  _jStr_evtEJ_mTjjjj = new std::vector<double>;
  _jStr_evtEJ_zetaPlus = new std::vector<double>;
  _jStr_evtEJ_zetaMinus = new std::vector<double>;
  _jStr_evtEJ_B12 = new std::vector<double>;
  _jStr_evtEJ_BT12 = new std::vector<double>;
  _jStr_evtEJ_Alpha = new std::vector<double>;
  _jStr_evtEJ_AlphaT = new std::vector<double>;
  _jStr_evtEJ_massDemocracy = new std::vector<double>;
  _jStr_evtEJ_betaflow = new std::vector<double>;
  _jStr_evtEJ_betaflow_1GeV = new std::vector<double>;
  _jStr_evtEJ_betaflow_5GeV = new std::vector<double>;
  _jStr_evtEJ_y23 = new std::vector<double>;
  _jStr_evtEJ_y23_1GeV = new std::vector<double>;
  _jStr_evtEJ_y23_5GeV = new std::vector<double>;
  _jStr_evtEJ_lny23 = new std::vector<double>;
  _jStr_evtEJ_lny23_1GeV = new std::vector<double>;
  _jStr_evtEJ_lny23_5GeV = new std::vector<double>;
  _jStr_evtEJ_theta = new std::vector<double>;
  _jStr_evtEJ_asym = new std::vector<double>;
  _jStr_evtEJ_yB = new std::vector<double>;
  _jStr_evtEJ_yStar = new std::vector<double>;
  _jStr_evtEJ_thetaStar = new std::vector<double>;
  _jStr_evtEJ_chi = new std::vector<double>;
  _jStr_evtEJ_deltaPhiJJ = new std::vector<double>;
  _jStr_evtEJ_deltaThetaJJ = new std::vector<double>;
  _jStr_evtEJ_deltaEtaJJ = new std::vector<double>;
  _jStr_evtEJ_deltaRapidityJJ = new std::vector<double>;
  _jStr_evtEJ_deltaRJJ = new std::vector<double>;
  _jStr_evtEJ_deltaRJJY = new std::vector<double>;
  _jStr_evtEJ_sigmaPhiJJ = new std::vector<double>;
  _jStr_evtEJ_sigmaThetaJJ = new std::vector<double>;
  _jStr_evtEJ_sigmaEtaJJ = new std::vector<double>;
  _jStr_evtEJ_sigmaRapidityJJ = new std::vector<double>;
  _jStr_evtEJ_sigmaPtJJ = new std::vector<double>;
  _jStr_evtEJ_sigmaEtJJ = new std::vector<double>;
  _jStr_evtEJ_sigmaEt12 = new std::vector<double>;
  _jStr_evtEJ_sigmaEt34 = new std::vector<double>;
  _jStr_evtEJ_A234 = new std::vector<double>;
  _jStr_evtEJ_asymPhiJJ = new std::vector<double>;
  _jStr_evtEJ_asymThetaJJ = new std::vector<double>;
  _jStr_evtEJ_asymEtaJJ = new std::vector<double>;
  _jStr_evtEJ_asymRapidityJJ = new std::vector<double>;
  _jStr_evtEJ_acoplanarity = new std::vector<double>;
  _jStr_evtEJ_twist = new std::vector<double>;
  _jStr_evtEJ_twistY = new std::vector<double>;
  _jStr_evtEJ_jetSumE = new std::vector<double>;
  _jStr_evtEJ_jetSumET = new std::vector<double>;
  _jStr_evtEJ_jetSumPT = new std::vector<double>;
  _jStr_evtEJ_jetSumM = new std::vector<double>;
  _jStr_evtEJ_jetSumMT = new std::vector<double>;
  _jStr_evtEJ_HTprime = new std::vector<double>;
  _jStr_evtEJ_centrality = new std::vector<double>;
  _jStr_evtEJ_centralityP = new std::vector<double>;
  _jStr_evtEJ_zminJ1J2 = new std::vector<double>;
  _jStr_evtEJ_zmaxJ1J2 = new std::vector<double>;
  _jStr_evtEJ_zminAllJets = new std::vector<double>;
  _jStr_evtEJ_zmaxAllJets = new std::vector<double>;
  _jStr_evtEJ_zminJ1J2Phi = new std::vector<double>;
  _jStr_evtEJ_zminJ1J2Theta = new std::vector<double>;
  _jStr_evtEJ_zminJ1J2Eta = new std::vector<double>;
  _jStr_evtEJ_zminJ1J2Rapidity = new std::vector<double>;
  _jStr_evtEJ_cosHelicityJ1 = new std::vector<double>;
  _jStr_evtEJ_helicityJ1 = new std::vector<double>;
  _jStr_evtEJ_azilicityJ1 = new std::vector<double>;
  _jStr_evtEJ_cosHelicityJ2 = new std::vector<double>;
  _jStr_evtEJ_helicityJ2 = new std::vector<double>;
  _jStr_evtEJ_azilicityJ2 = new std::vector<double>;
  _jStr_evtEJ_cosThetaJ1 = new std::vector<double>;
  _jStr_evtEJ_cosThetaJ2 = new std::vector<double>;
  _jStr_evtEJ_deltaRapidityXtoJ1CM = new std::vector<double>;
  _jStr_evtEJ_deltaRapidityXtoJ2CM = new std::vector<double>;
  _jStr_evtEJ_deltaRapidityXtoJ1 = new std::vector<double>;
  _jStr_evtEJ_deltaRapidityXtoJ2 = new std::vector<double>;
  _jStr_evtEJ_cosTheta1 = new std::vector<double>;
  _jStr_evtEJ_cosTheta2 = new std::vector<double>;
  _jStr_evtEJ_cosThetaStar1 = new std::vector<double>;
  _jStr_evtEJ_cosThetaStar2 = new std::vector<double>;
  _jStr_evtEJ_cosPhiTilde1 = new std::vector<double>;
  _jStr_evtEJ_cosPhiTilde2 = new std::vector<double>;
  _jStr_evtEJ_cosPhi = new std::vector<double>;
  _jStr_evtEJ_PhiTilde1 = new std::vector<double>;
  _jStr_evtEJ_PhiTilde2 = new std::vector<double>;
  _jStr_evtEJ_Phi = new std::vector<double>;
  _jStr_evtEJ_M = new std::vector<double>;
  _jStr_evtEJ_DeltaM = new std::vector<double>;
  _jStr_evtEJ_asymM = new std::vector<double>;
  _jStr_evtEJ_Q = new std::vector<double>;
  _jStr_evtEJ_SCDF = new std::vector<double>;
  _jStr_evtEJ_dPhiIJCDF = new std::vector<double>;
  _jStr_evtEJ_dPhiKLCDF = new std::vector<double>;
  _jStr_evtEJ_PtIPtJCDF = new std::vector<double>;
  _jStr_evtEJ_PtKPtLCDF = new std::vector<double>;
  _jStr_evtEJ_SD0 = new std::vector<double>;
  _jStr_evtEJ_dPhiIJD0 = new std::vector<double>;
  _jStr_evtEJ_dPhiKLD0 = new std::vector<double>;
  _jStr_evtEJ_PtIPtJD0 = new std::vector<double>;
  _jStr_evtEJ_PtKPtLD0 = new std::vector<double>;
  _jStr_evtEJ_DeltaSCDF = new std::vector<double>;
  _jStr_evtEJ_DeltaSD0 = new std::vector<double>;
  _jStr_evtEJ_Delta12 = new std::vector<double>;
  _jStr_evtEJ_Delta12n = new std::vector<double>;
  _jStr_evtEJ_M14 = new std::vector<double>;
  _jStr_evtEJ_y1y2 = new std::vector<double>;

  _jStr_evtEJ_Mprime4 = new std::vector<double>;
  _jStr_evtEJ_dMprime4 = new std::vector<double>;
  _jStr_evtEJ_MprimeAvg4 = new std::vector<double>;
  _jStr_evtEJ_DeltaMin4 = new std::vector<double>;
  _jStr_evtEJ_DeltaMax4 = new std::vector<double>;
  _jStr_evtEJ_DeltaPhiXX4 = new std::vector<double>;
  _jStr_evtEJ_DeltaYXX4 = new std::vector<double>;
  _jStr_evtEJ_DeltaRXX4 = new std::vector<double>;
  _jStr_evtEJ_TwistXX4 = new std::vector<double>;
  _jStr_evtEJ_separated4 = new std::vector<int>;
  _jStr_evtEJ_Mprime = new std::vector<double>;
  _jStr_evtEJ_dMprime = new std::vector<double>;
  _jStr_evtEJ_MprimeAvg = new std::vector<double>;
  _jStr_evtEJ_DeltaMin = new std::vector<double>;
  _jStr_evtEJ_DeltaMax = new std::vector<double>;
  _jStr_evtEJ_DeltaPhiXX = new std::vector<double>;
  _jStr_evtEJ_DeltaYXX = new std::vector<double>;
  _jStr_evtEJ_DeltaRXX = new std::vector<double>;
  _jStr_evtEJ_TwistXX = new std::vector<double>;
  _jStr_evtEJ_separated = new std::vector<int>;
  
  // EventShape, mini event (jet in own CM frame), subjets as input:
  _jStr_evtJj_etaJ1 = new std::vector<double>;
  _jStr_evtJj_etaJ2 = new std::vector<double>;
  _jStr_evtJj_etaJ3 = new std::vector<double>;
  _jStr_evtJj_etaJ4 = new std::vector<double>;
  _jStr_evtJj_pTJ1 = new std::vector<double>;
  _jStr_evtJj_pTJ2 = new std::vector<double>;
  _jStr_evtJj_pTJ3 = new std::vector<double>;
  _jStr_evtJj_pTJ4 = new std::vector<double>;

  _jStr_evtJj_fE1 = new std::vector<double>;
  _jStr_evtJj_fE2 = new std::vector<double>;
  _jStr_evtJj_fE3 = new std::vector<double>;
  _jStr_evtJj_fET1 = new std::vector<double>;
  _jStr_evtJj_fET2 = new std::vector<double>;
  _jStr_evtJj_fET3 = new std::vector<double>;
  _jStr_evtJj_mct = new std::vector<double>;
  _jStr_evtJj_q = new std::vector<double>;
  _jStr_evtJj_mjj = new std::vector<double>;
  _jStr_evtJj_mTjj = new std::vector<double>;
  _jStr_evtJj_mjjj = new std::vector<double>;
  _jStr_evtJj_mTjjj = new std::vector<double>;
  _jStr_evtJj_mjjjj = new std::vector<double>;
  _jStr_evtJj_mTjjjj = new std::vector<double>;
  _jStr_evtJj_zetaPlus = new std::vector<double>;
  _jStr_evtJj_zetaMinus = new std::vector<double>;
  _jStr_evtJj_B12 = new std::vector<double>;
  _jStr_evtJj_BT12 = new std::vector<double>;
  _jStr_evtJj_Alpha = new std::vector<double>;
  _jStr_evtJj_AlphaT = new std::vector<double>;
  _jStr_evtJj_massDemocracy = new std::vector<double>;
  _jStr_evtJj_betaflow = new std::vector<double>;
  _jStr_evtJj_betaflow_1GeV = new std::vector<double>;
  _jStr_evtJj_betaflow_5GeV = new std::vector<double>;
  _jStr_evtJj_y23 = new std::vector<double>;
  _jStr_evtJj_y23_1GeV = new std::vector<double>;
  _jStr_evtJj_y23_5GeV = new std::vector<double>;
  _jStr_evtJj_lny23 = new std::vector<double>;
  _jStr_evtJj_lny23_1GeV = new std::vector<double>;
  _jStr_evtJj_lny23_5GeV = new std::vector<double>;
  _jStr_evtJj_theta = new std::vector<double>;
  _jStr_evtJj_asym = new std::vector<double>;
  _jStr_evtJj_yB = new std::vector<double>;
  _jStr_evtJj_yStar = new std::vector<double>;
  _jStr_evtJj_thetaStar = new std::vector<double>;
  _jStr_evtJj_chi = new std::vector<double>;
  _jStr_evtJj_deltaPhiJJ = new std::vector<double>;
  _jStr_evtJj_deltaThetaJJ = new std::vector<double>;
  _jStr_evtJj_deltaEtaJJ = new std::vector<double>;
  _jStr_evtJj_deltaRapidityJJ = new std::vector<double>;
  _jStr_evtJj_deltaRJJ = new std::vector<double>;
  _jStr_evtJj_deltaRJJY = new std::vector<double>;
  _jStr_evtJj_sigmaPhiJJ = new std::vector<double>;
  _jStr_evtJj_sigmaThetaJJ = new std::vector<double>;
  _jStr_evtJj_sigmaEtaJJ = new std::vector<double>;
  _jStr_evtJj_sigmaRapidityJJ = new std::vector<double>;
  _jStr_evtJj_sigmaPtJJ = new std::vector<double>;
  _jStr_evtJj_sigmaEtJJ = new std::vector<double>;
  _jStr_evtJj_sigmaEt12 = new std::vector<double>;
  _jStr_evtJj_sigmaEt34 = new std::vector<double>;
  _jStr_evtJj_A234 = new std::vector<double>;
  _jStr_evtJj_asymPhiJJ = new std::vector<double>;
  _jStr_evtJj_asymThetaJJ = new std::vector<double>;
  _jStr_evtJj_asymEtaJJ = new std::vector<double>;
  _jStr_evtJj_asymRapidityJJ = new std::vector<double>;
  _jStr_evtJj_acoplanarity = new std::vector<double>;
  _jStr_evtJj_twist = new std::vector<double>;
  _jStr_evtJj_twistY = new std::vector<double>;
  _jStr_evtJj_jetSumE = new std::vector<double>;
  _jStr_evtJj_jetSumET = new std::vector<double>;
  _jStr_evtJj_jetSumPT = new std::vector<double>;
  _jStr_evtJj_jetSumM = new std::vector<double>;
  _jStr_evtJj_jetSumMT = new std::vector<double>;
  _jStr_evtJj_HTprime = new std::vector<double>;
  _jStr_evtJj_centrality = new std::vector<double>;
  _jStr_evtJj_centralityP = new std::vector<double>;
  _jStr_evtJj_zminJ1J2 = new std::vector<double>;
  _jStr_evtJj_zmaxJ1J2 = new std::vector<double>;
  _jStr_evtJj_zminAllJets = new std::vector<double>;
  _jStr_evtJj_zmaxAllJets = new std::vector<double>;
  _jStr_evtJj_zminJ1J2Phi = new std::vector<double>;
  _jStr_evtJj_zminJ1J2Theta = new std::vector<double>;
  _jStr_evtJj_zminJ1J2Eta = new std::vector<double>;
  _jStr_evtJj_zminJ1J2Rapidity = new std::vector<double>;
  _jStr_evtJj_cosHelicityJ1 = new std::vector<double>;
  _jStr_evtJj_helicityJ1 = new std::vector<double>;
  _jStr_evtJj_azilicityJ1 = new std::vector<double>;
  _jStr_evtJj_cosHelicityJ2 = new std::vector<double>;
  _jStr_evtJj_helicityJ2 = new std::vector<double>;
  _jStr_evtJj_azilicityJ2 = new std::vector<double>;
  _jStr_evtJj_cosThetaJ1 = new std::vector<double>;
  _jStr_evtJj_cosThetaJ2 = new std::vector<double>;
  _jStr_evtJj_deltaRapidityXtoJ1CM = new std::vector<double>;
  _jStr_evtJj_deltaRapidityXtoJ2CM = new std::vector<double>;
  _jStr_evtJj_deltaRapidityXtoJ1 = new std::vector<double>;
  _jStr_evtJj_deltaRapidityXtoJ2 = new std::vector<double>;
  _jStr_evtJj_cosTheta1 = new std::vector<double>;
  _jStr_evtJj_cosTheta2 = new std::vector<double>;
  _jStr_evtJj_cosThetaStar1 = new std::vector<double>;
  _jStr_evtJj_cosThetaStar2 = new std::vector<double>;
  _jStr_evtJj_cosPhiTilde1 = new std::vector<double>;
  _jStr_evtJj_cosPhiTilde2 = new std::vector<double>;
  _jStr_evtJj_cosPhi = new std::vector<double>;
  _jStr_evtJj_PhiTilde1 = new std::vector<double>;
  _jStr_evtJj_PhiTilde2 = new std::vector<double>;
  _jStr_evtJj_Phi = new std::vector<double>;
  _jStr_evtJj_M = new std::vector<double>;
  _jStr_evtJj_DeltaM = new std::vector<double>;
  _jStr_evtJj_asymM = new std::vector<double>;
  _jStr_evtJj_Q = new std::vector<double>;
  _jStr_evtJj_SCDF = new std::vector<double>;
  _jStr_evtJj_dPhiIJCDF = new std::vector<double>;
  _jStr_evtJj_dPhiKLCDF = new std::vector<double>;
  _jStr_evtJj_PtIPtJCDF = new std::vector<double>;
  _jStr_evtJj_PtKPtLCDF = new std::vector<double>;
  _jStr_evtJj_SD0 = new std::vector<double>;
  _jStr_evtJj_dPhiIJD0 = new std::vector<double>;
  _jStr_evtJj_dPhiKLD0 = new std::vector<double>;
  _jStr_evtJj_PtIPtJD0 = new std::vector<double>;
  _jStr_evtJj_PtKPtLD0 = new std::vector<double>;
  _jStr_evtJj_DeltaSCDF = new std::vector<double>;
  _jStr_evtJj_DeltaSD0 = new std::vector<double>;
  _jStr_evtJj_Delta12 = new std::vector<double>;
  _jStr_evtJj_Delta12n = new std::vector<double>;
  _jStr_evtJj_M14 = new std::vector<double>;
  _jStr_evtJj_y1y2 = new std::vector<double>;

  _jStr_evtJj_Mprime4 = new std::vector<double>;
  _jStr_evtJj_dMprime4 = new std::vector<double>;
  _jStr_evtJj_MprimeAvg4 = new std::vector<double>;
  _jStr_evtJj_DeltaMin4 = new std::vector<double>;
  _jStr_evtJj_DeltaMax4 = new std::vector<double>;
  _jStr_evtJj_DeltaPhiXX4 = new std::vector<double>;
  _jStr_evtJj_DeltaYXX4 = new std::vector<double>;
  _jStr_evtJj_DeltaRXX4 = new std::vector<double>;
  _jStr_evtJj_TwistXX4 = new std::vector<double>;
  _jStr_evtJj_separated4 = new std::vector<int>;
  _jStr_evtJj_Mprime = new std::vector<double>;
  _jStr_evtJj_dMprime = new std::vector<double>;
  _jStr_evtJj_MprimeAvg = new std::vector<double>;
  _jStr_evtJj_DeltaMin = new std::vector<double>;
  _jStr_evtJj_DeltaMax = new std::vector<double>;
  _jStr_evtJj_DeltaPhiXX = new std::vector<double>;
  _jStr_evtJj_DeltaYXX = new std::vector<double>;
  _jStr_evtJj_DeltaRXX = new std::vector<double>;
  _jStr_evtJj_TwistXX = new std::vector<double>;
  _jStr_evtJj_separated = new std::vector<int>;

  _jStr_fwmEJ_xiPlus = new std::vector<double>;
  _jStr_fwmEJ_xiMinus = new std::vector<double>;

  _jStr_fwmEJ_xPlus = new std::vector<double>;
  _jStr_fwmEJ_xMinus = new std::vector<double>;

  // Fox-Wolfram, whole event, jets as input:
  _jStr_fwmEJ_Psi1 = new std::vector<double>;
  _jStr_fwmEJ_B0 = new std::vector<double>;
  _jStr_fwmEJ_B1 = new std::vector<double>;
  _jStr_fwmEJ_B2 = new std::vector<double>;
  _jStr_fwmEJ_B3 = new std::vector<double>;
  _jStr_fwmEJ_B4 = new std::vector<double>;
  _jStr_fwmEJ_B5 = new std::vector<double>;
  _jStr_fwmEJ_B6 = new std::vector<double>;
  _jStr_fwmEJ_B7 = new std::vector<double>;
  _jStr_fwmEJ_B8 = new std::vector<double>;
  _jStr_fwmEJ_C0 = new std::vector<double>;
  _jStr_fwmEJ_C1 = new std::vector<double>;
  _jStr_fwmEJ_C2 = new std::vector<double>;
  _jStr_fwmEJ_C3 = new std::vector<double>;
  _jStr_fwmEJ_C4 = new std::vector<double>;
  _jStr_fwmEJ_C5 = new std::vector<double>;
  _jStr_fwmEJ_C6 = new std::vector<double>;
  _jStr_fwmEJ_C7 = new std::vector<double>;
  _jStr_fwmEJ_C8 = new std::vector<double>;
  _jStr_fwmEJ_K0 = new std::vector<double>;
  _jStr_fwmEJ_K1 = new std::vector<double>;
  _jStr_fwmEJ_K2 = new std::vector<double>;
  _jStr_fwmEJ_K3 = new std::vector<double>;
  _jStr_fwmEJ_K4 = new std::vector<double>;
  _jStr_fwmEJ_K5 = new std::vector<double>;
  _jStr_fwmEJ_K6 = new std::vector<double>;
  _jStr_fwmEJ_K7 = new std::vector<double>;
  _jStr_fwmEJ_K8 = new std::vector<double>;
  _jStr_fwmEJ_D0 = new std::vector<double>;
  _jStr_fwmEJ_D1 = new std::vector<double>;
  _jStr_fwmEJ_D2 = new std::vector<double>;
  _jStr_fwmEJ_D3 = new std::vector<double>;
  _jStr_fwmEJ_D4 = new std::vector<double>;
  _jStr_fwmEJ_D5 = new std::vector<double>;
  _jStr_fwmEJ_D6 = new std::vector<double>;
  _jStr_fwmEJ_D7 = new std::vector<double>;
  _jStr_fwmEJ_D8 = new std::vector<double>;
  _jStr_fwmEJ_H0 = new std::vector<double>;
  _jStr_fwmEJ_H1 = new std::vector<double>;
  _jStr_fwmEJ_H2 = new std::vector<double>;
  _jStr_fwmEJ_H3 = new std::vector<double>;
  _jStr_fwmEJ_H4 = new std::vector<double>;
  _jStr_fwmEJ_H5 = new std::vector<double>;
  _jStr_fwmEJ_H6 = new std::vector<double>;
  _jStr_fwmEJ_H7 = new std::vector<double>;
  _jStr_fwmEJ_H8 = new std::vector<double>;
  _jStr_fwmEJ_Q0 = new std::vector<double>;
  _jStr_fwmEJ_Q1 = new std::vector<double>;
  _jStr_fwmEJ_Q2 = new std::vector<double>;
  _jStr_fwmEJ_Q3 = new std::vector<double>;
  _jStr_fwmEJ_Q4 = new std::vector<double>;
  _jStr_fwmEJ_Q5 = new std::vector<double>;
  _jStr_fwmEJ_Q6 = new std::vector<double>;
  _jStr_fwmEJ_Q7 = new std::vector<double>;
  _jStr_fwmEJ_Q8 = new std::vector<double>;
  _jStr_fwmEJ_Pi1 = new std::vector<double>;
  _jStr_fwmEJ_Pi2 = new std::vector<double>;
  _jStr_fwmEJ_Pi3 = new std::vector<double>;
  _jStr_fwmEJ_Pi4 = new std::vector<double>;
  _jStr_fwmEJ_B10 = new std::vector<double>;
  _jStr_fwmEJ_B20 = new std::vector<double>;
  _jStr_fwmEJ_B30 = new std::vector<double>;
  _jStr_fwmEJ_B40 = new std::vector<double>;
  _jStr_fwmEJ_B50 = new std::vector<double>;
  _jStr_fwmEJ_B60 = new std::vector<double>;
  _jStr_fwmEJ_B70 = new std::vector<double>;
  _jStr_fwmEJ_B80 = new std::vector<double>;
  _jStr_fwmEJ_C10 = new std::vector<double>;
  _jStr_fwmEJ_C20 = new std::vector<double>;
  _jStr_fwmEJ_C30 = new std::vector<double>;
  _jStr_fwmEJ_C40 = new std::vector<double>;
  _jStr_fwmEJ_C50 = new std::vector<double>;
  _jStr_fwmEJ_C60 = new std::vector<double>;
  _jStr_fwmEJ_C70 = new std::vector<double>;
  _jStr_fwmEJ_C80 = new std::vector<double>;
  _jStr_fwmEJ_K10 = new std::vector<double>;
  _jStr_fwmEJ_K20 = new std::vector<double>;
  _jStr_fwmEJ_K30 = new std::vector<double>;
  _jStr_fwmEJ_K40 = new std::vector<double>;
  _jStr_fwmEJ_K50 = new std::vector<double>;
  _jStr_fwmEJ_K60 = new std::vector<double>;
  _jStr_fwmEJ_K70 = new std::vector<double>;
  _jStr_fwmEJ_K80 = new std::vector<double>;
  _jStr_fwmEJ_D10 = new std::vector<double>;
  _jStr_fwmEJ_D20 = new std::vector<double>;
  _jStr_fwmEJ_D30 = new std::vector<double>;
  _jStr_fwmEJ_D40 = new std::vector<double>;
  _jStr_fwmEJ_D50 = new std::vector<double>;
  _jStr_fwmEJ_D60 = new std::vector<double>;
  _jStr_fwmEJ_D70 = new std::vector<double>;
  _jStr_fwmEJ_D80 = new std::vector<double>;
  _jStr_fwmEJ_H10 = new std::vector<double>;
  _jStr_fwmEJ_H20 = new std::vector<double>;
  _jStr_fwmEJ_H30 = new std::vector<double>;
  _jStr_fwmEJ_H40 = new std::vector<double>;
  _jStr_fwmEJ_H50 = new std::vector<double>;
  _jStr_fwmEJ_H60 = new std::vector<double>;
  _jStr_fwmEJ_H70 = new std::vector<double>;
  _jStr_fwmEJ_H80 = new std::vector<double>;
  _jStr_fwmEJ_Q10 = new std::vector<double>;
  _jStr_fwmEJ_Q20 = new std::vector<double>;
  _jStr_fwmEJ_Q30 = new std::vector<double>;
  _jStr_fwmEJ_Q40 = new std::vector<double>;
  _jStr_fwmEJ_Q50 = new std::vector<double>;
  _jStr_fwmEJ_Q60 = new std::vector<double>;
  _jStr_fwmEJ_Q70 = new std::vector<double>;
  _jStr_fwmEJ_Q80 = new std::vector<double>;

  // // Fox-Wolfram, whole event, jet constituents as input:
  // _jStr_fwmEi_Psi1 = new std::vector<double>;
  // _jStr_fwmEi_B0 = new std::vector<double>;
  // _jStr_fwmEi_B1 = new std::vector<double>;
  // _jStr_fwmEi_B2 = new std::vector<double>;
  // _jStr_fwmEi_B3 = new std::vector<double>;
  // _jStr_fwmEi_B4 = new std::vector<double>;
  // _jStr_fwmEi_B5 = new std::vector<double>;
  // _jStr_fwmEi_B6 = new std::vector<double>;
  // _jStr_fwmEi_B7 = new std::vector<double>;
  // _jStr_fwmEi_B8 = new std::vector<double>;
  // _jStr_fwmEi_C0 = new std::vector<double>;
  // _jStr_fwmEi_C1 = new std::vector<double>;
  // _jStr_fwmEi_C2 = new std::vector<double>;
  // _jStr_fwmEi_C3 = new std::vector<double>;
  // _jStr_fwmEi_C4 = new std::vector<double>;
  // _jStr_fwmEi_C5 = new std::vector<double>;
  // _jStr_fwmEi_C6 = new std::vector<double>;
  // _jStr_fwmEi_C7 = new std::vector<double>;
  // _jStr_fwmEi_C8 = new std::vector<double>;
  // _jStr_fwmEi_K0 = new std::vector<double>;
  // _jStr_fwmEi_K1 = new std::vector<double>;
  // _jStr_fwmEi_K2 = new std::vector<double>;
  // _jStr_fwmEi_K3 = new std::vector<double>;
  // _jStr_fwmEi_K4 = new std::vector<double>;
  // _jStr_fwmEi_K5 = new std::vector<double>;
  // _jStr_fwmEi_K6 = new std::vector<double>;
  // _jStr_fwmEi_K7 = new std::vector<double>;
  // _jStr_fwmEi_K8 = new std::vector<double>;
  // _jStr_fwmEi_D0 = new std::vector<double>;
  // _jStr_fwmEi_D1 = new std::vector<double>;
  // _jStr_fwmEi_D2 = new std::vector<double>;
  // _jStr_fwmEi_D3 = new std::vector<double>;
  // _jStr_fwmEi_D4 = new std::vector<double>;
  // _jStr_fwmEi_D5 = new std::vector<double>;
  // _jStr_fwmEi_D6 = new std::vector<double>;
  // _jStr_fwmEi_D7 = new std::vector<double>;
  // _jStr_fwmEi_D8 = new std::vector<double>;
  // _jStr_fwmEi_H0 = new std::vector<double>;
  // _jStr_fwmEi_H1 = new std::vector<double>;
  // _jStr_fwmEi_H2 = new std::vector<double>;
  // _jStr_fwmEi_H3 = new std::vector<double>;
  // _jStr_fwmEi_H4 = new std::vector<double>;
  // _jStr_fwmEi_H5 = new std::vector<double>;
  // _jStr_fwmEi_H6 = new std::vector<double>;
  // _jStr_fwmEi_H7 = new std::vector<double>;
  // _jStr_fwmEi_H8 = new std::vector<double>;
  // _jStr_fwmEi_Q0 = new std::vector<double>;
  // _jStr_fwmEi_Q1 = new std::vector<double>;
  // _jStr_fwmEi_Q2 = new std::vector<double>;
  // _jStr_fwmEi_Q3 = new std::vector<double>;
  // _jStr_fwmEi_Q4 = new std::vector<double>;
  // _jStr_fwmEi_Q5 = new std::vector<double>;
  // _jStr_fwmEi_Q6 = new std::vector<double>;
  // _jStr_fwmEi_Q7 = new std::vector<double>;
  // _jStr_fwmEi_Q8 = new std::vector<double>;
  // _jStr_fwmEi_Pi1 = new std::vector<double>;
  // _jStr_fwmEi_Pi2 = new std::vector<double>;
  // _jStr_fwmEi_Pi3 = new std::vector<double>;
  // _jStr_fwmEi_Pi4 = new std::vector<double>;
  // _jStr_fwmEi_B10 = new std::vector<double>;
  // _jStr_fwmEi_B20 = new std::vector<double>;
  // _jStr_fwmEi_B30 = new std::vector<double>;
  // _jStr_fwmEi_B40 = new std::vector<double>;
  // _jStr_fwmEi_B50 = new std::vector<double>;
  // _jStr_fwmEi_B60 = new std::vector<double>;
  // _jStr_fwmEi_B70 = new std::vector<double>;
  // _jStr_fwmEi_B80 = new std::vector<double>;
  // _jStr_fwmEi_C10 = new std::vector<double>;
  // _jStr_fwmEi_C20 = new std::vector<double>;
  // _jStr_fwmEi_C30 = new std::vector<double>;
  // _jStr_fwmEi_C40 = new std::vector<double>;
  // _jStr_fwmEi_C50 = new std::vector<double>;
  // _jStr_fwmEi_C60 = new std::vector<double>;
  // _jStr_fwmEi_C70 = new std::vector<double>;
  // _jStr_fwmEi_C80 = new std::vector<double>;
  // _jStr_fwmEi_K10 = new std::vector<double>;
  // _jStr_fwmEi_K20 = new std::vector<double>;
  // _jStr_fwmEi_K30 = new std::vector<double>;
  // _jStr_fwmEi_K40 = new std::vector<double>;
  // _jStr_fwmEi_K50 = new std::vector<double>;
  // _jStr_fwmEi_K60 = new std::vector<double>;
  // _jStr_fwmEi_K70 = new std::vector<double>;
  // _jStr_fwmEi_K80 = new std::vector<double>;
  // _jStr_fwmEi_D10 = new std::vector<double>;
  // _jStr_fwmEi_D20 = new std::vector<double>;
  // _jStr_fwmEi_D30 = new std::vector<double>;
  // _jStr_fwmEi_D40 = new std::vector<double>;
  // _jStr_fwmEi_D50 = new std::vector<double>;
  // _jStr_fwmEi_D60 = new std::vector<double>;
  // _jStr_fwmEi_D70 = new std::vector<double>;
  // _jStr_fwmEi_D80 = new std::vector<double>;
  // _jStr_fwmEi_H10 = new std::vector<double>;
  // _jStr_fwmEi_H20 = new std::vector<double>;
  // _jStr_fwmEi_H30 = new std::vector<double>;
  // _jStr_fwmEi_H40 = new std::vector<double>;
  // _jStr_fwmEi_H50 = new std::vector<double>;
  // _jStr_fwmEi_H60 = new std::vector<double>;
  // _jStr_fwmEi_H70 = new std::vector<double>;
  // _jStr_fwmEi_H80 = new std::vector<double>;
  // _jStr_fwmEi_Q10 = new std::vector<double>;
  // _jStr_fwmEi_Q20 = new std::vector<double>;
  // _jStr_fwmEi_Q30 = new std::vector<double>;
  // _jStr_fwmEi_Q40 = new std::vector<double>;
  // _jStr_fwmEi_Q50 = new std::vector<double>;
  // _jStr_fwmEi_Q60 = new std::vector<double>;
  // _jStr_fwmEi_Q70 = new std::vector<double>;
  // _jStr_fwmEi_Q80 = new std::vector<double>;

  _jStr_fwmJj_xiPlus = new std::vector<double>;
  _jStr_fwmJj_xiMinus = new std::vector<double>;

  _jStr_fwmJj_xPlus = new std::vector<double>;
  _jStr_fwmJj_xMinus = new std::vector<double>;

  // Fox-Wolfram, mini event, subjets as input:
  _jStr_fwmJj_Psi1 = new std::vector<double>;
  _jStr_fwmJj_B0 = new std::vector<double>;
  _jStr_fwmJj_B1 = new std::vector<double>;
  _jStr_fwmJj_B2 = new std::vector<double>;
  _jStr_fwmJj_B3 = new std::vector<double>;
  _jStr_fwmJj_B4 = new std::vector<double>;
  _jStr_fwmJj_B5 = new std::vector<double>;
  _jStr_fwmJj_B6 = new std::vector<double>;
  _jStr_fwmJj_B7 = new std::vector<double>;
  _jStr_fwmJj_B8 = new std::vector<double>;
  _jStr_fwmJj_C0 = new std::vector<double>;
  _jStr_fwmJj_C1 = new std::vector<double>;
  _jStr_fwmJj_C2 = new std::vector<double>;
  _jStr_fwmJj_C3 = new std::vector<double>;
  _jStr_fwmJj_C4 = new std::vector<double>;
  _jStr_fwmJj_C5 = new std::vector<double>;
  _jStr_fwmJj_C6 = new std::vector<double>;
  _jStr_fwmJj_C7 = new std::vector<double>;
  _jStr_fwmJj_C8 = new std::vector<double>;
  _jStr_fwmJj_K0 = new std::vector<double>;
  _jStr_fwmJj_K1 = new std::vector<double>;
  _jStr_fwmJj_K2 = new std::vector<double>;
  _jStr_fwmJj_K3 = new std::vector<double>;
  _jStr_fwmJj_K4 = new std::vector<double>;
  _jStr_fwmJj_K5 = new std::vector<double>;
  _jStr_fwmJj_K6 = new std::vector<double>;
  _jStr_fwmJj_K7 = new std::vector<double>;
  _jStr_fwmJj_K8 = new std::vector<double>;
  _jStr_fwmJj_D0 = new std::vector<double>;
  _jStr_fwmJj_D1 = new std::vector<double>;
  _jStr_fwmJj_D2 = new std::vector<double>;
  _jStr_fwmJj_D3 = new std::vector<double>;
  _jStr_fwmJj_D4 = new std::vector<double>;
  _jStr_fwmJj_D5 = new std::vector<double>;
  _jStr_fwmJj_D6 = new std::vector<double>;
  _jStr_fwmJj_D7 = new std::vector<double>;
  _jStr_fwmJj_D8 = new std::vector<double>;
  _jStr_fwmJj_H0 = new std::vector<double>;
  _jStr_fwmJj_H1 = new std::vector<double>;
  _jStr_fwmJj_H2 = new std::vector<double>;
  _jStr_fwmJj_H3 = new std::vector<double>;
  _jStr_fwmJj_H4 = new std::vector<double>;
  _jStr_fwmJj_H5 = new std::vector<double>;
  _jStr_fwmJj_H6 = new std::vector<double>;
  _jStr_fwmJj_H7 = new std::vector<double>;
  _jStr_fwmJj_H8 = new std::vector<double>;
  _jStr_fwmJj_Q0 = new std::vector<double>;
  _jStr_fwmJj_Q1 = new std::vector<double>;
  _jStr_fwmJj_Q2 = new std::vector<double>;
  _jStr_fwmJj_Q3 = new std::vector<double>;
  _jStr_fwmJj_Q4 = new std::vector<double>;
  _jStr_fwmJj_Q5 = new std::vector<double>;
  _jStr_fwmJj_Q6 = new std::vector<double>;
  _jStr_fwmJj_Q7 = new std::vector<double>;
  _jStr_fwmJj_Q8 = new std::vector<double>;
  _jStr_fwmJj_Pi1 = new std::vector<double>;
  _jStr_fwmJj_Pi2 = new std::vector<double>;
  _jStr_fwmJj_Pi3 = new std::vector<double>;
  _jStr_fwmJj_Pi4 = new std::vector<double>;
  _jStr_fwmJj_B10 = new std::vector<double>;
  _jStr_fwmJj_B20 = new std::vector<double>;
  _jStr_fwmJj_B30 = new std::vector<double>;
  _jStr_fwmJj_B40 = new std::vector<double>;
  _jStr_fwmJj_B50 = new std::vector<double>;
  _jStr_fwmJj_B60 = new std::vector<double>;
  _jStr_fwmJj_B70 = new std::vector<double>;
  _jStr_fwmJj_B80 = new std::vector<double>;
  _jStr_fwmJj_C10 = new std::vector<double>;
  _jStr_fwmJj_C20 = new std::vector<double>;
  _jStr_fwmJj_C30 = new std::vector<double>;
  _jStr_fwmJj_C40 = new std::vector<double>;
  _jStr_fwmJj_C50 = new std::vector<double>;
  _jStr_fwmJj_C60 = new std::vector<double>;
  _jStr_fwmJj_C70 = new std::vector<double>;
  _jStr_fwmJj_C80 = new std::vector<double>;
  _jStr_fwmJj_K10 = new std::vector<double>;
  _jStr_fwmJj_K20 = new std::vector<double>;
  _jStr_fwmJj_K30 = new std::vector<double>;
  _jStr_fwmJj_K40 = new std::vector<double>;
  _jStr_fwmJj_K50 = new std::vector<double>;
  _jStr_fwmJj_K60 = new std::vector<double>;
  _jStr_fwmJj_K70 = new std::vector<double>;
  _jStr_fwmJj_K80 = new std::vector<double>;
  _jStr_fwmJj_D10 = new std::vector<double>;
  _jStr_fwmJj_D20 = new std::vector<double>;
  _jStr_fwmJj_D30 = new std::vector<double>;
  _jStr_fwmJj_D40 = new std::vector<double>;
  _jStr_fwmJj_D50 = new std::vector<double>;
  _jStr_fwmJj_D60 = new std::vector<double>;
  _jStr_fwmJj_D70 = new std::vector<double>;
  _jStr_fwmJj_D80 = new std::vector<double>;
  _jStr_fwmJj_H10 = new std::vector<double>;
  _jStr_fwmJj_H20 = new std::vector<double>;
  _jStr_fwmJj_H30 = new std::vector<double>;
  _jStr_fwmJj_H40 = new std::vector<double>;
  _jStr_fwmJj_H50 = new std::vector<double>;
  _jStr_fwmJj_H60 = new std::vector<double>;
  _jStr_fwmJj_H70 = new std::vector<double>;
  _jStr_fwmJj_H80 = new std::vector<double>;
  _jStr_fwmJj_Q10 = new std::vector<double>;
  _jStr_fwmJj_Q20 = new std::vector<double>;
  _jStr_fwmJj_Q30 = new std::vector<double>;
  _jStr_fwmJj_Q40 = new std::vector<double>;
  _jStr_fwmJj_Q50 = new std::vector<double>;
  _jStr_fwmJj_Q60 = new std::vector<double>;
  _jStr_fwmJj_Q70 = new std::vector<double>;
  _jStr_fwmJj_Q80 = new std::vector<double>;

  _jStr_fwmJi_xiPlus = new std::vector<double>;
  _jStr_fwmJi_xiMinus = new std::vector<double>;

  _jStr_fwmJi_xPlus = new std::vector<double>;
  _jStr_fwmJi_xMinus = new std::vector<double>;

  // Fox-Wolfram, mini event, jet constituents as input:
  _jStr_fwmJi_Psi1 = new std::vector<double>;
  _jStr_fwmJi_B0 = new std::vector<double>;
  _jStr_fwmJi_B1 = new std::vector<double>;
  _jStr_fwmJi_B2 = new std::vector<double>;
  _jStr_fwmJi_B3 = new std::vector<double>;
  _jStr_fwmJi_B4 = new std::vector<double>;
  _jStr_fwmJi_B5 = new std::vector<double>;
  _jStr_fwmJi_B6 = new std::vector<double>;
  _jStr_fwmJi_B7 = new std::vector<double>;
  _jStr_fwmJi_B8 = new std::vector<double>;
  _jStr_fwmJi_C0 = new std::vector<double>;
  _jStr_fwmJi_C1 = new std::vector<double>;
  _jStr_fwmJi_C2 = new std::vector<double>;
  _jStr_fwmJi_C3 = new std::vector<double>;
  _jStr_fwmJi_C4 = new std::vector<double>;
  _jStr_fwmJi_C5 = new std::vector<double>;
  _jStr_fwmJi_C6 = new std::vector<double>;
  _jStr_fwmJi_C7 = new std::vector<double>;
  _jStr_fwmJi_C8 = new std::vector<double>;
  _jStr_fwmJi_K0 = new std::vector<double>;
  _jStr_fwmJi_K1 = new std::vector<double>;
  _jStr_fwmJi_K2 = new std::vector<double>;
  _jStr_fwmJi_K3 = new std::vector<double>;
  _jStr_fwmJi_K4 = new std::vector<double>;
  _jStr_fwmJi_K5 = new std::vector<double>;
  _jStr_fwmJi_K6 = new std::vector<double>;
  _jStr_fwmJi_K7 = new std::vector<double>;
  _jStr_fwmJi_K8 = new std::vector<double>;
  _jStr_fwmJi_D0 = new std::vector<double>;
  _jStr_fwmJi_D1 = new std::vector<double>;
  _jStr_fwmJi_D2 = new std::vector<double>;
  _jStr_fwmJi_D3 = new std::vector<double>;
  _jStr_fwmJi_D4 = new std::vector<double>;
  _jStr_fwmJi_D5 = new std::vector<double>;
  _jStr_fwmJi_D6 = new std::vector<double>;
  _jStr_fwmJi_D7 = new std::vector<double>;
  _jStr_fwmJi_D8 = new std::vector<double>;
  _jStr_fwmJi_H0 = new std::vector<double>;
  _jStr_fwmJi_H1 = new std::vector<double>;
  _jStr_fwmJi_H2 = new std::vector<double>;
  _jStr_fwmJi_H3 = new std::vector<double>;
  _jStr_fwmJi_H4 = new std::vector<double>;
  _jStr_fwmJi_H5 = new std::vector<double>;
  _jStr_fwmJi_H6 = new std::vector<double>;
  _jStr_fwmJi_H7 = new std::vector<double>;
  _jStr_fwmJi_H8 = new std::vector<double>;
  _jStr_fwmJi_Q0 = new std::vector<double>;
  _jStr_fwmJi_Q1 = new std::vector<double>;
  _jStr_fwmJi_Q2 = new std::vector<double>;
  _jStr_fwmJi_Q3 = new std::vector<double>;
  _jStr_fwmJi_Q4 = new std::vector<double>;
  _jStr_fwmJi_Q5 = new std::vector<double>;
  _jStr_fwmJi_Q6 = new std::vector<double>;
  _jStr_fwmJi_Q7 = new std::vector<double>;
  _jStr_fwmJi_Q8 = new std::vector<double>;
  _jStr_fwmJi_Pi1 = new std::vector<double>;
  _jStr_fwmJi_Pi2 = new std::vector<double>;
  _jStr_fwmJi_Pi3 = new std::vector<double>;
  _jStr_fwmJi_Pi4 = new std::vector<double>;
  _jStr_fwmJi_B10 = new std::vector<double>;
  _jStr_fwmJi_B20 = new std::vector<double>;
  _jStr_fwmJi_B30 = new std::vector<double>;
  _jStr_fwmJi_B40 = new std::vector<double>;
  _jStr_fwmJi_B50 = new std::vector<double>;
  _jStr_fwmJi_B60 = new std::vector<double>;
  _jStr_fwmJi_B70 = new std::vector<double>;
  _jStr_fwmJi_B80 = new std::vector<double>;
  _jStr_fwmJi_C10 = new std::vector<double>;
  _jStr_fwmJi_C20 = new std::vector<double>;
  _jStr_fwmJi_C30 = new std::vector<double>;
  _jStr_fwmJi_C40 = new std::vector<double>;
  _jStr_fwmJi_C50 = new std::vector<double>;
  _jStr_fwmJi_C60 = new std::vector<double>;
  _jStr_fwmJi_C70 = new std::vector<double>;
  _jStr_fwmJi_C80 = new std::vector<double>;
  _jStr_fwmJi_K10 = new std::vector<double>;
  _jStr_fwmJi_K20 = new std::vector<double>;
  _jStr_fwmJi_K30 = new std::vector<double>;
  _jStr_fwmJi_K40 = new std::vector<double>;
  _jStr_fwmJi_K50 = new std::vector<double>;
  _jStr_fwmJi_K60 = new std::vector<double>;
  _jStr_fwmJi_K70 = new std::vector<double>;
  _jStr_fwmJi_K80 = new std::vector<double>;
  _jStr_fwmJi_D10 = new std::vector<double>;
  _jStr_fwmJi_D20 = new std::vector<double>;
  _jStr_fwmJi_D30 = new std::vector<double>;
  _jStr_fwmJi_D40 = new std::vector<double>;
  _jStr_fwmJi_D50 = new std::vector<double>;
  _jStr_fwmJi_D60 = new std::vector<double>;
  _jStr_fwmJi_D70 = new std::vector<double>;
  _jStr_fwmJi_D80 = new std::vector<double>;
  _jStr_fwmJi_H10 = new std::vector<double>;
  _jStr_fwmJi_H20 = new std::vector<double>;
  _jStr_fwmJi_H30 = new std::vector<double>;
  _jStr_fwmJi_H40 = new std::vector<double>;
  _jStr_fwmJi_H50 = new std::vector<double>;
  _jStr_fwmJi_H60 = new std::vector<double>;
  _jStr_fwmJi_H70 = new std::vector<double>;
  _jStr_fwmJi_H80 = new std::vector<double>;
  _jStr_fwmJi_Q10 = new std::vector<double>;
  _jStr_fwmJi_Q20 = new std::vector<double>;
  _jStr_fwmJi_Q30 = new std::vector<double>;
  _jStr_fwmJi_Q40 = new std::vector<double>;
  _jStr_fwmJi_Q50 = new std::vector<double>;
  _jStr_fwmJi_Q60 = new std::vector<double>;
  _jStr_fwmJi_Q70 = new std::vector<double>;
  _jStr_fwmJi_Q80 = new std::vector<double>;

  // Flags:
  _jStr_authorE = new std::vector<std::string>;
  _jStr_authorJ = new std::vector<std::string>;
  _jStr_radialParamE = new std::vector<int>;
  _jStr_radialParamJ = new std::vector<int>;
  _jStr_algE = new std::vector<int>;
  _jStr_algJ = new std::vector<int>;
  _jStr_inputE = new std::vector<int>;
  _jStr_inputJ = new std::vector<int>;
  _jStr_bjetE = new std::vector<int>;
  _jStr_bjetJ = new std::vector<int>;
  _jStr_xjetE = new std::vector<int>;
  _jStr_xjetJ = new std::vector<int>;
  _jStr_myJetsE = new std::vector<int>;
  _jStr_myJetsJ = new std::vector<int>;
  _jStr_indexJ = new std::vector<int>;

  // Misc jet properties: 
  _jStr_jetMiscJ_numConstituents = new std::vector<int>;
  _jStr_jetMiscJ_jetM = new std::vector<double>;
  _jStr_jetMiscJ_jetMt = new std::vector<double>;
  _jStr_jetMiscJ_jetE = new std::vector<double>;
  _jStr_jetMiscJ_jetP = new std::vector<double>;
  _jStr_jetMiscJ_jetEt = new std::vector<double>;
  _jStr_jetMiscJ_jetPt = new std::vector<double>;
  _jStr_jetMiscJ_jetPhi = new std::vector<double>;
  _jStr_jetMiscJ_jetEta = new std::vector<double>;
  _jStr_jetMiscJ_jetRapidity = new std::vector<double>;
  _jStr_jetMiscJ_xJ = new std::vector<double>;
  _jStr_jetMiscJ_gamma = new std::vector<double>;
  _jStr_jetMiscJ_R = new std::vector<double>;
  _jStr_jetMiscJ_Cbar = new std::vector<double>;
  _jStr_jetMiscJ_numSubjets = new std::vector<int>;

  // Jet pull vars: 
  _jStr_pullJ_det = new std::vector<double>;
  _jStr_pullJ_ratio = new std::vector<double>;
  _jStr_pullJ_pullPf = new std::vector<double>;
  _jStr_pullJ_angularEccentricity = new std::vector<double>;
  _jStr_pullJ_orientation = new std::vector<double>;
  _jStr_pullJ_girth = new std::vector<double>;
  _jStr_pullJ_Cbar = new std::vector<double>;
  _jStr_pullJ_g = new std::vector<double>;
  _jStr_pullJ_e = new std::vector<double>;
  _jStr_pullJ_B = new std::vector<double>;
  _jStr_pullJ_logB = new std::vector<double>;
  _jStr_pullJ_pullTheta = new std::vector<double>;
  _jStr_pullJ_pullMag = new std::vector<double>;

  // Unburied Higgs vars: 
  _jStr_ubhJ_z = new std::vector<double>;
  _jStr_ubhJ_z2 = new std::vector<double>;
  _jStr_ubhJ_a1 = new std::vector<double>;
  _jStr_ubhJ_a2 = new std::vector<double>;
  _jStr_ubhJ_a3 = new std::vector<double>;
  _jStr_ubhJ_meanpt = new std::vector<double>;
  _jStr_ubhJ_meanet = new std::vector<double>;
  _jStr_ubhJ_mbar = new std::vector<double>;
  _jStr_ubhJ_massDemocracy = new std::vector<double>;
  _jStr_ubhJ_fE1 = new std::vector<double>;
  _jStr_ubhJ_fE2 = new std::vector<double>;
  _jStr_ubhJ_fE3 = new std::vector<double>;
  _jStr_ubhJ_fET1 = new std::vector<double>;
  _jStr_ubhJ_fET2 = new std::vector<double>;
  _jStr_ubhJ_fET3 = new std::vector<double>;
  _jStr_ubhJ_Alpha = new std::vector<double>;
  _jStr_ubhJ_AlphaT = new std::vector<double>;
  _jStr_ubhJ_betaflow = new std::vector<double>;
  _jStr_ubhJ_betaflow_1GeV = new std::vector<double>;
  _jStr_ubhJ_betaflow_5GeV = new std::vector<double>;
  _jStr_ubhJ_y23 = new std::vector<double>;
  _jStr_ubhJ_y23_1GeV = new std::vector<double>;
  _jStr_ubhJ_y23_5GeV = new std::vector<double>;
  _jStr_ubhJ_lny23 = new std::vector<double>;
  _jStr_ubhJ_lny23_1GeV = new std::vector<double>;
  _jStr_ubhJ_lny23_5GeV = new std::vector<double>;
  _jStr_ubhJ_subjetAsymmetry = new std::vector<double>;
  _jStr_dipJ_dipolarity = new std::vector<double>;

  // N-subjettiness vars: 
  _jStr_nsubjnessJ_tau1 = new std::vector<double>;
  _jStr_nsubjnessJ_tau2 = new std::vector<double>;
  _jStr_nsubjnessJ_tau3 = new std::vector<double>;
  _jStr_nsubjnessJ_tau2tau1 = new std::vector<double>;
  _jStr_nsubjnessJ_tau3tau2 = new std::vector<double>;

  // psi jet shapes: 
  _jStr_psiJ2_psi = new std::vector<double>;
  _jStr_psiJ2_rho = new std::vector<double>;
  _jStr_psiRatiosJ_psi1 = new std::vector<double>;
  _jStr_psiRatiosJ_psi2 = new std::vector<double>;
  _jStr_psiRatiosJ_psi3 = new std::vector<double>;
  _jStr_psiRatiosJ_psi7 = new std::vector<double>;
  _jStr_psiRatiosJ_psi717 = new std::vector<double>;
  _jStr_psiRatiosJ_psi127 = new std::vector<double>;
  _jStr_psiRatiosJ_psi37 = new std::vector<double>;

  // BCF jet shapes: 
  _jStr_bcfJ_bcfVersion_v0a1 = new std::vector<int>;
  _jStr_bcfJ_a_v0a1 = new std::vector<int>;
  _jStr_bcfJ_bcfT_v0a1 = new std::vector<double>;
  _jStr_bcfJ_bcf_v0a1 = new std::vector<double>;
  _jStr_bcfJ_bcfAsymY_v0a1 = new std::vector<double>;
  _jStr_bcfJ_bcfAsymPhi_v0a1 = new std::vector<double>;
  _jStr_bcfJ_bcfAsymYPhi_v0a1 = new std::vector<double>;
  _jStr_bcfJ_bcfAsymYPhi2_v0a1 = new std::vector<double>;

  _jStr_bcfJ_bcfVersion_v0a2 = new std::vector<int>;
  _jStr_bcfJ_a_v0a2 = new std::vector<int>;
  _jStr_bcfJ_bcfT_v0a2 = new std::vector<double>;
  _jStr_bcfJ_bcf_v0a2 = new std::vector<double>;
  _jStr_bcfJ_bcfAsymY_v0a2 = new std::vector<double>;
  _jStr_bcfJ_bcfAsymPhi_v0a2 = new std::vector<double>;
  _jStr_bcfJ_bcfAsymYPhi_v0a2 = new std::vector<double>;
  _jStr_bcfJ_bcfAsymYPhi2_v0a2 = new std::vector<double>;

  // planar flow jet shapes:
  _jStr_pfJ_pf = new std::vector<double>;
  _jStr_pfJ_detST = new std::vector<double>;
  _jStr_pfJ_lambdaST = new std::vector<double>;

  // splittings:
  _jStr_zminJ = new std::vector<double>;
  _jStr_zmaxJ = new std::vector<double>;

  // angularities
  _jStr_tauJ09_a = new std::vector<double>;
  _jStr_tauJ09_tau = new std::vector<double>;
  _jStr_tauJ20_a = new std::vector<double>;
  _jStr_tauJ20_tau = new std::vector<double>;
  _jStr_tauJ40_a = new std::vector<double>;
  _jStr_tauJ40_tau = new std::vector<double>;

  _jStr_radEi_bcfJ1_v0a1 = new std::vector<double>;
  _jStr_radEi_bcfJ2_v0a1 = new std::vector<double>;
  _jStr_radEi_bcfJ3_v0a1 = new std::vector<double>;
  _jStr_radEi_bcfTJ1_v0a1 = new std::vector<double>;
  _jStr_radEi_bcfTJ2_v0a1 = new std::vector<double>;
  _jStr_radEi_bcfTJ3_v0a1 = new std::vector<double>;
  _jStr_radEi_bcfAsymYJ1_v0a1 = new std::vector<double>;
  _jStr_radEi_bcfAsymYJ2_v0a1 = new std::vector<double>;
  _jStr_radEi_bcfAsymYJ3_v0a1 = new std::vector<double>;
  _jStr_radEi_bcfAsymPhiJ1_v0a1 = new std::vector<double>;
  _jStr_radEi_bcfAsymPhiJ2_v0a1 = new std::vector<double>;
  _jStr_radEi_bcfAsymPhiJ3_v0a1 = new std::vector<double>;
  _jStr_radEi_bcfAsymYPhiJ1_v0a1 = new std::vector<double>;
  _jStr_radEi_bcfAsymYPhiJ2_v0a1 = new std::vector<double>;
  _jStr_radEi_bcfAsymYPhiJ3_v0a1 = new std::vector<double>;
  _jStr_radEi_bcfAsymYPhi2J1_v0a1 = new std::vector<double>;
  _jStr_radEi_bcfAsymYPhi2J2_v0a1 = new std::vector<double>;
  _jStr_radEi_bcfAsymYPhi2J3_v0a1 = new std::vector<double>;

  _jStr_radEi_bcfJ1_v0a2 = new std::vector<double>;
  _jStr_radEi_bcfJ2_v0a2 = new std::vector<double>;
  _jStr_radEi_bcfJ3_v0a2 = new std::vector<double>;
  _jStr_radEi_bcfTJ1_v0a2 = new std::vector<double>;
  _jStr_radEi_bcfTJ2_v0a2 = new std::vector<double>;
  _jStr_radEi_bcfTJ3_v0a2 = new std::vector<double>;
  _jStr_radEi_bcfAsymYJ1_v0a2 = new std::vector<double>;
  _jStr_radEi_bcfAsymYJ2_v0a2 = new std::vector<double>;
  _jStr_radEi_bcfAsymYJ3_v0a2 = new std::vector<double>;
  _jStr_radEi_bcfAsymPhiJ1_v0a2 = new std::vector<double>;
  _jStr_radEi_bcfAsymPhiJ2_v0a2 = new std::vector<double>;
  _jStr_radEi_bcfAsymPhiJ3_v0a2 = new std::vector<double>;
  _jStr_radEi_bcfAsymYPhiJ1_v0a2 = new std::vector<double>;
  _jStr_radEi_bcfAsymYPhiJ2_v0a2 = new std::vector<double>;
  _jStr_radEi_bcfAsymYPhiJ3_v0a2 = new std::vector<double>;
  _jStr_radEi_bcfAsymYPhi2J1_v0a2 = new std::vector<double>;
  _jStr_radEi_bcfAsymYPhi2J2_v0a2 = new std::vector<double>;
  _jStr_radEi_bcfAsymYPhi2J3_v0a2 = new std::vector<double>;

  _jStr_radEi_BJ1 = new std::vector<double>;
  _jStr_radEi_BJ2 = new std::vector<double>;
  _jStr_radEi_BJ3 = new std::vector<double>;
  _jStr_radEi_girthJ1 = new std::vector<double>;
  _jStr_radEi_girthJ2 = new std::vector<double>;
  _jStr_radEi_girthJ3 = new std::vector<double>;
  _jStr_radEi_girth32 = new std::vector<double>;
  _jStr_radEi_girth21 = new std::vector<double>;
  _jStr_radEi_girthAsymJ1J2 = new std::vector<double>;
  _jStr_radEi_alpha1 = new std::vector<double>;
  _jStr_radEi_alpha2 = new std::vector<double>;
  _jStr_radEi_alpha = new std::vector<double>;
  _jStr_radEi_beta1 = new std::vector<double>;
  _jStr_radEi_beta2 = new std::vector<double>;
  _jStr_radEi_beta = new std::vector<double>;
  _jStr_radEi_thetaJ1J2 = new std::vector<double>;
  _jStr_radEi_dipolarityInfLine = new std::vector<double>;
  _jStr_radEi_dipolarityLineSeg = new std::vector<double>;
  _jStr_radEi_BCF1 = new std::vector<double>;
  _jStr_radEi_BCF2 = new std::vector<double>;
  _jStr_radEi_BCF3 = new std::vector<double>;
  _jStr_radEi_dipolarity = new std::vector<double>;

  // whole event, jets as input:
  _jStr_sphEJ_detSphericity = new std::vector<double>;
  _jStr_sphEJ_detSpherocity = new std::vector<double>;
  _jStr_sphEJ_sphericityLambda1 = new std::vector<double>;
  _jStr_sphEJ_sphericityLambda2 = new std::vector<double>;
  _jStr_sphEJ_sphericityLambda3 = new std::vector<double>;
  _jStr_sphEJ_spherocityLambda1 = new std::vector<double>;
  _jStr_sphEJ_spherocityLambda2 = new std::vector<double>;
  _jStr_sphEJ_spherocityLambda3 = new std::vector<double>;
  _jStr_sphEJ_circularity = new std::vector<double>;
  _jStr_sphEJ_sphericity = new std::vector<double>;
  _jStr_sphEJ_spherocity = new std::vector<double>;
  _jStr_sphEJ_aplanarity = new std::vector<double>;
  _jStr_sphEJ_aplanority = new std::vector<double>;
  _jStr_sphEJ_Y = new std::vector<double>;
  _jStr_sphEJ_planarity = new std::vector<double>;
  _jStr_sphEJ_planority = new std::vector<double>;
  _jStr_sphEJ_Dshape = new std::vector<double>;
  _jStr_sphEJ_Cshape = new std::vector<double>;
  _jStr_sphEJ_H2 = new std::vector<double>;
  _jStr_sphEJ_Fshape = new std::vector<double>;
  _jStr_sphEJ_beamThrust = new std::vector<double>;
  _jStr_sphEJ_G = new std::vector<double>;
  _jStr_sphEJ_ST2D = new std::vector<double>;
  _jStr_sphEJ_detMlin = new std::vector<double>;
  _jStr_sphEJ_pft = new std::vector<double>;

  // whole event, jet constituents as input:
  _jStr_sphEi_detSphericity = new std::vector<double>;
  _jStr_sphEi_detSpherocity = new std::vector<double>;
  _jStr_sphEi_sphericityLambda1 = new std::vector<double>;
  _jStr_sphEi_sphericityLambda2 = new std::vector<double>;
  _jStr_sphEi_sphericityLambda3 = new std::vector<double>;
  _jStr_sphEi_spherocityLambda1 = new std::vector<double>;
  _jStr_sphEi_spherocityLambda2 = new std::vector<double>;
  _jStr_sphEi_spherocityLambda3 = new std::vector<double>;
  _jStr_sphEi_circularity = new std::vector<double>;
  _jStr_sphEi_sphericity = new std::vector<double>;
  _jStr_sphEi_spherocity = new std::vector<double>;
  _jStr_sphEi_aplanarity = new std::vector<double>;
  _jStr_sphEi_aplanority = new std::vector<double>;
  _jStr_sphEi_Y = new std::vector<double>;
  _jStr_sphEi_planarity = new std::vector<double>;
  _jStr_sphEi_planority = new std::vector<double>;
  _jStr_sphEi_Dshape = new std::vector<double>;
  _jStr_sphEi_Cshape = new std::vector<double>;
  _jStr_sphEi_H2 = new std::vector<double>;
  _jStr_sphEi_Fshape = new std::vector<double>;
  _jStr_sphEi_beamThrust = new std::vector<double>;
  _jStr_sphEi_G = new std::vector<double>;
  _jStr_sphEi_ST2D = new std::vector<double>;
  _jStr_sphEi_detMlin = new std::vector<double>;
  _jStr_sphEi_pft = new std::vector<double>;

  // mini event, subjets as input:
  _jStr_sphJj_detSphericity = new std::vector<double>;
  _jStr_sphJj_detSpherocity = new std::vector<double>;
  _jStr_sphJj_sphericityLambda1 = new std::vector<double>;
  _jStr_sphJj_sphericityLambda2 = new std::vector<double>;
  _jStr_sphJj_sphericityLambda3 = new std::vector<double>;
  _jStr_sphJj_spherocityLambda1 = new std::vector<double>;
  _jStr_sphJj_spherocityLambda2 = new std::vector<double>;
  _jStr_sphJj_spherocityLambda3 = new std::vector<double>;
  _jStr_sphJj_circularity = new std::vector<double>;
  _jStr_sphJj_sphericity = new std::vector<double>;
  _jStr_sphJj_spherocity = new std::vector<double>;
  _jStr_sphJj_aplanarity = new std::vector<double>;
  _jStr_sphJj_aplanority = new std::vector<double>;
  _jStr_sphJj_Y = new std::vector<double>;
  _jStr_sphJj_planarity = new std::vector<double>;
  _jStr_sphJj_planority = new std::vector<double>;
  _jStr_sphJj_Dshape = new std::vector<double>;
  _jStr_sphJj_Cshape = new std::vector<double>;
  _jStr_sphJj_H2 = new std::vector<double>;
  _jStr_sphJj_Fshape = new std::vector<double>;
  _jStr_sphJj_beamThrust = new std::vector<double>;
  _jStr_sphJj_G = new std::vector<double>;
  _jStr_sphJj_ST2D = new std::vector<double>;
  _jStr_sphJj_detMlin = new std::vector<double>;
  _jStr_sphJj_pft = new std::vector<double>;

  // mini event, jet constituents as input:
  _jStr_sphJi_detSphericity = new std::vector<double>;
  _jStr_sphJi_detSpherocity = new std::vector<double>;
  _jStr_sphJi_sphericityLambda1 = new std::vector<double>;
  _jStr_sphJi_sphericityLambda2 = new std::vector<double>;
  _jStr_sphJi_sphericityLambda3 = new std::vector<double>;
  _jStr_sphJi_spherocityLambda1 = new std::vector<double>;
  _jStr_sphJi_spherocityLambda2 = new std::vector<double>;
  _jStr_sphJi_spherocityLambda3 = new std::vector<double>;
  _jStr_sphJi_circularity = new std::vector<double>;
  _jStr_sphJi_sphericity = new std::vector<double>;
  _jStr_sphJi_spherocity = new std::vector<double>;
  _jStr_sphJi_aplanarity = new std::vector<double>;
  _jStr_sphJi_aplanority = new std::vector<double>;
  _jStr_sphJi_Y = new std::vector<double>;
  _jStr_sphJi_planarity = new std::vector<double>;
  _jStr_sphJi_planority = new std::vector<double>;
  _jStr_sphJi_Dshape = new std::vector<double>;
  _jStr_sphJi_Cshape = new std::vector<double>;
  _jStr_sphJi_H2 = new std::vector<double>;
  _jStr_sphJi_Fshape = new std::vector<double>;
  _jStr_sphJi_beamThrust = new std::vector<double>;
  _jStr_sphJi_G = new std::vector<double>;
  _jStr_sphJi_ST2D = new std::vector<double>;
  _jStr_sphJi_detMlin = new std::vector<double>;
  _jStr_sphJi_pft = new std::vector<double>;

  _nConst = new std::vector<int>;
  _const_et = new std::vector<std::vector<double> >;
  _const_rapidity = new std::vector<std::vector<double> >;
  _const_phi = new std::vector<std::vector<double> >;
  _const_Ret = new std::vector<std::vector<double> >;
  _const_Drapidity = new std::vector<std::vector<double> >;
  _const_Dphi = new std::vector<std::vector<double> >;

  _YFlip12 = new std::vector<double>;
  _YFlip23 = new std::vector<double>;
  _EM_FRACTION_CLUSTER = new std::vector<double>;
  _ELLIPTICAREA = new std::vector<double>;
  _AMAREA = new std::vector<double>;
  _HULL_LENGTH = new std::vector<double>;
  _HULL_AREA = new std::vector<double>;
  _EM_FRACTION_MCTRUTH = new std::vector<double>;
  _LowEtConstituentsFrac = new std::vector<double>;
  _JetEccentricity = new std::vector<double>;
  _DRtoReco = new std::vector<double>;
  _PtNearest = new std::vector<double>;
  _WIDTH = new std::vector<double>;
  _rbb = new std::vector<double>;
  _rfilt = new std::vector<double>;

  _jvf = new std::vector<double>;
  _ntrk = new std::vector<int>;
  _trkpt = new std::vector<double>;

  _QGtag = new std::vector<long long int>;
  _QGtag2 = new std::vector<long long int>;
  _jetTrueFlavour = new std::vector<int>;
  w_cmb = new std::vector<double>;
  w_TrackCounting2D = new std::vector<double>;
  w_JetProb = new std::vector<double>;
  w_IP1D = new std::vector<double>;
  w_IP2D = new std::vector<double>;
  w_IP3D = new std::vector<double>;
  w_SV0 = new std::vector<double>;
  w_SV1 = new std::vector<double>;
  w_SV2 = new std::vector<double>;
  w_BaselineTagger = new std::vector<double>;
  w_JetFitterTag = new std::vector<double>;
  w_JetFitterCOMB = new std::vector<double>;
  w_JetFitterTagNN = new std::vector<double>;
  w_JetFitterCOMBNN = new std::vector<double>;
  w_SoftMuonTag = new std::vector<double>;
  w_SoftElectronTag = new std::vector<double>;

  /** now add branches and leaves to the AAN tree */
  // the TTree
  m_tree_AS = new TTree("tree_AS","TTree of AnalysisSkleton");
  sc = m_thistSvc->regTree("/AANT/tree_AS", m_tree_AS);

  // first add Event info stuff
  m_tree_AS->Branch("Run",  &m_runNumber,   "Run/I");    // run number
  m_tree_AS->Branch("Event",&m_eventNumber, "Event/I");  // event number
  m_tree_AS->Branch("Time", &m_eventTime,   "Time/I");   // time stamp
  m_tree_AS->Branch("LumiBlock", &m_lumiBlock,"LumiBlock/I"); // lum block num 
  m_tree_AS->Branch("BCID", &m_bCID,"BCID/I"); // bunch crossing ID
  m_tree_AS->Branch("Weight", &m_eventWeight, "Weight/D"); // weight

  m_tree_AS->Branch("Nparticles", &_nparticles, "Nparticles/I");
  m_tree_AS->Branch("NFPE", &_NFPE, "NFPE/I");

  m_tree_AS->Branch("jStr_numPVE", &_jStr_numPVE, "jStr_numPVE/I");
  m_tree_AS->Branch("jStr_PVxE", &_jStr_PVxE, "jStr_PVxE/D");
  m_tree_AS->Branch("jStr_PVyE", &_jStr_PVyE, "jStr_PVyE/D");
  m_tree_AS->Branch("jStr_PVzE", &_jStr_PVzE, "jStr_PVzE/D");
  m_tree_AS->Branch("jStr_PVrE", &_jStr_PVrE, "jStr_PVrE/D");
  m_tree_AS->Branch("jStr_PVchiSqrE", &_jStr_PVchiSqrE, "jStr_PVchiSqrE/D");
  m_tree_AS->Branch("jStr_PVnDoFE", &_jStr_PVnDoFE, "jStr_PVnDoFE/D");
  m_tree_AS->Branch("jStr_PVfitE", &_jStr_PVfitE, "jStr_PVfitE/D");
  m_tree_AS->Branch("jStr_PVnTrkE", &_jStr_PVnTrkE, "jStr_PVnTrkE/I");
  m_tree_AS->Branch("jStr_PVtypeE", &_jStr_PVtypeE, "jStr_PVtypeE/I");

  m_tree_AS->Branch("jStr_numPVJ", &_jStr_numPVJ);
  m_tree_AS->Branch("jStr_PVxJ", &_jStr_PVxJ);
  m_tree_AS->Branch("jStr_PVyJ", &_jStr_PVyJ);
  m_tree_AS->Branch("jStr_PVzJ", &_jStr_PVzJ);
  m_tree_AS->Branch("jStr_PVrJ", &_jStr_PVrJ);
  m_tree_AS->Branch("jStr_PVchiSqrJ", &_jStr_PVchiSqrJ);
  m_tree_AS->Branch("jStr_PVnDoFJ", &_jStr_PVnDoFJ);
  m_tree_AS->Branch("jStr_PVfitJ", &_jStr_PVfitJ);
  m_tree_AS->Branch("jStr_PVnTrkJ", &_jStr_PVnTrkJ);
  m_tree_AS->Branch("jStr_PVtypeJ", &_jStr_PVtypeJ);

  m_tree_AS->Branch("jStr_numJets", &_jStr_numJets, "jStr_numJets/i");

  m_tree_AS->Branch("jStr_evtEJ_etaJ1", &_jStr_evtEJ_etaJ1);
  m_tree_AS->Branch("jStr_evtEJ_etaJ2", &_jStr_evtEJ_etaJ2);
  m_tree_AS->Branch("jStr_evtEJ_etaJ3", &_jStr_evtEJ_etaJ3);
  m_tree_AS->Branch("jStr_evtEJ_etaJ4", &_jStr_evtEJ_etaJ4);
  m_tree_AS->Branch("jStr_evtEJ_pTJ1", &_jStr_evtEJ_pTJ1);
  m_tree_AS->Branch("jStr_evtEJ_pTJ2", &_jStr_evtEJ_pTJ2);
  m_tree_AS->Branch("jStr_evtEJ_pTJ3", &_jStr_evtEJ_pTJ3);
  m_tree_AS->Branch("jStr_evtEJ_pTJ4", &_jStr_evtEJ_pTJ4);

  m_tree_AS->Branch("jStr_evtEJ_fE1", &_jStr_evtEJ_fE1);
  m_tree_AS->Branch("jStr_evtEJ_fE2", &_jStr_evtEJ_fE2);
  m_tree_AS->Branch("jStr_evtEJ_fE3", &_jStr_evtEJ_fE3);
  m_tree_AS->Branch("jStr_evtEJ_fET1", &_jStr_evtEJ_fET1);
  m_tree_AS->Branch("jStr_evtEJ_fET2", &_jStr_evtEJ_fET2);
  m_tree_AS->Branch("jStr_evtEJ_fET3", &_jStr_evtEJ_fET3);
  m_tree_AS->Branch("jStr_evtEJ_mct", &_jStr_evtEJ_mct);
  m_tree_AS->Branch("jStr_evtEJ_q", &_jStr_evtEJ_q);
  m_tree_AS->Branch("jStr_evtEJ_mjj", &_jStr_evtEJ_mjj);
  m_tree_AS->Branch("jStr_evtEJ_mTjj", &_jStr_evtEJ_mTjj);
  m_tree_AS->Branch("jStr_evtEJ_mjjj", &_jStr_evtEJ_mjjj);
  m_tree_AS->Branch("jStr_evtEJ_mTjjj", &_jStr_evtEJ_mTjjj);
  m_tree_AS->Branch("jStr_evtEJ_mjjjj", &_jStr_evtEJ_mjjjj);
  m_tree_AS->Branch("jStr_evtEJ_mTjjjj", &_jStr_evtEJ_mTjjjj);
  m_tree_AS->Branch("jStr_evtEJ_zetaPlus", &_jStr_evtEJ_zetaPlus);
  m_tree_AS->Branch("jStr_evtEJ_zetaMinus", &_jStr_evtEJ_zetaMinus);
  m_tree_AS->Branch("jStr_evtEJ_B12", &_jStr_evtEJ_B12);
  m_tree_AS->Branch("jStr_evtEJ_BT12", &_jStr_evtEJ_BT12);
  m_tree_AS->Branch("jStr_evtEJ_Alpha", &_jStr_evtEJ_Alpha);
  m_tree_AS->Branch("jStr_evtEJ_AlphaT", &_jStr_evtEJ_AlphaT);
  m_tree_AS->Branch("jStr_evtEJ_massDemocracy", &_jStr_evtEJ_massDemocracy);
  m_tree_AS->Branch("jStr_evtEJ_betaflow", &_jStr_evtEJ_betaflow);
  m_tree_AS->Branch("jStr_evtEJ_betaflow_1GeV", &_jStr_evtEJ_betaflow_1GeV);
  m_tree_AS->Branch("jStr_evtEJ_betaflow_5GeV", &_jStr_evtEJ_betaflow_5GeV);
  m_tree_AS->Branch("jStr_evtEJ_y23", &_jStr_evtEJ_y23);
  m_tree_AS->Branch("jStr_evtEJ_y23_1GeV", &_jStr_evtEJ_y23_1GeV);
  m_tree_AS->Branch("jStr_evtEJ_y23_5GeV", &_jStr_evtEJ_y23_5GeV);
  m_tree_AS->Branch("jStr_evtEJ_lny23", &_jStr_evtEJ_lny23);
  m_tree_AS->Branch("jStr_evtEJ_lny23_1GeV", &_jStr_evtEJ_lny23_1GeV);
  m_tree_AS->Branch("jStr_evtEJ_lny23_5GeV", &_jStr_evtEJ_lny23_5GeV);
  m_tree_AS->Branch("jStr_evtEJ_theta", &_jStr_evtEJ_theta);
  m_tree_AS->Branch("jStr_evtEJ_asym", &_jStr_evtEJ_asym);
  m_tree_AS->Branch("jStr_evtEJ_yB", &_jStr_evtEJ_yB);
  m_tree_AS->Branch("jStr_evtEJ_yStar", &_jStr_evtEJ_yStar);
  m_tree_AS->Branch("jStr_evtEJ_thetaStar", &_jStr_evtEJ_thetaStar);
  m_tree_AS->Branch("jStr_evtEJ_chi", &_jStr_evtEJ_chi);
  m_tree_AS->Branch("jStr_evtEJ_deltaPhiJJ", &_jStr_evtEJ_deltaPhiJJ);
  m_tree_AS->Branch("jStr_evtEJ_deltaThetaJJ", &_jStr_evtEJ_deltaThetaJJ);
  m_tree_AS->Branch("jStr_evtEJ_deltaEtaJJ", &_jStr_evtEJ_deltaEtaJJ);
  m_tree_AS->Branch("jStr_evtEJ_deltaRapidityJJ", &_jStr_evtEJ_deltaRapidityJJ);
  m_tree_AS->Branch("jStr_evtEJ_deltaRJJ", &_jStr_evtEJ_deltaRJJ);
  m_tree_AS->Branch("jStr_evtEJ_deltaRJJY", &_jStr_evtEJ_deltaRJJY);
  m_tree_AS->Branch("jStr_evtEJ_sigmaPhiJJ", &_jStr_evtEJ_sigmaPhiJJ);
  m_tree_AS->Branch("jStr_evtEJ_sigmaThetaJJ", &_jStr_evtEJ_sigmaThetaJJ);
  m_tree_AS->Branch("jStr_evtEJ_sigmaEtaJJ", &_jStr_evtEJ_sigmaEtaJJ);
  m_tree_AS->Branch("jStr_evtEJ_sigmaRapidityJJ", &_jStr_evtEJ_sigmaRapidityJJ);
  m_tree_AS->Branch("jStr_evtEJ_sigmaPtJJ", &_jStr_evtEJ_sigmaPtJJ);
  m_tree_AS->Branch("jStr_evtEJ_sigmaEtJJ", &_jStr_evtEJ_sigmaEtJJ);
  m_tree_AS->Branch("jStr_evtEJ_sigmaEt12", &_jStr_evtEJ_sigmaEt12);
  m_tree_AS->Branch("jStr_evtEJ_sigmaEt34", &_jStr_evtEJ_sigmaEt34);
  m_tree_AS->Branch("jStr_evtEJ_A234", &_jStr_evtEJ_A234);
  m_tree_AS->Branch("jStr_evtEJ_asymPhiJJ", &_jStr_evtEJ_asymPhiJJ);
  m_tree_AS->Branch("jStr_evtEJ_asymThetaJJ", &_jStr_evtEJ_asymThetaJJ);
  m_tree_AS->Branch("jStr_evtEJ_asymEtaJJ", &_jStr_evtEJ_asymEtaJJ);
  m_tree_AS->Branch("jStr_evtEJ_asymRapidityJJ", &_jStr_evtEJ_asymRapidityJJ);
  m_tree_AS->Branch("jStr_evtEJ_acoplanarity", &_jStr_evtEJ_acoplanarity);
  m_tree_AS->Branch("jStr_evtEJ_twist", &_jStr_evtEJ_twist);
  m_tree_AS->Branch("jStr_evtEJ_twistY", &_jStr_evtEJ_twistY);
  m_tree_AS->Branch("jStr_evtEJ_jetSumE", &_jStr_evtEJ_jetSumE);
  m_tree_AS->Branch("jStr_evtEJ_jetSumET", &_jStr_evtEJ_jetSumET);
  m_tree_AS->Branch("jStr_evtEJ_jetSumPT", &_jStr_evtEJ_jetSumPT);
  m_tree_AS->Branch("jStr_evtEJ_jetSumM", &_jStr_evtEJ_jetSumM);
  m_tree_AS->Branch("jStr_evtEJ_jetSumMT", &_jStr_evtEJ_jetSumMT);
  m_tree_AS->Branch("jStr_evtEJ_HTprime", &_jStr_evtEJ_HTprime);
  m_tree_AS->Branch("jStr_evtEJ_centrality", &_jStr_evtEJ_centrality);
  m_tree_AS->Branch("jStr_evtEJ_centralityP", &_jStr_evtEJ_centralityP);
  m_tree_AS->Branch("jStr_evtEJ_zminJ1J2", &_jStr_evtEJ_zminJ1J2);
  m_tree_AS->Branch("jStr_evtEJ_zmaxJ1J2", &_jStr_evtEJ_zmaxJ1J2);
  m_tree_AS->Branch("jStr_evtEJ_zminAllJets", &_jStr_evtEJ_zminAllJets);
  m_tree_AS->Branch("jStr_evtEJ_zmaxAllJets", &_jStr_evtEJ_zmaxAllJets);
  m_tree_AS->Branch("jStr_evtEJ_zminJ1J2Phi", &_jStr_evtEJ_zminJ1J2Phi);
  m_tree_AS->Branch("jStr_evtEJ_zminJ1J2Theta", &_jStr_evtEJ_zminJ1J2Theta);
  m_tree_AS->Branch("jStr_evtEJ_zminJ1J2Eta", &_jStr_evtEJ_zminJ1J2Eta);
  m_tree_AS->Branch("jStr_evtEJ_zminJ1J2Rapidity", &_jStr_evtEJ_zminJ1J2Rapidity);
  m_tree_AS->Branch("jStr_evtEJ_cosHelicityJ1", &_jStr_evtEJ_cosHelicityJ1);
  m_tree_AS->Branch("jStr_evtEJ_helicityJ1", &_jStr_evtEJ_helicityJ1);
  m_tree_AS->Branch("jStr_evtEJ_azilicityJ1", &_jStr_evtEJ_azilicityJ1);
  m_tree_AS->Branch("jStr_evtEJ_cosHelicityJ2", &_jStr_evtEJ_cosHelicityJ2);
  m_tree_AS->Branch("jStr_evtEJ_helicityJ2", &_jStr_evtEJ_helicityJ2);
  m_tree_AS->Branch("jStr_evtEJ_azilicityJ2", &_jStr_evtEJ_azilicityJ2);
  m_tree_AS->Branch("jStr_evtEJ_cosThetaJ1", &_jStr_evtEJ_cosThetaJ1);
  m_tree_AS->Branch("jStr_evtEJ_cosThetaJ2", &_jStr_evtEJ_cosThetaJ2);
  m_tree_AS->Branch("jStr_evtEJ_deltaRapidityXtoJ1CM", &_jStr_evtEJ_deltaRapidityXtoJ1CM);
  m_tree_AS->Branch("jStr_evtEJ_deltaRapidityXtoJ2CM", &_jStr_evtEJ_deltaRapidityXtoJ2CM);
  m_tree_AS->Branch("jStr_evtEJ_deltaRapidityXtoJ1", &_jStr_evtEJ_deltaRapidityXtoJ1);
  m_tree_AS->Branch("jStr_evtEJ_deltaRapidityXtoJ2", &_jStr_evtEJ_deltaRapidityXtoJ2);
  m_tree_AS->Branch("jStr_evtEJ_cosTheta1", &_jStr_evtEJ_cosTheta1);
  m_tree_AS->Branch("jStr_evtEJ_cosTheta2", &_jStr_evtEJ_cosTheta2);
  m_tree_AS->Branch("jStr_evtEJ_cosThetaStar1", &_jStr_evtEJ_cosThetaStar1);
  m_tree_AS->Branch("jStr_evtEJ_cosThetaStar2", &_jStr_evtEJ_cosThetaStar2);
  m_tree_AS->Branch("jStr_evtEJ_cosPhiTilde1", &_jStr_evtEJ_cosPhiTilde1);
  m_tree_AS->Branch("jStr_evtEJ_cosPhiTilde2", &_jStr_evtEJ_cosPhiTilde2);
  m_tree_AS->Branch("jStr_evtEJ_cosPhi", &_jStr_evtEJ_cosPhi);
  m_tree_AS->Branch("jStr_evtEJ_PhiTilde1", &_jStr_evtEJ_PhiTilde1);
  m_tree_AS->Branch("jStr_evtEJ_PhiTilde2", &_jStr_evtEJ_PhiTilde2);
  m_tree_AS->Branch("jStr_evtEJ_Phi", &_jStr_evtEJ_Phi);
  m_tree_AS->Branch("jStr_evtEJ_M", &_jStr_evtEJ_M);
  m_tree_AS->Branch("jStr_evtEJ_DeltaM", &_jStr_evtEJ_DeltaM);
  m_tree_AS->Branch("jStr_evtEJ_asymM", &_jStr_evtEJ_asymM);
  m_tree_AS->Branch("jStr_evtEJ_Q", &_jStr_evtEJ_Q);
  m_tree_AS->Branch("jStr_evtEJ_SCDF", &_jStr_evtEJ_SCDF);
  m_tree_AS->Branch("jStr_evtEJ_dPhiIJCDF", &_jStr_evtEJ_dPhiIJCDF);
  m_tree_AS->Branch("jStr_evtEJ_dPhiKLCDF", &_jStr_evtEJ_dPhiKLCDF);
  m_tree_AS->Branch("jStr_evtEJ_PtIPtJCDF", &_jStr_evtEJ_PtIPtJCDF);
  m_tree_AS->Branch("jStr_evtEJ_PtKPtLCDF", &_jStr_evtEJ_PtKPtLCDF);
  m_tree_AS->Branch("jStr_evtEJ_SD0", &_jStr_evtEJ_SD0);
  m_tree_AS->Branch("jStr_evtEJ_dPhiIJD0", &_jStr_evtEJ_dPhiIJD0);
  m_tree_AS->Branch("jStr_evtEJ_dPhiKLD0", &_jStr_evtEJ_dPhiKLD0);
  m_tree_AS->Branch("jStr_evtEJ_PtIPtJD0", &_jStr_evtEJ_PtIPtJD0);
  m_tree_AS->Branch("jStr_evtEJ_PtKPtLD0", &_jStr_evtEJ_PtKPtLD0);
  m_tree_AS->Branch("jStr_evtEJ_DeltaSCDF", &_jStr_evtEJ_DeltaSCDF);
  m_tree_AS->Branch("jStr_evtEJ_DeltaSD0", &_jStr_evtEJ_DeltaSD0);
  m_tree_AS->Branch("jStr_evtEJ_Delta12", &_jStr_evtEJ_Delta12);
  m_tree_AS->Branch("jStr_evtEJ_Delta12n", &_jStr_evtEJ_Delta12n);
  m_tree_AS->Branch("jStr_evtEJ_M14", &_jStr_evtEJ_M14);
  m_tree_AS->Branch("jStr_evtEJ_y1y2", &_jStr_evtEJ_y1y2);
  m_tree_AS->Branch("jStr_evtEJ_Mprime4", &_jStr_evtEJ_Mprime4);
  m_tree_AS->Branch("jStr_evtEJ_dMprime4", &_jStr_evtEJ_dMprime4);
  m_tree_AS->Branch("jStr_evtEJ_MprimeAvg4", &_jStr_evtEJ_MprimeAvg4);
  m_tree_AS->Branch("jStr_evtEJ_DeltaMin4", &_jStr_evtEJ_DeltaMin4);
  m_tree_AS->Branch("jStr_evtEJ_DeltaMax4", &_jStr_evtEJ_DeltaMax4);
  m_tree_AS->Branch("jStr_evtEJ_DeltaPhiXX4", &_jStr_evtEJ_DeltaPhiXX4);
  m_tree_AS->Branch("jStr_evtEJ_DeltaYXX4", &_jStr_evtEJ_DeltaYXX4);
  m_tree_AS->Branch("jStr_evtEJ_DeltaRXX4", &_jStr_evtEJ_DeltaRXX4);
  m_tree_AS->Branch("jStr_evtEJ_TwistXX4", &_jStr_evtEJ_TwistXX4);
  m_tree_AS->Branch("jStr_evtEJ_separated4", &_jStr_evtEJ_separated4);
  m_tree_AS->Branch("jStr_evtEJ_Mprime", &_jStr_evtEJ_Mprime);
  m_tree_AS->Branch("jStr_evtEJ_dMprime", &_jStr_evtEJ_dMprime);
  m_tree_AS->Branch("jStr_evtEJ_MprimeAvg", &_jStr_evtEJ_MprimeAvg);
  m_tree_AS->Branch("jStr_evtEJ_DeltaMin", &_jStr_evtEJ_DeltaMin);
  m_tree_AS->Branch("jStr_evtEJ_DeltaMax", &_jStr_evtEJ_DeltaMax);
  m_tree_AS->Branch("jStr_evtEJ_DeltaPhiXX", &_jStr_evtEJ_DeltaPhiXX);
  m_tree_AS->Branch("jStr_evtEJ_DeltaYXX", &_jStr_evtEJ_DeltaYXX);
  m_tree_AS->Branch("jStr_evtEJ_DeltaRXX", &_jStr_evtEJ_DeltaRXX);
  m_tree_AS->Branch("jStr_evtEJ_TwistXX", &_jStr_evtEJ_TwistXX);
  m_tree_AS->Branch("jStr_evtEJ_separated", &_jStr_evtEJ_separated);

  m_tree_AS->Branch("jStr_evtJj_etaJ1", &_jStr_evtJj_etaJ1);
  m_tree_AS->Branch("jStr_evtJj_etaJ2", &_jStr_evtJj_etaJ2);
  m_tree_AS->Branch("jStr_evtJj_etaJ3", &_jStr_evtJj_etaJ3);
  m_tree_AS->Branch("jStr_evtJj_etaJ4", &_jStr_evtJj_etaJ4);
  m_tree_AS->Branch("jStr_evtJj_pTJ1", &_jStr_evtJj_pTJ1);
  m_tree_AS->Branch("jStr_evtJj_pTJ2", &_jStr_evtJj_pTJ2);
  m_tree_AS->Branch("jStr_evtJj_pTJ3", &_jStr_evtJj_pTJ3);
  m_tree_AS->Branch("jStr_evtJj_pTJ4", &_jStr_evtJj_pTJ4);

  m_tree_AS->Branch("jStr_evtJj_fE1", &_jStr_evtJj_fE1);
  m_tree_AS->Branch("jStr_evtJj_fE2", &_jStr_evtJj_fE2);
  m_tree_AS->Branch("jStr_evtJj_fE3", &_jStr_evtJj_fE3);
  m_tree_AS->Branch("jStr_evtJj_fET1", &_jStr_evtJj_fET1);
  m_tree_AS->Branch("jStr_evtJj_fET2", &_jStr_evtJj_fET2);
  m_tree_AS->Branch("jStr_evtJj_fET3", &_jStr_evtJj_fET3);
  m_tree_AS->Branch("jStr_evtJj_mct", &_jStr_evtJj_mct);
  m_tree_AS->Branch("jStr_evtJj_q", &_jStr_evtJj_q);
  m_tree_AS->Branch("jStr_evtJj_mjj", &_jStr_evtJj_mjj);
  m_tree_AS->Branch("jStr_evtJj_mTjj", &_jStr_evtJj_mTjj);
  m_tree_AS->Branch("jStr_evtJj_mjjj", &_jStr_evtJj_mjjj);
  m_tree_AS->Branch("jStr_evtJj_mTjjj", &_jStr_evtJj_mTjjj);
  m_tree_AS->Branch("jStr_evtJj_mjjjj", &_jStr_evtJj_mjjjj);
  m_tree_AS->Branch("jStr_evtJj_mTjjjj", &_jStr_evtJj_mTjjjj);
  m_tree_AS->Branch("jStr_evtJj_zetaPlus", &_jStr_evtJj_zetaPlus);
  m_tree_AS->Branch("jStr_evtJj_zetaMinus", &_jStr_evtJj_zetaMinus);
  m_tree_AS->Branch("jStr_evtJj_B12", &_jStr_evtJj_B12);
  m_tree_AS->Branch("jStr_evtJj_BT12", &_jStr_evtJj_BT12);
  m_tree_AS->Branch("jStr_evtJj_Alpha", &_jStr_evtJj_Alpha);
  m_tree_AS->Branch("jStr_evtJj_AlphaT", &_jStr_evtJj_AlphaT);
  m_tree_AS->Branch("jStr_evtJj_massDemocracy", &_jStr_evtJj_massDemocracy);
  m_tree_AS->Branch("jStr_evtJj_betaflow", &_jStr_evtJj_betaflow);
  m_tree_AS->Branch("jStr_evtJj_betaflow_1GeV", &_jStr_evtJj_betaflow_1GeV);
  m_tree_AS->Branch("jStr_evtJj_betaflow_5GeV", &_jStr_evtJj_betaflow_5GeV);
  m_tree_AS->Branch("jStr_evtJj_y23", &_jStr_evtJj_y23);
  m_tree_AS->Branch("jStr_evtJj_y23_1GeV", &_jStr_evtJj_y23_1GeV);
  m_tree_AS->Branch("jStr_evtJj_y23_5GeV", &_jStr_evtJj_y23_5GeV);
  m_tree_AS->Branch("jStr_evtJj_lny23", &_jStr_evtJj_lny23);
  m_tree_AS->Branch("jStr_evtJj_lny23_1GeV", &_jStr_evtJj_lny23_1GeV);
  m_tree_AS->Branch("jStr_evtJj_lny23_5GeV", &_jStr_evtJj_lny23_5GeV);
  m_tree_AS->Branch("jStr_evtJj_theta", &_jStr_evtJj_theta);
  m_tree_AS->Branch("jStr_evtJj_asym", &_jStr_evtJj_asym);
  m_tree_AS->Branch("jStr_evtJj_yB", &_jStr_evtJj_yB);
  m_tree_AS->Branch("jStr_evtJj_yStar", &_jStr_evtJj_yStar);
  m_tree_AS->Branch("jStr_evtJj_thetaStar", &_jStr_evtJj_thetaStar);
  m_tree_AS->Branch("jStr_evtJj_chi", &_jStr_evtJj_chi);
  m_tree_AS->Branch("jStr_evtJj_deltaPhiJJ", &_jStr_evtJj_deltaPhiJJ);
  m_tree_AS->Branch("jStr_evtJj_deltaThetaJJ", &_jStr_evtJj_deltaThetaJJ);
  m_tree_AS->Branch("jStr_evtJj_deltaEtaJJ", &_jStr_evtJj_deltaEtaJJ);
  m_tree_AS->Branch("jStr_evtJj_deltaRapidityJJ", &_jStr_evtJj_deltaRapidityJJ);
  m_tree_AS->Branch("jStr_evtJj_deltaRJJ", &_jStr_evtJj_deltaRJJ);
  m_tree_AS->Branch("jStr_evtJj_deltaRJJY", &_jStr_evtJj_deltaRJJY);
  m_tree_AS->Branch("jStr_evtJj_sigmaPhiJJ", &_jStr_evtJj_sigmaPhiJJ);
  m_tree_AS->Branch("jStr_evtJj_sigmaThetaJJ", &_jStr_evtJj_sigmaThetaJJ);
  m_tree_AS->Branch("jStr_evtJj_sigmaEtaJJ", &_jStr_evtJj_sigmaEtaJJ);
  m_tree_AS->Branch("jStr_evtJj_sigmaRapidityJJ", &_jStr_evtJj_sigmaRapidityJJ);
  m_tree_AS->Branch("jStr_evtJj_sigmaPtJJ", &_jStr_evtJj_sigmaPtJJ);
  m_tree_AS->Branch("jStr_evtJj_sigmaEtJJ", &_jStr_evtJj_sigmaEtJJ);
  m_tree_AS->Branch("jStr_evtJj_sigmaEt12", &_jStr_evtJj_sigmaEt12);
  m_tree_AS->Branch("jStr_evtJj_sigmaEt34", &_jStr_evtJj_sigmaEt34);
  m_tree_AS->Branch("jStr_evtJj_A234", &_jStr_evtJj_A234);
  m_tree_AS->Branch("jStr_evtJj_asymPhiJJ", &_jStr_evtJj_asymPhiJJ);
  m_tree_AS->Branch("jStr_evtJj_asymThetaJJ", &_jStr_evtJj_asymThetaJJ);
  m_tree_AS->Branch("jStr_evtJj_asymEtaJJ", &_jStr_evtJj_asymEtaJJ);
  m_tree_AS->Branch("jStr_evtJj_asymRapidityJJ", &_jStr_evtJj_asymRapidityJJ);
  m_tree_AS->Branch("jStr_evtJj_acoplanarity", &_jStr_evtJj_acoplanarity);
  m_tree_AS->Branch("jStr_evtJj_twist", &_jStr_evtJj_twist);
  m_tree_AS->Branch("jStr_evtJj_twistY", &_jStr_evtJj_twistY);
  m_tree_AS->Branch("jStr_evtJj_jetSumE", &_jStr_evtJj_jetSumE);
  m_tree_AS->Branch("jStr_evtJj_jetSumET", &_jStr_evtJj_jetSumET);
  m_tree_AS->Branch("jStr_evtJj_jetSumPT", &_jStr_evtJj_jetSumPT);
  m_tree_AS->Branch("jStr_evtJj_jetSumM", &_jStr_evtJj_jetSumM);
  m_tree_AS->Branch("jStr_evtJj_jetSumMT", &_jStr_evtJj_jetSumMT);
  m_tree_AS->Branch("jStr_evtJj_HTprime", &_jStr_evtJj_HTprime);
  m_tree_AS->Branch("jStr_evtJj_centrality", &_jStr_evtJj_centrality);
  m_tree_AS->Branch("jStr_evtJj_centralityP", &_jStr_evtJj_centralityP);
  m_tree_AS->Branch("jStr_evtJj_zminJ1J2", &_jStr_evtJj_zminJ1J2);
  m_tree_AS->Branch("jStr_evtJj_zmaxJ1J2", &_jStr_evtJj_zmaxJ1J2);
  m_tree_AS->Branch("jStr_evtJj_zminAllJets", &_jStr_evtJj_zminAllJets);
  m_tree_AS->Branch("jStr_evtJj_zmaxAllJets", &_jStr_evtJj_zmaxAllJets);
  m_tree_AS->Branch("jStr_evtJj_zminJ1J2Phi", &_jStr_evtJj_zminJ1J2Phi);
  m_tree_AS->Branch("jStr_evtJj_zminJ1J2Theta", &_jStr_evtJj_zminJ1J2Theta);
  m_tree_AS->Branch("jStr_evtJj_zminJ1J2Eta", &_jStr_evtJj_zminJ1J2Eta);
  m_tree_AS->Branch("jStr_evtJj_zminJ1J2Rapidity", &_jStr_evtJj_zminJ1J2Rapidity);
  m_tree_AS->Branch("jStr_evtJj_cosHelicityJ1", &_jStr_evtJj_cosHelicityJ1);
  m_tree_AS->Branch("jStr_evtJj_helicityJ1", &_jStr_evtJj_helicityJ1);
  m_tree_AS->Branch("jStr_evtJj_azilicityJ1", &_jStr_evtJj_azilicityJ1);
  m_tree_AS->Branch("jStr_evtJj_cosHelicityJ2", &_jStr_evtJj_cosHelicityJ2);
  m_tree_AS->Branch("jStr_evtJj_helicityJ2", &_jStr_evtJj_helicityJ2);
  m_tree_AS->Branch("jStr_evtJj_azilicityJ2", &_jStr_evtJj_azilicityJ2);
  m_tree_AS->Branch("jStr_evtJj_cosThetaJ1", &_jStr_evtJj_cosThetaJ1);
  m_tree_AS->Branch("jStr_evtJj_cosThetaJ2", &_jStr_evtJj_cosThetaJ2);
  m_tree_AS->Branch("jStr_evtJj_deltaRapidityXtoJ1CM", &_jStr_evtJj_deltaRapidityXtoJ1CM);
  m_tree_AS->Branch("jStr_evtJj_deltaRapidityXtoJ2CM", &_jStr_evtJj_deltaRapidityXtoJ2CM);
  m_tree_AS->Branch("jStr_evtJj_deltaRapidityXtoJ1", &_jStr_evtJj_deltaRapidityXtoJ1);
  m_tree_AS->Branch("jStr_evtJj_deltaRapidityXtoJ2", &_jStr_evtJj_deltaRapidityXtoJ2);
  m_tree_AS->Branch("jStr_evtJj_cosTheta1", &_jStr_evtJj_cosTheta1);
  m_tree_AS->Branch("jStr_evtJj_cosTheta2", &_jStr_evtJj_cosTheta2);
  m_tree_AS->Branch("jStr_evtJj_cosThetaStar1", &_jStr_evtJj_cosThetaStar1);
  m_tree_AS->Branch("jStr_evtJj_cosThetaStar2", &_jStr_evtJj_cosThetaStar2);
  m_tree_AS->Branch("jStr_evtJj_cosPhiTilde1", &_jStr_evtJj_cosPhiTilde1);
  m_tree_AS->Branch("jStr_evtJj_cosPhiTilde2", &_jStr_evtJj_cosPhiTilde2);
  m_tree_AS->Branch("jStr_evtJj_cosPhi", &_jStr_evtJj_cosPhi);
  m_tree_AS->Branch("jStr_evtJj_PhiTilde1", &_jStr_evtJj_PhiTilde1);
  m_tree_AS->Branch("jStr_evtJj_PhiTilde2", &_jStr_evtJj_PhiTilde2);
  m_tree_AS->Branch("jStr_evtJj_Phi", &_jStr_evtJj_Phi);
  m_tree_AS->Branch("jStr_evtJj_M", &_jStr_evtJj_M);
  m_tree_AS->Branch("jStr_evtJj_DeltaM", &_jStr_evtJj_DeltaM);
  m_tree_AS->Branch("jStr_evtJj_asymM", &_jStr_evtJj_asymM);
  m_tree_AS->Branch("jStr_evtJj_Q", &_jStr_evtJj_Q);
  m_tree_AS->Branch("jStr_evtJj_SCDF", &_jStr_evtJj_SCDF);
  m_tree_AS->Branch("jStr_evtJj_dPhiIJCDF", &_jStr_evtJj_dPhiIJCDF);
  m_tree_AS->Branch("jStr_evtJj_dPhiKLCDF", &_jStr_evtJj_dPhiKLCDF);
  m_tree_AS->Branch("jStr_evtJj_PtIPtJCDF", &_jStr_evtJj_PtIPtJCDF);
  m_tree_AS->Branch("jStr_evtJj_PtKPtLCDF", &_jStr_evtJj_PtKPtLCDF);
  m_tree_AS->Branch("jStr_evtJj_SD0", &_jStr_evtJj_SD0);
  m_tree_AS->Branch("jStr_evtJj_dPhiIJD0", &_jStr_evtJj_dPhiIJD0);
  m_tree_AS->Branch("jStr_evtJj_dPhiKLD0", &_jStr_evtJj_dPhiKLD0);
  m_tree_AS->Branch("jStr_evtJj_PtIPtJD0", &_jStr_evtJj_PtIPtJD0);
  m_tree_AS->Branch("jStr_evtJj_PtKPtLD0", &_jStr_evtJj_PtKPtLD0);
  m_tree_AS->Branch("jStr_evtJj_DeltaSCDF", &_jStr_evtJj_DeltaSCDF);
  m_tree_AS->Branch("jStr_evtJj_DeltaSD0", &_jStr_evtJj_DeltaSD0);
  m_tree_AS->Branch("jStr_evtJj_Delta12", &_jStr_evtJj_Delta12);
  m_tree_AS->Branch("jStr_evtJj_Delta12n", &_jStr_evtJj_Delta12n);
  m_tree_AS->Branch("jStr_evtJj_M14", &_jStr_evtJj_M14);
  m_tree_AS->Branch("jStr_evtJj_y1y2", &_jStr_evtJj_y1y2);
  m_tree_AS->Branch("jStr_evtJj_Mprime4", &_jStr_evtJj_Mprime4);
  m_tree_AS->Branch("jStr_evtJj_dMprime4", &_jStr_evtJj_dMprime4);
  m_tree_AS->Branch("jStr_evtJj_MprimeAvg4", &_jStr_evtJj_MprimeAvg4);
  m_tree_AS->Branch("jStr_evtJj_DeltaMin4", &_jStr_evtJj_DeltaMin4);
  m_tree_AS->Branch("jStr_evtJj_DeltaMax4", &_jStr_evtJj_DeltaMax4);
  m_tree_AS->Branch("jStr_evtJj_DeltaPhiXX4", &_jStr_evtJj_DeltaPhiXX4);
  m_tree_AS->Branch("jStr_evtJj_DeltaYXX4", &_jStr_evtJj_DeltaYXX4);
  m_tree_AS->Branch("jStr_evtJj_DeltaRXX4", &_jStr_evtJj_DeltaRXX4);
  m_tree_AS->Branch("jStr_evtJj_TwistXX4", &_jStr_evtJj_TwistXX4);
  m_tree_AS->Branch("jStr_evtJj_separated4", &_jStr_evtJj_separated4);
  m_tree_AS->Branch("jStr_evtJj_Mprime", &_jStr_evtJj_Mprime);
  m_tree_AS->Branch("jStr_evtJj_dMprime", &_jStr_evtJj_dMprime);
  m_tree_AS->Branch("jStr_evtJj_MprimeAvg", &_jStr_evtJj_MprimeAvg);
  m_tree_AS->Branch("jStr_evtJj_DeltaMin", &_jStr_evtJj_DeltaMin);
  m_tree_AS->Branch("jStr_evtJj_DeltaMax", &_jStr_evtJj_DeltaMax);
  m_tree_AS->Branch("jStr_evtJj_DeltaPhiXX", &_jStr_evtJj_DeltaPhiXX);
  m_tree_AS->Branch("jStr_evtJj_DeltaYXX", &_jStr_evtJj_DeltaYXX);
  m_tree_AS->Branch("jStr_evtJj_DeltaRXX", &_jStr_evtJj_DeltaRXX);
  m_tree_AS->Branch("jStr_evtJj_TwistXX", &_jStr_evtJj_TwistXX);
  m_tree_AS->Branch("jStr_evtJj_separated", &_jStr_evtJj_separated);
  m_tree_AS->Branch("jStr_fwmEJ_xiPlus", &_jStr_fwmEJ_xiPlus);
  m_tree_AS->Branch("jStr_fwmEJ_xiMinus", &_jStr_fwmEJ_xiMinus);
  m_tree_AS->Branch("jStr_fwmEJ_xPlus", &_jStr_fwmEJ_xPlus);
  m_tree_AS->Branch("jStr_fwmEJ_xMinus", &_jStr_fwmEJ_xMinus);
  m_tree_AS->Branch("jStr_fwmEJ_Psi1", &_jStr_fwmEJ_Psi1);
  m_tree_AS->Branch("jStr_fwmEJ_B0", &_jStr_fwmEJ_B0);
  m_tree_AS->Branch("jStr_fwmEJ_B1", &_jStr_fwmEJ_B1);
  m_tree_AS->Branch("jStr_fwmEJ_B2", &_jStr_fwmEJ_B2);
  m_tree_AS->Branch("jStr_fwmEJ_B3", &_jStr_fwmEJ_B3);
  m_tree_AS->Branch("jStr_fwmEJ_B4", &_jStr_fwmEJ_B4);
  m_tree_AS->Branch("jStr_fwmEJ_B5", &_jStr_fwmEJ_B5);
  m_tree_AS->Branch("jStr_fwmEJ_B6", &_jStr_fwmEJ_B6);
  m_tree_AS->Branch("jStr_fwmEJ_B7", &_jStr_fwmEJ_B7);
  m_tree_AS->Branch("jStr_fwmEJ_B8", &_jStr_fwmEJ_B8);
  m_tree_AS->Branch("jStr_fwmEJ_C0", &_jStr_fwmEJ_C0);
  m_tree_AS->Branch("jStr_fwmEJ_C1", &_jStr_fwmEJ_C1);
  m_tree_AS->Branch("jStr_fwmEJ_C2", &_jStr_fwmEJ_C2);
  m_tree_AS->Branch("jStr_fwmEJ_C3", &_jStr_fwmEJ_C3);
  m_tree_AS->Branch("jStr_fwmEJ_C4", &_jStr_fwmEJ_C4);
  m_tree_AS->Branch("jStr_fwmEJ_C5", &_jStr_fwmEJ_C5);
  m_tree_AS->Branch("jStr_fwmEJ_C6", &_jStr_fwmEJ_C6);
  m_tree_AS->Branch("jStr_fwmEJ_C7", &_jStr_fwmEJ_C7);
  m_tree_AS->Branch("jStr_fwmEJ_C8", &_jStr_fwmEJ_C8);
  m_tree_AS->Branch("jStr_fwmEJ_K0", &_jStr_fwmEJ_K0);
  m_tree_AS->Branch("jStr_fwmEJ_K1", &_jStr_fwmEJ_K1);
  m_tree_AS->Branch("jStr_fwmEJ_K2", &_jStr_fwmEJ_K2);
  m_tree_AS->Branch("jStr_fwmEJ_K3", &_jStr_fwmEJ_K3);
  m_tree_AS->Branch("jStr_fwmEJ_K4", &_jStr_fwmEJ_K4);
  m_tree_AS->Branch("jStr_fwmEJ_K5", &_jStr_fwmEJ_K5);
  m_tree_AS->Branch("jStr_fwmEJ_K6", &_jStr_fwmEJ_K6);
  m_tree_AS->Branch("jStr_fwmEJ_K7", &_jStr_fwmEJ_K7);
  m_tree_AS->Branch("jStr_fwmEJ_K8", &_jStr_fwmEJ_K8);
  m_tree_AS->Branch("jStr_fwmEJ_D0", &_jStr_fwmEJ_D0);
  m_tree_AS->Branch("jStr_fwmEJ_D1", &_jStr_fwmEJ_D1);
  m_tree_AS->Branch("jStr_fwmEJ_D2", &_jStr_fwmEJ_D2);
  m_tree_AS->Branch("jStr_fwmEJ_D3", &_jStr_fwmEJ_D3);
  m_tree_AS->Branch("jStr_fwmEJ_D4", &_jStr_fwmEJ_D4);
  m_tree_AS->Branch("jStr_fwmEJ_D5", &_jStr_fwmEJ_D5);
  m_tree_AS->Branch("jStr_fwmEJ_D6", &_jStr_fwmEJ_D6);
  m_tree_AS->Branch("jStr_fwmEJ_D7", &_jStr_fwmEJ_D7);
  m_tree_AS->Branch("jStr_fwmEJ_D8", &_jStr_fwmEJ_D8);
  m_tree_AS->Branch("jStr_fwmEJ_H0", &_jStr_fwmEJ_H0);
  m_tree_AS->Branch("jStr_fwmEJ_H1", &_jStr_fwmEJ_H1);
  m_tree_AS->Branch("jStr_fwmEJ_H2", &_jStr_fwmEJ_H2);
  m_tree_AS->Branch("jStr_fwmEJ_H3", &_jStr_fwmEJ_H3);
  m_tree_AS->Branch("jStr_fwmEJ_H4", &_jStr_fwmEJ_H4);
  m_tree_AS->Branch("jStr_fwmEJ_H5", &_jStr_fwmEJ_H5);
  m_tree_AS->Branch("jStr_fwmEJ_H6", &_jStr_fwmEJ_H6);
  m_tree_AS->Branch("jStr_fwmEJ_H7", &_jStr_fwmEJ_H7);
  m_tree_AS->Branch("jStr_fwmEJ_H8", &_jStr_fwmEJ_H8);
  m_tree_AS->Branch("jStr_fwmEJ_Q0", &_jStr_fwmEJ_Q0);
  m_tree_AS->Branch("jStr_fwmEJ_Q1", &_jStr_fwmEJ_Q1);
  m_tree_AS->Branch("jStr_fwmEJ_Q2", &_jStr_fwmEJ_Q2);
  m_tree_AS->Branch("jStr_fwmEJ_Q3", &_jStr_fwmEJ_Q3);
  m_tree_AS->Branch("jStr_fwmEJ_Q4", &_jStr_fwmEJ_Q4);
  m_tree_AS->Branch("jStr_fwmEJ_Q5", &_jStr_fwmEJ_Q5);
  m_tree_AS->Branch("jStr_fwmEJ_Q6", &_jStr_fwmEJ_Q6);
  m_tree_AS->Branch("jStr_fwmEJ_Q7", &_jStr_fwmEJ_Q7);
  m_tree_AS->Branch("jStr_fwmEJ_Q8", &_jStr_fwmEJ_Q8);
  m_tree_AS->Branch("jStr_fwmEJ_Pi1", &_jStr_fwmEJ_Pi1);
  m_tree_AS->Branch("jStr_fwmEJ_Pi2", &_jStr_fwmEJ_Pi2);
  m_tree_AS->Branch("jStr_fwmEJ_Pi3", &_jStr_fwmEJ_Pi3);
  m_tree_AS->Branch("jStr_fwmEJ_Pi4", &_jStr_fwmEJ_Pi4);
  m_tree_AS->Branch("jStr_fwmEJ_B10", &_jStr_fwmEJ_B10);
  m_tree_AS->Branch("jStr_fwmEJ_B20", &_jStr_fwmEJ_B20);
  m_tree_AS->Branch("jStr_fwmEJ_B30", &_jStr_fwmEJ_B30);
  m_tree_AS->Branch("jStr_fwmEJ_B40", &_jStr_fwmEJ_B40);
  m_tree_AS->Branch("jStr_fwmEJ_B50", &_jStr_fwmEJ_B50);
  m_tree_AS->Branch("jStr_fwmEJ_B60", &_jStr_fwmEJ_B60);
  m_tree_AS->Branch("jStr_fwmEJ_B70", &_jStr_fwmEJ_B70);
  m_tree_AS->Branch("jStr_fwmEJ_B80", &_jStr_fwmEJ_B80);
  m_tree_AS->Branch("jStr_fwmEJ_C10", &_jStr_fwmEJ_C10);
  m_tree_AS->Branch("jStr_fwmEJ_C20", &_jStr_fwmEJ_C20);
  m_tree_AS->Branch("jStr_fwmEJ_C30", &_jStr_fwmEJ_C30);
  m_tree_AS->Branch("jStr_fwmEJ_C40", &_jStr_fwmEJ_C40);
  m_tree_AS->Branch("jStr_fwmEJ_C50", &_jStr_fwmEJ_C50);
  m_tree_AS->Branch("jStr_fwmEJ_C60", &_jStr_fwmEJ_C60);
  m_tree_AS->Branch("jStr_fwmEJ_C70", &_jStr_fwmEJ_C70);
  m_tree_AS->Branch("jStr_fwmEJ_C80", &_jStr_fwmEJ_C80);
  m_tree_AS->Branch("jStr_fwmEJ_K10", &_jStr_fwmEJ_K10);
  m_tree_AS->Branch("jStr_fwmEJ_K20", &_jStr_fwmEJ_K20);
  m_tree_AS->Branch("jStr_fwmEJ_K30", &_jStr_fwmEJ_K30);
  m_tree_AS->Branch("jStr_fwmEJ_K40", &_jStr_fwmEJ_K40);
  m_tree_AS->Branch("jStr_fwmEJ_K50", &_jStr_fwmEJ_K50);
  m_tree_AS->Branch("jStr_fwmEJ_K60", &_jStr_fwmEJ_K60);
  m_tree_AS->Branch("jStr_fwmEJ_K70", &_jStr_fwmEJ_K70);
  m_tree_AS->Branch("jStr_fwmEJ_K80", &_jStr_fwmEJ_K80);
  m_tree_AS->Branch("jStr_fwmEJ_D10", &_jStr_fwmEJ_D10);
  m_tree_AS->Branch("jStr_fwmEJ_D20", &_jStr_fwmEJ_D20);
  m_tree_AS->Branch("jStr_fwmEJ_D30", &_jStr_fwmEJ_D30);
  m_tree_AS->Branch("jStr_fwmEJ_D40", &_jStr_fwmEJ_D40);
  m_tree_AS->Branch("jStr_fwmEJ_D50", &_jStr_fwmEJ_D50);
  m_tree_AS->Branch("jStr_fwmEJ_D60", &_jStr_fwmEJ_D60);
  m_tree_AS->Branch("jStr_fwmEJ_D70", &_jStr_fwmEJ_D70);
  m_tree_AS->Branch("jStr_fwmEJ_D80", &_jStr_fwmEJ_D80);
  m_tree_AS->Branch("jStr_fwmEJ_H10", &_jStr_fwmEJ_H10);
  m_tree_AS->Branch("jStr_fwmEJ_H20", &_jStr_fwmEJ_H20);
  m_tree_AS->Branch("jStr_fwmEJ_H30", &_jStr_fwmEJ_H30);
  m_tree_AS->Branch("jStr_fwmEJ_H40", &_jStr_fwmEJ_H40);
  m_tree_AS->Branch("jStr_fwmEJ_H50", &_jStr_fwmEJ_H50);
  m_tree_AS->Branch("jStr_fwmEJ_H60", &_jStr_fwmEJ_H60);
  m_tree_AS->Branch("jStr_fwmEJ_H70", &_jStr_fwmEJ_H70);
  m_tree_AS->Branch("jStr_fwmEJ_H80", &_jStr_fwmEJ_H80);
  m_tree_AS->Branch("jStr_fwmEJ_Q10", &_jStr_fwmEJ_Q10);
  m_tree_AS->Branch("jStr_fwmEJ_Q20", &_jStr_fwmEJ_Q20);
  m_tree_AS->Branch("jStr_fwmEJ_Q30", &_jStr_fwmEJ_Q30);
  m_tree_AS->Branch("jStr_fwmEJ_Q40", &_jStr_fwmEJ_Q40);
  m_tree_AS->Branch("jStr_fwmEJ_Q50", &_jStr_fwmEJ_Q50);
  m_tree_AS->Branch("jStr_fwmEJ_Q60", &_jStr_fwmEJ_Q60);
  m_tree_AS->Branch("jStr_fwmEJ_Q70", &_jStr_fwmEJ_Q70);
  m_tree_AS->Branch("jStr_fwmEJ_Q80", &_jStr_fwmEJ_Q80);
  // m_tree_AS->Branch("jStr_fwmEi_Psi1", &_jStr_fwmEi_Psi1);
  // m_tree_AS->Branch("jStr_fwmEi_B0", &_jStr_fwmEi_B0);
  // m_tree_AS->Branch("jStr_fwmEi_B1", &_jStr_fwmEi_B1);
  // m_tree_AS->Branch("jStr_fwmEi_B2", &_jStr_fwmEi_B2);
  // m_tree_AS->Branch("jStr_fwmEi_B3", &_jStr_fwmEi_B3);
  // m_tree_AS->Branch("jStr_fwmEi_B4", &_jStr_fwmEi_B4);
  // m_tree_AS->Branch("jStr_fwmEi_B5", &_jStr_fwmEi_B5);
  // m_tree_AS->Branch("jStr_fwmEi_B6", &_jStr_fwmEi_B6);
  // m_tree_AS->Branch("jStr_fwmEi_B7", &_jStr_fwmEi_B7);
  // m_tree_AS->Branch("jStr_fwmEi_B8", &_jStr_fwmEi_B8);
  // m_tree_AS->Branch("jStr_fwmEi_C0", &_jStr_fwmEi_C0);
  // m_tree_AS->Branch("jStr_fwmEi_C1", &_jStr_fwmEi_C1);
  // m_tree_AS->Branch("jStr_fwmEi_C2", &_jStr_fwmEi_C2);
  // m_tree_AS->Branch("jStr_fwmEi_C3", &_jStr_fwmEi_C3);
  // m_tree_AS->Branch("jStr_fwmEi_C4", &_jStr_fwmEi_C4);
  // m_tree_AS->Branch("jStr_fwmEi_C5", &_jStr_fwmEi_C5);
  // m_tree_AS->Branch("jStr_fwmEi_C6", &_jStr_fwmEi_C6);
  // m_tree_AS->Branch("jStr_fwmEi_C7", &_jStr_fwmEi_C7);
  // m_tree_AS->Branch("jStr_fwmEi_C8", &_jStr_fwmEi_C8);
  // m_tree_AS->Branch("jStr_fwmEi_K0", &_jStr_fwmEi_K0);
  // m_tree_AS->Branch("jStr_fwmEi_K1", &_jStr_fwmEi_K1);
  // m_tree_AS->Branch("jStr_fwmEi_K2", &_jStr_fwmEi_K2);
  // m_tree_AS->Branch("jStr_fwmEi_K3", &_jStr_fwmEi_K3);
  // m_tree_AS->Branch("jStr_fwmEi_K4", &_jStr_fwmEi_K4);
  // m_tree_AS->Branch("jStr_fwmEi_K5", &_jStr_fwmEi_K5);
  // m_tree_AS->Branch("jStr_fwmEi_K6", &_jStr_fwmEi_K6);
  // m_tree_AS->Branch("jStr_fwmEi_K7", &_jStr_fwmEi_K7);
  // m_tree_AS->Branch("jStr_fwmEi_K8", &_jStr_fwmEi_K8);
  // m_tree_AS->Branch("jStr_fwmEi_D0", &_jStr_fwmEi_D0);
  // m_tree_AS->Branch("jStr_fwmEi_D1", &_jStr_fwmEi_D1);
  // m_tree_AS->Branch("jStr_fwmEi_D2", &_jStr_fwmEi_D2);
  // m_tree_AS->Branch("jStr_fwmEi_D3", &_jStr_fwmEi_D3);
  // m_tree_AS->Branch("jStr_fwmEi_D4", &_jStr_fwmEi_D4);
  // m_tree_AS->Branch("jStr_fwmEi_D5", &_jStr_fwmEi_D5);
  // m_tree_AS->Branch("jStr_fwmEi_D6", &_jStr_fwmEi_D6);
  // m_tree_AS->Branch("jStr_fwmEi_D7", &_jStr_fwmEi_D7);
  // m_tree_AS->Branch("jStr_fwmEi_D8", &_jStr_fwmEi_D8);
  // m_tree_AS->Branch("jStr_fwmEi_H0", &_jStr_fwmEi_H0);
  // m_tree_AS->Branch("jStr_fwmEi_H1", &_jStr_fwmEi_H1);
  // m_tree_AS->Branch("jStr_fwmEi_H2", &_jStr_fwmEi_H2);
  // m_tree_AS->Branch("jStr_fwmEi_H3", &_jStr_fwmEi_H3);
  // m_tree_AS->Branch("jStr_fwmEi_H4", &_jStr_fwmEi_H4);
  // m_tree_AS->Branch("jStr_fwmEi_H5", &_jStr_fwmEi_H5);
  // m_tree_AS->Branch("jStr_fwmEi_H6", &_jStr_fwmEi_H6);
  // m_tree_AS->Branch("jStr_fwmEi_H7", &_jStr_fwmEi_H7);
  // m_tree_AS->Branch("jStr_fwmEi_H8", &_jStr_fwmEi_H8);
  // m_tree_AS->Branch("jStr_fwmEi_Q0", &_jStr_fwmEi_Q0);
  // m_tree_AS->Branch("jStr_fwmEi_Q1", &_jStr_fwmEi_Q1);
  // m_tree_AS->Branch("jStr_fwmEi_Q2", &_jStr_fwmEi_Q2);
  // m_tree_AS->Branch("jStr_fwmEi_Q3", &_jStr_fwmEi_Q3);
  // m_tree_AS->Branch("jStr_fwmEi_Q4", &_jStr_fwmEi_Q4);
  // m_tree_AS->Branch("jStr_fwmEi_Q5", &_jStr_fwmEi_Q5);
  // m_tree_AS->Branch("jStr_fwmEi_Q6", &_jStr_fwmEi_Q6);
  // m_tree_AS->Branch("jStr_fwmEi_Q7", &_jStr_fwmEi_Q7);
  // m_tree_AS->Branch("jStr_fwmEi_Q8", &_jStr_fwmEi_Q8);
  // m_tree_AS->Branch("jStr_fwmEi_Pi1", &_jStr_fwmEi_Pi1);
  // m_tree_AS->Branch("jStr_fwmEi_Pi2", &_jStr_fwmEi_Pi2);
  // m_tree_AS->Branch("jStr_fwmEi_Pi3", &_jStr_fwmEi_Pi3);
  // m_tree_AS->Branch("jStr_fwmEi_Pi4", &_jStr_fwmEi_Pi4);
  // m_tree_AS->Branch("jStr_fwmEi_B10", &_jStr_fwmEi_B10);
  // m_tree_AS->Branch("jStr_fwmEi_B20", &_jStr_fwmEi_B20);
  // m_tree_AS->Branch("jStr_fwmEi_B30", &_jStr_fwmEi_B30);
  // m_tree_AS->Branch("jStr_fwmEi_B40", &_jStr_fwmEi_B40);
  // m_tree_AS->Branch("jStr_fwmEi_B50", &_jStr_fwmEi_B50);
  // m_tree_AS->Branch("jStr_fwmEi_B60", &_jStr_fwmEi_B60);
  // m_tree_AS->Branch("jStr_fwmEi_B70", &_jStr_fwmEi_B70);
  // m_tree_AS->Branch("jStr_fwmEi_B80", &_jStr_fwmEi_B80);
  // m_tree_AS->Branch("jStr_fwmEi_C10", &_jStr_fwmEi_C10);
  // m_tree_AS->Branch("jStr_fwmEi_C20", &_jStr_fwmEi_C20);
  // m_tree_AS->Branch("jStr_fwmEi_C30", &_jStr_fwmEi_C30);
  // m_tree_AS->Branch("jStr_fwmEi_C40", &_jStr_fwmEi_C40);
  // m_tree_AS->Branch("jStr_fwmEi_C50", &_jStr_fwmEi_C50);
  // m_tree_AS->Branch("jStr_fwmEi_C60", &_jStr_fwmEi_C60);
  // m_tree_AS->Branch("jStr_fwmEi_C70", &_jStr_fwmEi_C70);
  // m_tree_AS->Branch("jStr_fwmEi_C80", &_jStr_fwmEi_C80);
  // m_tree_AS->Branch("jStr_fwmEi_K10", &_jStr_fwmEi_K10);
  // m_tree_AS->Branch("jStr_fwmEi_K20", &_jStr_fwmEi_K20);
  // m_tree_AS->Branch("jStr_fwmEi_K30", &_jStr_fwmEi_K30);
  // m_tree_AS->Branch("jStr_fwmEi_K40", &_jStr_fwmEi_K40);
  // m_tree_AS->Branch("jStr_fwmEi_K50", &_jStr_fwmEi_K50);
  // m_tree_AS->Branch("jStr_fwmEi_K60", &_jStr_fwmEi_K60);
  // m_tree_AS->Branch("jStr_fwmEi_K70", &_jStr_fwmEi_K70);
  // m_tree_AS->Branch("jStr_fwmEi_K80", &_jStr_fwmEi_K80);
  // m_tree_AS->Branch("jStr_fwmEi_D10", &_jStr_fwmEi_D10);
  // m_tree_AS->Branch("jStr_fwmEi_D20", &_jStr_fwmEi_D20);
  // m_tree_AS->Branch("jStr_fwmEi_D30", &_jStr_fwmEi_D30);
  // m_tree_AS->Branch("jStr_fwmEi_D40", &_jStr_fwmEi_D40);
  // m_tree_AS->Branch("jStr_fwmEi_D50", &_jStr_fwmEi_D50);
  // m_tree_AS->Branch("jStr_fwmEi_D60", &_jStr_fwmEi_D60);
  // m_tree_AS->Branch("jStr_fwmEi_D70", &_jStr_fwmEi_D70);
  // m_tree_AS->Branch("jStr_fwmEi_D80", &_jStr_fwmEi_D80);
  // m_tree_AS->Branch("jStr_fwmEi_H10", &_jStr_fwmEi_H10);
  // m_tree_AS->Branch("jStr_fwmEi_H20", &_jStr_fwmEi_H20);
  // m_tree_AS->Branch("jStr_fwmEi_H30", &_jStr_fwmEi_H30);
  // m_tree_AS->Branch("jStr_fwmEi_H40", &_jStr_fwmEi_H40);
  // m_tree_AS->Branch("jStr_fwmEi_H50", &_jStr_fwmEi_H50);
  // m_tree_AS->Branch("jStr_fwmEi_H60", &_jStr_fwmEi_H60);
  // m_tree_AS->Branch("jStr_fwmEi_H70", &_jStr_fwmEi_H70);
  // m_tree_AS->Branch("jStr_fwmEi_H80", &_jStr_fwmEi_H80);
  // m_tree_AS->Branch("jStr_fwmEi_Q10", &_jStr_fwmEi_Q10);
  // m_tree_AS->Branch("jStr_fwmEi_Q20", &_jStr_fwmEi_Q20);
  // m_tree_AS->Branch("jStr_fwmEi_Q30", &_jStr_fwmEi_Q30);
  // m_tree_AS->Branch("jStr_fwmEi_Q40", &_jStr_fwmEi_Q40);
  // m_tree_AS->Branch("jStr_fwmEi_Q50", &_jStr_fwmEi_Q50);
  // m_tree_AS->Branch("jStr_fwmEi_Q60", &_jStr_fwmEi_Q60);
  // m_tree_AS->Branch("jStr_fwmEi_Q70", &_jStr_fwmEi_Q70);
  // m_tree_AS->Branch("jStr_fwmEi_Q80", &_jStr_fwmEi_Q80);
  m_tree_AS->Branch("jStr_fwmJj_xiPlus", &_jStr_fwmJj_xiPlus);
  m_tree_AS->Branch("jStr_fwmJj_xiMinus", &_jStr_fwmJj_xiMinus);
  m_tree_AS->Branch("jStr_fwmJj_xPlus", &_jStr_fwmJj_xPlus);
  m_tree_AS->Branch("jStr_fwmJj_xMinus", &_jStr_fwmJj_xMinus);
  m_tree_AS->Branch("jStr_fwmJj_Psi1", &_jStr_fwmJj_Psi1);
  m_tree_AS->Branch("jStr_fwmJj_B0", &_jStr_fwmJj_B0);
  m_tree_AS->Branch("jStr_fwmJj_B1", &_jStr_fwmJj_B1);
  m_tree_AS->Branch("jStr_fwmJj_B2", &_jStr_fwmJj_B2);
  m_tree_AS->Branch("jStr_fwmJj_B3", &_jStr_fwmJj_B3);
  m_tree_AS->Branch("jStr_fwmJj_B4", &_jStr_fwmJj_B4);
  m_tree_AS->Branch("jStr_fwmJj_B5", &_jStr_fwmJj_B5);
  m_tree_AS->Branch("jStr_fwmJj_B6", &_jStr_fwmJj_B6);
  m_tree_AS->Branch("jStr_fwmJj_B7", &_jStr_fwmJj_B7);
  m_tree_AS->Branch("jStr_fwmJj_B8", &_jStr_fwmJj_B8);
  m_tree_AS->Branch("jStr_fwmJj_C0", &_jStr_fwmJj_C0);
  m_tree_AS->Branch("jStr_fwmJj_C1", &_jStr_fwmJj_C1);
  m_tree_AS->Branch("jStr_fwmJj_C2", &_jStr_fwmJj_C2);
  m_tree_AS->Branch("jStr_fwmJj_C3", &_jStr_fwmJj_C3);
  m_tree_AS->Branch("jStr_fwmJj_C4", &_jStr_fwmJj_C4);
  m_tree_AS->Branch("jStr_fwmJj_C5", &_jStr_fwmJj_C5);
  m_tree_AS->Branch("jStr_fwmJj_C6", &_jStr_fwmJj_C6);
  m_tree_AS->Branch("jStr_fwmJj_C7", &_jStr_fwmJj_C7);
  m_tree_AS->Branch("jStr_fwmJj_C8", &_jStr_fwmJj_C8);
  m_tree_AS->Branch("jStr_fwmJj_K0", &_jStr_fwmJj_K0);
  m_tree_AS->Branch("jStr_fwmJj_K1", &_jStr_fwmJj_K1);
  m_tree_AS->Branch("jStr_fwmJj_K2", &_jStr_fwmJj_K2);
  m_tree_AS->Branch("jStr_fwmJj_K3", &_jStr_fwmJj_K3);
  m_tree_AS->Branch("jStr_fwmJj_K4", &_jStr_fwmJj_K4);
  m_tree_AS->Branch("jStr_fwmJj_K5", &_jStr_fwmJj_K5);
  m_tree_AS->Branch("jStr_fwmJj_K6", &_jStr_fwmJj_K6);
  m_tree_AS->Branch("jStr_fwmJj_K7", &_jStr_fwmJj_K7);
  m_tree_AS->Branch("jStr_fwmJj_K8", &_jStr_fwmJj_K8);
  m_tree_AS->Branch("jStr_fwmJj_D0", &_jStr_fwmJj_D0);
  m_tree_AS->Branch("jStr_fwmJj_D1", &_jStr_fwmJj_D1);
  m_tree_AS->Branch("jStr_fwmJj_D2", &_jStr_fwmJj_D2);
  m_tree_AS->Branch("jStr_fwmJj_D3", &_jStr_fwmJj_D3);
  m_tree_AS->Branch("jStr_fwmJj_D4", &_jStr_fwmJj_D4);
  m_tree_AS->Branch("jStr_fwmJj_D5", &_jStr_fwmJj_D5);
  m_tree_AS->Branch("jStr_fwmJj_D6", &_jStr_fwmJj_D6);
  m_tree_AS->Branch("jStr_fwmJj_D7", &_jStr_fwmJj_D7);
  m_tree_AS->Branch("jStr_fwmJj_D8", &_jStr_fwmJj_D8);
  m_tree_AS->Branch("jStr_fwmJj_H0", &_jStr_fwmJj_H0);
  m_tree_AS->Branch("jStr_fwmJj_H1", &_jStr_fwmJj_H1);
  m_tree_AS->Branch("jStr_fwmJj_H2", &_jStr_fwmJj_H2);
  m_tree_AS->Branch("jStr_fwmJj_H3", &_jStr_fwmJj_H3);
  m_tree_AS->Branch("jStr_fwmJj_H4", &_jStr_fwmJj_H4);
  m_tree_AS->Branch("jStr_fwmJj_H5", &_jStr_fwmJj_H5);
  m_tree_AS->Branch("jStr_fwmJj_H6", &_jStr_fwmJj_H6);
  m_tree_AS->Branch("jStr_fwmJj_H7", &_jStr_fwmJj_H7);
  m_tree_AS->Branch("jStr_fwmJj_H8", &_jStr_fwmJj_H8);
  m_tree_AS->Branch("jStr_fwmJj_Q0", &_jStr_fwmJj_Q0);
  m_tree_AS->Branch("jStr_fwmJj_Q1", &_jStr_fwmJj_Q1);
  m_tree_AS->Branch("jStr_fwmJj_Q2", &_jStr_fwmJj_Q2);
  m_tree_AS->Branch("jStr_fwmJj_Q3", &_jStr_fwmJj_Q3);
  m_tree_AS->Branch("jStr_fwmJj_Q4", &_jStr_fwmJj_Q4);
  m_tree_AS->Branch("jStr_fwmJj_Q5", &_jStr_fwmJj_Q5);
  m_tree_AS->Branch("jStr_fwmJj_Q6", &_jStr_fwmJj_Q6);
  m_tree_AS->Branch("jStr_fwmJj_Q7", &_jStr_fwmJj_Q7);
  m_tree_AS->Branch("jStr_fwmJj_Q8", &_jStr_fwmJj_Q8);
  m_tree_AS->Branch("jStr_fwmJj_Pi1", &_jStr_fwmJj_Pi1);
  m_tree_AS->Branch("jStr_fwmJj_Pi2", &_jStr_fwmJj_Pi2);
  m_tree_AS->Branch("jStr_fwmJj_Pi3", &_jStr_fwmJj_Pi3);
  m_tree_AS->Branch("jStr_fwmJj_Pi4", &_jStr_fwmJj_Pi4);
  m_tree_AS->Branch("jStr_fwmJj_B10", &_jStr_fwmJj_B10);
  m_tree_AS->Branch("jStr_fwmJj_B20", &_jStr_fwmJj_B20);
  m_tree_AS->Branch("jStr_fwmJj_B30", &_jStr_fwmJj_B30);
  m_tree_AS->Branch("jStr_fwmJj_B40", &_jStr_fwmJj_B40);
  m_tree_AS->Branch("jStr_fwmJj_B50", &_jStr_fwmJj_B50);
  m_tree_AS->Branch("jStr_fwmJj_B60", &_jStr_fwmJj_B60);
  m_tree_AS->Branch("jStr_fwmJj_B70", &_jStr_fwmJj_B70);
  m_tree_AS->Branch("jStr_fwmJj_B80", &_jStr_fwmJj_B80);
  m_tree_AS->Branch("jStr_fwmJj_C10", &_jStr_fwmJj_C10);
  m_tree_AS->Branch("jStr_fwmJj_C20", &_jStr_fwmJj_C20);
  m_tree_AS->Branch("jStr_fwmJj_C30", &_jStr_fwmJj_C30);
  m_tree_AS->Branch("jStr_fwmJj_C40", &_jStr_fwmJj_C40);
  m_tree_AS->Branch("jStr_fwmJj_C50", &_jStr_fwmJj_C50);
  m_tree_AS->Branch("jStr_fwmJj_C60", &_jStr_fwmJj_C60);
  m_tree_AS->Branch("jStr_fwmJj_C70", &_jStr_fwmJj_C70);
  m_tree_AS->Branch("jStr_fwmJj_C80", &_jStr_fwmJj_C80);
  m_tree_AS->Branch("jStr_fwmJj_K10", &_jStr_fwmJj_K10);
  m_tree_AS->Branch("jStr_fwmJj_K20", &_jStr_fwmJj_K20);
  m_tree_AS->Branch("jStr_fwmJj_K30", &_jStr_fwmJj_K30);
  m_tree_AS->Branch("jStr_fwmJj_K40", &_jStr_fwmJj_K40);
  m_tree_AS->Branch("jStr_fwmJj_K50", &_jStr_fwmJj_K50);
  m_tree_AS->Branch("jStr_fwmJj_K60", &_jStr_fwmJj_K60);
  m_tree_AS->Branch("jStr_fwmJj_K70", &_jStr_fwmJj_K70);
  m_tree_AS->Branch("jStr_fwmJj_K80", &_jStr_fwmJj_K80);
  m_tree_AS->Branch("jStr_fwmJj_D10", &_jStr_fwmJj_D10);
  m_tree_AS->Branch("jStr_fwmJj_D20", &_jStr_fwmJj_D20);
  m_tree_AS->Branch("jStr_fwmJj_D30", &_jStr_fwmJj_D30);
  m_tree_AS->Branch("jStr_fwmJj_D40", &_jStr_fwmJj_D40);
  m_tree_AS->Branch("jStr_fwmJj_D50", &_jStr_fwmJj_D50);
  m_tree_AS->Branch("jStr_fwmJj_D60", &_jStr_fwmJj_D60);
  m_tree_AS->Branch("jStr_fwmJj_D70", &_jStr_fwmJj_D70);
  m_tree_AS->Branch("jStr_fwmJj_D80", &_jStr_fwmJj_D80);
  m_tree_AS->Branch("jStr_fwmJj_H10", &_jStr_fwmJj_H10);
  m_tree_AS->Branch("jStr_fwmJj_H20", &_jStr_fwmJj_H20);
  m_tree_AS->Branch("jStr_fwmJj_H30", &_jStr_fwmJj_H30);
  m_tree_AS->Branch("jStr_fwmJj_H40", &_jStr_fwmJj_H40);
  m_tree_AS->Branch("jStr_fwmJj_H50", &_jStr_fwmJj_H50);
  m_tree_AS->Branch("jStr_fwmJj_H60", &_jStr_fwmJj_H60);
  m_tree_AS->Branch("jStr_fwmJj_H70", &_jStr_fwmJj_H70);
  m_tree_AS->Branch("jStr_fwmJj_H80", &_jStr_fwmJj_H80);
  m_tree_AS->Branch("jStr_fwmJj_Q10", &_jStr_fwmJj_Q10);
  m_tree_AS->Branch("jStr_fwmJj_Q20", &_jStr_fwmJj_Q20);
  m_tree_AS->Branch("jStr_fwmJj_Q30", &_jStr_fwmJj_Q30);
  m_tree_AS->Branch("jStr_fwmJj_Q40", &_jStr_fwmJj_Q40);
  m_tree_AS->Branch("jStr_fwmJj_Q50", &_jStr_fwmJj_Q50);
  m_tree_AS->Branch("jStr_fwmJj_Q60", &_jStr_fwmJj_Q60);
  m_tree_AS->Branch("jStr_fwmJj_Q70", &_jStr_fwmJj_Q70);
  m_tree_AS->Branch("jStr_fwmJj_Q80", &_jStr_fwmJj_Q80);
  m_tree_AS->Branch("jStr_fwmJi_xiPlus", &_jStr_fwmJi_xiPlus);
  m_tree_AS->Branch("jStr_fwmJi_xiMinus", &_jStr_fwmJi_xiMinus);
  m_tree_AS->Branch("jStr_fwmJi_xPlus", &_jStr_fwmJi_xPlus);
  m_tree_AS->Branch("jStr_fwmJi_xMinus", &_jStr_fwmJi_xMinus);
  m_tree_AS->Branch("jStr_fwmJi_Psi1", &_jStr_fwmJi_Psi1);
  m_tree_AS->Branch("jStr_fwmJi_B0", &_jStr_fwmJi_B0);
  m_tree_AS->Branch("jStr_fwmJi_B1", &_jStr_fwmJi_B1);
  m_tree_AS->Branch("jStr_fwmJi_B2", &_jStr_fwmJi_B2);
  m_tree_AS->Branch("jStr_fwmJi_B3", &_jStr_fwmJi_B3);
  m_tree_AS->Branch("jStr_fwmJi_B4", &_jStr_fwmJi_B4);
  m_tree_AS->Branch("jStr_fwmJi_B5", &_jStr_fwmJi_B5);
  m_tree_AS->Branch("jStr_fwmJi_B6", &_jStr_fwmJi_B6);
  m_tree_AS->Branch("jStr_fwmJi_B7", &_jStr_fwmJi_B7);
  m_tree_AS->Branch("jStr_fwmJi_B8", &_jStr_fwmJi_B8);
  m_tree_AS->Branch("jStr_fwmJi_C0", &_jStr_fwmJi_C0);
  m_tree_AS->Branch("jStr_fwmJi_C1", &_jStr_fwmJi_C1);
  m_tree_AS->Branch("jStr_fwmJi_C2", &_jStr_fwmJi_C2);
  m_tree_AS->Branch("jStr_fwmJi_C3", &_jStr_fwmJi_C3);
  m_tree_AS->Branch("jStr_fwmJi_C4", &_jStr_fwmJi_C4);
  m_tree_AS->Branch("jStr_fwmJi_C5", &_jStr_fwmJi_C5);
  m_tree_AS->Branch("jStr_fwmJi_C6", &_jStr_fwmJi_C6);
  m_tree_AS->Branch("jStr_fwmJi_C7", &_jStr_fwmJi_C7);
  m_tree_AS->Branch("jStr_fwmJi_C8", &_jStr_fwmJi_C8);
  m_tree_AS->Branch("jStr_fwmJi_K0", &_jStr_fwmJi_K0);
  m_tree_AS->Branch("jStr_fwmJi_K1", &_jStr_fwmJi_K1);
  m_tree_AS->Branch("jStr_fwmJi_K2", &_jStr_fwmJi_K2);
  m_tree_AS->Branch("jStr_fwmJi_K3", &_jStr_fwmJi_K3);
  m_tree_AS->Branch("jStr_fwmJi_K4", &_jStr_fwmJi_K4);
  m_tree_AS->Branch("jStr_fwmJi_K5", &_jStr_fwmJi_K5);
  m_tree_AS->Branch("jStr_fwmJi_K6", &_jStr_fwmJi_K6);
  m_tree_AS->Branch("jStr_fwmJi_K7", &_jStr_fwmJi_K7);
  m_tree_AS->Branch("jStr_fwmJi_K8", &_jStr_fwmJi_K8);
  m_tree_AS->Branch("jStr_fwmJi_D0", &_jStr_fwmJi_D0);
  m_tree_AS->Branch("jStr_fwmJi_D1", &_jStr_fwmJi_D1);
  m_tree_AS->Branch("jStr_fwmJi_D2", &_jStr_fwmJi_D2);
  m_tree_AS->Branch("jStr_fwmJi_D3", &_jStr_fwmJi_D3);
  m_tree_AS->Branch("jStr_fwmJi_D4", &_jStr_fwmJi_D4);
  m_tree_AS->Branch("jStr_fwmJi_D5", &_jStr_fwmJi_D5);
  m_tree_AS->Branch("jStr_fwmJi_D6", &_jStr_fwmJi_D6);
  m_tree_AS->Branch("jStr_fwmJi_D7", &_jStr_fwmJi_D7);
  m_tree_AS->Branch("jStr_fwmJi_D8", &_jStr_fwmJi_D8);
  m_tree_AS->Branch("jStr_fwmJi_H0", &_jStr_fwmJi_H0);
  m_tree_AS->Branch("jStr_fwmJi_H1", &_jStr_fwmJi_H1);
  m_tree_AS->Branch("jStr_fwmJi_H2", &_jStr_fwmJi_H2);
  m_tree_AS->Branch("jStr_fwmJi_H3", &_jStr_fwmJi_H3);
  m_tree_AS->Branch("jStr_fwmJi_H4", &_jStr_fwmJi_H4);
  m_tree_AS->Branch("jStr_fwmJi_H5", &_jStr_fwmJi_H5);
  m_tree_AS->Branch("jStr_fwmJi_H6", &_jStr_fwmJi_H6);
  m_tree_AS->Branch("jStr_fwmJi_H7", &_jStr_fwmJi_H7);
  m_tree_AS->Branch("jStr_fwmJi_H8", &_jStr_fwmJi_H8);
  m_tree_AS->Branch("jStr_fwmJi_Q0", &_jStr_fwmJi_Q0);
  m_tree_AS->Branch("jStr_fwmJi_Q1", &_jStr_fwmJi_Q1);
  m_tree_AS->Branch("jStr_fwmJi_Q2", &_jStr_fwmJi_Q2);
  m_tree_AS->Branch("jStr_fwmJi_Q3", &_jStr_fwmJi_Q3);
  m_tree_AS->Branch("jStr_fwmJi_Q4", &_jStr_fwmJi_Q4);
  m_tree_AS->Branch("jStr_fwmJi_Q5", &_jStr_fwmJi_Q5);
  m_tree_AS->Branch("jStr_fwmJi_Q6", &_jStr_fwmJi_Q6);
  m_tree_AS->Branch("jStr_fwmJi_Q7", &_jStr_fwmJi_Q7);
  m_tree_AS->Branch("jStr_fwmJi_Q8", &_jStr_fwmJi_Q8);
  m_tree_AS->Branch("jStr_fwmJi_Pi1", &_jStr_fwmJi_Pi1);
  m_tree_AS->Branch("jStr_fwmJi_Pi2", &_jStr_fwmJi_Pi2);
  m_tree_AS->Branch("jStr_fwmJi_Pi3", &_jStr_fwmJi_Pi3);
  m_tree_AS->Branch("jStr_fwmJi_Pi4", &_jStr_fwmJi_Pi4);
  m_tree_AS->Branch("jStr_fwmJi_B10", &_jStr_fwmJi_B10);
  m_tree_AS->Branch("jStr_fwmJi_B20", &_jStr_fwmJi_B20);
  m_tree_AS->Branch("jStr_fwmJi_B30", &_jStr_fwmJi_B30);
  m_tree_AS->Branch("jStr_fwmJi_B40", &_jStr_fwmJi_B40);
  m_tree_AS->Branch("jStr_fwmJi_B50", &_jStr_fwmJi_B50);
  m_tree_AS->Branch("jStr_fwmJi_B60", &_jStr_fwmJi_B60);
  m_tree_AS->Branch("jStr_fwmJi_B70", &_jStr_fwmJi_B70);
  m_tree_AS->Branch("jStr_fwmJi_B80", &_jStr_fwmJi_B80);
  m_tree_AS->Branch("jStr_fwmJi_C10", &_jStr_fwmJi_C10);
  m_tree_AS->Branch("jStr_fwmJi_C20", &_jStr_fwmJi_C20);
  m_tree_AS->Branch("jStr_fwmJi_C30", &_jStr_fwmJi_C30);
  m_tree_AS->Branch("jStr_fwmJi_C40", &_jStr_fwmJi_C40);
  m_tree_AS->Branch("jStr_fwmJi_C50", &_jStr_fwmJi_C50);
  m_tree_AS->Branch("jStr_fwmJi_C60", &_jStr_fwmJi_C60);
  m_tree_AS->Branch("jStr_fwmJi_C70", &_jStr_fwmJi_C70);
  m_tree_AS->Branch("jStr_fwmJi_C80", &_jStr_fwmJi_C80);
  m_tree_AS->Branch("jStr_fwmJi_K10", &_jStr_fwmJi_K10);
  m_tree_AS->Branch("jStr_fwmJi_K20", &_jStr_fwmJi_K20);
  m_tree_AS->Branch("jStr_fwmJi_K30", &_jStr_fwmJi_K30);
  m_tree_AS->Branch("jStr_fwmJi_K40", &_jStr_fwmJi_K40);
  m_tree_AS->Branch("jStr_fwmJi_K50", &_jStr_fwmJi_K50);
  m_tree_AS->Branch("jStr_fwmJi_K60", &_jStr_fwmJi_K60);
  m_tree_AS->Branch("jStr_fwmJi_K70", &_jStr_fwmJi_K70);
  m_tree_AS->Branch("jStr_fwmJi_K80", &_jStr_fwmJi_K80);
  m_tree_AS->Branch("jStr_fwmJi_D10", &_jStr_fwmJi_D10);
  m_tree_AS->Branch("jStr_fwmJi_D20", &_jStr_fwmJi_D20);
  m_tree_AS->Branch("jStr_fwmJi_D30", &_jStr_fwmJi_D30);
  m_tree_AS->Branch("jStr_fwmJi_D40", &_jStr_fwmJi_D40);
  m_tree_AS->Branch("jStr_fwmJi_D50", &_jStr_fwmJi_D50);
  m_tree_AS->Branch("jStr_fwmJi_D60", &_jStr_fwmJi_D60);
  m_tree_AS->Branch("jStr_fwmJi_D70", &_jStr_fwmJi_D70);
  m_tree_AS->Branch("jStr_fwmJi_D80", &_jStr_fwmJi_D80);
  m_tree_AS->Branch("jStr_fwmJi_H10", &_jStr_fwmJi_H10);
  m_tree_AS->Branch("jStr_fwmJi_H20", &_jStr_fwmJi_H20);
  m_tree_AS->Branch("jStr_fwmJi_H30", &_jStr_fwmJi_H30);
  m_tree_AS->Branch("jStr_fwmJi_H40", &_jStr_fwmJi_H40);
  m_tree_AS->Branch("jStr_fwmJi_H50", &_jStr_fwmJi_H50);
  m_tree_AS->Branch("jStr_fwmJi_H60", &_jStr_fwmJi_H60);
  m_tree_AS->Branch("jStr_fwmJi_H70", &_jStr_fwmJi_H70);
  m_tree_AS->Branch("jStr_fwmJi_H80", &_jStr_fwmJi_H80);
  m_tree_AS->Branch("jStr_fwmJi_Q10", &_jStr_fwmJi_Q10);
  m_tree_AS->Branch("jStr_fwmJi_Q20", &_jStr_fwmJi_Q20);
  m_tree_AS->Branch("jStr_fwmJi_Q30", &_jStr_fwmJi_Q30);
  m_tree_AS->Branch("jStr_fwmJi_Q40", &_jStr_fwmJi_Q40);
  m_tree_AS->Branch("jStr_fwmJi_Q50", &_jStr_fwmJi_Q50);
  m_tree_AS->Branch("jStr_fwmJi_Q60", &_jStr_fwmJi_Q60);
  m_tree_AS->Branch("jStr_fwmJi_Q70", &_jStr_fwmJi_Q70);
  m_tree_AS->Branch("jStr_fwmJi_Q80", &_jStr_fwmJi_Q80);
  m_tree_AS->Branch("jStr_authorE", &_jStr_authorE);
  m_tree_AS->Branch("jStr_authorJ", &_jStr_authorJ);
  m_tree_AS->Branch("jStr_radialParamE", &_jStr_radialParamE);
  m_tree_AS->Branch("jStr_radialParamJ", &_jStr_radialParamJ);
  m_tree_AS->Branch("jStr_algE", &_jStr_algE);
  m_tree_AS->Branch("jStr_algJ", &_jStr_algJ);
  m_tree_AS->Branch("jStr_inputE", &_jStr_inputE);
  m_tree_AS->Branch("jStr_inputJ", &_jStr_inputJ);
  m_tree_AS->Branch("jStr_bjetE", &_jStr_bjetE);
  m_tree_AS->Branch("jStr_bjetJ", &_jStr_bjetJ);
  m_tree_AS->Branch("jStr_xjetE", &_jStr_xjetE);
  m_tree_AS->Branch("jStr_xjetJ", &_jStr_xjetJ);
  m_tree_AS->Branch("jStr_myJetsE", &_jStr_myJetsE);
  m_tree_AS->Branch("jStr_myJetsJ", &_jStr_myJetsJ);
  m_tree_AS->Branch("jStr_indexJ", &_jStr_indexJ);
  m_tree_AS->Branch("jStr_jetMiscJ_numConstituents", &_jStr_jetMiscJ_numConstituents);
  m_tree_AS->Branch("jStr_jetMiscJ_jetM", &_jStr_jetMiscJ_jetM);
  m_tree_AS->Branch("jStr_jetMiscJ_jetMt", &_jStr_jetMiscJ_jetMt);
  m_tree_AS->Branch("jStr_jetMiscJ_jetE", &_jStr_jetMiscJ_jetE);
  m_tree_AS->Branch("jStr_jetMiscJ_jetP", &_jStr_jetMiscJ_jetP);
  m_tree_AS->Branch("jStr_jetMiscJ_jetEt", &_jStr_jetMiscJ_jetEt);
  m_tree_AS->Branch("jStr_jetMiscJ_jetPt", &_jStr_jetMiscJ_jetPt);
  m_tree_AS->Branch("jStr_jetMiscJ_jetPhi", &_jStr_jetMiscJ_jetPhi);
  m_tree_AS->Branch("jStr_jetMiscJ_jetEta", &_jStr_jetMiscJ_jetEta);
  m_tree_AS->Branch("jStr_jetMiscJ_jetRapidity", &_jStr_jetMiscJ_jetRapidity);
  m_tree_AS->Branch("jStr_jetMiscJ_xJ", &_jStr_jetMiscJ_xJ);
  m_tree_AS->Branch("jStr_jetMiscJ_gamma", &_jStr_jetMiscJ_gamma);
  m_tree_AS->Branch("jStr_jetMiscJ_R", &_jStr_jetMiscJ_R);
  m_tree_AS->Branch("jStr_jetMiscJ_Cbar", &_jStr_jetMiscJ_Cbar);
  m_tree_AS->Branch("jStr_jetMiscJ_numSubjets", &_jStr_jetMiscJ_numSubjets);
  m_tree_AS->Branch("jStr_pullJ_det", &_jStr_pullJ_det);
  m_tree_AS->Branch("jStr_pullJ_ratio", &_jStr_pullJ_ratio);
  m_tree_AS->Branch("jStr_pullJ_pullPf", &_jStr_pullJ_pullPf);
  m_tree_AS->Branch("jStr_pullJ_angularEccentricity", &_jStr_pullJ_angularEccentricity);
  m_tree_AS->Branch("jStr_pullJ_orientation", &_jStr_pullJ_orientation);
  m_tree_AS->Branch("jStr_pullJ_girth", &_jStr_pullJ_girth);
  m_tree_AS->Branch("jStr_pullJ_Cbar", &_jStr_pullJ_Cbar);
  m_tree_AS->Branch("jStr_pullJ_g", &_jStr_pullJ_g);
  m_tree_AS->Branch("jStr_pullJ_e", &_jStr_pullJ_e);
  m_tree_AS->Branch("jStr_pullJ_B", &_jStr_pullJ_B);
  m_tree_AS->Branch("jStr_pullJ_logB", &_jStr_pullJ_logB);
  m_tree_AS->Branch("jStr_pullJ_pullTheta", &_jStr_pullJ_pullTheta);
  m_tree_AS->Branch("jStr_pullJ_pullMag", &_jStr_pullJ_pullMag);
  m_tree_AS->Branch("jStr_ubhJ_z", &_jStr_ubhJ_z);
  m_tree_AS->Branch("jStr_ubhJ_z2", &_jStr_ubhJ_z2);
  m_tree_AS->Branch("jStr_ubhJ_a1", &_jStr_ubhJ_a1);
  m_tree_AS->Branch("jStr_ubhJ_a2", &_jStr_ubhJ_a2);
  m_tree_AS->Branch("jStr_ubhJ_a3", &_jStr_ubhJ_a3);
  m_tree_AS->Branch("jStr_ubhJ_meanpt", &_jStr_ubhJ_meanpt);
  m_tree_AS->Branch("jStr_ubhJ_meanet", &_jStr_ubhJ_meanet);
  m_tree_AS->Branch("jStr_ubhJ_mbar", &_jStr_ubhJ_mbar);
  m_tree_AS->Branch("jStr_ubhJ_massDemocracy", &_jStr_ubhJ_massDemocracy);
  m_tree_AS->Branch("jStr_ubhJ_fE1", &_jStr_ubhJ_fE1);
  m_tree_AS->Branch("jStr_ubhJ_fE2", &_jStr_ubhJ_fE2);
  m_tree_AS->Branch("jStr_ubhJ_fE3", &_jStr_ubhJ_fE3);
  m_tree_AS->Branch("jStr_ubhJ_fET1", &_jStr_ubhJ_fET1);
  m_tree_AS->Branch("jStr_ubhJ_fET2", &_jStr_ubhJ_fET2);
  m_tree_AS->Branch("jStr_ubhJ_fET3", &_jStr_ubhJ_fET3);
  m_tree_AS->Branch("jStr_ubhJ_Alpha", &_jStr_ubhJ_Alpha);
  m_tree_AS->Branch("jStr_ubhJ_AlphaT", &_jStr_ubhJ_AlphaT);
  m_tree_AS->Branch("jStr_ubhJ_betaflow", &_jStr_ubhJ_betaflow);
  m_tree_AS->Branch("jStr_ubhJ_betaflow_1GeV", &_jStr_ubhJ_betaflow_1GeV);
  m_tree_AS->Branch("jStr_ubhJ_betaflow_5GeV", &_jStr_ubhJ_betaflow_5GeV);
  m_tree_AS->Branch("jStr_ubhJ_y23", &_jStr_ubhJ_y23);
  m_tree_AS->Branch("jStr_ubhJ_y23_1GeV", &_jStr_ubhJ_y23_1GeV);
  m_tree_AS->Branch("jStr_ubhJ_y23_5GeV", &_jStr_ubhJ_y23_5GeV);
  m_tree_AS->Branch("jStr_ubhJ_lny23", &_jStr_ubhJ_lny23);
  m_tree_AS->Branch("jStr_ubhJ_lny23_1GeV", &_jStr_ubhJ_lny23_1GeV);
  m_tree_AS->Branch("jStr_ubhJ_lny23_5GeV", &_jStr_ubhJ_lny23_5GeV);
  m_tree_AS->Branch("jStr_ubhJ_subjetAsymmetry", &_jStr_ubhJ_subjetAsymmetry);
  m_tree_AS->Branch("jStr_dipJ_dipolarity", &_jStr_dipJ_dipolarity);
  m_tree_AS->Branch("jStr_nsubjnessJ_tau1", &_jStr_nsubjnessJ_tau1);
  m_tree_AS->Branch("jStr_nsubjnessJ_tau2", &_jStr_nsubjnessJ_tau2);
  m_tree_AS->Branch("jStr_nsubjnessJ_tau3", &_jStr_nsubjnessJ_tau3);
  m_tree_AS->Branch("jStr_nsubjnessJ_tau2tau1", &_jStr_nsubjnessJ_tau2tau1);
  m_tree_AS->Branch("jStr_nsubjnessJ_tau3tau2", &_jStr_nsubjnessJ_tau3tau2);
  m_tree_AS->Branch("jStr_psiJ2_psi", &_jStr_psiJ2_psi);
  m_tree_AS->Branch("jStr_psiJ2_rho", &_jStr_psiJ2_rho);
  m_tree_AS->Branch("jStr_psiRatiosJ_psi1", &_jStr_psiRatiosJ_psi1);
  m_tree_AS->Branch("jStr_psiRatiosJ_psi2", &_jStr_psiRatiosJ_psi2);
  m_tree_AS->Branch("jStr_psiRatiosJ_psi3", &_jStr_psiRatiosJ_psi3);
  m_tree_AS->Branch("jStr_psiRatiosJ_psi7", &_jStr_psiRatiosJ_psi7);
  m_tree_AS->Branch("jStr_psiRatiosJ_psi717", &_jStr_psiRatiosJ_psi717);
  m_tree_AS->Branch("jStr_psiRatiosJ_psi127", &_jStr_psiRatiosJ_psi127);
  m_tree_AS->Branch("jStr_psiRatiosJ_psi37", &_jStr_psiRatiosJ_psi37);
  m_tree_AS->Branch("jStr_bcfJ_bcfVersion_v0a1", &_jStr_bcfJ_bcfVersion_v0a1);
  m_tree_AS->Branch("jStr_bcfJ_a_v0a1", &_jStr_bcfJ_a_v0a1);
  m_tree_AS->Branch("jStr_bcfJ_bcfT_v0a1", &_jStr_bcfJ_bcfT_v0a1);
  m_tree_AS->Branch("jStr_bcfJ_bcf_v0a1", &_jStr_bcfJ_bcf_v0a1);
  m_tree_AS->Branch("jStr_bcfJ_bcfAsymY_v0a1", &_jStr_bcfJ_bcfAsymY_v0a1);
  m_tree_AS->Branch("jStr_bcfJ_bcfAsymPhi_v0a1", &_jStr_bcfJ_bcfAsymPhi_v0a1);
  m_tree_AS->Branch("jStr_bcfJ_bcfAsymYPhi_v0a1", &_jStr_bcfJ_bcfAsymYPhi_v0a1);
  m_tree_AS->Branch("jStr_bcfJ_bcfAsymYPhi2_v0a1", &_jStr_bcfJ_bcfAsymYPhi2_v0a1);
  m_tree_AS->Branch("jStr_bcfJ_bcfVersion_v0a2", &_jStr_bcfJ_bcfVersion_v0a2);
  m_tree_AS->Branch("jStr_bcfJ_a_v0a2", &_jStr_bcfJ_a_v0a2);
  m_tree_AS->Branch("jStr_bcfJ_bcfT_v0a2", &_jStr_bcfJ_bcfT_v0a2);
  m_tree_AS->Branch("jStr_bcfJ_bcf_v0a2", &_jStr_bcfJ_bcf_v0a2);
  m_tree_AS->Branch("jStr_bcfJ_bcfAsymY_v0a2", &_jStr_bcfJ_bcfAsymY_v0a2);
  m_tree_AS->Branch("jStr_bcfJ_bcfAsymPhi_v0a2", &_jStr_bcfJ_bcfAsymPhi_v0a2);
  m_tree_AS->Branch("jStr_bcfJ_bcfAsymYPhi_v0a2", &_jStr_bcfJ_bcfAsymYPhi_v0a2);
  m_tree_AS->Branch("jStr_bcfJ_bcfAsymYPhi2_v0a2", &_jStr_bcfJ_bcfAsymYPhi2_v0a2);
  m_tree_AS->Branch("jStr_pfJ_pf", &_jStr_pfJ_pf);
  m_tree_AS->Branch("jStr_pfJ_detST", &_jStr_pfJ_detST);
  m_tree_AS->Branch("jStr_pfJ_lambdaST", &_jStr_pfJ_lambdaST);
  m_tree_AS->Branch("jStr_zminJ", &_jStr_zminJ);
  m_tree_AS->Branch("jStr_zmaxJ", &_jStr_zmaxJ);
  m_tree_AS->Branch("jStr_tauJ09_a", &_jStr_tauJ09_a);
  m_tree_AS->Branch("jStr_tauJ09_tau", &_jStr_tauJ09_tau);
  m_tree_AS->Branch("jStr_tauJ20_a", &_jStr_tauJ20_a);
  m_tree_AS->Branch("jStr_tauJ20_tau", &_jStr_tauJ20_tau);
  m_tree_AS->Branch("jStr_tauJ40_a", &_jStr_tauJ40_a);
  m_tree_AS->Branch("jStr_tauJ40_tau", &_jStr_tauJ40_tau);
  m_tree_AS->Branch("jStr_radEi_bcfJ1_v0a1", &_jStr_radEi_bcfJ1_v0a1);
  m_tree_AS->Branch("jStr_radEi_bcfJ2_v0a1", &_jStr_radEi_bcfJ2_v0a1);
  m_tree_AS->Branch("jStr_radEi_bcfJ3_v0a1", &_jStr_radEi_bcfJ3_v0a1);
  m_tree_AS->Branch("jStr_radEi_bcfTJ1_v0a1", &_jStr_radEi_bcfTJ1_v0a1);
  m_tree_AS->Branch("jStr_radEi_bcfTJ2_v0a1", &_jStr_radEi_bcfTJ2_v0a1);
  m_tree_AS->Branch("jStr_radEi_bcfTJ3_v0a1", &_jStr_radEi_bcfTJ3_v0a1);
  m_tree_AS->Branch("jStr_radEi_bcfAsymYJ1_v0a1", &_jStr_radEi_bcfAsymYJ1_v0a1);
  m_tree_AS->Branch("jStr_radEi_bcfAsymYJ2_v0a1", &_jStr_radEi_bcfAsymYJ2_v0a1);
  m_tree_AS->Branch("jStr_radEi_bcfAsymYJ3_v0a1", &_jStr_radEi_bcfAsymYJ3_v0a1);
  m_tree_AS->Branch("jStr_radEi_bcfAsymPhiJ1_v0a1", &_jStr_radEi_bcfAsymPhiJ1_v0a1);
  m_tree_AS->Branch("jStr_radEi_bcfAsymPhiJ2_v0a1", &_jStr_radEi_bcfAsymPhiJ2_v0a1);
  m_tree_AS->Branch("jStr_radEi_bcfAsymPhiJ3_v0a1", &_jStr_radEi_bcfAsymPhiJ3_v0a1);
  m_tree_AS->Branch("jStr_radEi_bcfAsymYPhiJ1_v0a1", &_jStr_radEi_bcfAsymYPhiJ1_v0a1);
  m_tree_AS->Branch("jStr_radEi_bcfAsymYPhiJ2_v0a1", &_jStr_radEi_bcfAsymYPhiJ2_v0a1);
  m_tree_AS->Branch("jStr_radEi_bcfAsymYPhiJ3_v0a1", &_jStr_radEi_bcfAsymYPhiJ3_v0a1);
  m_tree_AS->Branch("jStr_radEi_bcfAsymYPhi2J1_v0a1", &_jStr_radEi_bcfAsymYPhi2J1_v0a1);
  m_tree_AS->Branch("jStr_radEi_bcfAsymYPhi2J2_v0a1", &_jStr_radEi_bcfAsymYPhi2J2_v0a1);
  m_tree_AS->Branch("jStr_radEi_bcfAsymYPhi2J3_v0a1", &_jStr_radEi_bcfAsymYPhi2J3_v0a1);
  m_tree_AS->Branch("jStr_radEi_bcfJ1_v0a2", &_jStr_radEi_bcfJ1_v0a2);
  m_tree_AS->Branch("jStr_radEi_bcfJ2_v0a2", &_jStr_radEi_bcfJ2_v0a2);
  m_tree_AS->Branch("jStr_radEi_bcfJ3_v0a2", &_jStr_radEi_bcfJ3_v0a2);
  m_tree_AS->Branch("jStr_radEi_bcfTJ1_v0a2", &_jStr_radEi_bcfTJ1_v0a2);
  m_tree_AS->Branch("jStr_radEi_bcfTJ2_v0a2", &_jStr_radEi_bcfTJ2_v0a2);
  m_tree_AS->Branch("jStr_radEi_bcfTJ3_v0a2", &_jStr_radEi_bcfTJ3_v0a2);
  m_tree_AS->Branch("jStr_radEi_bcfAsymYJ1_v0a2", &_jStr_radEi_bcfAsymYJ1_v0a2);
  m_tree_AS->Branch("jStr_radEi_bcfAsymYJ2_v0a2", &_jStr_radEi_bcfAsymYJ2_v0a2);
  m_tree_AS->Branch("jStr_radEi_bcfAsymYJ3_v0a2", &_jStr_radEi_bcfAsymYJ3_v0a2);
  m_tree_AS->Branch("jStr_radEi_bcfAsymPhiJ1_v0a2", &_jStr_radEi_bcfAsymPhiJ1_v0a2);
  m_tree_AS->Branch("jStr_radEi_bcfAsymPhiJ2_v0a2", &_jStr_radEi_bcfAsymPhiJ2_v0a2);
  m_tree_AS->Branch("jStr_radEi_bcfAsymPhiJ3_v0a2", &_jStr_radEi_bcfAsymPhiJ3_v0a2);
  m_tree_AS->Branch("jStr_radEi_bcfAsymYPhiJ1_v0a2", &_jStr_radEi_bcfAsymYPhiJ1_v0a2);
  m_tree_AS->Branch("jStr_radEi_bcfAsymYPhiJ2_v0a2", &_jStr_radEi_bcfAsymYPhiJ2_v0a2);
  m_tree_AS->Branch("jStr_radEi_bcfAsymYPhiJ3_v0a2", &_jStr_radEi_bcfAsymYPhiJ3_v0a2);
  m_tree_AS->Branch("jStr_radEi_bcfAsymYPhi2J1_v0a2", &_jStr_radEi_bcfAsymYPhi2J1_v0a2);
  m_tree_AS->Branch("jStr_radEi_bcfAsymYPhi2J2_v0a2", &_jStr_radEi_bcfAsymYPhi2J2_v0a2);
  m_tree_AS->Branch("jStr_radEi_bcfAsymYPhi2J3_v0a2", &_jStr_radEi_bcfAsymYPhi2J3_v0a2);
  m_tree_AS->Branch("jStr_radEi_BJ1", &_jStr_radEi_BJ1);
  m_tree_AS->Branch("jStr_radEi_BJ2", &_jStr_radEi_BJ2);
  m_tree_AS->Branch("jStr_radEi_BJ3", &_jStr_radEi_BJ3);
  m_tree_AS->Branch("jStr_radEi_girthJ1", &_jStr_radEi_girthJ1);
  m_tree_AS->Branch("jStr_radEi_girthJ2", &_jStr_radEi_girthJ2);
  m_tree_AS->Branch("jStr_radEi_girthJ3", &_jStr_radEi_girthJ3);
  m_tree_AS->Branch("jStr_radEi_girth32", &_jStr_radEi_girth32);
  m_tree_AS->Branch("jStr_radEi_girth21", &_jStr_radEi_girth21);
  m_tree_AS->Branch("jStr_radEi_girthAsymJ1J2", &_jStr_radEi_girthAsymJ1J2);
  m_tree_AS->Branch("jStr_radEi_alpha1", &_jStr_radEi_alpha1);
  m_tree_AS->Branch("jStr_radEi_alpha2", &_jStr_radEi_alpha2);
  m_tree_AS->Branch("jStr_radEi_alpha", &_jStr_radEi_alpha);
  m_tree_AS->Branch("jStr_radEi_beta1", &_jStr_radEi_beta1);
  m_tree_AS->Branch("jStr_radEi_beta2", &_jStr_radEi_beta2);
  m_tree_AS->Branch("jStr_radEi_beta", &_jStr_radEi_beta);
  m_tree_AS->Branch("jStr_radEi_thetaJ1J2", &_jStr_radEi_thetaJ1J2);
  m_tree_AS->Branch("jStr_radEi_dipolarityInfLine", &_jStr_radEi_dipolarityInfLine);
  m_tree_AS->Branch("jStr_radEi_dipolarityLineSeg", &_jStr_radEi_dipolarityLineSeg);
  m_tree_AS->Branch("jStr_radEi_BCF1", &_jStr_radEi_BCF1);
  m_tree_AS->Branch("jStr_radEi_BCF2", &_jStr_radEi_BCF2);
  m_tree_AS->Branch("jStr_radEi_BCF3", &_jStr_radEi_BCF3);
  m_tree_AS->Branch("jStr_radEi_dipolarity", &_jStr_radEi_dipolarity);
  m_tree_AS->Branch("jStr_sphEJ_detSphericity", &_jStr_sphEJ_detSphericity);
  m_tree_AS->Branch("jStr_sphEJ_detSpherocity", &_jStr_sphEJ_detSpherocity);
  m_tree_AS->Branch("jStr_sphEJ_sphericityLambda1", &_jStr_sphEJ_sphericityLambda1);
  m_tree_AS->Branch("jStr_sphEJ_sphericityLambda2", &_jStr_sphEJ_sphericityLambda2);
  m_tree_AS->Branch("jStr_sphEJ_sphericityLambda3", &_jStr_sphEJ_sphericityLambda3);
  m_tree_AS->Branch("jStr_sphEJ_spherocityLambda1", &_jStr_sphEJ_spherocityLambda1);
  m_tree_AS->Branch("jStr_sphEJ_spherocityLambda2", &_jStr_sphEJ_spherocityLambda2);
  m_tree_AS->Branch("jStr_sphEJ_spherocityLambda3", &_jStr_sphEJ_spherocityLambda3);
  m_tree_AS->Branch("jStr_sphEJ_circularity", &_jStr_sphEJ_circularity);
  m_tree_AS->Branch("jStr_sphEJ_sphericity", &_jStr_sphEJ_sphericity);
  m_tree_AS->Branch("jStr_sphEJ_spherocity", &_jStr_sphEJ_spherocity);
  m_tree_AS->Branch("jStr_sphEJ_aplanarity", &_jStr_sphEJ_aplanarity);
  m_tree_AS->Branch("jStr_sphEJ_aplanority", &_jStr_sphEJ_aplanority);
  m_tree_AS->Branch("jStr_sphEJ_Y", &_jStr_sphEJ_Y);
  m_tree_AS->Branch("jStr_sphEJ_planarity", &_jStr_sphEJ_planarity);
  m_tree_AS->Branch("jStr_sphEJ_planority", &_jStr_sphEJ_planority);
  m_tree_AS->Branch("jStr_sphEJ_Dshape", &_jStr_sphEJ_Dshape);
  m_tree_AS->Branch("jStr_sphEJ_Cshape", &_jStr_sphEJ_Cshape);
  m_tree_AS->Branch("jStr_sphEJ_H2", &_jStr_sphEJ_H2);
  m_tree_AS->Branch("jStr_sphEJ_Fshape", &_jStr_sphEJ_Fshape);
  m_tree_AS->Branch("jStr_sphEJ_beamThrust", &_jStr_sphEJ_beamThrust);
  m_tree_AS->Branch("jStr_sphEJ_G", &_jStr_sphEJ_G);
  m_tree_AS->Branch("jStr_sphEJ_ST2D", &_jStr_sphEJ_ST2D);
  m_tree_AS->Branch("jStr_sphEJ_detMlin", &_jStr_sphEJ_detMlin);
  m_tree_AS->Branch("jStr_sphEJ_pft", &_jStr_sphEJ_pft);
  m_tree_AS->Branch("jStr_sphEi_detSphericity", &_jStr_sphEi_detSphericity);
  m_tree_AS->Branch("jStr_sphEi_detSpherocity", &_jStr_sphEi_detSpherocity);
  m_tree_AS->Branch("jStr_sphEi_sphericityLambda1", &_jStr_sphEi_sphericityLambda1);
  m_tree_AS->Branch("jStr_sphEi_sphericityLambda2", &_jStr_sphEi_sphericityLambda2);
  m_tree_AS->Branch("jStr_sphEi_sphericityLambda3", &_jStr_sphEi_sphericityLambda3);
  m_tree_AS->Branch("jStr_sphEi_spherocityLambda1", &_jStr_sphEi_spherocityLambda1);
  m_tree_AS->Branch("jStr_sphEi_spherocityLambda2", &_jStr_sphEi_spherocityLambda2);
  m_tree_AS->Branch("jStr_sphEi_spherocityLambda3", &_jStr_sphEi_spherocityLambda3);
  m_tree_AS->Branch("jStr_sphEi_circularity", &_jStr_sphEi_circularity);
  m_tree_AS->Branch("jStr_sphEi_sphericity", &_jStr_sphEi_sphericity);
  m_tree_AS->Branch("jStr_sphEi_spherocity", &_jStr_sphEi_spherocity);
  m_tree_AS->Branch("jStr_sphEi_aplanarity", &_jStr_sphEi_aplanarity);
  m_tree_AS->Branch("jStr_sphEi_aplanority", &_jStr_sphEi_aplanority);
  m_tree_AS->Branch("jStr_sphEi_Y", &_jStr_sphEi_Y);
  m_tree_AS->Branch("jStr_sphEi_planarity", &_jStr_sphEi_planarity);
  m_tree_AS->Branch("jStr_sphEi_planority", &_jStr_sphEi_planority);
  m_tree_AS->Branch("jStr_sphEi_Dshape", &_jStr_sphEi_Dshape);
  m_tree_AS->Branch("jStr_sphEi_Cshape", &_jStr_sphEi_Cshape);
  m_tree_AS->Branch("jStr_sphEi_H2", &_jStr_sphEi_H2);
  m_tree_AS->Branch("jStr_sphEi_Fshape", &_jStr_sphEi_Fshape);
  m_tree_AS->Branch("jStr_sphEi_beamThrust", &_jStr_sphEi_beamThrust);
  m_tree_AS->Branch("jStr_sphEi_G", &_jStr_sphEi_G);
  m_tree_AS->Branch("jStr_sphEi_ST2D", &_jStr_sphEi_ST2D);
  m_tree_AS->Branch("jStr_sphEi_detMlin", &_jStr_sphEi_detMlin);
  m_tree_AS->Branch("jStr_sphEi_pft", &_jStr_sphEi_pft);
  m_tree_AS->Branch("jStr_sphJj_detSphericity", &_jStr_sphJj_detSphericity);
  m_tree_AS->Branch("jStr_sphJj_detSpherocity", &_jStr_sphJj_detSpherocity);
  m_tree_AS->Branch("jStr_sphJj_sphericityLambda1", &_jStr_sphJj_sphericityLambda1);
  m_tree_AS->Branch("jStr_sphJj_sphericityLambda2", &_jStr_sphJj_sphericityLambda2);
  m_tree_AS->Branch("jStr_sphJj_sphericityLambda3", &_jStr_sphJj_sphericityLambda3);
  m_tree_AS->Branch("jStr_sphJj_spherocityLambda1", &_jStr_sphJj_spherocityLambda1);
  m_tree_AS->Branch("jStr_sphJj_spherocityLambda2", &_jStr_sphJj_spherocityLambda2);
  m_tree_AS->Branch("jStr_sphJj_spherocityLambda3", &_jStr_sphJj_spherocityLambda3);
  m_tree_AS->Branch("jStr_sphJj_circularity", &_jStr_sphJj_circularity);
  m_tree_AS->Branch("jStr_sphJj_sphericity", &_jStr_sphJj_sphericity);
  m_tree_AS->Branch("jStr_sphJj_spherocity", &_jStr_sphJj_spherocity);
  m_tree_AS->Branch("jStr_sphJj_aplanarity", &_jStr_sphJj_aplanarity);
  m_tree_AS->Branch("jStr_sphJj_aplanority", &_jStr_sphJj_aplanority);
  m_tree_AS->Branch("jStr_sphJj_Y", &_jStr_sphJj_Y);
  m_tree_AS->Branch("jStr_sphJj_planarity", &_jStr_sphJj_planarity);
  m_tree_AS->Branch("jStr_sphJj_planority", &_jStr_sphJj_planority);
  m_tree_AS->Branch("jStr_sphJj_Dshape", &_jStr_sphJj_Dshape);
  m_tree_AS->Branch("jStr_sphJj_Cshape", &_jStr_sphJj_Cshape);
  m_tree_AS->Branch("jStr_sphJj_H2", &_jStr_sphJj_H2);
  m_tree_AS->Branch("jStr_sphJj_Fshape", &_jStr_sphJj_Fshape);
  m_tree_AS->Branch("jStr_sphJj_beamThrust", &_jStr_sphJj_beamThrust);
  m_tree_AS->Branch("jStr_sphJj_G", &_jStr_sphJj_G);
  m_tree_AS->Branch("jStr_sphJj_ST2D", &_jStr_sphJj_ST2D);
  m_tree_AS->Branch("jStr_sphJj_detMlin", &_jStr_sphJj_detMlin);
  m_tree_AS->Branch("jStr_sphJj_pft", &_jStr_sphJj_pft);
  m_tree_AS->Branch("jStr_sphJi_detSphericity", &_jStr_sphJi_detSphericity);
  m_tree_AS->Branch("jStr_sphJi_detSpherocity", &_jStr_sphJi_detSpherocity);
  m_tree_AS->Branch("jStr_sphJi_sphericityLambda1", &_jStr_sphJi_sphericityLambda1);
  m_tree_AS->Branch("jStr_sphJi_sphericityLambda2", &_jStr_sphJi_sphericityLambda2);
  m_tree_AS->Branch("jStr_sphJi_sphericityLambda3", &_jStr_sphJi_sphericityLambda3);
  m_tree_AS->Branch("jStr_sphJi_spherocityLambda1", &_jStr_sphJi_spherocityLambda1);
  m_tree_AS->Branch("jStr_sphJi_spherocityLambda2", &_jStr_sphJi_spherocityLambda2);
  m_tree_AS->Branch("jStr_sphJi_spherocityLambda3", &_jStr_sphJi_spherocityLambda3);
  m_tree_AS->Branch("jStr_sphJi_circularity", &_jStr_sphJi_circularity);
  m_tree_AS->Branch("jStr_sphJi_sphericity", &_jStr_sphJi_sphericity);
  m_tree_AS->Branch("jStr_sphJi_spherocity", &_jStr_sphJi_spherocity);
  m_tree_AS->Branch("jStr_sphJi_aplanarity", &_jStr_sphJi_aplanarity);
  m_tree_AS->Branch("jStr_sphJi_aplanority", &_jStr_sphJi_aplanority);
  m_tree_AS->Branch("jStr_sphJi_Y", &_jStr_sphJi_Y);
  m_tree_AS->Branch("jStr_sphJi_planarity", &_jStr_sphJi_planarity);
  m_tree_AS->Branch("jStr_sphJi_planority", &_jStr_sphJi_planority);
  m_tree_AS->Branch("jStr_sphJi_Dshape", &_jStr_sphJi_Dshape);
  m_tree_AS->Branch("jStr_sphJi_Cshape", &_jStr_sphJi_Cshape);
  m_tree_AS->Branch("jStr_sphJi_H2", &_jStr_sphJi_H2);
  m_tree_AS->Branch("jStr_sphJi_Fshape", &_jStr_sphJi_Fshape);
  m_tree_AS->Branch("jStr_sphJi_beamThrust", &_jStr_sphJi_beamThrust);
  m_tree_AS->Branch("jStr_sphJi_G", &_jStr_sphJi_G);
  m_tree_AS->Branch("jStr_sphJi_ST2D", &_jStr_sphJi_ST2D);
  m_tree_AS->Branch("jStr_sphJi_detMlin", &_jStr_sphJi_detMlin);
  m_tree_AS->Branch("jStr_sphJi_pft", &_jStr_sphJi_pft);

#ifdef OUTPUT_CONSTITUENTS
  ATH_MSG_WARNING("Writing out jet constituents; this will bloat output ntuple!");
  m_tree_AS->Branch("jStr_nConst", &_nConst);
  m_tree_AS->Branch("jStr_const_et", &_const_et);
  m_tree_AS->Branch("jStr_const_rapidity", &_const_rapidity);
  m_tree_AS->Branch("jStr_const_phi", &_const_phi);
  m_tree_AS->Branch("jStr_const_Ret", &_const_Ret);
  m_tree_AS->Branch("jStr_const_Drapidity", &_const_Drapidity);
  m_tree_AS->Branch("jStr_const_Dphi", &_const_Dphi);
#endif

  m_tree_AS->Branch("YFlip12", &_YFlip12);
  m_tree_AS->Branch("YFlip23", &_YFlip23);
  m_tree_AS->Branch("EM_FRACTION_CLUSTER", &_EM_FRACTION_CLUSTER);
  m_tree_AS->Branch("ELLIPTICAREA", &_ELLIPTICAREA);
  m_tree_AS->Branch("AMAREA", &_AMAREA);
  m_tree_AS->Branch("HULL_LENGTH", &_HULL_LENGTH);
  m_tree_AS->Branch("HULL_AREA", &_HULL_AREA);
  m_tree_AS->Branch("EM_FRACTION_MCTRUTH", &_EM_FRACTION_MCTRUTH);
  m_tree_AS->Branch("LowEtConstituentsFrac", &_LowEtConstituentsFrac);
  m_tree_AS->Branch("JetEccentricity", &_JetEccentricity);
  m_tree_AS->Branch("DRtoReco", &_DRtoReco);
  m_tree_AS->Branch("PtNearest", &_PtNearest);
  m_tree_AS->Branch("WIDTH", &_WIDTH);
  m_tree_AS->Branch("rbb", &_rbb);
  m_tree_AS->Branch("rfilt", &_rfilt);

  m_tree_AS->Branch("JVF", &_jvf);
  m_tree_AS->Branch("NTRK", &_ntrk);
  m_tree_AS->Branch("TRKPT", &_trkpt);

  m_tree_AS->Branch("QGtag", &_QGtag);
  m_tree_AS->Branch("QGtag2", &_QGtag2);
  m_tree_AS->Branch("jetTrueFlavour", &_jetTrueFlavour);
  m_tree_AS->Branch("IP3DSV1", &w_cmb);
  m_tree_AS->Branch("TrackCounting2D", &w_TrackCounting2D);
  m_tree_AS->Branch("JetProb", &w_JetProb);
  m_tree_AS->Branch("IP1D", &w_IP1D);
  m_tree_AS->Branch("IP2D", &w_IP2D);
  m_tree_AS->Branch("IP3D", &w_IP3D);
  m_tree_AS->Branch("SV0", &w_SV0);
  m_tree_AS->Branch("SV1", &w_SV1);
  m_tree_AS->Branch("SV2", &w_SV2);
  m_tree_AS->Branch("BaselineTagger", &w_BaselineTagger);
  m_tree_AS->Branch("JetFitterTag", &w_JetFitterTag);
  m_tree_AS->Branch("JetFitterCOMB", &w_JetFitterCOMB);
  m_tree_AS->Branch("JetFitterTagNN", &w_JetFitterTagNN);
  m_tree_AS->Branch("JetFitterCOMBNN", &w_JetFitterCOMBNN);
  m_tree_AS->Branch("SoftMuonTag", &w_SoftMuonTag);
  m_tree_AS->Branch("SoftElectronTag", &w_SoftElectronTag);

  if(sc.isFailure()) { 
    ATH_MSG_ERROR("ROOT Hist registration failed"); 
    return sc; 
  }
  
  m_eventNr=0;

  return StatusCode::SUCCESS;
}		 

StatusCode MyAnalysis::finalize() {
  return StatusCode::SUCCESS;
}

StatusCode MyAnalysis::initEvent() {
  /// For Athena-Aware NTuple

  _nparticles = 0;
  _NFPE = 0;

  _jStr_numPVE = 0;
  _jStr_PVxE = 0.;
  _jStr_PVyE = 0.;
  _jStr_PVzE = 0.;
  _jStr_PVrE = 0.;
  _jStr_PVchiSqrE = 0.;
  _jStr_PVnDoFE = 0.;
  _jStr_PVfitE = 0.;
  _jStr_PVnTrkE = 0;
  _jStr_PVtypeE = 0;

  _jStr_numPVJ->clear(); DEBUG(0);
  _jStr_PVxJ->clear(); DEBUG(0);
  _jStr_PVyJ->clear(); DEBUG(0);
  _jStr_PVzJ->clear(); DEBUG(0);
  _jStr_PVrJ->clear(); DEBUG(0);
  _jStr_PVchiSqrJ->clear(); DEBUG(0);
  _jStr_PVnDoFJ->clear(); DEBUG(0);
  _jStr_PVfitJ->clear(); DEBUG(0);
  _jStr_PVnTrkJ->clear(); DEBUG(0);
  _jStr_PVtypeJ->clear(); DEBUG(0);

  _jStr_numJets = 0; DEBUG(0);

  _jStr_evtEJ_etaJ1->clear(); DEBUG(0);
  _jStr_evtEJ_etaJ2->clear(); DEBUG(0);
  _jStr_evtEJ_etaJ3->clear(); DEBUG(0);
  _jStr_evtEJ_etaJ4->clear(); DEBUG(0);
  _jStr_evtEJ_pTJ1->clear(); DEBUG(0);
  _jStr_evtEJ_pTJ2->clear(); DEBUG(0);
  _jStr_evtEJ_pTJ3->clear(); DEBUG(0);
  _jStr_evtEJ_pTJ4->clear(); DEBUG(0);

  _jStr_evtEJ_fE1->clear(); DEBUG(0);
  _jStr_evtEJ_fE2->clear(); DEBUG(0);
  _jStr_evtEJ_fE3->clear(); DEBUG(0);
  _jStr_evtEJ_fET1->clear(); DEBUG(0);
  _jStr_evtEJ_fET2->clear(); DEBUG(0);
  _jStr_evtEJ_fET3->clear(); DEBUG(0);
  _jStr_evtEJ_mct->clear(); DEBUG(0);
  _jStr_evtEJ_q->clear(); DEBUG(0);
  _jStr_evtEJ_mjj->clear(); DEBUG(0);
  _jStr_evtEJ_mTjj->clear(); DEBUG(0);
  _jStr_evtEJ_mjjj->clear(); DEBUG(0);
  _jStr_evtEJ_mTjjj->clear(); DEBUG(0);
  _jStr_evtEJ_mjjjj->clear(); DEBUG(0);
  _jStr_evtEJ_mTjjjj->clear(); DEBUG(0);
  _jStr_evtEJ_zetaPlus->clear(); DEBUG(0);
  _jStr_evtEJ_zetaMinus->clear(); DEBUG(0);
  _jStr_evtEJ_B12->clear(); DEBUG(0);
  _jStr_evtEJ_BT12->clear(); DEBUG(0);
  _jStr_evtEJ_Alpha->clear(); DEBUG(0);
  _jStr_evtEJ_AlphaT->clear(); DEBUG(0);
  _jStr_evtEJ_massDemocracy->clear(); DEBUG(0);
  _jStr_evtEJ_betaflow->clear(); DEBUG(0);
  _jStr_evtEJ_betaflow_1GeV->clear(); DEBUG(0);
  _jStr_evtEJ_betaflow_5GeV->clear(); DEBUG(0);
  _jStr_evtEJ_y23->clear(); DEBUG(0);
  _jStr_evtEJ_y23_1GeV->clear(); DEBUG(0);
  _jStr_evtEJ_y23_5GeV->clear(); DEBUG(0);
  _jStr_evtEJ_lny23->clear(); DEBUG(0);
  _jStr_evtEJ_lny23_1GeV->clear(); DEBUG(0);
  _jStr_evtEJ_lny23_5GeV->clear(); DEBUG(0);
  _jStr_evtEJ_theta->clear(); DEBUG(0);
  _jStr_evtEJ_asym->clear(); DEBUG(0);
  _jStr_evtEJ_yB->clear(); DEBUG(0);
  _jStr_evtEJ_yStar->clear(); DEBUG(0);
  _jStr_evtEJ_thetaStar->clear(); DEBUG(0);
  _jStr_evtEJ_chi->clear(); DEBUG(0);
  _jStr_evtEJ_deltaPhiJJ->clear(); DEBUG(0);
  _jStr_evtEJ_deltaThetaJJ->clear(); DEBUG(0);
  _jStr_evtEJ_deltaEtaJJ->clear(); DEBUG(0);
  _jStr_evtEJ_deltaRapidityJJ->clear(); DEBUG(0);
  _jStr_evtEJ_deltaRJJ->clear(); DEBUG(0);
  _jStr_evtEJ_deltaRJJY->clear(); DEBUG(0);
  _jStr_evtEJ_sigmaPhiJJ->clear(); DEBUG(0);
  _jStr_evtEJ_sigmaThetaJJ->clear(); DEBUG(0);
  _jStr_evtEJ_sigmaEtaJJ->clear(); DEBUG(0);
  _jStr_evtEJ_sigmaRapidityJJ->clear(); DEBUG(0);
  _jStr_evtEJ_sigmaPtJJ->clear(); DEBUG(0);
  _jStr_evtEJ_sigmaEtJJ->clear(); DEBUG(0);
  _jStr_evtEJ_sigmaEt12->clear(); DEBUG(0);
  _jStr_evtEJ_sigmaEt34->clear(); DEBUG(0);
  _jStr_evtEJ_A234->clear(); DEBUG(0);
  _jStr_evtEJ_asymPhiJJ->clear(); DEBUG(0);
  _jStr_evtEJ_asymThetaJJ->clear(); DEBUG(0);
  _jStr_evtEJ_asymEtaJJ->clear(); DEBUG(0);
  _jStr_evtEJ_asymRapidityJJ->clear(); DEBUG(0);
  _jStr_evtEJ_acoplanarity->clear(); DEBUG(0);
  _jStr_evtEJ_twist->clear(); DEBUG(0);
  _jStr_evtEJ_twistY->clear(); DEBUG(0);
  _jStr_evtEJ_jetSumE->clear(); DEBUG(0);
  _jStr_evtEJ_jetSumET->clear(); DEBUG(0);
  _jStr_evtEJ_jetSumPT->clear(); DEBUG(0);
  _jStr_evtEJ_jetSumM->clear(); DEBUG(0);
  _jStr_evtEJ_jetSumMT->clear(); DEBUG(0);
  _jStr_evtEJ_HTprime->clear(); DEBUG(0);
  _jStr_evtEJ_centrality->clear(); DEBUG(0);
  _jStr_evtEJ_centralityP->clear(); DEBUG(0);
  _jStr_evtEJ_zminJ1J2->clear(); DEBUG(0);
  _jStr_evtEJ_zmaxJ1J2->clear(); DEBUG(0);
  _jStr_evtEJ_zminAllJets->clear(); DEBUG(0);
  _jStr_evtEJ_zmaxAllJets->clear(); DEBUG(0);
  _jStr_evtEJ_zminJ1J2Phi->clear(); DEBUG(0);
  _jStr_evtEJ_zminJ1J2Theta->clear(); DEBUG(0);
  _jStr_evtEJ_zminJ1J2Eta->clear(); DEBUG(0);
  _jStr_evtEJ_zminJ1J2Rapidity->clear(); DEBUG(0);
  _jStr_evtEJ_cosHelicityJ1->clear(); DEBUG(0);
  _jStr_evtEJ_helicityJ1->clear(); DEBUG(0);
  _jStr_evtEJ_azilicityJ1->clear(); DEBUG(0);
  _jStr_evtEJ_cosHelicityJ2->clear(); DEBUG(0);
  _jStr_evtEJ_helicityJ2->clear(); DEBUG(0);
  _jStr_evtEJ_azilicityJ2->clear(); DEBUG(0);
  _jStr_evtEJ_cosThetaJ1->clear(); DEBUG(0);
  _jStr_evtEJ_cosThetaJ2->clear(); DEBUG(0);
  _jStr_evtEJ_deltaRapidityXtoJ1CM->clear(); DEBUG(0);
  _jStr_evtEJ_deltaRapidityXtoJ2CM->clear(); DEBUG(0);
  _jStr_evtEJ_deltaRapidityXtoJ1->clear(); DEBUG(0);
  _jStr_evtEJ_deltaRapidityXtoJ2->clear(); DEBUG(0);
  _jStr_evtEJ_cosTheta1->clear(); DEBUG(0);
  _jStr_evtEJ_cosTheta2->clear(); DEBUG(0);
  _jStr_evtEJ_cosThetaStar1->clear(); DEBUG(0);
  _jStr_evtEJ_cosThetaStar2->clear(); DEBUG(0);
  _jStr_evtEJ_cosPhiTilde1->clear(); DEBUG(0);
  _jStr_evtEJ_cosPhiTilde2->clear(); DEBUG(0);
  _jStr_evtEJ_cosPhi->clear(); DEBUG(0);
  _jStr_evtEJ_PhiTilde1->clear(); DEBUG(0);
  _jStr_evtEJ_PhiTilde2->clear(); DEBUG(0);
  _jStr_evtEJ_Phi->clear(); DEBUG(0);
  _jStr_evtEJ_M->clear(); DEBUG(0);
  _jStr_evtEJ_DeltaM->clear(); DEBUG(0);
  _jStr_evtEJ_asymM->clear(); DEBUG(0);
  _jStr_evtEJ_Q->clear(); DEBUG(0);
  _jStr_evtEJ_SCDF->clear(); DEBUG(0);
  _jStr_evtEJ_dPhiIJCDF->clear(); DEBUG(0);
  _jStr_evtEJ_dPhiKLCDF->clear(); DEBUG(0);
  _jStr_evtEJ_PtIPtJCDF->clear(); DEBUG(0);
  _jStr_evtEJ_PtKPtLCDF->clear(); DEBUG(0);
  _jStr_evtEJ_SD0->clear(); DEBUG(0);
  _jStr_evtEJ_dPhiIJD0->clear(); DEBUG(0);
  _jStr_evtEJ_dPhiKLD0->clear(); DEBUG(0);
  _jStr_evtEJ_PtIPtJD0->clear(); DEBUG(0);
  _jStr_evtEJ_PtKPtLD0->clear(); DEBUG(0);
  _jStr_evtEJ_DeltaSCDF->clear(); DEBUG(0);
  _jStr_evtEJ_DeltaSD0->clear(); DEBUG(0);
  _jStr_evtEJ_Delta12->clear(); DEBUG(0);
  _jStr_evtEJ_Delta12n->clear(); DEBUG(0);
  _jStr_evtEJ_M14->clear(); DEBUG(0);
  _jStr_evtEJ_y1y2->clear(); DEBUG(0);

  _jStr_evtEJ_Mprime4->clear(); DEBUG(0);
  _jStr_evtEJ_dMprime4->clear(); DEBUG(0);
  _jStr_evtEJ_MprimeAvg4->clear(); DEBUG(0);
  _jStr_evtEJ_DeltaMin4->clear(); DEBUG(0);
  _jStr_evtEJ_DeltaMax4->clear(); DEBUG(0);
  _jStr_evtEJ_DeltaPhiXX4->clear(); DEBUG(0);
  _jStr_evtEJ_DeltaYXX4->clear(); DEBUG(0);
  _jStr_evtEJ_DeltaRXX4->clear(); DEBUG(0);
  _jStr_evtEJ_TwistXX4->clear(); DEBUG(0);
  _jStr_evtEJ_separated4->clear(); DEBUG(0);
  _jStr_evtEJ_Mprime->clear(); DEBUG(0);
  _jStr_evtEJ_dMprime->clear(); DEBUG(0);
  _jStr_evtEJ_MprimeAvg->clear(); DEBUG(0);
  _jStr_evtEJ_DeltaMin->clear(); DEBUG(0);
  _jStr_evtEJ_DeltaMax->clear(); DEBUG(0);
  _jStr_evtEJ_DeltaPhiXX->clear(); DEBUG(0);
  _jStr_evtEJ_DeltaYXX->clear(); DEBUG(0);
  _jStr_evtEJ_DeltaRXX->clear(); DEBUG(0);
  _jStr_evtEJ_TwistXX->clear(); DEBUG(0);
  _jStr_evtEJ_separated->clear(); DEBUG(0);

  _jStr_evtJj_etaJ1->clear(); DEBUG(0);
  _jStr_evtJj_etaJ2->clear(); DEBUG(0);
  _jStr_evtJj_etaJ3->clear(); DEBUG(0);
  _jStr_evtJj_etaJ4->clear(); DEBUG(0);
  _jStr_evtJj_pTJ1->clear(); DEBUG(0);
  _jStr_evtJj_pTJ2->clear(); DEBUG(0);
  _jStr_evtJj_pTJ3->clear(); DEBUG(0);
  _jStr_evtJj_pTJ4->clear(); DEBUG(0);

  _jStr_evtJj_fE1->clear(); DEBUG(0);
  _jStr_evtJj_fE2->clear(); DEBUG(0);
  _jStr_evtJj_fE3->clear(); DEBUG(0);
  _jStr_evtJj_fET1->clear(); DEBUG(0);
  _jStr_evtJj_fET2->clear(); DEBUG(0);
  _jStr_evtJj_fET3->clear(); DEBUG(0);
  _jStr_evtJj_mct->clear(); DEBUG(0);
  _jStr_evtJj_q->clear(); DEBUG(0);
  _jStr_evtJj_mjj->clear(); DEBUG(0);
  _jStr_evtJj_mTjj->clear(); DEBUG(0);
  _jStr_evtJj_mjjj->clear(); DEBUG(0);
  _jStr_evtJj_mTjjj->clear(); DEBUG(0);
  _jStr_evtJj_mjjjj->clear(); DEBUG(0);
  _jStr_evtJj_mTjjjj->clear(); DEBUG(0);
  _jStr_evtJj_zetaPlus->clear(); DEBUG(0);
  _jStr_evtJj_zetaMinus->clear(); DEBUG(0);
  _jStr_evtJj_B12->clear(); DEBUG(0);
  _jStr_evtJj_BT12->clear(); DEBUG(0);
  _jStr_evtJj_Alpha->clear(); DEBUG(0);
  _jStr_evtJj_AlphaT->clear(); DEBUG(0);
  _jStr_evtJj_massDemocracy->clear(); DEBUG(0);
  _jStr_evtJj_betaflow->clear(); DEBUG(0);
  _jStr_evtJj_betaflow_1GeV->clear(); DEBUG(0);
  _jStr_evtJj_betaflow_5GeV->clear(); DEBUG(0);
  _jStr_evtJj_y23->clear(); DEBUG(0);
  _jStr_evtJj_y23_1GeV->clear(); DEBUG(0);
  _jStr_evtJj_y23_5GeV->clear(); DEBUG(0);
  _jStr_evtJj_lny23->clear(); DEBUG(0);
  _jStr_evtJj_lny23_1GeV->clear(); DEBUG(0);
  _jStr_evtJj_lny23_5GeV->clear(); DEBUG(0);
  _jStr_evtJj_theta->clear(); DEBUG(0);
  _jStr_evtJj_asym->clear(); DEBUG(0);
  _jStr_evtJj_yB->clear(); DEBUG(0);
  _jStr_evtJj_yStar->clear(); DEBUG(0);
  _jStr_evtJj_thetaStar->clear(); DEBUG(0);
  _jStr_evtJj_chi->clear(); DEBUG(0);
  _jStr_evtJj_deltaPhiJJ->clear(); DEBUG(0);
  _jStr_evtJj_deltaThetaJJ->clear(); DEBUG(0);
  _jStr_evtJj_deltaEtaJJ->clear(); DEBUG(0);
  _jStr_evtJj_deltaRapidityJJ->clear(); DEBUG(0);
  _jStr_evtJj_deltaRJJ->clear(); DEBUG(0);
  _jStr_evtJj_deltaRJJY->clear(); DEBUG(0);
  _jStr_evtJj_sigmaPhiJJ->clear(); DEBUG(0);
  _jStr_evtJj_sigmaThetaJJ->clear(); DEBUG(0);
  _jStr_evtJj_sigmaEtaJJ->clear(); DEBUG(0);
  _jStr_evtJj_sigmaRapidityJJ->clear(); DEBUG(0);
  _jStr_evtJj_sigmaPtJJ->clear(); DEBUG(0);
  _jStr_evtJj_sigmaEtJJ->clear(); DEBUG(0);
  _jStr_evtJj_sigmaEt12->clear(); DEBUG(0);
  _jStr_evtJj_sigmaEt34->clear(); DEBUG(0);
  _jStr_evtJj_A234->clear(); DEBUG(0);
  _jStr_evtJj_asymPhiJJ->clear(); DEBUG(0);
  _jStr_evtJj_asymThetaJJ->clear(); DEBUG(0);
  _jStr_evtJj_asymEtaJJ->clear(); DEBUG(0);
  _jStr_evtJj_asymRapidityJJ->clear(); DEBUG(0);
  _jStr_evtJj_acoplanarity->clear(); DEBUG(0);
  _jStr_evtJj_twist->clear(); DEBUG(0);
  _jStr_evtJj_twistY->clear(); DEBUG(0);
  _jStr_evtJj_jetSumE->clear(); DEBUG(0);
  _jStr_evtJj_jetSumET->clear(); DEBUG(0);
  _jStr_evtJj_jetSumPT->clear(); DEBUG(0);
  _jStr_evtJj_jetSumM->clear(); DEBUG(0);
  _jStr_evtJj_jetSumMT->clear(); DEBUG(0);
  _jStr_evtJj_HTprime->clear(); DEBUG(0);
  _jStr_evtJj_centrality->clear(); DEBUG(0);
  _jStr_evtJj_centralityP->clear(); DEBUG(0);
  _jStr_evtJj_zminJ1J2->clear(); DEBUG(0);
  _jStr_evtJj_zmaxJ1J2->clear(); DEBUG(0);
  _jStr_evtJj_zminAllJets->clear(); DEBUG(0);
  _jStr_evtJj_zmaxAllJets->clear(); DEBUG(0);
  _jStr_evtJj_zminJ1J2Phi->clear(); DEBUG(0);
  _jStr_evtJj_zminJ1J2Theta->clear(); DEBUG(0);
  _jStr_evtJj_zminJ1J2Eta->clear(); DEBUG(0);
  _jStr_evtJj_zminJ1J2Rapidity->clear(); DEBUG(0);
  _jStr_evtJj_cosHelicityJ1->clear(); DEBUG(0);
  _jStr_evtJj_helicityJ1->clear(); DEBUG(0);
  _jStr_evtJj_azilicityJ1->clear(); DEBUG(0);
  _jStr_evtJj_cosHelicityJ2->clear(); DEBUG(0);
  _jStr_evtJj_helicityJ2->clear(); DEBUG(0);
  _jStr_evtJj_azilicityJ2->clear(); DEBUG(0);
  _jStr_evtJj_cosThetaJ1->clear(); DEBUG(0);
  _jStr_evtJj_cosThetaJ2->clear(); DEBUG(0);
  _jStr_evtJj_deltaRapidityXtoJ1CM->clear(); DEBUG(0);
  _jStr_evtJj_deltaRapidityXtoJ2CM->clear(); DEBUG(0);
  _jStr_evtJj_deltaRapidityXtoJ1->clear(); DEBUG(0);
  _jStr_evtJj_deltaRapidityXtoJ2->clear(); DEBUG(0);
  _jStr_evtJj_cosTheta1->clear(); DEBUG(0);
  _jStr_evtJj_cosTheta2->clear(); DEBUG(0);
  _jStr_evtJj_cosThetaStar1->clear(); DEBUG(0);
  _jStr_evtJj_cosThetaStar2->clear(); DEBUG(0);
  _jStr_evtJj_cosPhiTilde1->clear(); DEBUG(0);
  _jStr_evtJj_cosPhiTilde2->clear(); DEBUG(0);
  _jStr_evtJj_cosPhi->clear(); DEBUG(0);
  _jStr_evtJj_PhiTilde1->clear(); DEBUG(0);
  _jStr_evtJj_PhiTilde2->clear(); DEBUG(0);
  _jStr_evtJj_Phi->clear(); DEBUG(0);
  _jStr_evtJj_M->clear(); DEBUG(0);
  _jStr_evtJj_DeltaM->clear(); DEBUG(0);
  _jStr_evtJj_asymM->clear(); DEBUG(0);
  _jStr_evtJj_Q->clear(); DEBUG(0);
  _jStr_evtJj_SCDF->clear(); DEBUG(0);
  _jStr_evtJj_dPhiIJCDF->clear(); DEBUG(0);
  _jStr_evtJj_dPhiKLCDF->clear(); DEBUG(0);
  _jStr_evtJj_PtIPtJCDF->clear(); DEBUG(0);
  _jStr_evtJj_PtKPtLCDF->clear(); DEBUG(0);
  _jStr_evtJj_SD0->clear(); DEBUG(0);
  _jStr_evtJj_dPhiIJD0->clear(); DEBUG(0);
  _jStr_evtJj_dPhiKLD0->clear(); DEBUG(0);
  _jStr_evtJj_PtIPtJD0->clear(); DEBUG(0);
  _jStr_evtJj_PtKPtLD0->clear(); DEBUG(0);
  _jStr_evtJj_DeltaSCDF->clear(); DEBUG(0);
  _jStr_evtJj_DeltaSD0->clear(); DEBUG(0);
  _jStr_evtJj_Delta12->clear(); DEBUG(0);
  _jStr_evtJj_Delta12n->clear(); DEBUG(0);
  _jStr_evtJj_M14->clear(); DEBUG(0);
  _jStr_evtJj_y1y2->clear(); DEBUG(0);

  _jStr_evtJj_Mprime4->clear(); DEBUG(0);
  _jStr_evtJj_dMprime4->clear(); DEBUG(0);
  _jStr_evtJj_MprimeAvg4->clear(); DEBUG(0);
  _jStr_evtJj_DeltaMin4->clear(); DEBUG(0);
  _jStr_evtJj_DeltaMax4->clear(); DEBUG(0);
  _jStr_evtJj_DeltaPhiXX4->clear(); DEBUG(0);
  _jStr_evtJj_DeltaYXX4->clear(); DEBUG(0);
  _jStr_evtJj_DeltaRXX4->clear(); DEBUG(0);
  _jStr_evtJj_TwistXX4->clear(); DEBUG(0);
  _jStr_evtJj_separated4->clear(); DEBUG(0);
  _jStr_evtJj_Mprime->clear(); DEBUG(0);
  _jStr_evtJj_dMprime->clear(); DEBUG(0);
  _jStr_evtJj_MprimeAvg->clear(); DEBUG(0);
  _jStr_evtJj_DeltaMin->clear(); DEBUG(0);
  _jStr_evtJj_DeltaMax->clear(); DEBUG(0);
  _jStr_evtJj_DeltaPhiXX->clear(); DEBUG(0);
  _jStr_evtJj_DeltaYXX->clear(); DEBUG(0);
  _jStr_evtJj_DeltaRXX->clear(); DEBUG(0);
  _jStr_evtJj_TwistXX->clear(); DEBUG(0);
  _jStr_evtJj_separated->clear(); DEBUG(0);

  _jStr_fwmEJ_xiPlus->clear(); DEBUG(0);
  _jStr_fwmEJ_xiMinus->clear(); DEBUG(0);
  _jStr_fwmEJ_xPlus->clear(); DEBUG(0);
  _jStr_fwmEJ_xMinus->clear(); DEBUG(0);
  _jStr_fwmEJ_Psi1->clear(); DEBUG(0);
  _jStr_fwmEJ_B0->clear(); DEBUG(0);
  _jStr_fwmEJ_B1->clear(); DEBUG(0);
  _jStr_fwmEJ_B2->clear(); DEBUG(0);
  _jStr_fwmEJ_B3->clear(); DEBUG(0);
  _jStr_fwmEJ_B4->clear(); DEBUG(0);
  _jStr_fwmEJ_B5->clear(); DEBUG(0);
  _jStr_fwmEJ_B6->clear(); DEBUG(0);
  _jStr_fwmEJ_B7->clear(); DEBUG(0);
  _jStr_fwmEJ_B8->clear(); DEBUG(0);
  _jStr_fwmEJ_C0->clear(); DEBUG(0);
  _jStr_fwmEJ_C1->clear(); DEBUG(0);
  _jStr_fwmEJ_C2->clear(); DEBUG(0);
  _jStr_fwmEJ_C3->clear(); DEBUG(0);
  _jStr_fwmEJ_C4->clear(); DEBUG(0);
  _jStr_fwmEJ_C5->clear(); DEBUG(0);
  _jStr_fwmEJ_C6->clear(); DEBUG(0);
  _jStr_fwmEJ_C7->clear(); DEBUG(0);
  _jStr_fwmEJ_C8->clear(); DEBUG(0);
  _jStr_fwmEJ_K0->clear(); DEBUG(0);
  _jStr_fwmEJ_K1->clear(); DEBUG(0);
  _jStr_fwmEJ_K2->clear(); DEBUG(0);
  _jStr_fwmEJ_K3->clear(); DEBUG(0);
  _jStr_fwmEJ_K4->clear(); DEBUG(0);
  _jStr_fwmEJ_K5->clear(); DEBUG(0);
  _jStr_fwmEJ_K6->clear(); DEBUG(0);
  _jStr_fwmEJ_K7->clear(); DEBUG(0);
  _jStr_fwmEJ_K8->clear(); DEBUG(0);
  _jStr_fwmEJ_D0->clear(); DEBUG(0);
  _jStr_fwmEJ_D1->clear(); DEBUG(0);
  _jStr_fwmEJ_D2->clear(); DEBUG(0);
  _jStr_fwmEJ_D3->clear(); DEBUG(0);
  _jStr_fwmEJ_D4->clear(); DEBUG(0);
  _jStr_fwmEJ_D5->clear(); DEBUG(0);
  _jStr_fwmEJ_D6->clear(); DEBUG(0);
  _jStr_fwmEJ_D7->clear(); DEBUG(0);
  _jStr_fwmEJ_D8->clear(); DEBUG(0);
  _jStr_fwmEJ_H0->clear(); DEBUG(0);
  _jStr_fwmEJ_H1->clear(); DEBUG(0);
  _jStr_fwmEJ_H2->clear(); DEBUG(0);
  _jStr_fwmEJ_H3->clear(); DEBUG(0);
  _jStr_fwmEJ_H4->clear(); DEBUG(0);
  _jStr_fwmEJ_H5->clear(); DEBUG(0);
  _jStr_fwmEJ_H6->clear(); DEBUG(0);
  _jStr_fwmEJ_H7->clear(); DEBUG(0);
  _jStr_fwmEJ_H8->clear(); DEBUG(0);
  _jStr_fwmEJ_Q0->clear(); DEBUG(0);
  _jStr_fwmEJ_Q1->clear(); DEBUG(0);
  _jStr_fwmEJ_Q2->clear(); DEBUG(0);
  _jStr_fwmEJ_Q3->clear(); DEBUG(0);
  _jStr_fwmEJ_Q4->clear(); DEBUG(0);
  _jStr_fwmEJ_Q5->clear(); DEBUG(0);
  _jStr_fwmEJ_Q6->clear(); DEBUG(0);
  _jStr_fwmEJ_Q7->clear(); DEBUG(0);
  _jStr_fwmEJ_Q8->clear(); DEBUG(0);
  _jStr_fwmEJ_Pi1->clear(); DEBUG(0);
  _jStr_fwmEJ_Pi2->clear(); DEBUG(0);
  _jStr_fwmEJ_Pi3->clear(); DEBUG(0);
  _jStr_fwmEJ_Pi4->clear(); DEBUG(0);
  _jStr_fwmEJ_B10->clear(); DEBUG(0);
  _jStr_fwmEJ_B20->clear(); DEBUG(0);
  _jStr_fwmEJ_B30->clear(); DEBUG(0);
  _jStr_fwmEJ_B40->clear(); DEBUG(0);
  _jStr_fwmEJ_B50->clear(); DEBUG(0);
  _jStr_fwmEJ_B60->clear(); DEBUG(0);
  _jStr_fwmEJ_B70->clear(); DEBUG(0);
  _jStr_fwmEJ_B80->clear(); DEBUG(0);
  _jStr_fwmEJ_C10->clear(); DEBUG(0);
  _jStr_fwmEJ_C20->clear(); DEBUG(0);
  _jStr_fwmEJ_C30->clear(); DEBUG(0);
  _jStr_fwmEJ_C40->clear(); DEBUG(0);
  _jStr_fwmEJ_C50->clear(); DEBUG(0);
  _jStr_fwmEJ_C60->clear(); DEBUG(0);
  _jStr_fwmEJ_C70->clear(); DEBUG(0);
  _jStr_fwmEJ_C80->clear(); DEBUG(0);
  _jStr_fwmEJ_K10->clear(); DEBUG(0);
  _jStr_fwmEJ_K20->clear(); DEBUG(0);
  _jStr_fwmEJ_K30->clear(); DEBUG(0);
  _jStr_fwmEJ_K40->clear(); DEBUG(0);
  _jStr_fwmEJ_K50->clear(); DEBUG(0);
  _jStr_fwmEJ_K60->clear(); DEBUG(0);
  _jStr_fwmEJ_K70->clear(); DEBUG(0);
  _jStr_fwmEJ_K80->clear(); DEBUG(0);
  _jStr_fwmEJ_D10->clear(); DEBUG(0);
  _jStr_fwmEJ_D20->clear(); DEBUG(0);
  _jStr_fwmEJ_D30->clear(); DEBUG(0);
  _jStr_fwmEJ_D40->clear(); DEBUG(0);
  _jStr_fwmEJ_D50->clear(); DEBUG(0);
  _jStr_fwmEJ_D60->clear(); DEBUG(0);
  _jStr_fwmEJ_D70->clear(); DEBUG(0);
  _jStr_fwmEJ_D80->clear(); DEBUG(0);
  _jStr_fwmEJ_H10->clear(); DEBUG(0);
  _jStr_fwmEJ_H20->clear(); DEBUG(0);
  _jStr_fwmEJ_H30->clear(); DEBUG(0);
  _jStr_fwmEJ_H40->clear(); DEBUG(0);
  _jStr_fwmEJ_H50->clear(); DEBUG(0);
  _jStr_fwmEJ_H60->clear(); DEBUG(0);
  _jStr_fwmEJ_H70->clear(); DEBUG(0);
  _jStr_fwmEJ_H80->clear(); DEBUG(0);
  _jStr_fwmEJ_Q10->clear(); DEBUG(0);
  _jStr_fwmEJ_Q20->clear(); DEBUG(0);
  _jStr_fwmEJ_Q30->clear(); DEBUG(0);
  _jStr_fwmEJ_Q40->clear(); DEBUG(0);
  _jStr_fwmEJ_Q50->clear(); DEBUG(0);
  _jStr_fwmEJ_Q60->clear(); DEBUG(0);
  _jStr_fwmEJ_Q70->clear(); DEBUG(0);
  _jStr_fwmEJ_Q80->clear(); DEBUG(0);
  // _jStr_fwmEi_Psi1->clear(); DEBUG(0);
  // _jStr_fwmEi_B0->clear(); DEBUG(0);
  // _jStr_fwmEi_B1->clear(); DEBUG(0);
  // _jStr_fwmEi_B2->clear(); DEBUG(0);
  // _jStr_fwmEi_B3->clear(); DEBUG(0);
  // _jStr_fwmEi_B4->clear(); DEBUG(0);
  // _jStr_fwmEi_B5->clear(); DEBUG(0);
  // _jStr_fwmEi_B6->clear(); DEBUG(0);
  // _jStr_fwmEi_B7->clear(); DEBUG(0);
  // _jStr_fwmEi_B8->clear(); DEBUG(0);
  // _jStr_fwmEi_C0->clear(); DEBUG(0);
  // _jStr_fwmEi_C1->clear(); DEBUG(0);
  // _jStr_fwmEi_C2->clear(); DEBUG(0);
  // _jStr_fwmEi_C3->clear(); DEBUG(0);
  // _jStr_fwmEi_C4->clear(); DEBUG(0);
  // _jStr_fwmEi_C5->clear(); DEBUG(0);
  // _jStr_fwmEi_C6->clear(); DEBUG(0);
  // _jStr_fwmEi_C7->clear(); DEBUG(0);
  // _jStr_fwmEi_C8->clear(); DEBUG(0);
  // _jStr_fwmEi_K0->clear(); DEBUG(0);
  // _jStr_fwmEi_K1->clear(); DEBUG(0);
  // _jStr_fwmEi_K2->clear(); DEBUG(0);
  // _jStr_fwmEi_K3->clear(); DEBUG(0);
  // _jStr_fwmEi_K4->clear(); DEBUG(0);
  // _jStr_fwmEi_K5->clear(); DEBUG(0);
  // _jStr_fwmEi_K6->clear(); DEBUG(0);
  // _jStr_fwmEi_K7->clear(); DEBUG(0);
  // _jStr_fwmEi_K8->clear(); DEBUG(0);
  // _jStr_fwmEi_D0->clear(); DEBUG(0);
  // _jStr_fwmEi_D1->clear(); DEBUG(0);
  // _jStr_fwmEi_D2->clear(); DEBUG(0);
  // _jStr_fwmEi_D3->clear(); DEBUG(0);
  // _jStr_fwmEi_D4->clear(); DEBUG(0);
  // _jStr_fwmEi_D5->clear(); DEBUG(0);
  // _jStr_fwmEi_D6->clear(); DEBUG(0);
  // _jStr_fwmEi_D7->clear(); DEBUG(0);
  // _jStr_fwmEi_D8->clear(); DEBUG(0);
  // _jStr_fwmEi_H0->clear(); DEBUG(0);
  // _jStr_fwmEi_H1->clear(); DEBUG(0);
  // _jStr_fwmEi_H2->clear(); DEBUG(0);
  // _jStr_fwmEi_H3->clear(); DEBUG(0);
  // _jStr_fwmEi_H4->clear(); DEBUG(0);
  // _jStr_fwmEi_H5->clear(); DEBUG(0);
  // _jStr_fwmEi_H6->clear(); DEBUG(0);
  // _jStr_fwmEi_H7->clear(); DEBUG(0);
  // _jStr_fwmEi_H8->clear(); DEBUG(0);
  // _jStr_fwmEi_Q0->clear(); DEBUG(0);
  // _jStr_fwmEi_Q1->clear(); DEBUG(0);
  // _jStr_fwmEi_Q2->clear(); DEBUG(0);
  // _jStr_fwmEi_Q3->clear(); DEBUG(0);
  // _jStr_fwmEi_Q4->clear(); DEBUG(0);
  // _jStr_fwmEi_Q5->clear(); DEBUG(0);
  // _jStr_fwmEi_Q6->clear(); DEBUG(0);
  // _jStr_fwmEi_Q7->clear(); DEBUG(0);
  // _jStr_fwmEi_Q8->clear(); DEBUG(0);
  // _jStr_fwmEi_Pi1->clear(); DEBUG(0);
  // _jStr_fwmEi_Pi2->clear(); DEBUG(0);
  // _jStr_fwmEi_Pi3->clear(); DEBUG(0);
  // _jStr_fwmEi_Pi4->clear(); DEBUG(0);
  // _jStr_fwmEi_B10->clear(); DEBUG(0);
  // _jStr_fwmEi_B20->clear(); DEBUG(0);
  // _jStr_fwmEi_B30->clear(); DEBUG(0);
  // _jStr_fwmEi_B40->clear(); DEBUG(0);
  // _jStr_fwmEi_B50->clear(); DEBUG(0);
  // _jStr_fwmEi_B60->clear(); DEBUG(0);
  // _jStr_fwmEi_B70->clear(); DEBUG(0);
  // _jStr_fwmEi_B80->clear(); DEBUG(0);
  // _jStr_fwmEi_C10->clear(); DEBUG(0);
  // _jStr_fwmEi_C20->clear(); DEBUG(0);
  // _jStr_fwmEi_C30->clear(); DEBUG(0);
  // _jStr_fwmEi_C40->clear(); DEBUG(0);
  // _jStr_fwmEi_C50->clear(); DEBUG(0);
  // _jStr_fwmEi_C60->clear(); DEBUG(0);
  // _jStr_fwmEi_C70->clear(); DEBUG(0);
  // _jStr_fwmEi_C80->clear(); DEBUG(0);
  // _jStr_fwmEi_K10->clear(); DEBUG(0);
  // _jStr_fwmEi_K20->clear(); DEBUG(0);
  // _jStr_fwmEi_K30->clear(); DEBUG(0);
  // _jStr_fwmEi_K40->clear(); DEBUG(0);
  // _jStr_fwmEi_K50->clear(); DEBUG(0);
  // _jStr_fwmEi_K60->clear(); DEBUG(0);
  // _jStr_fwmEi_K70->clear(); DEBUG(0);
  // _jStr_fwmEi_K80->clear(); DEBUG(0);
  // _jStr_fwmEi_D10->clear(); DEBUG(0);
  // _jStr_fwmEi_D20->clear(); DEBUG(0);
  // _jStr_fwmEi_D30->clear(); DEBUG(0);
  // _jStr_fwmEi_D40->clear(); DEBUG(0);
  // _jStr_fwmEi_D50->clear(); DEBUG(0);
  // _jStr_fwmEi_D60->clear(); DEBUG(0);
  // _jStr_fwmEi_D70->clear(); DEBUG(0);
  // _jStr_fwmEi_D80->clear(); DEBUG(0);
  // _jStr_fwmEi_H10->clear(); DEBUG(0);
  // _jStr_fwmEi_H20->clear(); DEBUG(0);
  // _jStr_fwmEi_H30->clear(); DEBUG(0);
  // _jStr_fwmEi_H40->clear(); DEBUG(0);
  // _jStr_fwmEi_H50->clear(); DEBUG(0);
  // _jStr_fwmEi_H60->clear(); DEBUG(0);
  // _jStr_fwmEi_H70->clear(); DEBUG(0);
  // _jStr_fwmEi_H80->clear(); DEBUG(0);
  // _jStr_fwmEi_Q10->clear(); DEBUG(0);
  // _jStr_fwmEi_Q20->clear(); DEBUG(0);
  // _jStr_fwmEi_Q30->clear(); DEBUG(0);
  // _jStr_fwmEi_Q40->clear(); DEBUG(0);
  // _jStr_fwmEi_Q50->clear(); DEBUG(0);
  // _jStr_fwmEi_Q60->clear(); DEBUG(0);
  // _jStr_fwmEi_Q70->clear(); DEBUG(0);
  // _jStr_fwmEi_Q80->clear(); DEBUG(0);
  _jStr_fwmJj_xiPlus->clear(); DEBUG(0);
  _jStr_fwmJj_xiMinus->clear(); DEBUG(0);
  _jStr_fwmJj_xPlus->clear(); DEBUG(0);
  _jStr_fwmJj_xMinus->clear(); DEBUG(0);
  _jStr_fwmJj_Psi1->clear(); DEBUG(0);
  _jStr_fwmJj_B0->clear(); DEBUG(0);
  _jStr_fwmJj_B1->clear(); DEBUG(0);
  _jStr_fwmJj_B2->clear(); DEBUG(0);
  _jStr_fwmJj_B3->clear(); DEBUG(0);
  _jStr_fwmJj_B4->clear(); DEBUG(0);
  _jStr_fwmJj_B5->clear(); DEBUG(0);
  _jStr_fwmJj_B6->clear(); DEBUG(0);
  _jStr_fwmJj_B7->clear(); DEBUG(0);
  _jStr_fwmJj_B8->clear(); DEBUG(0);
  _jStr_fwmJj_C0->clear(); DEBUG(0);
  _jStr_fwmJj_C1->clear(); DEBUG(0);
  _jStr_fwmJj_C2->clear(); DEBUG(0);
  _jStr_fwmJj_C3->clear(); DEBUG(0);
  _jStr_fwmJj_C4->clear(); DEBUG(0);
  _jStr_fwmJj_C5->clear(); DEBUG(0);
  _jStr_fwmJj_C6->clear(); DEBUG(0);
  _jStr_fwmJj_C7->clear(); DEBUG(0);
  _jStr_fwmJj_C8->clear(); DEBUG(0);
  _jStr_fwmJj_K0->clear(); DEBUG(0);
  _jStr_fwmJj_K1->clear(); DEBUG(0);
  _jStr_fwmJj_K2->clear(); DEBUG(0);
  _jStr_fwmJj_K3->clear(); DEBUG(0);
  _jStr_fwmJj_K4->clear(); DEBUG(0);
  _jStr_fwmJj_K5->clear(); DEBUG(0);
  _jStr_fwmJj_K6->clear(); DEBUG(0);
  _jStr_fwmJj_K7->clear(); DEBUG(0);
  _jStr_fwmJj_K8->clear(); DEBUG(0);
  _jStr_fwmJj_D0->clear(); DEBUG(0);
  _jStr_fwmJj_D1->clear(); DEBUG(0);
  _jStr_fwmJj_D2->clear(); DEBUG(0);
  _jStr_fwmJj_D3->clear(); DEBUG(0);
  _jStr_fwmJj_D4->clear(); DEBUG(0);
  _jStr_fwmJj_D5->clear(); DEBUG(0);
  _jStr_fwmJj_D6->clear(); DEBUG(0);
  _jStr_fwmJj_D7->clear(); DEBUG(0);
  _jStr_fwmJj_D8->clear(); DEBUG(0);
  _jStr_fwmJj_H0->clear(); DEBUG(0);
  _jStr_fwmJj_H1->clear(); DEBUG(0);
  _jStr_fwmJj_H2->clear(); DEBUG(0);
  _jStr_fwmJj_H3->clear(); DEBUG(0);
  _jStr_fwmJj_H4->clear(); DEBUG(0);
  _jStr_fwmJj_H5->clear(); DEBUG(0);
  _jStr_fwmJj_H6->clear(); DEBUG(0);
  _jStr_fwmJj_H7->clear(); DEBUG(0);
  _jStr_fwmJj_H8->clear(); DEBUG(0);
  _jStr_fwmJj_Q0->clear(); DEBUG(0);
  _jStr_fwmJj_Q1->clear(); DEBUG(0);
  _jStr_fwmJj_Q2->clear(); DEBUG(0);
  _jStr_fwmJj_Q3->clear(); DEBUG(0);
  _jStr_fwmJj_Q4->clear(); DEBUG(0);
  _jStr_fwmJj_Q5->clear(); DEBUG(0);
  _jStr_fwmJj_Q6->clear(); DEBUG(0);
  _jStr_fwmJj_Q7->clear(); DEBUG(0);
  _jStr_fwmJj_Q8->clear(); DEBUG(0);
  _jStr_fwmJj_Pi1->clear(); DEBUG(0);
  _jStr_fwmJj_Pi2->clear(); DEBUG(0);
  _jStr_fwmJj_Pi3->clear(); DEBUG(0);
  _jStr_fwmJj_Pi4->clear(); DEBUG(0);
  _jStr_fwmJj_B10->clear(); DEBUG(0);
  _jStr_fwmJj_B20->clear(); DEBUG(0);
  _jStr_fwmJj_B30->clear(); DEBUG(0);
  _jStr_fwmJj_B40->clear(); DEBUG(0);
  _jStr_fwmJj_B50->clear(); DEBUG(0);
  _jStr_fwmJj_B60->clear(); DEBUG(0);
  _jStr_fwmJj_B70->clear(); DEBUG(0);
  _jStr_fwmJj_B80->clear(); DEBUG(0);
  _jStr_fwmJj_C10->clear(); DEBUG(0);
  _jStr_fwmJj_C20->clear(); DEBUG(0);
  _jStr_fwmJj_C30->clear(); DEBUG(0);
  _jStr_fwmJj_C40->clear(); DEBUG(0);
  _jStr_fwmJj_C50->clear(); DEBUG(0);
  _jStr_fwmJj_C60->clear(); DEBUG(0);
  _jStr_fwmJj_C70->clear(); DEBUG(0);
  _jStr_fwmJj_C80->clear(); DEBUG(0);
  _jStr_fwmJj_K10->clear(); DEBUG(0);
  _jStr_fwmJj_K20->clear(); DEBUG(0);
  _jStr_fwmJj_K30->clear(); DEBUG(0);
  _jStr_fwmJj_K40->clear(); DEBUG(0);
  _jStr_fwmJj_K50->clear(); DEBUG(0);
  _jStr_fwmJj_K60->clear(); DEBUG(0);
  _jStr_fwmJj_K70->clear(); DEBUG(0);
  _jStr_fwmJj_K80->clear(); DEBUG(0);
  _jStr_fwmJj_D10->clear(); DEBUG(0);
  _jStr_fwmJj_D20->clear(); DEBUG(0);
  _jStr_fwmJj_D30->clear(); DEBUG(0);
  _jStr_fwmJj_D40->clear(); DEBUG(0);
  _jStr_fwmJj_D50->clear(); DEBUG(0);
  _jStr_fwmJj_D60->clear(); DEBUG(0);
  _jStr_fwmJj_D70->clear(); DEBUG(0);
  _jStr_fwmJj_D80->clear(); DEBUG(0);
  _jStr_fwmJj_H10->clear(); DEBUG(0);
  _jStr_fwmJj_H20->clear(); DEBUG(0);
  _jStr_fwmJj_H30->clear(); DEBUG(0);
  _jStr_fwmJj_H40->clear(); DEBUG(0);
  _jStr_fwmJj_H50->clear(); DEBUG(0);
  _jStr_fwmJj_H60->clear(); DEBUG(0);
  _jStr_fwmJj_H70->clear(); DEBUG(0);
  _jStr_fwmJj_H80->clear(); DEBUG(0);
  _jStr_fwmJj_Q10->clear(); DEBUG(0);
  _jStr_fwmJj_Q20->clear(); DEBUG(0);
  _jStr_fwmJj_Q30->clear(); DEBUG(0);
  _jStr_fwmJj_Q40->clear(); DEBUG(0);
  _jStr_fwmJj_Q50->clear(); DEBUG(0);
  _jStr_fwmJj_Q60->clear(); DEBUG(0);
  _jStr_fwmJj_Q70->clear(); DEBUG(0);
  _jStr_fwmJj_Q80->clear(); DEBUG(0);
  _jStr_fwmJi_xiPlus->clear(); DEBUG(0);
  _jStr_fwmJi_xiMinus->clear(); DEBUG(0);
  _jStr_fwmJi_xPlus->clear(); DEBUG(0);
  _jStr_fwmJi_xMinus->clear(); DEBUG(0);
  _jStr_fwmJi_Psi1->clear(); DEBUG(0);
  _jStr_fwmJi_B0->clear(); DEBUG(0);
  _jStr_fwmJi_B1->clear(); DEBUG(0);
  _jStr_fwmJi_B2->clear(); DEBUG(0);
  _jStr_fwmJi_B3->clear(); DEBUG(0);
  _jStr_fwmJi_B4->clear(); DEBUG(0);
  _jStr_fwmJi_B5->clear(); DEBUG(0);
  _jStr_fwmJi_B6->clear(); DEBUG(0);
  _jStr_fwmJi_B7->clear(); DEBUG(0);
  _jStr_fwmJi_B8->clear(); DEBUG(0);
  _jStr_fwmJi_C0->clear(); DEBUG(0);
  _jStr_fwmJi_C1->clear(); DEBUG(0);
  _jStr_fwmJi_C2->clear(); DEBUG(0);
  _jStr_fwmJi_C3->clear(); DEBUG(0);
  _jStr_fwmJi_C4->clear(); DEBUG(0);
  _jStr_fwmJi_C5->clear(); DEBUG(0);
  _jStr_fwmJi_C6->clear(); DEBUG(0);
  _jStr_fwmJi_C7->clear(); DEBUG(0);
  _jStr_fwmJi_C8->clear(); DEBUG(0);
  _jStr_fwmJi_K0->clear(); DEBUG(0);
  _jStr_fwmJi_K1->clear(); DEBUG(0);
  _jStr_fwmJi_K2->clear(); DEBUG(0);
  _jStr_fwmJi_K3->clear(); DEBUG(0);
  _jStr_fwmJi_K4->clear(); DEBUG(0);
  _jStr_fwmJi_K5->clear(); DEBUG(0);
  _jStr_fwmJi_K6->clear(); DEBUG(0);
  _jStr_fwmJi_K7->clear(); DEBUG(0);
  _jStr_fwmJi_K8->clear(); DEBUG(0);
  _jStr_fwmJi_D0->clear(); DEBUG(0);
  _jStr_fwmJi_D1->clear(); DEBUG(0);
  _jStr_fwmJi_D2->clear(); DEBUG(0);
  _jStr_fwmJi_D3->clear(); DEBUG(0);
  _jStr_fwmJi_D4->clear(); DEBUG(0);
  _jStr_fwmJi_D5->clear(); DEBUG(0);
  _jStr_fwmJi_D6->clear(); DEBUG(0);
  _jStr_fwmJi_D7->clear(); DEBUG(0);
  _jStr_fwmJi_D8->clear(); DEBUG(0);
  _jStr_fwmJi_H0->clear(); DEBUG(0);
  _jStr_fwmJi_H1->clear(); DEBUG(0);
  _jStr_fwmJi_H2->clear(); DEBUG(0);
  _jStr_fwmJi_H3->clear(); DEBUG(0);
  _jStr_fwmJi_H4->clear(); DEBUG(0);
  _jStr_fwmJi_H5->clear(); DEBUG(0);
  _jStr_fwmJi_H6->clear(); DEBUG(0);
  _jStr_fwmJi_H7->clear(); DEBUG(0);
  _jStr_fwmJi_H8->clear(); DEBUG(0);
  _jStr_fwmJi_Q0->clear(); DEBUG(0);
  _jStr_fwmJi_Q1->clear(); DEBUG(0);
  _jStr_fwmJi_Q2->clear(); DEBUG(0);
  _jStr_fwmJi_Q3->clear(); DEBUG(0);
  _jStr_fwmJi_Q4->clear(); DEBUG(0);
  _jStr_fwmJi_Q5->clear(); DEBUG(0);
  _jStr_fwmJi_Q6->clear(); DEBUG(0);
  _jStr_fwmJi_Q7->clear(); DEBUG(0);
  _jStr_fwmJi_Q8->clear(); DEBUG(0);
  _jStr_fwmJi_Pi1->clear(); DEBUG(0);
  _jStr_fwmJi_Pi2->clear(); DEBUG(0);
  _jStr_fwmJi_Pi3->clear(); DEBUG(0);
  _jStr_fwmJi_Pi4->clear(); DEBUG(0);
  _jStr_fwmJi_B10->clear(); DEBUG(0);
  _jStr_fwmJi_B20->clear(); DEBUG(0);
  _jStr_fwmJi_B30->clear(); DEBUG(0);
  _jStr_fwmJi_B40->clear(); DEBUG(0);
  _jStr_fwmJi_B50->clear(); DEBUG(0);
  _jStr_fwmJi_B60->clear(); DEBUG(0);
  _jStr_fwmJi_B70->clear(); DEBUG(0);
  _jStr_fwmJi_B80->clear(); DEBUG(0);
  _jStr_fwmJi_C10->clear(); DEBUG(0);
  _jStr_fwmJi_C20->clear(); DEBUG(0);
  _jStr_fwmJi_C30->clear(); DEBUG(0);
  _jStr_fwmJi_C40->clear(); DEBUG(0);
  _jStr_fwmJi_C50->clear(); DEBUG(0);
  _jStr_fwmJi_C60->clear(); DEBUG(0);
  _jStr_fwmJi_C70->clear(); DEBUG(0);
  _jStr_fwmJi_C80->clear(); DEBUG(0);
  _jStr_fwmJi_K10->clear(); DEBUG(0);
  _jStr_fwmJi_K20->clear(); DEBUG(0);
  _jStr_fwmJi_K30->clear(); DEBUG(0);
  _jStr_fwmJi_K40->clear(); DEBUG(0);
  _jStr_fwmJi_K50->clear(); DEBUG(0);
  _jStr_fwmJi_K60->clear(); DEBUG(0);
  _jStr_fwmJi_K70->clear(); DEBUG(0);
  _jStr_fwmJi_K80->clear(); DEBUG(0);
  _jStr_fwmJi_D10->clear(); DEBUG(0);
  _jStr_fwmJi_D20->clear(); DEBUG(0);
  _jStr_fwmJi_D30->clear(); DEBUG(0);
  _jStr_fwmJi_D40->clear(); DEBUG(0);
  _jStr_fwmJi_D50->clear(); DEBUG(0);
  _jStr_fwmJi_D60->clear(); DEBUG(0);
  _jStr_fwmJi_D70->clear(); DEBUG(0);
  _jStr_fwmJi_D80->clear(); DEBUG(0);
  _jStr_fwmJi_H10->clear(); DEBUG(0);
  _jStr_fwmJi_H20->clear(); DEBUG(0);
  _jStr_fwmJi_H30->clear(); DEBUG(0);
  _jStr_fwmJi_H40->clear(); DEBUG(0);
  _jStr_fwmJi_H50->clear(); DEBUG(0);
  _jStr_fwmJi_H60->clear(); DEBUG(0);
  _jStr_fwmJi_H70->clear(); DEBUG(0);
  _jStr_fwmJi_H80->clear(); DEBUG(0);
  _jStr_fwmJi_Q10->clear(); DEBUG(0);
  _jStr_fwmJi_Q20->clear(); DEBUG(0);
  _jStr_fwmJi_Q30->clear(); DEBUG(0);
  _jStr_fwmJi_Q40->clear(); DEBUG(0);
  _jStr_fwmJi_Q50->clear(); DEBUG(0);
  _jStr_fwmJi_Q60->clear(); DEBUG(0);
  _jStr_fwmJi_Q70->clear(); DEBUG(0);
  _jStr_fwmJi_Q80->clear(); DEBUG(0);
  _jStr_authorE->clear(); DEBUG(0);
  _jStr_authorJ->clear(); DEBUG(0);
  _jStr_radialParamE->clear(); DEBUG(0);
  _jStr_radialParamJ->clear(); DEBUG(0);
  _jStr_algE->clear(); DEBUG(0);
  _jStr_algJ->clear(); DEBUG(0);
  _jStr_inputE->clear(); DEBUG(0);
  _jStr_inputJ->clear(); DEBUG(0);
  _jStr_bjetE->clear(); DEBUG(0);
  _jStr_bjetJ->clear(); DEBUG(0);
  _jStr_xjetE->clear(); DEBUG(0);
  _jStr_xjetJ->clear(); DEBUG(0);
  _jStr_myJetsE->clear(); DEBUG(0);
  _jStr_myJetsJ->clear(); DEBUG(0);
  _jStr_indexJ->clear(); DEBUG(0);
  _jStr_jetMiscJ_numConstituents->clear(); DEBUG(0);
  _jStr_jetMiscJ_jetM->clear(); DEBUG(0);
  _jStr_jetMiscJ_jetMt->clear(); DEBUG(0);
  _jStr_jetMiscJ_jetE->clear(); DEBUG(0);
  _jStr_jetMiscJ_jetP->clear(); DEBUG(0);
  _jStr_jetMiscJ_jetEt->clear(); DEBUG(0);
  _jStr_jetMiscJ_jetPt->clear(); DEBUG(0);
  _jStr_jetMiscJ_jetPhi->clear(); DEBUG(0);
  _jStr_jetMiscJ_jetEta->clear(); DEBUG(0);
  _jStr_jetMiscJ_jetRapidity->clear(); DEBUG(0);
  _jStr_jetMiscJ_xJ->clear(); DEBUG(0);
  _jStr_jetMiscJ_gamma->clear(); DEBUG(0);
  _jStr_jetMiscJ_R->clear(); DEBUG(0);
  _jStr_jetMiscJ_Cbar->clear(); DEBUG(0);
  _jStr_jetMiscJ_numSubjets->clear(); DEBUG(0);
  _jStr_pullJ_det->clear(); DEBUG(0);
  _jStr_pullJ_ratio->clear(); DEBUG(0);
  _jStr_pullJ_pullPf->clear(); DEBUG(0);
  _jStr_pullJ_angularEccentricity->clear(); DEBUG(0);
  _jStr_pullJ_orientation->clear(); DEBUG(0);
  _jStr_pullJ_girth->clear(); DEBUG(0);
  _jStr_pullJ_Cbar->clear(); DEBUG(0);
  _jStr_pullJ_g->clear(); DEBUG(0);
  _jStr_pullJ_e->clear(); DEBUG(0);
  _jStr_pullJ_B->clear(); DEBUG(0);
  _jStr_pullJ_logB->clear(); DEBUG(0);
  _jStr_pullJ_pullTheta->clear(); DEBUG(0);
  _jStr_pullJ_pullMag->clear(); DEBUG(0);
  _jStr_ubhJ_z->clear(); DEBUG(0);
  _jStr_ubhJ_z2->clear(); DEBUG(0);
  _jStr_ubhJ_a1->clear(); DEBUG(0);
  _jStr_ubhJ_a2->clear(); DEBUG(0);
  _jStr_ubhJ_a3->clear(); DEBUG(0);
  _jStr_ubhJ_meanpt->clear(); DEBUG(0);
  _jStr_ubhJ_meanet->clear(); DEBUG(0);
  _jStr_ubhJ_mbar->clear(); DEBUG(0);
  _jStr_ubhJ_massDemocracy->clear(); DEBUG(0);
  _jStr_ubhJ_fE1->clear(); DEBUG(0);
  _jStr_ubhJ_fE2->clear(); DEBUG(0);
  _jStr_ubhJ_fE3->clear(); DEBUG(0);
  _jStr_ubhJ_fET1->clear(); DEBUG(0);
  _jStr_ubhJ_fET2->clear(); DEBUG(0);
  _jStr_ubhJ_fET3->clear(); DEBUG(0);
  _jStr_ubhJ_Alpha->clear(); DEBUG(0);
  _jStr_ubhJ_AlphaT->clear(); DEBUG(0);
  _jStr_ubhJ_betaflow->clear(); DEBUG(0);
  _jStr_ubhJ_betaflow_1GeV->clear(); DEBUG(0);
  _jStr_ubhJ_betaflow_5GeV->clear(); DEBUG(0);
  _jStr_ubhJ_y23->clear(); DEBUG(0);
  _jStr_ubhJ_y23_1GeV->clear(); DEBUG(0);
  _jStr_ubhJ_y23_5GeV->clear(); DEBUG(0);
  _jStr_ubhJ_lny23->clear(); DEBUG(0);
  _jStr_ubhJ_lny23_1GeV->clear(); DEBUG(0);
  _jStr_ubhJ_lny23_5GeV->clear(); DEBUG(0);
  _jStr_ubhJ_subjetAsymmetry->clear(); DEBUG(0);
  _jStr_dipJ_dipolarity->clear(); DEBUG(0);
  _jStr_nsubjnessJ_tau1->clear(); DEBUG(0);
  _jStr_nsubjnessJ_tau2->clear(); DEBUG(0);
  _jStr_nsubjnessJ_tau3->clear(); DEBUG(0);
  _jStr_nsubjnessJ_tau2tau1->clear(); DEBUG(0);
  _jStr_nsubjnessJ_tau3tau2->clear(); DEBUG(0);
  _jStr_psiJ2_psi->clear(); DEBUG(0);
  _jStr_psiJ2_rho->clear(); DEBUG(0);
  _jStr_psiRatiosJ_psi1->clear(); DEBUG(0);
  _jStr_psiRatiosJ_psi2->clear(); DEBUG(0);
  _jStr_psiRatiosJ_psi3->clear(); DEBUG(0);
  _jStr_psiRatiosJ_psi7->clear(); DEBUG(0);
  _jStr_psiRatiosJ_psi717->clear(); DEBUG(0);
  _jStr_psiRatiosJ_psi127->clear(); DEBUG(0);
  _jStr_psiRatiosJ_psi37->clear(); DEBUG(0);
  _jStr_bcfJ_bcfVersion_v0a1->clear(); DEBUG(0);
  _jStr_bcfJ_a_v0a1->clear(); DEBUG(0);
  _jStr_bcfJ_bcfT_v0a1->clear(); DEBUG(0);
  _jStr_bcfJ_bcf_v0a1->clear(); DEBUG(0);
  _jStr_bcfJ_bcfAsymY_v0a1->clear(); DEBUG(0);
  _jStr_bcfJ_bcfAsymPhi_v0a1->clear(); DEBUG(0);
  _jStr_bcfJ_bcfAsymYPhi_v0a1->clear(); DEBUG(0);
  _jStr_bcfJ_bcfAsymYPhi2_v0a1->clear(); DEBUG(0);
  _jStr_bcfJ_bcfVersion_v0a2->clear(); DEBUG(0);
  _jStr_bcfJ_a_v0a2->clear(); DEBUG(0);
  _jStr_bcfJ_bcfT_v0a2->clear(); DEBUG(0);
  _jStr_bcfJ_bcf_v0a2->clear(); DEBUG(0);
  _jStr_bcfJ_bcfAsymY_v0a2->clear(); DEBUG(0);
  _jStr_bcfJ_bcfAsymPhi_v0a2->clear(); DEBUG(0);
  _jStr_bcfJ_bcfAsymYPhi_v0a2->clear(); DEBUG(0);
  _jStr_bcfJ_bcfAsymYPhi2_v0a2->clear(); DEBUG(0);
  _jStr_pfJ_pf->clear(); DEBUG(0);
  _jStr_pfJ_detST->clear(); DEBUG(0);
  _jStr_pfJ_lambdaST->clear(); DEBUG(0);
  _jStr_zminJ->clear(); DEBUG(0);
  _jStr_zmaxJ->clear(); DEBUG(0);
  _jStr_tauJ09_a->clear(); DEBUG(0);
  _jStr_tauJ09_tau->clear(); DEBUG(0);
  _jStr_tauJ20_a->clear(); DEBUG(0);
  _jStr_tauJ20_tau->clear(); DEBUG(0);
  _jStr_tauJ40_a->clear(); DEBUG(0);
  _jStr_tauJ40_tau->clear(); DEBUG(0);
  _jStr_radEi_bcfJ1_v0a1->clear(); DEBUG(0);
  _jStr_radEi_bcfJ2_v0a1->clear(); DEBUG(0);
  _jStr_radEi_bcfJ3_v0a1->clear(); DEBUG(0);
  _jStr_radEi_bcfTJ1_v0a1->clear(); DEBUG(0);
  _jStr_radEi_bcfTJ2_v0a1->clear(); DEBUG(0);
  _jStr_radEi_bcfTJ3_v0a1->clear(); DEBUG(0);
  _jStr_radEi_bcfAsymYJ1_v0a1->clear(); DEBUG(0);
  _jStr_radEi_bcfAsymYJ2_v0a1->clear(); DEBUG(0);
  _jStr_radEi_bcfAsymYJ3_v0a1->clear(); DEBUG(0);
  _jStr_radEi_bcfAsymPhiJ1_v0a1->clear(); DEBUG(0);
  _jStr_radEi_bcfAsymPhiJ2_v0a1->clear(); DEBUG(0);
  _jStr_radEi_bcfAsymPhiJ3_v0a1->clear(); DEBUG(0);
  _jStr_radEi_bcfAsymYPhiJ1_v0a1->clear(); DEBUG(0);
  _jStr_radEi_bcfAsymYPhiJ2_v0a1->clear(); DEBUG(0);
  _jStr_radEi_bcfAsymYPhiJ3_v0a1->clear(); DEBUG(0);
  _jStr_radEi_bcfAsymYPhi2J1_v0a1->clear(); DEBUG(0);
  _jStr_radEi_bcfAsymYPhi2J2_v0a1->clear(); DEBUG(0);
  _jStr_radEi_bcfAsymYPhi2J3_v0a1->clear(); DEBUG(0);
  _jStr_radEi_bcfJ1_v0a2->clear(); DEBUG(0);
  _jStr_radEi_bcfJ2_v0a2->clear(); DEBUG(0);
  _jStr_radEi_bcfJ3_v0a2->clear(); DEBUG(0);
  _jStr_radEi_bcfTJ1_v0a2->clear(); DEBUG(0);
  _jStr_radEi_bcfTJ2_v0a2->clear(); DEBUG(0);
  _jStr_radEi_bcfTJ3_v0a2->clear(); DEBUG(0);
  _jStr_radEi_bcfAsymYJ1_v0a2->clear(); DEBUG(0);
  _jStr_radEi_bcfAsymYJ2_v0a2->clear(); DEBUG(0);
  _jStr_radEi_bcfAsymYJ3_v0a2->clear(); DEBUG(0);
  _jStr_radEi_bcfAsymPhiJ1_v0a2->clear(); DEBUG(0);
  _jStr_radEi_bcfAsymPhiJ2_v0a2->clear(); DEBUG(0);
  _jStr_radEi_bcfAsymPhiJ3_v0a2->clear(); DEBUG(0);
  _jStr_radEi_bcfAsymYPhiJ1_v0a2->clear(); DEBUG(0);
  _jStr_radEi_bcfAsymYPhiJ2_v0a2->clear(); DEBUG(0);
  _jStr_radEi_bcfAsymYPhiJ3_v0a2->clear(); DEBUG(0);
  _jStr_radEi_bcfAsymYPhi2J1_v0a2->clear(); DEBUG(0);
  _jStr_radEi_bcfAsymYPhi2J2_v0a2->clear(); DEBUG(0);
  _jStr_radEi_bcfAsymYPhi2J3_v0a2->clear(); DEBUG(0);
  _jStr_radEi_BJ1->clear(); DEBUG(0);
  _jStr_radEi_BJ2->clear(); DEBUG(0);
  _jStr_radEi_BJ3->clear(); DEBUG(0);
  _jStr_radEi_girthJ1->clear(); DEBUG(0);
  _jStr_radEi_girthJ2->clear(); DEBUG(0);
  _jStr_radEi_girthJ3->clear(); DEBUG(0);
  _jStr_radEi_girth32->clear(); DEBUG(0);
  _jStr_radEi_girth21->clear(); DEBUG(0);
  _jStr_radEi_girthAsymJ1J2->clear(); DEBUG(0);
  _jStr_radEi_alpha1->clear(); DEBUG(0);
  _jStr_radEi_alpha2->clear(); DEBUG(0);
  _jStr_radEi_alpha->clear(); DEBUG(0);
  _jStr_radEi_beta1->clear(); DEBUG(0);
  _jStr_radEi_beta2->clear(); DEBUG(0);
  _jStr_radEi_beta->clear(); DEBUG(0);
  _jStr_radEi_thetaJ1J2->clear(); DEBUG(0);
  _jStr_radEi_dipolarityInfLine->clear(); DEBUG(0);
  _jStr_radEi_dipolarityLineSeg->clear(); DEBUG(0);
  _jStr_radEi_BCF1->clear(); DEBUG(0);
  _jStr_radEi_BCF2->clear(); DEBUG(0);
  _jStr_radEi_BCF3->clear(); DEBUG(0);
  _jStr_radEi_dipolarity->clear(); DEBUG(0);
  _jStr_sphEJ_detSphericity->clear(); DEBUG(0);
  _jStr_sphEJ_detSpherocity->clear(); DEBUG(0);
  _jStr_sphEJ_sphericityLambda1->clear(); DEBUG(0);
  _jStr_sphEJ_sphericityLambda2->clear(); DEBUG(0);
  _jStr_sphEJ_sphericityLambda3->clear(); DEBUG(0);
  _jStr_sphEJ_spherocityLambda1->clear(); DEBUG(0);
  _jStr_sphEJ_spherocityLambda2->clear(); DEBUG(0);
  _jStr_sphEJ_spherocityLambda3->clear(); DEBUG(0);
  _jStr_sphEJ_circularity->clear(); DEBUG(0);
  _jStr_sphEJ_sphericity->clear(); DEBUG(0);
  _jStr_sphEJ_spherocity->clear(); DEBUG(0);
  _jStr_sphEJ_aplanarity->clear(); DEBUG(0);
  _jStr_sphEJ_aplanority->clear(); DEBUG(0);
  _jStr_sphEJ_Y->clear(); DEBUG(0);
  _jStr_sphEJ_planarity->clear(); DEBUG(0);
  _jStr_sphEJ_planority->clear(); DEBUG(0);
  _jStr_sphEJ_Dshape->clear(); DEBUG(0);
  _jStr_sphEJ_Cshape->clear(); DEBUG(0);
  _jStr_sphEJ_H2->clear(); DEBUG(0);
  _jStr_sphEJ_Fshape->clear(); DEBUG(0);
  _jStr_sphEJ_beamThrust->clear(); DEBUG(0);
  _jStr_sphEJ_G->clear(); DEBUG(0);
  _jStr_sphEJ_ST2D->clear(); DEBUG(0);
  _jStr_sphEJ_detMlin->clear(); DEBUG(0);
  _jStr_sphEJ_pft->clear(); DEBUG(0);
  _jStr_sphEi_detSphericity->clear(); DEBUG(0);
  _jStr_sphEi_detSpherocity->clear(); DEBUG(0);
  _jStr_sphEi_sphericityLambda1->clear(); DEBUG(0);
  _jStr_sphEi_sphericityLambda2->clear(); DEBUG(0);
  _jStr_sphEi_sphericityLambda3->clear(); DEBUG(0);
  _jStr_sphEi_spherocityLambda1->clear(); DEBUG(0);
  _jStr_sphEi_spherocityLambda2->clear(); DEBUG(0);
  _jStr_sphEi_spherocityLambda3->clear(); DEBUG(0);
  _jStr_sphEi_circularity->clear(); DEBUG(0);
  _jStr_sphEi_sphericity->clear(); DEBUG(0);
  _jStr_sphEi_spherocity->clear(); DEBUG(0);
  _jStr_sphEi_aplanarity->clear(); DEBUG(0);
  _jStr_sphEi_aplanority->clear(); DEBUG(0);
  _jStr_sphEi_Y->clear(); DEBUG(0);
  _jStr_sphEi_planarity->clear(); DEBUG(0);
  _jStr_sphEi_planority->clear(); DEBUG(0);
  _jStr_sphEi_Dshape->clear(); DEBUG(0);
  _jStr_sphEi_Cshape->clear(); DEBUG(0);
  _jStr_sphEi_H2->clear(); DEBUG(0);
  _jStr_sphEi_Fshape->clear(); DEBUG(0);
  _jStr_sphEi_beamThrust->clear(); DEBUG(0);
  _jStr_sphEi_G->clear(); DEBUG(0);
  _jStr_sphEi_ST2D->clear(); DEBUG(0);
  _jStr_sphEi_detMlin->clear(); DEBUG(0);
  _jStr_sphEi_pft->clear(); DEBUG(0);
  _jStr_sphJj_detSphericity->clear(); DEBUG(0);
  _jStr_sphJj_detSpherocity->clear(); DEBUG(0);
  _jStr_sphJj_sphericityLambda1->clear(); DEBUG(0);
  _jStr_sphJj_sphericityLambda2->clear(); DEBUG(0);
  _jStr_sphJj_sphericityLambda3->clear(); DEBUG(0);
  _jStr_sphJj_spherocityLambda1->clear(); DEBUG(0);
  _jStr_sphJj_spherocityLambda2->clear(); DEBUG(0);
  _jStr_sphJj_spherocityLambda3->clear(); DEBUG(0);
  _jStr_sphJj_circularity->clear(); DEBUG(0);
  _jStr_sphJj_sphericity->clear(); DEBUG(0);
  _jStr_sphJj_spherocity->clear(); DEBUG(0);
  _jStr_sphJj_aplanarity->clear(); DEBUG(0);
  _jStr_sphJj_aplanority->clear(); DEBUG(0);
  _jStr_sphJj_Y->clear(); DEBUG(0);
  _jStr_sphJj_planarity->clear(); DEBUG(0);
  _jStr_sphJj_planority->clear(); DEBUG(0);
  _jStr_sphJj_Dshape->clear(); DEBUG(0);
  _jStr_sphJj_Cshape->clear(); DEBUG(0);
  _jStr_sphJj_H2->clear(); DEBUG(0);
  _jStr_sphJj_Fshape->clear(); DEBUG(0);
  _jStr_sphJj_beamThrust->clear(); DEBUG(0);
  _jStr_sphJj_G->clear(); DEBUG(0);
  _jStr_sphJj_ST2D->clear(); DEBUG(0);
  _jStr_sphJj_detMlin->clear(); DEBUG(0);
  _jStr_sphJj_pft->clear(); DEBUG(0);
  _jStr_sphJi_detSphericity->clear(); DEBUG(0);
  _jStr_sphJi_detSpherocity->clear(); DEBUG(0);
  _jStr_sphJi_sphericityLambda1->clear(); DEBUG(0);
  _jStr_sphJi_sphericityLambda2->clear(); DEBUG(0);
  _jStr_sphJi_sphericityLambda3->clear(); DEBUG(0);
  _jStr_sphJi_spherocityLambda1->clear(); DEBUG(0);
  _jStr_sphJi_spherocityLambda2->clear(); DEBUG(0);
  _jStr_sphJi_spherocityLambda3->clear(); DEBUG(0);
  _jStr_sphJi_circularity->clear(); DEBUG(0);
  _jStr_sphJi_sphericity->clear(); DEBUG(0);
  _jStr_sphJi_spherocity->clear(); DEBUG(0);
  _jStr_sphJi_aplanarity->clear(); DEBUG(0);
  _jStr_sphJi_aplanority->clear(); DEBUG(0);
  _jStr_sphJi_Y->clear(); DEBUG(0);
  _jStr_sphJi_planarity->clear(); DEBUG(0);
  _jStr_sphJi_planority->clear(); DEBUG(0);
  _jStr_sphJi_Dshape->clear(); DEBUG(0);
  _jStr_sphJi_Cshape->clear(); DEBUG(0);
  _jStr_sphJi_H2->clear(); DEBUG(0);
  _jStr_sphJi_Fshape->clear(); DEBUG(0);
  _jStr_sphJi_beamThrust->clear(); DEBUG(0);
  _jStr_sphJi_G->clear(); DEBUG(0);
  _jStr_sphJi_ST2D->clear(); DEBUG(0);
  _jStr_sphJi_detMlin->clear(); DEBUG(0);
  _jStr_sphJi_pft->clear(); DEBUG(0);

#ifdef OUTPUT_CONSTITUENTS
  _nConst->clear(); DEBUG(0);
  _const_et->clear(); DEBUG(0);
  _const_rapidity->clear(); DEBUG(0);
  _const_phi->clear(); DEBUG(0);
  _const_Ret->clear(); DEBUG(0);
  _const_Drapidity->clear(); DEBUG(0);
  _const_Dphi->clear(); DEBUG(0);
#endif

  _YFlip12->clear(); DEBUG(0);
  _YFlip23->clear(); DEBUG(0);
  _EM_FRACTION_CLUSTER->clear(); DEBUG(0);
  _ELLIPTICAREA->clear(); DEBUG(0);
  _AMAREA->clear(); DEBUG(0);
  _HULL_LENGTH->clear(); DEBUG(0);
  _HULL_AREA->clear(); DEBUG(0);
  _EM_FRACTION_MCTRUTH->clear(); DEBUG(0);
  _LowEtConstituentsFrac->clear(); DEBUG(0);
  _JetEccentricity->clear(); DEBUG(0);
  _DRtoReco->clear(); DEBUG(0);
  _PtNearest->clear(); DEBUG(0);
  _WIDTH->clear(); DEBUG(0);
  _rbb->clear(); DEBUG(0);
  _rfilt->clear(); DEBUG(0);

  _jvf->clear(); DEBUG(0);
  _ntrk->clear(); DEBUG(0);
  _trkpt->clear(); DEBUG(0);

  _QGtag->clear(); DEBUG(0);
  _QGtag2->clear(); DEBUG(0);
  _jetTrueFlavour->clear(); DEBUG(0);
  w_cmb->clear(); DEBUG(0);
  w_TrackCounting2D->clear(); DEBUG(0);
  w_JetProb->clear(); DEBUG(0);
  w_IP1D->clear(); DEBUG(0);
  w_IP2D->clear(); DEBUG(0);
  w_IP3D->clear(); DEBUG(0);
  w_SV0->clear(); DEBUG(0);
  w_SV1->clear(); DEBUG(0);
  w_SV2->clear(); DEBUG(0);
  w_BaselineTagger->clear(); DEBUG(0);
  w_JetFitterTag->clear(); DEBUG(0);
  w_JetFitterCOMB->clear(); DEBUG(0);
  w_JetFitterTagNN->clear(); DEBUG(0);
  w_JetFitterCOMBNN->clear(); DEBUG(0);
  w_SoftMuonTag->clear(); DEBUG(0);
  w_SoftElectronTag->clear(); DEBUG(0);

  return StatusCode::SUCCESS;
}

StatusCode MyAnalysis::execute() {
  // fedisableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
  // feclearexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
  // CALLGRIND_START_INSTRUMENTATION;
  m_eventNr++; 

  ATH_MSG_DEBUG(" in execute()");
  ATH_MSG_DEBUG("Events processed: " << m_eventNr);

  StatusCode sc;

  // initialize first before processing each event
  sc = initEvent();
  if ( sc.isFailure() ) {
    ATH_MSG_WARNING("initEvent failed. Continue");
  }

  sc = addEventInfo();
  if (sc.isFailure() ) {
    ATH_MSG_WARNING("Failure in getEventInfo() ");
    return StatusCode::FAILURE;
  }

  StatusCode m_sc=service("StoreGateSvc", m_storeGate);

  const TruthParticleContainer* mcpartTES = 0;
  sc = m_storeGate->retrieve(mcpartTES, m_truthParticleContainerName);
  if(sc.isFailure() || !mcpartTES) {
    ATH_MSG_WARNING("No AOD MC truth particle container found in TDS (maybe this is real data?)");
    // return StatusCode::SUCCESS;
  }
  else {
    ATH_MSG_DEBUG("MC Truth Container Successfully Retrieved");
    _nparticles = mcpartTES->size();
    ATH_MSG_INFO("Number of TruthParticles: " << _nparticles);
  }

  // retrieve VxFinder container from TES should be here
  const VxContainer* container(0);
  sc = m_storeGate->retrieve(container, m_vxCandidatesName);
  if(sc.isFailure() || container == 0) {
    ATH_MSG_WARNING("No VxCandidate container " << m_vxCandidatesName << " found in StoreGate!");
    return sc;
  }

  try { // Pokemon exception handling

    _jStr_numPVE = container->size(); DEBUG(_jStr_numPVE);
    if((*container)[0]) {
      _jStr_PVxE = (*container)[0]->recVertex().position().x(); DEBUG(_jStr_PVxE);
      _jStr_PVyE = (*container)[0]->recVertex().position().y(); DEBUG(_jStr_PVyE);
      _jStr_PVzE = (*container)[0]->recVertex().position().z(); DEBUG(_jStr_PVzE);
      _jStr_PVrE = sqrt((_jStr_PVxE * _jStr_PVxE) + (_jStr_PVyE * _jStr_PVyE)); DEBUG(_jStr_PVrE);
      _jStr_PVchiSqrE = (*container)[0]->recVertex().fitQuality().chiSquared(); DEBUG(_jStr_PVchiSqrE);
      _jStr_PVnDoFE = (*container)[0]->recVertex().fitQuality().doubleNumberDoF(); DEBUG(_jStr_PVnDoFE);
      _jStr_PVfitE = fabs(Divide(_jStr_PVchiSqrE, _jStr_PVnDoFE)); DEBUG(_jStr_PVfitE);
      _jStr_PVnTrkE = (*container)[0]->vxTrackAtVertex()->size(); DEBUG(_jStr_PVnTrkE);
      _jStr_PVtypeE = (*container)[0]->vertexType(); DEBUG(_jStr_PVtypeE);
    }
  
    ATH_MSG_INFO("Number of primary vertices: " << _jStr_numPVE);

    ATH_MSG_INFO("Get JetCollections");

    for(unsigned int i = 0; i < m_jetcollnames.size(); ++i) { // loop over jet collection names

      std::string jetcollname = m_jetcollnames[i]; DEBUG(jetcollname);
      const JetCollection* jetColl = 0;

      sc = m_storeGate->retrieve(jetColl, jetcollname); DEBUG(0); // jetcollname is a std::string containing the SG key
      if(sc.isFailure() || !jetColl) {
	ATH_MSG_WARNING("No Jet Collection " << jetcollname);
	return StatusCode::FAILURE; DEBUG(0);
      }

      ATH_MSG_INFO("Retrieved JetCollection " << jetcollname << " with " << jetColl->size() << " jets");

      JetCollection::const_iterator jetItr  = jetColl->begin(); DEBUG(0);
      JetCollection::const_iterator jetItrE = jetColl->end(); DEBUG(0);

      FilterBjets bjetFilter;
      if(Bjet()(jetColl)) { // bjets
	if(Input()(jetColl) == 5) { // == Truth
	  int numBjets = bjetFilter(jetColl);
	  ATH_MSG_INFO("JetCollection is truth with " << numBjets << " jets");
	}
	else {
	  int numBjets = bjetFilter(jetColl, m_cutIP3DSV1);
	  ATH_MSG_INFO("JetCollection is not truth with " << numBjets << " jets");
	}
	ATH_MSG_INFO("After filtering, JetCollection " << jetcollname << " has " << jetColl->size() << " b-jets");
      }

      FilterJets jetFilter;
      if(Xjet()(jetColl)) { // Jets to be filtered
	int numXjets = jetFilter(jetColl, 2.8, 80.*GeV);
	ATH_MSG_INFO("After eta/pT filtering, JetCollection " << jetcollname << " has " << jetColl->size() << " jets");
      }
    
      _jStr_numJets = jetColl->size(); DEBUG(_jStr_numJets);
    
      JetCollectionHelper::jetcollection_t theJets(jetColl->begin(), jetColl->end());
      JetCollectionHelper::sort_jets(theJets, JetSorters::sortJetByEtDown());
      JetCollectionHelper::jetcollection_t::iterator firstJet = theJets.begin();
      JetCollectionHelper::jetcollection_t::iterator lastJet  = theJets.end();

      unsigned int index = 1; DEBUG(index);
      for(; firstJet != lastJet; ++firstJet) { DEBUG(*firstJet); // loop over jets
	if((*firstJet) == 0) continue;
	const Jet* aJet = (*firstJet)->clone(true, true); DEBUG(aJet);
	if(aJet == 0) continue;

	_jStr_numPVJ->push_back(_jStr_numPVE);
	_jStr_PVxJ->push_back(_jStr_PVxE);
	_jStr_PVyJ->push_back(_jStr_PVyE);
	_jStr_PVzJ->push_back(_jStr_PVzE);
	_jStr_PVrJ->push_back(_jStr_PVrE);
	_jStr_PVchiSqrJ->push_back(_jStr_PVchiSqrE);
	_jStr_PVnDoFJ->push_back(_jStr_PVnDoFE);
	_jStr_PVfitJ->push_back(_jStr_PVfitE);
	_jStr_PVnTrkJ->push_back(_jStr_PVnTrkE);
	_jStr_PVtypeJ->push_back(_jStr_PVtypeE);
      
	// --- number of tags attached to the jet:
	int ntag = aJet->jetTagInfoVector().size(); DEBUG(ntag);
	// --- get btagging weights:
	w_cmb->push_back(aJet->getFlavourTagWeight()); DEBUG(0); // combination of IP3D and SV1, the recommanded tagger
      
	w_TrackCounting2D->push_back(aJet->getFlavourTagWeight("TrackCounting2D")); DEBUG(0);
	w_JetProb->push_back(aJet->getFlavourTagWeight("JetProb")); DEBUG(0);
	w_IP1D->push_back(aJet->getFlavourTagWeight("IP1D")); DEBUG(0);
	w_IP2D->push_back(aJet->getFlavourTagWeight("IP2D")); DEBUG(0);
	w_IP3D->push_back(aJet->getFlavourTagWeight("IP3D")); DEBUG(0);
	w_SV0->push_back(aJet->getFlavourTagWeight("SV0")); DEBUG(0);
	w_SV1->push_back(aJet->getFlavourTagWeight("SV1")); DEBUG(0);
	w_SV2->push_back(aJet->getFlavourTagWeight("SV2")); DEBUG(0);
	w_BaselineTagger->push_back(aJet->getFlavourTagWeight("BaselineTagger")); DEBUG(0);
	w_JetFitterTag->push_back(aJet->getFlavourTagWeight("JetFitterTag")); DEBUG(0);
	w_JetFitterCOMB->push_back(aJet->getFlavourTagWeight("JetFitterCOMB")); DEBUG(0);
	w_JetFitterTagNN->push_back(aJet->getFlavourTagWeight("JetFitterTagNN")); DEBUG(0);
	w_JetFitterCOMBNN->push_back(aJet->getFlavourTagWeight("JetFitterCOMBNN")); DEBUG(0);
	w_SoftMuonTag->push_back(aJet->getFlavourTagWeight("SoftMuonTag")); DEBUG(0);
	w_SoftElectronTag->push_back(aJet->getFlavourTagWeight("SoftElectronTag")); DEBUG(0);
      
	std::string jetname = JetAuthor()(aJet); DEBUG(jetname);

	// --- get the true label of the jet from MC Truth:
	std::string label("N/A"); DEBUG(label);
	const Analysis::TruthInfo* mcinfo = aJet->tagInfo<Analysis::TruthInfo>("TruthInfo"); DEBUG(0);
	if(mcinfo) {
	  label = mcinfo->jetTruthLabel(); DEBUG(label);
	}
	else {
	  ATH_MSG_WARNING("Could not find TruthInfo for jet " << jetname);
	}
	int iflav = 0;
	if(label=="B") iflav = 5; // b-jet
	if(label=="C") iflav = 4; // c-jet
	if(label=="T") iflav = 15; // tau
    
	ATH_MSG_INFO("Jet " << jetname << " with " << ntag
		     << " tags, IP3D+SV1 weight: "<< aJet->getFlavourTagWeight()
		     << ", truth label: " << iflav); DEBUG(0);
      
	_jetTrueFlavour->push_back(iflav); DEBUG(iflav);

	//_QGtag->push_back(QGtag(mcpartTES, aJet));
	_QGtag2->push_back(QGtag2(mcpartTES, aJet));

	// AngularStructure as(aJet); DEBUG(0);
	// AngularCorrelation acR2(aJet, 0.2); DEBUG(0);

	JetMisc jetMiscJ(aJet); DEBUG(0);
	BCF bcfJ01(aJet, bcfVersion, 1); DEBUG(0);
	BCF bcfJ02(aJet, bcfVersion, 2); DEBUG(0);
	JetDipolarity dipJ(aJet); DEBUG(0);
	EventShape evtJj(jetMiscJ.subjets(), aJet); DEBUG(0);
	FoxWolfram fwmJi(jetMiscJ.subjets(), aJet, true, maxTermsFW); DEBUG(0);
	FoxWolfram fwmJj(jetMiscJ.subjets(), aJet, false, maxTermsFW); DEBUG(0);
	Nsubjettiness nsubjnessJ(aJet); DEBUG(0);
	PlanarFlow pfJ(aJet); DEBUG(0);
	PsiRatios psiRatiosJ(aJet); DEBUG(0);
	Psi psiJ2(aJet, 0.2); DEBUG(0);
	Pull pullJ(aJet); DEBUG(0);
	SphericitySpherocity sphJi(jetMiscJ.subjets(), aJet, true); DEBUG(0);
	SphericitySpherocity sphJj(jetMiscJ.subjets(), aJet); DEBUG(0);
	Tau tauJ09(-0.9, aJet); DEBUG(0);
	Tau tauJ20(-2.0, aJet); DEBUG(0);
	Tau tauJ40(-4.0, aJet); DEBUG(0);
	UnburiedHiggs ubhJ(aJet); DEBUG(0);
	ZminJ zminJ; DEBUG(0);
	ZmaxJ zmaxJ; DEBUG(0);

	// Get the JVF, number of tracks, and sum track pT of this jet with respect to each primary vertex (vectors of doubles or ints) 
	std::vector<double> jvf; DEBUG(0);
	std::vector<int> ntrk_vector; DEBUG(0);
	std::vector<double> trkpt; DEBUG(0);

	double jvf0 = NaN; DEBUG(0);
	int ntrk0 = -1; DEBUG(0);
	double trkpt0 = NaN; DEBUG(0);
      
	CHECK(m_jvfTool->getJVF(aJet, jvf, ntrk_vector, trkpt)); DEBUG(0);
      
	ATH_MSG_INFO("Jet author: " << jetname
		     << ", JVF size: " << jvf.size()
		     << ", ntrk_vector size: " << ntrk_vector.size()
		     << ", trkpt size: " << trkpt.size()); DEBUG(0);
      
	if(!jvf.empty()) {
	  jvf0 = jvf[0];
	  ntrk0 = ntrk_vector[0];
	  trkpt0 = trkpt[0];
	  ATH_MSG_INFO("Jet author: " << jetname
		       << ", JVF0 = " << jvf0
		       << ", ntrk0 = " << ntrk0
		       << ", trkpt0 = " << trkpt0); DEBUG(0);
	}
	_jvf->push_back(jvf0); DEBUG(jvf0);
	_ntrk->push_back(ntrk0); DEBUG(ntrk0);
	_trkpt->push_back(trkpt0); DEBUG(trkpt0);
      
	std::vector<std::string> keys = aJet->getMomentKeys(); DEBUG(0);
	std::vector<double> values; DEBUG(0);
	int nKeys = keys.size(); DEBUG(nKeys);

	// "If the moment does not exist, the value 0 will be returned"
	for(int i = 0; i != nKeys; i++) {
	  values.push_back(aJet->getMoment(keys.at(i))); DEBUG(keys.at(i));
	}
      
	std::tr1::unordered_map<std::string, double> map; DEBUG(0);
	std::transform(keys.begin(), keys.end(),
		       values.begin(),
		       std::inserter(map, map.begin()),
		       std::make_pair<std::string,double>); DEBUG(0);

	// Now we replace all those unphysical values with NaN
	_YFlip12->push_back(map.count("YFlip12")? map.find("YFlip12")->second: NaN); DEBUG(0);
	_YFlip23->push_back(map.count("YFlip23")? map.find("YFlip23")->second: NaN); DEBUG(0);
	_EM_FRACTION_CLUSTER->push_back(map.count("EM_FRACTION_CLUSTER")? map.find("EM_FRACTION_CLUSTER")->second: NaN); DEBUG(0);
	_ELLIPTICAREA->push_back(map.count("ELLIPTICAREA")? map.find("ELLIPTICAREA")->second: NaN); DEBUG(0);
	_AMAREA->push_back(map.count("AMAREA")? map.find("AMAREA")->second: NaN); DEBUG(0);
	_HULL_LENGTH->push_back(map.count("HULL_LENGHT")? map.find("HULL_LENGHT")->second: NaN); DEBUG(0);
	_HULL_AREA->push_back(map.count("HULL_AREA")? map.find("HULL_AREA")->second: NaN); DEBUG(0);
	_EM_FRACTION_MCTRUTH->push_back(map.count("EM_FRACTION_MCTRUTH")? map.find("EM_FRACTION_MCTRUTH")->second: NaN); DEBUG(0);
	_LowEtConstituentsFrac->push_back(map.count("LowEtConstituentsFrac")? map.find("LowEtConstituentsFrac")->second: NaN); DEBUG(0);
	_JetEccentricity->push_back(map.count("JetEccentricity")? map.find("JetEccentricity")->second: NaN); DEBUG(0);
	_DRtoReco->push_back(map.count("DRtoReco")? map.find("DRtoReco")->second: NaN); DEBUG(0);
	_PtNearest->push_back(map.count("PtNearest")? map.find("PtNearest")->second: NaN); DEBUG(0);
	_WIDTH->push_back(map.count("WIDTH")? map.find("WIDTH")->second: NaN); DEBUG(0);
	_rbb->push_back(map.count("CORE_RBB")? map.find("CORE_RBB")->second: NaN); DEBUG(0);
	_rfilt->push_back(map.count("CORE_RFILT")? map.find("CORE_RFILT")->second: NaN); DEBUG(0);

	JetConstituentIterator firstConstituent = JetConstituentIterator::first(aJet);
	JetConstituentIterator lastConstituent = JetConstituentIterator::last(aJet);

#ifdef OUTPUT_CONSTITUENTS
	std::vector<double> et; DEBUG(0);
	et.clear(); DEBUG(0);
	std::vector<double> phi; DEBUG(0);
	phi.clear(); DEBUG(0);
	std::vector<double> rapidity; DEBUG(0);
	rapidity.clear(); DEBUG(0);
	std::vector<double> Ret; DEBUG(0);
	Ret.clear(); DEBUG(0);
	std::vector<double> Drapidity; DEBUG(0);
	Drapidity.clear(); DEBUG(0);
	std::vector<double> Dphi; DEBUG(0);
	Dphi.clear(); DEBUG(0);
#endif
	DEBUG(0);
	for(; firstConstituent != lastConstituent; ++firstConstituent) {
	  pullJ(firstConstituent);
	  psiJ2(firstConstituent);
	  psiRatiosJ(firstConstituent);
	  pfJ(firstConstituent);
	  tauJ09(firstConstituent);
	  tauJ20(firstConstituent);
	  tauJ40(firstConstituent);
	  bcfJ01(firstConstituent);
	  bcfJ02(firstConstituent);
	  nsubjnessJ(firstConstituent);
	  dipJ(firstConstituent);

#ifdef OUTPUT_CONSTITUENTS
	  et.push_back(firstConstituent.et()); DEBUG(0);
	  rapidity.push_back(firstConstituent.rapidity()); DEBUG(0);
	  phi.push_back(firstConstituent.phi()); DEBUG(0);
	  Ret.push_back(Divide(firstConstituent.et(), aJet->et())); DEBUG(0);
	  Drapidity.push_back(aJet->rapidity() - firstConstituent.rapidity()); DEBUG(0);
	  Dphi.push_back(DeltaPhi()(aJet->phi(), firstConstituent.phi())); DEBUG(0);
#endif
	}
	DEBUG(0);

#ifdef OUTPUT_CONSTITUENTS
	_const_et->push_back(et); DEBUG(0);
	_const_rapidity->push_back(rapidity); DEBUG(0);
	_const_phi->push_back(phi); DEBUG(0);
      
	_const_Ret->push_back(Ret); DEBUG(0);
	_const_Drapidity->push_back(Drapidity); DEBUG(0);
	_const_Dphi->push_back(Dphi); DEBUG(0);
#endif      

	_jStr_authorJ->push_back(JetAuthor()(aJet)); DEBUG(0);
	_jStr_radialParamJ->push_back(RadialParameter()(aJet)); DEBUG(0);
	_jStr_algJ->push_back(JetAlg()(aJet)); DEBUG(0);
	_jStr_inputJ->push_back(Input()(aJet)); DEBUG(0);
	_jStr_bjetJ->push_back(Bjet()(aJet)); DEBUG(0);
	_jStr_xjetJ->push_back(Xjet()(aJet)); DEBUG(0);
	_jStr_myJetsJ->push_back(MyJets()(aJet)); DEBUG(0);
	_jStr_indexJ->push_back(index++); DEBUG(0);
      
	_jStr_evtJj_etaJ1->push_back(evtJj.etaJ1()); DEBUG(0);
	_jStr_evtJj_etaJ2->push_back(evtJj.etaJ2()); DEBUG(0);
	_jStr_evtJj_etaJ3->push_back(evtJj.etaJ3()); DEBUG(0);
	_jStr_evtJj_etaJ4->push_back(evtJj.etaJ4()); DEBUG(0);
	_jStr_evtJj_pTJ1->push_back(evtJj.pTJ1()); DEBUG(0);
	_jStr_evtJj_pTJ2->push_back(evtJj.pTJ2()); DEBUG(0);
	_jStr_evtJj_pTJ3->push_back(evtJj.pTJ3()); DEBUG(0);
	_jStr_evtJj_pTJ4->push_back(evtJj.pTJ4()); DEBUG(0);

	_jStr_evtJj_fE1->push_back(evtJj.fE1()); DEBUG(0);
	_jStr_evtJj_fE2->push_back(evtJj.fE2()); DEBUG(0);
	_jStr_evtJj_fE3->push_back(evtJj.fE3()); DEBUG(0);
	_jStr_evtJj_fET1->push_back(evtJj.fET1()); DEBUG(0);
	_jStr_evtJj_fET2->push_back(evtJj.fET2()); DEBUG(0);
	_jStr_evtJj_fET3->push_back(evtJj.fET3()); DEBUG(0);
	_jStr_evtJj_mct->push_back(evtJj.mct()); DEBUG(0);
	_jStr_evtJj_q->push_back(evtJj.q()); DEBUG(0);
	_jStr_evtJj_mjj->push_back(evtJj.mjj()); DEBUG(0);
	_jStr_evtJj_mTjj->push_back(evtJj.mTjj()); DEBUG(0);
	_jStr_evtJj_mjjj->push_back(evtJj.mjjj()); DEBUG(0);
	_jStr_evtJj_mTjjj->push_back(evtJj.mTjjj()); DEBUG(0);
	_jStr_evtJj_mjjjj->push_back(evtJj.mjjjj()); DEBUG(0);
	_jStr_evtJj_mTjjjj->push_back(evtJj.mTjjjj()); DEBUG(0);
	_jStr_evtJj_zetaPlus->push_back(evtJj.zetaPlus()); DEBUG(0);
	_jStr_evtJj_zetaMinus->push_back(evtJj.zetaMinus()); DEBUG(0);
	_jStr_evtJj_B12->push_back(evtJj.B12()); DEBUG(0);
	_jStr_evtJj_BT12->push_back(evtJj.BT12()); DEBUG(0);
	_jStr_evtJj_Alpha->push_back(evtJj.Alpha()); DEBUG(0);
	_jStr_evtJj_AlphaT->push_back(evtJj.AlphaT()); DEBUG(0);
	_jStr_evtJj_massDemocracy->push_back(evtJj.massDemocracy()); DEBUG(0);
	_jStr_evtJj_betaflow->push_back(evtJj.betaflow()); DEBUG(0);
	_jStr_evtJj_betaflow_1GeV->push_back(evtJj.betaflow(1.*GeV)); DEBUG(0);
	_jStr_evtJj_betaflow_5GeV->push_back(evtJj.betaflow(5.*GeV)); DEBUG(0);
	_jStr_evtJj_y23->push_back(evtJj.y23()); DEBUG(0);
	_jStr_evtJj_y23_1GeV->push_back(evtJj.y23(1.*GeV)); DEBUG(0);
	_jStr_evtJj_y23_5GeV->push_back(evtJj.y23(5.*GeV)); DEBUG(0);
	_jStr_evtJj_lny23->push_back(evtJj.lny23()); DEBUG(0);
	_jStr_evtJj_lny23_1GeV->push_back(evtJj.lny23(1.*GeV)); DEBUG(0);
	_jStr_evtJj_lny23_5GeV->push_back(evtJj.lny23(5.*GeV)); DEBUG(0);
	_jStr_evtJj_theta->push_back(evtJj.theta()); DEBUG(0);
	_jStr_evtJj_asym->push_back(evtJj.asym()); DEBUG(0);
	_jStr_evtJj_yB->push_back(evtJj.yB()); DEBUG(0);
	_jStr_evtJj_yStar->push_back(evtJj.yStar()); DEBUG(0);
	_jStr_evtJj_thetaStar->push_back(evtJj.thetaStar()); DEBUG(0);
	_jStr_evtJj_chi->push_back(evtJj.chi()); DEBUG(0);
	_jStr_evtJj_deltaPhiJJ->push_back(evtJj.deltaPhiJJ()); DEBUG(0);
	_jStr_evtJj_deltaThetaJJ->push_back(evtJj.deltaThetaJJ()); DEBUG(0);
	_jStr_evtJj_deltaEtaJJ->push_back(evtJj.deltaEtaJJ()); DEBUG(0);
	_jStr_evtJj_deltaRapidityJJ->push_back(evtJj.deltaRapidityJJ()); DEBUG(0);
	_jStr_evtJj_deltaRJJ->push_back(evtJj.deltaRJJ()); DEBUG(0);
	_jStr_evtJj_deltaRJJY->push_back(evtJj.deltaRJJY()); DEBUG(0);
	_jStr_evtJj_sigmaPhiJJ->push_back(evtJj.sigmaPhiJJ()); DEBUG(0);
	_jStr_evtJj_sigmaThetaJJ->push_back(evtJj.sigmaThetaJJ()); DEBUG(0);
	_jStr_evtJj_sigmaEtaJJ->push_back(evtJj.sigmaEtaJJ()); DEBUG(0);
	_jStr_evtJj_sigmaRapidityJJ->push_back(evtJj.sigmaRapidityJJ()); DEBUG(0);
	_jStr_evtJj_sigmaPtJJ->push_back(evtJj.sigmaPtJJ()); DEBUG(0);
	_jStr_evtJj_sigmaEtJJ->push_back(evtJj.sigmaEtJJ()); DEBUG(0);
	_jStr_evtJj_sigmaEt12->push_back(evtJj.sigmaEt12()); DEBUG(0);
	_jStr_evtJj_sigmaEt34->push_back(evtJj.sigmaEt34()); DEBUG(0);
	_jStr_evtJj_A234->push_back(evtJj.A234()); DEBUG(0);
	_jStr_evtJj_asymPhiJJ->push_back(evtJj.asymPhiJJ()); DEBUG(0);
	_jStr_evtJj_asymThetaJJ->push_back(evtJj.asymThetaJJ()); DEBUG(0);
	_jStr_evtJj_asymEtaJJ->push_back(evtJj.asymEtaJJ()); DEBUG(0);
	_jStr_evtJj_asymRapidityJJ->push_back(evtJj.asymRapidityJJ()); DEBUG(0);
	_jStr_evtJj_acoplanarity->push_back(evtJj.acoplanarity()); DEBUG(0);
	_jStr_evtJj_twist->push_back(evtJj.twist()); DEBUG(0);
	_jStr_evtJj_twistY->push_back(evtJj.twistY()); DEBUG(0);
	_jStr_evtJj_jetSumE->push_back(evtJj.jetSumE()); DEBUG(0);
	_jStr_evtJj_jetSumET->push_back(evtJj.jetSumET()); DEBUG(0);
	_jStr_evtJj_jetSumPT->push_back(evtJj.jetSumPT()); DEBUG(0);
	_jStr_evtJj_jetSumM->push_back(evtJj.jetSumM()); DEBUG(0);
	_jStr_evtJj_jetSumMT->push_back(evtJj.jetSumMT()); DEBUG(0);
	_jStr_evtJj_HTprime->push_back(evtJj.HTprime()); DEBUG(0);
	_jStr_evtJj_centrality->push_back(evtJj.centrality()); DEBUG(0);
	_jStr_evtJj_centralityP->push_back(evtJj.centralityP()); DEBUG(0);
	_jStr_evtJj_zminJ1J2->push_back(evtJj.zminJ1J2()); DEBUG(0);
	_jStr_evtJj_zmaxJ1J2->push_back(evtJj.zmaxJ1J2()); DEBUG(0);
	_jStr_evtJj_zminAllJets->push_back(evtJj.zminAllJets()); DEBUG(0);
	_jStr_evtJj_zmaxAllJets->push_back(evtJj.zmaxAllJets()); DEBUG(0);
	_jStr_evtJj_zminJ1J2Phi->push_back(evtJj.zminJ1J2Phi()); DEBUG(0);
	_jStr_evtJj_zminJ1J2Theta->push_back(evtJj.zminJ1J2Theta()); DEBUG(0);
	_jStr_evtJj_zminJ1J2Eta->push_back(evtJj.zminJ1J2Eta()); DEBUG(0);
	_jStr_evtJj_zminJ1J2Rapidity->push_back(evtJj.zminJ1J2Rapidity()); DEBUG(0);
	_jStr_evtJj_cosHelicityJ1->push_back(evtJj.cosHelicityJ1()); DEBUG(0);
	_jStr_evtJj_helicityJ1->push_back(evtJj.helicityJ1()); DEBUG(0);
	_jStr_evtJj_azilicityJ1->push_back(evtJj.azilicityJ1()); DEBUG(0);
	_jStr_evtJj_cosHelicityJ2->push_back(evtJj.cosHelicityJ2()); DEBUG(0);
	_jStr_evtJj_helicityJ2->push_back(evtJj.helicityJ2()); DEBUG(0);
	_jStr_evtJj_azilicityJ2->push_back(evtJj.azilicityJ2()); DEBUG(0);
	_jStr_evtJj_cosThetaJ1->push_back(evtJj.cosThetaJ1()); DEBUG("the penultimate line");
	_jStr_evtJj_cosThetaJ2->push_back(evtJj.cosThetaJ2()); DEBUG("the last line");
	_jStr_evtJj_deltaRapidityXtoJ1CM->push_back(evtJj.deltaRapidityXtoJ1CM()); DEBUG(0);
	_jStr_evtJj_deltaRapidityXtoJ2CM->push_back(evtJj.deltaRapidityXtoJ2CM()); DEBUG(0);
	_jStr_evtJj_deltaRapidityXtoJ1->push_back(evtJj.deltaRapidityXtoJ1()); DEBUG(0);
	_jStr_evtJj_deltaRapidityXtoJ2->push_back(evtJj.deltaRapidityXtoJ2()); DEBUG(0);
	_jStr_evtJj_cosTheta1->push_back(evtJj.cosTheta1()); DEBUG(0);
	_jStr_evtJj_cosTheta2->push_back(evtJj.cosTheta2()); DEBUG(0);
	_jStr_evtJj_cosThetaStar1->push_back(evtJj.cosThetaStar1()); DEBUG(0);
	_jStr_evtJj_cosThetaStar2->push_back(evtJj.cosThetaStar2()); DEBUG(0);
	_jStr_evtJj_cosPhiTilde1->push_back(evtJj.cosPhiTilde1()); DEBUG(0);
	_jStr_evtJj_cosPhiTilde2->push_back(evtJj.cosPhiTilde2()); DEBUG(0);
	_jStr_evtJj_cosPhi->push_back(evtJj.cosPhi()); DEBUG(0);
	_jStr_evtJj_PhiTilde1->push_back(evtJj.PhiTilde1()); DEBUG(0);
	_jStr_evtJj_PhiTilde2->push_back(evtJj.PhiTilde2()); DEBUG(0);
	_jStr_evtJj_Phi->push_back(evtJj.Phi()); DEBUG(0);
	_jStr_evtJj_M->push_back(evtJj.M()); DEBUG(0);
	_jStr_evtJj_DeltaM->push_back(evtJj.DeltaM()); DEBUG(0);
	_jStr_evtJj_asymM->push_back(evtJj.asymM()); DEBUG(0);
	_jStr_evtJj_Q->push_back(evtJj.Q()); DEBUG(0);
	_jStr_evtJj_SCDF->push_back(evtJj.SCDF()); DEBUG(0);
	_jStr_evtJj_dPhiIJCDF->push_back(evtJj.dPhiIJCDF()); DEBUG(0);
	_jStr_evtJj_dPhiKLCDF->push_back(evtJj.dPhiKLCDF()); DEBUG(0);
	_jStr_evtJj_PtIPtJCDF->push_back(evtJj.PtIPtJCDF()); DEBUG(0);
	_jStr_evtJj_PtKPtLCDF->push_back(evtJj.PtKPtLCDF()); DEBUG(0);
	_jStr_evtJj_SD0->push_back(evtJj.SD0()); DEBUG(0);
	_jStr_evtJj_dPhiIJD0->push_back(evtJj.dPhiIJD0()); DEBUG(0);
	_jStr_evtJj_dPhiKLD0->push_back(evtJj.dPhiKLD0()); DEBUG(0);
	_jStr_evtJj_PtIPtJD0->push_back(evtJj.PtIPtJD0()); DEBUG(0);
	_jStr_evtJj_PtKPtLD0->push_back(evtJj.PtKPtLD0()); DEBUG(0);
	_jStr_evtJj_DeltaSCDF->push_back(evtJj.DeltaSCDF()); DEBUG(0);
	_jStr_evtJj_DeltaSD0->push_back(evtJj.DeltaSD0()); DEBUG(0);
	_jStr_evtJj_Delta12->push_back(evtJj.Delta12()); DEBUG(0);
	_jStr_evtJj_Delta12n->push_back(evtJj.Delta12n()); DEBUG(0);
	_jStr_evtJj_M14->push_back(evtJj.M14()); DEBUG(0);
	_jStr_evtJj_y1y2->push_back(evtJj.y1y2()); DEBUG(0);

	_jStr_evtJj_Mprime4->push_back(evtJj.Mprime4()); DEBUG(0);
	_jStr_evtJj_dMprime4->push_back(evtJj.dMprime4()); DEBUG(0);
	_jStr_evtJj_MprimeAvg4->push_back(evtJj.MprimeAvg4()); DEBUG(0);
	_jStr_evtJj_DeltaMin4->push_back(evtJj.DeltaMin4()); DEBUG(0);
	_jStr_evtJj_DeltaMax4->push_back(evtJj.DeltaMax4()); DEBUG(0);
	_jStr_evtJj_DeltaPhiXX4->push_back(evtJj.DeltaPhiXX4()); DEBUG(0);
	_jStr_evtJj_DeltaYXX4->push_back(evtJj.DeltaYXX4()); DEBUG(0);
	_jStr_evtJj_DeltaRXX4->push_back(evtJj.DeltaRXX4()); DEBUG(0);
	_jStr_evtJj_TwistXX4->push_back(evtJj.TwistXX4()); DEBUG(0);
	_jStr_evtJj_separated4->push_back(evtJj.separated4()); DEBUG(0);
	_jStr_evtJj_Mprime->push_back(evtJj.Mprime()); DEBUG(0);
	_jStr_evtJj_dMprime->push_back(evtJj.dMprime()); DEBUG(0);
	_jStr_evtJj_MprimeAvg->push_back(evtJj.MprimeAvg()); DEBUG(0);
	_jStr_evtJj_DeltaMin->push_back(evtJj.DeltaMin()); DEBUG(0);
	_jStr_evtJj_DeltaMax->push_back(evtJj.DeltaMax()); DEBUG(0);
	_jStr_evtJj_DeltaPhiXX->push_back(evtJj.DeltaPhiXX()); DEBUG(0);
	_jStr_evtJj_DeltaYXX->push_back(evtJj.DeltaYXX()); DEBUG(0);
	_jStr_evtJj_DeltaRXX->push_back(evtJj.DeltaRXX()); DEBUG(0);
	_jStr_evtJj_TwistXX->push_back(evtJj.TwistXX()); DEBUG(0);
	_jStr_evtJj_separated->push_back(evtJj.separated()); DEBUG(0);

	_jStr_fwmJj_xiPlus->push_back(fwmJj.xiPlus()); DEBUG(0);
	_jStr_fwmJj_xiMinus->push_back(fwmJj.xiMinus()); DEBUG(0);
	_jStr_fwmJj_xPlus->push_back(fwmJj.xPlus()); DEBUG(0);
	_jStr_fwmJj_xMinus->push_back(fwmJj.xMinus()); DEBUG(0);
	_jStr_fwmJj_Psi1->push_back(fwmJj.Psi1()); DEBUG(0);
	_jStr_fwmJj_B0->push_back(fwmJj.B0()); DEBUG(0);
	_jStr_fwmJj_B1->push_back(fwmJj.B1()); DEBUG(0);
	_jStr_fwmJj_B2->push_back(fwmJj.B2()); DEBUG(0);
	_jStr_fwmJj_B3->push_back(fwmJj.B3()); DEBUG(0);
	_jStr_fwmJj_B4->push_back(fwmJj.B4()); DEBUG(0);
	_jStr_fwmJj_B5->push_back(fwmJj.B5()); DEBUG(0);
	_jStr_fwmJj_B6->push_back(fwmJj.B6()); DEBUG(0);
	_jStr_fwmJj_B7->push_back(fwmJj.B7()); DEBUG(0);
	_jStr_fwmJj_B8->push_back(fwmJj.B8()); DEBUG(0);
	_jStr_fwmJj_C0->push_back(fwmJj.C0()); DEBUG(0);
	_jStr_fwmJj_C1->push_back(fwmJj.C1()); DEBUG(0);
	_jStr_fwmJj_C2->push_back(fwmJj.C2()); DEBUG(0);
	_jStr_fwmJj_C3->push_back(fwmJj.C3()); DEBUG(0);
	_jStr_fwmJj_C4->push_back(fwmJj.C4()); DEBUG(0);
	_jStr_fwmJj_C5->push_back(fwmJj.C5()); DEBUG(0);
	_jStr_fwmJj_C6->push_back(fwmJj.C6()); DEBUG(0);
	_jStr_fwmJj_C7->push_back(fwmJj.C7()); DEBUG(0);
	_jStr_fwmJj_C8->push_back(fwmJj.C8()); DEBUG(0);
	_jStr_fwmJj_K0->push_back(fwmJj.K0()); DEBUG(0);
	_jStr_fwmJj_K1->push_back(fwmJj.K1()); DEBUG(0);
	_jStr_fwmJj_K2->push_back(fwmJj.K2()); DEBUG(0);
	_jStr_fwmJj_K3->push_back(fwmJj.K3()); DEBUG(0);
	_jStr_fwmJj_K4->push_back(fwmJj.K4()); DEBUG(0);
	_jStr_fwmJj_K5->push_back(fwmJj.K5()); DEBUG(0);
	_jStr_fwmJj_K6->push_back(fwmJj.K6()); DEBUG(0);
	_jStr_fwmJj_K7->push_back(fwmJj.K7()); DEBUG(0);
	_jStr_fwmJj_K8->push_back(fwmJj.K8()); DEBUG(0);
	_jStr_fwmJj_D0->push_back(fwmJj.D0()); DEBUG(0);
	_jStr_fwmJj_D1->push_back(fwmJj.D1()); DEBUG(0);
	_jStr_fwmJj_D2->push_back(fwmJj.D2()); DEBUG(0);
	_jStr_fwmJj_D3->push_back(fwmJj.D3()); DEBUG(0);
	_jStr_fwmJj_D4->push_back(fwmJj.D4()); DEBUG(0);
	_jStr_fwmJj_D5->push_back(fwmJj.D5()); DEBUG(0);
	_jStr_fwmJj_D6->push_back(fwmJj.D6()); DEBUG(0);
	_jStr_fwmJj_D7->push_back(fwmJj.D7()); DEBUG(0);
	_jStr_fwmJj_D8->push_back(fwmJj.D8()); DEBUG(0);
	_jStr_fwmJj_H0->push_back(fwmJj.H0()); DEBUG(0);
	_jStr_fwmJj_H1->push_back(fwmJj.H1()); DEBUG(0);
	_jStr_fwmJj_H2->push_back(fwmJj.H2()); DEBUG(0);
	_jStr_fwmJj_H3->push_back(fwmJj.H3()); DEBUG(0);
	_jStr_fwmJj_H4->push_back(fwmJj.H4()); DEBUG(0);
	_jStr_fwmJj_H5->push_back(fwmJj.H5()); DEBUG(0);
	_jStr_fwmJj_H6->push_back(fwmJj.H6()); DEBUG(0);
	_jStr_fwmJj_H7->push_back(fwmJj.H7()); DEBUG(0);
	_jStr_fwmJj_H8->push_back(fwmJj.H8()); DEBUG(0);
	_jStr_fwmJj_Q0->push_back(fwmJj.Q0()); DEBUG(0);
	_jStr_fwmJj_Q1->push_back(fwmJj.Q1()); DEBUG(0);
	_jStr_fwmJj_Q2->push_back(fwmJj.Q2()); DEBUG(0);
	_jStr_fwmJj_Q3->push_back(fwmJj.Q3()); DEBUG(0);
	_jStr_fwmJj_Q4->push_back(fwmJj.Q4()); DEBUG(0);
	_jStr_fwmJj_Q5->push_back(fwmJj.Q5()); DEBUG(0);
	_jStr_fwmJj_Q6->push_back(fwmJj.Q6()); DEBUG(0);
	_jStr_fwmJj_Q7->push_back(fwmJj.Q7()); DEBUG(0);
	_jStr_fwmJj_Q8->push_back(fwmJj.Q8()); DEBUG(0);
	_jStr_fwmJj_Pi1->push_back(fwmJj.Pi1()); DEBUG(0);
	_jStr_fwmJj_Pi2->push_back(fwmJj.Pi2()); DEBUG(0);
	_jStr_fwmJj_Pi3->push_back(fwmJj.Pi3()); DEBUG(0);
	_jStr_fwmJj_Pi4->push_back(fwmJj.Pi4()); DEBUG(0);
	_jStr_fwmJj_B10->push_back(fwmJj.B10()); DEBUG(0);
	_jStr_fwmJj_B20->push_back(fwmJj.B20()); DEBUG(0);
	_jStr_fwmJj_B30->push_back(fwmJj.B30()); DEBUG(0);
	_jStr_fwmJj_B40->push_back(fwmJj.B40()); DEBUG(0);
	_jStr_fwmJj_B50->push_back(fwmJj.B50()); DEBUG(0);
	_jStr_fwmJj_B60->push_back(fwmJj.B60()); DEBUG(0);
	_jStr_fwmJj_B70->push_back(fwmJj.B70()); DEBUG(0);
	_jStr_fwmJj_B80->push_back(fwmJj.B80()); DEBUG(0);
	_jStr_fwmJj_C10->push_back(fwmJj.C10()); DEBUG(0);
	_jStr_fwmJj_C20->push_back(fwmJj.C20()); DEBUG(0);
	_jStr_fwmJj_C30->push_back(fwmJj.C30()); DEBUG(0);
	_jStr_fwmJj_C40->push_back(fwmJj.C40()); DEBUG(0);
	_jStr_fwmJj_C50->push_back(fwmJj.C50()); DEBUG(0);
	_jStr_fwmJj_C60->push_back(fwmJj.C60()); DEBUG(0);
	_jStr_fwmJj_C70->push_back(fwmJj.C70()); DEBUG(0);
	_jStr_fwmJj_C80->push_back(fwmJj.C80()); DEBUG(0);
	_jStr_fwmJj_K10->push_back(fwmJj.K10()); DEBUG(0);
	_jStr_fwmJj_K20->push_back(fwmJj.K20()); DEBUG(0);
	_jStr_fwmJj_K30->push_back(fwmJj.K30()); DEBUG(0);
	_jStr_fwmJj_K40->push_back(fwmJj.K40()); DEBUG(0);
	_jStr_fwmJj_K50->push_back(fwmJj.K50()); DEBUG(0);
	_jStr_fwmJj_K60->push_back(fwmJj.K60()); DEBUG(0);
	_jStr_fwmJj_K70->push_back(fwmJj.K70()); DEBUG(0);
	_jStr_fwmJj_K80->push_back(fwmJj.K80()); DEBUG(0);
	_jStr_fwmJj_D10->push_back(fwmJj.D10()); DEBUG(0);
	_jStr_fwmJj_D20->push_back(fwmJj.D20()); DEBUG(0);
	_jStr_fwmJj_D30->push_back(fwmJj.D30()); DEBUG(0);
	_jStr_fwmJj_D40->push_back(fwmJj.D40()); DEBUG(0);
	_jStr_fwmJj_D50->push_back(fwmJj.D50()); DEBUG(0);
	_jStr_fwmJj_D60->push_back(fwmJj.D60()); DEBUG(0);
	_jStr_fwmJj_D70->push_back(fwmJj.D70()); DEBUG(0);
	_jStr_fwmJj_D80->push_back(fwmJj.D80()); DEBUG(0);
	_jStr_fwmJj_H10->push_back(fwmJj.H10()); DEBUG(0);
	_jStr_fwmJj_H20->push_back(fwmJj.H20()); DEBUG(0);
	_jStr_fwmJj_H30->push_back(fwmJj.H30()); DEBUG(0);
	_jStr_fwmJj_H40->push_back(fwmJj.H40()); DEBUG(0);
	_jStr_fwmJj_H50->push_back(fwmJj.H50()); DEBUG(0);
	_jStr_fwmJj_H60->push_back(fwmJj.H60()); DEBUG(0);
	_jStr_fwmJj_H70->push_back(fwmJj.H70()); DEBUG(0);
	_jStr_fwmJj_H80->push_back(fwmJj.H80()); DEBUG(0);
	_jStr_fwmJj_Q10->push_back(fwmJj.Q10()); DEBUG(0);
	_jStr_fwmJj_Q20->push_back(fwmJj.Q20()); DEBUG(0);
	_jStr_fwmJj_Q30->push_back(fwmJj.Q30()); DEBUG(0);
	_jStr_fwmJj_Q40->push_back(fwmJj.Q40()); DEBUG(0);
	_jStr_fwmJj_Q50->push_back(fwmJj.Q50()); DEBUG(0);
	_jStr_fwmJj_Q60->push_back(fwmJj.Q60()); DEBUG(0);
	_jStr_fwmJj_Q70->push_back(fwmJj.Q70()); DEBUG(0);
	_jStr_fwmJj_Q80->push_back(fwmJj.Q80()); DEBUG(0);
	_jStr_fwmJi_xiPlus->push_back(fwmJi.xiPlus()); DEBUG(0);
	_jStr_fwmJi_xiMinus->push_back(fwmJi.xiMinus()); DEBUG(0);
	_jStr_fwmJi_xPlus->push_back(fwmJi.xPlus()); DEBUG(0);
	_jStr_fwmJi_xMinus->push_back(fwmJi.xMinus()); DEBUG(0);
	_jStr_fwmJi_Psi1->push_back(fwmJi.Psi1()); DEBUG(0);
	_jStr_fwmJi_B0->push_back(fwmJi.B0()); DEBUG(0);
	_jStr_fwmJi_B1->push_back(fwmJi.B1()); DEBUG(0);
	_jStr_fwmJi_B2->push_back(fwmJi.B2()); DEBUG(0);
	_jStr_fwmJi_B3->push_back(fwmJi.B3()); DEBUG(0);
	_jStr_fwmJi_B4->push_back(fwmJi.B4()); DEBUG(0);
	_jStr_fwmJi_B5->push_back(fwmJi.B5()); DEBUG(0);
	_jStr_fwmJi_B6->push_back(fwmJi.B6()); DEBUG(0);
	_jStr_fwmJi_B7->push_back(fwmJi.B7()); DEBUG(0);
	_jStr_fwmJi_B8->push_back(fwmJi.B8()); DEBUG(0);
	_jStr_fwmJi_C0->push_back(fwmJi.C0()); DEBUG(0);
	_jStr_fwmJi_C1->push_back(fwmJi.C1()); DEBUG(0);
	_jStr_fwmJi_C2->push_back(fwmJi.C2()); DEBUG(0);
	_jStr_fwmJi_C3->push_back(fwmJi.C3()); DEBUG(0);
	_jStr_fwmJi_C4->push_back(fwmJi.C4()); DEBUG(0);
	_jStr_fwmJi_C5->push_back(fwmJi.C5()); DEBUG(0);
	_jStr_fwmJi_C6->push_back(fwmJi.C6()); DEBUG(0);
	_jStr_fwmJi_C7->push_back(fwmJi.C7()); DEBUG(0);
	_jStr_fwmJi_C8->push_back(fwmJi.C8()); DEBUG(0);
	_jStr_fwmJi_K0->push_back(fwmJi.K0()); DEBUG(0);
	_jStr_fwmJi_K1->push_back(fwmJi.K1()); DEBUG(0);
	_jStr_fwmJi_K2->push_back(fwmJi.K2()); DEBUG(0);
	_jStr_fwmJi_K3->push_back(fwmJi.K3()); DEBUG(0);
	_jStr_fwmJi_K4->push_back(fwmJi.K4()); DEBUG(0);
	_jStr_fwmJi_K5->push_back(fwmJi.K5()); DEBUG(0);
	_jStr_fwmJi_K6->push_back(fwmJi.K6()); DEBUG(0);
	_jStr_fwmJi_K7->push_back(fwmJi.K7()); DEBUG(0);
	_jStr_fwmJi_K8->push_back(fwmJi.K8()); DEBUG(0);
	_jStr_fwmJi_D0->push_back(fwmJi.D0()); DEBUG(0);
	_jStr_fwmJi_D1->push_back(fwmJi.D1()); DEBUG(0);
	_jStr_fwmJi_D2->push_back(fwmJi.D2()); DEBUG(0);
	_jStr_fwmJi_D3->push_back(fwmJi.D3()); DEBUG(0);
	_jStr_fwmJi_D4->push_back(fwmJi.D4()); DEBUG(0);
	_jStr_fwmJi_D5->push_back(fwmJi.D5()); DEBUG(0);
	_jStr_fwmJi_D6->push_back(fwmJi.D6()); DEBUG(0);
	_jStr_fwmJi_D7->push_back(fwmJi.D7()); DEBUG(0);
	_jStr_fwmJi_D8->push_back(fwmJi.D8()); DEBUG(0);
	_jStr_fwmJi_H0->push_back(fwmJi.H0()); DEBUG(0);
	_jStr_fwmJi_H1->push_back(fwmJi.H1()); DEBUG(0);
	_jStr_fwmJi_H2->push_back(fwmJi.H2()); DEBUG(0);
	_jStr_fwmJi_H3->push_back(fwmJi.H3()); DEBUG(0);
	_jStr_fwmJi_H4->push_back(fwmJi.H4()); DEBUG(0);
	_jStr_fwmJi_H5->push_back(fwmJi.H5()); DEBUG(0);
	_jStr_fwmJi_H6->push_back(fwmJi.H6()); DEBUG(0);
	_jStr_fwmJi_H7->push_back(fwmJi.H7()); DEBUG(0);
	_jStr_fwmJi_H8->push_back(fwmJi.H8()); DEBUG(0);
	_jStr_fwmJi_Q0->push_back(fwmJi.Q0()); DEBUG(0);
	_jStr_fwmJi_Q1->push_back(fwmJi.Q1()); DEBUG(0);
	_jStr_fwmJi_Q2->push_back(fwmJi.Q2()); DEBUG(0);
	_jStr_fwmJi_Q3->push_back(fwmJi.Q3()); DEBUG(0);
	_jStr_fwmJi_Q4->push_back(fwmJi.Q4()); DEBUG(0);
	_jStr_fwmJi_Q5->push_back(fwmJi.Q5()); DEBUG(0);
	_jStr_fwmJi_Q6->push_back(fwmJi.Q6()); DEBUG(0);
	_jStr_fwmJi_Q7->push_back(fwmJi.Q7()); DEBUG(0);
	_jStr_fwmJi_Q8->push_back(fwmJi.Q8()); DEBUG(0);
	_jStr_fwmJi_Pi1->push_back(fwmJi.Pi1()); DEBUG(0);
	_jStr_fwmJi_Pi2->push_back(fwmJi.Pi2()); DEBUG(0);
	_jStr_fwmJi_Pi3->push_back(fwmJi.Pi3()); DEBUG(0);
	_jStr_fwmJi_Pi4->push_back(fwmJi.Pi4()); DEBUG(0);
	_jStr_fwmJi_B10->push_back(fwmJi.B10()); DEBUG(0);
	_jStr_fwmJi_B20->push_back(fwmJi.B20()); DEBUG(0);
	_jStr_fwmJi_B30->push_back(fwmJi.B30()); DEBUG(0);
	_jStr_fwmJi_B40->push_back(fwmJi.B40()); DEBUG(0);
	_jStr_fwmJi_B50->push_back(fwmJi.B50()); DEBUG(0);
	_jStr_fwmJi_B60->push_back(fwmJi.B60()); DEBUG(0);
	_jStr_fwmJi_B70->push_back(fwmJi.B70()); DEBUG(0);
	_jStr_fwmJi_B80->push_back(fwmJi.B80()); DEBUG(0);
	_jStr_fwmJi_C10->push_back(fwmJi.C10()); DEBUG(0);
	_jStr_fwmJi_C20->push_back(fwmJi.C20()); DEBUG(0);
	_jStr_fwmJi_C30->push_back(fwmJi.C30()); DEBUG(0);
	_jStr_fwmJi_C40->push_back(fwmJi.C40()); DEBUG(0);
	_jStr_fwmJi_C50->push_back(fwmJi.C50()); DEBUG(0);
	_jStr_fwmJi_C60->push_back(fwmJi.C60()); DEBUG(0);
	_jStr_fwmJi_C70->push_back(fwmJi.C70()); DEBUG(0);
	_jStr_fwmJi_C80->push_back(fwmJi.C80()); DEBUG(0);
	_jStr_fwmJi_K10->push_back(fwmJi.K10()); DEBUG(0);
	_jStr_fwmJi_K20->push_back(fwmJi.K20()); DEBUG(0);
	_jStr_fwmJi_K30->push_back(fwmJi.K30()); DEBUG(0);
	_jStr_fwmJi_K40->push_back(fwmJi.K40()); DEBUG(0);
	_jStr_fwmJi_K50->push_back(fwmJi.K50()); DEBUG(0);
	_jStr_fwmJi_K60->push_back(fwmJi.K60()); DEBUG(0);
	_jStr_fwmJi_K70->push_back(fwmJi.K70()); DEBUG(0);
	_jStr_fwmJi_K80->push_back(fwmJi.K80()); DEBUG(0);
	_jStr_fwmJi_D10->push_back(fwmJi.D10()); DEBUG(0);
	_jStr_fwmJi_D20->push_back(fwmJi.D20()); DEBUG(0);
	_jStr_fwmJi_D30->push_back(fwmJi.D30()); DEBUG(0);
	_jStr_fwmJi_D40->push_back(fwmJi.D40()); DEBUG(0);
	_jStr_fwmJi_D50->push_back(fwmJi.D50()); DEBUG(0);
	_jStr_fwmJi_D60->push_back(fwmJi.D60()); DEBUG(0);
	_jStr_fwmJi_D70->push_back(fwmJi.D70()); DEBUG(0);
	_jStr_fwmJi_D80->push_back(fwmJi.D80()); DEBUG(0);
	_jStr_fwmJi_H10->push_back(fwmJi.H10()); DEBUG(0);
	_jStr_fwmJi_H20->push_back(fwmJi.H20()); DEBUG(0);
	_jStr_fwmJi_H30->push_back(fwmJi.H30()); DEBUG(0);
	_jStr_fwmJi_H40->push_back(fwmJi.H40()); DEBUG(0);
	_jStr_fwmJi_H50->push_back(fwmJi.H50()); DEBUG(0);
	_jStr_fwmJi_H60->push_back(fwmJi.H60()); DEBUG(0);
	_jStr_fwmJi_H70->push_back(fwmJi.H70()); DEBUG(0);
	_jStr_fwmJi_H80->push_back(fwmJi.H80()); DEBUG(0);
	_jStr_fwmJi_Q10->push_back(fwmJi.Q10()); DEBUG(0);
	_jStr_fwmJi_Q20->push_back(fwmJi.Q20()); DEBUG(0);
	_jStr_fwmJi_Q30->push_back(fwmJi.Q30()); DEBUG(0);
	_jStr_fwmJi_Q40->push_back(fwmJi.Q40()); DEBUG(0);
	_jStr_fwmJi_Q50->push_back(fwmJi.Q50()); DEBUG(0);
	_jStr_fwmJi_Q60->push_back(fwmJi.Q60()); DEBUG(0);
	_jStr_fwmJi_Q70->push_back(fwmJi.Q70()); DEBUG(0);
	_jStr_fwmJi_Q80->push_back(fwmJi.Q80()); DEBUG(0);
      
	_jStr_jetMiscJ_numConstituents->push_back(jetMiscJ.numConstituents()); DEBUG(0);
	_jStr_jetMiscJ_jetM->push_back(jetMiscJ.jetM()); DEBUG(0);
	_jStr_jetMiscJ_jetMt->push_back(jetMiscJ.jetMt()); DEBUG(0);
	_jStr_jetMiscJ_jetE->push_back(jetMiscJ.jetE()); DEBUG(0);
	_jStr_jetMiscJ_jetP->push_back(jetMiscJ.jetP()); DEBUG(0);
	_jStr_jetMiscJ_jetEt->push_back(jetMiscJ.jetEt()); DEBUG(0);
	_jStr_jetMiscJ_jetPt->push_back(jetMiscJ.jetPt()); DEBUG(0);
	_jStr_jetMiscJ_jetPhi->push_back(jetMiscJ.jetPhi()); DEBUG(0);
	_jStr_jetMiscJ_jetEta->push_back(jetMiscJ.jetEta()); DEBUG(0);
	_jStr_jetMiscJ_jetRapidity->push_back(jetMiscJ.jetRapidity()); DEBUG(0);
	_jStr_jetMiscJ_xJ->push_back(jetMiscJ.xJ()); DEBUG(0);
	_jStr_jetMiscJ_gamma->push_back(jetMiscJ.gamma()); DEBUG(0);
	_jStr_jetMiscJ_R->push_back(jetMiscJ.R()); DEBUG(0);
	_jStr_jetMiscJ_Cbar->push_back(jetMiscJ.Cbar()); DEBUG(0);
	_jStr_jetMiscJ_numSubjets->push_back(jetMiscJ.numSubjets()); DEBUG(0);
	_jStr_pullJ_det->push_back(pullJ.det()); DEBUG(0);
	_jStr_pullJ_ratio->push_back(pullJ.ratio()); DEBUG(0);
	_jStr_pullJ_pullPf->push_back(pullJ.pullPf()); DEBUG(0);
	_jStr_pullJ_angularEccentricity->push_back(pullJ.angularEccentricity()); DEBUG(0);
	_jStr_pullJ_orientation->push_back(pullJ.orientation()); DEBUG(0);
	_jStr_pullJ_girth->push_back(pullJ.girth()); DEBUG(0);
	_jStr_pullJ_Cbar->push_back(pullJ.Cbar()); DEBUG(0);
	_jStr_pullJ_g->push_back(pullJ.g()); DEBUG(0);
	_jStr_pullJ_e->push_back(pullJ.e()); DEBUG(0);
	_jStr_pullJ_B->push_back(pullJ.B()); DEBUG(0);
	_jStr_pullJ_logB->push_back(pullJ.logB()); DEBUG(0);
	_jStr_pullJ_pullTheta->push_back(pullJ.pullTheta()); DEBUG(0);
	_jStr_pullJ_pullMag->push_back(pullJ.pullMag()); DEBUG(0);
	_jStr_ubhJ_z->push_back(ubhJ.z()); DEBUG(0);
	_jStr_ubhJ_z2->push_back(ubhJ.z2()); DEBUG(0);
	_jStr_ubhJ_a1->push_back(ubhJ.a1()); DEBUG(0);
	_jStr_ubhJ_a2->push_back(ubhJ.a2()); DEBUG(0);
	_jStr_ubhJ_a3->push_back(ubhJ.a3()); DEBUG(0);
	_jStr_ubhJ_meanpt->push_back(ubhJ.meanpt()); DEBUG(0);
	_jStr_ubhJ_meanet->push_back(ubhJ.meanet()); DEBUG(0);
	_jStr_ubhJ_mbar->push_back(ubhJ.mbar()); DEBUG(0);
	_jStr_ubhJ_massDemocracy->push_back(ubhJ.massDemocracy()); DEBUG(0);
	_jStr_ubhJ_fE1->push_back(ubhJ.fE1()); DEBUG(0);
	_jStr_ubhJ_fE2->push_back(ubhJ.fE2()); DEBUG(0);
	_jStr_ubhJ_fE3->push_back(ubhJ.fE3()); DEBUG(0);
	_jStr_ubhJ_fET1->push_back(ubhJ.fET1()); DEBUG(0);
	_jStr_ubhJ_fET2->push_back(ubhJ.fET2()); DEBUG(0);
	_jStr_ubhJ_fET3->push_back(ubhJ.fET3()); DEBUG(0);
	_jStr_ubhJ_Alpha->push_back(ubhJ.Alpha()); DEBUG(0);
	_jStr_ubhJ_AlphaT->push_back(ubhJ.AlphaT()); DEBUG(0);
	_jStr_ubhJ_betaflow->push_back(ubhJ.betaflow()); DEBUG(0);
	_jStr_ubhJ_betaflow_1GeV->push_back(ubhJ.betaflow(1.*GeV)); DEBUG(0);
	_jStr_ubhJ_betaflow_5GeV->push_back(ubhJ.betaflow(5.*GeV)); DEBUG(0);
	_jStr_ubhJ_y23->push_back(ubhJ.y23()); DEBUG(0);
	_jStr_ubhJ_y23_1GeV->push_back(ubhJ.y23(1.*GeV)); DEBUG(0);
	_jStr_ubhJ_y23_5GeV->push_back(ubhJ.y23(5.*GeV)); DEBUG(0);
	_jStr_ubhJ_lny23->push_back(ubhJ.lny23()); DEBUG(0);
	_jStr_ubhJ_lny23_1GeV->push_back(ubhJ.lny23(1.*GeV)); DEBUG(0);
	_jStr_ubhJ_lny23_5GeV->push_back(ubhJ.lny23(5.*GeV)); DEBUG(0);
	_jStr_ubhJ_subjetAsymmetry->push_back(ubhJ.subjetAsymmetry()); DEBUG(0);
	_jStr_dipJ_dipolarity->push_back(dipJ.dipolarity()); DEBUG(0);
	_jStr_nsubjnessJ_tau1->push_back(nsubjnessJ.tau1()); DEBUG(0);
	_jStr_nsubjnessJ_tau2->push_back(nsubjnessJ.tau2()); DEBUG(0);
	_jStr_nsubjnessJ_tau3->push_back(nsubjnessJ.tau3()); DEBUG(0);
	_jStr_nsubjnessJ_tau2tau1->push_back(nsubjnessJ.tau2tau1()); DEBUG(0);
	_jStr_nsubjnessJ_tau3tau2->push_back(nsubjnessJ.tau3tau2()); DEBUG(0);
	_jStr_psiJ2_psi->push_back(psiJ2.psi()); DEBUG(0);
	_jStr_psiJ2_rho->push_back(psiJ2.rho()); DEBUG(0);
	_jStr_psiRatiosJ_psi1->push_back(psiRatiosJ.psi1()); DEBUG(0);
	_jStr_psiRatiosJ_psi2->push_back(psiRatiosJ.psi2()); DEBUG(0);
	_jStr_psiRatiosJ_psi3->push_back(psiRatiosJ.psi3()); DEBUG(0);
	_jStr_psiRatiosJ_psi7->push_back(psiRatiosJ.psi7()); DEBUG(0);
	_jStr_psiRatiosJ_psi717->push_back(psiRatiosJ.psi717()); DEBUG(0);
	_jStr_psiRatiosJ_psi127->push_back(psiRatiosJ.psi127()); DEBUG(0);
	_jStr_psiRatiosJ_psi37->push_back(psiRatiosJ.psi37()); DEBUG(0);
	_jStr_bcfJ_bcfVersion_v0a1->push_back(bcfJ01.bcfVersion()); DEBUG(0);
	_jStr_bcfJ_a_v0a1->push_back(bcfJ01.a()); DEBUG(0);
	_jStr_bcfJ_bcfT_v0a1->push_back(bcfJ01.bcfT()); DEBUG(0);
	_jStr_bcfJ_bcf_v0a1->push_back(bcfJ01.bcf()); DEBUG(0);
	_jStr_bcfJ_bcfAsymY_v0a1->push_back(bcfJ01.bcfAsymY()); DEBUG(0);
	_jStr_bcfJ_bcfAsymPhi_v0a1->push_back(bcfJ01.bcfAsymPhi()); DEBUG(0);
	_jStr_bcfJ_bcfAsymYPhi_v0a1->push_back(bcfJ01.bcfAsymYPhi()); DEBUG(0);
	_jStr_bcfJ_bcfAsymYPhi2_v0a1->push_back(bcfJ01.bcfAsymYPhi2()); DEBUG(0);
	_jStr_bcfJ_bcfVersion_v0a2->push_back(bcfJ02.bcfVersion()); DEBUG(0);
	_jStr_bcfJ_a_v0a2->push_back(bcfJ02.a()); DEBUG(0);
	_jStr_bcfJ_bcfT_v0a2->push_back(bcfJ02.bcfT()); DEBUG(0);
	_jStr_bcfJ_bcf_v0a2->push_back(bcfJ02.bcf()); DEBUG(0);
	_jStr_bcfJ_bcfAsymY_v0a2->push_back(bcfJ02.bcfAsymY()); DEBUG(0);
	_jStr_bcfJ_bcfAsymPhi_v0a2->push_back(bcfJ02.bcfAsymPhi()); DEBUG(0);
	_jStr_bcfJ_bcfAsymYPhi_v0a2->push_back(bcfJ02.bcfAsymYPhi()); DEBUG(0);
	_jStr_bcfJ_bcfAsymYPhi2_v0a2->push_back(bcfJ02.bcfAsymYPhi2()); DEBUG(0);
	_jStr_pfJ_pf->push_back(pfJ.pf()); DEBUG(0);
	_jStr_pfJ_detST->push_back(pfJ.detST()); DEBUG(0);
	_jStr_pfJ_lambdaST->push_back(pfJ.lambdaST()); DEBUG(0);
	_jStr_zminJ->push_back(zminJ(aJet)); DEBUG(0);
	_jStr_zmaxJ->push_back(zmaxJ(aJet)); DEBUG(0);
	//std::cout << "zmaxJ: " << zmaxJ(aJet) << std::endl;
	//if(zmaxJ(aJet)>0.8) {
	std::cout << " zmaxJ: " << zmaxJ(aJet)
		  << " R: " << RadialParameter()(aJet) 
		  << " Alg: " << JetAlg()(aJet) 
		  << " Input: " << Input()(aJet)
		  << " Bjet: " << Bjet()(aJet)
		  << " gamma: " << jetMiscJ.gamma()
		  << " E: " << jetMiscJ.jetE()
		  << " P: " << jetMiscJ.jetP()
		  << " M: " << jetMiscJ.jetM(); _QGtag->push_back(QGtag(mcpartTES, aJet));
	//}
	_jStr_tauJ09_a->push_back(tauJ09.a()); DEBUG(0);
	_jStr_tauJ09_tau->push_back(tauJ09.tau()); DEBUG(0);
	_jStr_tauJ20_a->push_back(tauJ20.a()); DEBUG(0);
	_jStr_tauJ20_tau->push_back(tauJ20.tau()); DEBUG(0);
	_jStr_tauJ40_a->push_back(tauJ40.a()); DEBUG(0);
	_jStr_tauJ40_tau->push_back(tauJ40.tau()); DEBUG(0);
	_jStr_sphJj_detSphericity->push_back(sphJj.detSphericity()); DEBUG(0);
	_jStr_sphJj_detSpherocity->push_back(sphJj.detSpherocity()); DEBUG(0);
	_jStr_sphJj_sphericityLambda1->push_back(sphJj.sphericityLambda1()); DEBUG(0);
	_jStr_sphJj_sphericityLambda2->push_back(sphJj.sphericityLambda2()); DEBUG(0);
	_jStr_sphJj_sphericityLambda3->push_back(sphJj.sphericityLambda3()); DEBUG(0);
	_jStr_sphJj_spherocityLambda1->push_back(sphJj.spherocityLambda1()); DEBUG(0);
	_jStr_sphJj_spherocityLambda2->push_back(sphJj.spherocityLambda2()); DEBUG(0);
	_jStr_sphJj_spherocityLambda3->push_back(sphJj.spherocityLambda3()); DEBUG(0);
	_jStr_sphJj_circularity->push_back(sphJj.circularity()); DEBUG(0);
	_jStr_sphJj_sphericity->push_back(sphJj.sphericity()); DEBUG(0);
	_jStr_sphJj_spherocity->push_back(sphJj.spherocity()); DEBUG(0);
	_jStr_sphJj_aplanarity->push_back(sphJj.aplanarity()); DEBUG(0);
	_jStr_sphJj_aplanority->push_back(sphJj.aplanority()); DEBUG(0);
	_jStr_sphJj_Y->push_back(sphJj.Y()); DEBUG(0);
	_jStr_sphJj_planarity->push_back(sphJj.planarity()); DEBUG(0);
	_jStr_sphJj_planority->push_back(sphJj.planority()); DEBUG(0);
	_jStr_sphJj_Dshape->push_back(sphJj.Dshape()); DEBUG(0);
	_jStr_sphJj_Cshape->push_back(sphJj.Cshape()); DEBUG(0);
	_jStr_sphJj_H2->push_back(sphJj.H2()); DEBUG(0);
	_jStr_sphJj_Fshape->push_back(sphJj.Fshape()); DEBUG(0);
	_jStr_sphJj_beamThrust->push_back(sphJj.beamThrust()); DEBUG(0);
	_jStr_sphJj_G->push_back(sphJj.G()); DEBUG(0);
	_jStr_sphJj_ST2D->push_back(sphJj.ST2D()); DEBUG(0);
	_jStr_sphJj_detMlin->push_back(sphJj.detMlin()); DEBUG(0);
	_jStr_sphJj_pft->push_back(sphJj.pft()); DEBUG(0);
	_jStr_sphJi_detSphericity->push_back(sphJi.detSphericity()); DEBUG(0);
	_jStr_sphJi_detSpherocity->push_back(sphJi.detSpherocity()); DEBUG(0);
	_jStr_sphJi_sphericityLambda1->push_back(sphJi.sphericityLambda1()); DEBUG(0);
	_jStr_sphJi_sphericityLambda2->push_back(sphJi.sphericityLambda2()); DEBUG(0);
	_jStr_sphJi_sphericityLambda3->push_back(sphJi.sphericityLambda3()); DEBUG(0);
	_jStr_sphJi_spherocityLambda1->push_back(sphJi.spherocityLambda1()); DEBUG(0);
	_jStr_sphJi_spherocityLambda2->push_back(sphJi.spherocityLambda2()); DEBUG(0);
	_jStr_sphJi_spherocityLambda3->push_back(sphJi.spherocityLambda3()); DEBUG(0);
	_jStr_sphJi_circularity->push_back(sphJi.circularity()); DEBUG(0);
	_jStr_sphJi_sphericity->push_back(sphJi.sphericity()); DEBUG(0);
	_jStr_sphJi_spherocity->push_back(sphJi.spherocity()); DEBUG(0);
	_jStr_sphJi_aplanarity->push_back(sphJi.aplanarity()); DEBUG(0);
	_jStr_sphJi_aplanority->push_back(sphJi.aplanority()); DEBUG(0);
	_jStr_sphJi_Y->push_back(sphJi.Y()); DEBUG(0);
	_jStr_sphJi_planarity->push_back(sphJi.planarity()); DEBUG(0);
	_jStr_sphJi_planority->push_back(sphJi.planority()); DEBUG(0);
	_jStr_sphJi_Dshape->push_back(sphJi.Dshape()); DEBUG(0);
	_jStr_sphJi_Cshape->push_back(sphJi.Cshape()); DEBUG(0);
	_jStr_sphJi_H2->push_back(sphJi.H2()); DEBUG(0);
	_jStr_sphJi_Fshape->push_back(sphJi.Fshape()); DEBUG(0);
	_jStr_sphJi_beamThrust->push_back(sphJi.beamThrust()); DEBUG(0);
	_jStr_sphJi_G->push_back(sphJi.G()); DEBUG(0);
	_jStr_sphJi_ST2D->push_back(sphJi.ST2D()); DEBUG(0);
	_jStr_sphJi_detMlin->push_back(sphJi.detMlin()); DEBUG(0);
	_jStr_sphJi_pft->push_back(sphJi.pft()); DEBUG(0);
      
	delete aJet; DEBUG(0);
      } // loop over jets

      EventShape evtEJ(jetColl); DEBUG(0);
      RadiationVariables radEi(jetColl); DEBUG(0);

      FoxWolfram fwmEJ(jetColl, 0, false, maxTermsFW); DEBUG(0);
      //FoxWolfram fwmEi(jetColl, 0, true, 10); DEBUG(0);
      SphericitySpherocity sphEJ(jetColl); DEBUG(0);
      SphericitySpherocity sphEi(jetColl, 0, true); DEBUG(0);
    
      _jStr_authorE->push_back(JetAuthor()(jetColl)); DEBUG(0);
      _jStr_radialParamE->push_back(RadialParameter()(jetColl)); DEBUG(0);
      _jStr_algE->push_back(JetAlg()(jetColl)); DEBUG(0);
      _jStr_inputE->push_back(Input()(jetColl)); DEBUG(0);
      _jStr_bjetE->push_back(Bjet()(jetColl)); DEBUG(0);
      _jStr_xjetE->push_back(Xjet()(jetColl)); DEBUG(0);
      _jStr_myJetsE->push_back(MyJets()(jetColl)); DEBUG(0);

      _jStr_evtEJ_etaJ1->push_back(evtEJ.etaJ1()); DEBUG(0);
      _jStr_evtEJ_etaJ2->push_back(evtEJ.etaJ2()); DEBUG(0);
      _jStr_evtEJ_etaJ3->push_back(evtEJ.etaJ3()); DEBUG(0);
      _jStr_evtEJ_etaJ4->push_back(evtEJ.etaJ4()); DEBUG(0);
      _jStr_evtEJ_pTJ1->push_back(evtEJ.pTJ1()); DEBUG(0);
      _jStr_evtEJ_pTJ2->push_back(evtEJ.pTJ2()); DEBUG(0);
      _jStr_evtEJ_pTJ3->push_back(evtEJ.pTJ3()); DEBUG(0);
      _jStr_evtEJ_pTJ4->push_back(evtEJ.pTJ4()); DEBUG(0);

      _jStr_evtEJ_fE1->push_back(evtEJ.fE1()); DEBUG(0);
      _jStr_evtEJ_fE2->push_back(evtEJ.fE2()); DEBUG(0);
      _jStr_evtEJ_fE3->push_back(evtEJ.fE3()); DEBUG(0);
      _jStr_evtEJ_fET1->push_back(evtEJ.fET1()); DEBUG(0);
      _jStr_evtEJ_fET2->push_back(evtEJ.fET2()); DEBUG(0);
      _jStr_evtEJ_fET3->push_back(evtEJ.fET3()); DEBUG(0);
      _jStr_evtEJ_mct->push_back(evtEJ.mct()); DEBUG(0);
      _jStr_evtEJ_q->push_back(evtEJ.q()); DEBUG(0);
      _jStr_evtEJ_mjj->push_back(evtEJ.mjj()); DEBUG(0);
      _jStr_evtEJ_mTjj->push_back(evtEJ.mTjj()); DEBUG(0);
      _jStr_evtEJ_mjjj->push_back(evtEJ.mjjj()); DEBUG(0);
      _jStr_evtEJ_mTjjj->push_back(evtEJ.mTjjj()); DEBUG(0);
      _jStr_evtEJ_mjjjj->push_back(evtEJ.mjjjj()); DEBUG(0);
      _jStr_evtEJ_mTjjjj->push_back(evtEJ.mTjjjj()); DEBUG(0);
      _jStr_evtEJ_zetaPlus->push_back(evtEJ.zetaPlus()); DEBUG(0);
      _jStr_evtEJ_zetaMinus->push_back(evtEJ.zetaMinus()); DEBUG(0);
      _jStr_evtEJ_B12->push_back(evtEJ.B12()); DEBUG(0);
      _jStr_evtEJ_BT12->push_back(evtEJ.BT12()); DEBUG(0);
      _jStr_evtEJ_Alpha->push_back(evtEJ.Alpha()); DEBUG(0);
      _jStr_evtEJ_AlphaT->push_back(evtEJ.AlphaT()); DEBUG(0);
      _jStr_evtEJ_massDemocracy->push_back(evtEJ.massDemocracy()); DEBUG(0);
      _jStr_evtEJ_betaflow->push_back(evtEJ.betaflow()); DEBUG(0);
      _jStr_evtEJ_betaflow_1GeV->push_back(evtEJ.betaflow(1.*GeV)); DEBUG(0);
      _jStr_evtEJ_betaflow_5GeV->push_back(evtEJ.betaflow(5.*GeV)); DEBUG(0);
      _jStr_evtEJ_y23->push_back(evtEJ.y23()); DEBUG(0);
      _jStr_evtEJ_y23_1GeV->push_back(evtEJ.y23(1.*GeV)); DEBUG(0);
      _jStr_evtEJ_y23_5GeV->push_back(evtEJ.y23(5.*GeV)); DEBUG(0);
      _jStr_evtEJ_lny23->push_back(evtEJ.lny23()); DEBUG(0);
      _jStr_evtEJ_lny23_1GeV->push_back(evtEJ.lny23(1.*GeV)); DEBUG(0);
      _jStr_evtEJ_lny23_5GeV->push_back(evtEJ.lny23(5.*GeV)); DEBUG(0);
      _jStr_evtEJ_theta->push_back(evtEJ.theta()); DEBUG(0);
      _jStr_evtEJ_asym->push_back(evtEJ.asym()); DEBUG(0);
      _jStr_evtEJ_yB->push_back(evtEJ.yB()); DEBUG(0);
      _jStr_evtEJ_yStar->push_back(evtEJ.yStar()); DEBUG(0);
      _jStr_evtEJ_thetaStar->push_back(evtEJ.thetaStar()); DEBUG(0);
      _jStr_evtEJ_chi->push_back(evtEJ.chi()); DEBUG(0);
      _jStr_evtEJ_deltaPhiJJ->push_back(evtEJ.deltaPhiJJ()); DEBUG(0);
      _jStr_evtEJ_deltaThetaJJ->push_back(evtEJ.deltaThetaJJ()); DEBUG(0);
      _jStr_evtEJ_deltaEtaJJ->push_back(evtEJ.deltaEtaJJ()); DEBUG(0);
      _jStr_evtEJ_deltaRapidityJJ->push_back(evtEJ.deltaRapidityJJ()); DEBUG(0);
      _jStr_evtEJ_deltaRJJ->push_back(evtEJ.deltaRJJ()); DEBUG(0);
      _jStr_evtEJ_deltaRJJY->push_back(evtEJ.deltaRJJY()); DEBUG(0);
      _jStr_evtEJ_sigmaPhiJJ->push_back(evtEJ.sigmaPhiJJ()); DEBUG(0);
      _jStr_evtEJ_sigmaThetaJJ->push_back(evtEJ.sigmaThetaJJ()); DEBUG(0);
      _jStr_evtEJ_sigmaEtaJJ->push_back(evtEJ.sigmaEtaJJ()); DEBUG(0);
      _jStr_evtEJ_sigmaRapidityJJ->push_back(evtEJ.sigmaRapidityJJ()); DEBUG(0);
      _jStr_evtEJ_sigmaPtJJ->push_back(evtEJ.sigmaPtJJ()); DEBUG(0);
      _jStr_evtEJ_sigmaEtJJ->push_back(evtEJ.sigmaEtJJ()); DEBUG(0);
      _jStr_evtEJ_sigmaEt12->push_back(evtEJ.sigmaEt12()); DEBUG(0);
      _jStr_evtEJ_sigmaEt34->push_back(evtEJ.sigmaEt34()); DEBUG(0);
      _jStr_evtEJ_A234->push_back(evtEJ.A234()); DEBUG(0);
      _jStr_evtEJ_asymPhiJJ->push_back(evtEJ.asymPhiJJ()); DEBUG(0);
      _jStr_evtEJ_asymThetaJJ->push_back(evtEJ.asymThetaJJ()); DEBUG(0);
      _jStr_evtEJ_asymEtaJJ->push_back(evtEJ.asymEtaJJ()); DEBUG(0);
      _jStr_evtEJ_asymRapidityJJ->push_back(evtEJ.asymRapidityJJ()); DEBUG(0);
      _jStr_evtEJ_acoplanarity->push_back(evtEJ.acoplanarity()); DEBUG(0);
      _jStr_evtEJ_twist->push_back(evtEJ.twist()); DEBUG(0);
      _jStr_evtEJ_twistY->push_back(evtEJ.twistY()); DEBUG(0);
      _jStr_evtEJ_jetSumE->push_back(evtEJ.jetSumE()); DEBUG(0);
      _jStr_evtEJ_jetSumET->push_back(evtEJ.jetSumET()); DEBUG(0);
      _jStr_evtEJ_jetSumPT->push_back(evtEJ.jetSumPT()); DEBUG(0);
      _jStr_evtEJ_jetSumM->push_back(evtEJ.jetSumM()); DEBUG(0);
      _jStr_evtEJ_jetSumMT->push_back(evtEJ.jetSumMT()); DEBUG(0);
      _jStr_evtEJ_HTprime->push_back(evtEJ.HTprime()); DEBUG(0);
      _jStr_evtEJ_centrality->push_back(evtEJ.centrality()); DEBUG(0);
      _jStr_evtEJ_centralityP->push_back(evtEJ.centralityP()); DEBUG(0);
      _jStr_evtEJ_zminJ1J2->push_back(evtEJ.zminJ1J2()); DEBUG(0);
      _jStr_evtEJ_zmaxJ1J2->push_back(evtEJ.zmaxJ1J2()); DEBUG(0);
      _jStr_evtEJ_zminAllJets->push_back(evtEJ.zminAllJets()); DEBUG(0);
      _jStr_evtEJ_zmaxAllJets->push_back(evtEJ.zmaxAllJets()); DEBUG(0);
      _jStr_evtEJ_zminJ1J2Phi->push_back(evtEJ.zminJ1J2Phi()); DEBUG(0);
      _jStr_evtEJ_zminJ1J2Theta->push_back(evtEJ.zminJ1J2Theta()); DEBUG(0);
      _jStr_evtEJ_zminJ1J2Eta->push_back(evtEJ.zminJ1J2Eta()); DEBUG(0);
      _jStr_evtEJ_zminJ1J2Rapidity->push_back(evtEJ.zminJ1J2Rapidity()); DEBUG(0);
      _jStr_evtEJ_cosHelicityJ1->push_back(evtEJ.cosHelicityJ1()); DEBUG(0);
      _jStr_evtEJ_helicityJ1->push_back(evtEJ.helicityJ1()); DEBUG(0);
      _jStr_evtEJ_azilicityJ1->push_back(evtEJ.azilicityJ1()); DEBUG(0);
      _jStr_evtEJ_cosHelicityJ2->push_back(evtEJ.cosHelicityJ2()); DEBUG(0);
      _jStr_evtEJ_helicityJ2->push_back(evtEJ.helicityJ2()); DEBUG(0);
      _jStr_evtEJ_azilicityJ2->push_back(evtEJ.azilicityJ2()); DEBUG(0);
      _jStr_evtEJ_cosThetaJ1->push_back(evtEJ.cosThetaJ1()); DEBUG(0);
      _jStr_evtEJ_cosThetaJ2->push_back(evtEJ.cosThetaJ2()); DEBUG(0);
      _jStr_evtEJ_deltaRapidityXtoJ1CM->push_back(evtEJ.deltaRapidityXtoJ1CM()); DEBUG(0);
      _jStr_evtEJ_deltaRapidityXtoJ2CM->push_back(evtEJ.deltaRapidityXtoJ2CM()); DEBUG(0);
      _jStr_evtEJ_deltaRapidityXtoJ1->push_back(evtEJ.deltaRapidityXtoJ1()); DEBUG(0);
      _jStr_evtEJ_deltaRapidityXtoJ2->push_back(evtEJ.deltaRapidityXtoJ2()); DEBUG(0);
      _jStr_evtEJ_cosTheta1->push_back(evtEJ.cosTheta1()); DEBUG(0);
      _jStr_evtEJ_cosTheta2->push_back(evtEJ.cosTheta2()); DEBUG(0);
      _jStr_evtEJ_cosThetaStar1->push_back(evtEJ.cosThetaStar1()); DEBUG(0);
      _jStr_evtEJ_cosThetaStar2->push_back(evtEJ.cosThetaStar2()); DEBUG(0);
      _jStr_evtEJ_cosPhiTilde1->push_back(evtEJ.cosPhiTilde1()); DEBUG(0);
      _jStr_evtEJ_cosPhiTilde2->push_back(evtEJ.cosPhiTilde2()); DEBUG(0);
      _jStr_evtEJ_cosPhi->push_back(evtEJ.cosPhi()); DEBUG(0);
      _jStr_evtEJ_PhiTilde1->push_back(evtEJ.PhiTilde1()); DEBUG(0);
      _jStr_evtEJ_PhiTilde2->push_back(evtEJ.PhiTilde2()); DEBUG(0);
      _jStr_evtEJ_Phi->push_back(evtEJ.Phi()); DEBUG(0);
      _jStr_evtEJ_M->push_back(evtEJ.M()); DEBUG(0);
      _jStr_evtEJ_DeltaM->push_back(evtEJ.DeltaM()); DEBUG(0);
      _jStr_evtEJ_asymM->push_back(evtEJ.asymM()); DEBUG(0);
      _jStr_evtEJ_Q->push_back(evtEJ.Q()); DEBUG(0);
      _jStr_evtEJ_SCDF->push_back(evtEJ.SCDF()); DEBUG(0);
      _jStr_evtEJ_dPhiIJCDF->push_back(evtEJ.dPhiIJCDF()); DEBUG(0);
      _jStr_evtEJ_dPhiKLCDF->push_back(evtEJ.dPhiKLCDF()); DEBUG(0);
      _jStr_evtEJ_PtIPtJCDF->push_back(evtEJ.PtIPtJCDF()); DEBUG(0);
      _jStr_evtEJ_PtKPtLCDF->push_back(evtEJ.PtKPtLCDF()); DEBUG(0);
      _jStr_evtEJ_SD0->push_back(evtEJ.SD0()); DEBUG(0);
      _jStr_evtEJ_dPhiIJD0->push_back(evtEJ.dPhiIJD0()); DEBUG(0);
      _jStr_evtEJ_dPhiKLD0->push_back(evtEJ.dPhiKLD0()); DEBUG(0);
      _jStr_evtEJ_PtIPtJD0->push_back(evtEJ.PtIPtJD0()); DEBUG(0);
      _jStr_evtEJ_PtKPtLD0->push_back(evtEJ.PtKPtLD0()); DEBUG(0);
      _jStr_evtEJ_DeltaSCDF->push_back(evtEJ.DeltaSCDF()); DEBUG(0);
      _jStr_evtEJ_DeltaSD0->push_back(evtEJ.DeltaSD0()); DEBUG(0);
      _jStr_evtEJ_Delta12->push_back(evtEJ.Delta12()); DEBUG(0);
      _jStr_evtEJ_Delta12n->push_back(evtEJ.Delta12n()); DEBUG(0);
      _jStr_evtEJ_M14->push_back(evtEJ.M14()); DEBUG(0);
      _jStr_evtEJ_y1y2->push_back(evtEJ.y1y2()); DEBUG(0);
    
      _jStr_evtEJ_Mprime4->push_back(evtEJ.Mprime4()); DEBUG(0);
      _jStr_evtEJ_dMprime4->push_back(evtEJ.dMprime4()); DEBUG(0);
      _jStr_evtEJ_MprimeAvg4->push_back(evtEJ.MprimeAvg4()); DEBUG(0);
      _jStr_evtEJ_DeltaMin4->push_back(evtEJ.DeltaMin4()); DEBUG(0);
      _jStr_evtEJ_DeltaMax4->push_back(evtEJ.DeltaMax4()); DEBUG(0);
      _jStr_evtEJ_DeltaPhiXX4->push_back(evtEJ.DeltaPhiXX4()); DEBUG(0);
      _jStr_evtEJ_DeltaYXX4->push_back(evtEJ.DeltaYXX4()); DEBUG(0);
      _jStr_evtEJ_DeltaRXX4->push_back(evtEJ.DeltaRXX4()); DEBUG(0);
      _jStr_evtEJ_TwistXX4->push_back(evtEJ.TwistXX4()); DEBUG(0);
      _jStr_evtEJ_separated4->push_back(evtEJ.separated4()); DEBUG(0);
      _jStr_evtEJ_Mprime->push_back(evtEJ.Mprime()); DEBUG(0);
      _jStr_evtEJ_dMprime->push_back(evtEJ.dMprime()); DEBUG(0);
      _jStr_evtEJ_MprimeAvg->push_back(evtEJ.MprimeAvg()); DEBUG(0);
      _jStr_evtEJ_DeltaMin->push_back(evtEJ.DeltaMin()); DEBUG(0);
      _jStr_evtEJ_DeltaMax->push_back(evtEJ.DeltaMax()); DEBUG(0);
      _jStr_evtEJ_DeltaPhiXX->push_back(evtEJ.DeltaPhiXX()); DEBUG(0);
      _jStr_evtEJ_DeltaYXX->push_back(evtEJ.DeltaYXX()); DEBUG(0);
      _jStr_evtEJ_DeltaRXX->push_back(evtEJ.DeltaRXX()); DEBUG(0);
      _jStr_evtEJ_TwistXX->push_back(evtEJ.TwistXX()); DEBUG(0);
      _jStr_evtEJ_separated->push_back(evtEJ.separated()); DEBUG(0);
    
      _jStr_fwmEJ_xiPlus->push_back(fwmEJ.xiPlus()); DEBUG(0);
      _jStr_fwmEJ_xiMinus->push_back(fwmEJ.xiMinus()); DEBUG(0);
      _jStr_fwmEJ_xPlus->push_back(fwmEJ.xPlus()); DEBUG(0);
      _jStr_fwmEJ_xMinus->push_back(fwmEJ.xMinus()); DEBUG(0);

      _jStr_fwmEJ_Psi1->push_back(fwmEJ.Psi1()); DEBUG(0);
      _jStr_fwmEJ_B0->push_back(fwmEJ.B0()); DEBUG(0);
      _jStr_fwmEJ_B1->push_back(fwmEJ.B1()); DEBUG(0);
      _jStr_fwmEJ_B2->push_back(fwmEJ.B2()); DEBUG(0);
      _jStr_fwmEJ_B3->push_back(fwmEJ.B3()); DEBUG(0);
      _jStr_fwmEJ_B4->push_back(fwmEJ.B4()); DEBUG(0);
      _jStr_fwmEJ_B5->push_back(fwmEJ.B5()); DEBUG(0);
      _jStr_fwmEJ_B6->push_back(fwmEJ.B6()); DEBUG(0);
      _jStr_fwmEJ_B7->push_back(fwmEJ.B7()); DEBUG(0);
      _jStr_fwmEJ_B8->push_back(fwmEJ.B8()); DEBUG(0);
      _jStr_fwmEJ_C0->push_back(fwmEJ.C0()); DEBUG(0);
      _jStr_fwmEJ_C1->push_back(fwmEJ.C1()); DEBUG(0);
      _jStr_fwmEJ_C2->push_back(fwmEJ.C2()); DEBUG(0);
      _jStr_fwmEJ_C3->push_back(fwmEJ.C3()); DEBUG(0);
      _jStr_fwmEJ_C4->push_back(fwmEJ.C4()); DEBUG(0);
      _jStr_fwmEJ_C5->push_back(fwmEJ.C5()); DEBUG(0);
      _jStr_fwmEJ_C6->push_back(fwmEJ.C6()); DEBUG(0);
      _jStr_fwmEJ_C7->push_back(fwmEJ.C7()); DEBUG(0);
      _jStr_fwmEJ_C8->push_back(fwmEJ.C8()); DEBUG(0);
      _jStr_fwmEJ_K0->push_back(fwmEJ.K0()); DEBUG(0);
      _jStr_fwmEJ_K1->push_back(fwmEJ.K1()); DEBUG(0);
      _jStr_fwmEJ_K2->push_back(fwmEJ.K2()); DEBUG(0);
      _jStr_fwmEJ_K3->push_back(fwmEJ.K3()); DEBUG(0);
      _jStr_fwmEJ_K4->push_back(fwmEJ.K4()); DEBUG(0);
      _jStr_fwmEJ_K5->push_back(fwmEJ.K5()); DEBUG(0);
      _jStr_fwmEJ_K6->push_back(fwmEJ.K6()); DEBUG(0);
      _jStr_fwmEJ_K7->push_back(fwmEJ.K7()); DEBUG(0);
      _jStr_fwmEJ_K8->push_back(fwmEJ.K8()); DEBUG(0);
      _jStr_fwmEJ_D0->push_back(fwmEJ.D0()); DEBUG(0);
      _jStr_fwmEJ_D1->push_back(fwmEJ.D1()); DEBUG(0);
      _jStr_fwmEJ_D2->push_back(fwmEJ.D2()); DEBUG(0);
      _jStr_fwmEJ_D3->push_back(fwmEJ.D3()); DEBUG(0);
      _jStr_fwmEJ_D4->push_back(fwmEJ.D4()); DEBUG(0);
      _jStr_fwmEJ_D5->push_back(fwmEJ.D5()); DEBUG(0);
      _jStr_fwmEJ_D6->push_back(fwmEJ.D6()); DEBUG(0);
      _jStr_fwmEJ_D7->push_back(fwmEJ.D7()); DEBUG(0);
      _jStr_fwmEJ_D8->push_back(fwmEJ.D8()); DEBUG(0);
      _jStr_fwmEJ_H0->push_back(fwmEJ.H0()); DEBUG(0);
      _jStr_fwmEJ_H1->push_back(fwmEJ.H1()); DEBUG(0);
      _jStr_fwmEJ_H2->push_back(fwmEJ.H2()); DEBUG(0);
      _jStr_fwmEJ_H3->push_back(fwmEJ.H3()); DEBUG(0);
      _jStr_fwmEJ_H4->push_back(fwmEJ.H4()); DEBUG(0);
      _jStr_fwmEJ_H5->push_back(fwmEJ.H5()); DEBUG(0);
      _jStr_fwmEJ_H6->push_back(fwmEJ.H6()); DEBUG(0);
      _jStr_fwmEJ_H7->push_back(fwmEJ.H7()); DEBUG(0);
      _jStr_fwmEJ_H8->push_back(fwmEJ.H8()); DEBUG(0);
      _jStr_fwmEJ_Q0->push_back(fwmEJ.Q0()); DEBUG(0);
      _jStr_fwmEJ_Q1->push_back(fwmEJ.Q1()); DEBUG(0);
      _jStr_fwmEJ_Q2->push_back(fwmEJ.Q2()); DEBUG(0);
      _jStr_fwmEJ_Q3->push_back(fwmEJ.Q3()); DEBUG(0);
      _jStr_fwmEJ_Q4->push_back(fwmEJ.Q4()); DEBUG(0);
      _jStr_fwmEJ_Q5->push_back(fwmEJ.Q5()); DEBUG(0);
      _jStr_fwmEJ_Q6->push_back(fwmEJ.Q6()); DEBUG(0);
      _jStr_fwmEJ_Q7->push_back(fwmEJ.Q7()); DEBUG(0);
      _jStr_fwmEJ_Q8->push_back(fwmEJ.Q8()); DEBUG(0);
      _jStr_fwmEJ_Pi1->push_back(fwmEJ.Pi1()); DEBUG(0);
      _jStr_fwmEJ_Pi2->push_back(fwmEJ.Pi2()); DEBUG(0);
      _jStr_fwmEJ_Pi3->push_back(fwmEJ.Pi3()); DEBUG(0);
      _jStr_fwmEJ_Pi4->push_back(fwmEJ.Pi4()); DEBUG(0);
      _jStr_fwmEJ_B10->push_back(fwmEJ.B10()); DEBUG(0);
      _jStr_fwmEJ_B20->push_back(fwmEJ.B20()); DEBUG(0);
      _jStr_fwmEJ_B30->push_back(fwmEJ.B30()); DEBUG(0);
      _jStr_fwmEJ_B40->push_back(fwmEJ.B40()); DEBUG(0);
      _jStr_fwmEJ_B50->push_back(fwmEJ.B50()); DEBUG(0);
      _jStr_fwmEJ_B60->push_back(fwmEJ.B60()); DEBUG(0);
      _jStr_fwmEJ_B70->push_back(fwmEJ.B70()); DEBUG(0);
      _jStr_fwmEJ_B80->push_back(fwmEJ.B80()); DEBUG(0);
      _jStr_fwmEJ_C10->push_back(fwmEJ.C10()); DEBUG(0);
      _jStr_fwmEJ_C20->push_back(fwmEJ.C20()); DEBUG(0);
      _jStr_fwmEJ_C30->push_back(fwmEJ.C30()); DEBUG(0);
      _jStr_fwmEJ_C40->push_back(fwmEJ.C40()); DEBUG(0);
      _jStr_fwmEJ_C50->push_back(fwmEJ.C50()); DEBUG(0);
      _jStr_fwmEJ_C60->push_back(fwmEJ.C60()); DEBUG(0);
      _jStr_fwmEJ_C70->push_back(fwmEJ.C70()); DEBUG(0);
      _jStr_fwmEJ_C80->push_back(fwmEJ.C80()); DEBUG(0);
      _jStr_fwmEJ_K10->push_back(fwmEJ.K10()); DEBUG(0);
      _jStr_fwmEJ_K20->push_back(fwmEJ.K20()); DEBUG(0);
      _jStr_fwmEJ_K30->push_back(fwmEJ.K30()); DEBUG(0);
      _jStr_fwmEJ_K40->push_back(fwmEJ.K40()); DEBUG(0);
      _jStr_fwmEJ_K50->push_back(fwmEJ.K50()); DEBUG(0);
      _jStr_fwmEJ_K60->push_back(fwmEJ.K60()); DEBUG(0);
      _jStr_fwmEJ_K70->push_back(fwmEJ.K70()); DEBUG(0);
      _jStr_fwmEJ_K80->push_back(fwmEJ.K80()); DEBUG(0);
      _jStr_fwmEJ_D10->push_back(fwmEJ.D10()); DEBUG(0);
      _jStr_fwmEJ_D20->push_back(fwmEJ.D20()); DEBUG(0);
      _jStr_fwmEJ_D30->push_back(fwmEJ.D30()); DEBUG(0);
      _jStr_fwmEJ_D40->push_back(fwmEJ.D40()); DEBUG(0);
      _jStr_fwmEJ_D50->push_back(fwmEJ.D50()); DEBUG(0);
      _jStr_fwmEJ_D60->push_back(fwmEJ.D60()); DEBUG(0);
      _jStr_fwmEJ_D70->push_back(fwmEJ.D70()); DEBUG(0);
      _jStr_fwmEJ_D80->push_back(fwmEJ.D80()); DEBUG(0);
      _jStr_fwmEJ_H10->push_back(fwmEJ.H10()); DEBUG(0);
      _jStr_fwmEJ_H20->push_back(fwmEJ.H20()); DEBUG(0);
      _jStr_fwmEJ_H30->push_back(fwmEJ.H30()); DEBUG(0);
      _jStr_fwmEJ_H40->push_back(fwmEJ.H40()); DEBUG(0);
      _jStr_fwmEJ_H50->push_back(fwmEJ.H50()); DEBUG(0);
      _jStr_fwmEJ_H60->push_back(fwmEJ.H60()); DEBUG(0);
      _jStr_fwmEJ_H70->push_back(fwmEJ.H70()); DEBUG(0);
      _jStr_fwmEJ_H80->push_back(fwmEJ.H80()); DEBUG(0);
      _jStr_fwmEJ_Q10->push_back(fwmEJ.Q10()); DEBUG(0);
      _jStr_fwmEJ_Q20->push_back(fwmEJ.Q20()); DEBUG(0);
      _jStr_fwmEJ_Q30->push_back(fwmEJ.Q30()); DEBUG(0);
      _jStr_fwmEJ_Q40->push_back(fwmEJ.Q40()); DEBUG(0);
      _jStr_fwmEJ_Q50->push_back(fwmEJ.Q50()); DEBUG(0);
      _jStr_fwmEJ_Q60->push_back(fwmEJ.Q60()); DEBUG(0);
      _jStr_fwmEJ_Q70->push_back(fwmEJ.Q70()); DEBUG(0);
      _jStr_fwmEJ_Q80->push_back(fwmEJ.Q80()); DEBUG(0);
    
      _jStr_sphEJ_detSphericity->push_back(sphEJ.detSphericity()); DEBUG(0);
      _jStr_sphEJ_detSpherocity->push_back(sphEJ.detSpherocity()); DEBUG(0);
      _jStr_sphEJ_sphericityLambda1->push_back(sphEJ.sphericityLambda1()); DEBUG(0);
      _jStr_sphEJ_sphericityLambda2->push_back(sphEJ.sphericityLambda2()); DEBUG(0);
      _jStr_sphEJ_sphericityLambda3->push_back(sphEJ.sphericityLambda3()); DEBUG(0);
      _jStr_sphEJ_spherocityLambda1->push_back(sphEJ.spherocityLambda1()); DEBUG(0);
      _jStr_sphEJ_spherocityLambda2->push_back(sphEJ.spherocityLambda2()); DEBUG(0);
      _jStr_sphEJ_spherocityLambda3->push_back(sphEJ.spherocityLambda3()); DEBUG(0);
      _jStr_sphEJ_circularity->push_back(sphEJ.circularity()); DEBUG(0);
      _jStr_sphEJ_sphericity->push_back(sphEJ.sphericity()); DEBUG(0);
      _jStr_sphEJ_spherocity->push_back(sphEJ.spherocity()); DEBUG(0);
      _jStr_sphEJ_aplanarity->push_back(sphEJ.aplanarity()); DEBUG(0);
      _jStr_sphEJ_aplanority->push_back(sphEJ.aplanority()); DEBUG(0);
      _jStr_sphEJ_Y->push_back(sphEJ.Y()); DEBUG(0);
      _jStr_sphEJ_planarity->push_back(sphEJ.planarity()); DEBUG(0);
      _jStr_sphEJ_planority->push_back(sphEJ.planority()); DEBUG(0);
      _jStr_sphEJ_Dshape->push_back(sphEJ.Dshape()); DEBUG(0);
      _jStr_sphEJ_Cshape->push_back(sphEJ.Cshape()); DEBUG(0);
      _jStr_sphEJ_H2->push_back(sphEJ.H2()); DEBUG(0);
      _jStr_sphEJ_Fshape->push_back(sphEJ.Fshape()); DEBUG(0);
      _jStr_sphEJ_beamThrust->push_back(sphEJ.beamThrust()); DEBUG(0);
      _jStr_sphEJ_G->push_back(sphEJ.G()); DEBUG(0);
      _jStr_sphEJ_ST2D->push_back(sphEJ.ST2D()); DEBUG(0);
      _jStr_sphEJ_detMlin->push_back(sphEJ.detMlin()); DEBUG(0);
      _jStr_sphEJ_pft->push_back(sphEJ.pft()); DEBUG(0);
    
      // _jStr_fwmEi_Psi1->push_back(fwmEi.Psi1()); DEBUG(0);
      // _jStr_fwmEi_B0->push_back(fwmEi.B0()); DEBUG(0);
      // _jStr_fwmEi_B1->push_back(fwmEi.B1()); DEBUG(0);
      // _jStr_fwmEi_B2->push_back(fwmEi.B2()); DEBUG(0);
      // _jStr_fwmEi_B3->push_back(fwmEi.B3()); DEBUG(0);
      // _jStr_fwmEi_B4->push_back(fwmEi.B4()); DEBUG(0);
      // _jStr_fwmEi_B5->push_back(fwmEi.B5()); DEBUG(0);
      // _jStr_fwmEi_B6->push_back(fwmEi.B6()); DEBUG(0);
      // _jStr_fwmEi_B7->push_back(fwmEi.B7()); DEBUG(0);
      // _jStr_fwmEi_B8->push_back(fwmEi.B8()); DEBUG(0);
      // _jStr_fwmEi_C0->push_back(fwmEi.C0()); DEBUG(0);
      // _jStr_fwmEi_C1->push_back(fwmEi.C1()); DEBUG(0);
      // _jStr_fwmEi_C2->push_back(fwmEi.C2()); DEBUG(0);
      // _jStr_fwmEi_C3->push_back(fwmEi.C3()); DEBUG(0);
      // _jStr_fwmEi_C4->push_back(fwmEi.C4()); DEBUG(0);
      // _jStr_fwmEi_C5->push_back(fwmEi.C5()); DEBUG(0);
      // _jStr_fwmEi_C6->push_back(fwmEi.C6()); DEBUG(0);
      // _jStr_fwmEi_C7->push_back(fwmEi.C7()); DEBUG(0);
      // _jStr_fwmEi_C8->push_back(fwmEi.C8()); DEBUG(0);
      // _jStr_fwmEi_K0->push_back(fwmEi.K0()); DEBUG(0);
      // _jStr_fwmEi_K1->push_back(fwmEi.K1()); DEBUG(0);
      // _jStr_fwmEi_K2->push_back(fwmEi.K2()); DEBUG(0);
      // _jStr_fwmEi_K3->push_back(fwmEi.K3()); DEBUG(0);
      // _jStr_fwmEi_K4->push_back(fwmEi.K4()); DEBUG(0);
      // _jStr_fwmEi_K5->push_back(fwmEi.K5()); DEBUG(0);
      // _jStr_fwmEi_K6->push_back(fwmEi.K6()); DEBUG(0);
      // _jStr_fwmEi_K7->push_back(fwmEi.K7()); DEBUG(0);
      // _jStr_fwmEi_K8->push_back(fwmEi.K8()); DEBUG(0);
      // _jStr_fwmEi_D0->push_back(fwmEi.D0()); DEBUG(0);
      // _jStr_fwmEi_D1->push_back(fwmEi.D1()); DEBUG(0);
      // _jStr_fwmEi_D2->push_back(fwmEi.D2()); DEBUG(0);
      // _jStr_fwmEi_D3->push_back(fwmEi.D3()); DEBUG(0);
      // _jStr_fwmEi_D4->push_back(fwmEi.D4()); DEBUG(0);
      // _jStr_fwmEi_D5->push_back(fwmEi.D5()); DEBUG(0);
      // _jStr_fwmEi_D6->push_back(fwmEi.D6()); DEBUG(0);
      // _jStr_fwmEi_D7->push_back(fwmEi.D7()); DEBUG(0);
      // _jStr_fwmEi_D8->push_back(fwmEi.D8()); DEBUG(0);
      // _jStr_fwmEi_H0->push_back(fwmEi.H0()); DEBUG(0);
      // _jStr_fwmEi_H1->push_back(fwmEi.H1()); DEBUG(0);
      // _jStr_fwmEi_H2->push_back(fwmEi.H2()); DEBUG(0);
      // _jStr_fwmEi_H3->push_back(fwmEi.H3()); DEBUG(0);
      // _jStr_fwmEi_H4->push_back(fwmEi.H4()); DEBUG(0);
      // _jStr_fwmEi_H5->push_back(fwmEi.H5()); DEBUG(0);
      // _jStr_fwmEi_H6->push_back(fwmEi.H6()); DEBUG(0);
      // _jStr_fwmEi_H7->push_back(fwmEi.H7()); DEBUG(0);
      // _jStr_fwmEi_H8->push_back(fwmEi.H8()); DEBUG(0);
      // _jStr_fwmEi_Q0->push_back(fwmEi.Q0()); DEBUG(0);
      // _jStr_fwmEi_Q1->push_back(fwmEi.Q1()); DEBUG(0);
      // _jStr_fwmEi_Q2->push_back(fwmEi.Q2()); DEBUG(0);
      // _jStr_fwmEi_Q3->push_back(fwmEi.Q3()); DEBUG(0);
      // _jStr_fwmEi_Q4->push_back(fwmEi.Q4()); DEBUG(0);
      // _jStr_fwmEi_Q5->push_back(fwmEi.Q5()); DEBUG(0);
      // _jStr_fwmEi_Q6->push_back(fwmEi.Q6()); DEBUG(0);
      // _jStr_fwmEi_Q7->push_back(fwmEi.Q7()); DEBUG(0);
      // _jStr_fwmEi_Q8->push_back(fwmEi.Q8()); DEBUG(0);
      // _jStr_fwmEi_Pi1->push_back(fwmEi.Pi1()); DEBUG(0);
      // _jStr_fwmEi_Pi2->push_back(fwmEi.Pi2()); DEBUG(0);
      // _jStr_fwmEi_Pi3->push_back(fwmEi.Pi3()); DEBUG(0);
      // _jStr_fwmEi_Pi4->push_back(fwmEi.Pi4()); DEBUG(0);
      // _jStr_fwmEi_B10->push_back(fwmEi.B10()); DEBUG(0);
      // _jStr_fwmEi_B20->push_back(fwmEi.B20()); DEBUG(0);
      // _jStr_fwmEi_B30->push_back(fwmEi.B30()); DEBUG(0);
      // _jStr_fwmEi_B40->push_back(fwmEi.B40()); DEBUG(0);
      // _jStr_fwmEi_B50->push_back(fwmEi.B50()); DEBUG(0);
      // _jStr_fwmEi_B60->push_back(fwmEi.B60()); DEBUG(0);
      // _jStr_fwmEi_B70->push_back(fwmEi.B70()); DEBUG(0);
      // _jStr_fwmEi_B80->push_back(fwmEi.B80()); DEBUG(0);
      // _jStr_fwmEi_C10->push_back(fwmEi.C10()); DEBUG(0);
      // _jStr_fwmEi_C20->push_back(fwmEi.C20()); DEBUG(0);
      // _jStr_fwmEi_C30->push_back(fwmEi.C30()); DEBUG(0);
      // _jStr_fwmEi_C40->push_back(fwmEi.C40()); DEBUG(0);
      // _jStr_fwmEi_C50->push_back(fwmEi.C50()); DEBUG(0);
      // _jStr_fwmEi_C60->push_back(fwmEi.C60()); DEBUG(0);
      // _jStr_fwmEi_C70->push_back(fwmEi.C70()); DEBUG(0);
      // _jStr_fwmEi_C80->push_back(fwmEi.C80()); DEBUG(0);
      // _jStr_fwmEi_K10->push_back(fwmEi.K10()); DEBUG(0);
      // _jStr_fwmEi_K20->push_back(fwmEi.K20()); DEBUG(0);
      // _jStr_fwmEi_K30->push_back(fwmEi.K30()); DEBUG(0);
      // _jStr_fwmEi_K40->push_back(fwmEi.K40()); DEBUG(0);
      // _jStr_fwmEi_K50->push_back(fwmEi.K50()); DEBUG(0);
      // _jStr_fwmEi_K60->push_back(fwmEi.K60()); DEBUG(0);
      // _jStr_fwmEi_K70->push_back(fwmEi.K70()); DEBUG(0);
      // _jStr_fwmEi_K80->push_back(fwmEi.K80()); DEBUG(0);
      // _jStr_fwmEi_D10->push_back(fwmEi.D10()); DEBUG(0);
      // _jStr_fwmEi_D20->push_back(fwmEi.D20()); DEBUG(0);
      // _jStr_fwmEi_D30->push_back(fwmEi.D30()); DEBUG(0);
      // _jStr_fwmEi_D40->push_back(fwmEi.D40()); DEBUG(0);
      // _jStr_fwmEi_D50->push_back(fwmEi.D50()); DEBUG(0);
      // _jStr_fwmEi_D60->push_back(fwmEi.D60()); DEBUG(0);
      // _jStr_fwmEi_D70->push_back(fwmEi.D70()); DEBUG(0);
      // _jStr_fwmEi_D80->push_back(fwmEi.D80()); DEBUG(0);
      // _jStr_fwmEi_H10->push_back(fwmEi.H10()); DEBUG(0);
      // _jStr_fwmEi_H20->push_back(fwmEi.H20()); DEBUG(0);
      // _jStr_fwmEi_H30->push_back(fwmEi.H30()); DEBUG(0);
      // _jStr_fwmEi_H40->push_back(fwmEi.H40()); DEBUG(0);
      // _jStr_fwmEi_H50->push_back(fwmEi.H50()); DEBUG(0);
      // _jStr_fwmEi_H60->push_back(fwmEi.H60()); DEBUG(0);
      // _jStr_fwmEi_H70->push_back(fwmEi.H70()); DEBUG(0);
      // _jStr_fwmEi_H80->push_back(fwmEi.H80()); DEBUG(0);
      // _jStr_fwmEi_Q10->push_back(fwmEi.Q10()); DEBUG(0);
      // _jStr_fwmEi_Q20->push_back(fwmEi.Q20()); DEBUG(0);
      // _jStr_fwmEi_Q30->push_back(fwmEi.Q30()); DEBUG(0);
      // _jStr_fwmEi_Q40->push_back(fwmEi.Q40()); DEBUG(0);
      // _jStr_fwmEi_Q50->push_back(fwmEi.Q50()); DEBUG(0);
      // _jStr_fwmEi_Q60->push_back(fwmEi.Q60()); DEBUG(0);
      // _jStr_fwmEi_Q70->push_back(fwmEi.Q70()); DEBUG(0);
      // _jStr_fwmEi_Q80->push_back(fwmEi.Q80()); DEBUG(0);
    
      _jStr_radEi_bcfJ1_v0a1->push_back(radEi.bcfJ1(bcfVersion, 1)); DEBUG(0);
      _jStr_radEi_bcfJ2_v0a1->push_back(radEi.bcfJ2(bcfVersion, 1)); DEBUG(0);
      _jStr_radEi_bcfJ3_v0a1->push_back(radEi.bcfJ3(bcfVersion, 1)); DEBUG(0);
      _jStr_radEi_bcfTJ1_v0a1->push_back(radEi.bcfTJ1(bcfVersion, 1)); DEBUG(0);
      _jStr_radEi_bcfTJ2_v0a1->push_back(radEi.bcfTJ2(bcfVersion, 1)); DEBUG(0);
      _jStr_radEi_bcfTJ3_v0a1->push_back(radEi.bcfTJ3(bcfVersion, 1)); DEBUG(0);
      _jStr_radEi_bcfAsymYJ1_v0a1->push_back(radEi.bcfAsymYJ1(bcfVersion, 1)); DEBUG(0);
      _jStr_radEi_bcfAsymYJ2_v0a1->push_back(radEi.bcfAsymYJ2(bcfVersion, 1)); DEBUG(0);
      _jStr_radEi_bcfAsymYJ3_v0a1->push_back(radEi.bcfAsymYJ3(bcfVersion, 1)); DEBUG(0);
      _jStr_radEi_bcfAsymPhiJ1_v0a1->push_back(radEi.bcfAsymPhiJ1(bcfVersion, 1)); DEBUG(0);
      _jStr_radEi_bcfAsymPhiJ2_v0a1->push_back(radEi.bcfAsymPhiJ2(bcfVersion, 1)); DEBUG(0);
      _jStr_radEi_bcfAsymPhiJ3_v0a1->push_back(radEi.bcfAsymPhiJ3(bcfVersion, 1)); DEBUG(0);
      _jStr_radEi_bcfAsymYPhiJ1_v0a1->push_back(radEi.bcfAsymYPhiJ1(bcfVersion, 1)); DEBUG(0);
      _jStr_radEi_bcfAsymYPhiJ2_v0a1->push_back(radEi.bcfAsymYPhiJ2(bcfVersion, 1)); DEBUG(0);
      _jStr_radEi_bcfAsymYPhiJ3_v0a1->push_back(radEi.bcfAsymYPhiJ3(bcfVersion, 1)); DEBUG(0);
      _jStr_radEi_bcfAsymYPhi2J1_v0a1->push_back(radEi.bcfAsymYPhi2J1(bcfVersion, 1)); DEBUG(0);
      _jStr_radEi_bcfAsymYPhi2J2_v0a1->push_back(radEi.bcfAsymYPhi2J2(bcfVersion, 1)); DEBUG(0);
      _jStr_radEi_bcfAsymYPhi2J3_v0a1->push_back(radEi.bcfAsymYPhi2J3(bcfVersion, 1)); DEBUG(0);
      _jStr_radEi_bcfJ1_v0a2->push_back(radEi.bcfJ1(bcfVersion, 2)); DEBUG(0);
      _jStr_radEi_bcfJ2_v0a2->push_back(radEi.bcfJ2(bcfVersion, 2)); DEBUG(0);
      _jStr_radEi_bcfJ3_v0a2->push_back(radEi.bcfJ3(bcfVersion, 2)); DEBUG(0);
      _jStr_radEi_bcfTJ1_v0a2->push_back(radEi.bcfTJ1(bcfVersion, 2)); DEBUG(0);
      _jStr_radEi_bcfTJ2_v0a2->push_back(radEi.bcfTJ2(bcfVersion, 2)); DEBUG(0);
      _jStr_radEi_bcfTJ3_v0a2->push_back(radEi.bcfTJ3(bcfVersion, 2)); DEBUG(0);
      _jStr_radEi_bcfAsymYJ1_v0a2->push_back(radEi.bcfAsymYJ1(bcfVersion, 2)); DEBUG(0);
      _jStr_radEi_bcfAsymYJ2_v0a2->push_back(radEi.bcfAsymYJ2(bcfVersion, 2)); DEBUG(0);
      _jStr_radEi_bcfAsymYJ3_v0a2->push_back(radEi.bcfAsymYJ3(bcfVersion, 2)); DEBUG(0);
      _jStr_radEi_bcfAsymPhiJ1_v0a2->push_back(radEi.bcfAsymPhiJ1(bcfVersion, 2)); DEBUG(0);
      _jStr_radEi_bcfAsymPhiJ2_v0a2->push_back(radEi.bcfAsymPhiJ2(bcfVersion, 2)); DEBUG(0);
      _jStr_radEi_bcfAsymPhiJ3_v0a2->push_back(radEi.bcfAsymPhiJ3(bcfVersion, 2)); DEBUG(0);
      _jStr_radEi_bcfAsymYPhiJ1_v0a2->push_back(radEi.bcfAsymYPhiJ1(bcfVersion, 2)); DEBUG(0);
      _jStr_radEi_bcfAsymYPhiJ2_v0a2->push_back(radEi.bcfAsymYPhiJ2(bcfVersion, 2)); DEBUG(0);
      _jStr_radEi_bcfAsymYPhiJ3_v0a2->push_back(radEi.bcfAsymYPhiJ3(bcfVersion, 2)); DEBUG(0);
      _jStr_radEi_bcfAsymYPhi2J1_v0a2->push_back(radEi.bcfAsymYPhi2J1(bcfVersion, 2)); DEBUG(0);
      _jStr_radEi_bcfAsymYPhi2J2_v0a2->push_back(radEi.bcfAsymYPhi2J2(bcfVersion, 2)); DEBUG(0);
      _jStr_radEi_bcfAsymYPhi2J3_v0a2->push_back(radEi.bcfAsymYPhi2J3(bcfVersion, 2)); DEBUG(0);
      _jStr_radEi_BJ1->push_back(radEi.BJ1()); DEBUG(0);
      _jStr_radEi_BJ2->push_back(radEi.BJ2()); DEBUG(0);
      _jStr_radEi_BJ3->push_back(radEi.BJ3()); DEBUG(0);
      _jStr_radEi_girthJ1->push_back(radEi.girthJ1()); DEBUG(0);
      _jStr_radEi_girthJ2->push_back(radEi.girthJ2()); DEBUG(0);
      _jStr_radEi_girthJ3->push_back(radEi.girthJ3()); DEBUG(0);
      _jStr_radEi_girth32->push_back(radEi.girth32()); DEBUG(0);
      _jStr_radEi_girth21->push_back(radEi.girth21()); DEBUG(0);
      _jStr_radEi_girthAsymJ1J2->push_back(radEi.girthAsymJ1J2()); DEBUG(0);
      _jStr_radEi_alpha1->push_back(radEi.alpha1()); DEBUG(0);
      _jStr_radEi_alpha2->push_back(radEi.alpha2()); DEBUG(0);
      _jStr_radEi_alpha->push_back(radEi.alpha()); DEBUG(0);
      _jStr_radEi_beta1->push_back(radEi.beta1()); DEBUG(0);
      _jStr_radEi_beta2->push_back(radEi.beta2()); DEBUG(0);
      _jStr_radEi_beta->push_back(radEi.beta()); DEBUG(0);
      _jStr_radEi_thetaJ1J2->push_back(radEi.thetaJ1J2()); DEBUG(0);
      _jStr_radEi_dipolarityInfLine->push_back(radEi.dipolarityInfLine()); DEBUG(0);
      _jStr_radEi_dipolarityLineSeg->push_back(radEi.dipolarityLineSeg()); DEBUG(0);
      _jStr_radEi_BCF1->push_back(radEi.BCF1()); DEBUG(0);
      _jStr_radEi_BCF2->push_back(radEi.BCF2()); DEBUG(0);
      _jStr_radEi_BCF3->push_back(radEi.BCF3()); DEBUG(0);
      _jStr_radEi_dipolarity->push_back(radEi.dipolarity()); DEBUG(0);
      _jStr_sphEi_detSphericity->push_back(sphEi.detSphericity()); DEBUG(0);
      _jStr_sphEi_detSpherocity->push_back(sphEi.detSpherocity()); DEBUG(0);
      _jStr_sphEi_sphericityLambda1->push_back(sphEi.sphericityLambda1()); DEBUG(0);
      _jStr_sphEi_sphericityLambda2->push_back(sphEi.sphericityLambda2()); DEBUG(0);
      _jStr_sphEi_sphericityLambda3->push_back(sphEi.sphericityLambda3()); DEBUG(0);
      _jStr_sphEi_spherocityLambda1->push_back(sphEi.spherocityLambda1()); DEBUG(0);
      _jStr_sphEi_spherocityLambda2->push_back(sphEi.spherocityLambda2()); DEBUG(0);
      _jStr_sphEi_spherocityLambda3->push_back(sphEi.spherocityLambda3()); DEBUG(0);
      _jStr_sphEi_circularity->push_back(sphEi.circularity()); DEBUG(0);
      _jStr_sphEi_sphericity->push_back(sphEi.sphericity()); DEBUG(0);
      _jStr_sphEi_spherocity->push_back(sphEi.spherocity()); DEBUG(0);
      _jStr_sphEi_aplanarity->push_back(sphEi.aplanarity()); DEBUG(0);
      _jStr_sphEi_aplanority->push_back(sphEi.aplanority()); DEBUG(0);
      _jStr_sphEi_Y->push_back(sphEi.Y()); DEBUG(0);
      _jStr_sphEi_planarity->push_back(sphEi.planarity()); DEBUG(0);
      _jStr_sphEi_planority->push_back(sphEi.planority()); DEBUG(0);
      _jStr_sphEi_Dshape->push_back(sphEi.Dshape()); DEBUG(0);
      _jStr_sphEi_Cshape->push_back(sphEi.Cshape()); DEBUG(0);
      _jStr_sphEi_H2->push_back(sphEi.H2()); DEBUG(0);
      _jStr_sphEi_Fshape->push_back(sphEi.Fshape()); DEBUG(0);
      _jStr_sphEi_beamThrust->push_back(sphEi.beamThrust()); DEBUG(0);
      _jStr_sphEi_G->push_back(sphEi.G()); DEBUG(0);
      _jStr_sphEi_ST2D->push_back(sphEi.ST2D()); DEBUG(0);
      _jStr_sphEi_detMlin->push_back(sphEi.detMlin()); DEBUG(0);
      _jStr_sphEi_pft->push_back(sphEi.pft()); DEBUG(0);
    
      ATH_MSG_INFO("Clean-up jet collections...");
      //Deleter()(jetColl); DEBUG(0); // No longer needed if using new CoreSubJetsTool
    } // loop over different collections
    DEBUG(0);
    m_tree_AS->Fill(); DEBUG(0);

  } // try block; my first use of exceptions!

  catch(const std::exception& e) {
    ATH_MSG_ERROR("Exception thrown:\n" << e.what() << "\n" << "at line " << __LINE__ << " in file " << __FILE__ << " in function " << __FUNCTION__);
  }
  catch(...) {
    // unknown exception, should not happen
    ATH_MSG_ERROR("An exception occurred!");
  }
  
  // CALLGRIND_STOP_INSTRUMENTATION;
  // feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);

  _NFPE = fetestexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
  feclearexcept(FE_INVALID | FE_DIVBYZERO| FE_OVERFLOW);
  return StatusCode::SUCCESS;
}

int MyAnalysis::getQuarkJetFlavour(JetCollection::const_iterator jetItr) {
  /** flavour of quark that originated this jet */
  // --- get the true label of the jet from MC Truth

  std::string label("N/A");
    
  const Analysis::TruthInfo* mcinfo = (*jetItr)->tagInfo<Analysis::TruthInfo>("TruthInfo");
  if(mcinfo) {
    label = mcinfo->jetTruthLabel();
  }
  else {
    ATH_MSG_INFO("Could not find TruthInfo for jet");
  }
  int iflav(0);
  if(label == "B") {
    return iflav = 5;
  }
  if(label == "C") {
    return iflav = 4;
  }
  if(label == "T") {
    return iflav = 15;
  }
  return iflav;
}

StatusCode MyAnalysis::addEventInfo() {
  ATH_MSG_DEBUG( "in addEventInfo");
  //get EventInfo for run and event number

  const EventInfo* eventInfo;
  StatusCode sc = evtStore()->retrieve(eventInfo); DEBUG(eventInfo);
  
  if(sc.isFailure()) {
    ATH_MSG_WARNING("Could not retrieve event info");
    return sc;
  }
  
  const EventID* myEventID = eventInfo->event_ID(); DEBUG(myEventID);

  m_runNumber = myEventID->run_number();
  m_eventNumber = myEventID->event_number();
  ATH_MSG_DEBUG("event " << m_eventNumber);

  m_eventTime = myEventID->time_stamp(); 
  m_lumiBlock = myEventID->lumi_block();
  m_bCID = myEventID->bunch_crossing_id();

  const EventType* myEventType = eventInfo->event_type(); DEBUG(myEventType);
  if(myEventType) {
    m_eventWeight = myEventType->mc_event_weight();
  }
  else {
    m_eventWeight = -999;
  }
  return StatusCode::SUCCESS;
}

// http://cdsweb.cern.ch/record/1342550/files/ATLAS-CONF-2011-053.pdf
// http://pdg.lbl.gov/2011/reviews/rpp2011-rev-monte-carlo-numbering.pdf 
long long int MyAnalysis::QGtag(const TruthParticleContainer* mcpartTES, const Jet* j) {
  int result = -1;
  if(mcpartTES && j) {
    double e1 = 0.;
    double dR1 = -1.;    
    RadialParameter R;
    double radialParameter = 0.1 * R(j);

    TruthParticleContainer::const_iterator mcpItr  = mcpartTES->begin();
    TruthParticleContainer::const_iterator mcpItrE = mcpartTES->end();
    int mothid = -1;
    int gmothid = -1;
    int ggmothid = -1;
    int gggmothid = -1;
    int ggggmothid = -1;

    int rmothid = -1;
    int rgmothid = -1;
    int rggmothid = -1;
    int rgggmothid = -1;
    int rggggmothid = -1;

    for(; mcpItr != mcpItrE; ++mcpItr) {
      const HepMC::GenParticle* part =(*mcpItr)->genParticle();

      const TruthParticle* tpmoth = (*mcpItr)->mother();
      if(tpmoth) {
	const HepMC::GenParticle* gpmoth = tpmoth->genParticle();
	if(gpmoth) mothid = abs(gpmoth->pdg_id());
	const TruthParticle* tpgmoth = tpmoth->mother();
	if(tpgmoth) {
	  const HepMC::GenParticle* gpgmoth = tpgmoth->genParticle();
	  if(gpgmoth) gmothid = abs(gpgmoth->pdg_id());
	  const TruthParticle* tpggmoth = tpgmoth->mother();
	  if(tpggmoth) {
	    const HepMC::GenParticle* gpggmoth = tpggmoth->genParticle();
	    if(gpggmoth) ggmothid = abs(gpggmoth->pdg_id());
	    const TruthParticle* tpgggmoth = tpggmoth->mother();
	    if(tpgggmoth) {
	      const HepMC::GenParticle* gpgggmoth = tpgggmoth->genParticle();
	      if(gpgggmoth) gggmothid = abs(gpgggmoth->pdg_id());
	      const TruthParticle* tpggggmoth = tpgggmoth->mother();
	      if(tpggggmoth) {
		const HepMC::GenParticle* gpggggmoth = tpggggmoth->genParticle();
		if(gpggggmoth) ggggmothid = abs(gpggggmoth->pdg_id());
	      }
	    }
	  }
	}
      }
      
      int pdgid = abs(part->pdg_id());

      if(pdgid == 1 || // d
	 pdgid == 2 || // u
	 pdgid == 3 || // s
	 pdgid == 4 || // c
	 pdgid == 5 || // b
	 pdgid == 6 || // t
	 pdgid == 15 || // tau
	 pdgid == 21) { // g
	
	double eta1 = part->momentum().eta();
	double phi1 = part->momentum().phi();
	double eta2 = j->eta();
	double phi2 = j->phi();
	double dEta = fabs(eta1 - eta2);
	double dPhi = DeltaPhi()(phi1, phi2);
	double dR = sqrt((dEta * dEta) + (dPhi * dPhi));

	if(dR < radialParameter) {
	  if(part->momentum().e() > e1) {
	    e1 = part->momentum().e();
	    result = pdgid;
	    rmothid = mothid;
	    rgmothid = gmothid;
	    rggmothid = ggmothid;
	    rgggmothid = gggmothid;
	    rggggmothid = ggggmothid;
	    dR1 = dR;
	  } 
	} // inside cone
      } // have a parton
    } // loop over truth particles

    std::cout << " result so far: " << result;

    if(result > 0) std::cout << " part: " << result
			     << " mother: " << rmothid
			     << " grandmother: " << rgmothid
			     << " great grandmother: " << rggmothid
			     << " great great grandmother: " << rgggmothid
			     << " great great great grandmother: " << rggggmothid;

    if(result < 0) { // We didn't succeed to find a parton
      long long int tempResult = 0;
      e1 = 0.;
      dR1 = -1.;  
      mcpItr  = mcpartTES->begin();
      mcpItrE = mcpartTES->end();
      
      for(; mcpItr != mcpItrE; ++mcpItr) {
	const HepMC::GenParticle* part =(*mcpItr)->genParticle();
	int pdgid = abs(part->pdg_id());

	if(pdgid <= 25) { // What we're willing to accept

	  double eta1 = part->momentum().eta();
	  double phi1 = part->momentum().phi();
	  double eta2 = j->eta();
	  double phi2 = j->phi();
	  double dEta = fabs(eta1 - eta2);
	  double dPhi = DeltaPhi()(phi1, phi2);
	  double dR = sqrt((dEta * dEta) + (dPhi * dPhi));
	  
	  if(dR < radialParameter) {
	    if(part->momentum().e() > e1) {
	      e1 = part->momentum().e();
	      tempResult = pdgid*1000000;
	      dR1 = dR;
	    } 
	  } // inside cone
	} // have a parton
      } // loop over truth particles
      if(tempResult == 0) tempResult = 99*1000000;
      result = tempResult;

      tempResult = 0;
      e1 = 0.;
      dR1 = -1.;  
      mcpItr  = mcpartTES->begin();
      mcpItrE = mcpartTES->end();
      
      for(; mcpItr != mcpItrE; ++mcpItr) {
	const HepMC::GenParticle* part =(*mcpItr)->genParticle();
	int pdgid = abs(part->pdg_id());
	
	if(pdgid < 25) { // What we're willing to accept
	  
	  double eta1 = part->momentum().eta();
	  double phi1 = part->momentum().phi();
	  double eta2 = j->eta();
	  double phi2 = j->phi();
	  double dEta = fabs(eta1 - eta2);
	  double dPhi = DeltaPhi()(phi1, phi2);
	  double dR = sqrt((dEta * dEta) + (dPhi * dPhi));
	  
	  if(dR < radialParameter) {
	    if(part->momentum().e() > e1) {
	      e1 = part->momentum().e();
	      tempResult = pdgid*10000;
	      dR1 = dR;
	    } 
	  } // inside cone
	} // have a parton
      } // loop over truth particles
      if(tempResult == 0) tempResult = 99*10000;
      result += tempResult;

      tempResult = 0;
      e1 = 0.;
      dR1 = -1.;  
      mcpItr  = mcpartTES->begin();
      mcpItrE = mcpartTES->end();
      
      for(; mcpItr != mcpItrE; ++mcpItr) {
	const HepMC::GenParticle* part =(*mcpItr)->genParticle();
	int pdgid = abs(part->pdg_id());
	
	if(pdgid < 24) { // What we're willing to accept
	  
	  double eta1 = part->momentum().eta();
	  double phi1 = part->momentum().phi();
	  double eta2 = j->eta();
	  double phi2 = j->phi();
	  double dEta = fabs(eta1 - eta2);
	  double dPhi = DeltaPhi()(phi1, phi2);
	  double dR = sqrt((dEta * dEta) + (dPhi * dPhi));
	  
	  if(dR < radialParameter) {
	    if(part->momentum().e() > e1) {
	      e1 = part->momentum().e();
	      tempResult = pdgid*100;
	      dR1 = dR;
	    } 
	  } // inside cone
	} // have a parton
      } // loop over truth particles
      if(tempResult == 0) tempResult = 99*100;
      result += tempResult;

      tempResult = 0;
      e1 = 0.;
      dR1 = -1.;  
      mcpItr  = mcpartTES->begin();
      mcpItrE = mcpartTES->end();
      
      for(; mcpItr != mcpItrE; ++mcpItr) {
	const HepMC::GenParticle* part =(*mcpItr)->genParticle();
	
	const TruthParticle* tpmoth = (*mcpItr)->mother();
	if(tpmoth) {
	  const HepMC::GenParticle* gpmoth = tpmoth->genParticle();
	  if(gpmoth) mothid = abs(gpmoth->pdg_id());
	  const TruthParticle* tpgmoth = tpmoth->mother();
	  if(tpgmoth) {
	    const HepMC::GenParticle* gpgmoth = tpgmoth->genParticle();
	    if(gpgmoth) gmothid = abs(gpgmoth->pdg_id());
	    const TruthParticle* tpggmoth = tpgmoth->mother();
	    if(tpggmoth) {
	      const HepMC::GenParticle* gpggmoth = tpggmoth->genParticle();
	      if(gpggmoth) ggmothid = abs(gpggmoth->pdg_id());
	      const TruthParticle* tpgggmoth = tpggmoth->mother();
	      if(tpgggmoth) {
		const HepMC::GenParticle* gpgggmoth = tpgggmoth->genParticle();
		if(gpgggmoth) gggmothid = abs(gpgggmoth->pdg_id());
		const TruthParticle* tpggggmoth = tpgggmoth->mother();
		if(tpggggmoth) {
		  const HepMC::GenParticle* gpggggmoth = tpggggmoth->genParticle();
		  if(gpggggmoth) ggggmothid = abs(gpggggmoth->pdg_id());
		}
	      }
	    }
	  }
	}
	
	int pdgid = abs(part->pdg_id());
	
	if(pdgid < 23 &&
	   pdgid != 12 &&
	   pdgid != 14 &&
	   pdgid != 16) { // What we're willing to accept
	  
	  double eta1 = part->momentum().eta();
	  double phi1 = part->momentum().phi();
	  double eta2 = j->eta();
	  double phi2 = j->phi();
	  double dEta = fabs(eta1 - eta2);
	  double dPhi = DeltaPhi()(phi1, phi2);
	  double dR = sqrt((dEta * dEta) + (dPhi * dPhi));
	  
	  if(dR < radialParameter) {
	    if(part->momentum().e() > e1) {
	      e1 = part->momentum().e();
	      tempResult = pdgid*1;
	      rmothid = mothid;
	      rgmothid = gmothid;
	      rggmothid = ggmothid;
	      rgggmothid = gggmothid;
	      rggggmothid = ggggmothid;
	      dR1 = dR;
	    } 
	  } // inside cone
	} // have a parton
      } // loop over truth particles



      if(tempResult == 0) tempResult = 99*1;
      result += tempResult;

      result *= -1;

      std::cout << " part: " << tempResult
		<< " mother: " << rmothid
		<< " grandmother: " << rgmothid
		<< " great grandmother: " << rggmothid
		<< " great great grandmother: " << rgggmothid
		<< " great great great grandmother: " << rggggmothid;
    }
  } // have truth info

  std::cout << " FINAL RESULT: " << result << std::endl;
  return result;
}

long long int MyAnalysis::QGtag2(const TruthParticleContainer* mcpartTES, const Jet* j) {
  int result = -1;
  if(mcpartTES && j) {
    double e1 = 0.;
    double dR1 = -1.;    
    RadialParameter R;
    double radialParameter = 0.1 * R(j);

    TruthParticleContainer::const_iterator mcpItr  = mcpartTES->begin();
    TruthParticleContainer::const_iterator mcpItrE = mcpartTES->end();

    for(; mcpItr != mcpItrE; ++mcpItr) {
      const HepMC::GenParticle* part =(*mcpItr)->genParticle();
      int pdgid = abs(part->pdg_id());

      if(pdgid == 1 || // d
	 pdgid == 2 || // u
	 pdgid == 3 || // s
	 pdgid == 4 || // c
	 pdgid == 5 || // b
	 pdgid == 6 || // t
	 pdgid == 11 || // e
	 pdgid == 13 || // mu
	 pdgid == 15 || // tau
	 pdgid == 21 || // g
	 pdgid == 22) { // gamma
	
	double eta1 = part->momentum().eta();
	double phi1 = part->momentum().phi();
	double eta2 = j->eta();
	double phi2 = j->phi();
	double dEta = fabs(eta1 - eta2);
	double dPhi = DeltaPhi()(phi1, phi2);
	double dR = sqrt((dEta * dEta) + (dPhi * dPhi));

	if(dR < radialParameter) {
	  if(part->momentum().e() > e1) {
	    e1 = part->momentum().e();
	    result = pdgid;
	    dR1 = dR;
	  } 
	} // inside cone
      } // have a parton
    } // loop over truth particles
  } // have truth info
  return result;
}

