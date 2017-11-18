#include "UserAnalysis/RealDataSkeleton.h"

#include <sstream> 
#include <string>
#include <algorithm>
#include <fstream>

#include "GaudiKernel/MsgStream.h" 
#include "GaudiKernel/ISvcLocator.h"
#include "EventInfo/EventInfo.h"
#include "EventInfo/EventID.h"

#include "StoreGate/StoreGateSvc.h" 
#include "StoreGate/DataHandle.h"
#include "AthenaKernel/DefaultKey.h"

#include "VxVertex/VxContainer.h"
#include "TrkV0Vertex/V0Container.h"
#include "TrkVertexAnalysisUtils/V0Tools.h"
#include "muonEvent/MuonContainer.h"

/////////////////////////////////////////////////////////////

RealDataSkeleton::RealDataSkeleton(const std::string& name, ISvcLocator* pSvcLocator) :
   CBNT_AthenaAwareBase(name, pSvcLocator),
   m_trkSelector("InDet::InDetTrackSelectorTool")
{
 
  // Declare user-defined properties from jobOpts - cuts and vertexing methods etc
  declareProperty("TrackSelectorTool", m_trkSelector);
  declareProperty("TrackCollection", m_trackCollName);
  declareProperty("VertexCollection", m_vertexCollName);
  declareProperty("V0Collection", m_v0CollName);
  declareProperty("MuonCollection", m_muonCollectionKey);
 
}

StatusCode RealDataSkeleton::CBNT_initializeBeforeEventLoop() {

  ATH_MSG_DEBUG ("Initializing RealDataSkeleton (before eventloop)" );

  // retrieve trigger decision tool
  // needs to be done before the first run/event since a number of
  // BeginRun/BeginEvents are registered by dependent services

  if ( m_doTrigger ) {
     if ( m_trigDec.retrieve().isFailure() ){
        ATH_MSG_ERROR ( "Can't get handle on TrigDecisionTool" );
     } else {
       ATH_MSG_DEBUG ( "Got handle on TrigDecisionTool" );
     }
  }

  // Initialize the trigger passed counters
  // this can not be done before initialize, since the properties need to be set from job-options first
  std::vector<std::string>::const_iterator it;
  for(it = m_triggerChains.begin();it != m_triggerChains.end(); it++)
     m_triggersPassed[*it] = 0;

  return StatusCode::SUCCESS;
} 

StatusCode RealDataSkeleton::CBNT_clear() {

  m_runNum->clear();
  m_eventNum->clear();
  m_lumiBlock->clear();
  m_timeStamp->clear();

  return StatusCode::SUCCESS;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

StatusCode RealDataSkeleton::CBNT_initialize(){

  ATH_MSG_DEBUG("in CBNT_initialize()");

  // Get the track selector tool from ToolSvc
  if ( m_trkSelector.retrieve().isFailure() ) {
    ATH_MSG_FATAL("Failed to retrieve tool " << m_trkSelector);
    return StatusCode::FAILURE;
  } else {
    ATH_MSG_INFO("Retrieved tool " << m_trkSelector);
  }

  addBranch("RunNumber", m_runNum);
  addBranch("EventNumber", m_eventNum);
  addBranch("LumiBlockr", m_lumiBlock);
  addBranch("TimeStamp", m_timeStamp);

  return StatusCode::SUCCESS;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

StatusCode RealDataSkeleton::CBNT_execute() {

  ATH_MSG_DEBUG("in CBNT_execute()");

// NOW RETRIEVE EVENT INFO, TRACKS, MUONS AND PRIMARY VERTEX FROM STOREGATE

  int eventNum=0;
  int runNum=0;
  int lumiBlock=0;
  int timeStamp=0;

  const EventInfo*  p_evt = 0;
  StatusCode sc = evtStore()->retrieve(p_evt);
  
  if(sc.isSuccess() && p_evt!=0) {
    eventNum = p_evt->event_ID()->event_number();
    runNum = p_evt->event_ID()->run_number();
    lumiBlock = p_evt->event_ID()->lumi_block();
    timeStamp = p_evt->event_ID()->time_stamp();
    ATH_MSG_DEBUG ("Got run number = " << m_runNum << ", event number = " << m_eventNum);
  } else {
    ATH_MSG_WARNING ("Unable to retrieve EventInfo from StoreGate. Return failure.");
    return sc;
  }

  m_eventNum->push_back( eventNum );
  m_runNum->push_back( runNum );
  m_lumiBlock->push_back( lumiBlock );
  m_timeStamp->push_back( timeStamp ); 

  // Get the V0Container from StoreGate
  const V0Container* v0container(0);
  sc = evtStore()->retrieve(v0container, m_v0CollName);
  if (sc.isFailure() || 0 == v0container) {
    ATH_MSG_ERROR("No V0Candidate container " << m_v0CollName << " found in StoreGate!");
    return sc;
  }
  ATH_MSG_DEBUG("You have " << v0container->size() << " V0 candidates in the container");

  // Get the muons from StoreGate
  const Analysis::MuonContainer* importedMuonCollection;
  sc = evtStore()->retrieve(importedMuonCollection,m_muonCollectionKey);
  if(sc.isFailure()){
    ATH_MSG_ERROR("No muon collection with key " << m_muonCollectionKey << " found in StoreGate.");
      return sc;
  } else {
     ATH_MSG_DEBUG("Muon container size "<<importedMuonCollection->size());
  }
 
  return StatusCode::SUCCESS;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

StatusCode RealDataSkeleton::CBNT_finalize() {

  return StatusCode::SUCCESS;
}

