#include "DataFormats.h"

bool isKinTight(_MuonInfo& muon) {

  bool isKinTight=false;

 
   if (!muon.isTracker) return isKinTight;
   if (!muon.isGlobal)  return isKinTight;
  
  // acceptance cuts
  if (muon.pt < 20)         return isKinTight; // pt cut
  if (fabs(muon.eta) > 2.1) return isKinTight; // eta cut

  // kinematic cuts
  if (muon.numTrackerLayers < 9) return isKinTight; // # hits in tracker

  // iso cut
  if (muon.trackIsoSumPt/muon.pt>0.1) return isKinTight;

  // beam spot cut
  //if (fabs(muon.d0) > 0.2) return isKinTight;
  if (fabs(muon.d0_PV) > 0.2) return isKinTight;

  if ( muon.numValidMuonHits  < 1 ) return isKinTight;
  if ( muon.numValidPixelHits < 1 ) return isKinTight;
  if ( muon.numOfMatchedStations < 2 ) return isKinTight;

  isKinTight=true;
  return isKinTight;
}

bool isKinTight_2012(_MuonInfo& muon) {

  bool isKinTight_2012=false;

  
  if (!muon.isGlobal)  return isKinTight_2012;
  if (!muon.isPFMuon) return isKinTight_2012;

  // acceptance cuts
  if (muon.pt < 20)         return isKinTight_2012; // pt cut
  if (fabs(muon.eta) > 2.1) return isKinTight_2012; // eta cut

  // kinematic cuts
  if (muon.numTrackerLayers < 6) return isKinTight_2012; // # hits in tracker

  // iso cut
  if (muon.trackIsoSumPt/muon.pt>0.1) return isKinTight_2012;

  // beam spot cut
   if (fabs(muon.d0_PV) > 0.2) return isKinTight_2012;
   if (fabs(muon.dz_PV) > 0.5) return isKinTight_2012;
//   if (fabs(muon.d0) > 0.2) return isKinTight_2012;

  if ( muon.numValidMuonHits  < 1 ) return isKinTight_2012;
  if ( muon.numValidPixelHits < 1 ) return isKinTight_2012;
  if ( muon.numOfMatchedStations < 2 ) return isKinTight_2012;
  if ( muon.normChiSquare > 10)     return isKinTight_2012;

  isKinTight_2012=true;
  return isKinTight_2012;
}

  
bool isKinTight_2012WZ(_MuonInfo& muon) {

  bool isKinTight_2012WZ=false;


  if (!muon.isGlobal)  return isKinTight_2012WZ;
  if (!muon.isPFMuon) return isKinTight_2012WZ;

  // acceptance cuts
  if (muon.pt < 25)         return isKinTight_2012WZ; // pt cut
  if (fabs(muon.eta) > 2.1) return isKinTight_2012WZ; // eta cut

  // kinematic cuts
  if (muon.numTrackerLayers < 6) return isKinTight_2012WZ; // # hits in tracker

  // PF iso cut with deltaB correction
    if ( ( muon.sumChargedHadronPtR04 + max(0.0,muon.sumNeutralHadronEtR04 + muon.sumPhotonEtR04 - 0.5*muon.sumPUPtR04) )/muon.pt > 0.12) return isKinTight_2012WZ;
 
//   if (muon.trackIsoSumPt/muon.pt>0.1) return isKinTight_2012WZ;

  // beam spot cut
  if (fabs(muon.d0_PV) > 0.02) return isKinTight_2012WZ;
  if (fabs(muon.dz_PV) > 0.5) return isKinTight_2012WZ;
  
  if ( muon.numValidMuonHits  < 1 ) return isKinTight_2012WZ;
  if ( muon.numValidPixelHits < 1 ) return isKinTight_2012WZ;
  if ( muon.numOfMatchedStations < 2 ) return isKinTight_2012WZ;
  if ( muon.normChiSquare > 10)     return isKinTight_2012WZ;

  isKinTight_2012WZ=true;
  return isKinTight_2012WZ;
}


// of course pt and eta cuts have to be adjusted to the
// HLT requirements
bool isKinTight_2011A(_MuonInfo& muon) {
  bool isKinTight=false;

  if (!muon.isTracker) return isKinTight;  
  if (!muon.isGlobal)  return isKinTight;  

  // acceptance cuts
  if (muon.pt < 45)         return isKinTight; // pt cut
  if (fabs(muon.eta) > 2.1) return isKinTight; // eta cut
  
  // kinematic cuts
  if (muon.numValidTrackerHits < 11) return isKinTight; // # hits in tracker
  
  // iso cut
  if (muon.trackIsoSumPt>10) return isKinTight;

  // beam spot cut
  //if (fabs(muon.d0) > 0.2) return isKinTight;
  if (fabs(muon.d0_PV) > 0.2) return isKinTight;

  if ( muon.numValidMuonHits  < 1 ) return isKinTight;
  if ( muon.numValidPixelHits < 1 ) return isKinTight;
  if ( muon.numSegmentMatches < 2 ) return isKinTight;
  if ( muon.normChiSquare > 10)     return isKinTight;

  isKinTight=true;
  return isKinTight;
}

bool isKinTight_PF(_MuonInfo& muon) {

  bool isKinTight_PF=false;

//  if (!muon.isTracker) return isKinTight;  
//  if (!muon.isGlobal)  return isKinTight;  

  if (!muon.isPFMuon) return isKinTight_PF;

  // acceptance cuts
  if (muon.pfPt < 20)         return isKinTight_PF; // pt cut
  if (fabs(muon.pfEta) > 2.1) return isKinTight_PF; // eta cut

  // kinematic cuts
  if (muon.numTrackerLayers < 6) return isKinTight_PF; // # hits in tracker

  // iso cut
  if (muon.trackIsoSumPt/muon.pt>0.1) return isKinTight_PF;

  // beam spot cut
  //if (fabs(muon.d0) > 0.2) return isKinTight_PF;
  if (fabs(muon.d0_PV) > 0.2) return isKinTight_PF;

if ( muon.numValidMuonHits  < 1 ) return isKinTight_PF;
  if ( muon.numValidPixelHits < 1 ) return isKinTight_PF;
  if ( muon.numOfMatchedStations < 2 ) return isKinTight_PF;
  if ( muon.normChiSquare > 10)     return isKinTight_PF;

  isKinTight_PF=true;
  return isKinTight_PF;
}

bool isKinMedium(_MuonInfo& muon) {

  bool isKinMedium=false;

  // only being tracker is enough
  if (!muon.isTracker) return isKinMedium;  

  // acceptance cuts
  if (muon.pt < 45)         return isKinMedium; // pt cut
  if (fabs(muon.eta) > 2.1) return isKinMedium; // eta cut
  //if (muon.pt < 35)         return isKinMedium; // pt cut

  // kinematic cuts
  if (muon.numTrackerLayers < 9) return isKinMedium; // # hits in tracker

  // iso cut
  if (muon.trackIsoSumPt/muon.pt>0.1) return isKinMedium;

  // beam spot cut
  //if (fabs(muon.d0) > 0.2) return isKinMedium;
  if (fabs(muon.d0_PV) > 0.2) return isKinMedium;

  if ( muon.numValidPixelHits < 1 ) return isKinMedium;

  isKinMedium=true;
  return isKinMedium;
}

bool nMenOne_pt(_MuonInfo& muon) {

  bool passed=false;

  if (!muon.isPFMuon) return passed;
  if (!muon.isGlobal)  return passed;

  // acceptance cuts
  
  if (fabs(muon.eta) > 2.1) return passed; // eta cut

  // kinematic cuts
  if (muon.numTrackerLayers < 6) return passed; // # hits in tracker

  // iso cut
  if (muon.trackIsoSumPt/muon.pt>0.1) return passed;

 // beam spot cut
  //if (fabs(muon.d0) > 0.2) return passed;
  if (fabs(muon.d0_PV) > 0.2) return passed;

  if ( muon.numValidMuonHits  < 1 ) return passed;
  if ( muon.numValidPixelHits < 1 ) return passed;
  if ( muon.numOfMatchedStations < 2 ) return passed;
  if ( muon.normChiSquare > 10 ) return passed;

  passed=true;
  return passed;
}

bool nMenOne_eta(_MuonInfo& muon) {

  bool passed=false;

  if (!muon.isPFMuon) return passed;
  if (!muon.isGlobal)  return passed;

  // acceptance cuts

   if (muon.pt < 20)         return passed; // pt cut

  // kinematic cuts
  if (muon.numTrackerLayers < 6) return passed; // # hits in tracker

  // iso cut
  if (muon.trackIsoSumPt/muon.pt>0.1) return passed;

 // beam spot cut
  if (fabs(muon.d0_PV) > 0.2) return passed;
  if (fabs(muon.dz_PV) > 0.5) return passed;

  if ( muon.numValidMuonHits  < 1 ) return passed;
  if ( muon.numValidPixelHits < 1 ) return passed;
  if ( muon.numOfMatchedStations < 2 ) return passed;
  if ( muon.normChiSquare > 10 ) return passed;
 passed=true;
  return passed;
}


bool nMenOne_d0(_MuonInfo& muon) {

  bool passed=false;

  if (!muon.isPFMuon) return passed;  
  if (!muon.isGlobal)  return passed;  

  // acceptance cuts
  if (muon.pt < 20)         return passed; // pt cut
  if (fabs(muon.eta) > 2.1) return passed; // eta cut

  // kinematic cuts
  if (muon.numTrackerLayers < 6) return passed; // # hits in tracker

  if (fabs(muon.dz_PV) > 0.5) return passed;

  // iso cut
  if (muon.trackIsoSumPt/muon.pt>0.1) return passed;

  if ( muon.numValidMuonHits  < 1 ) return passed;
  if ( muon.numValidPixelHits < 1 ) return passed;
  if ( muon.numOfMatchedStations < 2 ) return passed;
  if ( muon.normChiSquare > 10 ) return passed;
  passed=true;
  return passed;
}

bool nMenOne_dz(_MuonInfo& muon) {

  bool passed=false;

  if (!muon.isPFMuon) return passed;
  if (!muon.isGlobal)  return passed;

  // acceptance cuts
  if (muon.pt < 20)         return passed; // pt cut
  if (fabs(muon.eta) > 2.1) return passed; // eta cut

  // kinematic cuts
  if (muon.numTrackerLayers < 6) return passed; // # hits in tracker

  if (fabs(muon.d0_PV) > 0.2) return passed;

  // iso cut
  if (muon.trackIsoSumPt/muon.pt>0.1) return passed;

  if ( muon.numValidMuonHits  < 1 ) return passed;
  if ( muon.numValidPixelHits < 1 ) return passed;
  if ( muon.numOfMatchedStations < 2 ) return passed;
  if ( muon.normChiSquare > 10 ) return passed;
  passed=true;
  return passed;
}

bool nMenOne_numValMuHits(_MuonInfo& muon) {

  bool passed=false;

  if (!muon.isPFMuon) return passed;  
  if (!muon.isGlobal)  return passed;  

  // acceptance cuts
  if (muon.pt < 20)         return passed; // pt cut
  if (fabs(muon.eta) > 2.1) return passed; // eta cut

  // kinematic cuts
  if (muon.numTrackerLayers < 6) return passed; // # hits in tracker

  // iso cut
  if (muon.trackIsoSumPt/muon.pt>0.1) return passed;

  // beam spot cut
  if (fabs(muon.d0_PV) > 0.2) return passed;
  if (fabs(muon.dz_PV) > 0.5) return passed;

  if ( muon.numValidPixelHits < 1 ) return passed;
  if ( muon.numOfMatchedStations < 2 ) return passed;
  if ( muon.normChiSquare > 10 ) return passed;

  passed=true;
  return passed;
}

bool nMenOne_numValPxlHits(_MuonInfo& muon) {

  bool passed=false;

  if (!muon.isPFMuon) return passed;  
  if (!muon.isGlobal)  return passed;  

  // acceptance cuts
  if (muon.pt < 20)         return passed; // pt cut
  if (fabs(muon.eta) > 2.1) return passed; // eta cut

  // kinematic cuts
  if (muon.numTrackerLayers < 6) return passed; // # hits in tracker

  // iso cut
  if (muon.trackIsoSumPt/muon.pt>0.1) return passed;

  // beam spot cut
  if (fabs(muon.d0_PV) > 0.2) return passed;
  if (fabs(muon.dz_PV) > 0.5) return passed;


  if ( muon.numValidMuonHits  < 1 ) return passed;
  if ( muon.numOfMatchedStations < 2 ) return passed;
  if ( muon.normChiSquare > 10 ) return passed;

  passed=true;
  return passed;
}

bool nMenOne_numStatMtc(_MuonInfo& muon) {

  bool passed=false;

  if (!muon.isPFMuon) return passed;  
  if (!muon.isGlobal)  return passed;  

  // acceptance cuts
  if (muon.pt < 20)         return passed; // pt cut
  if (fabs(muon.eta) > 2.1) return passed; // eta cut

  // kinematic cuts
  if (muon.numTrackerLayers < 6) return passed; // # hits in tracker

  // iso cut
  if (muon.trackIsoSumPt/muon.pt>0.1) return passed;

  // beam spot cut
  if (fabs(muon.d0_PV) > 0.2) return passed;
  if (fabs(muon.dz_PV) > 0.5) return passed;


  if ( muon.numValidMuonHits  < 1 ) return passed;
  if ( muon.numValidPixelHits < 1 ) return passed;
  if ( muon.normChiSquare > 10 ) return passed;


  passed=true;
  return passed;
}


bool nMenOne_numTrkLys(_MuonInfo& muon) {

  bool passed=false;

  if (!muon.isPFMuon) return passed;  
  if (!muon.isGlobal)  return passed;  

  // acceptance cuts
  if (muon.pt < 20)         return passed; // pt cut
  if (fabs(muon.eta) > 2.1) return passed; // eta cut

  // iso cut
  if (muon.trackIsoSumPt/muon.pt>0.1) return passed;

  // beam spot cut
  if (fabs(muon.d0_PV) > 0.2) return passed;
  if (fabs(muon.dz_PV) > 0.5) return passed;


  if ( muon.numValidMuonHits  < 1 ) return passed;
  if ( muon.numValidPixelHits < 1 ) return passed;
  if ( muon.numOfMatchedStations < 2 ) return passed;
  if ( muon.normChiSquare > 10 ) return passed;

  passed=true;
  return passed;
}

bool nMenOne_RelIso(_MuonInfo& muon) {

  bool passed=false;

  if (!muon.isPFMuon) return passed;  
  if (!muon.isGlobal)  return passed;  

  // acceptance cuts
  if (muon.pt < 20)         return passed; // pt cut
  if (fabs(muon.eta) > 2.1) return passed; // eta cut

  // kinematic cuts
  if (muon.numTrackerLayers < 6) return passed; // # hits in tracker

  // beam spot cut
  
  if (fabs(muon.d0_PV) > 0.2) return passed;
  if (fabs(muon.dz_PV) > 0.5) return passed;

  if ( muon.numValidMuonHits  < 1 ) return passed;
  if ( muon.numValidPixelHits < 1 ) return passed;
  if ( muon.numOfMatchedStations < 2 ) return passed;
  if ( muon.normChiSquare > 10 ) return passed;

  passed=true;
  return passed;
}

bool nMenOne_normChiSq(_MuonInfo& muon) {

  bool passed=false;

  if (!muon.isPFMuon) return passed; 
  if (!muon.isGlobal)  return passed;

  // acceptance cuts
  if (muon.pt < 20)         return passed; // pt cut
  if (fabs(muon.eta) > 2.1) return passed; // eta cut

  // kinematic cuts
  if (muon.numTrackerLayers < 6) return passed; // # hits in tracker

  // iso cut
  if (muon.trackIsoSumPt/muon.pt>0.1) return passed;


  // beam spot cut
  if (fabs(muon.d0_PV) > 0.2) return passed;
  if (fabs(muon.dz_PV) > 0.5) return passed;


  if ( muon.numValidMuonHits  < 1 ) return passed;
  if ( muon.numValidPixelHits < 1 ) return passed;
  if ( muon.numOfMatchedStations < 2 ) return passed;
 
passed=true;
  return passed;
}

