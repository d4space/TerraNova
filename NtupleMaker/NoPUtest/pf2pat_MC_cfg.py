import FWCore.ParameterSet.Config as cms
import os
# Options
produceTaus=True

runOnMC=True
postfix = "PFlow"
jetAlgo="AK5"

produceNoPUPFMET=True
produceMVAPFMET=True

process = cms.Process("PAT")

#PF2PAT
process.load("PhysicsTools.PatAlgos.patSequences_cff")
from PhysicsTools.PatAlgos.tools.pfTools import *
from KoSMP.WAnalyzer.pat_S10_cfg import *
from KoSMP.WAnalyzer.eventContent_cff import *
process.load("KoSMP.WAnalyzer.pf2pat_template_MC_MVaNoPuMEt_cfg")
#from KoSMP.WAnalyzer.tools import *

## Options and Output Report
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
#process.GlobalTag.globaltag = cms.string( 'START53_V18::All' )
process.GlobalTag.globaltag = myGlobaltag

## Source
from PhysicsTools.PatAlgos.tools.cmsswVersionTools import pickRelValInputFiles
#process.source = cms.Source("PoolSource",
#    fileNames = cms.untracked.vstring(
#    pickRelValInputFiles( cmsswVersion  = os.environ['CMSSW_VERSION']
#                        , relVal        = 'RelValTTbar'
#                        , globalTag     = process.GlobalTag.globaltag.value().split(':')[0]
#                        , numberOfFiles = 1
#                        )
#    )
#)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

# output
## Output Module Configuration (expects a path 'p')
process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('patTuple_skim.root'),
    # save only events passing the full path
    SelectEvents   = cms.untracked.PSet( SelectEvents = cms.vstring('p') ),
    # save PAT Layer 1 output; you need a '*' to
    # unpack the list of commands 'patEventContent'
    outputCommands = cms.untracked.vstring('drop *')
)

if runOnMC:
  usePF2PAT(process,runPF2PAT=True, jetAlgo=jetAlgo, runOnMC=runOnMC, postfix=postfix, typeIMetCorrections=True)
else :
  usePF2PAT(process,runPF2PAT=True, jetAlgo=jetAlgo, runOnMC=runOnMC, postfix=postfix, jetCorrections=('AK5PFchs',['L1FastJet','L2Relative', 'L3Absolute','L2L3Residual'] ), typeIMetCorrections=True)


# Pre Settings ######################
# 6(pt>10, tight, for DiMu) 7(pt>27, tight for SingleMu)
# -1(no cut) 0(check cut, isocut pset)
process.acceptedMuons.version = cms.untracked.int32(6)
# -1(no cut), 0(check cut, isocut pset), 1(WptCut)
# 13(medium pt 30) 14(medium pt 15)
process.acceptedElectrons.version = cms.untracked.int32(14)
### taus
process.load("RecoTauTag.Configuration.RecoPFTauTag_cff")

#process.load('RecoJets.JetProducers.PileupJetID_cfi')
#process.pileupJetIdProducer.jets = cms.InputTag('selectedPatJetsPFlow')
#process.pileupJetIdProducer.jets = cms.InputTag('selectedPatJets')

### ============= NoPU and MVA MET ===============###
#if produceMVAPFMET or produceNoPUPFMET :
  # Muon Cuts

  #electron

# No PU MET
process.load('JetMETCorrections.METPUSubtraction.noPileUpPFMET_cff')
if runOnMC:
  process.calibratedAK5PFJetsForNoPileUpPFMEt.correctors = cms.vstring('ak5PFL1FastL2L3')
else:
  process.calibratedAK5PFJetsForNoPileUpPFMEt.correctors = cms.vstring('ak5PFL1FastL2L3Residual')

process.noPileUpPFMEt.srcLeptons = cms.VInputTag(["acceptedMuons","acceptedElectrons","acceptedTaus"])

### MVA MET
process.load('JetMETCorrections.METPUSubtraction.mvaPFMET_cff')
if runOnMC:
  process.calibratedAK5PFJetsForPFMEtMVA.correctors = cms.vstring('ak5PFL1FastL2L3')
else:
  process.calibratedAK5PFJetsForPFMEtMVA.correctors = cms.vstring('ak5PFL1FastL2L3Residual')

process.pfMEtMVA.srcLeptons = cms.VInputTag( ["acceptedMuons","acceptedElectrons","acceptedTaus"]) #selectedPatMuons


#process.pfPileUpIsoPFlow.checkClosestZVertex = cms.bool(False)
#process.pfPileUpIso.checkClosestZVertex = cms.bool(False)

#change cone size
#changeConeSize(process,postfix)

#FastJet!
#applyFastJet(process,postfix)

#REMOVE ISOLATION FROM PF2PAT!!!
#addLooseLeptons(process)

# top projections in PF2PAT:
#topProjection(process,postfix)


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.patMuonFilter.minNumber = 1

## Source
process.load("KoSMP.WAnalyzer.Sources.SourceTemplate_cff")

process.out.outputCommands +=pf2patEventContent

process.outpath = cms.EndPath(process.out)
#process.load("KoPFA.CommonTools.recoPFCandCountFilter_cfi")
process.p = cms.Path(
    process.nEventsTotal*
    process.noscraping*
    process.nEventsNoscrap*
#    process.goodOfflinePrimaryVertices*
    process.HBHENoiseFilter *
    process.nEventsHBHE
#    process.nEventsClean
#    process.patElectronIDs*
#    process.eidCiCSequence
)

if produceTaus:
  process.p += process.recoTauClassicHPSSequence
  print "Using Tau Sequence =================================="
#process.p += process.hltHighLevelMuMuRD
process.p += process.nEventsHLT
process.p += getattr(process,"patPF2PATSequence"+postfix)
process.p += process.acceptedMuons
process.p += process.acceptedElectrons
process.p += process.acceptedTaus
process.p += process.patObjectCountFilter
#process.p += process.patMuonFilter
#process.p += process.patElectronFilter
#process.p += process.patTauFilter
process.p += process.nEventsFiltered
process.p += process.noPileUpPFMEtSequence
process.p += process.pfMEtMVAsequence

