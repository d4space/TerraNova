import FWCore.ParameterSet.Config as cms

process = cms.Process("TerraTuple")

process.load("TerraNova.NtupleMaker.Common_cfg")
process.load("TerraNova.NtupleMaker.NtupleMaker_MC_cff")

process.TTsemiLept.Channel = cms.untracked.string("Muon")
process.TTsemiLept.leptonLabel = cms.InputTag("Muons")


from TerraNova.NtupleMaker.pat_22Jan2013_MC_cfg import *
process.GlobalTag.globaltag = myGlobaltag

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = 1000


#process.load("FWCore.MessageLogger.MessageLogger_cfi")
#process.MessageLogger = cms.Service("MessageLogger")
#process.MessageLogger.destinations = ['cout']
#process.MessageLogger.cout = cms.untracked.PSet(
#    threshold = cms.untracked.string('INFO'),
#    FwkReport = cms.untracked.PSet(reportEvery=cms.untracked.int32(1000))
#)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring("file:../../prod/Proto/patTuple_skim.root")
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('TerraTuple.root')
)

#process.load("TerraNova.NtupleMaker.Sources.DYToMuMu_S8_Skim_cff")
#process.load("TerraNova.NtupleMaker.Sources.PatSkimTemplate_cff")
#process.load("TerraNova.NtupleMaker.Sources.WplusToMuNu_S8_8TeV_AODSIM_PatSkim_local_cff")

process.p = cms.Path(
    process.TTsemiLeptMuMCSequence
    #process.WMuNeuAnalysisMCSequence
)

