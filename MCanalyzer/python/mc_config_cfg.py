import FWCore.ParameterSet.Config as cms
from glob import glob

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

files = glob("/eos/user/h/hcrottel/PrivateMC-2020-b_kmumu_PHSPS/crab_PrivateMC-2020-b_kmumu_PHSPS-2020-09-14-06-57/200914_045749/0000/step3*")
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
	['file:'+file_ for file_ in files[:10]]
	#"root://cms-xrd-global.cern.ch//store/mc/Spring14miniaod/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/F6EDDC10-8DFC-E311-BC5D-0025905A60D6.root"
	
	)
)

process.demo = cms.EDAnalyzer('MiniAODGenPartAnalyzer',
packed = cms.InputTag("packedGenParticles"),
pruned = cms.InputTag("prunedGenParticles")
)


process.p = cms.Path(process.demo)
