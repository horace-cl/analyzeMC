import FWCore.ParameterSet.Config as cms
from glob import glob
import sys

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )


#Source
path1 = '/eos/user/h/hcrottel/PrivateMC-2020-b_kmumu_PHSPS/crab_PrivateMC-2020-b_kmumu_PHSPS-2020-09-14-07-03/200914_050331/0000/'
path2 = '/eos/user/h/hcrottel/PrivateMC-2020-b_kmumu_PHSPS/crab_PrivateMC-2020-b_kmumu_PHSPS-2020-09-14-06-57/200914_045749/0000/'
files = glob(path1+"step3*")
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
		['file:'+file_ for file_ in files[:10]]	
	)
)


#Analyzer
process.MiniSim = cms.EDAnalyzer('MiniAODGenPartAnalyzer',
	packed = cms.InputTag("packedGenParticles"),
	pruned = cms.InputTag("prunedGenParticles")
)



#OUTPUT
print(sys.argv)
if len(sys.argv)>2:
	file_name = sys.argv[2]+'.root'
else:
	file_name = 'MINI-SIM2.root'
print(file_name)
process.TFileService = cms.Service("TFileService",
        fileName = cms.string(file_name),                                  
)



process.p = cms.Path(process.MiniSim)
