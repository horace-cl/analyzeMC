import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.Config as cms
from glob import glob

process = cms.Process("tracks")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

path1= '/eos/user/h/hcrottel/PrivateMC-2020-b_kmumu_PHSPS/crab_PrivateMC-2020-b_kmumu_PHSPS-2020-09-14-06-57/200914_045749/0000/'
path2= '/eos/user/h/hcrottel/PrivateMC-2020-b_kmumu_PHSPS/crab_PrivateMC-2020-b_kmumu_PHSPS-2020-09-14-07-03/200914_050331/0000/'

files_ = glob(path2+'step0*')


process.source = cms.Source("PoolSource",
                                # replace 'myfile.root' with the source file you want to use
                                fileNames = cms.untracked.vstring(
		['file:'+file_ for file_ in files_]
            #"file:/eos/user/h/hcrottel/PrivateMC-2020-b_kmumu_PHSPS/crab_PrivateMC-2020-b_kmumu_PHSPS-2020-09-10-15-12/200910_131235/0000/step0-GS-b_kmumu_PHSPS_65.root"
	    #'file:/afs/cern.ch/work/h/hcrottel/private/GSnewFilter.root'
	    #'file:/afs/cern.ch/user/h/hcrottel/private/ONLYGEN/GSnewFilterOutput.root'
            #'file:/afs/cern.ch/cms/Tutorials/TWIKI_DATA/TTJets_8TeV_53X.root'
                )
                            )

process.TFileService = cms.Service("TFileService",
#        fileName = cms.string('Rootuple_BstoJpsiphi_2018_MiniAOD.root'),
         fileName = cms.string('GEN-SIM-tracks-07.root'),                                  
)

process.tracks = cms.EDAnalyzer('MCTracks'
                              )

process.p = cms.Path(process.tracks)
