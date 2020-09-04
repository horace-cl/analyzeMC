import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
                                # replace 'myfile.root' with the source file you want to use
                                fileNames = cms.untracked.vstring(
            'file:/afs/cern.ch/user/h/hcrottel/private/ONLYGEN/GS.root'
            #'file:/afs/cern.ch/cms/Tutorials/TWIKI_DATA/TTJets_8TeV_53X.root'
                )
                            )

process.TFileService = cms.Service("TFileService",
#        fileName = cms.string('Rootuple_BstoJpsiphi_2018_MiniAOD.root'),
         fileName = cms.string('Rootuple_ONLYGEN25K.root'),                                  
)

process.demo = cms.EDAnalyzer('MCanalyzer'
                              )

process.p = cms.Path(process.demo)
