import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.load('Configuration.Geometry.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')

process.load('TrackingTools.TransientTrack.TransientTrackBuilder_cfi')
process.load("Configuration.StandardSequences.MagneticField_38T_PostLS1_cff")

#needs to be commented out or else it complains about: 'two EventSetup Producers want to deliver type="CaloSubdetectorGeometry" label="ZDC"'

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '80X_dataRun2_Prompt_v8', '')
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_data', '')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100000) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
       # 'root://cms-xrd-global.cern.ch///store/data/Run2016B/SingleMuon/MINIAOD/PromptReco-v2/000/273/150/00000/34A57FB8-D819-E611-B0A4-02163E0144EE.root'
        '/store/data/Run2016B/SingleMuon/MINIAOD/PromptReco-v2/000/273/150/00000/34A57FB8-D819-E611-B0A4-02163E0144EE.root'
    )
)

from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
dataFormat = DataFormat.MiniAOD
switchOnVIDElectronIdProducer(process, dataFormat)
my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Spring15_25ns_V1_cff',
                 'RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV60_cff',
                 'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring15_25ns_nonTrig_V1_cff']
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)

process.TFileService = cms.Service("TFileService",
  fileName = cms.string("OutTree.root")
)


process.miniAOD = cms.EDAnalyzer('MiniAODAnalyzer',

                                 electronVetoIdMap   = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-veto"),
                                 electronLooseIdMap  = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-loose"),
                                 electronMediumIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-medium"),
                                 electronTightIdMap  = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-tight"),
                                 eleHEEPIdMap        = cms.InputTag("egmGsfElectronIDs:heepElectronID-HEEPV60"),

                                 #Trigger Paths
#                                 HLTPath1 =  cms.string( "HLT_IsoMu20_v4"),  #HLT_IsoMu24_eta2p1_v1" ),
         
                                 # input tags 
                                 muons               = cms.InputTag("slimmedMuons"),
                                 jets                = cms.InputTag("slimmedJets"),
                                 taus                = cms.InputTag("slimmedTaus"),
                                 electrons           = cms.InputTag("slimmedElectrons"),
                                 vertices            = cms.InputTag("offlineSlimmedPrimaryVertices"),
                                 mets                = cms.InputTag("slimmedMETs"),

                                 triggerBits                = cms.InputTag("TriggerResults","","HLT"),
                                 trigger_prescale           = cms.InputTag("patTrigger"),
                                 triggerObjects             = cms.InputTag("selectedPatTrigger"),			  
                                 l1min = cms.InputTag("patTrigger","l1min"),
                                 l1max = cms.InputTag("patTrigger","l1max")           
  

)


process.load('CommonTools.RecoAlgos.HBHENoiseFilterResultProducer_cfi')
process.HBHENoiseFilterResultProducer.minZeros = cms.int32(99999)
process.HBHENoiseFilterResultProducer.IgnoreTS4TS5ifJetInLowBVRegion=cms.bool(False) 
process.HBHENoiseFilterResultProducer.defaultDecision = cms.string("HBHENoiseFilterResultRun2Loose")
                              
process.ApplyBaselineHBHENoiseFilter = cms.EDFilter('BooleanFlagFilter',
   inputLabel = cms.InputTag('HBHENoiseFilterResultProducer','HBHENoiseFilterResult'),
   reverseDecision = cms.bool(False)
)                             
                              
process.ApplyBaselineHBHEIsoNoiseFilter = cms.EDFilter('BooleanFlagFilter',
   inputLabel = cms.InputTag('HBHENoiseFilterResultProducer','HBHEIsoNoiseFilterResult'),
   reverseDecision = cms.bool(False)
)                             

process.load("RecoMET/METProducers.METSignificance_cfi")
process.load("RecoMET/METProducers.METSignificanceParams_cfi")

process.p = cms.Path(
#                     process.HBHENoiseFilterResultProducer * 
#                     process.ApplyBaselineHBHENoiseFilter *    
                     process.egmGsfElectronIDSequence *
                     process.miniAOD
)


