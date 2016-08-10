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

#process.load("Configuration.StandardSequences.GeometryConf")

#needs to be commented out or else it complains about: 'two EventSetup Producers want to deliver type="CaloSubdetectorGeometry" label="ZDC"'

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

#process.load("Geometry.CMSCommonData.cmsIdealGeometryXML_cfi");
#process.load("Geometry.CaloEventSetup.CaloGeometry_cfi");
#process.load("Geometry.CaloEventSetup.CaloTopology_cfi");
#process.load("Configuration.Geometry.GeometryECALHCAL_cff")

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '80X_dataRun2_Prompt_v8', '')
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_data', '')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
       # 'root://cms-xrd-global.cern.ch///store/data/Run2016B/SingleMuon/MINIAOD/PromptReco-v2/000/273/150/00000/34A57FB8-D819-E611-B0A4-02163E0144EE.root'
        '/store/mc/RunIISpring16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/50000/000FF6AC-9F2A-E611-A063-0CC47A4C8EB0.root' # DY MC
#        '/store/data/Run2016B/SingleMuon/MINIAOD/PromptReco-v2/000/273/150/00000/34A57FB8-D819-E611-B0A4-02163E0144EE.root' # 2016B DATA
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
                                 HLTPath1 =  cms.string( "HLT_IsoMu18_v2" ),#"HLT_IsoMu20_v4"),  #HLT_IsoMu24_eta2p1_v1" ),
                                 HLTFilter1a= cms.string( "hltL3crIsoL1sSingleMu16erL1f0L2f10QL3f17QL3trkIsoFiltered0p09"),#hltL3crIsoL1sMu16L1f0L2f10QL3f18QL3trkIsoFiltered0p09"),
                                 HLTFilter1b= cms.string( "hltL3crIsoL1sSingleMu16erL1f0L2f10QL3f17QL3trkIsoFiltered0p09"),#hltL3crIsoL1sMu16L1f0L2f10QL3f18QL3trkIsoFiltered0p09"),
                                 HLTPath2 =  cms.string( "HLT_IsoMu19_eta2p1_MediumIsoPFTau32_Trk1_eta2p1_Reg_v2"),#HLT_IsoMu19_eta2p1_MediumIsoPFTau32_Trk1_eta2p1_Reg_v1)",#HLT_IsoMu19_eta2p1_MediumIsoPFTau32_Trk1_eta2p1_Reg_v2" ), #"HLT_IsoMu17_eta2p1_MediumIsoPFTau40_Trk1_eta2p1_Reg_v3"),
                                 HLTFilter2a = cms.string( "hltL3crIsoL1sMu18erIsoTau26erL1f0L2f10QL3f19QL3trkIsoFiltered0p09"),#hltL3crIsoL1sMu16erTauJet20erL1f0L2f10QL3f17QL3trkIsoFiltered0p09" ),
                                 HLTFilter2b = cms.string( "hltOverlapFilterIsoMu19L2IsoTau26"),#hltOverlapFilterIsoMu17MediumIsoPFTau40Reg" ),

                                 # samples input tag
                                 isMC = cms.bool(True),
                                 isZtau  = cms.bool(True),
                                 isZprime  = cms.bool(False),
                                 GenReq = cms.bool(False),

                                 # objects input tags 
                                 muons               = cms.InputTag("slimmedMuons"),
                                 jets                = cms.InputTag("slimmedJets"),
                                 taus                = cms.InputTag("slimmedTaus"),
                                 electrons           = cms.InputTag("slimmedElectrons"),
                                 vertices            = cms.InputTag("offlineSlimmedPrimaryVertices"),
                                 mets                = cms.InputTag("slimmedMETs"),

                                 #Trigger Input Tag
                                 triggerBits                = cms.InputTag("TriggerResults","","HLT"),
                                 trigger_prescale           = cms.InputTag("patTrigger"),
                                 triggerObjects             = cms.InputTag("selectedPatTrigger"),			  
                                 l1min = cms.InputTag("patTrigger","l1min"),
                                 l1max = cms.InputTag("patTrigger","l1max"),           
                                 
                                 #MC
                                 genparts= cms.InputTag("prunedGenParticles"),
                                 genweights= cms.InputTag("generator"),
                                 
                                 #Cut values
                                 isOSCharge = cms.bool(True),     
                                 TauPtCut = cms.double(16.0),
                                 TauEtaCut = cms.double(2.1),
                                 TauDMF = cms.string("decayModeFinding"),
                                 TauEleVeto = cms.string("againstElectronVLooseMVA6"),
                                 TauMuVeto = cms.string("againstMuonTight3"),
                                 TauIsoString = cms.string("byLooseCombinedIsolationDeltaBetaCorr3Hits"),
                                 DYOthersBG = cms.bool(False),
                                 TauIsoCutMax = cms.double(99999),
                                 TauIsoCutMin = cms.double(0.5),
                                 MuonPtCut = cms.double(18.0),
                                 MuonEtaCut = cms.double(2.1),
                                 IsoMuonMax = cms.double(0.1),
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


