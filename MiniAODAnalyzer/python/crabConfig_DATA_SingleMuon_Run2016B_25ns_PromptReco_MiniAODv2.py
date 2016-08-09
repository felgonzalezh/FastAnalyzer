#from CRABClient.UserUtilities import config, getUsernameFromSiteDB
#config = config()

#config.General.requestName = 'SingleMuon_Run2016B_25ns_PromptReco_MiniAODv2'
#config.General.transferOutputs = True
#config.General.transferLogs = True

#config.JobType.pluginName = 'Analysis'
#config.JobType.psetName = 'ConfFile_cfg.py'
#config.JobType.outputFiles = ['OutTree.root']

#config.Data.inputDataset = '/SingleMuon/Run2016B-PromptReco-v2/MINIAOD'
#config.Data.lumiMask = 'Cert_271036-273730_13TeV_PromptReco_Collisions16_JSON.txt'
#config.Data.runRange = '273150-273292' # '193093-194075'
#config.Data.inputDBS = 'global'
#config.Data.splitting = 'LumiBased'
#config.Data.unitsPerJob = 10
#config.Data.outLFNDirBase = '/store/user/cfgonzal/'
#config.Data.publication = True
#config.Data.Data.outputDatasetTag = 'SingleMuon_Run2016B_25ns_PromptReco_MiniAODv2'

#config.Site.storageSite = 'T3_US_FNALLPC'


from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = 'SingleMuon_Run2016B_25ns_PromptReco_MiniAODv2_total'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'ConfFile_cfg.py'

config.Data.inputDataset = '/SingleMuon/Run2016B-PromptReco-v2/MINIAOD'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 10
config.Data.outLFNDirBase = '/store/user/cfgonzal/'
config.Data.publication = True
config.Data.outputDatasetTag = 'SingleMuon_Run2016B_25ns_PromptReco_MiniAODv2_total'

config.Site.storageSite = 'T3_US_FNALLPC'
