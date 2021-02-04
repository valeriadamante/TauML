import FWCore.ParameterSet.Config as cms
def update(process):
    #print process.dumpPython()
    process.HLT_DoubleMediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_v5 = cms.Path(process.HLTBeginSequence + process.hltL1sDoubleTauBigOR + process.hltPreDoubleMediumChargedIsoPFTauHPS35Trk1eta2p1Reg + process.HLTL2TauJetsL1TauSeededSequence + process.hltDoubleL2Tau26eta2p2 + process.HLTL2p5IsoTauL1TauSeededSequence + process.hltDoubleL2IsoTau26eta2p2)
    process.schedule = cms.Schedule(*[ process.HLTriggerFirstPath, process.HLT_DoubleMediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_v5, process.HLTriggerFinalPath, process.endjob_step ], tasks=[process.patAlgosToolsTask])
    '''
    #my_sequence=cms.Task(*(process.HLTriggerFirstPath, process.HLT_DoubleMediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_v5, process.HLTriggerFinalPath ))
    process.schedule.append(process.HLTriggerFirstPath)
    process.schedule.append(process.HLT_DoubleMediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_v5)
    process.schedule.append(process.HLTriggerFinalPath)
    '''
    process.options.wantSummary = cms.untracked.bool(True)
    return process
