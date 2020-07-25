import PhysicsTools.HeppyCore.framework.config as cfg
from PhysicsTools.Heppy.analyzers.core.all import *
from PhysicsTools.Heppy.analyzers.objects.all import *
from PhysicsTools.Heppy.analyzers.gen.all import *
import os

from CMGTools.VVResonances.analyzers.Skimmer  import *
from CMGTools.DisplacedDiPhotons.analyzers.VHGGBuilder  import *
from CMGTools.RootTools.samples.triggers_13TeV_DATA2018 import *

# Pick individual events (normally not in the path)
eventSelector = cfg.Analyzer(
    EventSelector,name="EventSelector",
    toSelect = []  # here put the event numbers (actual event numbers from CMSSW)
    )

skimAnalyzer = cfg.Analyzer(
    SkimAnalyzerCount, name='skimAnalyzerCount',
    useLumiBlocks = False,
    )

# Apply json file (if the dataset has one)
jsonAna = cfg.Analyzer(
    JSONAnalyzer, name="JSONAnalyzer",
)


# Filter using the 'triggers' and 'vetoTriggers' specified in the dataset
triggerAna = cfg.Analyzer(
    TriggerBitFilter, name="TriggerBitFilter",
)

# Create flags for trigger bits
triggerFlagsAna = cfg.Analyzer(
    TriggerBitAnalyzer, name="TriggerFlags",

    processName = 'HLT',
    fallbackProcessName = 'HLT2',
    unrollbits = True,
    triggerBits = {
    "ISOMU":triggers_1mu_iso,
    "MU":triggers_1mu_noniso,
    "ISOMUMU":triggers_mumu_iso,
    "MUMU":triggers_1mu_noniso,
    "ELE":triggers_1e_noniso,
    "ISOELE":triggers_1e_iso,
    "EE": triggers_ee
    }
    )



# Create flags for MET filter bits
eventFlagsAna = cfg.Analyzer(
    TriggerBitAnalyzer, name="EventFlags",
    processName = 'PAT',
    fallbackProcessName = 'RECO',
    outprefix   = 'Flag',
    triggerBits = {
        "HBHENoiseFilter" : [ "Flag_HBHENoiseFilter" ],
        "HBHENoiseIsoFilter" : [ "Flag_HBHENoiseIsoFilter" ],
        "CSCTightHaloFilter" : [ "Flag_CSCTightHaloFilter" ],
        "globalTightHalo2016Filter" : [ "Flag_globalTightHalo2016Filter" ],
        "hcalLaserEventFilter" : [ "Flag_hcalLaserEventFilter" ],
        "EcalDeadCellTriggerPrimitiveFilter" : [ "Flag_EcalDeadCellTriggerPrimitiveFilter" ],
        "goodVertices" : [ "Flag_goodVertices" ],
        "trackingFailureFilter" : [ "Flag_trackingFailureFilter" ],
        "eeBadScFilter" : [ "Flag_eeBadScFilter" ],
        "ecalLaserCorrFilter" : [ "Flag_ecalLaserCorrFilter" ],
        "trkPOGFilters" : [ "Flag_trkPOGFilters" ],
        "trkPOG_manystripclus53X" : [ "Flag_trkPOG_manystripclus53X" ],
        "trkPOG_toomanystripclus53X" : [ "Flag_trkPOG_toomanystripclus53X" ],
        "trkPOG_logErrorTooManyClusters" : [ "Flag_trkPOG_logErrorTooManyClusters" ],
        "METFilters" : [ "Flag_METFilters" ],
    }
    )

from CMGTools.TTHAnalysis.analyzers.badChargedHadronAnalyzer import badChargedHadronAnalyzer
badChargedHadronAna = cfg.Analyzer(
    badChargedHadronAnalyzer, name = 'badChargedHadronAna',
    muons='slimmedMuons',
    packedCandidates = 'packedPFCandidates',
)

from CMGTools.TTHAnalysis.analyzers.badMuonAnalyzer import badMuonAnalyzer
badMuonAna = cfg.Analyzer(
    badMuonAnalyzer, name = 'badMuonAna',
    muons='slimmedMuons',
    packedCandidates = 'packedPFCandidates',
)



# Select a list of good primary vertices (generic)
vertexAna = cfg.Analyzer(
    VertexAnalyzer, name="VertexAnalyzer",
    vertexWeight = None,
    fixedWeight = 1,
    verbose = False
    )



# This analyzer actually does the pile-up reweighting (generic)
pileUpAna = cfg.Analyzer(
    PileUpAnalyzer, name="PileUpAnalyzer",
    true = True,  # use number of true interactions for reweighting
    autoPU = True,
    makeHists=False
)


## Gen Info Analyzer (generic, but should be revised)
genAna = cfg.Analyzer(
    GeneratorAnalyzer, name="GeneratorAnalyzer",
    # BSM particles that can appear with status <= 2 and should be kept
    stableBSMParticleIds = [ 1000022 ],
    # Particles of which we want to save the pre-FSR momentum (a la status 3).
    # Note that for quarks and gluons the post-FSR doesn't make sense,
    # so those should always be in the list
    savePreFSRParticleIds = [ 1,2,3,4,5, 11,12,13,14,15,16, 21 ],
    # Make also the list of all genParticles, for other analyzers to handle
    makeAllGenParticles = True,
    # Make also the splitted lists
    makeSplittedGenLists = True,
    allGenTaus = False,
    # Save LHE weights from LHEEventProduct
    makeLHEweights = True,
    # Print out debug information
    verbose = False,
    )

lepAna = cfg.Analyzer(
    LeptonAnalyzer, name="leptonAnalyzer",
    # input collections
    muons='slimmedMuons',
    electrons='slimmedElectrons',
    rhoMuon= 'fixedGridRhoFastjetAll',
    rhoElectron = 'fixedGridRhoFastjetAll',
    # energy scale corrections and ghost muon suppression (off by default)
    doMuonScaleCorrections=False,
    doElectronScaleCorrections=False, # "embedded" in 5.18 for regression
    doSegmentBasedMuonCleaning=False,
    # inclusive very loose muon selection
    inclusive_muon_id  = "",
    inclusive_muon_pt  = 20.0,
    inclusive_muon_eta = 2.4,
    inclusive_muon_dxy = 0.3,
    inclusive_muon_dz  = 20,
    muon_dxydz_track = "innerTrack",
    # loose muon selection
    loose_muon_id     = "",
    loose_muon_pt     = 20.0,
    loose_muon_eta    = 2.4,
    loose_muon_dxy    = 0.3,
    loose_muon_dz     = 20.0,
    loose_muon_isoCut = lambda x:True,
    # inclusive very loose electron selection
    inclusive_electron_id  = "",
    inclusive_electron_pt  = 10.0,
    inclusive_electron_eta = 2.5,
    inclusive_electron_dxy = 0.2,
    inclusive_electron_dz  = 0.2,
    inclusive_electron_lostHits = 1.0,
    # loose electron selection
    loose_electron_id     = "",
    loose_electron_pt     = 10.0,
    loose_electron_eta    = 2.5,
    loose_electron_dxy    = 0.2,
    loose_electron_dz     = 0.2,
    loose_electron_lostHits = 1.0,
    loose_electron_isoCut = lambda x:True,

    # muon isolation correction method (can be "rhoArea" or "deltaBeta")
    mu_isoCorr = "deltaBeta",
    mu_effectiveAreas = "Spring15_25ns_v1", #(can be 'Data2012' or 'Phys14_25ns_v1')
    # electron isolation correction method (can be "rhoArea" or "deltaBeta")
    ele_isoCorr = "rhoArea" ,
    el_effectiveAreas = "Spring15_25ns_v1" , #(can be 'Data2012' or 'Phys14_25ns_v1')
    ele_tightId = "" ,
    # Mini-isolation, with pT dependent cone: will fill in the miniRelIso, miniRelIsoCharged, miniRelIsoNeutral variables of the leptons (see https://indico.cern.ch/event/368826/ )
    doMiniIsolation = False, # off by default since it requires access to all PFCandidates
    packedCandidates = 'packedPFCandidates',
    miniIsolationPUCorr = 'deltaBeta', # Allowed options: 'rhoArea' (EAs for 03 cone scaled by R^2), 'deltaBeta', 'raw' (uncorrected), 'weights' (delta beta weights; not validated)
    miniIsolationVetoLeptons = 'inclusive', # use 'inclusive' to veto inclusive leptons and their footprint in all isolation cones
    # minimum deltaR between a loose electron and a loose muon (on overlaps, discard the electron)
    min_dr_electron_muon = 0.0,
    # do MC matching
    do_mc_match = True, # note: it will in any case try it only on MC, not on data
    match_inclusiveLeptons = False, # match to all inclusive leptons
    )



metAna = cfg.Analyzer(
    METAnalyzer, name="metAnalyzer",
    metCollection     = "slimmedMETs",
    noPUMetCollection = "slimmedMETs",
    copyMETsByValue = False,
    storePuppiExtra=False,
    doTkMet = False,
    doPuppiMet = False,
    doMetNoPU = True,
    doMetNoMu = False,
    doMetNoEle = False,
    doMetNoPhoton = False,
    recalibrate = False, # or "type1", or True
    applyJetSmearing = False, # does nothing unless the jet smearing is turned on in the jet analyzer
    old74XMiniAODs = False, # set to True to get the correct Raw MET when running on old 74X MiniAODs
    jetAnalyzerPostFix = "",
    candidates='packedPFCandidates',
    candidatesTypes='std::vector<pat::PackedCandidate>',
    dzMax = 0.1,
    collectionPostFix = "",
    )



photonAna = cfg.Analyzer(
    PhotonAnalyzer,
    photons='slimmedPhotons',
    ptMin = 10,
    etaMax = 2.5,
    # energy scale corrections (off by default)
    doPhotonScaleCorrections=False, 
    gammaID = "PhotonCutBasedIDLoose_CSA14",
    rhoPhoton = 'fixedGridRhoFastjetAll',
    gamma_isoCorr = 'rhoArea',
    # Footprint-removed isolation, removing all the footprint of the photon
    doFootprintRemovedIsolation = True, # off by default since it requires access to all PFCandidates
    packedCandidates = 'packedPFCandidates',
    footprintRemovedIsolationPUCorr = 'rhoArea', # Allowed options: 'rhoArea', 'raw' (uncorrected)
    conversionSafe_eleVeto = False,
    do_mc_match = True,
    do_randomCone = False,
)



#vvAna = cfg.Analyzer(
#    VVBuilder,name='vvAna',
#    suffix = '',
#    fDiscriminatorB = "pfDeepCSVJetTags:probb",
#    fDiscriminatorBB = "pfDeepCSVJetTags:probbb",
#    fDiscriminatorC = "pfDeepCSVJetTags:probc",
#    fDiscriminatorL = "pfDeepCSVJetTags:probudsg",
#    btagCSVFile = "${CMSSW_BASE}/src/CMGTools/VVResonances/data/DeepCSV_94XSF_V2_B_F.csv",
#    subjetBtagCSVFile = "${CMSSW_BASE}/src/CMGTools/VVResonances/data/subjet_DeepCSV_94XSF_V4_B_F.csv",
#    puppiJecCorrFile = "${CMSSW_BASE}/src/CMGTools/VVResonances/data/puppiCorr.root"
#
#)


vhGGSkimmer = cfg.Analyzer(
    Skimmer,
    name='vhGGSkimmer',
    required = ['ZX','WX','ZXX','WXX','looseZX','looseWX','looseZXX','looseWXX']
)

vhGGAna = cfg.Analyzer(
    VHGGBuilder,
    name='vhGGBuilder'
)



coreSequence = [
   #eventSelector,
    skimAnalyzer,
    jsonAna,
    triggerAna,
    pileUpAna,
    genAna,
    vertexAna,
    lepAna,
    metAna,
    photonAna,
    eventFlagsAna,
    triggerFlagsAna,
    badMuonAna,
    badChargedHadronAna,
    vhGGAna,
    vhGGSkimmer

]
