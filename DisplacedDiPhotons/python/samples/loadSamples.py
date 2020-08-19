import PhysicsTools.HeppyCore.framework.config as cfg
import os

# Load backgrounds from common place
from CMGTools.RootTools.samples.samples_13TeV_RunIIAutumn18MiniAOD import *
from CMGTools.RootTools.samples.samples_13TeV_DATA2018_MiniAOD import *
from CMGTools.RootTools.samples.triggers_13TeV_DATA2018 import *
from CMGTools.RootTools.samples.ComponentCreator import ComponentCreator

kreator = ComponentCreator()

# Signal samples
WHtoLNuGGdd0 = kreator.makeMyPrivateMCComponent("WHtoLNuGGdd0", "/WPlusH_HToSS_SToddgg_MH125_MS15_CTAU0/bachtis-STEP3_Wplus_S15_ctau0-3ee3afd6b5a1410aea6d0b4d52723d06/USER", "PRIVATE", ".*root", "phys03")
WHtoLNuGGdd3000 = kreator.makeMyPrivateMCComponent("WHtoLNuGGdd3000", "/WPlusH_HToSS_SToddgg_MH125_MS15_CTAU3000/bachtis-STEP3_Wplus_S15_ctau3000-3ee3afd6b5a1410aea6d0b4d52723d06/USER", "PRIVATE", ".*root", "phys03")

signalSamples = [WHtoLNuGGdd0, WHtoLNuGGdd3000]

#mcSamples = [DYJetsToLL_M50_LO, WW,WZ,ZZ,WJetsToLNu_LO,TTLep_pow,TTSemi_pow,WGtoLNuG,WGG] + QCD_EMs + QCD_Mus + DYJetsToLLM50HT + WJetsToLNuHT

#mcSamples = [DYJetsToLL_M50_LO]  #For testing
#mcSamples = [WW,ZZ]
mcSamples = signalSamples

dataSamples_EGamma = [EGamma_Run2018A_17Sep2018, EGamma_Run2018B_17Sep2018, EGamma_Run2018C_17Sep2018, EGamma_Run2018D_PromptReco_v2]

dataSamples_SingleMuon = [SingleMuon_Run2018A_17Sep2018,SingleMuon_Run2018B_17Sep2018,SingleMuon_Run2018C_17Sep2018,SingleMuon_Run2018D_PromptReco_v2]

dataSamples_DoubleMuon = [DoubleMuon_Run2018A_17Sep2018,DoubleMuon_Run2018B_17Sep2018,DoubleMuon_Run2018C_17Sep2018,DoubleMuon_Run2018D_PromptReco_v2]

for s in dataSamples_SingleMuon:
    s.triggers = triggers_1mu_noniso+triggers_1mu_iso
    s.vetoTriggers = []

for s in dataSamples_DoubleMuon:
    s.triggers = triggers_mumu_iso+triggers_mumu_noniso
    s.vetoTriggers = triggers_1mu_noniso+triggers_1mu_iso

for s in dataSamples_EGamma:
    s.triggers = triggers_1e_iso+triggers_1e_noniso+triggers_ee
    s.vetoTriggers = triggers_1mu_noniso+triggers_1mu_iso+triggers_mumu_iso+triggers_mumu_noniso



dataSamples= dataSamples_EGamma + dataSamples_SingleMuon + dataSamples_DoubleMuon




# Define splitting
from CMGTools.TTHAnalysis.setup.Efficiencies import *
dataDir = "$CMSSW_BASE/src/CMGTools/VVResonances/data"
for comp in mcSamples:
    comp.isMC = True
    comp.isData = False
    comp.splitFactor = 300
    comp.puFileMC=dataDir+"/pileup_MC.root"
    comp.puFileData=dataDir+"/pileup_DATA.root"
    comp.efficiency = eff2012
    comp.triggers=[]


for comp in dataSamples:
    comp.splitFactor = 500
    comp.isMC = False
    comp.isData = True

