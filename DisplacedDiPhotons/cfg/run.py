##########################################################
##       GENERIC SUSY TREE PRODUCTION CONFIG.           ##
## no skim is applied in this configuration, as it is   ##
## meant only to check that all common modules run ok   ##
##########################################################


#AAA
###
def autoAAA(selectedComponents,runOnlyRemoteSamples=False,forceAAA=False):
    newComp=[]
    import re
    from CMGTools.Production import changeComponentAccessMode
    from CMGTools.Production.localityChecker import LocalityChecker
    tier2Checker = LocalityChecker("T2_CH_CERN", datasets="/*/*/MINIAOD*")
    for comp in selectedComponents:
        if len(comp.files)==0:
            continue
        if not hasattr(comp,'dataset'): continue
        if not re.match("/[^/]+/[^/]+/MINIAOD(SIM)?", comp.dataset): continue
        if "/store/" not in comp.files[0]: continue
        if re.search("/store/(group|user|cmst3)/", comp.files[0]): continue
        if (not tier2Checker.available(comp.dataset)) or forceAAA:
            print "Dataset %s is not available, will use AAA" % comp.dataset
            changeComponentAccessMode.convertComponent(comp, "root://cms-xrd-global.cern.ch/%s")
            if 'X509_USER_PROXY' not in os.environ or "/afs/" not in os.environ['X509_USER_PROXY']:
                raise RuntimeError, "X509_USER_PROXY not defined or not pointing to /afs"
            newComp.append(comp)
    if runOnlyRemoteSamples:
        return newComp
    else:
        return selectedComponents


def runOnFNAL(selectedComponents,JSON):
    newComp=[]
    import re
    from CMGTools.Production import changeComponentAccessMode
    from CMGTools.Production.localityChecker import LocalityChecker

    for comp in selectedComponents:
        print comp.name
        if comp.isData:
            comp.json=JSON 
        if len(comp.files)==0:
            continue
        #if it is data attach the local JSON file
        if not hasattr(comp,'dataset'): continue
        if not re.match("/[^/]+/[^/]+/MINIAOD(SIM)?", comp.dataset): continue
        if "/store/" not in comp.files[0]: continue
        if re.search("/store/(group|user|cmst3)/", comp.files[0]): continue
        changeComponentAccessMode.convertComponent(comp, "root://cms-xrd-global.cern.ch/%s")
        if 'X509_USER_PROXY' not in os.environ:
            raise RuntimeError, "X509_USER_PROXY not defined or not pointing to /afs"

    return selectedComponents


def autoConfig(selectedComponents,sequence,services=[],xrd_aggressive=2):
    import PhysicsTools.HeppyCore.framework.config as cfg
    from PhysicsTools.HeppyCore.framework.eventsfwlite import Events
    from CMGTools.TTHAnalysis.tools.EOSEventsWithDownload import EOSEventsWithDownload
    event_class = EOSEventsWithDownload
    EOSEventsWithDownload.aggressive = xrd_aggressive
    if getHeppyOption("nofetch") or getHeppyOption("isCrab"):
        event_class = Events
    return cfg.Config( components = selectedComponents,
                     sequence = sequence,
                     services = services,
                     events_class = event_class)


###



import CMGTools.RootTools.fwlite.Config as cfg
from CMGTools.RootTools.fwlite.Config import printComps
from CMGTools.RootTools.RootTools import *
from PhysicsTools.HeppyCore.framework.heppy_loop import getHeppyOption

#Load all common analyzers
from CMGTools.DisplacedDiPhotons.analyzers.core_cff import *







#import pdb;pdb.set_trace()

#-------- Analyzer
from CMGTools.DisplacedDiPhotons.analyzers.tree_cff import *

#-------- SEQUENCE

sequence = cfg.Sequence(coreSequence+[vhGGTreeProducer])



#-------- HOW TO RUN
#runCommand values
#-1 -> CRAB
# 0 -> Condor LPC
# 1 ->test

runCommand = 1

if runCommand<0:
    import pickle
    f=open("component.pck")
    component=pickle.load(f)
    f.close()
    selectedComponents = [component]
    selectedComponents[0].name="Output"
else:
    from CMGTools.DisplacedDiPhotons.samples.loadSamples import *
    ####
    #selectedComponents = mcSamples
    selectedComponents = dataSamples
    ####
    if runCommand==1:
        selectedComponents = [selectedComponents[0]]
        for c in selectedComponents:
            c.files = c.files[:1]
            c.splitFactor = 1
    else:    # full scale production
        # split samples in a smarter way
        from CMGTools.RootTools.samples.configTools import configureSplittingFromTime, printSummary
        configureSplittingFromTime(selectedComponents, 60, 1)  # means 70 ms per event, job to last 12h
        # prnt summary of components to process
        printSummary(selectedComponents)

selectedComponents=runOnFNAL(selectedComponents,"$CMSSW_BASE/src/CMGTools/DisplacedDiPhotons/data/JSON2018.txt")
config=autoConfig(selectedComponents,sequence)
