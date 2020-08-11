#!/usr/bin/env python
from __future__ import print_function
import os
# probably easier to fetch everything without subdirs, but that's up to user preferences
import PhysicsTools.HeppyCore.framework.config as cfg
cfg.Analyzer.nosubdir=True

import PSet
import sys
import re
print("ARGV:",sys.argv)
JobNumber=sys.argv[1]
crabFiles=PSet.process.source.fileNames
print(crabFiles)
firstInput = crabFiles[0]
print("--------------- using edmFileUtil to convert PFN to LFN -------------------------")
for i in xrange(0,len(crabFiles)) :
     pfn=os.popen("edmFileUtil -d %s"%(crabFiles[i])).read() 
     pfn=re.sub("\n","",pfn)
     print(crabFiles[i],"->",pfn)
     #crabFiles[i]=pfn
     crabFiles[i]="root://cms-xrd-global.cern.ch/"+crabFiles[i]

import imp
handle = open("heppy_config.py", 'r')
cfo = imp.load_source("heppy_config", "heppy_config.py", handle)
config = cfo.config
handle.close()


config.components[0].files=crabFiles

#Use a simple self configured looper so that we know where the output goes
from PhysicsTools.HeppyCore.framework.looper import Looper
looper = Looper( 'Output', config, nPrint = 1)
looper.loop()
looper.write()


#create bare minimum FJR
fwkreport='''
<FrameworkJobReport>

</FrameworkJobReport>
'''

f1=open('./FrameworkJobReport.xml', 'w+')
f1.write(fwkreport)

#place the file in the main folder
os.rename("Output/tree.root", "tree.root")
os.rename("Output/SkimReport.pck", "SkimReport.pck")
os.rename("Output/fakerate.root", "fakerate.root")
