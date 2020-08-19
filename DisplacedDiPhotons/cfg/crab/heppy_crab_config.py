from CRABClient.UserUtilities import config
from CRABAPI.RawCommand import crabCommand
from CRABClient.ClientExceptions import ClientException
from httplib import HTTPException
import pickle
import optparse
import os
import commands



from CMGTools.DisplacedDiPhotons.samples.loadSamples import *
#selectedComponents = mcSamples
#selectedComponents = dataSamples
selectedComponents = [EGamma_Run2018A_17Sep2018, EGamma_Run2018B_17Sep2018, EGamma_Run2018C_17Sep2018,SingleMuon_Run2018A_17Sep2018,SingleMuon_Run2018B_17Sep2018,SingleMuon_Run2018C_17Sep2018,DoubleMuon_Run2018A_17Sep2018,DoubleMuon_Run2018B_17Sep2018,DoubleMuon_Run2018C_17Sep2018]

parser = optparse.OptionParser()
parser.add_option("-p","--production",dest="prod",default='DDP',help="Name Of Production")
parser.add_option("-u","--username",dest="username",default='bachtis',help="user name")
parser.add_option("-d","--database", dest="database", default='global', help="DBS")
(options,args) = parser.parse_args()


statusVector={}

comps = []

if args[0]=="merge":
    os.system("mkdir merged_"+options.prod)
            
for component in selectedComponents:
    comps.append("merged_"+options.prod+"/"+component.name + "_fakerate.root")

    print("Processing Component",component.name)
    configu = config()
    configu.General.requestName = component.name
    configu.General.workArea =options.prod
    configu.General.transferOutputs = True
    configu.General.transferLogs = True
    configu.JobType.pluginName = 'Analysis'
    configu.JobType.scriptExe = 'heppy_crab_script.sh'
    configu.JobType.psetName = 'heppy_crab_fake_pset.py'
    configu.JobType.maxMemoryMB = 2500
    configu.Data.inputDataset = component.dataset
    configu.Data.inputDBS = options.database
    #configu.Data.inputDBS = 'global'
    #configu.Data.inputDBS = 'phys03'
    configu.JobType.inputFiles = ['heppy_config.py','heppy_crab_script.py','component.pck']
    configu.JobType.outputFiles = ['tree.root','fakerate.root','SkimReport.pck']
    configu.JobType.sendPythonFolder = True
    configu.Data.splitting = 'FileBased'
    if component.isMC:
        configu.Data.unitsPerJob = 3
    if component.isData:
        configu.Data.unitsPerJob = 4
    configu.Data.outLFNDirBase = '/store/user/'+options.username+'/'+options.prod
    configu.Data.publication = False
    configu.Data.ignoreLocality =True
    configu.Site.whitelist = ['T2_US_*']
    configu.Data.outputDatasetTag = component.name
    configu.Site.storageSite = 'T3_US_FNALLPC'
    if args[0]=='list':
        print(component.name)


    if args[0]=='submit':
        if os.path.isdir(options.prod+'/crab_'+component.name):
            print("Directory already exists -doing nothing")
            continue
        f=open("component.pck","w")
        pickle.dump(component,f)
        f.close()

        try:
            crabCommand('submit', config = configu)
        except HTTPException as hte:
            print "Submission  failed: %s" % (hte.headers)
            break;
        except ClientException as cle:
            print "Submission for input dataset %s failed"
            break;
    if args[0]=='status':
        if os.path.isdir(options.prod+'/crab_'+component.name):

            result = crabCommand('status', dir = '/'.join([options.prod,'crab_'+component.name]))
            print('---------------------------------')
            print(component.name+':'+result['status'])
            print('---------------------------------')
            print(result['jobsPerStatus'])
            print('---------------------------------')
            statusVector[component.name] = result['jobsPerStatus']
            if result['status']=='FAILED':
                print("RESUBMITING")
                r = crabCommand('resubmit', dir = '/'.join([options.prod,'crab_'+component.name]))


    if args[0]=='merge':
        if not os.path.isdir(options.prod+'/crab_'+component.name):
            continue
        primaryDataset=component.dataset.split('/')[1]
        directory='/'.join([configu.Data.outLFNDirBase,primaryDataset,component.name])
        #timestamp-get latest
        status,output = commands.getstatusoutput("eos root://cmseos.fnal.gov ls  "+directory)
        directory=directory+'/'+output.split('\n')[-1]
        print("Directory and timestamp used:",directory)
        rootFiles = []
        pckFiles=[]
        #enumeration
        status,output = commands.getstatusoutput("eos root://cmseos.fnal.gov ls  "+directory)
        subdirs=output.split('\n')
        premerged={}
        premerged['tree'] = []
        premerged['fakerate'] = []
        if not os.path.isdir("/uscmst1b_scratch/lpc1/3DayLifetime/"+options.username+"/merged_"+options.prod):
            os.system("mkdir /uscmst1b_scratch/lpc1/3DayLifetime/"+options.username+"/merged_"+options.prod)

        '''
        for d in subdirs: 
            status,output = commands.getstatusoutput("xrdfs root://cmseos.fnal.gov ls -u "+directory+'/'+d)
            rootFiles=list(filter(lambda x: x.find("tree")!=-1,output.split('\n')))
            pckFiles=pckFiles+list(filter(lambda x: x.find(".pck")!=-1,output.split('\n')))
            os.system("hadd merged_"+options.prod+'/'+component.name+'_'+d+'.root '+' '.join(rootFiles))
            premerged['tree'].append("merged_"+options.prod+'/'+component.name+'_'+d+'.root')

            rateFiles=list(filter(lambda x: x.find("fakerate")!=-1,output.split('\n')))
            os.system("hadd merged_"+options.prod+'/'+'fakerate_'+d+'.root '+' '.join(rateFiles))
            premerged['fakerate'].append("merged_"+options.prod+'/'+'fakerate_'+d+'.root')
        '''
        for d in subdirs:
            status,output = commands.getstatusoutput("xrdfs root://cmseos.fnal.gov ls -u "+directory+'/'+d)
            rootFiles=list(filter(lambda x: x.find("tree")!=-1,output.split('\n')))
            pckFiles=pckFiles+list(filter(lambda x: x.find(".pck")!=-1,output.split('\n')))
            print(rootFiles)
            os.system("mkdir tmp")
            for f in rootFiles:
                os.system("xrdcp "+f+" tmp/")
            os.system("hadd /uscmst1b_scratch/lpc1/3DayLifetime/"+options.username+"/merged_"+options.prod+'/'+component.name+'_'+d+'.root  tmp/tree*')
            os.system("rm -rf tmp")
                #            os.system("hadd /uscmst1b_scratch/lpc1/3DayLifetime/"+options.username+"/merged_"+options.prod+'/'+component.name+'_'+d+'.root '+' '.join(rootFiles))                                                                                                               
            premerged['tree'].append('/uscmst1b_scratch/lpc1/3DayLifetime/'+options.username+"/merged_"+options.prod+'/'+component.name+'_'+d+'.root')
            
            os.system("mkdir tmp")

            rateFiles=list(filter(lambda x: x.find("fakerate")!=-1,output.split('\n')))
            print(rateFiles)
            for f in rateFiles:
                os.system("xrdcp "+f+" tmp/")
            os.system("hadd /uscmst1b_scratch/lpc1/3DayLifetime/"+options.username+"/merged_"+options.prod+'/'+'fakerate_'+d+'.root tmp/fakerate*')
            premerged['fakerate'].append("/uscmst1b_scratch/lpc1/3DayLifetime/"+options.username+"/merged_"+options.prod+'/'+'fakerate_'+d+'.root')
            os.system("rm -rf tmp")


        print "PRINTOUT"
        print(pckFiles)
        print(premerged['tree'])
        print(premerged['fakerate'])

        os.system("hadd merged_"+options.prod+'/'+component.name+'.root '+' '.join(premerged['tree']))
        os.system("rm "+' '.join(premerged['tree']))
        os.system("hadd merged_"+options.prod+'/'+component.name+'_fakerate.root '+' '.join(premerged['fakerate']))
        os.system("rm "+' '.join(premerged['fakerate']))

        os.system("mkdir tmp")
        for f in pckFiles:
            os.system("xrdcp "+f+" tmp/")
            
        status,output = commands.getstatusoutput("ls tmp/*.pck")
        pckFiles=output.split('\n')
        sumWeights=0
        for p in pckFiles:
            f=open(p)
            counter=pickle.load(f)
            if len(counter)>1 and counter[1][0]=='Sum Weights':
                sumWeights=sumWeights+counter[1][1]
            else:    
                sumWeights=sumWeights+counter[0][1]
            f.close()
        os.system("rm -rf tmp")    
        output=dict()
        if component.isMC:
            output['weight'] = 1.0/sumWeights
            output['events'] = sumWeights
        else:
            output['weight'] = 1.0
            output['events'] = sumWeights

        output['sigma']  = 1.0

        f=open("merged_"+options.prod+'/'+component.name+'.pck',"w")
        pickle.dump(output,f)
        f.close()

if args[0] == 'merge':
    os.system("hadd merged_"+options.prod+'/fakerate.root '+' '.join(comps))

if args[0]=='status':
    print statusVector
    print('---------------------------------')
    status = {}

    nJobs = 0

    for comp in statusVector:
        for s in statusVector[comp]:
            if s in status:
                status[s] += statusVector[comp][s]
            else:
                status[s] = statusVector[comp][s]
            nJobs += statusVector[comp][s]

    print "Status for all jobs: {}".format(status)
    print "Total jobs = {}".format(nJobs)

