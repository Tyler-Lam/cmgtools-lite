import ROOT

from CMGTools.VVResonances.plotting.TreePlotter import TreePlotter
from CMGTools.VVResonances.plotting.MergedPlotter import MergedPlotter
from CMGTools.VVResonances.plotting.StackPlotter import StackPlotter
from CMGTools.VVResonances.plotting.FakePlotter import FakePlotter
ROOT.gROOT.SetBatch(True)

def normalizeStack(stack):
    scale = 0
    newStack = ROOT.THStack("stack", "")
    for hist in stack.GetHists():
        scale += hist.Integral()
    if scale > 0:
        for hist in stack.GetHists():
            hist.Scale(1.0/scale)
            newStack.Add(hist)
    return newStack


cuts = {}
cuts['default'] = {} #Default Signal Region cuts

# Cuts on isolation for the leptons/photons
cuts['relIso'] = {}
cuts['relIso']['W'] = "WX_W_l1_relIso03<.1&&WX_X_g1_relIso<.1&&WX_X_g2_relIso<.1"
cuts['relIso']['Z'] = "ZX_Z_l1_relIso03<.1&&ZX_Z_l2_relIso03<.1&&ZX_X_g1_relIso<.1&&ZX_X_g2_relIso<.1"

# Cuts on MVA
cuts['mva'] = {}
for v in ['W', 'Z']:
    temp = {}
    for n in [1,2]:
        tag = "{}X_X_g{}".format(v,n)
        temp[n] =  "({t}_isEE&&{t}_MVANonTrigV1Values>-0.26)||({t}_isEB&&{t}_MVANonTrigV1Values>-0.02)".format(t=tag)
    cuts['mva'][v] = "("+temp[1]+")*("+temp[2]+")"

# Cuts for FSR on Z->mu+mu and misID for W->e+nu
cuts['fsr'] = {}
cuts['fsr']['Z'] = {}
cuts['fsr']['Z']['MU'] = "ZX_hasFSR==0"
cuts['fsr']['Z']['ELE'] = "1"
cuts['fsr']['W'] = {}
cuts['fsr']['W']['MU'] = "1"
cuts['fsr']['W']['ELE'] = "WX_misID!=1&&WX_misID!=3"

# Cuts for HLT
cuts['HLT'] = {}
cuts['HLT']['MU'] = "(HLT_ISOMU&&!HLT_ISOELE)"
cuts['HLT']['ELE'] = "(HLT_ISOELE&&!HLT_ISOMU)"

#Defining Signal region default cuts: iso, mva, fsr/misID
for v in ['W', 'Z']:
    cuts['default'][v] = {}
    for l in ['ELE', 'MU']:
        cuts['default'][v][l] = "("+cuts['relIso'][v]+")*("+cuts['mva'][v]+")*("+cuts['fsr'][v][l]+")"
cuts['default']['Z']['ELE'] =cuts['default']['Z']['ELE']+"*(ZX_Z_l1_pt>35&&ZX_Z_l2_pt>35)"


#Cuts for MC based matching
cuts['pi0'] = {}
cuts['isMCgamma'] = {}
cuts['mcFSR'] = {}
cuts['mcISR'] = {}
cuts['mcReal'] = {}
for v in ['Z', 'W']:
    cuts['pi0'][v] = {}
    cuts['isMCgamma'][v] = "abs({}X_X_g1_mcMatchId)==22&&abs({}X_X_g2_mcMatchId)==22".format(v,v)
    cuts['mcISR'][v] = {}
    cuts['mcFSR'][v] = {}
    cuts['mcReal'][v] = {}
    for g in ['g1', 'g2']:
        cuts['pi0'][v][g] = "abs({}X_X_{}_mcMotherId)==111".format(v,g)    
        cuts['mcISR'][v][g] = "(abs({v}X_X_{g}_mcMatchId)==22&&abs({v}X_X_{g}_mcMotherId)<9)".format(v=v, g=g)
        cuts['mcFSR'][v][g] = "(abs({v}X_X_{g}_mcMatchId)==22&&(abs({v}X_X_{g}_mcMotherId)==11||abs({v}X_X_{g}_mcMotherId)==13))".format(v=v, g=g)
        cuts['mcReal'][v][g] = "("+cuts['mcISR'][v][g]+"||"+cuts['mcFSR'][v][g]+")"



date = "08_10_20"
# Load MC Samples
dyjSamples = ['DYJetsToLL_M50_LO']
wjSamples = ['WJetsToLNu_LO']
wgSamples = ['WGG', 'WGToLNuG']
vvSamples = ['WW', 
             'WZ',
             'ZZ']
qcdSamples = ['QCD_Mu15',
              'QCD_Pt15to20_EMEnriched',
              'QCD_Pt15to20_Mu5',
              'QCD_Pt20to30_EMEnriched',
              'QCD_Pt20to30_Mu5',
              'QCD_Pt30to50_EMEnriched',
              'QCD_Pt30to50_Mu5',
              'QCD_Pt50to80_EMEnriched',
              'QCD_Pt50to80_Mu5',
              'QCD_Pt80to120_EMEnriched',
              'QCD_Pt80to120_Mu5',
              'QCD_Pt120to170_EMEnriched',
              'QCD_Pt120to170_Mu5',
              'QCD_Pt170to300_EMEnriched',
              'QCD_Pt170to300_Mu5',
              'QCD_Pt300toInf_EMEnriched',
              'QCD_Pt300to470_Mu5',
              'QCD_Pt470to600_Mu5',
              'QCD_Pt600to800_Mu5',
              'QCD_Pt800to1000_Mu5',
              'QCD_Pt1000toInf_Mu5']
ttSamples = ['TTLep_pow',
             'TTSemi_pow']

mcSamples = dyjSamples + wjSamples + wgSamples + vvSamples + qcdSamples + ttSamples
twoFakesPlotter = {}
oneFakePlotter = {}
zeroFakesPlotter = {}


dyjPlotters = []
wjPlotters = []
wgPlotters = []
vvPlotters = []
qcdPlotters = []
ttPlotters = []

# Load Data Samples

egammaSamples = ['EGamma_Run2018A_17Sep2018',
                 'EGamma_Run2018B_17Sep2018',
                 'EGamma_Run2018C_17Sep2018',
                 'EGamma_Run2018D_PromptReco_v2']
singleMuSamples = ['SingleMuon_Run2018A_17Sep2018',
                   'SingleMuon_Run2018B_17Sep2018',
                   'SingleMuon_Run2018C_17Sep2018',
                   'SingleMuon_Run2018D_PromptReco_v2']
doubleMuSamples=['DoubleMuon_Run2018A_17Sep2018',
                 'DoubleMuon_Run2018B_17Sep2018',
                 'DoubleMuon_Run2018C_17Sep2018',
                 'DoubleMuon_Run2018D_PromptReco_v2']
dataSamples = egammaSamples + singleMuSamples + doubleMuSamples

egammaPlotters = []
singleMuPlotters = []
doubleMuPlotters = []
fakePlotters = []

####################################################################################################

#MC
for sample in dyjSamples:
    dyjPlotters.append(TreePlotter('/scratch2/Tyler/samples/DDP/MC_{}/{}.root'.format(date,sample), 'tree'))
    dyjPlotters[-1].setupFromFile('/scratch2/Tyler/samples/DDP/MC_{}/{}.pck'.format(date,sample))
    dyjPlotters[-1].addCorrectionFactor('xsec', 'tree')
    dyjPlotters[-1].addCorrectionFactor('genWeight','tree')
    dyjPlotters[-1].addCorrectionFactor('puWeight','tree')

DYJets = MergedPlotter(dyjPlotters)


for sample in wjSamples:
    wjPlotters.append(TreePlotter('/scratch2/Tyler/samples/DDP/MC_{}/{}.root'.format(date,sample), 'tree'))
    wjPlotters[-1].setupFromFile('/scratch2/Tyler/samples/DDP/MC_{}/{}.pck'.format(date,sample))
    wjPlotters[-1].addCorrectionFactor('xsec', 'tree')
    wjPlotters[-1].addCorrectionFactor('genWeight','tree')
    wjPlotters[-1].addCorrectionFactor('puWeight','tree')

WJets = MergedPlotter(wjPlotters)

for sample in wgSamples:
    wgPlotters.append(TreePlotter('/scratch2/Tyler/samples/DDP/MC_{}/{}.root'.format(date,sample), 'tree'))
    wgPlotters[-1].setupFromFile('/scratch2/Tyler/samples/DDP/MC_{}/{}.pck'.format(date,sample))
    wgPlotters[-1].addCorrectionFactor('xsec', 'tree')
    wgPlotters[-1].addCorrectionFactor('genWeight','tree')
    wgPlotters[-1].addCorrectionFactor('puWeight','tree')

WGs = MergedPlotter(wgPlotters)



for sample in vvSamples:
    vvPlotters.append(TreePlotter('/scratch2/Tyler/samples/DDP/MC_{}/{}.root'.format(date,sample), 'tree'))
    vvPlotters[-1].setupFromFile('/scratch2/Tyler/samples/DDP/MC_{}/{}.pck'.format(date,sample))
    vvPlotters[-1].addCorrectionFactor('xsec', 'tree')
    vvPlotters[-1].addCorrectionFactor('genWeight','tree')
    vvPlotters[-1].addCorrectionFactor('puWeight','tree')

VVs = MergedPlotter(vvPlotters)


for sample in qcdSamples:

    qcdPlotters.append(TreePlotter('/scratch2/Tyler/samples/DDP/MC_{}/{}.root'.format(date,sample), 'tree'))
    qcdPlotters[-1].setupFromFile('/scratch2/Tyler/samples/DDP/MC_{}/{}.pck'.format(date,sample))
    qcdPlotters[-1].addCorrectionFactor('xsec', 'tree')
    qcdPlotters[-1].addCorrectionFactor('genWeight','tree')
    qcdPlotters[-1].addCorrectionFactor('puWeight','tree')

QCDs = MergedPlotter(qcdPlotters)


for sample in ttSamples:
    ttPlotters.append(TreePlotter('/scratch2/Tyler/samples/DDP/MC_{}/{}.root'.format(date,sample), 'tree'))
    ttPlotters[-1].setupFromFile('/scratch2/Tyler/samples/DDP/MC_{}/{}.pck'.format(date,sample))
    ttPlotters[-1].addCorrectionFactor('xsec', 'tree')
    ttPlotters[-1].addCorrectionFactor('genWeight','tree')
    ttPlotters[-1].addCorrectionFactor('puWeight','tree')

TTs = MergedPlotter(ttPlotters)


twoFakes = {}
oneFake = {}
zeroFakes = {}
for v in ['W', 'Z']:
    twoFakesPlotter[v] = []
    oneFakePlotter[v] = []
    zeroFakesPlotter[v] = []
    for sample in mcSamples:

        twoFakesPlotter[v].append(TreePlotter('/scratch2/Tyler/samples/DDP/MC_{}/{}.root'.format(date, sample), 'tree'))
        twoFakesPlotter[v][-1].setupFromFile('/scratch2/Tyler/samples/DDP/MC_{}/{}.pck'.format(date, sample))
        twoFakesPlotter[v][-1].addCorrectionFactor('xsec', 'tree')
        twoFakesPlotter[v][-1].addCorrectionFactor('genWeight', 'tree')
        twoFakesPlotter[v][-1].addCorrectionFactor('puWeight', 'tree')
        twoFakesPlotter[v][-1].addCorrectionFactor("(!("+cuts['mcReal'][v]['g1']+"||"+cuts['mcReal'][v]['g2']+"))", '')

        zeroFakesPlotter[v].append(TreePlotter('/scratch2/Tyler/samples/DDP/MC_{}/{}.root'.format(date, sample), 'tree'))
        zeroFakesPlotter[v][-1].setupFromFile('/scratch2/Tyler/samples/DDP/MC_{}/{}.pck'.format(date, sample))
        zeroFakesPlotter[v][-1].addCorrectionFactor('xsec', 'tree')
        zeroFakesPlotter[v][-1].addCorrectionFactor('genWeight', 'tree')
        zeroFakesPlotter[v][-1].addCorrectionFactor('puWeight', 'tree')
        zeroFakesPlotter[v][-1].addCorrectionFactor("("+cuts['mcReal'][v]['g1']+"&&"+cuts['mcReal'][v]['g2']+")", '')

        oneFakePlotter[v].append(TreePlotter('/scratch2/Tyler/samples/DDP/MC_{}/{}.root'.format(date, sample), 'tree'))
        oneFakePlotter[v][-1].setupFromFile('/scratch2/Tyler/samples/DDP/MC_{}/{}.pck'.format(date, sample))
        oneFakePlotter[v][-1].addCorrectionFactor('xsec', 'tree')
        oneFakePlotter[v][-1].addCorrectionFactor('genWeight', 'tree')
        oneFakePlotter[v][-1].addCorrectionFactor('puWeight', 'tree')
        oneFakePlotter[v][-1].addCorrectionFactor("(("+cuts['mcReal'][v]['g1']+"&&!"+cuts['mcReal'][v]['g2']+")||("+cuts['mcReal'][v]['g2']+"&&!"+cuts['mcReal'][v]['g1']+"))", '')

    twoFakes[v] = MergedPlotter(twoFakesPlotter[v])
    oneFake[v] = MergedPlotter(oneFakePlotter[v])
    zeroFakes[v] = MergedPlotter(zeroFakesPlotter[v])

#DATA
for sample in egammaSamples:
    egammaPlotters.append(TreePlotter('/scratch2/Tyler/samples/DDP/DATA_{}/{}.root'.format(date,sample),'tree'))

for sample in singleMuSamples:
    singleMuPlotters.append(TreePlotter('/scratch2/Tyler/samples/DDP/DATA_{}/{}.root'.format(date,sample),'tree'))
    
for sample in doubleMuSamples:
    doubleMuPlotters.append(TreePlotter('/scratch2/Tyler/samples/DDP/DATA_{}/{}.root'.format(date,sample),'tree'))

for sample in dataSamples:
    fakePlotters.append(FakePlotter('/scratch2/Tyler/samples/DDP/DATA_{}/{}.root'.format(date, sample), 'tree', '/scratch2/Tyler/samples/DDP/DATA_{}/{}_fakerate.root'.format(date, sample)))

dataEMU = MergedPlotter(egammaPlotters + singleMuPlotters + doubleMuPlotters)
fakesEMU = MergedPlotter(fakePlotters)

#SIGNAL
WHtoLNuGG = TreePlotter("/scratch2/Tyler/samples/DDP/SIGNAL_{}/WHtoLNuGGdd0.root".format(date), 'tree')
WHtoLNuGG.setupFromFile("/scratch2/Tyler/samples/DDP/SIGNAL_{}/WHtoLNuGGdd0.pck".format(date))
#WH Xsec * W->lep br * H->SS BR 
WHtoLNuGG.addCorrectionFactor('1.504*0.3*0.1', 'tree')
#WHtoLNuGG.addCorrectionFactor('1.504*0.3', 'tree')
WHtoLNuGG.setFillProperties(0, ROOT.kWhite)
WHtoLNuGG.setLineProperties(1, ROOT.kRed, 3)


####################################################################################################
DYJets.setFillProperties(1001, ROOT.kAzure+5)
WJets.setFillProperties(1001, ROOT.kAzure-9)
VVs.setFillProperties(1001, ROOT.kGreen-5)
WGs.setFillProperties(1001, ROOT.kMagenta-4)
QCDs.setFillProperties(1001, ROOT.kGray)
TTs.setFillProperties(1001, ROOT.kRed-6)

bkgStack = StackPlotter()
bkgStack.addPlotter(DYJets, "DYJets", "DY+Jets", "background")
bkgStack.addPlotter(WJets, "WJets", "W+Jets", "background")
bkgStack.addPlotter(VVs, "VV", "VV", "background")
bkgStack.addPlotter(WGs, "WGs", "WG", "background")
bkgStack.addPlotter(QCDs, "QCD", "QCD", "background")
bkgStack.addPlotter(TTs, "tt", "t#bar{t}", "background")
#bkgStack.addPlotter(WHtoLNuGG, "WHToLNuGG", "Signal", "signal")

XStack = StackPlotter()
XStack.addPlotter(DYJets, "DYJets", "DY+Jets", "background")
XStack.addPlotter(WJets, "WJets", "W+Jets", "background")
XStack.addPlotter(VVs, "VV", "VV", "background")
XStack.addPlotter(WGs, "WGs", "WG", "background")
XStack.addPlotter(QCDs, "QCD", "QCD", "background")
XStack.addPlotter(TTs, "tt", "t#bar{t}", "background")
XStack.addPlotter(fakesEMU,"data_obs", "Data", "data")

fakeBkgStack = {}
for v in ['W', 'Z']:
    twoFakes[v].setFillProperties(1001, ROOT.kAzure+5)
    oneFake[v].setFillProperties(1001, ROOT.kAzure-9)
    zeroFakes[v].setFillProperties(1001, ROOT.kGreen-5)
    fakeBkgStack[v] = StackPlotter()
    fakeBkgStack[v].addPlotter(twoFakes[v], "twoFakes", "V + 2 Fakes", "background")
    fakeBkgStack[v].addPlotter(oneFake[v], "oneFake", "V + 1 Fake + 1 Real", "background")
    fakeBkgStack[v].addPlotter(zeroFakes[v], "noFakes", "V + 2 Real", "background")
    fakeBgkStack[v].addPlotter(fakesEMU, "fakesEMU", "V + 2 Fakes (fakerate)", "data")
###############################################################################################

f = ROOT.TFile("fakerate_test.root", "RECREATE")

for l in ['MU', 'ELE']:
    zStack = fakeBkgStack['Z'].drawStack("ZX_Z_mass", cuts['HLT'][l]+"*"+cuts['default']['Z'][l], "59740", 45, 50, 140)
    zStack['canvas'].SaveAs("plots_08_10_20/ZX_Z_mass_{}_fakeBkg.pdf".format(l))
    zStack['canvas'].Write("ZX_Z_mass_{}_fakeBkg".format(l))

    wStack = fakeBkgStack['W'].drawStack("WX_W_mt", cuts['HLT'][l]+"*"+cuts['default']['W'][l], "59740", 20, 0, 200)
    wStack['canvas'].SaveAs("plots_08_10_20/WX_W_mt_{}_fakerate.pdf".format(l))
    wStack['canvas'].Write("WX_W_mt_{}_fakerate".format(l))

    zStack = bkgStack.drawStack("ZX_Z_mass", cuts['HLT'][l]+"*"+cuts['default']['Z'][l], "59740", 45, 50, 140)
    zStack['canvas'].SaveAs("plots_08_10_20/ZX_Z_mass_{}.pdf".format(l))
    zStack['canvas'].Write("ZX_Z_mass_{}".format(l))

    wStack = bkgStack.drawStack("WX_W_mt", cuts['HLT'][l]+"*"+cuts['default']['W'][l], "59740", 20, 0, 200)
    wStack['canvas'].SaveAs("plots_08_10_20/WX_W_mt_{}.pdf".format(l))
    wStack['canvas'].Write("WX_W_mt_{}".format(l))
