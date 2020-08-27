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


def convertToPoisson(h):
    graph = ROOT.TGraphAsymmErrors()
    q = (1-0.6827)/2.

    for i in range(1,h.GetNbinsX()+1):
        x=h.GetXaxis().GetBinCenter(i)
        xLow =h.GetXaxis().GetBinLowEdge(i) 
        xHigh =h.GetXaxis().GetBinUpEdge(i) 
        y=h.GetBinContent(i)
        yLow=0
        yHigh=0
        if y !=0.0:
            yLow = y-ROOT.Math.chisquared_quantile_c(1-q,2*y)/2.
            yHigh = ROOT.Math.chisquared_quantile_c(q,2*(y+1))/2.-y
            graph.SetPoint(i-1,x,y)
            graph.SetPointEYlow(i-1,yLow)
            graph.SetPointEYhigh(i-1,yHigh)
            graph.SetPointEXlow(i-1,0.0)
            graph.SetPointEXhigh(i-1,0.0)


    graph.SetMarkerStyle(20)
    graph.SetLineWidth(2)
    graph.SetMarkerSize(1.)
    graph.SetMarkerColor(ROOT.kBlack)
    

    return graph    


cuts = {}
# Cuts on isolation for the leptons/photons
cuts['relIso'] = {}
cuts['relIso']['W'] = "WX_W_l1_relIso03<.1&&WX_X_g1_relIso<.1&&WX_X_g2_relIso<.1"
cuts['relIso']['Z'] = "ZX_Z_l1_relIso03<.1&&ZX_Z_l2_relIso03<.1&&ZX_X_g1_relIso<.1&&ZX_X_g2_relIso<.1"
# Cuts on MVA
cuts['mva'] = {}
cuts['default'] = {}
for v in ['W', 'Z']:
    temp = {}
    for n in [1,2]:
        tag = "{}X_X_g{}".format(v,n)
        temp[n] =  "({t}_isEE&&{t}_MVANonTrigV1Values>-0.26)||({t}_isEB&&{t}_MVANonTrigV1Values>-0.02)".format(t=tag)
    cuts['mva'][v] = temp[1]+"&&"+temp[2]

# Cuts for FSR on Z->mumu and misID for W->e+nu
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

#Defining Signal region default cuts
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
    for g in ['g1', 'g2']:
        cuts['pi0'][v][g] = "abs({}X_X_{}_mcMotherId)==111".format(v,g)    
        cuts['mcISR'][v][g] = "(abs({v}X_X_{g}_mcMatchId)==22&&abs({v}X_X_{g}_mcMotherId)<9)".format(v=v, g=g)
        cuts['mcFSR'][v][g] = "(abs({v}X_X_{g}_mcMatchId)==22&&(abs({v}X_X_{g}_mcMotherId)==11||abs({v}X_X_{g}_mcMotherId)==13))".format(v=v, g=g)
    cuts['mcReal'][v] = "(("+cuts['mcISR'][v]['g1']+"||"+cuts['mcFSR'][v]['g1']+")&&("+cuts['mcISR'][v]['g2']+"||"+cuts['mcFSR'][v]['g2']+"))"


date = "08_10_20"
# Load MC Samples
'''
dyjSamples = ["DYJetsToLL_M50_HT100to200",
              "DYJetsToLL_M50_HT200to400",
              "DYJetsToLL_M50_HT400to600",
              "DYJetsToLL_M50_HT600to800",
              "DYJetsToLL_M50_HT800to1200",
              "DYJetsToLL_M50_HT2500toInf"]
wjSamples = ['WJetsToLNu_HT100to200',
             'WJetsToLNu_HT200to400',
             'WJetsToLNu_HT400to600',
             'WJetsToLNu_HT600to800',
             'WJetsToLNu_HT800to1200',
             'WJetsToLNu_HT1200to2500',
             'WJetsToLNu_HT2500toInf']
'''
dyjSamples = ['DYJetsToLL_M50_LO']
wjSamples = ['WJetsToLNu_LO']
wgSamples = ['WGG', 'WGToLNuG']
#wgSamples = []
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
#qcdSamples = ['QCD_Pt80to120',
#              'QCD_Pt120to170',
#              'QCD_Pt170to300',
#              'QCD_Pt300to470',
#              'QCD_Pt470to600',
#              'QCD_Pt470to600_ext1',
#              'QCD_Pt600to800',
#              'QCD_Pt800to1000',
#              'QCD_Pt1000to1400',
#              'QCD_Pt1400to1800',
#              'QCD_Pt1400to1800_ext1',
#              'QCD_Pt1800to2400',
#              'QCD_Pt1800to2400_ext1',
#              'QCD_Pt2400to3200',
#              'QCD_Pt2400to3200_ext1',
#              'QCD_Pt3200toInf']


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
#WHtoLNuGG.addCorrectionFactor('1.504*0.3*0.1', 'tree')
WHtoLNuGG.addCorrectionFactor('1.504*0.3', 'tree')
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
#bkgStack.addPlotter(QCDs, "QCD", "QCD", "background")
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
###############################################################################################

f = ROOT.TFile("fakerate_test.root", "RECREATE")

for l in ['MU', 'ELE']:
    zStack = XStack.drawStack("ZX_Z_mass", cuts['HLT'][l]+"*"+cuts['default']['Z'][l], "59740", 45, 50, 140)
    zStack['canvas'].SaveAs("ZX_Z_mass_{}_fakerate.pdf".format(l))
    zStack['canvas'].Write("ZX_Z_mass_{}_fakerate".format(l))

    wStack = XStack.drawStack("WX_W_mt", cuts['HLT'][l]+"*"+cuts['default']['W'][l], "59740", 20, 0, 200)
    wStack['canvas'].SaveAs("WX_W_mt_{}_fakerate.pdf".format(l))
    wStack['canvas'].Write("WX_W_mt_{}_fakerate".format(l))


#To test new MC: Plotting Z mass and W mt with default cuts
'''
f = ROOT.TFile("plots_dxy.root", "RECREATE")
masses = [10,20,30,40,50,60]
for v in ['W', 'Z']:
    for l in ['ELE', 'MU']:
        for m in masses:
            hData = WHtoLNuGG.drawTH1("{}X_X_vertex{}_d0".format(v,m), "HLT_ISO{}&&n{}XX<1&&{}X_X_vertex{}_valid&&".format(l,v,v,m)+cuts['relIso'][v]+"&&"+cuts['fsr'][v][l]+"&&!("+cuts['mva'][v]+")", "59740", 20, 0, 300)
            hData.Write("data_{}_{}X_X_vertex{}_d0".format(l,v,m))
            hData.SetLineColor(ROOT.kBlack)
            if hData.Integral() > 0:
                hData.Scale(1.0/hData.Integral())
            temp = bkgStack.drawStack("{}X_X_vertex{}_d0".format(v,m), "HLT_ISO{}&&n{}XX<1&&{}X_X_vertex{}_valid&&".format(l,v,v,m)+cuts['relIso'][v]+"&&"+cuts['fsr'][v][l]+"&&!("+cuts['mva'][v]+")", "59740", 20, 0, 300)
            temp['stack'].Write("mcInverted_{}_{}X_X_vertex{}_d0".format(l,v,m))
            hMC = normalizeStack(temp['stack'])
            c = ROOT.TCanvas("c", "")
            hData.Draw("P")
            hMC.Draw("HIST,same")
            temp['legend'].AddEntry(hData, "Data", "p")
            temp['legend'].Draw("same")
            c.SaveAs("{}X_{}_vertex{}_relIso_invertedMVA_misID.png".format(v,l,m))

            temp = bkgStack.drawStack("{}X_X_vertex{}_d0".format(v,m), "HLT_ISO{}&&n{}XX<1&&{}X_X_vertex{}_valid&&".format(l,v,v,m)+cuts['default'][v][l], "59740", 20, 0, 300)
            temp['stack'].Write("mc_{}_{}X_X_vertex{}_d0".format(l,v,m))
            hMC = normalizeStack(temp['stack'])
            c = ROOT.TCanvas("c", "")
            hData.Draw("P")
            hMC.Draw("HIST,same")
            temp['legend'].AddEntry(hData, "Data", "p")
            temp['legend'].Draw("same")
            c.SaveAs("{}X_{}_vertex{}_relIso_MVA_misID.png".format(v,l,m))


f = ROOT.TFile("plots_fakeRate.root", "RECREATE")

for l in ['MU', 'ELE']:
#    zstack = XStack.drawStack("ZX_Z_mass", "HLT_ISO{}&&".format(l)+cuts['default']['Z'][l], "59740", 45, 50, 140)
#    zstack['canvas'].SaveAs("plots_07_27_20/ZX_Z_mass_{}_default.pdf".format(l))
#    wstack = XStack.drawStack("WX_W_mt", "HLT_ISO{}&&".format(l)+cuts['default']['W'][l], "59740", 25, 0, 200)
#    wstack['canvas'].SaveAs("plots_07_27_20/WX_W_mt_{}_default.pdf".format(l))
'''        
