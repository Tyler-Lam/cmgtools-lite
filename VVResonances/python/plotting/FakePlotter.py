import ROOT
import sys
from array import array
import pickle
from TreePlotter import TreePlotter
from array import array


class FakePlotter(TreePlotter):

    def __init__(self, file, tree, fakerate):
        super(FakePlotter,self).__init__(file, tree)
        self.fakefile = ROOT.TFile(fakerate)
        self.name = fakerate
        self.num = self.fakefile.Get("num")
        self.denom = self.fakefile.Get("denom")
        fakerate = self.num.Clone("fakerate")
        fakerate.Divide(self.denom)
        self.fakerate = fakerate

    # NOTE: This only works if cuts use branches that exist for reco AND pf photons
    def drawTH1(self, var, cuts, lumi, bins, min, max, titlex = "", units = "", drawStyle="HIST"):
        print "Drawing fakerate for {}".format(self.name)
        h = ROOT.TH1D("tmpTH1", "", bins, min, max)
        h.Sumw2()
        h.SetLineStyle(self.linestyle)
        h.SetLineColor(self.linecolor)
        h.SetLineWidth(self.linewidth)
        h.SetFillStyle(self.fillstyle)
        h.SetFillColor(self.fillcolor)
        h.SetMarkerStyle(self.markerstyle)
        if units=="":
            h.GetXaxis().SetTitle(titlex)
        else:
            h.GetXaxis().SetTitle(titlex+ " ["+units+"]")
        
        corrString = '1'
        for corr in self.corrFactors:
            corrString = corrString + "*("+str(corr['value'])+")"
            
        # Convert reco photons to loose photons
        var = var.replace("WX", "looseWX")
        var = var.replace("ZX", "looseZX")
        split = cuts.split("*")
        branches = self.tree.GetListOfBranches()
        branchNames = []
        for b in branches:
            branchNames.append(b.GetName())
        for n,x in enumerate(split):
            for b in branchNames:
                if "loose" in b:
                    continue                
                temp = b
                temp = temp.replace("WX", "looseWX")
                temp = temp.replace("ZX", "looseZX")
                if b in x:
                    if temp not in branchNames:
                        split[n] = "(1)"
                        break
        cuts = "*".join(split)
        cuts = cuts.replace("WX", "looseWX")
        cuts = cuts.replace("ZX", "looseZX")           
        pt1 = ROOT.TTreeFormula("pt1", var[:7]+"_X_g1_pt", self.tree)
        pt2 = ROOT.TTreeFormula("pt2", var[:7]+"_X_g2_pt", self.tree)
        eta1 = ROOT.TTreeFormula("eta1", var[:7]+"_X_g1_eta", self.tree)
        eta2 = ROOT.TTreeFormula("eta2", var[:7]+"_X_g2_eta", self.tree)
        weights = ROOT.TTreeFormula("weights", "("+cuts+")*"+lumi+"*"+self.weight+"*("+corrString+")", self.tree)

        fakerate = self.fakerate
        nEntries = -1
        for e in self.tree:
            nEntries += 1
            if nEntries%10000==0:
                print "Processed {} entries".format(nEntries)
            if len(getattr(e, var))==0:
                continue
            v = getattr(e, var)[0]
            #Get fakerate stuff
            fr1 = fakerate.GetBinContent(fakerate.GetXaxis().FindBin(pt1.EvalInstance()), fakerate.GetYaxis().FindBin(eta1.EvalInstance()))
            fr2 = fakerate.GetBinContent(fakerate.GetXaxis().FindBin(pt2.EvalInstance()), fakerate.GetYaxis().FindBin(eta2.EvalInstance()))
            h.Fill(v, fr1*fr2*weights.EvalInstance())
            '''
            print "fr1={} fr2={} weights={}".format(fr1, fr2, weights.EvalInstance())
            print "var = {}".format(var)
            print "weights = {}".format(weights)
            print "corrFactors={}".format("("+cuts+")*"+lumi+"*"+self.weight+"*("+corrString+")")
            import pdb
            pdb.set_trace()
            '''
        return h

