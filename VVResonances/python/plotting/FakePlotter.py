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

        # Separating g1 and g2 histograms
        self.num_g1 = self.fakefile.Get("num_g1")
        self.denom_g1 = self.fakefile.Get("denom_g1")
        fakerate_g1 = self.num_g1.Clone("fakerate_g1")
        fakerate_g1.Divide(self.denom_g1)
        self.fakerate_g1 = fakerate_g1
        
        self.num_g2 = self.fakefile.Get("num_g2")
        self.denom_g2 = self.fakefile.Get("denom_g2")
        fakerate_g2 = self.num_g2.Clone("fakerate_g2")
        fakerate_g2.Divide(self.denom_g2)
        self.fakerate_g2 = fakerate_g2



    # NOTE: This only works for plotting VX variables
    def drawTH1(self, var, cuts, lumi, bins, minx, maxx, titlex = "", units = "", drawStyle="HIST"):
        print "Drawing fakerate for {}".format(self.name)
        h = ROOT.TH1D("tmpTH1", "", bins, minx, maxx)
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

        #Check to remove cuts if not common
        for n,x in enumerate(split):
            for b in branchNames:
                if "loose" in b:
                    continue                
                temp = b
                temp = temp.replace("WX", "looseWX")
                temp = temp.replace("ZX", "looseZX")
                stripped = x
                # Inelegant way to fix parenthesis for weird cases
                # e.g. (cut1)*((cut2)*(cut3))
                #      (1)*(1)*(cut3)) <- has extra ')'
                while stripped[0]=="(" and stripped[-1]==")":
                    stripped = stripped[1:-1]
                if stripped[0] == "(" and stripped[-1]!=")":
                    stripped = stripped[1:]
                if stripped[0] != "(" and stripped[-1]==")":
                    stripped = stripped[0:-1]
                if b in x:
                    if temp not in branchNames:
                        split[n] = split[n].replace(stripped, "1")
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
        fakerate_g1 = self.fakerate_g1
        fakerate_g2 = self.fakerate_g2
        for e in self.tree:
            if len(getattr(e, var))==0:
                continue
            v = getattr(e, var)[0]
            #Get fakerate

            binPt1 = fakerate.GetXaxis().FindBin(min(90,pt1.EvalInstance()))
            binEta1 = fakerate.GetYaxis().FindBin(eta1.EvalInstance())
            binPt2 = fakerate.GetXaxis().FindBin(min(90,pt2.EvalInstance()))
            binEta2 = fakerate.GetYaxis().FindBin(eta2.EvalInstance())
            
            fr1 = fakerate.GetBinContent(binPt1, binEta1)
            fr2 = fakerate.GetBinContent(binPt2, binEta2)
            
            #binPt1 = fakerate_g1.GetXaxis().FindBin(min(90,pt1.EvalInstance()))
            #binEta1 = fakerate_g1.GetYaxis().FindBin(eta1.EvalInstance())
            #binPt2 = fakerate_g2.GetXaxis().FindBin(min(90,pt2.EvalInstance()))
            #binEta2 = fakerate_g2.GetYaxis().FindBin(eta2.EvalInstance())


            #fr1 = fakerate_g1.GetBinContent(binPt1, binEta1)
            #fr2 = fakerate_g2.GetBinContent(binPt2, binEta2)
            
#            h.Fill(v, weights.EvalInstance())
            h.Fill(v, (fr1+fr2-fr1*fr2)*weights.EvalInstance())

        return h

