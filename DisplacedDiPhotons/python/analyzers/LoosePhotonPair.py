import math
import ROOT
from PhysicsTools.HeppyCore.utils.deltar import deltaR, deltaPhi
from CMGTools.DisplacedDiPhotons.analyzers.xVertex import *

class LoosePhotonPair(object):
    def __init__(self, leg1, leg2, pdg = 0):
        self.leg1 = leg1
        self.leg2 = leg2
        self.pdg = pdg
        self.LV = leg1.p4() + leg2.p4()

        self.vertex10 = self.vertex(10)
        self.vertex15 = self.vertex(15)
        self.vertex20 = self.vertex(20)
        self.vertex30 = self.vertex(30)
        self.vertex40 = self.vertex(40)
        self.vertex50 = self.vertex(50)
        self.vertex60 = self.vertex(60)

    def p4(self):
        return self.LV

    def m(self):
        return self.LV.mass()

    def pdgId(self):
        return self.pdg

    def deltaPhi(self):
        return abs(deltaPhi(self.leg1.phi(), self.leg2.phi()))

    def deltaR(self):
        return abs(deltaR(self.leg1.eta(), self.leg1.phi(), self.leg2.eta(), self.leg2.phi()))

    def vertex(self, mass):
        vertexCalculator = ROOT.cmg.VertexCalculator()
        vertexCalculator.run(self.leg1.caloPosition, self.leg2.caloPosition, self.leg1.energy(), self.leg2.energy(), mass)
        return xVertex(vertexCalculator.vertex(), vertexCalculator.pt(), vertexCalculator.phi(), vertexCalculator.d0(), vertexCalculator.ip3d(), vertexCalculator.valid())

    def __getattr__(self, name):
        return getattr(self.LV, name)
