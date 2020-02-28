import math
import ROOT
from PhysicsTools.HeppyCore.utils.deltar import deltaR, deltaPhi
from CMGTools.DisplacedDiPhotons.analyzers.xVertex import *

class PhotonPair(object):
    def __init__(self,leg1,leg2,pdg = 0):
        self.leg1 = leg1
        self.leg2 = leg2
        self.pdg = pdg
        self.LV = leg1.p4(2)+leg2.p4(2)
        et1 = math.sqrt(leg1.mass()*leg1.mass()+leg1.pt()*leg1.pt())
        et2 = math.sqrt(leg2.mass()*leg2.mass()+leg2.pt()*leg2.pt())
        #self.MT2  =self.leg1.p4(2).mass()*self.leg1.p4(2).mass()+\
        #                    self.leg2.p4(2).mass()*self.leg2.p4(2).mass()+2*(et1*et2-self.leg1.p4(2).px()*self.leg2.p4(2).px()-self.leg1.p4(2).py()*self.leg2.p4(2).py())
        self.vertex10 = xVertex(ROOT.TVector3(leg1.caloPosition().x(), leg1.caloPosition().y(),leg1.caloPosition().z()),
                           ROOT.TVector3(leg2.caloPosition().x(), leg2.caloPosition().y(),leg2.caloPosition().z()),
                           leg1.p4(2).energy(),
                           leg2.p4(2).energy(),
                           10)
        self.vertex20 = xVertex(ROOT.TVector3(leg1.caloPosition().x(), leg1.caloPosition().y(),leg1.caloPosition().z()),
                           ROOT.TVector3(leg2.caloPosition().x(), leg2.caloPosition().y(),leg2.caloPosition().z()),
                           leg1.p4(2).energy(),
                           leg2.p4(2).energy(),
                           20)
        self.vertex30 = xVertex(ROOT.TVector3(leg1.caloPosition().x(), leg1.caloPosition().y(),leg1.caloPosition().z()),
                           ROOT.TVector3(leg2.caloPosition().x(), leg2.caloPosition().y(),leg2.caloPosition().z()),
                           leg1.p4(2).energy(),
                           leg2.p4(2).energy(),
                           30)

    def rawP4(self):
        return self.leg1.p4(2)+self.leg2.p4(2)

    def p4(self):
        return self.LV
    
    def m(self):
        return self.LV.mass()
    
    def pdgId(self):
        return self.pdg
   
#    def mt2(self):
#        return self.MT2

#    def mt(self):
#        return self.MT

    def deltaPhi(self):
        return abs(deltaPhi(self.leg1.phi(),self.leg2.phi()))

    def deltaR(self):
        return abs(deltaR(self.leg1.eta(),self.leg1.phi(),self.leg2.eta(),self.leg2.phi()))

    # Vertix objects for various X mass values
    def vertex(self, mass):
        v1 = ROOT.TVector3(self.leg1.caloPosition().x(), self.leg1.caloPosition().y(),self.leg1.caloPosition().z())
        v2 = ROOT.TVector3(self.leg2.caloPosition().x(), self.leg2.caloPosition().y(),self.leg2.caloPosition().z())
        e1 = self.leg1.p4(2).energy()
        e2 = self.leg2.p4(2).energy()
        vertex = xVertex(v1, v2, e1, e2, mass)
        return vertex

    def __getattr__(self, name):
        return getattr(self.LV,name)
                       
