from PhysicsTools.Heppy.analyzers.core.Analyzer import Analyzer
from PhysicsTools.Heppy.analyzers.core.AutoHandle import AutoHandle
from PhysicsTools.HeppyCore.utils.deltar import *


import ROOT
import os
import math

import  itertools
from CMGTools.VVResonances.tools.Pair import *
from CMGTools.DisplacedDiPhotons.analyzers.PhotonPair import *
from CMGTools.DisplacedDiPhotons.analyzers.XPair import *
from CMGTools.DisplacedDiPhotons.analyzers.xVertex import *
from CMGTools.DisplacedDiPhotons.analyzers.PFPhoton import *
from CMGTools.DisplacedDiPhotons.analyzers.LoosePhotonPair import *
from CMGTools.DisplacedDiPhotons.analyzers.LooseXPair import *

debug = False

class VHGGBuilder(Analyzer):

    def __init__(self, cfg_ana, cfg_comp, looperName):
        super(VHGGBuilder, self).__init__(cfg_ana, cfg_comp, looperName)


    def declareHandles(self):
        super(VHGGBuilder, self).declareHandles()
        self.handles['packed'] = AutoHandle('packedPFCandidates', 'std::vector<pat::PackedCandidate>')

    def makeW(self,leptons,met):
        Ws=[]
        # Pair leptons to met
        for l in leptons:
            W=Pair(l,met,l.charge()*24)
            Ws.append(W)
        return Ws

    def makeZ(self,leptons):
        Zs=[]
        # Pair leptons if near Z mass
        for l1,l2 in itertools.combinations(leptons,2):
            if (l1.charge()+l2.charge())!=0:
                continue
            # Did this need to be added? Wasn't in and couldn't find any instances where it wasn't true
            if l1.pdgId()+l2.pdgId()!=0:
                continue
            Z = Pair(l1,l2,23)
            mass = Z.p4().mass()
            if mass<50 or mass>140:
                continue
            Zs.append(Z)
        return Zs

    def makeX(self,photons):
        Xs=[]
        # Pair photons (TODO: Conversions)
        for g1,g2 in itertools.combinations(photons,2):
            X = PhotonPair(g1,g2,1995)
            Xs.append(X)
            
        return Xs

    def makeXPair(self, photons):
        XXs = []
        # Loop over 4 photons to make XXs
        for g1,g2,g3,g4 in itertools.combinations(photons, 4):
            X12 = PhotonPair(g1, g2)
            X34 = PhotonPair(g3, g4)
            X13 = PhotonPair(g1, g3)
            X24 = PhotonPair(g2, g4)
            X23 = PhotonPair(g2, g3)
            X14 = PhotonPair(g1, g4)
            
            X1234 = XPair(X12, X34)
            X1324 = XPair(X13, X24)
            X1423 = XPair(X14, X23)
            XXs.append(X1234)
            XXs.append(X1324)
            XXs.append(X1423)
        return XXs

    def makeLooseX(self,photons):
        Xs=[]
        # Pair photons (TODO: Conversions)
        if len(photons) < 2:
            return Xs
        # Pick highest pt pf photons
        sortedPhotons = sorted(photons, key = lambda x: x.pt(), reverse = True)
        Xs.append(LoosePhotonPair(sortedPhotons[0], sortedPhotons[1]))
        return Xs

    def makeLooseXPair(self, photons):
        XXs = []
        # Loop over 4 photons to make XXs
        if len(photons) < 4:
            return XXs
        sortedPhotons = sorted(photons, key = lambda x: x.pt(), reverse = True)
        # Pick highest 4 pt pf photons
        for g1,g2,g3,g4 in itertools.combinations(sortedPhotons[0:3], 4):
            X12 = LoosePhotonPair(g1, g2)
            X34 = LoosePhotonPair(g3, g4)
            X13 = LoosePhotonPair(g1, g3)
            X24 = LoosePhotonPair(g2, g4)
            X23 = LoosePhotonPair(g2, g3)
            X14 = LoosePhotonPair(g1, g4)
            
            X1234 = LooseXPair(X12, X34)
            X1324 = LooseXPair(X13, X24)
            X1423 = LooseXPair(X14, X23)
            XXs.append(X1234)
            XXs.append(X1324)
            XXs.append(X1423)
        return XXs

    #Check fsr for ZX candidates
    def checkFSR_ZX(self, Z,X):
        flag = 0
        mz = Z.p4().mass()
        for g in [X.leg1, X.leg2]:
            p4New = Z.p4() + g.p4(2)
            if abs(p4New.mass()-91) < abs(mz-91) and abs(p4New.mass()-91) < 15:
                flag = 1

        p4Total = Z.p4() + X.p4()
        if abs(p4Total.mass()-91) < abs(mz-91) and abs(p4New.mass()-91) < 15:
            flag = 1

        return flag
        
    # Check fsr for loose (pf) ZX candidates
    def checkFSR_LooseZX(self, Z, X):
        flag = 0
        mz = Z.p4().mass()
        for g in [X.leg1, X.leg2]:
            p4New = Z.p4() + g.p4()
            if abs(p4New.mass()-91) < abs(mz-91) and abs(p4New.mass()-91) < 15:
                flag = 1

        p4Total = Z.p4() + X.p4()
        if abs(p4Total.mass()-91) < abs(mz-91) and abs(p4New.mass()-91) < 15:
            flag = 1

        return flag


    def checkFSR_ZXX(self, Z,XX):
        flag = 0
        mz = Z.p4().mass()
        for g in [XX.x1.leg1, XX.x1.leg2, XX.x2.leg1, XX.x2.leg2]:
            p4New = Z.p4() + g.p4(2)
            if abs(p4New.mass()-91) < abs(mz-91) and abs(p4New.mass()-91) < 15:
                flag =  1
        for g1, g2 in itertools.combinations([XX.x1.leg1, XX.x1.leg2, XX.x2.leg1, XX.x2.leg2],2):
            p4New = Z.p4() + g1.p4(2) + g2.p4(2)
            if abs(p4New.mass()-91) < abs(mz-91) and abs(p4New.mass()-91) < 15:
                flag = 1
        for g1, g2, g3 in itertools.combinations([XX.x1.leg1, XX.x1.leg2, XX.x2.leg1, XX.x2.leg2],3):
            p4New = Z.p4() + g1.p4(2) + g2.p4(2) + g3.p4(2)
            if abs(p4New.mass()-91) < abs(mz-91) and abs(p4New.mass()-91) < 15:
                flag = 1
        for g1, g2, g3, g4 in itertools.combinations([XX.x1.leg1, XX.x1.leg2, XX.x2.leg1, XX.x2.leg2],4):
            p4New = Z.p4() + g1.p4(2) + g2.p4(2) + g3.p4(2) + g4.p4(2)
            if abs(p4New.mass()-91) < abs(mz-91) and abs(p4New.mass()-91) < 15:
                flag = 1

        return flag


    def checkFSR_LooseZXX(self, Z,XX):
        flag = 0
        mz = Z.p4().mass()
        for g in [XX.x1.leg1, XX.x1.leg2, XX.x2.leg1, XX.x2.leg2]:
            p4New = Z.p4() + g.p4()
            if abs(p4New.mass()-91) < abs(mz-91) and abs(p4New.mass()-91) < 15:
                flag =  1
        for g1, g2 in itertools.combinations([XX.x1.leg1, XX.x1.leg2, XX.x2.leg1, XX.x2.leg2],2):
            p4New = Z.p4() + g1.p4() + g2.p4()
            if abs(p4New.mass()-91) < abs(mz-91) and abs(p4New.mass()-91) < 15:
                flag = 1
        for g1, g2, g3 in itertools.combinations([XX.x1.leg1, XX.x1.leg2, XX.x2.leg1, XX.x2.leg2],3):
            p4New = Z.p4() + g1.p4() + g2.p4() + g3.p4()
            if abs(p4New.mass()-91) < abs(mz-91) and abs(p4New.mass()-91) < 15:
                flag = 1
        for g1, g2, g3, g4 in itertools.combinations([XX.x1.leg1, XX.x1.leg2, XX.x2.leg1, XX.x2.leg2],4):
            p4New = Z.p4() + g1.p4() + g2.p4() + g3.p4() + g4.p4()
            if abs(p4New.mass()-91) < abs(mz-91) and abs(p4New.mass()-91) < 15:
                flag = 1

        return flag


    def makeZXs(self, bestZ, Xs):
        ZXs = []
        #Pair best Z to best X:
        if len(Xs) > 0:
            #Select best X by balancing momentum with lepton
            bestX = min(Xs, key = lambda x: (x.p4() + bestZ.p4()).pt())
            bestZX = Pair(bestZ, bestX)
            bestZX.hasFSR = self.checkFSR_ZX(bestZ, bestX)
            bestZX.deltaPhi_g1 = min(abs(deltaPhi(bestZ.leg1.phi(), bestX.leg1.phi())), abs(deltaPhi(bestZ.leg2.phi(), bestX.leg1.phi())))
            bestZX.deltaPhi_g2 = min(abs(deltaPhi(bestZ.leg1.phi(), bestX.leg2.phi())), abs(deltaPhi(bestZ.leg2.phi(), bestX.leg2.phi())))
            ZXs.append(bestZX)
            
        return ZXs

    def makeLooseZXs(self, bestZ, Xs):
        ZXs = []
        #Pair best Z to best X:
        if len(Xs) > 0:
            #Select best X by balancing momentum with lepton
            bestX = min(Xs, key = lambda x: (x.p4() + bestZ.p4()).pt())
            bestZX = Pair(bestZ, bestX)
            bestZX.hasFSR = self.checkFSR_LooseZX(bestZ, bestX)
            bestZX.deltaPhi_g1 = min(abs(deltaPhi(bestZ.leg1.phi(), bestX.leg1.phi())), abs(deltaPhi(bestZ.leg2.phi(), bestX.leg1.phi())))
            bestZX.deltaPhi_g2 = min(abs(deltaPhi(bestZ.leg1.phi(), bestX.leg2.phi())), abs(deltaPhi(bestZ.leg2.phi(), bestX.leg2.phi())))
            ZXs.append(bestZX)
            
        return ZXs


    def makeZXXs(self, bestZ, XXs):
        ZXXs = []
        if len(XXs) > 0:
            # Look at z + 4gamma combinations that minimize pt
            sortedXXs = sorted(XXs, key = lambda x: (x.p4()+bestZ.p4()).pt())[0:3]
            # Pick best XX by minimizing deltaR between photon pairs
            bestXX = min(sortedXXs, key = lambda x: deltaR(x.x1.leg1.eta(),x.x1.leg1.phi(),x.x1.leg2.eta(),x.x1.leg2.phi()) + deltaR(x.x2.leg1.eta(), x.x2.leg1.phi(), x.x2.leg2.eta(), x.x2.leg2.phi()))
            bestZXX = ZXX(bestZ, bestXX)
            bestZXX.hasFSR = self.checkFSR_ZXX(bestZ, bestXX)
            ZXXs.append(bestZXX)
        return ZXXs


    def makeLooseZXXs(self, bestZ, XXs):
        ZXXs = []
        if len(XXs) > 0:
            sortedXXs = sorted(XXs, key = lambda x: (x.p4()+bestZ.p4()).pt())[0:3]
            bestXX = min(sortedXXs, key = lambda x: deltaR(x.x1.leg1.eta(),x.x1.leg1.phi(),x.x1.leg2.eta(),x.x1.leg2.phi()) + deltaR(x.x2.leg1.eta(), x.x2.leg1.phi(), x.x2.leg2.eta(), x.x2.leg2.phi()))
            bestZXX = ZXX(bestZ, bestXX)
            bestZXX.hasFSR = self.checkFSR_LooseZXX(bestZ, bestXX)
            ZXXs.append(bestZXX)
        return ZXXs


    def makeWXs(self, bestW, Xs):
        WXs = []
        if len(Xs) > 0:
            bestX = min(Xs, key = lambda x: (x.p4() + bestW.p4()).pt())
            bestWX = Pair(bestW, bestX)
            bestWX.deltaPhi_g1 = deltaPhi(bestW.leg1.phi(), bestX.leg1.phi())
            bestWX.deltaPhi_g2 = deltaPhi(bestW.leg1.phi(), bestW.leg2.phi())
            misID = 0
            masses = {}
            if abs(bestW.leg1.pdgId()) == 11:
                p4_1 = bestW.leg1.p4() + bestX.leg1.p4(2)
                masses[1] = p4_1.mass()
                p4_2 = bestW.leg1.p4() + bestX.leg2.p4(2)
                masses[2] = p4_2.mass()
                p4_3 = bestW.leg1.p4() + bestX.leg1.p4(2) + bestX.leg2.p4(2)
                masses[3] = p4_3.mass()
                misID = min(masses, key = lambda x: abs(masses[x] - 90))
                if abs(masses[misID] - 90) > 10:
                        misID = 0
            bestWX.misID = misID
            WXs.append(bestWX)
        return WXs

    def makeLooseWXs(self, bestW, Xs):
        WXs = []
        if len(Xs) > 0:
            bestX = min(Xs, key = lambda x: (x.p4() + bestW.p4()).pt())
            bestWX = Pair(bestW, bestX)
            bestWX.deltaPhi_g1 = deltaPhi(bestW.leg1.phi(), bestX.leg1.phi())
            bestWX.deltaPhi_g2 = deltaPhi(bestW.leg1.phi(), bestW.leg2.phi())
            misID = 0
            masses = {}
            if abs(bestW.leg1.pdgId()) == 11:
                p4_1 = bestW.leg1.p4() + bestX.leg1.p4()
                masses[1] = p4_1.mass()
                p4_2 = bestW.leg1.p4() + bestX.leg2.p4()
                masses[2] = p4_2.mass()
                p4_3 = bestW.leg1.p4() + bestX.leg1.p4() + bestX.leg2.p4()
                masses[3] = p4_3.mass()
                misID = min(masses, key = lambda x: abs(masses[x] - 90))
                if abs(masses[misID] - 90) > 10:
                    misID = 0
            bestWX.misID = misID
            WXs.append(bestWX)
        return WXs



    def makeWXXs(self, bestW, XXs):
        WXXs = []
        if len(XXs) > 0:
            sortedXXs = sorted(XXs, key = lambda x: (x.p4()+bestW.p4()).pt())[0:3]
            bestXX = min(sortedXXs,  key = lambda x: deltaR(x.x1.leg1.eta(),x.x1.leg1.phi(),x.x1.leg2.eta(),x.x1.leg2.phi()) + deltaR(x.x2.leg1.eta(), x.x2.leg1.phi(), x.x2.leg2.eta(), x.x2.leg2.phi()))
            bestWXX = WXX(bestW, bestXX)
            bestWXX.deltaPhi_X1_g1 = deltaPhi(bestW.leg1.phi(), bestXX.x1.leg1.phi())
            bestWXX.deltaPhi_X1_g2 = deltaPhi(bestW.leg1.phi(), bestXX.x1.leg2.phi())
            bestWXX.deltaPhi_X2_g1 = deltaPhi(bestW.leg1.phi(), bestXX.x2.leg1.phi())
            bestWXX.deltaPhi_X2_g2 = deltaPhi(bestW.leg1.phi(), bestXX.x2.leg2.phi())
            WXXs.append(bestWXX)
        return WXXs

        
    def log(self, signal, genLeptons, genPhotons, recoLeptons, recoPhotons):
        for S in signal:
            print "S vertex=({:.3f}, {:.3f}, {:.3f}), pt={:.3f}, eta={:.3f}, phi={:.3f}, ID={}".format(S.vx(),S.vy(),S.vz(),S.pt(),S.eta(),S.phi(), S.pdgId())
            for d in range(S.numberOfDaughters()):
                print "    Daughter vertex=({:.3f}, {:.3f}, {:.3f}), pt={:.3f}, eta={:.3f}, phi={:.3f}, ID={}".format(S.daughter(d).vx(),S.daughter(d).vy(),S.daughter(d).vz(),S.daughter(d).pt(),S.daughter(d).eta(),S.daughter(d).phi(), S.daughter(d).pdgId())
        print "Gen Leptons"
        for l in genLeptons:
            print "    pt={:.3f}, eta={:.3f}, phi={:.3f}, ID={:.3f}, status={}".format(l.pt(), l.eta(), l.phi(), l.pdgId(), l.status())
        print "Gen Photons"
        for g in genPhotons:
            print "    pt={:.3f}, eta={:.3f}, phi={:.3f}, mother={}, status={}, vertex=({:.3f},{:.3f},{:.3f})".format(g.pt(),g.eta(),g.phi(),g.mother().pdgId(),g.status(),g.vx(),g.vy(),g.vz())
        print "Reco photons"
        for p in recoPhotons:
            print "   ",
            print p,
            print "mcMatchId={}, mcMotherId={}, MVANonTrigV1Values={:.3f}, hasPixelSeed={}, passElectronVeto={}".format(p.mcMatchId,p.mcMotherId,p.userFloat("PhotonMVAEstimatorRun2Spring16NonTrigV1Values"),p.hasPixelSeed(), p.passElectronVeto())
        print "Reco leptons"
        for l in recoLeptons:
            print "   ",
            print l

        
    def process(self, event):

        self.readCollections(event.input)

        pfCands = self.handles['packed'].product()
        pfPhotons = filter(lambda x: x.pdgId()==22 and x.pt() > 10 and abs(x.eta()) < 2.5, pfCands)
        
        event.ZX=[]
        event.WX=[]
        event.ZXX = []
        event.WXX = []
        event.looseZX = []
        event.looseWX = []
        event.looseZXX = []
        event.looseWXX = []
        leptons = filter(lambda x: x.pt() > 0, event.selectedLeptons)
        goodLeptons = filter(lambda x: x.relIso03 < .1, leptons)
        photons = filter(lambda x: x.pt()>0,event.selectedPhotons)
        
        #Get Gen info
        gen = []
        if self.cfg_comp.isMC:
            gen = event.genParticles
        genPhotons = filter(lambda x: x.pdgId() == 22 and x.status()==1, gen)
        genWs = filter(lambda x: x.pdgId() == 24 and x.numberOfDaughters()==2, gen)
        genLeptons = filter(lambda x: abs(x.pdgId())==11 or abs(x.pdgId())==13, gen)
        signal = filter(lambda x: abs(x.pdgId())==9000006, gen)
        signalPhotons = filter(lambda x: abs(x.mother().pdgId())==9000006, genPhotons)
        event.GenPhoton = genPhotons

        goodPhotons = []
        for x in photons:
            overlap = False
            for l in goodLeptons:
                if deltaR(l.eta(), l.phi(), x.eta(), x.phi()) < 0.3:
                    overlap = True
            if not overlap:
                goodPhotons.append(x)

        loosePhotons = []
        for x in pfPhotons:
            g = PFPhoton(x)
            overlap = False
            for l in goodLeptons:
                if deltaR(l.eta(), l.phi(),  x.eta(), x.phi()) < 0.3:
                    overlap = True
                    break
            if not overlap:
                loosePhotons.append(g)
        if self.cfg_comp.isMC:
            for g in loosePhotons:
                for gamma in genPhotons:
                    if deltaR(gamma.eta(), gamma.phi(), g.eta(), g.phi()) < 0.1:
                        g.mcMatch = 1
                        if gamma.numberOfMothers() > 0:
                            g.mcMotherId = gamma.mother().pdgId()
                        continue
            

        # Z to ee, used for photon misID cuts
        Zees = filter(lambda x: abs(x.leg1.pdgId())==11, self.makeZ(leptons))

        Zs = self.makeZ(goodLeptons)
        Ws = self.makeW(goodLeptons,event.met)
        Xs = self.makeX(goodPhotons)
        XXs = self.makeXPair(goodPhotons)
        LooseXs = self.makeLooseX(loosePhotons)
        LooseXXs = self.makeLooseXPair(loosePhotons)

        # Make ZX/ZXX Pairs first
        if len(Zs) > 0:
            bestZ = max(Zs, key = lambda x: x.leg1.pt() + x.leg2.pt())

            ZXs = self.makeZXs(bestZ, Xs)
            for zx in ZXs:
                zx.otherLeptons = len(leptons) - 2
                event.ZX.append(zx)
            
            ZXXs = self.makeZXXs(bestZ, XXs)
            for zxx in ZXXs:
                zxx.otherLeptons = len(leptons) - 2
                event.ZXX.append(zxx)

            LooseZXs = self.makeLooseZXs(bestZ, LooseXs)
            for zx in LooseZXs:
                zx.otherLeptons = len(leptons) - 2
                event.looseZX.append(zx)

            LooseZXXs = self.makeLooseZXXs(bestZ, LooseXXs)
            for zxx in LooseZXXs:
                zxx.otherLeptons = len(leptons) - 2
                event.looseZXX.append(zxx)

        elif len(Ws) > 0:
            bestW = max(Ws, key = lambda x: x.leg1.pt())
            WXs = self.makeWXs(bestW, Xs)
            for wx in WXs:
                wx.otherLeptons = len(leptons) - 1
                wx.hasZee = (len(Zees) > 0)
                event.WX.append(wx)
            
            WXXs = self.makeWXXs(bestW, XXs)
            for wxx in WXXs:
                wxx.otherLeptons = len(leptons) - 1
                wxx.hasZee = (len(Zees) > 0)
                event.WXX.append(wxx)

            LooseWXs = self.makeLooseWXs(bestW, LooseXs)
            for wx in LooseWXs:
                wx.otherLeptons = len(leptons) - 1
                wx.hasZee = (len(Zees) > 0)
                event.looseWX.append(wx)

            LooseWXXs = self.makeWXXs(bestW, LooseXXs)
            for wxx in LooseWXXs:
                wx.otherLeptons = len(leptons) - 1
                wx.hasZee = (len(Zees) > 0)
                event.looseWXXs.append(wxx)
