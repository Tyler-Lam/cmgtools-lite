from PhysicsTools.Heppy.analyzers.core.AutoFillTreeProducer  import * 

import ROOT
dummyLV=ROOT.math.XYZTLorentzVector(0.0,0.0,0.0001,0.0001)

conversionType = NTupleObjectType("conversion",baseObjectTypes = [], variables = [
    NTupleVariable("convX", lambda x: x.conversionVertex().position().x(), float , help = "Reco Conversion Vertex X"),
    NTupleVariable("convY", lambda x: x.conversionVertex().position().y(), float , help = "Reco Conversion Vertex Y"),
    NTupleVariable("convZ", lambda x: x.conversionVertex().position().z(), float , help = "Reco Conversion Vertex Z"),
])


ZType = NTupleObjectType("PairType", baseObjectTypes=[fourVectorType], variables = [
    NTupleVariable("deltaPhi",   lambda x : x.deltaPhi(), float),       
    NTupleVariable("deltaR",   lambda x : x.deltaR(), float),       
    NTupleVariable("mt", lambda x: x.mt(), float),
    NTupleSubObject("l1",  lambda x : x.leg1,leptonType),
    NTupleSubObject("l2",  lambda x : x.leg2,leptonType),    
])

XType = NTupleObjectType("PhotonPair", baseObjectTypes=[fourVectorType], variables = [
    NTupleVariable("deltaPhi",   lambda x : x.deltaPhi(), float),       
    NTupleVariable("deltaR",   lambda x : x.deltaR(), float),       
    NTupleSubObject("g1",  lambda x : x.leg1,photonType),
    NTupleSubObject("g2",  lambda x : x.leg2,photonType),    
])

WType = NTupleObjectType("PairType", baseObjectTypes=[fourVectorType], variables = [
    NTupleVariable("deltaPhi",   lambda x : x.deltaPhi(), float),       
    NTupleVariable("deltaR",   lambda x : x.deltaR(), float),
    NTupleVariable("mt", lambda x: x.mt(), float),
    NTupleSubObject("l1",  lambda x : x.leg1,leptonType),
    NTupleSubObject("l2",  lambda x : x.leg2,metType),    
])


ZXType = NTupleObjectType("ZXType", baseObjectTypes=[], variables = [
    NTupleSubObject("Z",  lambda x : x.leg1,ZType),
    NTupleSubObject("X",  lambda x : x.leg2,XType),
    NTupleVariable("deltaPhi_g1", lambda x: x.deltaPhi_g1, float),
    NTupleVariable("deltaPhi_g2", lambda x: x.deltaPhi_g2, float),
    NTupleVariable("mass", lambda x: x.m(), float),
    NTupleVariable("mass_Zg1", lambda x: (x.leg1.p4()+x.leg2.leg1.p4(2)).mass(), float),
    NTupleVariable("mass_Zg2", lambda x: (x.leg1.p4()+x.leg2.leg2.p4(2)).mass(), float),
    NTupleVariable("hasFSR", lambda x: x.hasFSR, int),
    NTupleVariable("otherLeptons",   lambda x : x.otherLeptons, int),          
])

XXType = NTupleObjectType("XXType", baseObjectTypes=[], variables = [
    NTupleSubObject("X1", lambda x : x.x1, XType),
    NTupleSubObject("X2", lambda x : x.x2, XType),
    NTupleVariable("m12", lambda x: x.m12(), float),
    NTupleVariable("m13", lambda x: x.m13(), float),
    NTupleVariable("m14", lambda x: x.m14(), float),
    NTupleVariable("m23", lambda x: x.m23(), float),
    NTupleVariable("m24", lambda x: x.m24(), float),
    NTupleVariable("m34", lambda x: x.m34(), float),
])

WXType = NTupleObjectType("WXType", baseObjectTypes=[], variables = [
    NTupleSubObject("W",  lambda x : x.leg1,WType),
    NTupleSubObject("X",  lambda x : x.leg2,XType),
    NTupleVariable("deltaPhi_g1", lambda x: x.deltaPhi_g1, float),
    NTupleVariable("deltaPhi_g2", lambda x: x.deltaPhi_g2, float),
    NTupleVariable("otherLeptons", lambda x: x.otherLeptons, int),
])

ZXXType = NTupleObjectType("ZXXType", baseObjectTypes=[], variables = [
    NTupleSubObject("Z", lambda x: x.Z, ZType),
    NTupleSubObject("XX", lambda x: x.XX, XXType),
    NTupleVariable("otherLeptons", lambda x: x.otherLeptons, int),
    NTupleVariable("hasFSR", lambda x: x.hasFSR, int),
])

WXXType = NTupleObjectType("WXXType", baseObjectTypes=[], variables = [
    NTupleSubObject("W", lambda x: x.W, WType),
    NTupleSubObject("XX", lambda x: x.XX, XXType),
    NTupleSubObject("deltaPhi_X1_g1", lambda x: x.deltaPhi_X1_g1, float),
    NTupleSubObject("deltaPhi_X1_g2", lambda x: x.deltaPhi_X1_g1, float),
    NTupleSubObject("deltaPhi_X2_g1", lambda x: x.deltaPhi_X1_g1, float),
    NTupleSubObject("deltaPhi_X2_g2", lambda x: x.deltaPhi_X1_g1, float),
    NTupleVariable("otherLeptons", lambda x: x.otherLeptons, int),
])
