import math
import ROOT

class xVertex(object):
    def __init__(self, vertex, pt, phi, d0, ip3d, valid):
        self.vertex = vertex
        self.pt = pt
        self.phi = phi
        self.d0 = d0
        self.ip3d = ip3d
        self.valid = valid

