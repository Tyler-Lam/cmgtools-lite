import math
import ROOT

class xVertex(object):
    def __init__(self, v1, v2, e1, e2, mass):
        self.v1 = v1
        self.v2 = v2
        self.e1 = e1
        self.e2 = e2
        self.mass = mass
        self.valid = 1
        self.vertex = ROOT.TVector3(-999,-999,-999)
        self.pt = 999
        self.setVertex()

    # Given two ecal photon positions, return angle between photon plane and xy plane
    def getRotAngle(self):
        v1 = self.v1
        v3 = v1.Cross(self.v2)
        v3.SetMag(1.0)
        return math.acos(v3[2])

    # Get axis of rotation for photon - xy plane rotation
    def getRotAxis(self):
        v1 = self.v1
        v3 = v1.Cross(self.v2)
        v3.SetMag(1.0)
        ez = ROOT.TVector3(0,0,1.0)
        return v3.Cross(ez)

    # Given mass and two ecal energies, find inscribed angle phi
    # If kinematically impossible, set valid to 0 and return pi
    def getPhi(self):
        m = self.mass
        e1 = self.e1
        e2 = self.e2
        if abs(1.0-m**2/(2*e1*e2))>=1:
            self.valid = 0
            return math.pi
        return math.acos(1.0-m**2/(2*e1*e2))

    ########FOLLOWING METHODS ARE FOR 2D WHEN GETTING VERTEX, AFTER WE ROTATE THE ECAL VECTORS###################

    # Get radius of circles given inscribed angle and two points
    def getRadius(self, x1,x2,y1,y2,phi):
        return math.sqrt(((x1-x2)/2)**2+((y1-y2)/2)**2)/math.sin(phi)

    # Get center of circles given inscribed angle and two points
    def getCenters(self, x1,x2,y1,y2,phi):
        r =self.getRadius(x1,x2,y1,y2,phi)
        q = math.sqrt((x2-x1)**2+(y2-y1)**2)
        x3 = (x1+x2)/2
        y3 = (y1+y2)/2
        xp = x3+math.sqrt(r**2-(q/2)**2)*(y1-y2)/q
        yp = y3+math.sqrt(r**2-(q/2)**2)*(x2-x1)/q
        xm = x3-math.sqrt(r**2-(q/2)**2)*(y1-y2)/q
        ym = y3-math.sqrt(r**2-(q/2)**2)*(x2-x1)/q
        return [[xp,yp],[xm,ym]]

    # Return pt of two points
    def getPt(self, x,y,x1,x2,y1,y2):
        e1 = self.e1
        e2 = self.e2
        v1 = ROOT.TVector3(x1,y1,0)
        v2 = ROOT.TVector3(x2,y2,0)
        v = ROOT.TVector3(x,y,0)
        v2 = v2-v
        v1 = v1-v
        sinTheta1 = (v.Cross(v1).Mag())/(v.Mag()*v1.Mag())*math.copysign(1,v.Cross(v1)[2])
        sinTheta2 = (v.Cross(v2).Mag())/(v.Mag()*v2.Mag())*math.copysign(1,v.Cross(v2)[2])
        pt1 = e1*sinTheta1
        pt2 = e2*sinTheta2
        return pt1+pt2


    # Get phi from vertex of (x,y) going to (x1,y1) (x2,y2) - to remove points directly between ecal hits
    def getPhiFromPoints(self, x,y,x1,x2,y1,y2):
        cosPhi = ((x1-x)*(x2-x)+(y1-y)*(y2-y))/(math.sqrt((x1-x)**2+(y1-y)**2)*math.sqrt((x2-x)**2+(y2-y)**2))
        phi = math.acos(cosPhi)
        return phi

    # Check if a vertex is between beamline and ecal hits
    def checkValid(self, vx, vy, x1, y1, x2, y2):
        x = (x1+x2)/2.
        y = (y1+y2)/2.
        return ((x**2+y**2) - (vx**2+vy**2)>-1 and (x**2+y**2) - ((vx-x)**2+(vy-y)**2))>-1


    # Put everything together, set the vertex, pt, and valid
    def setVertex(self):
        if self.v1[0]==self.v2[0] and self.v1[1]==self.v2[1] and self.v1[2]==self.v2[2]:
            self.valid = 0
            return None
        v1 = self.v1
        v2 = self.v2
        axis = self.getRotAxis()
        theta = self.getRotAngle()
        v1.Rotate(theta, axis)
        v2.Rotate(theta, axis)
        phi = self.getPhi()
        if not self.valid:
            return None
        radius = self.getRadius(v1[0], v2[0], v1[1], v2[1], phi)
        centers = self.getCenters(v1[0], v2[0], v1[1], v2[1], phi)
        stepsPhi = 1000
        deltaPhi = 2*math.pi/stepsPhi
        points = []
        c = min(centers, key = lambda x: x[0]**2+x[1]**2)
        for i in range(stepsPhi):
            angle = i*deltaPhi
            x = c[0]+radius*math.cos(angle)
            y = c[1]+radius*math.sin(angle)
            if abs(self.getPhiFromPoints(x,y,v1[0],v2[0],v1[1],v2[1]) - phi) > 0.001:
                continue
            points.append([x,y])
        goodPoints = filter(lambda x: self.checkValid(x[0], x[1], v1[0], v1[1], v2[0], v2[1]), points)
        if len(goodPoints) == 0:
            self.valid = 0
            return None
        best = min(goodPoints, key = lambda x: abs(self.getPt(x[0], x[1],v1[0],v2[0],v1[1],v2[1])))
        coord = ROOT.TVector3(best[0], best[1], 0)
        coord.Rotate(-theta, axis)
        self.vertex = coord
        self.pt = self.getPt(best[0], best[1], v1[0], v2[0], v1[1], v2[1])

    
