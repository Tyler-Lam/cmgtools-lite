#include "CMGTools/DisplacedDiPhotons/interface/VertexCalculator.h"

using namespace cmg;

VertexCalculator::VertexCalculator(){
  //  vertex_ = TVector3(-999,-999,-999);
  // pt_ = -999;
  //phi_ = TMath::Pi();
}

VertexCalculator::~VertexCalculator() {}


void VertexCalculator::setVertex(const TVector3& vertex){vertex_ = vertex;}
void VertexCalculator::setPt(double pt){pt_ = pt;}  
void VertexCalculator::setPhi(double phi){phi_ = phi;}
void VertexCalculator::setD0(double d0){d0_ = d0;}
void VertexCalculator::setIp3d(double ip3d){ip3d_ = ip3d;}
void VertexCalculator::setValid(bool valid){valid_ = valid;}
TVector3 VertexCalculator::vertex(){return vertex_;}    
float VertexCalculator::pt(){return pt_;}
float VertexCalculator::phi(){return phi_;}
float VertexCalculator::d0(){return d0_;}
float VertexCalculator::ip3d(){return ip3d_;}
bool VertexCalculator::valid(){return valid_;}

float VertexCalculator::getRotAngle(const TVector3& v1, const TVector3& v2){
  TVector3 v3 = v1.Cross(v2);
  v3.SetMag(1.0);
  return TMath::ACos(v3[2]);
}

TVector3 VertexCalculator::getRotAxis(const TVector3& v1, const TVector3& v2){
  if (v1[2]==0 && v2[2]==0)
    return TVector3(0,0,1.0);
  TVector3 v3 = v1.Cross(v2);
  v3.SetMag(1.0);
  TVector3 ez(0,0,1.0);
  return v3.Cross(ez);
}

double VertexCalculator::getPhi(const double mass, const double e1, const double e2){
  double cos = 1.0-mass*mass/(2.0*e1*e2);
  if (cos>=1 || cos<=-1)
    return TMath::Pi();
  return TMath::ACos(cos);
}

double VertexCalculator::getRadius(const TVector3& v1, const TVector3& v2, const double phi){
  double x1 = v1[0];
  double y1 = v1[1];
  double x2 = v2[0];
  double y2 = v2[1];
  return TMath::Sqrt(TMath::Power((x1-x2)/2.,2)+TMath::Power((y1-y2)/2.,2))/TMath::Sin(phi);
}

std::pair<double, double> VertexCalculator::getCenter(const TVector3& v1, const TVector3& v2, const double phi){
  double rad = getRadius(v1, v2, phi);
  double x1 = v1[0];
  double y1 = v1[1];
  double x2 = v2[0];
  double y2 = v2[1];
  double q = TMath::Sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
  double x3 = (x1+x2)/2.;
  double y3 = (y1+y2)/2.;
  double xp = x3 + TMath::Sqrt(rad*rad-(q*q/4.))*(y1-y2)/q;
  double yp = y3 + TMath::Sqrt(rad*rad-(q*q/4.))*(x2-x1)/q;
  double xm = x3 - TMath::Sqrt(rad*rad-(q*q/4.))*(y1-y2)/q;
  double ym = y3 - TMath::Sqrt(rad*rad-(q*q/4.))*(x2-x1)/q;
  std::pair<double, double> close, far;
  if (xp*xp+yp*yp < xm*xm+ym*ym){
    close.first = xp;
    close.second = yp;
    far.first = xm;
    far.second = ym;
  }
  else{
    close.first = xm;
    close.second = ym;
    far.first = xp;
    far.second = yp;
  }

  if (phi <= TMath::Pi()/2)
    return close;
  else
    return far;  
}

double VertexCalculator::getPt(const TVector3& vertex, const TVector3& v1, const TVector3& v2, const double e1, const double e2){
  TVector3 temp1 = v1-vertex;
  TVector3 temp2 = v2-vertex;
  double sinTheta1 = (vertex.Cross(temp1).Mag())/(vertex.Mag()*temp1.Mag())*std::copysign(1.0, vertex.Cross(temp1)[2]);
  double sinTheta2 = (vertex.Cross(v2).Mag())/(vertex.Mag()*temp2.Mag())*std::copysign(1.0, vertex.Cross(v2)[2]);
  double pt1 = e1 * sinTheta1;
  double pt2 = e2 * sinTheta2;
  return pt1 + pt2;
  
}
bool VertexCalculator::checkValid2D(const TVector3& vertex, const TVector3& v1, const TVector3& v2){
  double x = (v1[0]+v2[0])/2.;
  double y = (v1[1]+v2[1])/2.;
  double x1 = v1[0];
  double y1 = v1[1];
  double x2 = v2[0];
  double y2 = v2[1];

  double hypo = TMath::Sqrt(x*x+y*y);
  double leg1 = TMath::Sqrt(x1*x1+y1*y1);
  double leg2 = TMath::Sqrt(x2*x2+y2*y2);
  return ( (hypo-vertex.Mag())>0 && (leg1-vertex.Mag())>0 && (leg2-vertex.Mag())>0 );
  
}

    // Set the vertex, pt, and phi
void VertexCalculator::run(const TVector3& v1, const TVector3& v2, const double e1, const double e2, const double mass){


  double phi = TMath::Pi();
  double pt = 999;
  bool valid = false;
  TVector3 vertex(-999,-999,-999);
  


  setVertex(vertex);
  setPhi(phi);
  setPt(pt);
  setD0(TMath::Sqrt(vertex[0]*vertex[0]+vertex[1]*vertex[1]));
  setIp3d(TMath::Sqrt(vertex[0]*vertex[0]+vertex[1]*vertex[1]+vertex[2]*vertex[2]));
  setValid(false);

  TVector3 axis = getRotAxis(v1, v2);
  double theta = getRotAngle(v1, v2);
  phi = getPhi(mass, e1, e2);
  if (phi >= TMath::Pi()){
    return;
  }

  TVector3 temp1 = v1;
  TVector3 temp2 = v2;
  temp1.Rotate(theta, axis);
  temp2.Rotate(theta, axis);
  double radius = getRadius(temp1, temp2, phi);
  std::pair<double, double> center = getCenter(temp1, temp2, phi);

  if (radius >= TMath::Sqrt(center.first*center.first+center.second*center.second)+5){
    return;
  }

  int stepsPhi = 200;
  double deltaPhi = 2*TMath::Pi()/stepsPhi;

  std::pair<double, double> best;
  for (int i=0; i<stepsPhi; i++){
    double angle = i * deltaPhi;
    double x = center.first+radius*TMath::Cos(angle);
    double y = center.second+radius*TMath::Sin(angle);
    TVector3 temp(x,y,0.0);
    if (checkValid2D(temp, temp1, temp2)){
      double tempPt = getPt(temp, temp1, temp2, e1, e2);
      if (abs(tempPt) <= abs(pt)){
	pt = tempPt;
	vertex = temp;
	valid = true;
      }
    }
  }


  if (valid){
    vertex.Rotate(-theta, axis);
    setVertex(vertex);
    setPt(pt);
    setPhi(phi);
    setD0(TMath::Sqrt(vertex[0]*vertex[0]+vertex[1]*vertex[1]));
    setIp3d(TMath::Sqrt(vertex[0]*vertex[0]+vertex[1]*vertex[1]+vertex[2]*vertex[2]));
    setValid(valid);
  }
  return;
}

