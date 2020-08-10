#ifndef CMGTools_DisplacedDiPhotons_VertexCalculator_h
#define CMGTools_DisplacedDiPhotons_VertexCalculator_h

#include "TMath.h"
#include "TVector3.h"
#include <cmath>
#include <map>

namespace cmg {
   
  class VertexCalculator {
  private:
    TVector3 vertex_;
    float pt_;
    float phi_;
    float d0_;
    float ip3d_;
    bool valid_;
   
  public:
    VertexCalculator();
    ~VertexCalculator();

    void setVertex(const TVector3& vertex);
    void setPt(double);
    void setPhi(double);
    void setD0(double);
    void setIp3d(double);
    void setValid(bool);

    TVector3 vertex();
    float pt();
    float phi();
    float d0();
    float ip3d();
    bool valid();

    // Return angle between photon plane and x-y plane
    float getRotAngle(const TVector3& v1, const TVector3& v2);

    // Return axis of rotation between photon plane and x-y plane
    TVector3 getRotAxis(const TVector3& v1, const TVector3& v2);

    // Get phi between photons (inscribed angle for circle)
    double getPhi(const double mass, const double e1, const double e2);

    // Following methods are for 2D after rotating ECAL vectors

    //Get radius for circle of possible vertices
    double getRadius(const TVector3& v1, const TVector3& v2, const double phi);

    //Get center for circle of possible vertices
    std::pair<double, double> getCenter(const TVector3& v1, const TVector3& v2, const double phi);

    //Calculate unbalanced Pt for X->gg (should be ~0 for momentum conservation)
    double getPt(const TVector3& vertex, const TVector3& v1, const TVector3& v2, const double e1, const double e2);

    //In 2d (rotated) coords, check if point is valid
    bool checkValid2D(const TVector3& vertex, const TVector3& v1, const TVector3& v2);

    // Set the vertex, pt, and phi
    void run(const TVector3& v1, const TVector3& v2, const double e1, const double e2, const double mass);

  };
}
#endif
