#ifndef EIGENSOLVERSYM3X3_H
#define EIGENSOLVERSYM3X3_H

#include "CLHEP/Matrix/Matrix.h"
#include <algorithm> // for sorting eigenvalues
#include <limits> // For NaN

// debugging macros so we can pin down message provenance at a glance
#include <iostream>
#define DEBUG(x)              \
  std::cout << "DEBUG:\t" << "\t" << #x << "\t" << x << "\t"    \
  << __FUNCTION__ << "\t" << __FILE__ << "\t" << __LINE__ << std::endl;
#define DEBUG(x)

class EigenSolverSym3x3 {
public:
  EigenSolverSym3x3(CLHEP::HepMatrix A): _A(A) {DEBUG(0);}

  bool operator()(double root[3]) { DEBUG(0);
    double max = computeMaximumMagnitudeEntry(); DEBUG(max);
    if(max > 1) _A /= max; DEBUG(_A);
    bool status = computeRoots(root); DEBUG(status);
    if(max > 1) for(int i = 0; i != 3; i++) root[i] *= max; DEBUG(root[0]);DEBUG(root[1]);DEBUG(root[2]);
    return status;
  }

private:
  double computeMaximumMagnitudeEntry() { DEBUG(0);
    double max = 0.; DEBUG(max);
    for(int i = 0; i != _A.num_row(); i++) 
      for(int j = 0; j != _A.num_col(); j++) 
	if(fabs(_A[i][j]) > max) max = fabs(_A[i][j]);
    return max;
  }

  bool computeRoots(double root[3]) { DEBUG(0);
    const double inv3 = 1./3.; DEBUG(inv3);
    const double root3 = sqrt(3.); DEBUG(root3);
    double a00 = _A[0][0]; DEBUG(0);
    double a01 = _A[0][1]; DEBUG(0);
    double a02 = _A[0][2]; DEBUG(0);
    double a11 = _A[1][1]; DEBUG(0);
    double a12 = _A[1][2]; DEBUG(0);
    double a22 = _A[2][2]; DEBUG(0);
    double c0 = a00*a11*a22 + 2.0*a01*a02*a12 - a00*a12*a12 - a11*a02*a02 - a22*a01*a01; DEBUG(0);
    double c1 = a00*a11 - a01*a01 + a00*a22 - a02*a02 + a11*a22 - a12*a12; DEBUG(0);
    double c2 = a00 + a11 + a22; DEBUG(0);
    double c2Div3 = c2*inv3; DEBUG(0);
    double aDiv3 = (c1 - c2*c2Div3)*inv3; DEBUG(0);
    if(aDiv3 > 0.0) aDiv3 = 0.0; DEBUG(0);
    double mbDiv2 = 0.5*(c0 + c2Div3*(2.0*c2Div3*c2Div3 - c1)); DEBUG(0);
    double q = mbDiv2*mbDiv2 + aDiv3*aDiv3*aDiv3; DEBUG(0);
    if(q > 0.0) q = 0.0; DEBUG(0);
    double magnitude = -aDiv3 >= 0.? sqrt(-aDiv3): std::numeric_limits<double>::quiet_NaN(); DEBUG(0);
    double angle = -q >= 0.? atan2(sqrt(-q), mbDiv2)*inv3: std::numeric_limits<double>::quiet_NaN(); DEBUG(0);
    double cs = cos(angle); DEBUG(0);
    double sn = sin(angle); DEBUG(0);
    root[0] = c2Div3 + 2.0*magnitude*cs; DEBUG(0);
    root[1] = c2Div3 - magnitude*(cs + root3*sn); DEBUG(0);
    root[2] = c2Div3 - magnitude*(cs - root3*sn); DEBUG(0);

    // Sort the roots here to obtain root[0] > root[1] > root[2].
    std::sort(root, root+3, comp); DEBUG(0);
    return -aDiv3 >= 0. && -q >= 0.;
  }
  CLHEP::HepMatrix _A;

  struct Compare {
    bool operator()(double i, double j) { return i > j;}
  } comp;
};

#endif
